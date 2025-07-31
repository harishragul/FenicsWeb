from django.shortcuts import render, HttpResponse
from django.http import FileResponse
from .forms import HeatSolverForm, MeshForm
from SteadyStateThermal.conduction import solve_conduction, generate_solution_plot
from mesh.mesh import show_mesh, generate_mesh_plot
import io
import os
from datetime import datetime
from reportlab.lib.pagesizes import letter, A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY

# Create your views here.
def conduction(request):
    if request.method == "POST":
        form = HeatSolverForm(request.POST)
        if form.is_valid():
            dimension = form.cleaned_data['dimension']
            f = form.cleaned_data.get('f', 0)
            k = form.cleaned_data.get('k', 1)
            
            # Build boundary condition data structure
            bc_data = {}
            
            # Define boundary names based on dimension
            if dimension == '1D':
                boundary_names = ['left', 'right']
            elif dimension == '2D':
                boundary_names = ['left', 'right', 'top', 'bottom']
            elif dimension == '3D':
                boundary_names = ['left', 'right', 'top', 'bottom', 'front', 'back']
            
            # Process each boundary
            for boundary in boundary_names:
                bc_type = form.cleaned_data.get(f'{boundary}_bc_type')
                if bc_type:
                    bc_info = {'type': bc_type}
                    
                    if bc_type == 'dirichlet':
                        bc_value = form.cleaned_data.get(f'{boundary}_bc')
                        if bc_value is not None:
                            bc_info['value'] = bc_value
                        else:
                            # Skip this boundary if no value provided for Dirichlet BC
                            continue
                    elif bc_type == 'convection':
                        h = form.cleaned_data.get(f'{boundary}_h')
                        t_amb = form.cleaned_data.get(f'{boundary}_t_amb')
                        if h is not None and t_amb is not None:
                            bc_info['h'] = h
                            bc_info['t_amb'] = t_amb
                        else:
                            # Skip this boundary if convection parameters are missing
                            continue
                    
                    bc_data[boundary] = bc_info

            mesh_str = request.session['mesh_str']
            mesh = show_mesh(mesh_str)
            mesh_plot = generate_mesh_plot(mesh)
            
            try:
                solution, solution_list = solve_conduction(dimension, mesh, bc_data, f, k)
                solution_plots = generate_solution_plot(solution)

                # Store analysis data in session for PDF generation
                request.session['analysis_data'] = {
                    'dimension': dimension,
                    'f': f,
                    'k': k,
                    'bc_data': bc_data,
                    'solution_list': solution_list,
                    'mesh_plot': mesh_plot,
                    'solution_plots': solution_plots,
                    'timestamp': datetime.now().isoformat(),
                    'num_nodes': len(solution_list) if solution_list else 0,
                }

                context = {
                    'mesh_plot': mesh_plot,
                    'solution_plot_interactive': solution_plots['interactive'],
                    'solution_plot_static': solution_plots['static'],
                    'message': 'Problem Solved Successfully',
                    'solution_list': solution_list,
                }
                return render(request, 'conduction.html', context)
            except Exception as e:
                context = {
                    'mesh_plot': mesh_plot,
                    'message': f'Error solving problem: {str(e)}',
                    'form': HeatSolverForm(request.POST),
                }
                return render(request, 'conduction.html', context)
        else:
            form = MeshForm(request.POST)
            if form.is_valid():
                mesh = form.cleaned_data['mesh']
                request.session['mesh_str'] = mesh.mesh
                mesh = show_mesh(mesh.mesh)
                mesh_plot = generate_mesh_plot(mesh)

                context = {
                    'mesh_plot': mesh_plot,
                    'form': HeatSolverForm(),
                    'message': 'Enter the Boundary Conditions',
                }
                return render(request, 'conduction.html', context)
    else:
        form = MeshForm()
        context = {
            'form': form,
            'message': 'Select the Mesh to Solve',
        }
        return render(request, 'conduction.html', context)


def generate_pdf_report(request):
    """Generate a detailed PDF report of the thermal analysis."""
    
    # Check if analysis data exists in session
    analysis_data = request.session.get('analysis_data')
    if not analysis_data:
        return HttpResponse("No analysis data available. Please run an analysis first.", status=400)
    
    # Create a file-like buffer to receive PDF data
    buffer = io.BytesIO()
    
    # Create the PDF object with buffer as its 'file'
    doc = SimpleDocTemplate(buffer, pagesize=A4, 
                           rightMargin=72, leftMargin=72,
                           topMargin=72, bottomMargin=18)
    
    # Container for the 'Flowable' objects
    story = []
    
    # Define styles
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        'CustomTitle',
        parent=styles['Heading1'],
        fontSize=24,
        textColor=colors.HexColor('#2c3e50'),
        spaceAfter=30,
        alignment=TA_CENTER
    )
    
    heading_style = ParagraphStyle(
        'CustomHeading',
        parent=styles['Heading2'],
        fontSize=16,
        textColor=colors.HexColor('#34495e'),
        spaceAfter=12,
        spaceBefore=20
    )
    
    normal_style = ParagraphStyle(
        'CustomNormal',
        parent=styles['Normal'],
        fontSize=11,
        spaceAfter=12,
        alignment=TA_JUSTIFY
    )
    
    # Add title
    story.append(Paragraph("FenicsWeb - Thermal Analysis Report", title_style))
    story.append(Spacer(1, 12))
    
    # Add timestamp
    timestamp = datetime.fromisoformat(analysis_data['timestamp'])
    story.append(Paragraph(f"<b>Generated on:</b> {timestamp.strftime('%B %d, %Y at %I:%M %p')}", normal_style))
    story.append(Spacer(1, 20))
    
    # Problem Description
    story.append(Paragraph("Problem Description", heading_style))
    story.append(Paragraph(
        f"This report presents the results of a {analysis_data['dimension']} steady-state thermal conduction analysis "
        f"performed using the FEniCS finite element library. The analysis was conducted on a computational mesh "
        f"with {analysis_data['num_nodes']} nodes.",
        normal_style
    ))
    
    # Analysis Parameters
    story.append(Paragraph("Analysis Parameters", heading_style))
    
    param_data = [
        ['Parameter', 'Value', 'Units'],
        ['Dimension', analysis_data['dimension'], '-'],
        ['Heat Source (f)', str(analysis_data['f']), 'W/m³'],
        ['Thermal Conductivity (k)', str(analysis_data['k']), 'W/(m·K)'],
        ['Number of Nodes', str(analysis_data['num_nodes']), '-']
    ]
    
    param_table = Table(param_data, colWidths=[2.5*inch, 1.5*inch, 1.5*inch])
    param_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#3498db')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 12),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 1), (-1, -1), 10),
    ]))
    
    story.append(param_table)
    story.append(Spacer(1, 20))
    
    # Boundary Conditions
    story.append(Paragraph("Boundary Conditions", heading_style))
    
    bc_data_formatted = []
    bc_data_formatted.append(['Boundary', 'Type', 'Value', 'Additional Parameters'])
    
    for boundary, bc_info in analysis_data['bc_data'].items():
        bc_type = bc_info['type']
        if bc_type == 'dirichlet':
            bc_data_formatted.append([
                boundary.title(),
                'Dirichlet (Fixed Temperature)',
                f"{bc_info['value']} K",
                '-'
            ])
        elif bc_type == 'convection':
            bc_data_formatted.append([
                boundary.title(),
                'Convection',
                f"h = {bc_info['h']} W/(m²·K)",
                f"T_amb = {bc_info['t_amb']} K"
            ])
    
    bc_table = Table(bc_data_formatted, colWidths=[1.2*inch, 1.8*inch, 1.5*inch, 1.5*inch])
    bc_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#e74c3c')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 11),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 1), (-1, -1), 9),
    ]))
    
    story.append(bc_table)
    story.append(Spacer(1, 20))
    
    # Add mesh plot if available
    mesh_plot_path = os.path.join('static', analysis_data['mesh_plot'])
    if os.path.exists(mesh_plot_path):
        story.append(Paragraph("Computational Mesh", heading_style))
        story.append(Paragraph("The following figure shows the computational mesh used for the analysis:", normal_style))
        
        try:
            mesh_img = Image(mesh_plot_path, width=4*inch, height=3*inch)
            mesh_img.hAlign = 'CENTER'
            story.append(mesh_img)
            story.append(Spacer(1, 20))
        except:
            story.append(Paragraph("Mesh plot could not be loaded.", normal_style))
    
    # Add solution plot if available
    solution_plot_path = os.path.join('static', analysis_data['solution_plots']['static'])
    if os.path.exists(solution_plot_path):
        story.append(Paragraph("Solution Visualization", heading_style))
        story.append(Paragraph("The temperature distribution across the domain is shown below:", normal_style))
        
        try:
            solution_img = Image(solution_plot_path, width=4.5*inch, height=3.5*inch)
            solution_img.hAlign = 'CENTER'
            story.append(solution_img)
            story.append(Spacer(1, 20))
        except:
            story.append(Paragraph("Solution plot could not be loaded.", normal_style))
    
    # Results Summary
    story.append(Paragraph("Results Summary", heading_style))
    
    solution_list = analysis_data['solution_list']
    if solution_list:
        # Calculate statistics
        temperatures = [float(row[1]) for row in solution_list]
        min_temp = min(temperatures)
        max_temp = max(temperatures)
        avg_temp = sum(temperatures) / len(temperatures)
        
        summary_data = [
            ['Statistic', 'Value', 'Units'],
            ['Minimum Temperature', f"{min_temp:.3f}", 'K'],
            ['Maximum Temperature', f"{max_temp:.3f}", 'K'],
            ['Average Temperature', f"{avg_temp:.3f}", 'K'],
            ['Temperature Range', f"{max_temp - min_temp:.3f}", 'K'],
        ]
        
        summary_table = Table(summary_data, colWidths=[2.5*inch, 1.5*inch, 1.5*inch])
        summary_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#27ae60')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 12),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 1), (-1, -1), 10),
        ]))
        
        story.append(summary_table)
        story.append(Spacer(1, 20))
    
    # Detailed Solution Data (first 20 points)
    story.append(PageBreak())
    story.append(Paragraph("Detailed Solution Data", heading_style))
    story.append(Paragraph(
        f"The following table shows the temperature values at selected nodal coordinates. "
        f"Complete data for all {len(solution_list)} nodes is available through the web interface.",
        normal_style
    ))
    
    # Show first 20 data points
    display_data = [['Node', 'Coordinates', 'Temperature (K)']]
    for i, row in enumerate(solution_list[:20]):
        display_data.append([str(i+1), str(row[0]), f"{float(row[1]):.6f}"])
    
    if len(solution_list) > 20:
        display_data.append(['...', '...', '...'])
        display_data.append([str(len(solution_list)), str(solution_list[-1][0]), f"{float(solution_list[-1][1]):.6f}"])
    
    data_table = Table(display_data, colWidths=[0.8*inch, 2.5*inch, 1.7*inch])
    data_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#9b59b6')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 11),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 1), (-1, -1), 8),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
    ]))
    
    story.append(data_table)
    story.append(Spacer(1, 30))
    
    # Footer
    story.append(Paragraph("About FenicsWeb", heading_style))
    story.append(Paragraph(
        "FenicsWeb is a web-based finite element analysis platform powered by the FEniCS Project. "
        "This report was automatically generated to provide a comprehensive summary of the thermal analysis results. "
        "For more information or to perform additional analyses, please visit the FenicsWeb platform.",
        normal_style
    ))
    
    # Build PDF
    doc.build(story)
    
    # Get the value of the BytesIO buffer and write it to the response
    pdf = buffer.getvalue()
    buffer.close()
    
    # Create the HttpResponse object with the appropriate PDF headers
    response = HttpResponse(pdf, content_type='application/pdf')
    filename = f"thermal_analysis_report_{timestamp.strftime('%Y%m%d_%H%M%S')}.pdf"
    response['Content-Disposition'] = f'attachment; filename="{filename}"'
    
    return response