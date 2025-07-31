# FenicsWeb: Steady State Thermal Solver Web Application

## Introduction
FenicsWeb is a web-based application designed to provide an interactive platform for solving steady-state thermal conduction problems using the FEniCS finite element library. The project integrates scientific computing with a user-friendly interface, enabling users to generate meshes, solve PDEs, and visualize results directly from their browser. This report details the architecture, features, and implementation of FenicsWeb, developed as a college project.

## Objectives
- To create a web application for solving and visualizing steady-state thermal conduction problems.
- To provide interactive and static visualizations of solutions.
- To allow users to download solution data and reports for further analysis.
- To demonstrate the integration of Python scientific libraries with Django web framework.

## Technologies Used
- **Python**: Core programming language for backend logic and scientific computation.
- **Django**: Web framework for building the application, handling routing, forms, and templates.
- **FEniCS**: Open-source computing platform for solving partial differential equations (PDEs).
- **Bootstrap**: Frontend framework for responsive UI design.
- **HTML/CSS/JavaScript**: For templates, interactivity, and client-side features.

## Project Structure
- `FenicsWeb/`: Django project configuration and settings.
- `main/`, `mesh/`, `SteadyStateThermal/`: Django apps for core logic, mesh generation, and thermal solver.
- `static/`: Static files (CSS, JS, images) for frontend assets.
- `templates/`: HTML templates for rendering web pages.
- `db.sqlite3`: SQLite database for storing user and solution data.

## Key Features
### 1. Mesh Generation
- Users can generate computational meshes for the domain of interest.
- Meshes are visualized and displayed in the dashboard.

### 2. Steady-State Thermal Solver
- Users input boundary conditions and parameters via web forms.
- The backend uses FEniCS to solve the steady-state heat equation.
- Solution data (coordinates and temperature values) is computed and displayed.

### 3. Visualization
- **Interactive Plot**: Users can explore solution data interactively (zoom, pan, hover for values).
- **Static Plot**: High-quality static images for traditional visualization.
- **Plot Switching**: Users can toggle between interactive, static, or both plots side-by-side.

### 4. Solution Data Table
- Tabular display of solution data (coordinates and values).
- Search functionality to filter table rows by keywords.
- Download solution data as CSV for offline analysis.
- Download detailed PDF report for documentation.

## Implementation Details
### Django Backend
- Views handle user requests, form submissions, and solution computation.
- Models store mesh and solution data.
- Templates render dynamic content using Django's template language.

### Frontend Interactivity
- JavaScript functions for table filtering and CSV download:
  - `filterTable()`: Filters solution table rows based on user input.
  - `downloadCSV()`: Exports visible table data to a CSV file.
- Plot switching logic for interactive/static/both views.
- Bootstrap for responsive layout and styling.

### Data Export
- Users can download solution data in CSV format.
- PDF report generation is available for detailed documentation.

## Sample Code Snippets
### Table Filtering (JavaScript)
```javascript
function filterTable() {
    var input = document.getElementById("searchInput");
    var filter = input.value.toLowerCase();
    var table = document.getElementById("solutionTable");
    var tr = table.getElementsByTagName("tr");
    for (var i = 1; i < tr.length; i++) {
        tr[i].style.display = "none";
        var td = tr[i].getElementsByTagName("td");
        for (var j = 0; j < td.length; j++) {
            if (td[j]) {
                var txtValue = td[j].textContent || td[j].innerText;
                if (txtValue.toLowerCase().indexOf(filter) > -1) {
                    tr[i].style.display = "";
                    break;
                }
            }
        }
    }
}
```

### CSV Download (JavaScript)
```javascript
function downloadCSV() {
    var table = document.getElementById("solutionTable");
    var rows = table.getElementsByTagName("tr");
    var csvContent = "";
    for (var i = 0; i < rows.length; i++) {
        if (rows[i].style.display !== "none") {
            var cols = rows[i].getElementsByTagName("td");
            var rowContent = [];
            for (var j = 0; j < cols.length; j++) {
                rowContent.push(cols[j].innerText);
            }
            csvContent += rowContent.join(",") + "\n";
        }
    }
    var blob = new Blob([csvContent], { type: 'text/csv' });
    var url = URL.createObjectURL(blob);
    var a = document.createElement("a");
    a.setAttribute("hidden", "");
    a.setAttribute("href", url);
    a.setAttribute("download", "solution_data.csv");
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
}
```

## How to Use
1. Access the web application via browser.
2. Generate mesh and input problem parameters.
3. View solution plots and data table.
4. Use search and download features for data analysis and reporting.

## Conclusion
FenicsWeb demonstrates the power of combining scientific computing with modern web technologies. It provides an accessible interface for solving and visualizing thermal conduction problems, making it a valuable educational and research tool. The project showcases skills in Python, Django, frontend development, and scientific programming.

## References
- [FEniCS Project](https://fenicsproject.org/)
- [Django Documentation](https://docs.djangoproject.com/)
- [Bootstrap Documentation](https://getbootstrap.com/)

---
*Prepared by Harish Ragul for college project submission.*
