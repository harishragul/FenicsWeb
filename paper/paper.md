---
title: 'An Open-Source Application for Thermal Conduction Analysis: Integrating FEniCS with Django'
tags:
  - Python
  - FEniCS
  - FEM
  - Django
  - Thermal
authors:
  - name: Harish Ragul Rajaramaduraikarthik
    orcid: 0009-0007-5058-3995
    affiliation: 1
affiliations:
 - name: JJ College of Engineering and Technology, Trichy, India
   index: 1
date: 19 January 2025
---

# Summary

This project presents a powerful computation initiative for solving numerical analysis problems (like Finite Element Analysis and Computational Fluid Dynamics) using FEniCS integrated with Django. The application provides a web-based interface for engineers, researchers, and students to define, generate and analyze finite element models of thermal systems. For the initial phase of this work, I created a solver for 1D, 2D, and 3D Steady State Thermal Conduction Problems. The solver module employs the finite element method (FEM) to solve the governing heat conduction equation, incorporating user-defined thermal conductivity, boundary conditions, and heat source terms. Static visualization of both mesh structures and computed temperature distributions is provided using Matplotlib to get detailed insights into the thermal behavior of the system. The application stores all the data like mesh and solutions, which allows users for future reuse and retrieval. The application is designed to be open-source to serve as a foundational tool for educational purposes and research. The combination of FEniCS and Django bridges the gap between advanced numerical analysis and accessibility for non-experts.

# Statement of need

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References