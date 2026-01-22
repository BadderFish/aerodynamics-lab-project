# Computational Aerodynamics Laboratory Project

**Institution:** Politecnico di Milano
**Program:** Aerospace Engineering (Master's)
**Date:** December 2025

## Overview

This repository contains a comprehensive computational aerodynamics study implementing classical panel methods and lifting-line theory for aerodynamic analysis. The project demonstrates proficiency in numerical methods, validation against industry-standard tools, and practical applications to airfoil and wing design.

## Project Structure

```
aerodynamics-lab-project/
├── panel_method_airfoil_analysis/     # 2D airfoil analysis using Hess-Smith method
├── weissinger_lifting_line_theory/    # 3D wing analysis using Weissinger method
├── additional_methods/                # Additional method implementations
├── results/                           # Generated plots and visualizations
├── Aerodynamics_Lab.pdf               # Complete laboratory report
└── README.md
```

## Technical Content

### 1. Panel Method for Airfoil Analysis

**Objective:** Implement and validate the Hess-Smith panel method for 2D incompressible flow analysis over airfoils.

#### Key Features:
- MATLAB implementation of the Hess-Smith panel method
- Integration with XFOIL for validation and comparison
- Analysis of NACA symmetric airfoils (0008, 0012)
- AG26 Bubble Dancer airfoil performance comparison
- Pressure coefficient (Cp) distribution calculations
- Boundary layer transition and separation point analysis

#### Results:
- **Validation accuracy:** <2% relative error compared to XFOIL
- **Reynolds number effects:** Investigated Re = 3×10⁵ to 5×10⁵
- **Panel discretization:** Up to 301 panels for high-fidelity results

#### Tools & Methods:
- **XFOIL:** Industry-standard airfoil analysis tool
- **Python:** Visualization and post-processing (matplotlib)
- **Bash scripting:** Automated XFOIL batch processing

### 2. Weissinger Lifting-Line Theory for Wings

**Objective:** Implement the Weissinger method for 3D finite-wing aerodynamic analysis with induced drag calculations.

#### Key Features:
- MATLAB implementation of iterative Weissinger lifting-line theory
- Validation against XFLR5 (3D panel/VLM code)
- Analysis of Cessna 172 Skyhawk wing configuration
- Comparison with custom Vimana aircraft design
- Lift and induced drag coefficient calculations
- Wing-tail interference effects

#### Results:
- **CL curve validation:** Excellent agreement with XFLR5 at low angles of attack
- **CDi predictions:** Sufficiently accurate for preliminary design
- **Aircraft comparison:** Quantified lift and drag differences between configurations

#### Tools & Methods:
- **XFLR5:** 3D wing design and analysis software
- **MATLAB:** Vortex lattice and control point calculations
- **Iterative solver:** Converged solutions for wing circulation distribution

## Key Findings

### Airfoil Analysis (Panel Method)
1. The Hess-Smith panel method achieves <2% error versus XFOIL for thin symmetric airfoils
2. Increasing Reynolds number advances transition location upstream (improves laminar flow extent)
3. AG26 Bubble Dancer generates ~73% more lift than NACA 0008 at α=3° due to cambered geometry

### Wing Analysis (Weissinger Method)
1. Weissinger method accuracy decreases at higher angles of attack (pre-stall nonlinearities)
2. Vimana aircraft has 18% higher CL but 67% higher CDi compared to Cessna 172
3. Tail downforce significantly affects total lift and induced drag distribution

## Computational Methods

### Hess-Smith Panel Method
- **Type:** Source + vortex panel method
- **Governing equations:** Laplace equation (incompressible, inviscid, irrotational flow)
- **Boundary condition:** Zero normal velocity at panel surfaces
- **Kutta condition:** Enforced at trailing edge for circulation determination

### Weissinger Lifting-Line Theory
- **Type:** Vortex lattice method with horseshoe vortices
- **Approach:** Iterative solution for spanwise circulation distribution
- **Induced drag:** Trefftz plane analysis
- **Limitations:** Valid for high aspect ratio wings, low angles of attack

## Usage

### Panel Method (Airfoil Analysis)

#### MATLAB Implementation:
```matlab
cd panel_method_airfoil_analysis/starting_point_HessSmith
main.m  % Run Hess-Smith solver
```

#### XFOIL Batch Processing:
```bash
cd panel_method_airfoil_analysis/example_script_bash
./run_xfoil.sh  # Automated XFOIL analysis
python plot_cp.py  # Generate Cp distribution plots
```

### Weissinger Method (Wing Analysis)

```matlab
cd weissinger_lifting_line_theory/Original\ Weissinger/Iterative\ Weissinger\ Script
Iterative_Weissinger.m  % Run lifting-line solver
```

## Requirements

### Software
- MATLAB R2020b or later
- Python 3.8+ with matplotlib, numpy
- XFOIL (included in example scripts)
- XFLR5 (optional, for comparison)

### MATLAB Toolboxes
- None required (base MATLAB only)

### Python Dependencies
```bash
pip install matplotlib numpy
```

## Results & Visualizations

See the `results/` directory for generated plots:
- Pressure coefficient distributions
- CL vs α curves
- CDi vs α curves
- MATLAB vs XFLR5 comparisons

Full analysis report: [`Aerodynamics_Lab.pdf`](Aerodynamics_Lab.pdf)

## Authors

- Arjun Nair
- Alessandro Pegoraro
- Alessandro Piana
- Francesco Poggi
- Kovindaraj Venkatesh

Politecnico di Milano - December 2025

## Acknowledgments

- Professor's starter code and validation datasets
- XFOIL developed by Mark Drela (MIT)
- XFLR5 developed by André Deperrois

## References

1. Katz, J., & Plotkin, A. (2001). *Low-Speed Aerodynamics* (2nd ed.). Cambridge University Press.
2. Drela, M. (1989). XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils.
3. Anderson, J. D. (2017). *Fundamentals of Aerodynamics* (6th ed.). McGraw-Hill Education.

---

*This project demonstrates practical application of computational aerodynamics methods for airfoil and wing analysis, with validation against industry-standard tools.*

