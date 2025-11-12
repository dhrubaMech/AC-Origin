**Title: Fluid flow induced biomechanical origin of collagen architecture in articular cartilage**  
**Authors: Dhruba Jyoti Mech and Mohd Suhail Rizvi**  
**Affliation: Indian Institute of Technology Hyderabad**  
**Article: https://doi.org/10.1101/2025.07.13.664559**  

# About this project
This project comprises the computational models and code used to investigate how the distinct, layered structure of collagen fibers in articular cartilage emerges. The central hypothesis is that joint movement-driven fluid flow within the cartilage guides the alignment of collagen fibers during their formation. The research employs two complementary computational approaches: a high-level continuum model to track fiber orientation patterns and a detailed 3D discrete fiber network model to simulate individual fiber behavior and deposition. These models were used to simulate cartilage development from embryogenesis to adulthood, its adaptation to different physical activities, and its degeneration in osteoarthritis, providing a unified biomechanical explanation for the tissue's structure, function, and pathology.

## Key features:
- **Dual-Model Approach**: Implements both a continuum orientation field model and a discrete 3D fiber network model for comprehensive multi-scale analysis.
- **Fluid-Flow Driven Mechanism**: The core of the simulation is a mechanobiological rule where synovial fluid flow, induced by joint movement, directs the preferential deposition and alignment of new collagen fibers.
- **Multi-Species & Multi-Joint Validation**: The code is designed to be validated against existing experimental data from various species, joint types, and developmental stages.
- **Activity-Dependent Adaptation**: Allows for systematic variation of mechanical loading (shear and compression durations) to simulate different physical activities and predict their impact on collagen architecture, tissue stiffness, and fluid viscosity.
- **Disease Modeling**: Includes functionality to simulate osteoarthritic remodeling by introducing localized disruptions to the normal fluid flow, demonstrating how this leads to progressive collagen disorganization.
- **Reproduction of Classic Architecture**: Successfully generates the well-known "Benninghoff" collagen architecture observed in mature, healthy articular cartilage.

## Goal:
The primary goal of this project is to establish a unifying biomechanical framework that explains:
- **Development**: How the zonal collagen architecture of articular cartilage emerges during embryogenesis and postnatal growth.
- **Function & Adaptation**: How the mature tissue's structure and mechanical properties (stiffness, fluid viscosity) adapt to specific, long-term mechanical loading patterns.
- **Degeneration**: How disruptions to the mechanical-fluid environment, as in osteoarthritis, lead to the progressive breakdown of the collagen network.

## Requirements
- [Numpy](https://numpy.org/)
- [Matplotlib](https://matplotlib.org/)
- [Numba](https://numba.pydata.org/)
- [PyVista](https://docs.pyvista.org/index.html)

You can install them using pip:
```bash
pip install numpy matplotlib numba pyvista[all]
```

## Usage
Make sure you have Python and dependent libraries installed. Then run the scripts as
```bash
python main.py
```
Use params.py file to change any simulation parameters.
The description of the key variables necesary for the simulation are provided in parameters file.

## License
This code is released for academic and research purposes only.  
Please cite the following paper if you use this code in your work:
> Modeling the role of ATP metabolism in articular cartilage and osteoarthritis\
> Dhruba Jyoti Mech and Mohd Suhail Rizvi\
> [Visit Pre-print in BioRxiv](https://doi.org/10.1101/2025.07.13.664559)
