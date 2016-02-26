tumor-hypoxia-simulation
========================

Matlab code for simulating tumor hypoxia and other cell-population-level phenomena, as featured in the manuscript "Simulating Heterogeneous Tumor Cell Populations".

There are four functions, corresponding to the four featured phenomena -- the filenames are self-explanatory:

    simulate_metabolic_symbiosis_2d.m
    simulate_stable_local_hypoxia_with_many_vessels_2d.m
    simulate_stable_local_hypoxia_with_many_vessels_3d.m
    simulate_tumor_stroma_signaling_2d.m
    
These functions have the following dependencies:

    1. David Legland's imMinkowski module available from the Mathworks File Exchange:
       http://www.mathworks.com/matlabcentral/fileexchange/33690-geometric-measures-in-2d3d-images/all_files
       functions used: imEuler2d, imEuler3d
    
    2. John Iversen's freezeColors module available from the Mathworks File Exchange:
       http://www.mathworks.com/matlabcentral/fileexchange/7943-freezecolors---unfreezecolors
       functions used: freezeColors
    
    3. Carlos Adrian Vargas Aguilera's COLORMAP and COLORBAR utilities (Jul 2014) module available from the Mathworks File Exchange:
       http://www.mathworks.com/matlabcentral/fileexchange/24371-colormap-and-colorbar-utilities--jul-2014-
       functions used: cbfreeze, cbhandle
    
For convenience, I've placed a copy of the imEuler2d, imEuler3d, freezeColors, cbfreeze, and cbhandle functions in the code path.
