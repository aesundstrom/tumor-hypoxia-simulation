tumor-hypoxia-simulation
========================

Matlab code for simulating tumor hypoxia and other cell-population-level phenomena, as featured in the manuscript "Simulating Heterogeneous Tumor Cell Populations".

There are four simulation configurations, corresponding to the four featured phenomena, and they are stored in the following subdirectories:

    1. metabolic_symbiosis_2d/
    2. stable_local_hypoxia_with_many_vessels_2d/
    3. stable_local_hypoxia_with_many_vessels_3d/
    4. tumor_stroma_signaling_2d/

Each of these subdirectories contains two configuration files:

    1. initialize_occupation_pattern.m    -- Initializes the cell occupation matrix with desired initial cell population.
    2. initialize_operating_parameters.m  -- Initializes the operating parameters for the simulation with desired values.

To run a simulation:

    1. In Matlab, change to this repository directory, or add it to your Matlab path.
    2. At the Matlab command line, type "simulate( '<simulation name>' )", where <simulation name> is one of the four subdirectory names above.
       e.g. "simulate( 'stable_local_hypoxia_with_many_vessels_3d' )"

Runtime explanation:

The simulator uses functions stored in the util/ subdirectory; in addition, the simulator has the following external dependencies:

    1. David Legland's imMinkowski module available from the Mathworks File Exchange:
       http://www.mathworks.com/matlabcentral/fileexchange/33690-geometric-measures-in-2d3d-images/all_files
       functions used: imEuler2d, imEuler3d
    
    2. John Iversen's freezeColors module available from the Mathworks File Exchange:
       http://www.mathworks.com/matlabcentral/fileexchange/7943-freezecolors---unfreezecolors
       functions used: freezeColors
    
    3. Carlos Adrian Vargas Aguilera's COLORMAP and COLORBAR utilities (Jul 2014) module available from the Mathworks File Exchange:
       http://www.mathworks.com/matlabcentral/fileexchange/24371-colormap-and-colorbar-utilities--jul-2014-
       functions used: cbfreeze, cbhandle
    
For convenience, I've placed a copy of the imEuler2d, imEuler3d, freezeColors, cbfreeze, and cbhandle functions in the external/ subdirectory.

The simulate.m function binds the simulator (simulator.m) to the util/ and external/ functions, and to the appropriate configuration data, by adding the util/ and external/ subdirectories, and the subdirectory named by the argument passed to simulate.m, respectively, to your Matlab path for the duration of the simulation; after the simulation terminates (by completion, ctrl-c interrupt, or an error) simulate.m removes the three subdirectories from your Matlab path.
