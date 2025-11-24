This code implements the methodology described in Section 3.1 of the manuscript “A Numerical Alternative to MR Thermometry for Safety Validation of Multi-Channel RF Transmit Coils.” To run the MATLAB script, three input files are required:

1. experimental_B1data;      -->    per-channel experimental B1+ data; 4D complex data with the format of [x,y,z,channels]
2. simulation_B1data;        -->    per-channel simulation B1+ data; 4D complex data with the format of [x,y,z,channels]
3. simulation_Q_matrices;    -->    10g-averaged local SAR data;  3D complex data with the size of [# of channels x # of channels x # of voxels]
