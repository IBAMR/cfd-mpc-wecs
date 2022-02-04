# cfd-mpc-wecs
A model predictive control (MPC) integrated multiphase immersed boundary (IB) framework for simulating wave energy converters (WEC) 

# Building IBAMR code
Refer https://ibamr.github.io/

# Steps to simulate the (CFD+MPC) solver
1) **Generate BEM_data.mat file:** Run the script bemio.m located in directory cfd-mpc-wecs -> BEMIO to generate BEM_data.mat file.
2) **Input file:** For the CFD solver the input file is located in cfd-mpc-wecs directory with the name input3d.cyl for 1st order regular water waves and input3d_irregwave.cyl for irregular waves. To load the BEM data and to input MPC parameters use load_mpc_paramters.m file. Currently, the input files are set to simulate the regular wave case with AR model given in Sec. 9.2 in the paper.
3) **Building the CFD code:** Use the command **make main3d** to compile the CFD code. (Set IBAMR directory path in the provided Makefile. For more details see the IBAMR website)
4) **Run the simulation:** Use the command **mpirun -np 32 ./main3d input3d.cyl** to run the simulation with 32 processors (number of processors can be adjusted here).
5) **Output data:** Results contain the files
                    * rbd.curve -> Contains the simulation time, WEC z displacement, and velocity.
                    * hydro_force_torque.curve -> Contains the simulation time and hydrodynamic forces acting on the WEC device.
                    * mpc_control_force.curve -> Contains the simulation time and the control force applied.

# Steps to simulate the (BEM+MPC) solver
1) **Generate BEM_data.mat file:** Method remains same as above.
2) **Input file:** Details remain same as above.
3) **Running the simulation:** The driving script for the BEM solver is vcyl_driving_script.m. Set the solver and wave parameters in this file and run the script.
4) **output Data:** Results are written in BEM_results.mat file. It includes simulation time, WEC displacement and velocity, wave excitation force and the control force applied.
