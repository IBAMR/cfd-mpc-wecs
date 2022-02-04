# cfd-mpc-wecs
A model predictive control (MPC) integrated multiphase immersed boundary (IB) framework for simulating wave energy converters (WEC) 

### Clone the repository

```
git clone https://github.com/IBAMR/cfd-mpc-wecs.git
```

# Build IBAMR 
Refer https://ibamr.github.io/

# Build PETSc with MATLAB
Refer https://petsc.org/main/docs/manual/matlab/ 

# Steps to execute the CFD+MPC solver:

1) **Generate BEM_data.mat file:** First excute the shell script `BEM_data/get_bemio_files.sh` to download some mising code files. It assumes the UNIX command `wget`is available on the system.  Next, execute the MATLAB script `BEM_data/process_bem_data.m`  to generate the BEM_data.mat file. The two scripts are located in the `BEM_data` directory. The directory also contains the output from ANSYS AQWA software for the 1:20 scaled cylinder: `AQWA_ANALYSIS.AH` and `AQWA_ANALYSIS.LIS` 

2) **Input files:** The CFD solver requires an input file. Two sample input files are provided in the main directory with the names `input3d.cyl` for first-order order regular waves and `input3d_irregwave.cyl` for irregular waves. The input files are set to simulate the AR-enabled CFD cases given in Sec. 9.2 of the paper.

3) **C++ driver code (main.cpp):** The BEM data and the MPC parameters are loaded by the `load_mpc_paramters.m` MATLAB script called within `main.cpp`. Similarly, the MATLAB MPC routines contained in the `MPC_matlab_code` directory are also called within the C++ driver code.   

3) **Building and linking the executable:** Use the command `make main3d` to compile the CFD code and link it with IBAMR. Modify the IBAMR source and build directory paths in the provided Makefile. For more details refer to https://ibamr.github.io/linking

4) **Run the simulation:** Use the command `mpirun -np 128 ./main3d input3d.cyl` to run the simulation with 128 processors (number of processors can be adjusted here). 

5) **Output data:** Results contain the files
                    * rbd.curve -> Contains the simulation time, WEC z displacement, and velocity.
                    * hydro_force_torque.curve -> Contains the simulation time and hydrodynamic forces acting on the WEC device.
                    * mpc_control_force.curve -> Contains the simulation time and the control force.



# Steps to execute the BEM+MPC solver:

1) **Generate BEM_data.mat file:** Method remains same as above.

2) **Input file:** Details remain same as above.

3) **Running the simulation:** The MATLAB driver script for the BEM solver is `vcyl_driving_script.m`. Set the solver and wave parameters in this file and run the script.

4) **output Data:** Results are saved in `BEM_results.mat` file. This includes the simulation time, WEC displacement and velocity, wave excitation and control force.
