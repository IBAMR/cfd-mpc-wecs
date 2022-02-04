// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
// All rights reserved.
//
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>

#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>
#include <petscvec.h>

// Header for PETSC MATLAB engine functions
#include <petscmatlab.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <LocationIndexRobinBcCoefs.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvDiffPredictorCorrectorHierarchyIntegrator.h>
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/BrinkmanPenalizationRigidBodyDynamics.h>
#include <ibamr/FESurfaceDistanceEvaluator.h>
#include <ibamr/FirstOrderStokesWaveGenerator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/IBInterpolantHierarchyIntegrator.h>
#include <ibamr/IBInterpolantMethod.h>
#include <ibamr/IBLevelSetMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>
#include <ibamr/IrregularWaveBcCoef.h>
#include <ibamr/IrregularWaveGenerator.h>
#include <ibamr/RelaxationLSMethod.h>
#include <ibamr/SpongeLayerForceFunction.h>
#include <ibamr/StokesFirstOrderWaveBcCoef.h>
#include <ibamr/SurfaceTensionForceFunction.h>
#include <ibamr/WaveDampingFunctions.h>
#include <ibamr/WaveGenerationFunctions.h>
#include <ibamr/app_namespaces.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Application specific includes.
#include "FlowGravityForcing.h"
#include "GravityForcing.h"
#include "LSLocateGasInterface.h"
#include "LevelSetInitialCondition.h"
#include "SetFluidGasSolidDensity.h"
#include "SetFluidGasSolidViscosity.h"
#include "SetLSProperties.h"
#include "TagLSRefinementCells.h"

int coarsest_ln, max_finest_ln;
double dx, ds;

CylinderInterface cylinder;

// Struct to maintain the properties of the MPC interface
struct MPC_Interface
{
    int sample_size, sampling_interval;
    double dt_sampling;

    std::list<double> t_past;
    std::list<double> eta_A_past;

    double m_plus_Ainf = 0.0;
    double next_tuning_time = 0.0, old_tuning_time = 0.0, dt_controller, dt_cfd_current, mpc_start_time;
    double initial_position, current_position = 0.0, velocity = 0.0;
    double eta_A_current = 0.0, ls_val_initial = 0.0;

    PetscScalar F_control = 0.0, F_control_old = 0.0, deltaU = 0.0;

    PetscMatlabEngine mengine;
};
MPC_Interface* mpc;

// Struct to maintain the properties of the AR model interface
struct AR_Interface
{
    int sample_size, sampling_interval;
    double AR_dt_sampling;

    double mpc_start_time;
    double eta_A_current = 0.0;

    std::list<double> t_past;
    std::list<double> eta_A_past;
    std::list<double> Z_past;
};
AR_Interface* AR_interface;

// Struct to reset solid level set
struct SolidLevelSetResetter
{
    Pointer<IBInterpolantMethod> ib_interp_ops;
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
    Pointer<CellVariable<NDIM, double> > ls_solid_var;
    Pointer<BrinkmanPenalizationRigidBodyDynamics> bp_rbd;
};

void
reset_solid_level_set_callback_fcn(double current_time, double new_time, int /*cycle_num*/, void* ctx)
{
    SolidLevelSetResetter* resetter = static_cast<SolidLevelSetResetter*>(ctx);
    resetter->ib_interp_ops->copyEulerianDataToIntegrator(new_time);

    // Get the new centroid of the body
    const double dt = new_time - current_time;
    Eigen::Vector3d XCOM_current = resetter->bp_rbd->getCurrentCOMPosn();
    Eigen::Vector3d XCOM_new = XCOM_current + dt * (resetter->bp_rbd->getNewCOMTransVelocity());

    double distance[3]; // 3D cylinder has three surfaces.

    // Set a large value away from the solid body.
    Pointer<PatchHierarchy<NDIM> > patch_hier = resetter->adv_diff_integrator->getPatchHierarchy();
    const int hier_finest_ln = patch_hier->getFinestLevelNumber();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_solid_idx =
        var_db->mapVariableAndContextToIndex(resetter->ls_solid_var, resetter->adv_diff_integrator->getNewContext());

    for (int ln = 0; ln <= hier_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > patch_level = patch_hier->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* patch_X_lower = patch_geom->getXLower();
            const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
            const double* const patch_dx = patch_geom->getDx();

            Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(ls_solid_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                const hier::Index<NDIM>& ci = it();
                Eigen::Vector3d coord = Eigen::Vector3d::Zero();
                for (int d = 0; d < NDIM; ++d)
                {
                    coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                }

                // Distance from the circular surface.
                distance[0] =
                    std::sqrt(std::pow((coord[0] - XCOM_new(0)), 2.0) + std::pow((coord[1] - XCOM_new(1)), 2.0)) -
                    cylinder.R;

                distance[1] = coord[2] - (XCOM_new(2) + cylinder.L / 2);
                distance[2] = (XCOM_new(2) - cylinder.L / 2) - coord[2];

                (*ls_solid_data)(ci) = *std::max_element(distance, distance + 3);
            }
        }
    }

    return;
}

void
generate_interp_mesh(const unsigned int& /*strct_num*/,
                     const int& /*ln*/,
                     int& /*num_vertices*/,
                     std::vector<IBTK::Point>& /*vertex_posn*/)
{
    return;

} // generate_interp_mesh

void
imposed_kinematics(double /*data_time*/,
                   int /*cycle_num*/,
                   Eigen::Vector3d& U_com,
                   Eigen::Vector3d& W_com,
                   void* /*ctx*/)
{
    U_com.setZero();
    W_com.setZero();
    return;
} // imposed_kinematics

void
external_force_torque(double data_time, int cycle_num, Eigen::Vector3d& F, Eigen::Vector3d& T, void* /*ctx*/)
{
    double control_force = 0.0;

    if (SAMRAI_MPI::getRank() == 0)
    {
        PetscMatlabEngineEvaluate(mpc->mengine,
                                  "time_current = %f; dt_cfd = %f; z = %f; z_dot = %f;",
                                  (data_time - mpc->dt_cfd_current),
                                  mpc->dt_cfd_current,
                                  (mpc->current_position - mpc->initial_position),
                                  mpc->velocity);

        std::cout << "\n(z - z0) = " << (mpc->current_position - mpc->initial_position) << std::endl;

        if (cycle_num == 0)
        {
            if ((data_time - mpc->dt_cfd_current) >= mpc->next_tuning_time)
            {
                PetscMatlabEngineEvaluate(mpc->mengine, "calculate_mpc_matrices; get_radiation_damping_xr;");

                std::cout << "\nCalculating the control force using MPC." << std::endl;
                PetscMatlabEngineEvaluate(mpc->mengine, "F_control = %f;", mpc->F_control);
                mpc->F_control_old = mpc->F_control;

                // Store mpc past data in MATLAB workspace.
                std::vector<double> mpc_t_past(mpc->t_past.size());
                std::vector<double> mpc_eta_A(mpc->eta_A_past.size());

                std::copy(mpc->t_past.begin(), mpc->t_past.end(), mpc_t_past.begin());
                std::copy(mpc->eta_A_past.begin(), mpc->eta_A_past.end(), mpc_eta_A.begin());

                mpc_t_past.pop_back();
                mpc_eta_A.pop_back();

                mpc_t_past.push_back(data_time - mpc->dt_cfd_current);
                mpc_eta_A.push_back(mpc->eta_A_current);

                PetscMatlabEnginePutArray(mpc->mengine, mpc_t_past.size(), 1, &(mpc_t_past[0]), "t_past");
                PetscMatlabEnginePutArray(mpc->mengine, mpc_eta_A.size(), 1, &(mpc_eta_A[0]), "eta_A_past");

                // Store AR model past data in MATLAB workspace.
                std::vector<double> AR_t_past(AR_interface->t_past.size());
                std::vector<double> AR_eta_A(AR_interface->eta_A_past.size());
                std::vector<double> AR_Z_past(AR_interface->Z_past.size());

                std::copy(AR_interface->t_past.begin(), AR_interface->t_past.end(), AR_t_past.begin());
                std::copy(AR_interface->eta_A_past.begin(), AR_interface->eta_A_past.end(), AR_eta_A.begin());
                std::copy(AR_interface->Z_past.begin(), AR_interface->Z_past.end(), AR_Z_past.begin());

                AR_t_past.pop_back();
                AR_eta_A.pop_back();
                AR_Z_past.pop_back();

                AR_t_past.push_back(data_time - mpc->dt_cfd_current);
                AR_eta_A.push_back(mpc->eta_A_current);
                AR_Z_past.push_back(mpc->current_position);

                PetscMatlabEnginePutArray(mpc->mengine, AR_t_past.size(), 1, &(AR_t_past[0]), "AR_t_past");
                PetscMatlabEnginePutArray(mpc->mengine, AR_eta_A.size(), 1, &(AR_eta_A[0]), "AR_eta_A_past");
                PetscMatlabEnginePutArray(mpc->mengine, AR_Z_past.size(), 1, &(AR_Z_past[0]), "AR_z_past");

                // Calculate and get the control force.
                PetscMatlabEngineEvaluate(mpc->mengine, "get_control_force;");
                PetscMatlabEngineGetArray(mpc->mengine, 1, 1, &(mpc->deltaU), "deltaU_first");

                mpc->old_tuning_time = mpc->next_tuning_time;
                mpc->next_tuning_time = mpc->old_tuning_time + mpc->dt_controller;
            }
            mpc->F_control = mpc->F_control_old +
                             mpc->m_plus_Ainf * mpc->deltaU * (data_time - mpc->old_tuning_time) / mpc->dt_controller;
            control_force = mpc->F_control;
        }
        else
        {
            control_force = mpc->F_control;
        }
        std::cout << "control_force = " << control_force << std::endl;
        std::cout << "Next tuning time = " << mpc->next_tuning_time << std::endl;
    }

    int n = 1;
    SAMRAI_MPI::sumReduction(&control_force, n);

    F.setZero();
    F[2] = cylinder.mass * cylinder.g_z + control_force;
    T.setZero();

    return;
} // external_force_torque

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    // Increase maximum patch data component indices
    SAMRAIManager::setMaxNumberPatchDataEntries(2500);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IBLevelSet.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
        if (dump_restart_data && (restart_dump_interval > 0) && !restart_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(restart_dump_dirname);
        }

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Setup solid information
        cylinder.R = input_db->getDouble("R");
        cylinder.L = input_db->getDouble("L");
        cylinder.mass = input_db->getDouble("MASS");
        cylinder.X0[0] = input_db->getDouble("XCOM");
        cylinder.X0[1] = input_db->getDouble("YCOM");
#if (NDIM == 3)
        cylinder.X0[2] = input_db->getDouble("ZCOM");
#endif

        // Model Predictive Control information
        mpc = new MPC_Interface;
        AR_interface = new AR_Interface;

        if (SAMRAI_MPI::getRank() == 0)
        {
            std::cout << "\nStarting Matlab engine...\n" << std::endl;
            PetscMatlabEngineCreate(PETSC_COMM_SELF, "master", &(mpc->mengine));
            std::cout << "Matlab engine started" << std::endl;
            PetscMatlabEngineEvaluate(mpc->mengine, "clc;  clear all;  close all;  addpath([cd,'/MPC_matlab_code'])");

            std::cout << "Storing variables in Matlab workspace...\n" << std::endl;
            PetscMatlabEngineEvaluate(mpc->mengine,
                                      "global hydro wave mpc_interface; wave.H_wave = %f; wave.Tp = %f; wave.omega = "
                                      "%f; wave.wave_no = %f; wave.depth = %f; mass = %f; r_cyl = %f; L_cyl = %f;",
                                      input_db->getDouble("HEIGHT"),
                                      input_db->getDouble("TIME_PERIOD"),
                                      input_db->getDouble("OMEGA"),
                                      input_db->getDouble("WAVENUMBER"),
                                      input_db->getDouble("DEPTH"),
                                      cylinder.mass,
                                      cylinder.R,
                                      cylinder.L);

            PetscMatlabEngineEvaluate(mpc->mengine, "dt_solver = %f;", input_db->getDouble("DT_MAX"));

            PetscMatlabEngineEvaluate(mpc->mengine, "load_mpc_parameters;");

            mpc->initial_position = input_db->getDouble("Z_MEAN");
            mpc->current_position = mpc->initial_position;

            PetscScalar sample_size, dt_sampling, sampling_interval;
            PetscMatlabEngineGetArray(mpc->mengine, 1, 1, &sample_size, "sample_size");
            PetscMatlabEngineGetArray(mpc->mengine, 1, 1, &dt_sampling, "dt_sampling");
            //PetscMatlabEngineGetArray(mpc->mengine, 1, 1, &sampling_interval, "sampling_interval");
            mpc->sample_size = int(sample_size);
            mpc->dt_sampling = double(dt_sampling);
            //mpc->sampling_interval = int(sampling_interval);
            std::cout << "\n mpc_sample_size = " << mpc->sample_size << std::endl;
            std::cout << "\n mpc_dt_sampling = " << mpc->dt_sampling << std::endl;
            //std::cout << "mpc_sampling_interval = " << mpc->sampling_interval << std::endl;

            mpc->t_past.resize(mpc->sample_size);
            mpc->eta_A_past.resize(mpc->sample_size);

            PetscScalar AR_sample_size, AR_dt_sampling, AR_sampling_interval;
            PetscMatlabEngineGetArray(mpc->mengine, 1, 1, &AR_sample_size, "AR_sample_size");
            PetscMatlabEngineGetArray(mpc->mengine, 1, 1, &AR_dt_sampling, "AR_dt_sampling");
            //PetscMatlabEngineGetArray(mpc->mengine, 1, 1, &AR_sampling_interval, "AR_sampling_interval");
            AR_interface->sample_size = int(AR_sample_size);
            AR_interface->AR_dt_sampling = double(AR_dt_sampling);
            //AR_interface->sampling_interval = int(AR_sampling_interval);
            std::cout << "\n AR_sample_size = " << AR_interface->sample_size << std::endl;
            std::cout << "\n AR_dt_sampling = " << AR_interface->AR_dt_sampling << std::endl;
            //std::cout << "AR_sampling_interval = " << AR_interface->sampling_interval << std::endl;

            AR_interface->t_past.resize(AR_interface->sample_size);
            AR_interface->eta_A_past.resize(AR_interface->sample_size);
            AR_interface->Z_past.resize(AR_interface->sample_size);

            PetscMatlabEngineEvaluate(mpc->mengine, "hydro.x_B = %f; hydro.x_A = hydro.x_B - wave.df", cylinder.X0[0]);

            PetscMatlabEngineGetArray(mpc->mengine, 1, 1, &(mpc->m_plus_Ainf), "a");

            PetscMatlabEngineGetArray(mpc->mengine, 1, 1, &(mpc->dt_controller), "dt_controller");

            PetscMatlabEngineGetArray(mpc->mengine, 1, 1, &(mpc->mpc_start_time), "mpc_start_time");
            mpc->next_tuning_time = mpc->mpc_start_time;

            PetscMatlabEngineEvaluate(mpc->mengine, "z_initial = %f", mpc->initial_position);
            
            PetscMatlabEngineEvaluate(mpc->mengine, "X_0 = [  %f  ;  %f  ;  %f  ]", cylinder.X0[0], cylinder.X0[1], mpc->initial_position);
        }

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSVCStaggeredHierarchyIntegrator> navier_stokes_integrator;
        const string discretization_form =
            app_initializer->getComponentDatabase("Main")->getString("discretization_form");
        const bool conservative_form = (discretization_form == "CONSERVATIVE");
        if (conservative_form)
        {
            navier_stokes_integrator = new INSVCStaggeredConservativeHierarchyIntegrator(
                "INSVCStaggeredConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredConservativeHierarchyIntegrator"));
        }
        else if (!conservative_form)
        {
            navier_stokes_integrator = new INSVCStaggeredNonConservativeHierarchyIntegrator(
                "INSVCStaggeredNonConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredNonConservativeHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << discretization_form << "\n"
                                                   << "Valid options are: CONSERVATIVE, NON_CONSERVATIVE");
        }

        // Set up the advection diffusion hierarchy integrator
        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
        const string adv_diff_solver_type = app_initializer->getComponentDatabase("Main")->getStringWithDefault(
            "adv_diff_solver_type", "PREDICTOR_CORRECTOR");
        if (adv_diff_solver_type == "PREDICTOR_CORRECTOR")
        {
            Pointer<AdvectorExplicitPredictorPatchOps> predictor = new AdvectorExplicitPredictorPatchOps(
                "AdvectorExplicitPredictorPatchOps",
                app_initializer->getComponentDatabase("AdvectorExplicitPredictorPatchOps"));
            adv_diff_integrator = new AdvDiffPredictorCorrectorHierarchyIntegrator(
                "AdvDiffPredictorCorrectorHierarchyIntegrator",
                app_initializer->getComponentDatabase("AdvDiffPredictorCorrectorHierarchyIntegrator"),
                predictor);
        }
        else if (adv_diff_solver_type == "SEMI_IMPLICIT")
        {
            adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
                "AdvDiffSemiImplicitHierarchyIntegrator",
                app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << adv_diff_solver_type << "\n"
                                                   << "Valid options are: PREDICTOR_CORRECTOR, SEMI_IMPLICIT");
        }
        navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);

        Pointer<IBFEMethod> ibfe_method_ops = nullptr;
        Pointer<IBInterpolantMethod> ib_interpolant_method_ops = new IBInterpolantMethod(
            "IBInterpolantMethod", app_initializer->getComponentDatabase("IBInterpolantMethod"));
        Pointer<IBLevelSetMethod> ib_level_set_method_ops =
            new IBLevelSetMethod(ib_interpolant_method_ops, ibfe_method_ops);

        Pointer<IBHierarchyIntegrator> time_integrator = new IBInterpolantHierarchyIntegrator(
            "IBInterpolantHierarchyIntegrator",
            app_initializer->getComponentDatabase("IBInterpolantHierarchyIntegrator"),
            ib_level_set_method_ops,
            navier_stokes_integrator);

        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Create level sets for solid interface.
        const string& ls_name_solid = "level_set_solid";
        Pointer<CellVariable<NDIM, double> > phi_var_solid = new CellVariable<NDIM, double>(ls_name_solid);

        // Create level sets for gas/liquid interface.
        const double fluid_height = input_db->getDouble("GAS_LS_INIT");
        const string& ls_name_gas = "level_set_gas";
        Pointer<CellVariable<NDIM, double> > phi_var_gas = new CellVariable<NDIM, double>(ls_name_gas);
        Pointer<RelaxationLSMethod> level_set_gas_ops =
            new RelaxationLSMethod(ls_name_gas, app_initializer->getComponentDatabase("LevelSet_Gas"));
        LSLocateGasInterface* ptr_LSLocateGasInterface =
            new LSLocateGasInterface("LSLocateGasInterface", adv_diff_integrator, phi_var_gas, fluid_height);
        level_set_gas_ops->registerInterfaceNeighborhoodLocatingFcn(&callLSLocateGasInterfaceCallbackFunction,
                                                                    static_cast<void*>(ptr_LSLocateGasInterface));

        // Register the level sets with advection diffusion integrator.
        adv_diff_integrator->registerTransportedQuantity(phi_var_solid);
        adv_diff_integrator->setDiffusionCoefficient(phi_var_solid, 0.0);
        adv_diff_integrator->setAdvectionVelocity(phi_var_solid,
                                                  navier_stokes_integrator->getAdvectionVelocityVariable());

        adv_diff_integrator->registerTransportedQuantity(phi_var_gas);
        adv_diff_integrator->setDiffusionCoefficient(phi_var_gas, 0.0);
        adv_diff_integrator->setAdvectionVelocity(phi_var_gas,
                                                  navier_stokes_integrator->getAdvectionVelocityVariable());

        // Register the reinitialization functions for the level set variables
        SetLSProperties* ptr_setSetLSProperties = new SetLSProperties("SetLSProperties", level_set_gas_ops);
        adv_diff_integrator->registerResetFunction(
            phi_var_gas, &callSetGasLSCallbackFunction, static_cast<void*>(ptr_setSetLSProperties));

        // Solid level set initial conditions
        Pointer<CartGridFunction> phi_solid_init = new LevelSetInitialCondition("solid_ls_init", cylinder);
        adv_diff_integrator->setInitialConditions(phi_var_solid, phi_solid_init);

        SolidLevelSetResetter solid_level_set_resetter;
        solid_level_set_resetter.ib_interp_ops = ib_interpolant_method_ops;
        solid_level_set_resetter.adv_diff_integrator = adv_diff_integrator;
        solid_level_set_resetter.ls_solid_var = phi_var_solid;
        adv_diff_integrator->registerIntegrateHierarchyCallback(&reset_solid_level_set_callback_fcn,
                                                                static_cast<void*>(&solid_level_set_resetter));

        // Setup the advected and diffused fluid quantities.
        Pointer<CellVariable<NDIM, double> > mu_var = new CellVariable<NDIM, double>("mu");
        Pointer<hier::Variable<NDIM> > rho_var;
        if (conservative_form)
        {
            rho_var = new SideVariable<NDIM, double>("rho");
        }
        else
        {
            rho_var = new CellVariable<NDIM, double>("rho");
        }
        navier_stokes_integrator->registerMassDensityVariable(rho_var);
        navier_stokes_integrator->registerViscosityVariable(mu_var);

        // Array for input into callback function
        const int ls_reinit_interval = input_db->getInteger("LS_REINIT_INTERVAL");
        const double rho_fluid = input_db->getDouble("RHO_F");
        const double rho_solid = input_db->getDouble("RHO_S");
        const double rho_gas = input_db->getDouble("RHO_G");
        const int num_solid_interface_cells = input_db->getDouble("NUM_SOLID_INTERFACE_CELLS");
        const int num_gas_interface_cells = input_db->getDouble("NUM_GAS_INTERFACE_CELLS");
        cylinder.rho_solid = rho_solid;
        SetFluidGasSolidDensity* ptr_setFluidGasSolidDensity = new SetFluidGasSolidDensity("SetFluidGasSolidDensity",
                                                                                           adv_diff_integrator,
                                                                                           phi_var_solid,
                                                                                           phi_var_gas,
                                                                                           rho_fluid,
                                                                                           rho_gas,
                                                                                           rho_solid,
                                                                                           ls_reinit_interval,
                                                                                           num_solid_interface_cells,
                                                                                           num_gas_interface_cells);
        navier_stokes_integrator->registerResetFluidDensityFcn(&callSetFluidGasSolidDensityCallbackFunction,
                                                               static_cast<void*>(ptr_setFluidGasSolidDensity));

        const double mu_fluid = input_db->getDouble("MU_F");
        const double mu_gas = input_db->getDouble("MU_G");
        const double mu_solid = input_db->getDoubleWithDefault("MU_S", std::numeric_limits<double>::quiet_NaN());
        const bool set_mu_solid = input_db->getBool("SET_MU_S");
        SetFluidGasSolidViscosity* ptr_setFluidGasSolidViscosity =
            new SetFluidGasSolidViscosity("SetFluidGasSolidViscosity",
                                          adv_diff_integrator,
                                          phi_var_solid,
                                          phi_var_gas,
                                          mu_fluid,
                                          mu_gas,
                                          mu_solid,
                                          ls_reinit_interval,
                                          num_solid_interface_cells,
                                          num_gas_interface_cells,
                                          set_mu_solid);
        navier_stokes_integrator->registerResetFluidViscosityFcn(&callSetFluidGasSolidViscosityCallbackFunction,
                                                                 static_cast<void*>(ptr_setFluidGasSolidViscosity));

        // Register callback function for tagging refined cells for level set data
        const double tag_value = input_db->getDouble("LS_TAG_VALUE");
        const double tag_thresh = input_db->getDouble("LS_TAG_ABS_THRESH");
        TagLSRefinementCells ls_gas_tagger;
        ls_gas_tagger.d_ls_var = phi_var_gas;
        ls_gas_tagger.d_tag_value = tag_value;
        ls_gas_tagger.d_tag_abs_thresh = tag_thresh;
        ls_gas_tagger.d_adv_diff_solver = adv_diff_integrator;
        TagLSRefinementCells ls_solid_tagger;
        ls_solid_tagger.d_ls_var = phi_var_solid;
        ls_solid_tagger.d_tag_value = tag_value;
        ls_solid_tagger.d_tag_abs_thresh = tag_thresh;
        ls_solid_tagger.d_adv_diff_solver = adv_diff_integrator;
        time_integrator->registerApplyGradientDetectorCallback(&callTagGasLSRefinementCellsCallbackFunction,
                                                               static_cast<void*>(&ls_gas_tagger));
        time_integrator->registerApplyGradientDetectorCallback(&callTagSolidLSRefinementCellsCallbackFunction,
                                                               static_cast<void*>(&ls_solid_tagger));

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        const string& wave_type = input_db->getString("WAVE_TYPE");
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                // u_bc_coefs[d] = new muParserRobinBcCoefs(
                //    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
                if (wave_type == "FIRST_ORDER_STOKES")
                {
                    u_bc_coefs[d] = new StokesFirstOrderWaveBcCoef(
                        bc_coefs_name, d, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
                }
                else if (wave_type == "IRREGULAR")
                {
                    u_bc_coefs[d] = new IrregularWaveBcCoef(
                        bc_coefs_name, d, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
                }
                else
                {
                    TBOX_ERROR("Unknown WAVE_TYPE = " << wave_type << " specified in the input file" << std::endl);
                }
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        RobinBcCoefStrategy<NDIM>* rho_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("DensityBcCoefs"))
        {
            rho_bc_coef = new muParserRobinBcCoefs(
                "rho_bc_coef", app_initializer->getComponentDatabase("DensityBcCoefs"), grid_geometry);
            navier_stokes_integrator->registerMassDensityBoundaryConditions(rho_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* mu_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("ViscosityBcCoefs"))
        {
            mu_bc_coef = new muParserRobinBcCoefs(
                "mu_bc_coef", app_initializer->getComponentDatabase("ViscosityBcCoefs"), grid_geometry);
            navier_stokes_integrator->registerViscosityBoundaryConditions(mu_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* phi_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("PhiBcCoefs"))
        {
            phi_bc_coef = new muParserRobinBcCoefs(
                "phi_bc_coef", app_initializer->getComponentDatabase("PhiBcCoefs"), grid_geometry);
        }
        adv_diff_integrator->setPhysicalBcCoef(phi_var_gas, phi_bc_coef);
        adv_diff_integrator->setPhysicalBcCoef(phi_var_solid, phi_bc_coef);

        // LS reinit boundary conditions, which is set to be the same as the BCs
        // for advection
        RobinBcCoefStrategy<NDIM>* ls_reinit_bcs = phi_bc_coef;
        level_set_gas_ops->registerPhysicalBoundaryCondition(ls_reinit_bcs);

        // Create a damping zone near the channel outlet to absorb water waves via time-splitting approach.
        // This method uses a time splitting approach and modifies the fluid momentum in the
        // post processing step.
        // Register callback function for tagging refined cells for level set data
        const double outlet_zone_start = input_db->getDouble("OUTLET_ZONE_START");
        const double outlet_zone_end = input_db->getDouble("OUTLET_ZONE_END");
        const double depth = input_db->getDouble("DEPTH");
        const double alpha = input_db->getDouble("ALPHA");
        WaveDampingData wave_damper;
        wave_damper.d_x_zone_start = outlet_zone_start;
        wave_damper.d_x_zone_end = outlet_zone_end;
        wave_damper.d_depth = depth;
        wave_damper.d_alpha = alpha;
        wave_damper.d_sign_gas_phase = -1;
        wave_damper.d_ins_hier_integrator = navier_stokes_integrator;
        wave_damper.d_adv_diff_hier_integrator = adv_diff_integrator;
        wave_damper.d_phi_var = phi_var_gas;

        const string wave_damping_method = input_db->getStringWithDefault("WAVE_DAMPING_METHOD", "RELAXATION");
        if (wave_damping_method == "RELAXATION")
        {
            time_integrator->registerPostprocessIntegrateHierarchyCallback(
                &WaveDampingFunctions::callRelaxationZoneCallbackFunction, static_cast<void*>(&wave_damper));
        }
        else if (wave_damping_method == "CONSERVING")
        {
            time_integrator->registerPostprocessIntegrateHierarchyCallback(
                &WaveDampingFunctions::callConservedWaveAbsorbingCallbackFunction, static_cast<void*>(&wave_damper));
        }
        else
        {
            pout << "WARNING: Unknown WAVE_DAMPING_METHOD = " << wave_damping_method << " specified in the input file"
                 << std::endl;
            pout << "WARNING: Running without any wave damping method." << std::endl;
        }

        // Create a generating zone near the channel inlet to absorb water waves reflected by the
        // WEC towards the channel inlet. This method also uses a time splitting approach and modifies
        // the fluid momentum in the post processing step.
        // Register callback function for tagging refined cells for level set data
        const string& wave_generator_type = input_db->getStringWithDefault("WAVE_GENERATOR_TYPE", "");
        Pointer<Database> wave_db =
            app_initializer->getComponentDatabase("VelocityBcCoefs_0")->getDatabase("wave_parameters_db");
        StokesWaveGeneratorStrategy* wave_generator = nullptr;
        if (wave_generator_type == "FIRST_ORDER_STOKES")
        {
            wave_generator = new FirstOrderStokesWaveGenerator("FIRST_ORDER_STOKES", wave_db);
        }
        if (wave_generator_type == "IRREGULAR")
        {
            wave_generator = new IrregularWaveGenerator("IRREGULAR", wave_db);
        }
        if (wave_generator)
        {
            const double inlet_zone_start = input_db->getDouble("INLET_ZONE_START");
            const double inlet_zone_end = input_db->getDouble("INLET_ZONE_END");
            wave_generator->d_wave_gen_data.d_x_zone_start = inlet_zone_start;
            wave_generator->d_wave_gen_data.d_x_zone_end = inlet_zone_end;
            wave_generator->d_wave_gen_data.d_alpha = alpha;
            wave_generator->d_wave_gen_data.d_sign_gas_phase = -1;
            wave_generator->d_wave_gen_data.d_ins_hier_integrator = navier_stokes_integrator;
            wave_generator->d_wave_gen_data.d_adv_diff_hier_integrator = adv_diff_integrator;
            wave_generator->d_wave_gen_data.d_phi_var = phi_var_gas;

            time_integrator->registerPostprocessIntegrateHierarchyCallback(
                &WaveGenerationFunctions::callStokesWaveRelaxationCallbackFunction, static_cast<void*>(wave_generator));
        }

        // Body forces.
        std::vector<double> grav_const(NDIM);
        input_db->getDoubleArray("GRAV_CONST", &grav_const[0], NDIM);
        cylinder.g_z = grav_const[2];
        Pointer<CartGridFunction> grav_force;
        const string grav_type = input_db->getString("GRAV_TYPE");
        if (grav_type == "FULL")
        {
            grav_force = new GravityForcing("GravityForcing", navier_stokes_integrator, grav_const);
        }
        else if (grav_type == "FLOW")
        {
            grav_force = new FlowGravityForcing("FlowGravityForcing",
                                                app_initializer->getComponentDatabase("FlowGravityForcing"),
                                                adv_diff_integrator,
                                                phi_var_gas,
                                                grav_const);
        }

        Pointer<SurfaceTensionForceFunction> surface_tension_force =
            new SurfaceTensionForceFunction("SurfaceTensionForceFunction",
                                            app_initializer->getComponentDatabase("SurfaceTensionForceFunction"),
                                            adv_diff_integrator,
                                            phi_var_gas);

        Pointer<SpongeLayerForceFunction> sponge_fcn =
            new SpongeLayerForceFunction("SpongeLayerForceFunction",
                                         app_initializer->getComponentDatabase("SpongeLayerForceFunction"),
                                         navier_stokes_integrator,
                                         grid_geometry);

        Pointer<CartGridFunctionSet> eul_forces = new CartGridFunctionSet("eulerian_forces");
        eul_forces->addFunction(grav_force);
        eul_forces->addFunction(surface_tension_force);
        eul_forces->addFunction(sponge_fcn);
        time_integrator->registerBodyForceFunction(eul_forces);

        // Configure the IBInterpolant solver.
        Pointer<IBRedundantInitializer> ib_initializer = new IBRedundantInitializer(
            "IBRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
        std::vector<std::string> struct_list_vec(1, "InterpolationMesh");
        coarsest_ln = 0;
        max_finest_ln = input_db->getInteger("MAX_LEVELS") - 1;
        ib_initializer->setStructureNamesOnLevel(max_finest_ln, struct_list_vec);
        ib_initializer->registerInitStructureFunction(generate_interp_mesh);
        ib_interpolant_method_ops->registerLInitStrategy(ib_initializer);
        ib_interpolant_method_ops->registerVariableAndHierarchyIntegrator(
            ls_name_solid, /*depth*/ 1, phi_var_solid, adv_diff_integrator);

        // Configure the Brinkman penalization object to do the rigid body dynamics.
        Pointer<BrinkmanPenalizationRigidBodyDynamics> bp_rbd =
            new BrinkmanPenalizationRigidBodyDynamics("Brinkman Body",
                                                      phi_var_solid,
                                                      adv_diff_integrator,
                                                      navier_stokes_integrator,
                                                      app_initializer->getComponentDatabase("BrinkmanPenalization"),
                                                      /*register_for_restart*/ true);
        FreeRigidDOFVector free_dofs;
        free_dofs << 0, 0, 1, 0, 0, 0;        // surge, sway, heave, roll, pitch, yaw
        Eigen::Vector3d U_i = Eigen::Vector3d::Zero();
        const double mass = cylinder.mass;
        bp_rbd->setSolveRigidBodyVelocity(free_dofs);
        bp_rbd->registerKinematicsFunction(&imposed_kinematics);
        bp_rbd->registerExternalForceTorqueFunction(&external_force_torque);
        bp_rbd->setInitialConditions(cylinder.X0, U_i, U_i, mass);
        navier_stokes_integrator->registerBrinkmanPenalizationStrategy(bp_rbd);
        solid_level_set_resetter.bp_rbd = bp_rbd;

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        if (uses_visit)
        {
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
            ib_interpolant_method_ops->registerLSiloDataWriter(silo_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        ib_interpolant_method_ops->freeLInitStrategy();
        ib_initializer.setNull();
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();

        if (!RestartManager::getManager()->isFromRestart())
        {
            navier_stokes_integrator->initializeCompositeHierarchyData(loop_time, /*initial_time*/ true);
        }

        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(iteration_num, loop_time);
            }
        }

        // Get the probe points from the input file
        Pointer<Database> probe_db = input_db->getDatabase("ProbePoints");
        const int num_probes = (probe_db->getAllKeys()).getSize();
        std::vector<std::vector<double> > probe_points;
        std::vector<std::ofstream*> probe_streams;
        probe_points.resize(num_probes);
        probe_streams.resize(num_probes);
        for (int i = 0; i < num_probes; ++i)
        {
            std::string probe_name = "probe_" + Utilities::intToString(i);
            probe_points[i].resize(NDIM);
            probe_db->getDoubleArray(probe_name, &probe_points[i][0], NDIM);

            if (SAMRAI_MPI::getRank() == 0)
            {
                probe_streams[i] = new std::ofstream(probe_name.c_str(), std::fstream::out);

                *probe_streams[i] << "Printing level set at cell center closest to point (" << probe_points[i][0]
                                  << ", " << probe_points[i][1]
#if (NDIM == 3)
                                  << ", " << probe_points[i][2]
#endif
                                  << ") " << std::endl;
                probe_streams[i]->precision(10);
            }
        }

        // Open streams to save position and velocity of the structure.
        ofstream rbd_stream, ft_stream, CF_stream;
        if (SAMRAI_MPI::getRank() == 0)
        {
            rbd_stream.open("rbd.curve", ios_base::out | ios_base::app);
            ft_stream.open("hydro_force_torque.curve", ios_base::out | std::ios_base::app);
            CF_stream.open("mpc_control_force.curve", ios_base::out | std::ios_base::app);
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            mpc->dt_cfd_current = dt;
            pout << "Advancing hierarchy with timestep size dt = " << dt << "\n";

            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // Print out the level set values at probe locations
            // Note that it will print the nearest cell center
            // Max reduction over ls_val array to ensure that only processor 0 prints
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            int phi_idx = var_db->mapVariableAndContextToIndex(phi_var_gas, adv_diff_integrator->getCurrentContext());
            std::vector<double> ls_val(num_probes);
            for (int i = 0; i < num_probes; ++i)
            {
                ls_val[i] = -std::numeric_limits<double>::max();
                bool found_point_in_patch = false;
                for (int ln = max_finest_ln; ln >= coarsest_ln; --ln)
                {
                    // Get the cell index for this point
                    Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                    const CellIndex<NDIM> cell_idx =
                        IndexUtilities::getCellIndex(&probe_points[i][0], level->getGridGeometry(), level->getRatio());
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM> > patch = level->getPatch(p());
                        const Box<NDIM>& patch_box = patch->getBox();
                        const bool contains_probe = patch_box.contains(cell_idx);
                        if (!contains_probe) continue;

                        // Get the level set value at this particular cell and print to stream
                        Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(phi_idx);
                        ls_val[i] = (*phi_data)(cell_idx);
                        found_point_in_patch = true;
                        if (found_point_in_patch) break;
                    }
                    if (found_point_in_patch) break;
                }
            }
            SAMRAI_MPI::maxReduction(&ls_val[0], num_probes);
            if (SAMRAI_MPI::getRank() == 0)
            {
                for (int i = 0; i < num_probes; ++i)
                {
                    *probe_streams[i] << loop_time << '\t' << ls_val[i] << std::endl;
                }
                if (loop_time == dt) 
                {
                    mpc->ls_val_initial = ls_val[0];
                }
                // Storing the level set value at x_A (upwave) point.
                mpc->eta_A_current = -ls_val[0] + mpc->ls_val_initial;
            }

            // At specified intervals, write visualization and restart files,
            // and print out timer data.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "Writing visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(iteration_num, loop_time);
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "Writing restart files...\n\nn";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "Writing timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }

            if (SAMRAI_MPI::getRank() == 0)
            {
                const Eigen::Vector3d& rbd_posn = bp_rbd->getCurrentCOMPosn();
                const Eigen::Vector3d& rbd_trans_vel = bp_rbd->getCurrentCOMTransVelocity();

                rbd_stream.precision(12);
                rbd_stream.setf(ios::fixed, ios::floatfield);
                rbd_stream << loop_time << "\t" << rbd_posn[2] << "\t" << rbd_trans_vel[2] << std::endl;

                Eigen::Vector3d hydro_force_pressure, hydro_force_viscous, hydro_torque_pressure, hydro_torque_viscous;
                bp_rbd->getHydrodynamicForceTorque(
                    hydro_force_pressure, hydro_force_viscous, hydro_torque_pressure, hydro_torque_viscous);
                ft_stream.precision(12);
                ft_stream.setf(ios::fixed, ios::floatfield);
                ft_stream << loop_time << "\t" << hydro_force_viscous[0] << "\t" << hydro_force_viscous[1] << "\t"
                          << hydro_force_viscous[2] << "\t" << hydro_force_pressure[0] << "\t"
                          << hydro_force_pressure[1] << "\t" << hydro_force_pressure[2] << std::endl;

                CF_stream.precision(12);
                CF_stream.setf(ios::fixed, ios::floatfield);
                CF_stream << loop_time << "\t" << mpc->F_control << std::endl;

                // At specified mpc sampling dt write the elevation data in t_past and eta_A_past.
                if ((mpc->t_past.back() + mpc->dt_sampling) <= loop_time)
                {
                    mpc->t_past.pop_front();
                    mpc->eta_A_past.pop_front();

                    mpc->t_past.push_back(loop_time);
                    mpc->eta_A_past.push_back(mpc->eta_A_current);
                }

                // At specified AR sampling dt write the elevation data in t_past and eta_A_past.
                if ((AR_interface->t_past.back() + AR_interface->AR_dt_sampling) <= loop_time)
                {
                    AR_interface->t_past.pop_front();
                    AR_interface->eta_A_past.pop_front();
                    AR_interface->Z_past.pop_front();

                    AR_interface->t_past.push_back(loop_time);
                    AR_interface->eta_A_past.push_back(mpc->eta_A_current);
                    AR_interface->Z_past.push_back(rbd_posn[2]);
                }

                mpc->current_position = rbd_posn[2];
                mpc->velocity = rbd_trans_vel[2];
            }
        }

        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
            rbd_stream.close();
            ft_stream.close();
            CF_stream.close();
        }

        // Delete dumb pointers.
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
        delete ptr_setFluidGasSolidDensity;
        delete ptr_setFluidGasSolidViscosity;
        delete rho_bc_coef;
        delete mu_bc_coef;
        delete phi_bc_coef;
        delete wave_generator;
        delete mpc;
        delete AR_interface;
        for (int i = 0; i < num_probes; ++i) delete probe_streams[i];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
} // main
