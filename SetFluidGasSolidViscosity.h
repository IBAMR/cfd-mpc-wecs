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


/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_SetFluidGasSolidViscosity
#define included_SetFluidGasSolidViscosity

///////////////////////////// INCLUDES ///////////////////////////////////

#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <Variable.h>

namespace IBTK
{
class HierarchyMathOps;
}

/*!
 * Pre processing call back function to be hooked into IBAMR::VCINSStaggeredHierarchyIntegrator class.
 *
 * \param rho_idx a patch data index for the current density variable maintained by the integrator.
 * \param ctx is the pointer to SetFluidGasSolidViscosity class object.
 *
 * \TODO: Let's move this out of the global namespace and use "snake case" for static function names.
 */
void callSetFluidGasSolidViscosityCallbackFunction(int mu_idx,
                                                   SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var,
                                                   SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                   const int cycle_num,
                                                   const double time,
                                                   const double current_time,
                                                   const double new_time,
                                                   void* ctx);

/*!
 * \brief Class SetFluidGasSolidViscosity is a utility class which sets the fluid and
 * solid Eulerian density based on the current level set information
 */
class SetFluidGasSolidViscosity
{
public:
    /*!
     * The only constructor of this class.
     */
    SetFluidGasSolidViscosity(const std::string& object_name,
                              SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_solid_var,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_gas_var,
                              const double mu_fluid,
                              const double mu_gas,
                              const double mu_solid,
                              const int ls_reinit_interval,
                              const double num_solid_interface_cells,
                              const double num_gas_interface_cells,
                              const bool set_mu_solid);

    /*!
     * Destructor for this class.
     */
    ~SetFluidGasSolidViscosity() = default;

    /*!
     * Set the density based on the current level set information
     */
    void setViscosityPatchData(int mu_idx,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var,
                               SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                               const int cycle_num,
                               const double time,
                               const double current_time,
                               const double new_time);

private:
    /*!
     * Default constructor is not implemented and should not be used.
     */
    SetFluidGasSolidViscosity() = delete;

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    SetFluidGasSolidViscosity& operator=(const SetFluidGasSolidViscosity& that) = delete;

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    SetFluidGasSolidViscosity(const SetFluidGasSolidViscosity& from) = delete;

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Pointer to advection-diffusion solver.
     */
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;

    /*!
     * Level set variables
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_solid_var, d_ls_gas_var;

    /*!
     * Density of the fluid and solid.
     */
    double d_mu_fluid, d_mu_gas, d_mu_solid;

    /*!
     * Level set reinitialization interval
     */
    int d_ls_reinit_interval;

    /*!
     * Number of cells over which to transition between values
     */
    double d_num_solid_interface_cells, d_num_gas_interface_cells;

    /*!
     * Whether or not to set viscosity in the solid region
     */
    bool d_set_mu_solid;

}; // SetFluidGasSolidViscosity

#endif // #ifndef included_SetFluidGasSolidViscosity
