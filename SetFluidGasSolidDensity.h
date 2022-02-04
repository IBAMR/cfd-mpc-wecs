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

#ifndef included_SetFluidGasSolidDensity
#define included_SetFluidGasSolidDensity

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
 * Pre processing call back function to be hooked into IBAMR::VCINSStaggeredHierarchyIntegratorclass.
 *
 * \param rho_idx a patch data index for the current density variable maintained by the integrator.
 * \param ctx is the pointer to SetFluidGasSolidDensity class object.
 *
 * \TODO: Let's move this out of the global namespace and use "snake case" for static function names.
 */
void callSetFluidGasSolidDensityCallbackFunction(int rho_idx,
                                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
                                                 SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                 const int cycle_num,
                                                 const double time,
                                                 const double current_time,
                                                 const double new_time,
                                                 void* ctx);

/*!
 * \brief Class SetFluidGasSolidDensity is a utility class which sets the fluid and
 * solid Eulerian density based on the current level set information
 */
class SetFluidGasSolidDensity
{
public:
    /*!
     * The only constructor of this class.
     */
    SetFluidGasSolidDensity(const std::string& object_name,
                            SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_solid_var,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_gas_var,
                            const double rho_fluid,
                            const double rho_gas,
                            const double rho_solid,
                            const int ls_reinit_interval,
                            const double num_solid_interface_cells,
                            const double num_gas_interface_cells);

    /*!
     * Destructor for this class.
     */
    ~SetFluidGasSolidDensity() = default;

    /*!
     * Set the density based on the current level set information
     */
    void setDensityPatchData(int rho_idx,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
                             SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                             const int cycle_num,
                             const double time,
                             const double current_time,
                             const double new_time);

private:
    /*!
     * Default constructor is not implemented and should not be used.
     */
    SetFluidGasSolidDensity() = delete;

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    SetFluidGasSolidDensity& operator=(const SetFluidGasSolidDensity& that) = delete;

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    SetFluidGasSolidDensity(const SetFluidGasSolidDensity& from) = delete;

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
    double d_rho_fluid, d_rho_gas, d_rho_solid;

    /*!
     * Level set reinitialization interval
     */
    int d_ls_reinit_interval;

    /*!
     * Number of interface cells
     */
    double d_num_solid_interface_cells, d_num_gas_interface_cells;

}; // SetFluidGasSolidDensity

#endif // #ifndef included_SetFluidGasSolidDensity
