!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief all DAMASK files without solver
!> @details List of files needed by MSC.Marc
!--------------------------------------------------------------------------------------------------
#include "parallelization.f90"
#include "IO.f90"
#include "YAML_types.f90"
#include "YAML_parse.f90"
#include "config.f90"
#include "LAPACK_interface.f90"
#include "math.f90"
#include "rotations.f90"
#include "element.f90"
#include "HDF5_utilities.f90"
#include "results.f90"
#include "geometry_plastic_nonlocal.f90"
#include "discretization.f90"
#include "Marc/discretization_Marc.f90"
#include "material.f90"
#include "lattice.f90"
#include "phase.f90"
#include "phase_mechanical.f90"
#include "phase_mechanical_plastic.f90"
#include "phase_mechanical_plastic_none.f90"
#include "phase_mechanical_plastic_isotropic.f90"
#include "phase_mechanical_plastic_phenopowerlaw.f90"
#include "phase_mechanical_plastic_kinehardening.f90"
#include "phase_mechanical_plastic_dislotwin.f90"
#include "phase_mechanical_plastic_dislotungsten.f90"
#include "phase_mechanical_plastic_nonlocal.f90"
#include "phase_mechanical_eigen.f90"
#include "phase_mechanical_eigen_cleavageopening.f90"
#include "phase_mechanical_eigen_slipplaneopening.f90"
#include "phase_mechanical_eigen_thermalexpansion.f90"
#include "phase_thermal.f90"
#include "phase_thermal_dissipation.f90"
#include "phase_thermal_externalheat.f90"
#include "phase_damage.f90"
#include "phase_damage_isobrittle.f90"
#include "phase_damage_isoductile.f90"
#include "phase_damage_anisobrittle.f90"
#include "homogenization.f90"
#include "homogenization_mechanical.f90"
#include "homogenization_mechanical_pass.f90"
#include "homogenization_mechanical_isostrain.f90"
#include "homogenization_mechanical_RGC.f90"
#include "homogenization_thermal.f90"
#include "homogenization_thermal_pass.f90"
#include "homogenization_thermal_isotemperature.f90"
#include "homogenization_damage.f90"
#include "homogenization_damage_pass.f90"
#include "CPFEM.f90"
