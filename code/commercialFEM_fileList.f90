!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief all DAMASK files without solver
!> @details List of files needed by MSC.Marc, Abaqus/Explicit, and Abaqus/Standard
!--------------------------------------------------------------------------------------------------
#include "IO.f90"
#include "libs.f90"
#include "numerics.f90"
#include "debug.f90"
#include "math.f90"
#include "FEsolving.f90"
#include "mesh.f90"
#include "material.f90"
#include "lattice.f90"
#include "source_thermal_dissipation.f90"
#include "source_thermal_externalheat.f90"
#include "source_damage_isoBrittle.f90"
#include "source_damage_isoDuctile.f90"
#include "source_damage_anisoBrittle.f90"
#include "source_damage_anisoDuctile.f90"
#include "source_vacancy_phenoplasticity.f90"
#include "source_vacancy_irradiation.f90"
#include "source_vacancy_thermalfluc.f90"
#include "kinematics_cleavage_opening.f90"
#include "kinematics_slipplane_opening.f90"
#include "kinematics_thermal_expansion.f90"
#include "kinematics_vacancy_strain.f90"
#include "kinematics_hydrogen_strain.f90"
#include "plastic_none.f90"
#include "plastic_j2.f90"
#include "plastic_phenopowerlaw.f90"
#include "plastic_titanmod.f90"
#include "plastic_dislotwin.f90"
#include "plastic_disloKMC.f90"
#include "plastic_disloUCLA.f90"
#include "plastic_nonlocal.f90"
#include "constitutive.f90"
#include "crystallite.f90"
#include "homogenization_none.f90"
#include "homogenization_isostrain.f90"
#include "homogenization_RGC.f90"
#include "thermal_isothermal.f90"
#include "thermal_adiabatic.f90"
#include "thermal_conduction.f90"
#include "damage_none.f90"
#include "damage_local.f90"
#include "damage_nonlocal.f90"
#include "vacancyflux_isoconc.f90"
#include "vacancyflux_isochempot.f90"
#include "vacancyflux_cahnhilliard.f90"
#include "porosity_none.f90"
#include "porosity_phasefield.f90"
#include "hydrogenflux_isoconc.f90"
#include "hydrogenflux_cahnhilliard.f90"
#include "homogenization.f90"
#include "CPFEM.f90"
