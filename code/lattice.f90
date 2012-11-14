! Copyright 2011 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced MAterial Simulation Kit.
!
! DAMASK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DAMASK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DAMASK. If not, see <http://www.gnu.org/licenses/>.
!
!##############################################################
!* $Id$
!************************************
!*         Module: LATTICE          *
!************************************
!* contains:                        *
!* - Lattice structure definition   *
!* - Slip system definition         *
!* - Schmid matrices calculation    *
!************************************

module lattice

 use prec, only: pReal, &
                 pInt

 implicit none
 private
!************************************
!*      Lattice structures          *
!************************************

 integer(pInt), parameter, public :: &
   lattice_maxNslipFamily  =  5_pInt, &                                                             !< max # of slip system families over lattice structures
   lattice_maxNtwinFamily  =  4_pInt, &                                                             !< max # of twin system families over lattice structures
   lattice_maxNslip        = 54_pInt, &                                                             !< max # of slip systems over lattice structures
   lattice_maxNtwin        = 24_pInt, &                                                             !< max # of twin systems over lattice structures
   lattice_maxNinteraction = 30_pInt                                                                !< max # of interaction types (in hardening matrix part)
 
 integer(pInt), allocatable, dimension(:,:), protected, public :: &
   lattice_NslipSystem, &                                                                           !< # of slip systems in each family
   lattice_NtwinSystem                                                                              !< # of twin systems in each family

 integer(pInt), allocatable, dimension(:,:,:), protected, public :: &
   lattice_interactionSlipSlip, &                                                                   !< interaction type between slip/slip
   lattice_interactionSlipTwin, &                                                                   !< interaction type between slip/twin
   lattice_interactionTwinSlip, &                                                                   !< interaction type between twin/slip
   lattice_interactionTwinTwin                                                                      !< interaction type between twin/twin


 real(pReal), allocatable, dimension(:,:,:,:), protected, public :: &
   lattice_Sslip                                                                                    !< Schmid matrices, normal, shear direction and d x n of slip systems                                          
  
 real(pReal), allocatable, dimension(:,:,:), protected, public :: &
   lattice_Sslip_v, &
   lattice_sn, &
   lattice_sd, &
   lattice_st

! rotation and Schmid matrices, normal, shear direction and d x n of twin systems
 real(pReal), allocatable, dimension(:,:,:,:), protected, public  :: &
   lattice_Stwin, &
   lattice_Qtwin

 real(pReal), allocatable, dimension(:,:,:), protected, public :: &
   lattice_Stwin_v, &
   lattice_tn, &
   lattice_td, &
   lattice_tt

 real(pReal), allocatable, dimension(:,:), protected, public :: &
   lattice_shearTwin                                                                                !< characteristic twin shear

 integer(pInt), private :: &
   lattice_Nhexagonal, &                                                                            !< # of hexagonal lattice structure (from tag CoverA_ratio)
   lattice_Nstructure                                                                               !< # of lattice structures (1: fcc,2: bcc,3+: hexagonal)

 integer(pInt), dimension(:,:), pointer, private :: &
   interactionSlipSlip, &
   interactionSlipTwin, &
   interactionTwinSlip, &
   interactionTwinTwin

!--------------------------------------------------------------------------------------------------
! fcc (1)

 integer(pInt), dimension(lattice_maxNslipFamily), parameter, private :: & 
   lattice_fcc_NslipSystem = int([12, 0, 0, 0, 0],pInt)
   
 integer(pInt), dimension(lattice_maxNtwinFamily), parameter, private :: &
   lattice_fcc_NtwinSystem = int([12, 0, 0, 0],pInt)
   
 integer(pInt), parameter, private  :: &
   lattice_fcc_Nslip = 12_pInt, &                                                                   ! sum(lattice_fcc_NslipSystem)
   lattice_fcc_Ntwin = 12_pInt                                                                      ! sum(lattice_fcc_NtwinSystem)
 
 integer(pInt), private :: &
   lattice_fcc_Nstructure = 0_pInt

 real(pReal), dimension(3+3,lattice_fcc_Nslip), parameter, private :: &
   lattice_fcc_systemSlip = reshape(real([&
      0, 1,-1,     1, 1, 1, &
     -1, 0, 1,     1, 1, 1, &
      1,-1, 0,     1, 1, 1, &
      0,-1,-1,    -1,-1, 1, &
      1, 0, 1,    -1,-1, 1, &
     -1, 1, 0,    -1,-1, 1, &
      0,-1, 1,     1,-1,-1, &
     -1, 0,-1,     1,-1,-1, &
      1, 1, 0,     1,-1,-1, &
      0, 1, 1,    -1, 1,-1, &
      1, 0,-1,    -1, 1,-1, &
     -1,-1, 0,    -1, 1,-1  &
     ],pReal),[ 3_pInt + 3_pInt,lattice_fcc_Nslip])                                                 !< Slip system <110>{111}  Sorted according to Eisenlohr & Hantcherli

 real(pReal), dimension(3+3,lattice_fcc_Ntwin), parameter, private :: &
   lattice_fcc_systemTwin = reshape(real( [&
     -2, 1, 1,     1, 1, 1, &
      1,-2, 1,     1, 1, 1, &
      1, 1,-2,     1, 1, 1, &
      2,-1, 1,    -1,-1, 1, &
     -1, 2, 1,    -1,-1, 1, &
     -1,-1,-2,    -1,-1, 1, &
     -2,-1,-1,     1,-1,-1, &
      1, 2,-1,     1,-1,-1, &
      1,-1, 2,     1,-1,-1, &
      2, 1,-1,    -1, 1,-1, &
     -1,-2,-1,    -1, 1,-1, &
     -1, 1, 2,    -1, 1,-1  &
     ],pReal),[ 3_pInt + 3_pInt ,lattice_fcc_Ntwin])                                                !< Twin system <112>{111}  Sorted according to Eisenlohr & Hantcherli

 real(pReal), dimension(lattice_fcc_Ntwin), parameter, private :: &
   lattice_fcc_shearTwin = reshape([&
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal  &
     ],[lattice_fcc_Ntwin])                                                                         !< Twin system <112>{111}  Sorted according to Eisenlohr & Hantcherli

 integer(pInt), dimension(lattice_fcc_Nslip,lattice_fcc_Nslip), target, private :: &
   lattice_fcc_interactionSlipSlip = reshape(int( [&
     1,2,2,4,6,5,3,5,5,4,5,6, &  ! ---> slip
     2,1,2,6,4,5,5,4,6,5,3,5, &  ! |
     2,2,1,5,5,3,5,6,4,6,5,4, &  ! |
     4,6,5,1,2,2,4,5,6,3,5,5, &  ! v slip
     6,4,5,2,1,2,5,3,5,5,4,6, &
     5,5,3,2,2,1,6,5,4,5,6,4, &
     3,5,5,4,5,6,1,2,2,4,6,5, &
     5,4,6,5,3,5,2,1,2,6,4,5, &
     5,6,4,6,5,4,2,2,1,5,5,3, &
     4,5,6,3,5,5,4,6,5,1,2,2, &
     5,3,5,5,4,6,6,4,5,2,1,2, &
     6,5,4,5,6,4,5,5,3,2,2,1  &
     ],pInt),[lattice_fcc_Nslip,lattice_fcc_Nslip],order=[2,1])
    !< Interaction types
    !< 1 --- self interaction
    !< 2 --- coplanar interaction
    !< 3 --- colinear interaction
    !< 4 --- Hirth locks
    !< 5 --- glissile junctions
    !< 6 --- Lomer locks
    
 integer(pInt), dimension(lattice_fcc_Nslip,lattice_fcc_Ntwin), target, private :: &
   lattice_fcc_interactionSlipTwin = reshape(int( [&
     1,1,1,3,3,3,2,2,2,3,3,3, & ! ---> twin
     1,1,1,3,3,3,3,3,3,2,2,2, & ! |
     1,1,1,2,2,2,3,3,3,3,3,3, & ! |
     3,3,3,1,1,1,3,3,3,2,2,2, & ! v slip
     3,3,3,1,1,1,2,2,2,3,3,3, &
     2,2,2,1,1,1,3,3,3,3,3,3, &
     2,2,2,3,3,3,1,1,1,3,3,3, &
     3,3,3,2,2,2,1,1,1,3,3,3, &
     3,3,3,3,3,3,1,1,1,2,2,2, &
     3,3,3,2,2,2,3,3,3,1,1,1, &
     2,2,2,3,3,3,3,3,3,1,1,1, &
     3,3,3,3,3,3,2,2,2,1,1,1  &
     ],pInt),[lattice_fcc_Nslip,lattice_fcc_Ntwin],order=[2,1])
    !< Interaction types
    !< 1 --- coplanar interaction
    !< 2 --- colinear interaction
    !< 3 --- hardened interaction

 integer(pInt), dimension(lattice_fcc_Ntwin,lattice_fcc_Nslip), target, private :: &
   lattice_fcc_interactionTwinSlip = 0_pInt

 integer(pInt), dimension(lattice_fcc_Ntwin,lattice_fcc_Ntwin), target, private :: &
   lattice_fcc_interactionTwinTwin = reshape(int( [&
     1,1,1,2,2,2,2,2,2,2,2,2, &  ! ---> twin
     1,1,1,2,2,2,2,2,2,2,2,2, &  ! |
     1,1,1,2,2,2,2,2,2,2,2,2, &  ! |
     2,2,2,1,1,1,2,2,2,2,2,2, &  ! v twin
     2,2,2,1,1,1,2,2,2,2,2,2, &
     2,2,2,1,1,1,2,2,2,2,2,2, &
     2,2,2,2,2,2,1,1,1,2,2,2, &
     2,2,2,2,2,2,1,1,1,2,2,2, &
     2,2,2,2,2,2,1,1,1,2,2,2, &
     2,2,2,2,2,2,2,2,2,1,1,1, &
     2,2,2,2,2,2,2,2,2,1,1,1, &
     2,2,2,2,2,2,2,2,2,1,1,1  &
     ],pInt),[lattice_fcc_Ntwin,lattice_fcc_Ntwin],order=[2,1])
 
 
!--------------------------------------------------------------------------------------------------
! bcc (2)

 integer(pInt), dimension(lattice_maxNslipFamily), parameter, private :: &
   lattice_bcc_NslipSystem = int([ 12, 12, 0, 0, 0], pInt)
   
 integer(pInt), dimension(lattice_maxNtwinFamily), parameter, private :: &
   lattice_bcc_NtwinSystem = int([ 12, 0, 0, 0], pInt)
   
 integer(pInt), parameter, private :: &
   lattice_bcc_Nslip = 24_pInt                                       ! sum(lattice_bcc_NslipSystem)
   
 integer(pInt), parameter, private :: &
   lattice_bcc_Ntwin = 12_pInt                                      ! sum(lattice_bcc_NtwinSystem)
   
 integer(pInt), private :: &
   lattice_bcc_Nstructure = 0_pInt

 real(pReal), dimension(3+3,lattice_bcc_Nslip), parameter, private :: &
   lattice_bcc_systemSlip = reshape(real([&
    ! Slip system <111>{110} 
      1,-1, 1,     0, 1, 1, &
     -1,-1, 1,     0, 1, 1, &
      1, 1, 1,     0,-1, 1, &
     -1, 1, 1,     0,-1, 1, &
     -1, 1, 1,     1, 0, 1, &
     -1,-1, 1,     1, 0, 1, &
      1, 1, 1,    -1, 0, 1, &
      1,-1, 1,    -1, 0, 1, &
     -1, 1, 1,     1, 1, 0, &
     -1, 1,-1,     1, 1, 0, &
      1, 1, 1,    -1, 1, 0, &
      1, 1,-1,    -1, 1, 0, &
    ! Slip system <111>{112}
     -1, 1, 1,     2, 1, 1, &
      1, 1, 1,    -2, 1, 1, &
      1, 1,-1,     2,-1, 1, &
      1,-1, 1,     2, 1,-1, &
      1,-1, 1,     1, 2, 1, &
      1, 1,-1,    -1, 2, 1, &
      1, 1, 1,     1,-2, 1, &
     -1, 1, 1,     1, 2,-1, &
      1, 1,-1,     1, 1, 2, &
      1,-1, 1,    -1, 1, 2, &
     -1, 1, 1,     1,-1, 2, &
      1, 1, 1,     1, 1,-2 &
    ! Slip system <111>{123}
      ! 1, 1,-1,     1, 2, 3, &
      ! 1,-1, 1,    -1, 2, 3, &
     ! -1, 1, 1,     1,-2, 3, &
      ! 1, 1, 1,     1, 2,-3, &
      ! 1,-1, 1,     1, 3, 2, &
      ! 1, 1,-1,    -1, 3, 2, &
      ! 1, 1, 1,     1,-3, 2, &
     ! -1, 1, 1,     1, 3,-2, &
      ! 1, 1,-1,     2, 1, 3, &
      ! 1,-1, 1,    -2, 1, 3, &
     ! -1, 1, 1,     2,-1, 3, &
      ! 1, 1, 1,     2, 1,-3, &
      ! 1,-1, 1,     2, 3, 1, &
      ! 1, 1,-1,    -2, 3, 1, &
      ! 1, 1, 1,     2,-3, 1, &
     ! -1, 1, 1,     2, 3,-1, &
     ! -1, 1, 1,     3, 1, 2, &
      ! 1, 1, 1,    -3, 1, 2, &
      ! 1, 1,-1,     3,-1, 2, &
      ! 1,-1, 1,     3, 1,-2, &
     ! -1, 1, 1,     3, 2, 1, &
      ! 1, 1, 1,    -3, 2, 1, &
      ! 1, 1,-1,     3,-2, 1, &
      ! 1,-1, 1,     3, 2,-1  &
     ],pReal),[ 3_pInt + 3_pInt ,lattice_bcc_Nslip])

 real(pReal), dimension(3+3,lattice_bcc_Ntwin), parameter, private :: &
   lattice_bcc_systemTwin = reshape(real([&
    ! Twin system <111>{112}
     -1, 1, 1,     2, 1, 1, &
      1, 1, 1,    -2, 1, 1, &
      1, 1,-1,     2,-1, 1, &
      1,-1, 1,     2, 1,-1, &
      1,-1, 1,     1, 2, 1, &
      1, 1,-1,    -1, 2, 1, &
      1, 1, 1,     1,-2, 1, &
     -1, 1, 1,     1, 2,-1, &
      1, 1,-1,     1, 1, 2, &
      1,-1, 1,    -1, 1, 2, &
     -1, 1, 1,     1,-1, 2, &
      1, 1, 1,     1, 1,-2  &
     ],pReal),[ 3_pInt + 3_pInt,lattice_bcc_Ntwin])

 real(pReal), dimension(lattice_bcc_Ntwin), parameter, private :: &
   lattice_bcc_shearTwin = reshape([&
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal, &
     0.7071067812_pReal  &
     ],[lattice_bcc_Ntwin])

!> slip--slip interactions for BCC structures (2)   from Lee et-al. Int J of Plast. (v15) 1999 pp. 625-645
 integer(pInt), dimension(lattice_bcc_Nslip,lattice_bcc_Nslip), target, private :: &
   lattice_bcc_interactionSlipSlip = reshape(int( [&
     1,3,6,6,5,4,4,2,4,2,5,4,6,6,4,2,2,4,6,6,4,2,6,6, &  ! ---> slip  
     3,1,6,6,4,2,5,4,5,4,4,2,6,6,2,4,4,2,6,6,2,4,6,6, &  ! |
     6,6,1,3,4,5,2,4,4,5,2,4,4,2,6,6,6,6,2,4,6,6,4,2, &  ! |
     6,6,3,1,2,4,4,5,2,4,4,5,2,4,6,6,6,6,4,2,6,6,2,4, &  ! v slip
     5,4,4,2,1,3,6,6,2,4,5,4,2,6,4,6,6,4,6,2,4,6,2,6, &
     4,2,5,4,3,1,6,6,4,5,4,2,4,6,2,6,6,2,6,4,2,6,4,6, &
     4,5,2,4,6,6,1,3,5,4,2,4,6,2,6,4,4,6,2,6,6,4,6,2, &
     2,4,4,5,6,6,3,1,4,2,4,5,6,4,6,2,2,6,4,6,6,2,6,4, &
     4,5,4,2,2,4,5,4,1,3,6,6,2,6,6,4,4,6,6,2,6,4,2,6, &
     2,4,5,4,4,5,4,2,3,1,6,6,4,6,6,2,2,6,6,4,6,2,4,6, &
     5,4,2,4,5,4,2,4,6,6,1,3,6,2,4,6,6,4,2,6,4,6,6,2, &
     4,2,4,5,4,2,4,5,6,6,3,1,6,4,2,6,6,2,4,6,2,6,6,4, &
     6,6,4,2,2,4,6,6,2,4,6,6,1,5,6,6,5,6,6,2,5,6,2,6, &
     6,6,2,4,6,6,2,4,6,6,2,4,5,1,6,6,6,5,2,6,6,5,6,2, &
     4,2,6,6,4,2,6,6,6,6,4,2,6,6,1,5,6,2,5,6,2,6,5,6, &
     2,4,6,6,6,6,4,2,4,2,6,6,6,6,5,1,2,6,6,5,6,2,6,5, &
     2,4,6,6,6,6,4,2,4,2,6,6,5,6,6,2,1,6,5,6,5,2,6,6, &
     4,2,6,6,4,2,6,6,6,6,4,2,6,5,2,6,6,1,6,5,2,5,6,6, &
     6,6,2,4,6,6,2,4,6,6,2,4,6,2,5,6,5,6,1,6,6,6,5,2, &
     6,6,4,2,2,4,6,6,2,4,6,6,2,6,6,5,6,5,6,1,6,6,2,5, &
     4,2,6,6,4,2,6,6,6,6,4,2,5,6,2,6,5,2,6,6,1,6,6,5, &
     2,4,6,6,6,6,4,2,4,2,6,6,6,5,6,2,2,5,6,6,6,1,5,6, &
     6,6,4,2,2,4,6,6,2,4,6,6,2,6,5,6,6,6,5,2,6,5,1,6, &
     6,6,2,4,6,6,2,4,6,6,2,4,6,2,6,5,6,6,2,5,5,6,6,1  &
     ],pInt),[lattice_bcc_Nslip,lattice_bcc_Nslip],order=[2,1])
    !< Interaction types
    !< 1 --- self interaction
    !< 2 --- no interaction
    !< 3 --- coplanar interaction
    !< 4 --- glissile interaction
    !< 5 --- weak sessile interaction
    !< 6 --- strong sessile interaction
 
 integer(pInt), dimension(lattice_bcc_Nslip,lattice_bcc_Ntwin), target, private :: &
  lattice_bcc_interactionSlipTwin = reshape(int( [&
     3,3,3,2,2,3,3,3,3,2,3,3, &  ! ---> twin
     3,3,2,3,3,2,3,3,2,3,3,3, &  ! |
     3,2,3,3,3,3,2,3,3,3,3,2, &  ! |
     2,3,3,3,3,3,3,2,3,3,2,3, &  ! v slip 
     2,3,3,3,3,3,3,2,3,3,2,3, &
     3,3,2,3,3,2,3,3,2,3,3,3, &
     3,2,3,3,3,3,2,3,3,3,3,2, &
     3,3,3,2,2,3,3,3,3,2,3,3, &
     2,3,3,3,3,3,3,2,3,3,2,3, &
     3,3,3,2,2,3,3,3,3,2,3,3, &
     3,2,3,3,3,3,2,3,3,3,3,2, &
     3,3,2,3,3,2,3,3,2,3,3,3, &
     1,3,3,3,3,3,3,2,3,3,2,3, &
     3,1,3,3,3,3,2,3,3,3,3,2, &
     3,3,1,3,3,2,3,3,2,3,3,3, &
     3,3,3,1,2,3,3,3,3,2,3,3, &
     3,3,3,2,1,3,3,3,3,2,3,3, &
     3,3,2,3,3,1,3,3,2,3,3,3, &
     3,2,3,3,3,3,1,3,3,3,3,2, &
     2,3,3,3,3,3,3,1,3,3,2,3, &
     3,3,2,3,3,2,3,3,1,3,3,3, &
     3,3,3,2,2,3,3,3,3,1,3,3, &
     2,3,3,3,3,3,3,2,3,3,1,3, &
     3,2,3,3,3,3,2,3,3,3,3,1 &
     ],pInt),[lattice_bcc_Nslip,lattice_bcc_Ntwin],order=[2,1])
    !< Interaction types
    !< 1 --- coplanar interaction
    !< 2 --- colinear interaction
    !< 3 --- hardened interaction

!>twin--slip interactions for BCC structures (2) MISSING: not implemented yet
 integer(pInt), dimension(lattice_bcc_Ntwin,lattice_bcc_Nslip), target, private :: &
   lattice_bcc_interactionTwinSlip = 0_pInt

!> twin-twin interactions for BCC structures (2) MISSING: not implemented yet
 integer(pInt), dimension(lattice_bcc_Ntwin,lattice_bcc_Ntwin), target, private :: &
   lattice_bcc_interactionTwinTwin = 0_pInt
 

!--------------------------------------------------------------------------------------------------
! hex (3+)

 integer(pInt), dimension(lattice_maxNslipFamily), parameter, private :: &
   lattice_hex_NslipSystem = int([ 3, 3, 6,12, 6],pInt)
   
 integer(pInt), dimension(lattice_maxNtwinFamily), parameter, private :: &
   lattice_hex_NtwinSystem = int([ 6, 6, 6, 6],pInt)
   
 integer(pInt), parameter , private :: &
   lattice_hex_Nslip = 30_pInt                                       ! sum(lattice_hex_NslipSystem)
 
 integer(pInt), parameter, private :: &
   lattice_hex_Ntwin = 24_pInt                                       ! sum(lattice_hex_NtwinSystem)
   
 integer(pInt), private :: &
   lattice_hex_Nstructure = 0_pInt

 !* sorted by A. Alankar & P. Eisenlohr
 real(pReal), dimension(4+4,lattice_hex_Nslip), parameter, private :: &
   lattice_hex_systemSlip = reshape(real([&
    ! Basal systems <1120>{0001} (independent of c/a-ratio, Bravais notation (4 coordinate base))
      2, -1, -1,  0,     0,  0,  0,  1, &
     -1,  2, -1,  0,     0,  0,  0,  1, &
     -1, -1,  2,  0,     0,  0,  0,  1, &
    ! 1st type prismatic systems <11.0>{10.0}  (independent of c/a-ratio)
      2, -1, -1,  0,     0,  1, -1,  0, &
     -1,  2, -1,  0,    -1,  0,  1,  0, &
     -1, -1,  2,  0,     1, -1,  0,  0, &
    ! 1st type 1st order pyramidal systems <11.0>{10.1} -- plane normals depend on the c/a-ratio
      2, -1, -1,  0,     0,  1, -1,  1, &
      1,  1, -2,  0,    -1,  1,  0,  1, &
     -1,  2, -1,  0,    -1,  0,  1,  1, &
     -2,  1,  1,  0,     0, -1,  1,  1, &
     -1, -1,  2,  0,     1, -1,  0,  1, &
      1, -2,  1,  0,     1,  0, -1,  1, &
    ! pyramidal system: c+a slip <-21.-3>{-10.1} -- plane normals depend on the c/a-ratio
     -1,  2, -1, -3,     0,  1, -1,  1, &
      1,  1, -2, -3,     0,  1, -1,  1, &
     -2,  1,  1, -3,    -1,  1,  0,  1, &
     -1,  2, -1, -3,    -1,  1,  0,  1, &
     -1, -1,  2, -3,    -1,  0,  1,  1, &
     -2,  1,  1, -3,    -1,  0,  1,  1, &
      1, -2,  1, -3,     0, -1,  1,  1, &
     -1, -1,  2, -3,     0, -1,  1,  1, &
      2, -1, -1, -3,     1, -1,  0,  1, &
      1, -2,  1, -3,     1, -1,  0,  1, &
      1,  1, -2, -3,     1,  0, -1,  1, &
      2, -1, -1, -3,     1,  0, -1,  1, &
    ! pyramidal system: c+a slip <11.-3>{11.2} -- as for hexagonal Ice (Castelnau et al 1996, similar to twin system found below) 
      2, -1, -1, -3,     2, -1, -1,  2, & ! <11.-3>{11.2} shear = 2((c/a)^2-2)/(3 c/a)
      1,  1, -2, -3,     1,  1, -2,  2, & ! not sorted, just copied from twin system
     -1,  2, -1, -3,    -1,  2, -1,  2, &
     -2,  1,  1, -3,    -2,  1,  1,  2, &
     -1, -1,  2, -3,    -1, -1,  2,  2, &
      1, -2,  1, -3,     1, -2,  1,  2  &
     ],pReal),[ 4_pInt + 4_pInt,lattice_hex_Nslip])

 real(pReal), dimension(4+4,lattice_hex_Ntwin), parameter, private :: &
   lattice_hex_systemTwin =  reshape(real([&
      0,  1, -1,  1,     0, -1,  1,  2, & ! <-10.1>{10.2} shear = (3-(c/a)^2)/(sqrt(3) c/a)
     -1,  1,  0,  1,     1, -1,  0,  2, &
     -1,  0,  1,  1,     1,  0, -1,  2, & !!
      0, -1,  1,  1,     0,  1, -1,  2, &
      1, -1,  0,  1,    -1,  1,  0,  2, &
      1,  0, -1,  1,    -1,  0,  1,  2, &
      2, -1, -1, -3,     2, -1, -1,  2, & ! <11.-3>{11.2} shear = 2((c/a)^2-2)/(3 c/a)
      1,  1, -2, -3,     1,  1, -2,  2, & !!
     -1,  2, -1, -3,    -1,  2, -1,  2, &
     -2,  1,  1, -3,    -2,  1,  1,  2, &
     -1, -1,  2, -3,    -1, -1,  2,  2, &
      1, -2,  1, -3,     1, -2,  1,  2, &
     -2,  1,  1,  6,     2, -1, -1,  1, & ! <-1-1.6>{11.1} shear = 1/(c/a)
     -1, -1,  2,  6,     1,  1, -2,  1, & !!
      1, -2,  1,  6,    -1,  2, -1,  1, &
      2, -1, -1,  6,    -2,  1,  1,  1, &
      1,  1, -2,  6,    -1, -1,  2,  1, &
     -1,  2, -1,  6,     1, -2,  1,  1, &
      1,  0, -1, -2,     1,  0, -1,  1, & !! <10.-2>{10.1} shear = (4(c/a)^2-9)/(4 sqrt(3) c/a)
     -1,  0,  1, -2,    -1,  0,  1,  1, &
      0,  1, -1, -2,     0,  1, -1,  1, &
      0, -1,  1, -2,     0, -1,  1,  1, &
      1, -1,  0, -2,     1, -1,  0,  1, &
     -1,  1,  0, -2,    -1,  1,  0,  1  &
     ],pReal),[ 4_pInt + 4_pInt ,lattice_hex_Ntwin])                 !* Sort? Numbering of twin system follows Prof. Tom Bieler's scheme (to be consistent with his work); but numbering in data was restarted from 1 &

 integer(pInt), dimension(lattice_hex_Ntwin), parameter, private :: &
   lattice_hex_shearTwin = reshape(int( [&   ! indicator to formula further below
     1, &  ! {10.2}<-10.1>
     1, &
     1, &
     1, &
     1, &
     1, &
     2, &  ! {11.2}<11.-3>
     2, &
     2, &
     2, &
     2, &
     2, &
     3, &  ! {11.1}<-1-1.6>
     3, &
     3, &
     3, &
     3, &
     3, &
     4, &  ! {10.1}<10.-2>
     4, &
     4, &
     4, &
     4, &
     4  &
     ],pInt),[lattice_hex_Ntwin])

!* four different interaction type matrix
 !* 1. slip-slip interaction - 30 types
 !* 2. slip-twin interaction - 20 types
 !* 3. twin-twin interaction - 20 types
 !* 4. twin-slip interaction - 16 types
   
 integer(pInt), dimension(lattice_hex_Nslip,lattice_hex_Nslip), target, private :: &
   lattice_hex_interactionSlipSlip = reshape(int( [&
      1, 6, 6,  11,11,11,  15,15,15,15,15,15,  18,18,18,18,18,18,18,18,18,18,18,18,  20,20,20,20,20,20,  &  ! ---> slip
      6, 1, 6,  11,11,11,  15,15,15,15,15,15,  18,18,18,18,18,18,18,18,18,18,18,18,  20,20,20,20,20,20,  &  ! |
      6, 6, 1,  11,11,11,  15,15,15,15,15,15,  18,18,18,18,18,18,18,18,18,18,18,18,  20,20,20,20,20,20,  &  ! |
    !                                                                                                         v slip
     21,21,21,   2, 7, 7,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16,  19,19,19,19,19,19,  &
     21,21,21,   7, 2, 7,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16,  19,19,19,19,19,19,  &
     21,21,21,   7, 7, 2,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16,  19,19,19,19,19,19,  &
    !
     25,25,25,  22,22,22,   3, 8, 8, 8, 8, 8,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17,  &
     25,25,25,  22,22,22,   8, 3, 8, 8, 8, 8,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17,  &
     25,25,25,  22,22,22,   8, 8, 3, 8, 8, 8,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17,  &
     25,25,25,  22,22,22,   8, 8, 8, 3, 8, 8,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17,  &
     25,25,25,  22,22,22,   8, 8, 8, 8, 3, 8,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17,  &
     25,25,25,  22,22,22,   8, 8, 8, 8, 8, 3,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17,  &
    !
     28,28,28,  26,26,26,  23,23,23,23,23,23,   4, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,  14,14,14,14,14,14,  &
     28,28,28,  26,26,26,  23,23,23,23,23,23,   9, 4, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,  14,14,14,14,14,14,  &
     28,28,28,  26,26,26,  23,23,23,23,23,23,   9, 9, 4, 9, 9, 9, 9, 9, 9, 9, 9, 9,  14,14,14,14,14,14,  &
     28,28,28,  26,26,26,  23,23,23,23,23,23,   9, 9, 9, 4, 9, 9, 9, 9, 9, 9, 9, 9,  14,14,14,14,14,14,  &
     28,28,28,  26,26,26,  23,23,23,23,23,23,   9, 9, 9, 9, 4, 9, 9, 9, 9, 9, 9, 9,  14,14,14,14,14,14,  &
     28,28,28,  26,26,26,  23,23,23,23,23,23,   9, 9, 9, 9, 9, 4, 9, 9, 9, 9, 9, 9,  14,14,14,14,14,14,  &
     28,28,28,  26,26,26,  23,23,23,23,23,23,   9, 9, 9, 9, 9, 9, 4, 9, 9, 9, 9, 9,  14,14,14,14,14,14,  &
     28,28,28,  26,26,26,  23,23,23,23,23,23,   9, 9, 9, 9, 9, 9, 9, 4, 9, 9, 9, 9,  14,14,14,14,14,14,  &
     28,28,28,  26,26,26,  23,23,23,23,23,23,   9, 9, 9, 9, 9, 9, 9, 9, 4, 9, 9, 9,  14,14,14,14,14,14,  &
     28,28,28,  26,26,26,  23,23,23,23,23,23,   9, 9, 9, 9, 9, 9, 9, 9, 9, 4, 9, 9,  14,14,14,14,14,14,  &
     28,28,28,  26,26,26,  23,23,23,23,23,23,   9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 4, 9,  14,14,14,14,14,14,  &
     28,28,28,  26,26,26,  23,23,23,23,23,23,   9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 4,  14,14,14,14,14,14,  &
    !
     30,30,30,  29,29,29,  27,27,27,27,27,27,  24,24,24,24,24,24,24,24,24,24,24,24,   5,10,10,10,10,10,  &
     30,30,30,  29,29,29,  27,27,27,27,27,27,  24,24,24,24,24,24,24,24,24,24,24,24,  10, 5,10,10,10,10,  &
     30,30,30,  29,29,29,  27,27,27,27,27,27,  24,24,24,24,24,24,24,24,24,24,24,24,  10,10, 5,10,10,10,  &
     30,30,30,  29,29,29,  27,27,27,27,27,27,  24,24,24,24,24,24,24,24,24,24,24,24,  10,10,10, 5,10,10,  &
     30,30,30,  29,29,29,  27,27,27,27,27,27,  24,24,24,24,24,24,24,24,24,24,24,24,  10,10,10,10, 5,10,  &
     30,30,30,  29,29,29,  27,27,27,27,27,27,  24,24,24,24,24,24,24,24,24,24,24,24,  10,10,10,10,10, 5   &
     ],pInt),[lattice_hex_Nslip,lattice_hex_Nslip],order=[2,1])
  
!* isotropic interaction at the moment
 integer(pInt), dimension(lattice_hex_Nslip,lattice_hex_Ntwin), target, private :: &
   lattice_hex_interactionSlipTwin = reshape(int( [&
      1, 1, 1, 1, 1, 1,   2, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   4, 4, 4, 4, 4, 4, & ! --> twin
      1, 1, 1, 1, 1, 1,   2, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   4, 4, 4, 4, 4, 4, & ! |
      1, 1, 1, 1, 1, 1,   2, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   4, 4, 4, 4, 4, 4, & ! |
    !                                                                                   v
      5, 5, 5, 5, 5, 5,   6, 6, 6, 6, 6, 6,   7, 7, 7, 7, 7, 7,   8, 8, 8, 8, 8, 8, & ! slip
      5, 5, 5, 5, 5, 5,   6, 6, 6, 6, 6, 6,   7, 7, 7, 7, 7, 7,   8, 8, 8, 8, 8, 8, & 
      5, 5, 5, 5, 5, 5,   6, 6, 6, 6, 6, 6,   7, 7, 7, 7, 7, 7,   8, 8, 8, 8, 8, 8, &
    !
      9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, &
      9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, &
      9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, &
      9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, &
      9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, &
      9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, &
    !
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
    !
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20  &
     ],pInt),[lattice_hex_Nslip,lattice_hex_Ntwin],order=[2,1])

 !* isotropic interaction at the moment
 integer(pInt), dimension(lattice_hex_Ntwin,lattice_hex_Nslip), target, private :: &
   lattice_hex_interactionTwinSlip = reshape(int( [&
      1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17, & ! --> slip
      1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17, & ! |
      1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17, & ! |
      1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17, & ! v
      1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17, & ! twin
      1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17, &
    !
      2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14,  18,18,18,18,18,18, &
      2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14,  18,18,18,18,18,18, &
      2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14,  18,18,18,18,18,18, &
      2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14,  18,18,18,18,18,18, &
      2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14,  18,18,18,18,18,18, &
      2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14,  18,18,18,18,18,18, &
    !
      3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15,  19,19,19,19,19,19, &
      3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15,  19,19,19,19,19,19, &
      3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15,  19,19,19,19,19,19, &
      3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15,  19,19,19,19,19,19, &
      3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15,  19,19,19,19,19,19, &
      3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15,  19,19,19,19,19,19, &
    !
      4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16,  20,20,20,20,20,20, &
      4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16,  20,20,20,20,20,20, &
      4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16,  20,20,20,20,20,20, &
      4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16,  20,20,20,20,20,20, &
      4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16,  20,20,20,20,20,20, &
      4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16,  20,20,20,20,20,20  &
     ],pInt),[lattice_hex_Ntwin,lattice_hex_Nslip],order=[2,1])


 integer(pInt), dimension(lattice_hex_Ntwin,lattice_hex_Ntwin), target, private :: &
   lattice_hex_interactionTwinTwin = reshape(int( [&
      1, 5, 5, 5, 5, 5,   9, 9, 9, 9, 9, 9,  12,12,12,12,12,12,  14,14,14,14,14,14, &  ! ---> twin
      5, 1, 5, 5, 5, 5,   9, 9, 9, 9, 9, 9,  12,12,12,12,12,12,  14,14,14,14,14,14, &  ! |
      5, 5, 1, 5, 5, 5,   9, 9, 9, 9, 9, 9,  12,12,12,12,12,12,  14,14,14,14,14,14, &  ! |
      5, 5, 5, 1, 5, 5,   9, 9, 9, 9, 9, 9,  12,12,12,12,12,12,  14,14,14,14,14,14, &  ! v twin
      5, 5, 5, 5, 1, 5,   9, 9, 9, 9, 9, 9,  12,12,12,12,12,12,  14,14,14,14,14,14, &
      5, 5, 5, 5, 5, 1,   9, 9, 9, 9, 9, 9,  12,12,12,12,12,12,  14,14,14,14,14,14, &
    !
     15,15,15,15,15,15,   2, 6, 6, 6, 6, 6,  10,10,10,10,10,10,  13,13,13,13,13,13, &
     15,15,15,15,15,15,   6, 2, 6, 6, 6, 6,  10,10,10,10,10,10,  13,13,13,13,13,13, &
     15,15,15,15,15,15,   6, 6, 2, 6, 6, 6,  10,10,10,10,10,10,  13,13,13,13,13,13, &
     15,15,15,15,15,15,   6, 6, 6, 2, 6, 6,  10,10,10,10,10,10,  13,13,13,13,13,13, &
     15,15,15,15,15,15,   6, 6, 6, 6, 2, 6,  10,10,10,10,10,10,  13,13,13,13,13,13, &
     15,15,15,15,15,15,   6, 6, 6, 6, 6, 2,  10,10,10,10,10,10,  13,13,13,13,13,13, &
    !
     18,18,18,18,18,18,  16,16,16,16,16,16,   3, 7, 7, 7, 7, 7,  11,11,11,11,11,11, &
     18,18,18,18,18,18,  16,16,16,16,16,16,   7, 3, 7, 7, 7, 7,  11,11,11,11,11,11, &
     18,18,18,18,18,18,  16,16,16,16,16,16,   7, 7, 3, 7, 7, 7,  11,11,11,11,11,11, &
     18,18,18,18,18,18,  16,16,16,16,16,16,   7, 7, 7, 3, 7, 7,  11,11,11,11,11,11, &
     18,18,18,18,18,18,  16,16,16,16,16,16,   7, 7, 7, 7, 3, 7,  11,11,11,11,11,11, &
     18,18,18,18,18,18,  16,16,16,16,16,16,   7, 7, 7, 7, 7, 3,  11,11,11,11,11,11, &
    !
     20,20,20,20,20,20,  19,19,19,19,19,19,  17,17,17,17,17,17,   4, 8, 8, 8, 8, 8, &
     20,20,20,20,20,20,  19,19,19,19,19,19,  17,17,17,17,17,17,   8, 4, 8, 8, 8, 8, &
     20,20,20,20,20,20,  19,19,19,19,19,19,  17,17,17,17,17,17,   8, 8, 4, 8, 8, 8, &
     20,20,20,20,20,20,  19,19,19,19,19,19,  17,17,17,17,17,17,   8, 8, 8, 4, 8, 8, &
     20,20,20,20,20,20,  19,19,19,19,19,19,  17,17,17,17,17,17,   8, 8, 8, 8, 4, 8, &
     20,20,20,20,20,20,  19,19,19,19,19,19,  17,17,17,17,17,17,   8, 8, 8, 8, 8, 4  &
     ],pInt),[lattice_hex_Ntwin,lattice_hex_Ntwin],order=[2,1])

 public :: &
  lattice_init, &
  lattice_initializeStructure, &
  lattice_symmetryType

contains


!--------------------------------------------------------------------------------------------------
!> @brief Maps structure to symmetry type 
!> @details fcc(1) and bcc(2) are cubic(1) hex(3+) is hexagonal(2)
!--------------------------------------------------------------------------------------------------
integer(pInt) pure function lattice_symmetryType(structID)

 implicit none
 integer(pInt), intent(in) :: structID

 select case(structID)
   case (1_pInt,2_pInt)
     lattice_symmetryType = 1_pInt
   case (3_pInt:)
     lattice_symmetryType = 2_pInt
   case default
     lattice_symmetryType = 0_pInt
  end select

 return
 
end function lattice_symmetryType


!--------------------------------------------------------------------------------------------------
!> @brief Module initialization
!--------------------------------------------------------------------------------------------------
subroutine lattice_init

 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use IO,       only: IO_open_file,&
                     IO_open_jobFile_stat, &
                     IO_countSections, &
                     IO_countTagInPart, &
                     IO_error
 use material, only: material_configfile, &
                     material_localFileExt, &
                     material_partPhase
 use debug,    only: debug_level, &
                     debug_lattice, &
                     debug_levelBasic

 implicit none
 integer(pInt), parameter :: fileunit = 200_pInt
 integer(pInt) :: Nsections

 !$OMP CRITICAL (write2out)
   write(6,*)
   write(6,*) '<<<+-  lattice init  -+>>>'
   write(6,*) '$Id$'
#include "compilation_info.f90"
 !$OMP END CRITICAL (write2out)

 if (.not. IO_open_jobFile_stat(fileunit,material_localFileExt)) then                               ! no local material configuration present...
   call IO_open_file(fileunit,material_configFile)                                                  ! ... open material.config file
 endif
 Nsections = IO_countSections(fileunit,material_partPhase)
 lattice_Nstructure = 2_pInt + sum(IO_countTagInPart(fileunit,material_partPhase,'covera_ratio',Nsections)) ! fcc + bcc + all hex
! lattice_Nstructure = Nsections + 2_pInt                                                           ! most conservative assumption
 close(fileunit)

 if (iand(debug_level(debug_lattice),debug_levelBasic) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
     write(6,'(a16,1x,i5)') '# phases:',Nsections
     write(6,'(a16,1x,i5)') '# structures:',lattice_Nstructure
     write(6,*)
   !$OMP END CRITICAL (write2out)
 endif

 allocate(lattice_Sslip(3,3,lattice_maxNslip,lattice_Nstructure)); lattice_Sslip   = 0.0_pReal
 allocate(lattice_Sslip_v(6,lattice_maxNslip,lattice_Nstructure)); lattice_Sslip_v = 0.0_pReal
 allocate(lattice_sd(3,lattice_maxNslip,lattice_Nstructure));      lattice_sd      = 0.0_pReal
 allocate(lattice_st(3,lattice_maxNslip,lattice_Nstructure));      lattice_st      = 0.0_pReal
 allocate(lattice_sn(3,lattice_maxNslip,lattice_Nstructure));      lattice_sn      = 0.0_pReal

 allocate(lattice_Qtwin(3,3,lattice_maxNtwin,lattice_Nstructure)); lattice_Qtwin   = 0.0_pReal
 allocate(lattice_Stwin(3,3,lattice_maxNtwin,lattice_Nstructure)); lattice_Stwin   = 0.0_pReal
 allocate(lattice_Stwin_v(6,lattice_maxNtwin,lattice_Nstructure)); lattice_Stwin_v = 0.0_pReal
 allocate(lattice_td(3,lattice_maxNtwin,lattice_Nstructure));      lattice_td      = 0.0_pReal
 allocate(lattice_tt(3,lattice_maxNtwin,lattice_Nstructure));      lattice_tt      = 0.0_pReal
 allocate(lattice_tn(3,lattice_maxNtwin,lattice_Nstructure));      lattice_tn      = 0.0_pReal

 allocate(lattice_shearTwin(lattice_maxNtwin,lattice_Nstructure)); lattice_shearTwin = 0.0_pReal

 allocate(lattice_NslipSystem(lattice_maxNslipFamily,lattice_Nstructure)); lattice_NslipSystem = 0_pInt
 allocate(lattice_NtwinSystem(lattice_maxNtwinFamily,lattice_Nstructure)); lattice_NtwinSystem = 0_pInt

 allocate(lattice_interactionSlipSlip(lattice_maxNslip,lattice_maxNslip,lattice_Nstructure))
          lattice_interactionSlipSlip = 0_pInt ! other:me
 allocate(lattice_interactionSlipTwin(lattice_maxNslip,lattice_maxNtwin,lattice_Nstructure))
          lattice_interactionSlipTwin = 0_pInt ! other:me
 allocate(lattice_interactionTwinSlip(lattice_maxNtwin,lattice_maxNslip,lattice_Nstructure))
          lattice_interactionTwinSlip = 0_pInt ! other:me
 allocate(lattice_interactionTwinTwin(lattice_maxNtwin,lattice_maxNtwin,lattice_Nstructure))
          lattice_interactionTwinTwin = 0_pInt ! other:me

end subroutine lattice_init


!--------------------------------------------------------------------------------------------------
!> @brief   Calculation of Schmid matrices, etc.
!--------------------------------------------------------------------------------------------------
integer(pInt) function lattice_initializeStructure(struct,CoverA)
 use prec, only: pReal,pInt
 use math, only: math_vectorproduct, &
                 math_tensorproduct, &
                 math_norm3, &
                 math_trace33, &
                 math_symmetric33, &
                 math_Mandel33to6, &
                 math_axisAngleToR, &
                 INRAD
 use IO, only: IO_error
 
 implicit none
 character(len=*) struct
 real(pReal) CoverA
 real(pReal), dimension(3,lattice_maxNslip) :: sd = 0.0_pReal, &
                                               sn = 0.0_pReal
 real(pReal), dimension(3,lattice_maxNtwin) :: td = 0.0_pReal, &
                                               tn = 0.0_pReal
 real(pReal), dimension(lattice_maxNtwin) ::   ts = 0.0_pReal
 integer(pInt), dimension(lattice_maxNslipFamily) :: myNslipSystem = 0_pInt
 integer(pInt), dimension(lattice_maxNtwinFamily) :: myNtwinSystem = 0_pInt
 integer(pInt) :: i,myNslip,myNtwin,myStructure = 0_pInt
 logical :: processMe

 processMe = .false.

 select case(struct(1:3))                          ! check first three chars of structure name
   case ('fcc')
     myStructure = 1_pInt
     myNslipSystem = lattice_fcc_NslipSystem       ! size of slip system families
     myNtwinSystem = lattice_fcc_NtwinSystem       ! size of twin system families
     myNslip = lattice_fcc_Nslip                   ! overall number of slip systems
     myNtwin = lattice_fcc_Ntwin                   ! overall number of twin systems
     lattice_fcc_Nstructure = lattice_fcc_Nstructure + 1_pInt    ! count fcc instances
     if (lattice_fcc_Nstructure == 1_pInt) then    ! me is first fcc structure
       processMe = .true.
       do i = 1_pInt,myNslip                            ! assign slip system vectors
         sd(1:3,i) = lattice_fcc_systemSlip(1:3,i)
         sn(1:3,i) = lattice_fcc_systemSlip(4:6,i)
       enddo
       do i = 1_pInt,myNtwin                            ! assign twin system vectors and shears
         td(1:3,i) = lattice_fcc_systemTwin(1:3,i)
         tn(1:3,i) = lattice_fcc_systemTwin(4:6,i)
         ts(i)     = lattice_fcc_shearTwin(i)
       enddo
       interactionSlipSlip => lattice_fcc_interactionSlipSlip
       interactionSlipTwin => lattice_fcc_interactionSlipTwin
       interactionTwinSlip => lattice_fcc_interactionTwinSlip
       interactionTwinTwin => lattice_fcc_interactionTwinTwin
     endif
     
   case ('bcc')
     myStructure = 2_pInt
     myNslipSystem = lattice_bcc_NslipSystem       ! size of slip system families
     myNtwinSystem = lattice_bcc_NtwinSystem       ! size of twin system families
     myNslip = lattice_bcc_Nslip                   ! overall number of slip systems
     myNtwin = lattice_bcc_Ntwin                   ! overall number of twin systems
     lattice_bcc_Nstructure = lattice_bcc_Nstructure + 1_pInt    ! count bcc instances
     if (lattice_bcc_Nstructure == 1_pInt) then    ! me is first bcc structure
       processMe = .true.
       do i = 1_pInt,myNslip                            ! assign slip system vectors
         sd(1:3,i) = lattice_bcc_systemSlip(1:3,i)
         sn(1:3,i) = lattice_bcc_systemSlip(4:6,i)
       enddo
       do i = 1_pInt,myNtwin                            ! assign twin system vectors and shears
         td(1:3,i) = lattice_bcc_systemTwin(1:3,i)
         tn(1:3,i) = lattice_bcc_systemTwin(4:6,i)
         ts(i)     = lattice_bcc_shearTwin(i)
       enddo
       interactionSlipSlip => lattice_bcc_interactionSlipSlip
       interactionSlipTwin => lattice_bcc_interactionSlipTwin
       interactionTwinSlip => lattice_bcc_interactionTwinSlip
       interactionTwinTwin => lattice_bcc_interactionTwinTwin
     endif
     
   case ('hex')
     if (CoverA >= 1.0_pReal) then                 ! checking physical significance of c/a
       lattice_hex_Nstructure = lattice_hex_Nstructure + 1_pInt  ! count instances of hex structures
       myStructure = 2_pInt + lattice_hex_Nstructure             ! 3,4,5,.. for hex
       myNslipSystem = lattice_hex_NslipSystem     ! size of slip system families
       myNtwinSystem = lattice_hex_NtwinSystem     ! size of twin system families
       myNslip = lattice_hex_Nslip                 ! overall number of slip systems
       myNtwin = lattice_hex_Ntwin                 ! overall number of twin systems
       processMe = .true.
! converting from 4 axes coordinate system (a1=a2=a3=c) to ortho-hexgonal system (a, b, c)
       do i = 1_pInt,myNslip
         sd(1,i) =  lattice_hex_systemSlip(1,i)*1.5_pReal ! direction [uvtw]->[3u/2 (u+2v)*sqrt(3)/2 w*(c/a)]
         sd(2,i) = (lattice_hex_systemSlip(1,i)+2.0_pReal*lattice_hex_systemSlip(2,i))*(0.5_pReal*sqrt(3.0_pReal))
         sd(3,i) =  lattice_hex_systemSlip(4,i)*CoverA
         sn(1,i) =  lattice_hex_systemSlip(5,i)           ! plane (hkil)->(h (h+2k)/sqrt(3) l/(c/a))
         sn(2,i) = (lattice_hex_systemSlip(5,i)+2.0_pReal*lattice_hex_systemSlip(6,i))/sqrt(3.0_pReal)
         sn(3,i) =  lattice_hex_systemSlip(8,i)/CoverA
       enddo
       do i = 1_pInt,myNtwin
         td(1,i) =  lattice_hex_systemTwin(1,i)*1.5_pReal
         td(2,i) = (lattice_hex_systemTwin(1,i)+2.0_pReal*lattice_hex_systemTwin(2,i))*(0.5_pReal*sqrt(3.0_pReal))
         td(3,i) =  lattice_hex_systemTwin(4,i)*CoverA
         tn(1,i) =  lattice_hex_systemTwin(5,i)
         tn(2,i) = (lattice_hex_systemTwin(5,i)+2.0_pReal*lattice_hex_systemTwin(6,i))/sqrt(3.0_pReal)
         tn(3,i) =  lattice_hex_systemTwin(8,i)/CoverA

         select case(lattice_hex_shearTwin(i))                                          ! from Christian & Mahajan 1995 p.29
           case (1_pInt)                                                                ! {10.2}<-10.1>
                    ts(i) = (3.0_pReal-CoverA*CoverA)/sqrt(3.0_pReal)/CoverA
           case (2_pInt)                                                                ! {11.2}<11.-3>
                    ts(i) = 2.0_pReal*(CoverA*CoverA-2.0_pReal)/3.0_pReal/CoverA
           case (3_pInt)                                                                ! {11.1}<-1-1.6>
                    ts(i) = 1.0_pReal/CoverA
           case (4_pInt)                                                                ! {10.1}<10.-2>
                    ts(i) = (4.0_pReal*CoverA*CoverA-9.0_pReal)/4.0_pReal/sqrt(3.0_pReal)/CoverA
         end select

       enddo
       interactionSlipSlip => lattice_hex_interactionSlipSlip
       interactionSlipTwin => lattice_hex_interactionSlipTwin
       interactionTwinSlip => lattice_hex_interactionTwinSlip
       interactionTwinTwin => lattice_hex_interactionTwinTwin
     endif
 end select

 if (processMe) then
   if  (myStructure > lattice_Nstructure) &
     call IO_error(666_pInt,myStructure,ext_msg = 'structure index out of bounds')     ! check for memory leakage
   do i = 1_pInt,myNslip                                                               ! store slip system vectors and Schmid matrix for my structure
     lattice_sd(1:3,i,myStructure) = sd(1:3,i)/math_norm3(sd(1:3,i))                   ! make unit vector
     lattice_sn(1:3,i,myStructure) = sn(1:3,i)/math_norm3(sn(1:3,i))                   ! make unit vector
     lattice_st(1:3,i,myStructure) = math_vectorproduct(lattice_sd(1:3,i,myStructure), &
                                                        lattice_sn(1:3,i,myStructure))
     lattice_Sslip(1:3,1:3,i,myStructure) = math_tensorproduct(lattice_sd(1:3,i,myStructure), &
                                                               lattice_sn(1:3,i,myStructure))
     lattice_Sslip_v(1:6,i,myStructure) = math_Mandel33to6(math_symmetric33(lattice_Sslip(1:3,1:3,i,myStructure)))
     if (abs(math_trace33(lattice_Sslip(1:3,1:3,i,myStructure))) > 1.0e-8) &
       call IO_error(0_pInt,myStructure,i,0_pInt,ext_msg = 'dilatational slip Schmid matrix')
   enddo
   do i = 1_pInt,myNtwin                                              ! store twin system vectors and Schmid plus rotation matrix for my structure
     lattice_td(1:3,i,myStructure) = td(1:3,i)/math_norm3(sd(1:3,i))                   ! make unit vector
     lattice_tn(1:3,i,myStructure) = tn(1:3,i)/math_norm3(sn(1:3,i))                   ! make unit vector
     lattice_tt(1:3,i,myStructure) = math_vectorproduct(lattice_td(1:3,i,myStructure), &
                                                        lattice_tn(1:3,i,myStructure))
     lattice_Stwin(1:3,1:3,i,myStructure) = math_tensorproduct(lattice_td(1:3,i,myStructure), &
                                                               lattice_tn(1:3,i,myStructure))
     lattice_Stwin_v(1:6,i,myStructure) = math_Mandel33to6(math_symmetric33(lattice_Stwin(1:3,1:3,i,myStructure)))
     lattice_Qtwin(1:3,1:3,i,myStructure) = math_AxisAngleToR(tn(1:3,i),180.0_pReal*INRAD)
     lattice_shearTwin(i,myStructure) = ts(i)
     if (abs(math_trace33(lattice_Stwin(1:3,1:3,i,myStructure))) > 1.0e-8) &
       call IO_error(0_pInt,myStructure,i,0_pInt,ext_msg = 'dilatational twin Schmid matrix')
   enddo
   lattice_NslipSystem(1:lattice_maxNslipFamily,myStructure) = myNslipSystem                                ! number of slip systems in each family
   lattice_NtwinSystem(1:lattice_maxNtwinFamily,myStructure) = myNtwinSystem                                ! number of twin systems in each family
   lattice_interactionSlipSlip(1:myNslip,1:myNslip,myStructure) = interactionSlipSlip(1:myNslip,1:myNslip)
   lattice_interactionSlipTwin(1:myNslip,1:myNtwin,myStructure) = interactionSlipTwin(1:myNslip,1:myNtwin)
   lattice_interactionTwinSlip(1:myNtwin,1:myNslip,myStructure) = interactionTwinSlip(1:myNtwin,1:myNslip)
   lattice_interactionTwinTwin(1:myNtwin,1:myNtwin,myStructure) = interactionTwinTwin(1:myNtwin,1:myNtwin)
 endif

 lattice_initializeStructure = myStructure        ! report my structure index back

end function lattice_initializeStructure

end module lattice
