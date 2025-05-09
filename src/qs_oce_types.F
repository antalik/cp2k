!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
MODULE qs_oce_types

   USE block_p_types,                   ONLY: block_p_type
   USE kinds,                           ONLY: dp
   USE sap_kind_types,                  ONLY: release_sap_int,&
                                              sap_int_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters (only in this module)

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_oce_types'

! *** Define a oce matrix type ***

! **************************************************************************************************
   TYPE qlist_type
      TYPE(block_p_type), DIMENSION(:), POINTER :: sac => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER      :: r2 => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER   :: r => NULL()
      INTEGER                                   :: n = -1
      INTEGER, DIMENSION(:), POINTER            :: index_list => NULL()
      INTEGER, DIMENSION(:), POINTER            :: list => NULL()
   END TYPE qlist_type

! **************************************************************************************************
   TYPE oce_matrix_type
      TYPE(sap_int_type), DIMENSION(:), POINTER        :: intac => NULL()
   END TYPE

! *** Public data types ***

   PUBLIC :: oce_matrix_type

! *** Public subroutines ***

   PUBLIC :: allocate_oce_set, &
             create_oce_set, &
             deallocate_oce_set

CONTAINS

! **************************************************************************************************
!> \brief   Allocate and initialize the matrix set of oce coefficients.
!> \param oce_set ...
!> \param nkind ...
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE allocate_oce_set(oce_set, nkind)
      TYPE(oce_matrix_type), POINTER                     :: oce_set
      INTEGER, INTENT(IN)                                :: nkind

      INTEGER                                            :: i

      ALLOCATE (oce_set%intac(nkind*nkind))
      DO i = 1, nkind*nkind
         NULLIFY (oce_set%intac(i)%alist)
         NULLIFY (oce_set%intac(i)%asort)
         NULLIFY (oce_set%intac(i)%aindex)
      END DO

   END SUBROUTINE allocate_oce_set

! **************************************************************************************************
!> \brief ...
!> \param oce_set ...
! **************************************************************************************************
   SUBROUTINE create_oce_set(oce_set)

      TYPE(oce_matrix_type), POINTER                     :: oce_set

      IF (ASSOCIATED(oce_set)) CALL deallocate_oce_set(oce_set)

      ALLOCATE (oce_set)

      NULLIFY (oce_set%intac)

   END SUBROUTINE create_oce_set

! **************************************************************************************************
!> \brief  Deallocate the matrix set of oce coefficients
!> \param oce_set ...
!> \date
!> \author
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE deallocate_oce_set(oce_set)
      TYPE(oce_matrix_type), POINTER                     :: oce_set

      IF (.NOT. ASSOCIATED(oce_set)) RETURN

      IF (ASSOCIATED(oce_set%intac)) CALL release_sap_int(oce_set%intac)

      DEALLOCATE (oce_set)

   END SUBROUTINE deallocate_oce_set

! **************************************************************************************************

END MODULE qs_oce_types
