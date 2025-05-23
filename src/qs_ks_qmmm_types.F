!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      05.2004 [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
MODULE qs_ks_qmmm_types
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cube_utils,                      ONLY: cube_info_type,&
                                              destroy_cube_info
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_release,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_qmmm_types'

   PUBLIC :: qs_ks_qmmm_env_type
   PUBLIC :: qs_ks_qmmm_release

! **************************************************************************************************
!> \brief calculation environment to calculate the ks_qmmm matrix,
!>      holds the QM/MM potential and all the needed variables to
!>      compute the QM/MM electrostatic 1-electron ks matrix
!>      assumes that the core hamiltonian and energy are up to date.
!>      v_metal_rspace is the potential at the metal sites within the image
!>      charge approach
!> \par History
!>      05.2004 created [tlaino]
!>      01.2012 added v_metal_rspace [dgolze]
!> \author Teodoro Laino
! **************************************************************************************************
   TYPE qs_ks_qmmm_env_type
      INTEGER :: n_evals = -1
      REAL(KIND=dp)                               :: pc_ener = -1.0_dp
      TYPE(pw_env_type), POINTER                  :: pw_env => NULL()
      TYPE(pw_r3d_rs_type)                             :: v_qmmm_rspace = pw_r3d_rs_type()
      TYPE(pw_r3d_rs_type), POINTER                             :: v_metal_rspace => NULL()
      TYPE(cube_info_type), DIMENSION(:), POINTER  :: cube_info => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), &
         POINTER                                :: matrix_h => NULL()
   END TYPE qs_ks_qmmm_env_type

CONTAINS

! **************************************************************************************************
!> \brief releases the ks_qmmm_env (see doc/ReferenceCounting.html)
!> \param ks_qmmm_env the ks_qmmm_env to be released
!> \par History
!>      05.2004 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE qs_ks_qmmm_release(ks_qmmm_env)
      TYPE(qs_ks_qmmm_env_type), INTENT(INOUT)           :: ks_qmmm_env

      INTEGER                                            :: i
      TYPE(pw_pool_type), POINTER                        :: pool

      CALL pw_env_get(ks_qmmm_env%pw_env, auxbas_pw_pool=pool)
      CALL pool%give_back_pw(ks_qmmm_env%v_qmmm_rspace)
      CALL pw_env_release(ks_qmmm_env%pw_env)
      IF (ASSOCIATED(ks_qmmm_env%cube_info)) THEN
         DO i = 1, SIZE(ks_qmmm_env%cube_info)
            CALL destroy_cube_info(ks_qmmm_env%cube_info(i))
         END DO
         DEALLOCATE (ks_qmmm_env%cube_info)
      END IF
      IF (ASSOCIATED(ks_qmmm_env%matrix_h)) THEN
         CALL dbcsr_deallocate_matrix_set(ks_qmmm_env%matrix_h)
      END IF

   END SUBROUTINE qs_ks_qmmm_release

END MODULE qs_ks_qmmm_types
