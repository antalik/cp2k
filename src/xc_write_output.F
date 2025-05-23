!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Writes information on XC functionals to output
! **************************************************************************************************
MODULE xc_write_output

   USE input_constants,                 ONLY: xc_none
   USE input_cp2k_check,                ONLY: xc_functionals_expand
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_get_subs_vals2,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length
   USE xc_derivatives,                  ONLY: xc_functional_get_info
   USE xc_libxc,                        ONLY: libxc_check_existence_in_libxc,&
                                              libxc_get_reference_length
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_write_output'

   PUBLIC :: xc_write

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param iounit ...
!> \param xc_section ...
!> \param lsd ...
! **************************************************************************************************
   SUBROUTINE xc_write(iounit, xc_section, lsd)
      INTEGER, INTENT(IN)                                :: iounit
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(IN)                                :: lsd

      CHARACTER(LEN=2*default_string_length)             :: shortform
      CHARACTER(LEN=:), ALLOCATABLE                      :: reference
      INTEGER                                            :: ifun, il, myfun
      TYPE(section_vals_type), POINTER                   :: xc_fun, xc_fun_section

      IF (iounit > 0) THEN

         xc_fun_section => section_vals_get_subs_vals(xc_section, &
                                                      "XC_FUNCTIONAL")
         CALL section_vals_val_get(xc_fun_section, "_SECTION_PARAMETERS_", i_val=myfun)
         IF (myfun /= xc_none) THEN

            !check if FUNCTIONAL_ROUTINE keyword present

            CALL xc_functionals_expand(xc_fun_section, xc_section)
            ifun = 0
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (libxc_check_existence_in_libxc(xc_fun)) THEN
                  ALLOCATE (CHARACTER(LEN=libxc_get_reference_length(xc_fun, lsd)) :: reference)
               ELSE
                  ALLOCATE (CHARACTER(LEN=20*default_string_length) :: reference)
               END IF
               CALL xc_functional_get_info(xc_fun, lsd=lsd, reference=reference, shortform=shortform)
               WRITE (iounit, fmt="(' FUNCTIONAL| ',a,':')") &
                  TRIM(xc_fun%section%name)
               DO il = 1, LEN_TRIM(reference), 67
                  WRITE (iounit, fmt="(' FUNCTIONAL| ',a67)") reference(il:)
               END DO
               DEALLOCATE (reference)
            END DO
         ELSE
            WRITE (iounit, fmt="(' FUNCTIONAL| NO EXCHANGE-CORRELATION FUNCTIONAL USED.')")
         END IF
      END IF

   END SUBROUTINE xc_write

END MODULE xc_write_output
