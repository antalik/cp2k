!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Module containing a few utility routines for a MiMiC run 
!> \par History
!>      03.2023 Created [AA]
!> \author Andrej Antalik
! **************************************************************************************************

MODULE mimic_utils

   USE force_env_types,                 ONLY: force_env_get, &
                                              force_env_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals, &
                                              section_vals_type, &
                                              section_vals_val_get
   USE kinds,                           ONLY: dp

#include "../base/base_uses.f90"

   IMPLICIT NONE
   
   PRIVATE

   PUBLIC :: mimic_get_input_val

CONTAINS

! **************************************************************************************************
!> \brief Fetches a value correponding to a keyword in the MIMIC input section
!> \param force_env ...
!> \param keyword ...
!> \param var [output] return value
! **************************************************************************************************
   SUBROUTINE mimic_get_input_val(force_env, keyword, var)
      TYPE(force_env_type), POINTER                      :: force_env
      CHARACTER(LEN=*)                                   :: keyword
      INTEGER, INTENT(OUT)                               :: var
      
      TYPE(section_vals_type), POINTER                   :: input, mimic_input

      CALL force_env_get(force_env, input=input)
      mimic_input => section_vals_get_subs_vals(input, "MIMIC")
      CALL section_vals_val_get(mimic_input, keyword, i_val=var)

   END SUBROUTINE mimic_get_input_val

END MODULE mimic_utils
