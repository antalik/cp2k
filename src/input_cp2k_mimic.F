!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Module for MiMiC input
!> \par History
!>      03.2023 created [AA]
!> \author Andrej Antalik
! **************************************************************************************************
MODULE input_cp2k_mimic
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_val_types,                 ONLY: char_t,&
                                              integer_t,&
                                              lchar_t,&
                                              logical_t,&
                                              real_t

#include "./base/base_uses.f90"


   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_mimic'

   PUBLIC :: create_mimic_section

CONTAINS

! **************************************************************************************************
!> \brief Creates the MiMiC section of the input
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_mimic_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="MIMIC", &
                          description="Input for QM/MM calculations with the MiMiC framework.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.) !, &
!                          citations=(/Laino2005, Laino2006/))

      NULLIFY(keyword, subsection)


      CALL keyword_create(keyword, __LOCATION__, name="MULTIPOLE_ORDER", &
                          description="TBD", &
                          usage="MULTIPOLE_ORDER {Integer}", default_i_val=5, &
                          n_var=1, type_of_var=integer_t, repeats=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)





   END SUBROUTINE create_mimic_section

END MODULE input_cp2k_mimic
