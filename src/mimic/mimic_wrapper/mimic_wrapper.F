!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wrapper module for MiMiC and MiMiC Communication Library (MCL) routines
!> \par History
!>       03.2023 Created [AA]
!> \author Andrej Antalik
! **************************************************************************************************

MODULE mimic_wrapper

   USE mimic_send_recv
#if defined(__MIMIC)
   USE mcl_fortran,                     ONLY: mcl_destroy, &
                                              mcl_handshake, &
                                              mcl_init
#endif

#include "../../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: mimic_init, mimic_handshake, mimic_destroy, mimic_receive, mimic_send

   INTERFACE mimic_receive
      MODULE PROCEDURE recv_int, recv_int_array1d, &
                       recv_dp,  recv_dp_array1d 
   END INTERFACE

   INTERFACE mimic_send
      MODULE PROCEDURE send_int, send_int_array1d, &
                       send_dp,  send_dp_array1d 
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Initialize MiMiC Communicator
!> \param communicator in: world communicator out: split communicator
!> \note Should be called right after MPI_init()
! **************************************************************************************************
   SUBROUTINE mimic_init(communicator)
      INTEGER, INTENT(INOUT)                           :: communicator

#if defined(__MIMIC)
      CALL mcl_init(communicator)
#else
      CPABORT("This CP2K executable has not been linked against MiMiC.")
#endif

   END SUBROUTINE mimic_init


! **************************************************************************************************
!> \brief Client server handshake function
!> \param path path to a working folder
! **************************************************************************************************
   SUBROUTINE mimic_handshake(path)
      CHARACTER(LEN=*), INTENT(IN)                     :: path

#if defined(__MIMIC)
      CALL mcl_handshake(path, ";", .false.)
#else
      CPABORT("This CP2K executable has not been linked against MiMiC.")
#endif

   END SUBROUTINE mimic_handshake


! **************************************************************************************************
!> \brief Destroys the endpoint
! **************************************************************************************************
   SUBROUTINE mimic_destroy()

#if defined(__MIMIC)
      CALL mcl_destroy()
#else
      CPABORT("This CP2K executable has not been linked against MiMiC.")
#endif

   END SUBROUTINE mimic_destroy

END MODULE mimic_wrapper

