!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Module containing a MiMiC communicator class
!> \par History
!>      03.2023 Created [AA]
!> \author Andrej Antalik
! **************************************************************************************************

MODULE mimic_communicator

   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get, &
                                              cp_subsys_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE force_env_types,                 ONLY: force_env_get, &
                                              force_env_type
   USE kinds,                           ONLY: default_path_length, &
                                              dp
   USE machine,                         ONLY: m_getcwd
   USE message_passing,                 ONLY: mp_para_env_type
   USE mimic_utils,                     ONLY: mimic_get_input_val
   USE particle_list_types,             ONLY: particle_list_type
   USE qs_environment_types,            ONLY: qs_environment_type, &
                                              get_qs_env
   USE qs_energy_types,                 ONLY: qs_energy_type

   USE mimic_wrapper,                   ONLY: mimic_destroy, &
                                              mimic_handshake, &
                                              mimic_receive, &
                                              mimic_send
#if defined(__MIMIC)
   USE mcl_requests
#endif

#include "../base/base_uses.f90"

   IMPLICIT NONE
   
   PRIVATE

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = 'mimic_communicator'

! **************************************************************************************************
!> \brief MiMiC communicator class that facilitates MiMiC client-server data exchange 
!> \par History
!>      03.2023 Created [AA]
! **************************************************************************************************
   TYPE, PUBLIC :: mimic_communicator_type
      PRIVATE
         !> communication
         TYPE(mp_para_env_type), POINTER                 :: para_env => Null()
         LOGICAL                                         :: is_ionode
         INTEGER                                         :: command = -1
         INTEGER                                         :: MCL_DATA = -1, &
                                                            MCL_COMMAND = -1
         !> data
         TYPE(particle_list_type), POINTER               :: atoms => Null()
         TYPE(atomic_kind_list_type), POINTER            :: kinds => Null()
         TYPE(qs_energy_type), POINTER                   :: energy => Null()
         INTEGER                                         :: multipole_order = -1, &
                                                            n_atoms = -1, &
                                                            n_kinds = -1
      CONTAINS

         PROCEDURE          :: initialize
         PROCEDURE          :: handshake
         PROCEDURE          :: finalize
         PROCEDURE          :: receive_command
         PROCEDURE          :: execute_command
         PROCEDURE          :: is_cmd
         PROCEDURE, PRIVATE :: send_value
         PROCEDURE, PRIVATE :: send_atom_kinds
         PROCEDURE, PRIVATE :: send_multipoles
         PROCEDURE, PRIVATE :: send_kind_masses
         PROCEDURE, PRIVATE :: send_atom_ids
         PROCEDURE, PRIVATE :: send_kind_elements
         PROCEDURE, PRIVATE :: send_coords
         PROCEDURE, PRIVATE :: receive_coords
         PROCEDURE, PRIVATE :: send_forces

   END TYPE mimic_communicator_type

CONTAINS

! **************************************************************************************************
!> \brief Initializes MiMiC communicator loading input and by saving pointers to relevant data
!> \par History
!>      03.2023 Created [AA]
! **************************************************************************************************
   SUBROUTINE initialize(this, force_env)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      TYPE(force_env_type), POINTER                      :: force_env
      
      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':initialize'
      
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      NULLIFY (subsys, qs_env)

      CPASSERT(ASSOCIATED(force_env))

      CALL force_env_get(force_env, subsys=subsys, para_env=this%para_env)
      CALL cp_subsys_get(subsys, natom=this%n_atoms, particles=this%atoms, &
                                 nkind=this%n_kinds, atomic_kinds=this%kinds)
      CALL force_env_get(force_env, qs_env=qs_env)
      CALL get_qs_env(qs_env, energy=this%energy)
      
      this%is_ionode = this%para_env%is_source()

      CALL mimic_get_input_val(force_env, "MULTIPOLE_ORDER", this%multipole_order)
      !> placeholder value for now
      this%multipole_order = 0
     
#if defined(__MIMIC)
      this%MCL_DATA = MCL_DATA
      this%MCL_COMMAND = MCL_COMMAND
#endif

      CALL timestop(handle)

   END SUBROUTINE initialize


! **************************************************************************************************
!> \brief Performs a handshake with the MiMiC server
! **************************************************************************************************
   SUBROUTINE handshake(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':handshake'

      CHARACTER(LEN=default_path_length)                 :: path
      INTEGER                                            :: handle

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      CALL m_getcwd(path)
      CALL mimic_handshake(trim(path))

      CALL timestop(handle)

   END SUBROUTINE handshake


! **************************************************************************************************
!> \brief Destroys the endpoint
! **************************************************************************************************
   SUBROUTINE finalize(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':finalize'

      INTEGER                                            :: handle

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      CALL mimic_destroy()

      CALL timestop(handle)

   END SUBROUTINE finalize


! **************************************************************************************************
!> \brief Receives a MiMiC command
! **************************************************************************************************
   SUBROUTINE receive_command(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':receive_command'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      this%command = -1
      IF (this%is_ionode) THEN
         CALL mimic_receive(this%command, 1, this%MCL_COMMAND, 0)
      END IF
      CALL this%para_env%bcast(this%command)

      CALL timestop(handle)

   END SUBROUTINE receive_command


! **************************************************************************************************
!> \brief Checks whether 'command' corresponds to the command received from the server
!> \param command command to compare with the received command
!> \retval is_match true if 'command' matches the received command 'this%command'
! **************************************************************************************************
   LOGICAL FUNCTION is_cmd(this, command) RESULT(is_match)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      INTEGER, INTENT(IN)                                :: command

      is_match = (this%command == command)

   END FUNCTION is_cmd


! **************************************************************************************************
!> \brief Executes an action corresponding to one of the communication commands
!> \param executed [output] returns true if one of the listed actions has been executed
!> \par History
!>      03.2023 Created [AA]
! **************************************************************************************************
   SUBROUTINE execute_command(this, executed)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      LOGICAL, INTENT(INOUT)                             :: executed
          
      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':execute_command'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      executed = .TRUE.

#if defined(__MIMIC)
      IF (this%is_cmd(MCL_SEND_NUM_ATOMS)) THEN
         CALL this%send_value("nAtoms")
      ELSE IF (this%is_cmd(MCL_SEND_NUM_ATOM_TYPES)) THEN
         CALL this%send_value("nKinds")
      ELSE IF (this%is_cmd(MCL_SEND_ATOM_SPECIES)) THEN
         CALL this%send_atom_kinds()
      ELSE IF (this%is_cmd(MCL_SEND_MULTIPOLE_ORDER)) THEN
         CALL this%send_value("multipoleOrder")
      ELSE IF (this%is_cmd(MCL_SEND_NUM_FRAGMENTS)) THEN
         CALL this%send_value("nFragments")
      ELSE IF (this%is_cmd(MCL_SEND_NUM_ATOMS_IN_FRAGMENTS)) THEN
         CALL this%send_value("nAtomsInFragments")
      ELSE IF (this%is_cmd(MCL_SEND_NUM_BONDS)) THEN
         CALL this%send_value("nBonds")
      ELSE IF (this%is_cmd(MCL_SEND_NUM_ANGLES)) THEN
         CALL this%send_value("nAngles")
      ELSE IF (this%is_cmd(MCL_SEND_ATOM_MULTIPOLES)) THEN
         CALL this%send_multipoles()  ! fake placeholder (these are intended only for MM)
      ELSE IF (this%is_cmd(MCL_SEND_ATOM_MASSES)) THEN
         CALL this%send_kind_masses()
      ELSE IF (this%is_cmd(MCL_SEND_ATOMS_IN_FRAGMENTS)) THEN
         CALL this%send_atom_ids()
      ELSE IF (this%is_cmd(MCL_SEND_ATOM_ELEMENTS)) THEN ! element per kind
         CALL this%send_kind_elements()
      ELSE IF (this%is_cmd(MCL_SEND_ATOM_COORDS)) THEN
         CALL this%send_coords()
      ELSE IF (this%is_cmd(MCL_RECV_ATOM_COORDS)) THEN
         CALL this%receive_coords()
      ELSE IF (this%is_cmd(MCL_SEND_ATOM_FORCES)) THEN
         CALL this%send_forces()
      ELSE IF (this%is_cmd(MCL_SEND_ENERGY)) THEN
         CALL this%send_value("Energy")
      ELSE
         executed = .FALSE.
      END IF
#endif

      CALL timestop(handle)

   END SUBROUTINE execute_command


! **************************************************************************************************
!> \brief Sends one of the single value data to MiMiC server
!> \param option word corresponding to available options 
!> \note for now hardcoded num_frag = 1  =>  nAtInFrag = nAt
! **************************************************************************************************
   SUBROUTINE send_value(this, option)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      CHARACTER(LEN=*)                                   :: option

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_value'

      INTEGER                                            :: handle

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      SELECT CASE (option)
         CASE ("nAtoms", "nAtomsInFragments") 
            CALL mimic_send(this%n_atoms, 1, this%MCL_DATA, 0)
         CASE ("nKinds")
            CALL mimic_send(this%n_kinds, 1, this%MCL_DATA, 0)
         CASE ("multipoleOrder")
            CALL mimic_send(this%multipole_order, 1, this%MCL_DATA, 0)
         CASE ("nFragments")
            CALL mimic_send(1, 1, this%MCL_DATA, 0)
         CASE ("nBonds") ! later use to communicate constraints
            CALL mimic_send(0, 1, this%MCL_DATA, 0)
         CASE ("nAngles") ! later use to communicate constraints
            CALL mimic_send(0, 1, this%MCL_DATA, 0)
         CASE ("Energy")
            CALL mimic_send(this%energy%total, 1, this%MCL_DATA, 0)
         CASE DEFAULT
            CPABORT("The value chosen in "//routineN//" is not implemented.")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE send_value


! **************************************************************************************************
!> \brief Send a list of kind ids for each atom to the server
! **************************************************************************************************
   SUBROUTINE send_atom_kinds(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_atom_kinds'

      INTEGER, DIMENSION(:), ALLOCATABLE                 :: buffer
      INTEGER                                            :: handle, i

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      ALLOCATE(buffer(this%n_atoms))
      DO i=1,this%n_atoms
         buffer(i) = this%atoms%els(i)%atomic_kind%kind_number
      END DO
      CALL mimic_send(buffer, this%n_atoms, this%MCL_DATA, 0)
      DEALLOCATE(buffer)

      CALL timestop(handle)

   END SUBROUTINE send_atom_kinds


! **************************************************************************************************
!> \brief TBD
! **************************************************************************************************
   SUBROUTINE send_multipoles(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_multipoles'

      REAL(dp), DIMENSION(:), ALLOCATABLE                :: buffer
      INTEGER                                            :: handle


      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      ALLOCATE(buffer(this%n_atoms))
      buffer = 0.0_dp ! placeholder
      CALL mimic_send(buffer, this%n_atoms, this%MCL_DATA, 0)
      DEALLOCATE(buffer)

      CALL timestop(handle)

   END SUBROUTINE send_multipoles


! **************************************************************************************************
!> \brief Send a list of masses for each kind to the server
! **************************************************************************************************
   SUBROUTINE send_kind_masses(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_kind_masses'

      REAL(dp), DIMENSION(:), ALLOCATABLE                :: buffer
      INTEGER                                            :: handle, i

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      ALLOCATE(buffer(this%n_kinds))
      DO i=1,this%n_kinds
         buffer(i) = cp_unit_from_cp2k(this%kinds%els(i)%mass, "AMU")
      END DO
      CALL mimic_send(buffer, this%n_kinds, this%MCL_DATA, 0)
      DEALLOCATE(buffer)

      CALL timestop(handle)

   END SUBROUTINE send_kind_masses


! **************************************************************************************************
!> \brief Send a list of atom ids to the server
! **************************************************************************************************
   SUBROUTINE send_atom_ids(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_atom_ids'

      INTEGER, DIMENSION(:), ALLOCATABLE                 :: buffer
      INTEGER                                            :: handle, i

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      ALLOCATE(buffer(this%n_atoms))
      DO i=1,this%n_atoms
         buffer(i) = this%atoms%els(i)%atom_index
      END DO
      CALL mimic_send(buffer, this%n_atoms, this%MCL_DATA, 0)
      DEALLOCATE(buffer)

      CALL timestop(handle)

   END SUBROUTINE send_atom_ids


! **************************************************************************************************
!> \brief Send a list of atomic numbers for each kind to the server
! **************************************************************************************************
   SUBROUTINE send_kind_elements(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_kind_elements'

      INTEGER, DIMENSION(:), ALLOCATABLE                 :: buffer
      INTEGER                                            :: handle, i

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      ALLOCATE(buffer(this%n_kinds))
      DO i=1,this%n_kinds
         CALL get_atomic_kind(this%kinds%els(i), z=buffer(i))
      END DO
      CALL mimic_send(buffer, this%n_kinds, this%MCL_DATA, 0)
      DEALLOCATE(buffer)

      CALL timestop(handle)

   END SUBROUTINE send_kind_elements


! **************************************************************************************************
!> \brief Send coordinates of all atoms to the server
! **************************************************************************************************
   SUBROUTINE send_coords(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_coords'

      REAL(dp), DIMENSION(:), ALLOCATABLE                :: buffer
      INTEGER                                            :: handle, i

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      ALLOCATE(buffer(3*this%n_atoms))
      DO i=1,this%n_atoms
         buffer((3*i-2):(3*i)) = this%atoms%els(i)%r
      END DO
      CALL mimic_send(buffer, 3*this%n_atoms, this%MCL_DATA, 0)
      DEALLOCATE(buffer)

      CALL timestop(handle)

   END SUBROUTINE send_coords


! **************************************************************************************************
!> \brief Receive coordinates of all atoms from the server
! **************************************************************************************************
   SUBROUTINE receive_coords(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':receive_coords'

      REAL(dp), DIMENSION(:), ALLOCATABLE                :: buffer
      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      ALLOCATE(buffer(3*this%n_atoms))
      IF (this%is_ionode) THEN
         CALL mimic_receive(buffer, 3*this%n_atoms, this%MCL_DATA, 0)
      END IF
      CALL this%para_env%bcast(buffer)
      DO i=1,this%n_atoms
         this%atoms%els(i)%r = buffer((3*i-2):(3*i))
      END DO
      DEALLOCATE(buffer)

      CALL timestop(handle)

   END SUBROUTINE receive_coords


! **************************************************************************************************
!> \brief Send QM forces of all atoms to the server
! **************************************************************************************************
   SUBROUTINE send_forces(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_forces'

      REAL(dp), DIMENSION(:), ALLOCATABLE                :: buffer
      INTEGER                                            :: handle, i

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      ALLOCATE(buffer(3*this%n_atoms))
      DO i=1,this%n_atoms
         buffer((3*i-2):(3*i)) = this%atoms%els(i)%f
      END DO
      CALL mimic_send(buffer, 3*this%n_atoms, this%MCL_DATA, 0)
      DEALLOCATE(buffer)

      CALL timestop(handle)

   END SUBROUTINE send_forces

END MODULE mimic_communicator

