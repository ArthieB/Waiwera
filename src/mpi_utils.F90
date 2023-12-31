!   Copyright 2016 University of Auckland.

!   This file is part of Waiwera.

!   Waiwera is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.

!   Waiwera is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.

!   You should have received a copy of the GNU Lesser General Public License
!   along with Waiwera.  If not, see <http://www.gnu.org/licenses/>.

module mpi_utils_module
  !! MPI utilities.

#include <petsc/finclude/petscsys.h>

  use petscsys

  implicit none
  private

  public :: mpi_broadcast_error_flag, mpi_broadcast_logical
  public :: get_mpi_int_gather_array

contains

!------------------------------------------------------------------------

  subroutine mpi_broadcast_error_flag(err)
    !! If integer error flag err is greater than zero on any rank,
    !! broadcast value 1 to all ranks.

    PetscErrorCode, intent(in out) :: err
    ! Locals:
    PetscBool :: any_err
    PetscErrorCode :: ierr

    call MPI_allreduce(err > 0, any_err, 1, MPI_LOGICAL, MPI_LOR, &
         PETSC_COMM_WORLD, ierr)
    if (any_err) then
       err = 1
    else
       err = 0
    end if

  end subroutine mpi_broadcast_error_flag

!------------------------------------------------------------------------

  subroutine mpi_broadcast_logical(val)
    !! If logical val is true on any rank, broadcast to all ranks.

    PetscBool, intent(in out) :: val
    ! Locals:
    PetscBool :: any_val
    PetscErrorCode :: ierr

    call MPI_allreduce(val, any_val, 1, MPI_LOGICAL, MPI_LOR, &
         PETSC_COMM_WORLD, ierr)
    val = any_val

  end subroutine mpi_broadcast_logical

!------------------------------------------------------------------------

  function get_mpi_int_gather_array() result(array)
    !! Returns integer array for use in MPI gather call. This is of
    !! size equal to the number of processes on the root rank, and
    !! size 1 on other ranks (needs to be allocated even though it is
    !! not actually used.)

    PetscInt, allocatable :: array(:)
    ! Locals:
    PetscMPIInt :: rank, num_procs
    PetscInt :: size
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(PETSC_COMM_WORLD, num_procs, ierr)

    if (rank == 0) then
       size = num_procs
    else ! have to allocate non-zero size, even if not actually used:
       size = 1
    end if
    allocate(array(size))

  end function get_mpi_int_gather_array

!------------------------------------------------------------------------

end module mpi_utils_module
