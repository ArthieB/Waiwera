module nacl_thermodynamics_test 
#include <petsc/finclude/petscsys.h>

        use petscsys
        use kinds_module
        use nacl_thermodynamics_module
        use zofu
        
        implicit none
        private
        public :: setup, teardown
        public :: test_halite_liquidus_fraction
            
    contains
        subroutine setup()
            use profiling_module, only: init_profiling
        
            PetscErrorCode :: ierr
        
            call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
            call init_profiling()
        
        end subroutine setup
        
        subroutine teardown()
        
            PetscErrorCode :: ierr
        
            call PetscFinalize(ierr); CHKERRQ(ierr)
        
        end subroutine teardown
            
        subroutine test_halite_liquidus_fraction(test)
                
            class(unit_test_type), intent(inout) :: test
            type(halite_liquidus) :: liquidus
            PetscReal :: temperature, pressure, mole_fraction, weight_fraction, expected
            PetscErrorCode :: err
            character(34) :: s = "Halite liquidus fraction, pressure"
            PetscMPIInt :: rank
            PetscInt :: ierr
        
            call liquidus%init()
        
            call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
        
            if (rank == 0) then
    
                pressure = 698_dp
                temperature = 623_dp
                expected = 0.7829_dp
                call liquidus%fraction(temperature, pressure, mole_fraction, err)
                call liquidus%weight_fraction(mole_fraction, weight_fraction)
                call test%assert(0, err, trim(s) // " 698 bar, temp. 623 deg C error")
                call test%assert(expected, weight_fraction, trim(s) // " 698 bar, temp. 623 deg C", 2E-2_dp)
    
                pressure = 467_dp
                temperature = 624_dp
                expected = 0.7829_dp
    
                call liquidus%fraction(temperature, pressure, mole_fraction, err)
                call liquidus%weight_fraction(mole_fraction, weight_fraction)
                call test%assert(0, err, trim(s) // " 467 bar, temp. 624 deg C error")
                call test%assert(expected, weight_fraction, trim(s) // " 467 bar, temp. 624 deg C", 2E-2_dp)
    
                pressure = 394_dp
                temperature = 638_dp
                expected = 0.7829_dp
                    
                call liquidus%fraction(temperature, pressure, mole_fraction, err)
                call liquidus%weight_fraction(mole_fraction, weight_fraction)
                call test%assert(0, err, trim(s) // " 394 bar, temp. 638 deg C error")
                call test%assert(expected, weight_fraction, trim(s) // " 394 bar, temp. 638 deg C", 2E-2_dp)
    
                pressure = 803_dp
                temperature = 685_dp
                expected = 0.8647_dp
                    
                call liquidus%fraction(temperature, pressure, mole_fraction, err)
                call liquidus%weight_fraction(mole_fraction, weight_fraction)
                call test%assert(0, err, trim(s) // " 803 bar, temp. 685 deg C error")
                call test%assert(expected, weight_fraction, trim(s) // " 803 bar, temp. 685 deg C", 2E-2_dp)
    
                pressure = 380_dp
                temperature = 684_dp
                expected = 0.8647_dp
                    
                call liquidus%fraction(temperature, pressure, mole_fraction, err)
                call liquidus%weight_fraction(mole_fraction, weight_fraction)
                call test%assert(0, err, trim(s) // " 380 bar, temp. 684 deg C error")
                call test%assert(expected, weight_fraction, trim(s) // " 380 bar, temp. 684 deg C", 2E-2_dp)
            end if
        
        end subroutine test_halite_liquidus_fraction

        subroutine test_vapor_liquid_halite(test)
            
            class(unit_test_type), intent(inout) :: test
            type(vapor_liquid_halite) :: vlh
            PetscReal :: temperature, pressure, expected
            PetscErrorCode :: err
            character(34) :: s = "Vapor liquid halite, pressure"
            PetscMPIInt :: rank
            PetscInt :: ierr 
            
            call vlh%init()
            
            call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
            
            if (rank == 0) then
        
                pressure = 78.8_dp
                temperature = 323.7_dp
                expected = 0.0071_dp
                call vlh%fraction(temperature, pressure, err)
                call test%assert(0, err, trim(s) // " 78.8 bar, temp. 323.7 deg C error")
                call test%assert(expected, pressure, trim(s) // " 78.8 bar, temp. 323.7 deg C", 2E-2_dp)
        
                pressure = 105.3_dp
                temperature = 348.5_dp
                expected = 0.0014_dp
        
                call vlh%fraction(temperature, pressure, err)
                call test%assert(0, err, trim(s) // " 105.3 bar, temp. 348.5 deg C error")
                call test%assert(expected, pressure, trim(s) // " 105.3 bar, temp. 348.5 deg C", 2E-2_dp)
        
                pressure = 137.4_dp
                temperature = 375.1_dp
                expected = 0.0032_dp
                        
                call vlh%fraction(temperature, pressure, err)
                call test%assert(0, err, trim(s) // " 137.4 bar, temp. 375.1 deg C error")
                call test%assert(expected, pressure, trim(s) // " 137.4 bar, temp. 375.1 deg C", 2E-2_dp)
        
                pressure = 138.4_dp
                temperature = 375.5_dp
                expected = 0.0026_dp
                    
                call vlh%fraction(temperature, pressure, err)
                call test%assert(0, err, trim(s) // " 138.4 bar, temp. 375.5 deg C error")
                call test%assert(expected, pressure, trim(s) // " 138.4 bar, temp. 375.5 deg C", 2E-2_dp)
        
                pressure = 172.7_dp
                temperature = 400.2_dp
                expected = 0.0065_dp
                      
                call vlh%fraction(temperature, pressure, err)
                call test%assert(0, err, trim(s) // " 172.7 bar, temp. 400.2 deg C error")
                call test%assert(expected, pressure, trim(s) // " 172.7 bar, temp. 400.2 deg C", 2E-2_dp)
            end if
            
        end subroutine test_vapor_liquid_halite
    end module nacl_thermodynamics_test