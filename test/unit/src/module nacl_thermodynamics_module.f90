module nacl_thermodynamics_module_2

#include <petsc/finclude/petscsys.h>
        
        use petscsys
        use kinds_module
        use powertable_module
        use thermodynamics_module
        use pure_nacl
                
        implicit none
        private

        type, public :: vapor_liquid
            private
            PetscReal :: critical_temperature_water = 373.976_dp
            PetscReal :: critical_pressure_water = 220.54915_dp            
            contains
                private
                procedure, public :: composition => liquid_composition
        end type vapor_liquid
        
        type, public :: critical_fraction
            private
            character(len = 6) : xcrit
            character(len = 6) : pcrit
            contains
                private
                procedure, public :: xcrit => critical_composition
                procedure, public :: pcrit => critical_pressure
                end type critical_fraction

        type, public :: vapor_liquid_halite
            private 
            class(pure_nacl), allocatable, private :: nacl
            PetscReal :: f(10) = (/4.64E-3, 5E-7, 16.9078, -2.69148E2, 7.63204E3, &
                -4.95636E4, 2.33119E5, -5.13556E5, 5.49708E5, -2.84628E5/) 
            PetscInt :: i(11) = (/ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 /)
            type (powertable_type) :: ti
            contains
                private
                procedure, public :: init => vapor_liquid_halite_init
                procedure, public :: destroy => vapor_liquid_halite_destroy
                procedure, public :: pressure => vapor_liquid_halite_pressure
        end type vapor_liquid_halite

        type, public :: liquid_branch
            private
            class(pure_nacl), allocatable, private :: nacl
            PetscReal :: f(11) = (/1.68486E-3, 2.19379E-4, 4.3858E2, 18.4508, -5.6765E-10, &
                6.73704E-6, 1.44951E-7, 3.84904E2, 7.07477, 6.06896E-5, 7.62859E-3/)
            PetscInt :: i(12) = (/ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11/)
            type (powertable_type) :: ti
            type (variable)
            contains
            private
            procedure, public :: init => liquid_branch_init
            procedure, public :: destroy => liquid_branch_destroy
            procedure, public :: fraction => liquid_branch_fraction
        end type liquid_branch

    contains
    
        subroutine vapor_liquid_halite_pressure(self, temperature, pressure, err)
            use utils_module, only: polynomial

            class(vapor_liquid_halite),intent(inout) :: self
            PetscReal, intent(in) :: temperature
            PetscReal, intent(out) :: pressure
            PetscErrorCode, intent(out) :: err

            err = 0

            f10 = self%nacl%triple_point_pressure - sum(self%f)
    
            call self%ti%compute(temperature / self%nacl%triple_point_temperature)

            pressure = sum(self%f * self%ti%power(self%i(:10)))     
        end subroutine

        subroutine vapor_liquid_halite_init(self)
            class(vapor_liquid_halite), intent(inout) :: self
    
            allocate(pure_nacl :: self%nacl)

            call self%ti%configure(self%i)
        end subroutine vapor_liquid_halite_init
    
        subroutine vapor_liquid_halite_destroy(self)
            class(vapor_liquid_halite), intent(inout) :: self
    
            deallocate(self%nacl)

            call self%ti%destroy()
        end subroutine vapor_liquid_halite_destroy

        subroutine liquid_branch_fraction(self, pressure, fraction, err)
            use utils_module, only: polynomial

            class(liquid_branch), intent(inout) :: self
            PetscReal, intent(in) :: pressure
            PetscReal, intent(out) :: fraction
            PetscErrorCode, intent(out) :: error
            PetscReal :: critical_pressure_water
            PetscReal :: critical_temperature_water
           
            liquid_composition = xcrit + (g0 * sqrt(critical_pressure_water - pressure)) &
                + (g1 * (critical_pressure_water - pressure)) & 
                + (g2 * (critical_pressure_water - pressure)E2)

            if ( pressure =  ) then
                
            end if
                if ( T<critical_temperature_water ) then pcrit
                
            end if    

        end subroutine

        subroutine liquid_branch_init(Self)
            class(liquid_branch), intent(inout) :: self

            allocate (pure_nacl :: self%nacl)

            call self%ti%configure(self%i)
        end subroutine liquid_branch_init

        subroutine liquid_branch_destroy
            class(liquid_branch), intent(inout) :: self

            deallocate(self%nacl)

            call self%ti%destroy
        end subroutine liquid_branch_destroy
    end module nacl_thermodynamics_module