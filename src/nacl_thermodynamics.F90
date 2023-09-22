module nacl_thermodynamics_module

#include <petsc/finclude/petscsys.h>
    
        use petscsys
        use kinds_module
        use powertable_module
        use thermodynamics_module
    
        implicit none
        private

        type, public :: pure_nacl
            private
            PetscReal :: triple_point_pressure = 5e-4_dp
            PetscReal :: triple_point_temperature = 800.7_dp
            PetscReal :: a = 2.47260e-2_dp
            contains
                private
                procedure, public :: melting => pure_nacl_melting_temperature
        end type pure_nacl

        type, public :: halite_liquidus
            private
            class(pure_nacl), allocatable, private :: nacl
            PetscReal :: e(5, 3) = reshape([ &
                0.0989944_dp, 0.00947257_dp, 0.610863_dp, -1.64994_dp, 3.36474_dp, &
                3.30796e-6_dp, -8.66460e-6_dp, -1.51716e-5_dp, 2.03441e-4_dp, -1.54023e-4_dp, &
                -4.71759e-10_dp, 1.69417e-9_dp, 1.19290e-8_dp, -6.46015e-8_dp, 8.17048e-8_dp &
                ], &
                [5, 3])
            PetscInt :: i(6) = (/ 0, 1, 2, 3, 4, 5 /)
            type(powertable_type) :: ti
            contains
                private
                procedure, public :: init => halite_liquidus_init
                procedure, public :: destroy => halite_liquidus_destroy
                procedure, public :: weight_fraction => nacl_weight_fraction
                procedure, public :: fraction => halite_liquidus_fraction
        end type halite_liquidus
       
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

    contains

        ! TODO: WHERE TO PUT THIS?
        subroutine nacl_weight_fraction(self, mole_fraction, weight_fraction)
            class(halite_liquidus), intent(inout) :: self
            PetscReal, intent(in) :: mole_fraction
            PetscReal, intent(out) :: weight_fraction
            ! Local variables:
            PetscReal :: nacl_molecular_weight = 58.443_dp
            PetscReal :: h2o_molecular_weight = 18.015_dp
            PetscReal :: weight

            weight = (1 - mole_fraction) * h2o_molecular_weight + mole_fraction * nacl_molecular_weight
            weight_fraction = mole_fraction * (nacl_molecular_weight / weight)
        end subroutine

        subroutine halite_liquidus_init(self)
            class(halite_liquidus), intent(inout) :: self
    
            allocate(pure_nacl :: self%nacl)

            call self%ti%configure(self%i)
        end subroutine halite_liquidus_init
    
        subroutine halite_liquidus_destroy(self)
            class(halite_liquidus), intent(inout) :: self
    
            deallocate(self%nacl)

            call self%ti%destroy()
        end subroutine halite_liquidus_destroy

        subroutine pure_nacl_melting_temperature(self, pressure, temperature, err)
            class(pure_nacl), intent(inout) :: self
            PetscReal, intent(in) :: pressure
            PetscReal, intent(out) :: temperature
            PetscErrorCode, intent(out) :: err

            err = 0

            temperature = self%triple_point_temperature + self%a * (pressure - self%triple_point_pressure)
        end subroutine
    
        subroutine halite_liquidus_fraction(self, temperature, pressure, fraction, err)
            use utils_module, only: polynomial
    
            class(halite_liquidus), intent(inout) :: self
            PetscReal, intent(in) :: temperature
            PetscReal, intent(in) :: pressure
            PetscReal, intent(out) :: fraction
            PetscErrorCode, intent(out) :: err !! Error code
            ! Local variable:
            PetscReal :: melting_temperature
            PetscReal :: e5
    
            e5 = 1.0_dp - sum(polynomial(self%e, pressure))

            call self%nacl%melting(pressure, melting_temperature, err)
            call self%ti%compute(temperature / melting_temperature)

            fraction = sum(polynomial(self%e, pressure) * self%ti%power(self%i(:5)))
            fraction = fraction + e5 * self%ti%power(self%i(6))
        end subroutine halite_liquidus_fraction 
        
        subroutine vapor_liquid_halite_pressure(self, temperature, pressure, err)
            use utils_module, only: polynomial

            class(vapor_liquid_halite),intent(inout) :: self
            PetscReal, intent(in) :: temperature
            PetscReal, intent(out) :: pressure
            PetscErrorCode, intent(out) :: err
            PetscReal :: f(10)
            
            err = 0 

            f(10) = self%nacl%triple_point_pressure - sum(self%f)
    
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
    end module nacl_thermodynamics_module