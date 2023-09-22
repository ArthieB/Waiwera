module eos_nacl_module
#include <petsc/finclude/petscsys.h>
    
    use petscsys
    use kinds_module
    use eos_module

    implicit none
    private
    
    type, public, extends(eos_type) :: eos_nacl_type
        private 
    contains
        private
        procedure, public :: init => eos_nacl_init
        procedure, public :: destroy => eos_nacl_destroy
        procedure, public :: transition => eos_nacl_transition
        procedure, public :: bulk_properties => eos_nacl_bulk_properties
        procedure, public :: phase_properties => eos_nacl_phase_properties
        procedure, public :: primary_variables => eos_nacl_primary_variables
        procedure, public :: check_primary_variables => eos_nacl_check_primary_variables
    end type eos_nacl_type

contains
    
    subroutine eos_nacl_init(self, json, thermo, logfile)
        use fson
        use fson_mpi_module, only: fson_get_mpi
        use logfile_module
        use thermodynamics_module

        class(eos_nacl_type), intent(inout) :: self
        type(fson_value), pointer, intent(in) :: json
        class(thermodynamics_type), intent(in), target :: thermo
        type(logfile_type), intent(inout), optional :: logfile

        PetscReal, parameter :: default_pressure = 1.0e5_dp
        PetscReal, parameter :: default_temperature = 20._dp
        PetscReal, parameter :: default_nacl_fraction = 0._dp
        PetscReal, parameter :: default_pressure_scale = 1.e6_dp
        PetscReal, parameter :: default_temperature_scale = 1.e2_dp

        self%name = 'nacl'
        self%description = 'H2O-NaCl'
        self%primary_variable_names = ['Pressure', 'Temperature', 'NaCl fraction']
        self%num_primary_variables = size(self%primary_variable_names)
        self%num_phases = 3
        self%phase_names = ['liquid', 'vapour', 'halite']
        self%num_components = 2
        self%num_components = 2
        self%component_names = ['water', 'nacl']
        self%required_output_fluid_fields = [ &
            'pressure     ', 'temperature  ', &
            'nacl_fraction']
        self%default_output_fluid_fields = [ &
        'pressure     ', 'temperature  ', &
        'nacl_fraction']
        
        call fson_get_mpi(json, 'eos.primary.scale.pressure', default_pressure_scale, pressure_scale, logfile)
        call fson_get_mpi(json, 'eos.primary.scale.temperature', default_temperature_scale, temperature_scale, logfile)
        

        self%thermo => thermo
    end subroutine eos_nacl_init

    subroutine eos_nacl_destroy(self)
        class(eos_nacl_type), intent(inout) :: self
        deallocate(self%primary_variable_names)
        deallocate(self%phase_names, self%component_names)
        
        self%thermo => null()
        
    end subroutine eos_nacl_destroy

    subroutine eos_nacl_transition(self, old_primary, primary, old_fluid, fluid, transition, err)
        use fluid_module, only: fluid_type

        class(eos_nacl_type), intent(inout) :: self
        PetscReal, intent(in) :: old_primary(self%num_primary_variables)
        PetscReal, intent(inout) :: primary(self%num_primary_variables)
        type(fluid_type), intent(in) :: old_fluid
        type(fluid_type), intent(inout) :: fluid
        PetscBool, intent(out) :: transition
        PetscErrorCode, intent(out) :: err

    end subroutine eos_nacl_transition

    subroutine eos_nacl_bulk_properties(self, primary, fluid, err)
        use fluid_module, only: fluid_type

        class(eos_nacl_type), intent(inout) :: self
        PetscReal, intent(in) :: primary(self%num_primary_variables)
        type(fluid_type), intent(inout) :: fluid
        PetscErrorCode, intent(out) :: err

    end subroutine eos_nacl_bulk_properties

    subroutine eos_nacl_phase_properties(self, primary, rock, fluid, err)
        use rock_module, only: rock_type
        use fluid_module, only: fluid_type

        class(eos_nacl_type), intent(inout) :: self
        PetscReal, intent(in) :: primary(self%num_primary_variables)
        type(rock_type), intent(inout) :: rock
        type(fluid_type), intent(inout) :: fluid
        PetscErrorCode, intent(out) :: err

    end subroutine eos_nacl_phase_properties

    subroutine eos_nacl_primary_variables(self, fluid, primary)
        use fluid_module, only: fluid_type

        class(eos_nacl_type), intent(in) :: self
        type(fluid_type), intent(in) :: fluid
        PetscReal, intent(out) :: primary(self%num_primary_variables)
    end subroutine eos_nacl_primary_variables

    subroutine eos_nacl_check_primary_variables(self, fluid, primary, changed, err)
        use fluid_module, only: fluid_type

        class(eos_nacl_type), intent(in) :: self
        type(fluid_type), intent(in) :: fluid
        PetscReal, intent(inout) :: primary(self%num_primary_variables)
        PetscBool, intent(out) :: changed
        PetscErrorCode, intent(out) :: err
    end subroutine eos_nacl_check_primary_variables
end module eos_nacl_module