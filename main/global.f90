    !********************************************************************************
    !
    !  MODULE:  global
    !
    !  PURPOSE: Keeps all global variables and constants
    !
    !  SUBROUTINES:
    !               allocate_problem: allocates space
    !               free_memory: deallocates space
    !
    !******************************************************************************** 
    MODULE global
        
        USE sensitivity, only: sensitivity_type
        USE montecarlo, only: nd_type
        USE io, only: inputoutput_type, limits_type, LSSS_type
        
        IMPLICIT NONE
        SAVE
        
        !================================================================================
        ! Declare Parameters
        !================================================================================
        REAL(8), PARAMETER :: E=2.71828, PI=3.141592, gravity=9.81, T_REF=500.0

        !================================================================================
        ! Declare Variables
        !================================================================================
        type(sensitivity_type)          :: sens         ! sensitivities of flibe properties
        type(nd_type)                   :: nd           ! norm dist for flibe properties
        type(inputoutput_type)          :: inputoutput  ! inputs and outputs
        type(limits_type)               :: limits       ! LSSS limits
        type(LSSS_type)                 :: LSSS         ! LSSS results

        !type(tally_type),allocatable    :: tal(:)       ! tally
        !integer(8)                      :: npart        ! number of particles/runs
        
  
        contains

        !==============================================================================  
        ! allocate
        !==============================================================================  
        subroutine allocate_problem()
            implicit none
    
            ! allocate tal for number of slabs
            !allocate(tal(geo%nsubslabs))
    
        end subroutine allocate_problem
  
        !==============================================================================  
        ! deallocate
        !==============================================================================  
        subroutine free_memory()
            implicit none
    
            ! deallocate tal from memory
            !deallocate(tal)
  
        end subroutine free_memory
    
END MODULE global