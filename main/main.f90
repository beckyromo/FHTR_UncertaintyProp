!****************************************************************************
!
!  PROGRAM: main
!
!  PURPOSE: Entry point for the main application.
!
!**************************************************************************** 
PROGRAM main    
    
    use global
    use io
    use sensitivity, ONLY: sensitivity_init, sensitivity_study
    USE pebblebed, ONLY: pebblebedLSSSloop,pebblebedcore
    USE prismatic, ONLY: prismaticLSSSloop,prismaticcore
    USE flibeprop, ONLY: flibe_rho
    USE trisoprop, ONLY: k_print
    !use PDFs, ONLY: get_normdist_value
    use random_numbers, ONLY: normal
        
            
    IMPLICIT NONE
        
        ! DECLARE VARIABLES
        integer             :: i  
       
        
        
        !================================================================================
        ! Print Thermal Conductivities
        !================================================================================
        !call k_print()
        

        !================================================================================
        ! Print Normal Distribution
        !================================================================================
        open(40,FILE="norm.txt") 
        do i=1,1000000
            write(40,'(F9.2)') normal(2386.0d0,71.58d0)  
        end do
        close(40)
      
        
        !================================================================================
        ! CALL INPUT SUBROUTINES IN MODULE IO
        !================================================================================
        call sensitivity_init(sens)
        call read_inputs(input)
        call read_limits(limits)
        call outputs_init(output)
        call LSSS_init(LSSS)
        

        !================================================================================
        ! Start Calculations
        !================================================================================
        write(*,*)
        write(*,*) "To Begin Program, hit enter:"
        read(*,*)
        write(*,*) "Calculating..."

        
     
        open(10,FILE="output.txt")
        open(20,FILE="LSSS.txt")
        
        
        
        !================================================================================
        ! Sensitivity Study Calculations
        !================================================================================ 
        !open(60,FILE="sensitivity.txt")
        !write(60,*) "Sensitivity Study"
        ! Initialize sensitivities to 1.0
        !call sensitivity_init(sens)
        ! Run sensitivity loop
        !call sensitivity_study(0.9_8, 0.01_8, 1.1_8, 4, sens, input, limits, LSSS) !2=cp, 3=k, 4=rho, 5=mu
        !close(60)
        

        !================================================================================
        ! Calculates LSSS for prismatic core
        !================================================================================
        ! reset sensitivities
        !call sensitivity_init(sens)
        !call prismaticLSSSloop(input,limits,LSSS)
        !call print_LSSS(LSSS,input%Q_core,limits)
        

        
        
        !================================================================================
        ! Calculates LSSS for single point for prismatic core in avg and hot channels
        !================================================================================
        !input%T_in=600.0_8
        !input%POWER=LSSS%minPOWER*1.0E6
        !write(10,'(A,F6.2)') "SINGLE POWER LEVEL OUTPUT OF ", input%POWER/1.0E6
        !write(*,*)
        !write(*,*) "Calculating..."
        !write(*,'(A,F6.2)') "SINGLE POWER LEVEL OUTPUT OF ", input%POWER/1.0E6
        !write(*,*)
        !write(*,*) "               Q       In      Out     Cool     Fuel"
        !call prismaticcore(input%POWER,input%Q_core,input%T_in,output,1) 
        !write(*,'(A,F9.4,F9.3,F9.3,F9.3,F9.3)') "AVG CH: ", input%Q_core, input%T_in, output         
        !call prismaticcore(input%POWER,input%Q_core,input%T_in,output,2)
        !write(*,'(A, F9.4,F9.3,F9.3,F9.3,F9.3)') "HOT CH: ", input%Q_core, input%T_in, output   
        
                                    
        close(10)
        close(20)
        
        write(*,*)
        write(*,*) 'Program Complete'
               
        
        !================================================================================
        ! Wait for user
        !================================================================================
        ! read(*,*)


    CONTAINS


    END PROGRAM main

    
        !! ===============================================================================
        !!  PEBBLE BED PEBBLE BED PEBBLE BED PEBBLE BED PEBBLE BED PEBBLE BED PEBBLE BED
        !! -------------------------------------------------------------------------------
        !! -------------------------------------------------------------------------------
        !
        !!================================================================================
        !! INPUTS
        !!================================================================================
        !POWER=20.0E6_8              ! [W]
        !Q_core=84.6_8/flibe_rho(600.0_8)*0.9              ! [m^3/s]
        !T_in=500.0_8                ! [Celcius]
        !T_fuel_limit=1300.0_8       ! [Celcius]
        !T_coolant_limit=1200.0_8    ! [Celcius]
        !T_out_limit=700.0_8         ! [Celcius]
        !
        !
        !!================================================================================
        !! Start Calculations
        !!================================================================================
        !write(*,*)
        !write(*,*) 'To Begin Program, hit enter:'
        !read(*,*)
        !
        !open(10,FILE="output.txt")
        !
        !open(20,FILE="LSSS.txt")
        !
        !!================================================================================
        !! Calculates LSSS for pebblebed core
        !!================================================================================
        !call pebblebedLSSSloop(Q_core,T_in,T_fuel_limit,T_coolant_limit,T_out_limit,T_out,T_coolant_max,T_core_max)
        !
        !!================================================================================
        !! Calculates LSSS for single point for pebblebed core
        !!================================================================================
        !POWER=44.28E6_8              ! [W]
        !Q_core=84.6_8/flibe_rho(600.0_8)              ! [m^3/s]
        !write(*,*) Q_core*3600.0
        !T_in=466.0_8                ! [Celcius]
        !write(10,*) "SINGLE POWER LEVEL OUTPUT OF ", POWER
        !write(*,*) "SINGLE POWER LEVEL OUTPUT OF ", POWER
        !write(*,*)
        !write(*,*) "               Q       In      Out     Cool     Fuel"
        !call pebblebedcore(POWER,Q_core,T_in,T_out,T_coolant_max,T_core_max,1)
        !write(*,'(A,F9.2,F9.3,F9.3,F9.3,F9.3)') "AVG CH: ", Q_core*3600.0, T_in, T_out, T_coolant_max, T_core_max         
        !call pebblebedcore(POWER,Q_core,T_in,T_out,T_coolant_max,T_core_max,2)
        !write(*,'(A, F9.2,F9.3,F9.3,F9.3,F9.3)') "HOT CH: ", Q_core*3600.0, T_in, T_out, T_coolant_max, T_core_max