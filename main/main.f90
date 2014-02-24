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
    use sensitivity, ONLY: sensitivity_init!, sensitivity_study
    use prismatic, ONLY: prismaticcoreTin,prismaticcoreTout,prismaticLSSS
    use flibeprop, ONLY: flibe_rho
    use trisoprop, ONLY: k_print
    use random_numbers, ONLY: rn_normal
        
            
    IMPLICIT NONE
        
        ! DECLARE VARIABLES
        integer                 :: i            ! Step counter
        real                    :: t1           ! Time 1
        real                    :: t2           ! Time 2 
        integer                 :: t            ! Time for system clock
        integer                 :: tt           ! Time for system clock
        integer                 :: clock_rate   ! clock rate
        integer                 :: clock_max    ! clock max
       
        
        
        !================================================================================
        ! Start Calculations
        !================================================================================
        write(*,*)
        write(*,*) "To Begin Program, hit enter:"
        read(*,*)
        write(*,*) "Calculating..."

        !================================================================================
        ! Open files for output
        !================================================================================       
        open(10,FILE="output.txt")
        !open(20,FILE="LSSS.txt")
        
        
        !================================================================================
        ! CALL INPUT SUBROUTINES IN MODULE IO
        !================================================================================
        call sensitivity_init(sens)
        call read_limits(limits)
        call LSSS_init(LSSS)
        
        
        !call read_inputs(input)
        !call outputs_init(output)

        
        !================================================================================
        ! Sensitivity Study Calculations
        !================================================================================ 
        !!!open(60,FILE="sensitivity.txt")
        !!!write(60,*) "Sensitivity Study"
        !!! Initialize sensitivities to 1.0
        !!!call sensitivity_init(sens)
        !!! Run sensitivity loop
        !!!call sensitivity_study(0.9_8, 0.01_8, 1.1_8, 4, sens, input, limits, LSSS) !2=cp, 3=k, 4=rho, 5=mu
        !!!close(60)
        

        !================================================================================
        ! Calculates LSSS for prismatic core
        !================================================================================
        call CPU_TIME(t1)
        call SYSTEM_CLOCK(t,clock_rate,clock_max)
        !!!open(50,FILE="histories.txt")
        !!!do i=1,10
        ! reset sensitivities
        call sensitivity_init(sens)
        call prismaticLSSS(inputoutput,limits,LSSS)
        call print_LSSS(LSSS,inputoutput%Q_core,limits)
        !!!end do
        !!!close(50)
        call CPU_TIME(t2)
        write ( *, * ) 'Elapsed CPU time = ', t2 - t1
        call SYSTEM_CLOCK(tt,clock_rate,clock_max)
        write ( *, * ) 'Elapsed real time = ', real ( tt - t ) / real ( clock_rate )
        write(*,*)
        write(*,*)
        write(*,*)
        
        
        !================================================================================
        ! Calculate LSSS for single pt for in avg and hot channels given Tin
        !================================================================================
        call CPU_TIME(t1)
        call SYSTEM_CLOCK(t,clock_rate,clock_max)
        ! Initialize inputs
        call inputoutput_init_Tin(inputoutput)
        ! Overide intialized arguments
        inputoutput%T_in=600.0_8
        !write(10,'(A,F6.2)') "SINGLE POWER LEVEL OUTPUT OF ", inputoutput%POWER/1.0E6
        write(*,*)
        write(*,*) "Calculating..."
        write(*,*)
        write(*,'(A,F6.2)') "SINGLE POWER LEVEL OUTPUT OF ", inputoutput%POWER/1.0E6
        write(*,*)
        write(*,*) "           POWER        W        Q       In      Out     Cool     Fuel"
        call prismaticcoreTin(inputoutput,1,1) 
        write(*,'(A,E9.2,F9.2,F9.4,F9.3,F9.3,F9.3,F9.3)') "AVG CH: ", inputoutput         
        call prismaticcoreTin(inputoutput,2,0)
        write(*,'(A,E9.2,F9.2,F9.4,F9.3,F9.3,F9.3,F9.3)') "HOT CH: ", inputoutput  
        call CPU_TIME(t2)
        write ( *, * ) 'Single Point: Elapsed CPU time = ', t2 - t1
        call SYSTEM_CLOCK(tt,clock_rate,clock_max)
        write ( *, * ) 'Single Point: Elapsed real time = ', real ( tt - t ) / real ( clock_rate )
        
        
        !================================================================================
        ! Calculate LSSS for single pt for in avg and hot channels given Toutavg
        !================================================================================
        call CPU_TIME(t1)
        call SYSTEM_CLOCK(t,clock_rate,clock_max)
        ! Initialize inputs
        call inputoutput_init_Tout(inputoutput)
        ! Overide intialized arguments
        inputoutput%T_out=700.002_8
        !write(10,'(A,F6.2)') "SINGLE POWER LEVEL OUTPUT OF ", inputoutput%POWER/1.0E6
        write(*,*)
        write(*,*) "Calculating..."
        write(*,*)
        write(*,'(A,F6.2)') "SINGLE POWER LEVEL OUTPUT OF ", inputoutput%POWER/1.0E6
        write(*,*)
        write(*,*) "           POWER        W        Q       In      Out     Cool     Fuel"
        call prismaticcoreTout(inputoutput,1,1) 
        write(*,'(A,E9.2,F9.2,F9.4,F9.3,F9.3,F9.3,F9.3)') "AVG CH: ", inputoutput 
        call prismaticcoreTin(inputoutput,2,1)
        write(*,'(A,E9.2,F9.2,F9.4,F9.3,F9.3,F9.3,F9.3)') "HOT CH: ", inputoutput  
        call CPU_TIME(t2)
        write ( *, * ) 'Single Point: Elapsed CPU time = ', t2 - t1
        call SYSTEM_CLOCK(tt,clock_rate,clock_max)
        write ( *, * ) 'Single Point: Elapsed real time = ', real ( tt - t ) / real ( clock_rate )
        
        
 
        
        !================================================================================
        ! Close files for output
        !================================================================================                            
        close(10)
        !close(20)
        
        
        !================================================================================
        ! Wait for user
        !================================================================================
        write(*,*)
        write(*,*) 'Program Complete'
        read(*,*)


    CONTAINS


END PROGRAM main

    
        !================================================================================
        ! Print Thermal Conductivities
        !================================================================================
        !call k_print()
        

        !================================================================================
        ! Print Normal Distribution
        !================================================================================
        !open(40,FILE="norm.txt") 
        !do i=1,1000000
        !    write(40,'(F9.2)') rn_normal(2386.0d0,71.58d0)  
        !end do
        !close(40)
    
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