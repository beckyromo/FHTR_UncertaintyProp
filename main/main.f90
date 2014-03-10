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
    use sensitivity, ONLY: sensitivity_init,sensitivity_study
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
        ! CALL INPUT SUBROUTINES IN MODULE IO FOR INTIALIZATION
        !================================================================================
        call sensitivity_init(sens)
        call read_limits(limits)
        call LSSS_init(LSSS)
        call inputoutput_init_inputs(inputoutput)
        call inputoutput_init_Tin(inputoutput)

        
        !================================================================================
        ! Sensitivity Study Calculations
        !================================================================================ 
        open(60,FILE="sensitivity.txt")
        write(60,*) "Sensitivity Study"
        ! Initialize sensitivities to 1.0
        call sensitivity_init(sens)
        ! Run sensitivity loop
        call sensitivity_study(0.9_8, 0.01_8, 1.1_8, 4, sens, inputoutput, limits, LSSS) !2=cp, 3=k, 4=rho, 5=mu
        close(60)
        

        !================================================================================
        ! Calculates LSSS for prismatic core
        !================================================================================
        call CPU_TIME(t1)
        call SYSTEM_CLOCK(t,clock_rate,clock_max)
        !!!open(50,FILE="histories.txt")
        !!!do i=1,10
            ! reset sensitivities
            call sensitivity_init(sens)
            ! initialize inputs and outputs
            call inputoutput_init_Tin(inputoutput)
            ! Run LSSS
            call prismaticLSSS(inputoutput,limits,LSSS)
            ! Print LSSS
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
        ! Initialize inputs and outputs
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
        ! Initialize inputs and outputs
        call inputoutput_init_Tout(inputoutput)
        ! Overide intialized arguments
        inputoutput%T_out=720.002_8
        inputoutput%POWER=21.13E6
        sens%s_rho=0.9_8
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
    
