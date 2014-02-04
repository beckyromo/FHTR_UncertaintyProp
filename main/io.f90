module io
    use flibeprop, only:flibe_rho
    implicit none
    private
    public :: read_limits, read_inputs, outputs_init, LSSS_init, print_LSSS
  
    ! declare variables
    type,public :: input_type
        real(8)             :: POWER                ! Reactor power [W]
        !real(8)             :: DECAY_HEAT           ! Decay heat of reactor [W]
        !real(8)             :: POWER_NORM           ! Norminal Fisison Power [W]
        real(8)             :: W_core               ! Mass flow rate in the core [kg/s]
        real(8)             :: Q_core               ! Volumetric flow rate in the core [m^3/s]
        real(8)             :: T_in                 ! Inlet temperature of the core [Celcius] 
    end type input_type
  
    type,public :: limits_type
        real(8)             :: T_in_limit           ! LSSS temperature limit for T_in (min T_in temp)
        real(8)             :: T_fuel_limit         ! LSSS temoerature limit for fuel (max fuel temp in hot channel) 
        real(8)             :: T_coolant_limit      ! LSSS temperature limit for coolant (max coolant temp in hot channel)
        real(8)             :: T_out_limit          ! LSSS temperature limit for T_out (max T_out temp in average channel)
    end type limits_type
    
    type,public :: output_type
        real(8)             :: T_out                ! Maximum coolant temperature [Celcius]
        real(8)             :: T_coolant_max        ! Maximum core temperature  [Celcius]
        real(8)             :: T_core_max           ! Outlet temperature of the core [Celcius] 
    end type output_type  
    
    type,public :: LSSS_type
        real(8)             :: INmaxcoolPOWER
        real(8)             :: INmaxcoolTin
        real(8)             :: INmaxcoolTout
        real(8)             :: INmaxcoolTmax
        real(8)             :: INmaxcoolToutavg
        real(8)             :: INmaxfuelPOWER
        real(8)             :: INmaxfuelTin
        real(8)             :: INmaxfuelTout
        real(8)             :: INmaxfuelTmax
        real(8)             :: INmaxfuelToutavg
        real(8)             :: INmaxoutPOWER
        real(8)             :: INmaxoutTin
        real(8)             :: INmaxoutTout
        real(8)             :: OUTmaxcoolPOWER
        real(8)             :: OUTmaxcoolTin
        real(8)             :: OUTmaxcoolTout
        real(8)             :: OUTmaxcoolTmax
        real(8)             :: OUTmaxcoolToutavg
        real(8)             :: OUTmaxfuelPOWER
        real(8)             :: OUTmaxfuelTin
        real(8)             :: OUTmaxfuelTout
        real(8)             :: OUTmaxfuelTmax
        real(8)             :: OUTmaxfuelToutavg
        real(8)             :: minPOWER
    end type LSSS_type 
    
contains

!==============================================================================
! read_limits
!==============================================================================
    subroutine read_limits(this)
        ! declare arguments
        type(limits_type)   :: this
    
        this%T_fuel_limit=1300.0_8               ! [Celcius]
        this%T_coolant_limit=1200.0_8            ! [Celcius]
        this%T_in_limit=470.0_8                  ! [Celcius]
        this%T_out_limit=720.0_8                 ! [Celcius]
        
    end subroutine read_limits
  
    
!==============================================================================
! read_inputs
!==============================================================================
    subroutine read_inputs(this)
        ! declare arguments
        type(input_type)    :: this
        
        this%POWER=20.0E6_8                      ! [W]
        this%W_core=83.82_8                      ! [kg/s]
        this%Q_core=this%W_core/flibe_rho(600.0_8)    ! [m^3/s]
        this%T_in=470.0_8                        ! [Celcius]
    
    end subroutine read_inputs
    
    
!==============================================================================
! outputs_init
!==============================================================================    
    subroutine outputs_init(this)
        !declare arguemtns
        type(output_type)   :: this
        
        this%T_out=0.0_8
        this%T_coolant_max=0.0_8
        this%T_core_max=0.0_8
        
    end subroutine outputs_init
    
    
!==============================================================================
! LSSS_init
!==============================================================================   
    subroutine LSSS_init(this)
        ! declare arguments
        type(LSSS_type)     :: this
        
        this%INmaxcoolPOWER=0.0_8
        this%INmaxcoolTin=0.0_8
        this%INmaxcoolTout=0.0_8
        this%INmaxcoolTmax=0.0_8
        this%INmaxcoolToutavg=0.0_8
        this%INmaxfuelPOWER=0.0_8
        this%INmaxfuelTin=0.0_8
        this%INmaxfuelTout=0.0_8
        this%INmaxfuelTmax=0.0_8
        this%INmaxfuelToutavg=0.0_8
        this%INmaxoutPOWER=0.0_8
        this%INmaxoutTin=0.0_8
        this%INmaxoutTout=0.0_8
        this%OUTmaxcoolPOWER=0.0_8
        this%OUTmaxcoolTin=0.0_8
        this%OUTmaxcoolTout=0.0_8
        this%OUTmaxcoolTmax=0.0_8
        this%OUTmaxcoolToutavg=0.0_8
        this%OUTmaxfuelPOWER=0.0_8
        this%OUTmaxfuelTin=0.0_8
        this%OUTmaxfuelTout=0.0_8
        this%OUTmaxfuelTmax=0.0_8
        this%OUTmaxfuelToutavg=0.0_8
        this%minPOWER=0.0_8
        
    end subroutine LSSS_init
    

!==============================================================================
! print_LSSS
!==============================================================================
    subroutine print_LSSS(LSSS,Q_core,limits)
        ! declare arguments
        type(LSSS_type)     :: LSSS
        real(8)             :: Q_core
        type(limits_type)   :: limits
        
        write(*,*)
        write(*,'(A)')  "======================================================================"
        write(*,*)      "LSSS Calculation Results"
        write(*,'(A)')  "======================================================================"
        
        write(*,*)
        write(*,'(A,F6.3,A,F6.1,A)') "For primary mass flow rate of ", Q_core, " and  inlet temperature of ", limits%T_in_limit, ":"
        write(*,'(A)') "----------------------------------------------------------------------"
        
        ! #3 Min inlet channel temperature and maximum outlet average channel temperature limits
        write(*,*)
        write(*,'(A,F7.2,A)') "The maximum outlet average channel temperature exceeds ", limits%T_out_limit, " at:"
        write(*,'(F5.2,A,F7.2,A,F7.2)') LSSS%INmaxoutPOWER, " MW    T_in = ", LSSS%INmaxoutTin, "   T_out = ", LSSS%INmaxoutTout
        
        ! #1 Min inlet channel temperature and maximum coolant hot channel temperature limits
        write(*,*) 
        write(*,'(A,F7.2,A)') "The maximum  coolant hot channel   temperature exceeds ", limits%T_coolant_limit, " at:"
        write(*,'(F5.2,A,F7.2,A,F7.2,A,F7.2)') LSSS%INmaxcoolPOWER, " MW    T_in = ", LSSS%INmaxcoolTin, &
            &   "   T_out = ", LSSS%INmaxcoolTout, "   TCMAX = ", LSSS%INmaxcoolTmax
        write(*,'(A,F7.2)') "T_out_avg = ", LSSS%INmaxcoolToutavg
        
        ! #2 Min inlet channel temperature and maximum fuel hot channel temperature limits
        write(*,*) 
        write(*,'(A,F7.2,A)') "The maximum    fuel hot channel    temperature exceeds ", limits%T_fuel_limit, " at:"
        write(*,'(F5.2,A,F7.2,A,F7.2,A,F7.2)') LSSS%INmaxfuelPOWER, " MW    T_in = ", LSSS%INmaxfuelTin, &
            &   "   T_out = ", LSSS%INmaxfuelTout, "   TFMAX = ", LSSS%INmaxfuelTmax
        write(*,'(A,F7.2)') "T_out_avg = ", LSSS%INmaxfuelToutavg
        write(*,*) 
        
        write(*,*)
        write(*,'(A,F6.3,A,F6.1,A)') "For primary mass flow rate of ", Q_core, " and outlet temperature of ", limits%T_out_limit, ":"
        write(*,'(A)') "----------------------------------------------------------------------"
        
        ! #4 Max outlet average channel temperature and maximum coolant hot channel temperature limits
        write(*,'(A,F7.2,A)') "The maximum   coolant hot channel  temperature exceeds ", limits%T_coolant_limit, " at:"
        write(*,'(F5.2,A,F7.2,A,F7.2,A,F7.2)') LSSS%OUTmaxcoolPOWER, " MW    T_in = ", LSSS%OUTmaxcoolTin, &
            &   "   T_out = ", LSSS%OUTmaxcoolToutavg, "   TFMAX = ", LSSS%OUTmaxcoolTmax
        
        ! #5 Max outlet average channel temperature and maximum fuel hot channel temperature limits    
        write(*,*)
        write(*,'(A,F7.2,A)') "The maximum    fuel hot channel    temperature exceeds ", limits%T_fuel_limit, " at:"
        write(*,'(F5.2,A,F7.2,A,F7.2,A,F7.2)') LSSS%OUTmaxfuelPOWER, " MW    T_in = ", LSSS%OUTmaxfuelTin, &
            &   "   T_out = ", LSSS%OUTmaxfuelToutavg, "   TFMAX = ", LSSS%OUTmaxfuelTmax
        
        write(*,*)
        write(*,*)
        write(*,*) "The minimum LSSS power is:"
        write(*,'(A)') "----------------------------------------------------------------------"
        write(*,'(F7.2)') LSSS%minPOWER
        
        write(*,*)
        write(*,'(A)') "======================================================================"
        write(*,*)
        
    end subroutine print_LSSS
    
    
    
!===============================================================================
end module io
