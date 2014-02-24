module sensitivity
  implicit none
  private
  public :: sensitivity_init!, sensitivity_study
  
  ! declare variables
  type,public :: sensitivity_type
    real(8)   :: s_cp   ! sensitivity parameter of heat capacity of flibe
    real(8)   :: s_k    ! sensitivity parameter of thermal conductivity of flibe
    real(8)   :: s_rho  ! sensitivity parameter of density of flibe
    real(8)   :: s_mu   ! sensitivity parameter of viscosity of flibe
  end type sensitivity_type
  
contains

!==============================================================================
! sensitivity_init
!==============================================================================
    subroutine sensitivity_init(this)
        ! declare arguments
        type(sensitivity_type) :: this
    
        this%s_cp=1.0_8
        this%s_k=1.0_8
        this%s_rho=1.0_8
        this%s_mu=1.0_8

    end subroutine sensitivity_init
  
!==============================================================================
! sensitivity_study
!==============================================================================
    !subroutine sensitivity_study(s_min,s_step,s_max,this,sens,inputoutput,limits,LSSS)
    !
    !    use io, only: inputoutput_type, limits_type, LSSS_type, LSSS_init
    !    use prismatic!, only: prismaticLSSSloop
    !    use flibeprop, only: flibe_cp, flibe_mu, flibe_rho, flibe_k
        
        !! declare arguments
        !real(8)                 :: s_min
        !real(8)                 :: s_step
        !real(8)                 :: s_max
        !integer                 :: this     ! variable selector: 1=all, 2=cp, 3=k, 4=rho, 5=mu
        !type(sensitivity_type)  :: sens
        !real(8)                 :: s
        !integer                 :: nsteps
        !type(inputoutput_type)  :: inputoutput
        !type(limits_type)       :: limits
        !type(LSSS_type)         :: LSSS
        
        !!write(*,*)
        !!write(*,'(A)')  "======================================================================"
        !!write(*,*)      "Sensitivity Study Begin"
        !!write(*,'(A)')  "======================================================================"
        !!
        !!
        !!
        !!s=s_min
        !!
        !!if (this==1) then
        !!    write(*,*) "Looping through all parameters:"
        !!end if
        !!
        !!
        !!do while (s<=s_max)
        !!
        !!    
        !!    if (this==2) then
        !!        sens%s_cp=s
        !!        
        !!    else if (this==3) then
        !!        sens%s_k=s
        !!    else if (this==4) then
        !!        sens%s_rho=s
        !!    else if (this==5) then
        !!        sens%s_mu=s
        !!        write(*,*) flibe_mu(700.0_8)
        !!    end if
        !!
        !!    call LSSS_init(LSSS)
        !!    call prismaticLSSSloop(input,limits,LSSS)
        !!    
        !!    s=s+s_step
        !!    
        !!    write(*,'(A,F7.2,A,F6.3,A,F8.2,A,E11.4,A,F7.2,A)') "For cp=", flibe_cp(input%T_in), ", k=", flibe_k(input%T_in), &
        !!        & ", rho=", flibe_rho(input%T_in), ", and mu=", flibe_mu(input%T_in), " at Tin=", input%T_in, ":"
        !!    write(*,*) LSSS
        !!    
        !!end do
        !!
        !!write(*,*)
        !!write(*,'(A)')  "======================================================================"
        !!write(*,*)      "Sensitivity Study Complete"
        !!write(*,'(A)')  "======================================================================"
        !!write(*,*)
        
    !end subroutine sensitivity_study
    
!===============================================================================
end module sensitivity
