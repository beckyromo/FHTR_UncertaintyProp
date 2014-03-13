module sensitivity
  implicit none
  private
  public :: sensitivity_init,sensitivity_study,MCsensitivity_study
  
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
    subroutine sensitivity_study(s_min,s_step,s_max,this,sens,inputoutput,limits,LSSS)
    
        use io, only: inputoutput_type, limits_type, LSSS_type, LSSS_init, inputoutput_init_Tin, inputoutput_init_outputs
        use prismatic, only: prismaticLSSS
        use flibeprop, only: flibe_cp, flibe_mu, flibe_rho, flibe_k
        
        ! declare arguments
        real(8)                 :: s_min
        real(8)                 :: s_step
        real(8)                 :: s_max
        integer                 :: this     ! variable selector: 1=all, 2=cp, 3=k, 4=rho, 5=mu
        type(sensitivity_type)  :: sens
        real(8)                 :: s
        type(inputoutput_type)  :: inputoutput
        type(limits_type)       :: limits
        type(LSSS_type)         :: LSSS
        
        write(*,*)
        write(*,'(A)')  "======================================================================"
        write(*,*)      "Sensitivity Study Begin"
        write(*,'(A)')  "======================================================================"
                
        
        s=s_min
        
        if (this==1) then
            write(*,*) "Looping through all parameters:"
        end if
        
        do while (s<=s_max)
            
            if (this==2) then
                sens%s_cp=s               
            else if (this==3) then
                sens%s_k=s
            else if (this==4) then
                sens%s_rho=s
            else if (this==5) then
                sens%s_mu=s
            end if

            ! Reset Outputs
            call LSSS_init(LSSS)
            call inputoutput_init_Tin(inputoutput)
            !call inputoutput_init_outputs(inputoutput)
            
            
            write(*,*)
            write(*,'(A,F7.2,A,F6.3,A,F8.2,A,E11.4,A,F6.1,A)') "cp=", flibe_cp(inputoutput%T_in), &
                & ", k=", flibe_k(inputoutput%T_in), ", rho_ref=", flibe_rho(inputoutput%T_in), &
                & ", and mu_ref=", flibe_mu(inputoutput%T_in), " at T_ref=", inputoutput%T_in, ":"
            

            ! Run LSSS
            call prismaticLSSS(inputoutput,limits,LSSS)
            
            s=s+s_step

            write(*,'(5(F8.2),/5(F8.2),/2(F8.2),F24.2,/5(F8.2),/5(F8.2),1(F24.2))') LSSS
            
        end do
        
        write(*,*)
        write(*,'(A)')  "======================================================================"
        write(*,*)      "Sensitivity Study Complete"
        write(*,'(A)')  "======================================================================"
        write(*,*)
        
    end subroutine sensitivity_study
    
!===============================================================================
    
    
!==============================================================================
! MCsensitivity_study
!==============================================================================
    subroutine MCsensitivity_study(this,nhist,sens,inputoutput,limits,LSSS)
    
        use io, only: inputoutput_type, limits_type, LSSS_type, LSSS_init, inputoutput_init_Tin, inputoutput_init_outputs
        use random_numbers, only: rn_normal
        use prismatic, only: prismaticLSSS
        use flibeprop, only: flibe_cp, flibe_mu, flibe_rho, flibe_k
        
        ! declare arguments
        integer                 :: this     ! variable selector: 1=all, 2=cp, 3=k, 4=rho, 5=mu
        type(sensitivity_type)  :: sens
        integer                 :: i        ! count
        integer                 :: nhist    ! number of histories
        type(inputoutput_type)  :: inputoutput
        type(limits_type)       :: limits
        type(LSSS_type)         :: LSSS
        
        write(*,*)
        write(*,'(A)')  "======================================================================"
        write(*,*)      "Monte Carlo Sensitivity Study Begin"
        write(*,'(A)')  "======================================================================"
                
        open(50,FILE="histories.txt")
        
        if (this==1) then
            write(*,*) "All parameters:"
        end if
        
        do i=1,nhist
            
            if (this==1) then
                sens%s_cp=rn_normal(1.0d0,0.03d0)
                sens%s_k=rn_normal(1.0d0,0.1d0)
                sens%s_rho=rn_normal(1.0d0,0.02d0)
                sens%s_mu=rn_normal(1.0d0,0.2d0)
            else if (this==2) then
                sens%s_cp=rn_normal(1.0d0,0.03d0)
            else if (this==3) then
                sens%s_k=rn_normal(1.0d0,0.1d0)
            else if (this==4) then
                sens%s_rho=rn_normal(1.0d0,0.02d0)
            else if (this==5) then
                sens%s_mu=rn_normal(1.0d0,0.2d0)
            end if

            ! Reset Outputs
            call LSSS_init(LSSS)
            call inputoutput_init_Tin(inputoutput)
            !call inputoutput_init_outputs(inputoutput)
                       
            write(*,*)
            write(*,'(I5,4(F9.6))') nhist-i+1, sens          

            ! Run LSSS
            call prismaticLSSS(inputoutput,limits,LSSS)

            !!write(*,'(5(F8.2),/5(F8.2),/2(F8.2),F24.2,/5(F8.2),/5(F8.2),1(F24.2))') LSSS
            !!write(*,*)
            !!write(*,'(A)')  "----------------------------------------------------------------------"
            
            ! Print LSSS
            call print_histories(LSSS)
            
        end do
        
        close(50)
        
        write(*,*)
        write(*,'(A)')  "======================================================================"
        write(*,*)      "Monte Carlo Sensitivity Study Complete"
        write(*,'(A)')  "======================================================================"
        write(*,*)
        
    end subroutine MCsensitivity_study   
!===============================================================================
    
!==============================================================================
! print_histories
!==============================================================================    
    subroutine print_histories(LSSS)
        
        use io, only: LSSS_type
        
        ! declare arguments
        type(LSSS_type)         :: LSSS
        
            write(50,'(F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2, &
            & F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2)') &
            & LSSS%minPOWER, LSSS%INmaxoutPOWER, LSSS%INmaxoutTin, LSSS%INmaxoutTout, &
            & LSSS%INmaxcoolPOWER, LSSS%INmaxcoolTin, LSSS%INmaxcoolTout, LSSS%INmaxcoolTmax, LSSS%INmaxcoolToutavg, &
            & LSSS%INmaxfuelPOWER, LSSS%INmaxfuelTin, LSSS%INmaxfuelTout, LSSS%INmaxfuelTmax, LSSS%INmaxfuelToutavg, &
            & LSSS%OUTmaxcoolPOWER, LSSS%OUTmaxcoolTin, LSSS%OUTmaxcoolToutavg, LSSS%OUTmaxcoolTmax, &
            & LSSS%OUTmaxfuelPOWER, LSSS%OUTmaxfuelTin, LSSS%OUTmaxfuelToutavg, LSSS%OUTmaxfuelTmax
    end subroutine print_histories
!===============================================================================  
end module sensitivity
