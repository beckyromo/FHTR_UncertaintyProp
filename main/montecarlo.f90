module montecarlo
  implicit none
  private
  public :: nd_init_off,mcnd_study
  
  ! declare variables
  type,public :: nd_type
      real(8) :: cp   ! sensitivity parameter of heat capacity of flibe
      real(8) :: k    ! sensitivity parameter of thermal conductivity of flibe
      real(8) :: rho  ! sensitivity parameter of density of flibe
      real(8) :: mu   ! sensitivity parameter of viscosity of flibe
  end type nd_type
  
contains

!==============================================================================
! nd_init_off
!==============================================================================
    subroutine nd_init_off(this)
        ! declare arguments
        type(nd_type) :: this
    
        this%cp=0_8
        this%k=0_8
        this%rho=0_8
        this%mu=0_8

    end subroutine nd_init_off   
  
!==============================================================================
! mcnd_study
!==============================================================================
    subroutine mcnd_study(nhist,nd,inputoutput,limits,LSSS)
    
        use io, only: inputoutput_type, limits_type, LSSS_type, LSSS_init, inputoutput_init_Tin, inputoutput_init_outputs
        use sensitivity, only: sensitivity_init
        use prismatic, only: prismaticLSSS
        
        ! declare arguments
        integer                 :: nhist        ! number of histories
        integer                 :: i            ! counter
        type(nd_type)           :: nd
        type(inputoutput_type)  :: inputoutput
        type(limits_type)       :: limits
        type(LSSS_type)         :: LSSS
        
        write(*,*)
        write(*,'(A)')  "======================================================================"
        write(*,*)      "MC Sampling of ND Study Begin"
        write(*,'(A)')  "======================================================================"
                
        

        open(50,FILE="histories.txt")
            
        do i=1,nhist
            
            write(*,*) i

            ! Reset Outputs
            call LSSS_init(LSSS)
            call inputoutput_init_Tin(inputoutput)
            ! Run LSSS
            call prismaticLSSS(inputoutput,limits,LSSS)
            ! Print LSSS
            call print_histories(LSSS)
            
        end do
        
        close(50)
       

        !    write(*,'(5(F8.2),/5(F8.2),/2(F8.2),F24.2,/5(F8.2),/5(F8.2),1(F24.2))') LSSS
            
        
        write(*,*)
        write(*,'(A)')  "======================================================================"
        write(*,*)      "MC Sampling of ND Study Complete"
        write(*,'(A)')  "======================================================================"
        write(*,*)
        
    end subroutine mcnd_study
!===============================================================================
   
 
    
!==============================================================================
! print_histories
!==============================================================================    
    subroutine print_histories(LSSS)
        
        use io, only: LSSS_type
        
        ! declare arguments
        type(LSSS_type)         :: LSSS
        
            !!write(50,'(F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2, &
            !!& F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2)') &
            !!& LSSS%minPOWER, LSSS%INmaxoutPOWER, LSSS%INmaxoutTin, LSSS%INmaxoutTout, &
            !!& LSSS%INmaxcoolPOWER, LSSS%INmaxcoolTin, LSSS%INmaxcoolTout, LSSS%INmaxcoolTmax, LSSS%INmaxcoolToutavg, &
            !!& LSSS%INmaxfuelPOWER, LSSS%INmaxfuelTin, LSSS%INmaxfuelTout, LSSS%INmaxfuelTmax, LSSS%INmaxfuelToutavg, &
            !!& LSSS%OUTmaxcoolPOWER, LSSS%OUTmaxcoolTin, LSSS%OUTmaxcoolToutavg, LSSS%OUTmaxcoolTmax, &
            !!& LSSS%OUTmaxfuelPOWER, LSSS%OUTmaxfuelTin, LSSS%OUTmaxfuelToutavg, LSSS%OUTmaxfuelTmax
    end subroutine print_histories
!===============================================================================    
    
end module montecarlo