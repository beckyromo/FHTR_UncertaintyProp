!==============================================================================!
! MODULE: pdfs
!==============================================================================!
module pdfs

  use random_number_generator, only: rand
  use global, only: PI

  implicit none

contains

    function get_normdist_value(mean,stddev)
  
        implicit none
  
        double precision, intent(in)  :: mean
        double precision, intent(in)  :: stddev
        double precision              :: rn1
        double precision              :: rn2
        double precision              :: get_normdist_value
!      double precision,allocatable  :: CDF(:)
!      double precision,allocatable  :: PDF(:) 
 !     integer                       :: i
!      integer                       :: N
        double precision              :: x
        double precision              :: y
        double precision              :: z
        double precision              :: a
        double precision              :: b
        double precision              :: c
        
  
        call RANDOM_SEED()
        call RANDOM_NUMBER(rn1)
        call RANDOM_NUMBER(rn2)
        !write(*,*) rn1, rn2
        
          ! get rn
          !rn = rand()
          z=2.0
          do while (z>1.0)
          x=2.0*rn1-1
          y=2.0*rn2-1
          z=x**2.0+y**2.0
            call RANDOM_NUMBER(rn1)
            call RANDOM_NUMBER(rn2)
          end do
          if (z<1.0d0) then
              c=sqrt(-2.0*log(z)/z)
              a=x*c
              b=y*c
          end if
    
          get_normdist_value = b*stddev+mean
      
    end function get_normdist_value

  
end module pdfs

        !!!!!!!    !================================================================================
        !!!!!!!! Testing get_normdist_value 
        !!!!!!!!================================================================================                
        !!!!!!!open(30,FILE="norm.txt") 
        !!!!!!!do i=1,1000000
        !!!!!!!    write(30,'(F9.2)') get_normdist_value(1.0d0,0.5d0)
        !!!!!!!    !write(30,'(F9.2)') exp(-(i-1.0-0.0)**2.0/sqrt(2.0*1.0**2.0))/(1.0*sqrt(2.0*PI))+exp(-(i-0.0)**2.0/sqrt(2.0*1.0**2.0))/(1.0*sqrt(2.0*PI))*1.0/1000
        !!!!!!!end do
        !!!!!!!close(30)
        
        
      
  !!=============================================================================
  !!> @brief Get a particle position in the slab
  !!>
  !!> @param[in]     L         Total width of domain.
  !!> @return                  A random location in the domain
  !!=============================================================================
  !function get_particle_pos(L)
  !
  !  double precision              :: get_particle_pos
  !  double precision, intent(in)  :: L
  !
  !  get_particle_pos = rand()*L
  !
  !end function get_particle_pos
  
  
        !N=1000
      !allocate(CDF(N))
      !allocate(PDF(N))
      !
      !! Probability Distribution Function (Normal Distribution)
      !do i=1,N
      !  PDF=exp(-(i-mean)**2.0/sqrt(2.0*stddev**2.0))/(stddev*sqrt(2.0*PI))
      !end do 
      !  
      !! Cumulative Distribution Function (Normal Distribution)
      !CDF(1)=0.0
      !do i=2,N
      !    !CDF=0.5*(1+erf( (rn-mean)/(stddev*sqrt(2.0)) ) )
      !    CDF=CDF(i-1)+PDF(i)/N
      !end do
      !
      !write(*,*) CDF
      !
      !get_normdist_value=1.0
      !
      !deallocate(PDF)
      !deallocate(CDF)