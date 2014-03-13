!********************************************************************************
!
!  MODULE:  random_numbers
!
!  PURPOSE: Contains functions random number generation and normal distribution
!
!  FUNCTIONS:
!   rn          generates a random number
!   rnspread    generates a random number between a min and max value
!   rn_normal   generates a normal distribution with a mean and stddev (sigma)
!
!******************************************************************************** 
MODULE random_numbers
    
    CONTAINS
    
    ! FUNCTION rn
    function rn()
        implicit none
        double precision :: rn, x
        call random_number(x)
        rn=x
    end function rn
    
    ! FUNCTION rnspread
    function rnspread(min,max)
        implicit none
        double precision :: rnspread, min, max
        rnspread = (max-min) * rn() + min
    end function rnspread
    
    ! FUNCTION rn_normal
    function rn_normal(mean,sigma)
        implicit none
        double precision :: rn_normal, mean, sigma
        double precision :: tmp, fac, gsave, rsq, r1, r2
        integer flag
        save flag, gsave
        data flag /0/
        if (flag == 0) then
            rsq=2.0d0
            do while (rsq >= 1.0d0 .OR. rsq == 0.0d0)
                r1=2.0d0*rn()-1.0d0
                r2=2.0d0*rn()-1.0d0
                rsq=r1*r1+r2*r2
            end do
            fac=sqrt(-2.0d0*log(rsq)/rsq)
            gsave=r1*fac
            tmp=r2*fac
            flag=1
        else
            tmp=gsave
            flag=0
        end if
        rn_normal=tmp*sigma+mean
    end function rn_normal
    
!!!!!!! http://sukhbinder.wordpress.com/fortran-random-number-generation/ !also has random number for t distribution
!!!!!!! Random Sample from normal (Gaussian) distribution
!!!!!!!
!!!!!!      FUNCTION rand_normal(mean,stdev) RESULT(c)
!!!!!!       DOUBLE PRECISION :: mean,stdev,c,temp(2)
!!!!!!      IF(stdev <= 0.0d0) THEN
!!!!!!
!!!!!!        WRITE(*,*) "Standard Deviation must be +ve"
!!!!!!      ELSE
!!!!!!        CALL RANDOM_NUMBER(temp) ! random_number is my rnspread
!!!!!!        r=(-2.0d0*log(temp(1)))**0.5
!!!!!!        theta = 2.0d0*PI*temp(2)
!!!!!!        c= mean+stdev*r*sin(theta)
!!!!!!      END IF
!!!!!!      END FUNCTION
    
END MODULE random_numbers