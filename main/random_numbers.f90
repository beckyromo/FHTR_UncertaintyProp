!********************************************************************************
!
!  MODULE:  random_numbers
!
!  PURPOSE: Contains functions random number generation and normal distribution
!
!  FUNCTIONS:
!   rn          generates a random number
!   rnspread    generates a random number between a min and max value
!   normal      generates a normal distribution with a mean and stddev (sigma)
!
!******************************************************************************** 
MODULE random_numbers
    
    CONTAINS
    
    function rn()
        implicit none
        double precision :: rn, x
        call random_number(x)
        rn=x
    end function rn
    
    function rnspread(min,max)
        implicit none
        double precision :: rnspread, min, max
        rnspread = (max-min) * rn() + min
    end function rnspread
    
    function normal(mean,sigma)
        implicit none
        double precision :: normal, mean, sigma
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
        normal=tmp*sigma+mean
    end function normal
    
END MODULE random_numbers