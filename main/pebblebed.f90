!********************************************************************************
!
!  MODULE:  pebblebed
!
!  PURPOSE: Contains functions for pebble bed core calculations and subroutines 
!           to do LSSS calculations for a pebble bed core
!
!  FUNCTIONS:
!  porosity_PB
!  HTC_PB
!  FF_PB
!
!  SUBROUTINES:
!  pebblebedcore        calculates core temperatures for given power, mass flow 
!                       rate, and inlet temperature for either average or hot 
!                       channel, prints core temperature results and outputs the 
!                       outlet temperature, maximum coolant temperature, and 
!                       maximum core/fuel temperature
!  pebblebedLSSSloop    calls subroutine pebblebedcore iterating on power to 
!                       print LSSS for given fuel, coolant, and outlet 
!                       temperature limits
!
!******************************************************************************** 
MODULE pebblebed
            
    USE global
    USE flibeprop, ONLY: flibe_cp, flibe_k, flibe_rho, flibe_mu
        
    IMPLICIT NONE


CONTAINS
    
    
    !================================================================================
    !  SUBROUTINE: pebblebedcore
    !================================================================================
    ! INPUTS
    !   POWER           :: Power [W]
    !   Q_core          :: Volumetric flow rate [m^3/s]
    !   T_in            :: Inlet temperature [Celcius]
    !   channel         :: Selects for average (1) or hot (2) channel calculation
    ! OUTPUTS    
    !   T_out           :: Outlet temperature [Celcius]
    !   T_coolant_max   :: Maximum coolant temperature [Celcius]
    !   T_core_max      :: Maximum core/fuel temperature [Celcius]
    ! REFERENCE
    !   Yao's Code, which uses power distribution data
    !================================================================================
    SUBROUTINE pebblebedcore(POWER,Q_core,T_in,T_out,T_coolant_max,T_core_max,channel)
    
        USE global
        USE flibeprop, ONLY: flibe_cp, flibe_k, flibe_rho, flibe_mu, flibe_enthalpy, flibe_temperature
        USE trisoprop, ONLY: k_TRISO, k_SiC, K_densePyC, k_PyC, k_UO2, k_TRISOlayer
        USE pipes, ONLY: FF_PIPE, dP_PIPE
        
        IMPLICIT NONE

        !================================================================================
        ! Declare Variables
        !--> Declare/Initialize, Declare arrays, initialize array inputs
        !================================================================================
        
        integer,intent(in)  :: channel          ! 1 for average, 2 for hot channel
        
        integer             :: I                ! Loop counter
        integer             :: m                ! Loop counter
        integer             :: N_core           ! Number of nodes(CVs) the core is split into
        ! INPUTS
        real(8)             :: AFP(21)
        real(8)             :: LAMBDA(6)
        real(8)             :: BETA(6)
        ! UO2,PyC,DPyC,SiC,DPyC,TRISO+Graphite,Graphite
        real(8)             :: RPF(7)           ! Fuel radii [m]
        real(8)             :: KPF(7)           ! Fuel thermal conductivity [W/m-C]
        real(8)             :: RHOPF(7)         ! Fuel density [kg/m^3]
        real(8)             :: CPPF(7)          ! Fuel heat capacity [J/kg-C]
        real(8)             :: D_pebble         ! Diameter of pebble [m]                     
        real(8)             :: DD1              ! Face to face distance 1 of hexagonal core
        real(8)             :: DD2              ! Face to face distance 2 of hexagonal core
        real(8)             :: H_core           ! Height of the core [m]
        real(8)             :: FH               ! HCF - enthalpy rise hot channel factor (HCF)
        real(8)             :: FDTW             ! HCF - film/wall temperature rise
        real(8)             :: FDTF             ! HCF - fuel temperature rise
        real(8)             :: FCORE            ! HCF - 
        real(8)             :: FFUEL            ! HCF - 
        real(8)             :: FFDF             ! HCF - channel flow disparity factor
        real(8)             :: FKZ              ! Radial power peaking factor
        real(8)             :: FKTRISO          ! TRISO power peaking factor          
        real(8),intent(in)  :: POWER            ! Reactor power [W]
        real(8)             :: DECAY_HEAT       ! Decay heat of reactor [W]
        real(8)             :: POWER_NORM       ! Norminal Fisison Power [W]
        real(8),intent(in)  :: Q_core           ! Volumetric flow rate in the core [m^3/s]
        real(8)             :: Q                ! Volumetric flow rate in the core [m^3/s]
        real(8)             :: W                ! Mass flow rate in the core [kg/s]
        real(8),intent(in)  :: T_in             ! Inlet temperature of the core [Celcius]
        
        real(8)             :: NTRISO           ! Number of TRISO particles per pebble
        real(8)             :: VolPF1           ! Volume of first layer of TRISO (UO2) [m^3]
        real(8)             :: VolPF5           ! Volume of first 5 layers of TRISO [m^3]
        real(8)             :: VolPF6           ! Volume of fuel matrix (TRISOs and graphite) of pebble [m^3]
        real(8)             :: VolPF7           ! Volume of pebble (fuel matrix plus graphite cladding) [m^3]
        real(8)             :: A_core           ! Area of the core [m^2]
        real(8)             :: D_core           ! Hydraulic diameter of the core [m]
        
        real(8)             :: FF_core          ! Friction factor in core 
        real(8)             :: HTC_core         ! Heat transfer coefficient in core [W/m^2-C]
        real(8)             :: A_HTC            ! Area of core with convection heat transfer [m^2]
        real(8)             :: ACTUAL_POWER     ! Actual Power [W]
        real(8)             :: dP_core          ! Core pressure drop [Pa] 
        real(8),allocatable :: enthalpy_core(:) ! Enthalpy in each node [J?/kg]
        real(8),allocatable :: T_w_core(:)      ! Pebble surface temperature at each node  [Celcius]
        real(8),allocatable :: T_CL_core(:)     ! Centerline temperature at each node  [Celcius]
        real(8),allocatable :: T_CL_TRISO(:)    ! Centerline temperature of centerline TRISO at each node [C]
        real(8),intent(out) :: T_coolant_max    ! Maximum coolant temperature [Celcius]
        real(8),intent(out) :: T_core_max       ! Maximum core temperature  [Celcius]
        real(8),intent(out) :: T_out            ! Outlet temperature of the core [Celcius]
        real(8)             :: TC               ! Temperature of coolant in core [Celcius]
        real(8)             :: TG               ! Temperature at inner edge of graphite cladding [Celcius]
        real(8)             :: TW               ! Temperature at the pebble surface [Celcius]
        real(8)             :: T_temp           ! Temperature place holder [Celcius]
        real(8)             :: NPPN             ! Number of pebbles in each node(CV)
        real(8)             :: PPP              ! Power per pebble [W]
        real(8)             :: PPT              ! Power per TRISO particle [W]
        real(8)             :: TOL              ! While loop error tolerance
        real(8)             :: U                ! Thermal resistance

        
        
        !================================================================================
        ! Initialize Variables / Inputs
        !================================================================================
        
        ! POWER DISTRIBUTION FROM SINAP REPORT
        DATA AFP/1.3522E-002,4.7941E-002,4.6349E-002,4.7030E-002,4.8748E-002,5.0615E-002,5.2291E-002,&
                 5.3749E-002,5.4762E-002,5.5298E-002,5.5540E-002,5.5100E-002,5.4286E-002,5.3050E-002,&
                 5.1321E-002,4.9228E-002,4.7067E-002,4.4800E-002,4.2871E-002,4.2380E-002,3.4045E-002/
        
        ! LAMBDA AND BETA FOR POINT KINETICS
        !DATA LAMBDA/1.271E-2,3.174E-2,1.16E-1,3.11E-1,1.4E0,3.87E0/
        !DATA BETA/2.3794E-4, 1.5831E-3,1.4183E-3, 2.8550E-3,8.3273E-4,3.0192E-4/
        
        ! ENGINERRING HOT CHANNEL FACTORS -------------------------------------
        FFDF=0.8_8
        FH=1.172_8      ! 1.173 is only MITR value, 1.172 also uses some HTC-10 pebblebed
        FDTW=1.268_8    ! 1.275 is only MITR value, 1.168 also uses some HTC-10 pebblebed
        FDTF=1.114_8    ! 1.123 is only MITR value, 1.172 also uses some HTC-10 pebblebed
        ! See ICONE paper of Yao's for MITR values, see Yao's report for HTC-10 values
        FCORE=0.965_8       ! NOT USED
        FFUEL=0.940_8       ! NOT USED
        FKZ=1.4_8
        FKTRISO=1.0_8       ! NOT USED

        ! CORE
        DECAY_HEAT=0.0_8  ! [W]
        POWER_NORM=1.0_8  ! [W]
        DD1=134.69E-2_8   ! Face to face distance 1 [m]
        DD2=139.0E-2_8    ! Face to face distance 2 [m]
        H_core=138.6_8    ! Height of the core [m]
        D_pebble=6.0E-2_8 ! Diameter of fuel pebble [m]
        N_core=21         ! Number of nodes(CVs) in the core

        
        ! FUEL: UO2,PyC,DPyC,SiC,DPyC,TRISO+Graphite,Graphite (
        DATA	RPF/0.25E-3,0.345E-3,0.385E-3,0.42E-3,0.46E-3,25.0E-3,30.0E-3/  ! Radius [m]
        DATA	KPF/2.5,0.5,4.0,13.9,4.0,30.0,30.0/                             ! Thermal conductivity [W/m-C] NOTE: OVERWRITTEN
        !DATA	RHOPF/10.5E3,1.05E3,1.9E3,3.18E3,1.9E3,1.73E3,1.73E3/           ! Density [kg/m^3]
        !DATA	CPPF/332,1.5,1.5,1.5,1.5,710,710/                               ! Heat capacity [J/kg-K]
        
        ! Calculated inputs
        VolPF5=4.0/3.0*PI*RPF(5)**3.0 ! TRISO volume
        VolPF6=4.0/3.0*PI*RPF(6)**3.0
        VolPF7=4.0/3.0*PI*RPF(7)**3.0
        VolPF1=4.0/3.0*PI*RPF(1)**3.0 ! Fuel kernel volume
        NTRISO=VolPF6*7.5/100.0/VolPF5
        A_core=DD1**2-((DD1*2.0**0.5-DD2)/2.0**0.5)**2*2.0-0.35*0.35
        D_core=(4.0*A_core/PI)**0.5                                     ! HYDRAULIC RADIUS OF THE CORE [m]
        ACTUAL_POWER=POWER_NORM*POWER+DECAY_HEAT       
        

        
        !================================================================================
        ! Set criteria for average (1) or hot (2) channel
        !--------------------------------------------------------------------------------
        !
        ! channel==1 --> average channel:
        ! Override hot channel factors to 1.0_8 such that the hot channel factors are 
        ! not applied for the average channel.
        !
        ! channel==2 --> hot channel:
        ! Adjust flow rate for hot channel using the channel flow disparity factor.
        !
        !================================================================================
        IF (channel==1) THEN
            W=Q_core*flibe_rho(T_in)
            FH=1.0_8
            FKZ=1.0_8
            FDTW=1.0_8
            FDTF=1.0_8
            FFDF=1.0_8
        ELSE IF (channel==2) THEN
            W=Q_core*flibe_rho(T_in)*FFDF
        END IF
            
                
       
        !================================================================================       
        ! Core Loop Calculations - Average channel
        !================================================================================         

        ! Initialize variables for loop
        dP_core=0.0
        T_coolant_max=0.0
        T_core_max=0.0        
        ! Allocate arrays
        allocate(enthalpy_core(N_core)) ! Array of enthalpies in each node(CV) of the core
        allocate(T_w_core(N_core))      ! Array of surface tempertures in each node(CV) of the core
        allocate(T_CL_core(N_core))     ! Array of central temperature in each node(CV) of the core
        allocate(T_CL_TRISO(N_core))    ! Array of central temperature in TRISO in each node(CV) of the core
        
        ! Loop through each core CV node
        DO I=1,N_core
            
            IF (I==1) THEN
                enthalpy_core(I)=flibe_enthalpy(T_in)+ACTUAL_POWER*AFP(I)*FH*FKZ/W
                TC=( flibe_temperature(enthalpy_core(I)) + T_in ) / 2.0       ! Average coolant temperature in node
            ELSE
                enthalpy_core(I)=enthalpy_core(I-1)+ACTUAL_POWER*AFP(I)*FH*FKZ/W
                TC=( flibe_temperature(enthalpy_core(I)) + flibe_temperature(enthalpy_core(I-1)) ) / 2.0
            END IF   
            
            ! Calculate pressure drop in node
            IF (channel==1) THEN
                FF_core=FF_PB(W,D_core,D_pebble,TC,1)
                dP_core=dP_core+0.5*FF_core*(H_core/N_core)*(W/A_core)**2.0/(flibe_rho(TC)*D_pebble) &
                        + flibe_rho(TC)*GRAVITY*H_core/N_core
            END IF
            
            ! Calculate heat transfer coefficient in node
            HTC_core=HTC_PB(W,D_core,D_pebble,TC,1)
            
            ! Find wall temperature
            IF (I==1) THEN
                NPPN=140.0                      ! Number of pebbles
            ELSE IF (I==21) THEN
                NPPN=416.0                      ! Number of pebbles
            ELSE
                NPPN=556.0                      ! Number of pebbles
            END  IF
            A_HTC=NPPN*PI*D_pebble**2.0
            TW=ACTUAL_POWER*AFP(I)*FKZ/(HTC_core*A_HTC) + TC
            TW=(TW-TC)*FDTW+TC
            T_w_core(I)=TW
            
            ! Check if temperature is maximum
            IF (TW>T_coolant_max) THEN
                T_coolant_max=TW
            END IF
            
            ! Calculate power in pebble and centerline TRISO particle
            PPP=ACTUAL_POWER*AFP(I)*FKZ/NPPN        ! Power per pebble
            PPT=PPP/NTRISO                      ! Power per TRISO particle

            ! Calculate thermal conductivties and find temperatures
            TOL=1.0e-5           
                        
            ! GRAPHITE CLADDING
            T_temp=TW
            TG=0.0
            do while (abs(T_temp-TG)>TOL)
                T_temp=TG
                KPF(7)=k_TRISO( (TW+T_temp)/2.0 )
                U=1.0/(4.0*PI*KPF(7)) * (1.0/RPF(6) - 1.0/RPF(7))
                TG=PPP*U + TW
            end do
            
            ! FUEL MATRIX (TRISOs IN GRAPHITE MATRIX)
            T_temp=TG
            T_CL_core(I)=0.0
            do while (abs(T_temp-T_CL_core(I))>TOL)
                T_temp=T_CL_core(I)
                KPF(6)=k_TRISO( (TG+T_temp)/2.0 )
                T_CL_core(I)=PPP*(RPF(6)**2.0/6.0/KPF(6)/VolPF6) + TG
            end do
            
            ! CENTERLINE TRISO PARTICLE
            T_temp=T_CL_core(I)
            U=0.0
            do m=2,5
                U=1.0/(4.0*PI*k_trisolayer(m,T_temp)) * (1.0/RPF(m-1) - 1.0/RPF(m)) + U
            end do
            T_CL_TRISO(I)=PPT*(RPF(1)**2.0/6.0/k_trisolayer(1,T_temp)/VolPF1+U) + T_CL_core(I)
            
            ! Apply hot channel factors
            T_CL_TRISO(I)=(T_CL_TRISO(I)-TW)*FDTF+TW
            T_CL_core(I)=(T_CL_core(I)-TW)*FDTF+TW
            
            ! Check if temperature is max
            if (T_CL_core(I)>T_core_max) then
                T_core_max=T_CL_core(I)
            end if
            
            ! Print results to file output.txt
            if (I==1) then
                write(10,*)
                if (channel==1) then
                    write(10,*) "      TEMPERATURES IN AVERAGE CHANNEL       "
                else if (channel==2) then
                    write(10,*) "        TEMPERATURES IN HOT CHANNEL         "
                end if
                write(10,*) "       W       TC       TW       TG      TCL     TCLT"
            end if
            write(10,'(F9.3,F9.3,F9.3,F9.3,F9.3,F9.3)') W, TC, TW, TG, T_CL_core(I), T_CL_TRISO(I)
            ! write to command line        
            !write(*,*) PPP, PPT, TC, TW, TG, T_CL_core(I), T_CL_TRISO(I)
            
        
        end do
        
        
        T_out=flibe_temperature(enthalpy_core(N_core))    
        
        
        !================================================================================
        ! Print LSSS
        !================================================================================  
        !write(*,*)
        !write(*,*) "                           LSSS Results                           "
        !write(*,*) "   POWER      W      T_IN    T_OUT    TC_MAX   TF_MAX"
        !write(*,'(F9.2,F9.3,F9.3,F9.3,F9.3,F9.3)') &
        !            POWER/1.0E6,W,T_in,T_out,T_coolant_max,T_core_max
        !write(*,*)
        
        write(10,*)
        write(10,*) "                           LSSS Results                           "
        write(10,*) "   POWER      Q      T_IN    T_OUT    TC_MAX   TF_MAX"
        write(10,'(F9.2,F9.2,F9.3,F9.3,F9.3,F9.3)') &
                    POWER/1.0E6,Q_core*3600.0,T_in,T_out,T_coolant_max,T_core_max
        write(10,*)
        write(10,*) "=================================================================="
        write(10,*)
                
        
        !================================================================================
        ! Deallocate all arrays
        !================================================================================
        deallocate(enthalpy_core)  
        deallocate(T_w_core)       
        deallocate(T_CL_core)      
        deallocate(T_CL_TRISO)
        
    END SUBROUTINE pebblebedcore
    
    
    
    
    
    
    !================================================================================
    !  SUBROUTINE: pebblebedLSSSloop  [W/m-K]
    !================================================================================
    ! INPUTS    
    ! OUTPUTS   
    ! REFERENCE 
    !================================================================================
    SUBROUTINE pebblebedLSSSloop(Q_core,T_in,T_fuel_limit,T_coolant_limit,T_out_limit,T_out,T_coolant_max,T_core_max)
        
        IMPLICIT NONE
        
        integer             :: I                    ! Loop counter
        logical             :: fuelflag             ! Flag if max fuel temperature has been met in hot channel
        logical             :: coolantflag          ! Flag if max coolant temperature has been met in hot channel
        logical             :: Toutflag             ! Flag if max coolant outlet temperature has been met in average channel
        real(8)             :: T_fuel_limit         ! LSSS temoerature limit for fuel (max fuel temp in hot channel) 
        real(8)             :: T_coolant_limit      ! LSSS temperature limit for coolant (max coolant temp in hot channel)
        real(8)             :: T_out_limit          ! LSSS temperature limit for T_out (max T_out temp in average channel)
        real(8)             :: TC                   ! Temperature [Celcius]
        real(8)             :: POWER                ! Reactor power [W]
        real(8)             :: DECAY_HEAT           ! Decay heat of reactor [W]
        real(8)             :: POWER_NORM           ! Norminal Fisison Power [W]
        real(8)             :: Q_core               ! Volumetric flow rate in the core [kg/s]
        real(8)             :: T_in                 ! Inlet temperature of the core [Celcius]
        real(8)             :: T_out                ! Maximum coolant temperature [Celcius]
        real(8)             :: T_out_avg            ! Holder for T_out in average channel
        real(8)             :: T_coolant_max        ! Maximum core temperature  [Celcius]
        real(8)             :: T_core_max           ! Outlet temperature of the core [Celcius]
        real(8)             :: step                 ! Step size for power do loop
        real(8)             :: T_in_coolantrange    
        real(8)             :: POWER_coolantrange   
        real(8)             :: T_in_fuelrange       
        real(8)             :: POWER_fuelrange

            
        write(20,*) "                                               LSSS Results"
        write(20,*)
        write(20,'(A54)',ADVANCE='NO') "                           AVG CHANNEL                             "    
        write(20,'(A54)')              "                           HOT CHANNEL                             "
        write(20,'(A54)',ADVANCE='NO') " =================================================================="    
        write(20,'(A54)')              " =================================================================="
        write(20,'(A54)',ADVANCE='NO') "   POWER      Q      T_IN    T_OUT    TC_MAX   TF_MAX"    
        write(20,'(A54)')              "   POWER      Q      T_IN    T_OUT    TC_MAX   TF_MAX"
        write(20,'(A54)',ADVANCE='NO') " ------------------------------------------------------------------"    
        write(20,'(A54)')              " ------------------------------------------------------------------"
            
        !================================================================================
        ! Find power level for LSSS limits for given T_in and mass flow rate 
        !================================================================================
        ! Reset flags to .false.
        fuelflag=.false. 
        coolantflag=.false.
        Toutflag=.false.
        write(*,*)
        write(*,'(A,F7.2,A,F7.2,A)') "For primary volumetric flow rate of ", Q_core*3600, " m^3/hr and  inlet temperature of ", T_in, ":"
        write(*,'(A)') "======================================================================"
        write(*,*)
        do POWER=0.0E6,60.0E6,0.01E6
            !================================================================================
            ! Calculate core temperatures for average channel
            call pebblebedcore(POWER,Q_core,T_in,T_out,T_coolant_max,T_core_max,1)
            write(20,'(F9.2,F9.3,F9.3,F9.3,F9.3,F9.3)', ADVANCE='NO') &
                        POWER/1.0E6,Q_core*3600.0,T_in,T_out,T_coolant_max,T_core_max
            T_out_avg=T_out
            ! Check for exceedance of LSSS T_out temperature limit
            if (T_out>=T_out_limit .AND. Toutflag==.false.) then
                write(*,'(A,F7.2,A)') "The maximum outlet average channel temperature exceeds ", T_out_limit, " at:"
                write(*,'(F5.2,A,F7.2,A,F7.2)') POWER/1.0E6, " MW    T_in = ", T_in, "   T_out = ", T_out
                write(*,*) 
                Toutflag=.true.
            end if
            !================================================================================
            ! Calculate core temperatures for hot channel
            call pebblebedcore(POWER,Q_core,T_in,T_out,T_coolant_max,T_core_max,2)
            write(20,'(F9.2,F9.3,F9.3,F9.3,F9.3,F9.3)') &
                        POWER/1.0E6,Q_core*3600,T_in,T_out,T_coolant_max,T_core_max
            ! Check for exceedance of LSSS coolant temperature limit   
            if (T_coolant_max>=T_coolant_limit .AND. coolantflag==.false.) then
                write(*,'(A,F7.2,A)') "The maximum  coolant hot channel   temperature exceeds ", T_coolant_limit, " at:"
                write(*,'(F5.2,A,F7.2,A,F7.2,A,F7.2)') POWER/1.0E6, " MW    T_in = ", T_in, "   T_out = ", T_out, "   TCMAX = ", T_coolant_max
                write(*,'(A,F7.2)') "T_out_avg = ", T_out_avg
                write(*,*) 
                coolantflag=.true.
            end if
			! Check for exceedance of LSSS fuel temperature limit
            if (T_core_max>=T_fuel_limit .AND. fuelflag==.false.) then
                write(*,'(A,F7.2,A)') "The maximum    fuel hot channel    temperature exceeds ", T_fuel_limit, " at:"
                write(*,'(F5.2,A,F7.2,A,F7.2,A,F7.2)') POWER/1.0E6, " MW    T_in = ", T_in, "   T_out = ", T_out, "   TFMAX = ", T_core_max
                write(*,'(A,F7.2)') "T_out_avg = ", T_out_avg
                write(*,*) 
                fuelflag=.true.
            end if
            
        end do

        
        
        !================================================================================
        ! Find power level for LSSS limits for given T_out and mass flow rate :: Find T_in and POWER range
        !================================================================================
        ! Reset flags to .false.
        fuelflag=.false. 
        coolantflag=.false.
        Toutflag=.false.
        write(*,*)
        write(*,'(A,F7.2,A,F7.2,A)') "For primary volumtric flow rate of ", Q_core*3600, " m^3/hr and outlet temperature of ", T_out_limit, ":"
        write(*,'(A)') "======================================================================"
        do T_in=400.0_8,T_out_limit,1.0_8 
            ! Check if all LSSS limits met
            if (Toutflag==.true. .AND. coolantflag==.true. .AND. fuelflag==.true.) then
                write(*,*) "LSSS Limits Range Found"
                write(*,*)
                EXIT
            else          
                ! Loop through power level to find T_out=T_out_limit
                ! Reset T_out flags to .false.
                Toutflag=.false.
                do POWER=0.0E6,50.0E6,1.0E6
                    if (Toutflag==.false.) then
                        ! Calculate core temperatures for average channel
                        call pebblebedcore(POWER,Q_core,T_in,T_out,T_coolant_max,T_core_max,1)
                        T_out_avg=T_out
                        ! Check for exceedance of LSSS T_out temperature limit
                        if (T_out>=T_out_limit .AND. Toutflag==.false.) then
                            Toutflag=.true.
                            ! Calculate core temperatures for hot channel
                            call pebblebedcore(POWER,Q_core,T_in,T_out,T_coolant_max,T_core_max,2)
                            if (T_coolant_max <= T_coolant_limit .AND. coolantflag==.false.) then
                                !write(*,*) T_in, T_out_avg
                                T_in_coolantrange=T_in-5.0
                                POWER_coolantrange=POWER
                                !write(*,'(A,F7.2,A,F5.2,A)') & 
                                !"For coolant LSSS start T_in at ", T_in_coolantrange, " and the power level at ", POWER/1.0E6, " MW" 
                                coolantflag=.true.
                            end if
                            if (T_core_max <= T_fuel_limit .AND. fuelflag==.false.) then
                                !write(*,*) T_in,  T_out_avg
                                T_in_fuelrange=T_in-7.0
                                POWER_fuelrange=POWER
                                !write(*,'(A,F7.2,A,F5.2,A)') & 
                                !"For fuel LSSS start T_in at ", T_in_fuelrange, " and the power level at ", POWER/1.0E6, " MW" 
                                !write(*,*) T_core_max
                                fuelflag=.true.
                            end if
                        end if
                    end if
                end do
            end if
        end do
        !================================================================================
        ! Find power level for LSSS limits for given T_out and mass flow rate :: use range found for T_in and POWER
        !================================================================================
        ! Reset flags to .false.
        fuelflag=.false. 
        coolantflag=.false.
        Toutflag=.false.
        !write(*,*)
        !write(*,'(A,F5.2,A,F7.2,A)') "For primary mass flow rate of ", W_core, " and outlet temperature of ", T_out_limit, ":"
        do T_in=T_in_coolantrange,T_out_limit,0.01_8 
            ! Check if all LSSS limits met
            if (Toutflag==.true. .AND. coolantflag==.true.) then
                !write(*,*) "LSSS Coolant Limit Found"
                EXIT
            else          
                ! Loop through power level to find T_out=T_out_limit
                ! Reset T_out flags to .false.
                Toutflag=.false.
                do POWER=POWER_coolantrange,50.0E6,0.01E6_8
                    if (Toutflag==.false.) then
                        ! Calculate core temperatures for average channel
                        call pebblebedcore(POWER,Q_core,T_in,T_out,T_coolant_max,T_core_max,1)
                        T_out_avg=T_out
                        ! Check for exceedance of LSSS T_out temperature limit
                        if (T_out>=T_out_limit .AND. Toutflag==.false.) then
                            Toutflag=.true.
                            ! Calculate core temperatures for hot channel
                            call pebblebedcore(POWER,Q_core,T_in,T_out,T_coolant_max,T_core_max,2)
                            if (T_coolant_max <= T_coolant_limit .AND. coolantflag==.false.) then
                                write(*,'(A,F7.2,A)') "The maximum   coolant hot channel  temperature exceeds ", T_coolant_limit, " at:"
                                write(*,'(F5.2,A,F7.2,A,F7.2,A,F7.2)') POWER/1.0E6, " MW    T_in = ", T_in, "   T_out = ", T_out_avg, "   TCMAX = ", T_coolant_max
                                write(*,*)
                                coolantflag=.true.
                            end if
                        end if
                    end if
                end do
            end if
        end do
        !================================================================================
        ! Find power level for LSSS limits for given T_out and mass flow rate 
        !================================================================================
        ! Reset flags to .false.
        fuelflag=.false. 
        !coolantflag=.false.
        Toutflag=.false.
        !write(*,*)
        !write(*,'(A,F5.2,A,F7.2,A)') "For primary mass flow rate of ", W_core, " and outlet temperature of ", T_out_limit, ":"
        do T_in=T_in_fuelrange,T_out_limit,0.01_8
            ! Check if all LSSS limits met
            if (Toutflag==.true. .AND. fuelflag==.true.) then
                !write(*,*) "LSSS Fuel Limit Found"
                EXIT
            else          
                ! Loop through power level to find T_out=T_out_limit
                ! Reset flags to .false.
                Toutflag=.false.
                do POWER=POWER_fuelrange,50.0E6,0.01E6
                    if (Toutflag==.false.) then
                        ! Calculate core temperatures for average channel
                        call pebblebedcore(POWER,Q_core,T_in,T_out,T_coolant_max,T_core_max,1)
                        T_out_avg=T_out
                        ! Check for exceedance of LSSS T_out temperature limit
                        if (T_out>=T_out_limit .AND. Toutflag==.false.) then
                            Toutflag=.true.
                            ! Calculate core temperatures for hot channel
                            call pebblebedcore(POWER,Q_core,T_in,T_out,T_coolant_max,T_core_max,2)
                            if (T_core_max <= T_fuel_limit .AND. fuelflag==.false.) then
                                write(*,'(A,F7.2,A)') "The maximum    fuel hot channel    temperature exceeds ", T_fuel_limit, " at:"
                                write(*,'(F5.2,A,F7.2,A,F7.2,A,F7.2)') POWER/1.0E6, " MW    T_in = ", T_in, "   T_out = ", T_out_avg, "   TFMAX = ", T_core_max 
                                write(*,*)
                                fuelflag=.true.
                            end if
                        end if
                    end if
                end do
            end if
        end do
        
    END SUBROUTINE pebblebedLSSSloop
    
    
    
    
    
    
    
    !================================================================================
    !  FUNCTION: POROSITY OF PEBBLE PED
    !================================================================================
    ! INPUT     ::  
    ! OUTPUT    ::  Porosity of the pebble bed
    ! REFERENCE ::  
    !================================================================================
    REAL(8) FUNCTION porosity_PB(R)
        IMPLICIT NONE
        REAL(8), INTENT(IN)     :: R     ! RADIUS OF AREA WHERE PEBBLES ARE [m]   
                
        IF (R>1.0.AND.R<=1.886) THEN
            porosity_PB=1.0-2.0/3.0*R**(-3.0)/(2.0/R-1.0)**0.5
        ELSEIF (R>1.886.AND.R<2.033) THEN
            porosity_PB=1.8578-0.6649*R
        ELSEIF (R>=2.033) THEN
            porosity_PB=0.151/(R-1.0)+0.36
        ENDIF
                
        porosity_PB=1.0-0.6806 !REVISE duty ratio 0.32 !! SET POROSITY BY SINAP DESIGN !!
                
    END FUNCTION porosity_PB
            
            
    !================================================================================
    !  FUNCTION: CALCULATE HEAT TRANSFER COEFFICIENT IN PEBBLE BED
    !================================================================================
    ! INPUTS
    !   W       :: MASS FLOW RATE [kg/s]
    !   DC      :: CHANNEL DAIMETER [m]
    !   DP      :: PEBBLE DAIMETER [m]
    !   T       :: COOLANT TEMPERATURE [Celcius]
    !   N       :: LOGICAL TRIP
    !       N=1 Wakao Correlation 
    !       N=2 Kunii and Levenspiel Correlation
    ! OUTPUT    
    !   HTC_PB  ::  Heat transfer coefficient of the pebble bed [W/m^2-K]            
    ! REFERENCE
    !
    !================================================================================
    REAL(8) FUNCTION HTC_PB(W,DC,DP,T,N)
        IMPLICIT NONE
        REAL(8), INTENT(IN)     :: W    ! Mass flow rate [kg/s]
        REAL(8), INTENT(IN)     :: DC   ! Channel diameter [m]
        REAL(8), INTENT(IN)     :: DP   ! Pebble diameter [m]
        REAL(8), INTENT(IN)     :: T    ! Temperature [Celcius]
        INTEGER, INTENT(IN)     :: N    ! Correlation Selector
        REAL(8)                 :: G    ! Mass flux = mass flow rate / area
        REAL(8)                 :: Re   ! Reynolds number
        REAL(8)                 :: Pr   ! Prandtl number
        REAL(8)                 :: Nu   ! Nusselt number
                
        ! AVERAGE FLUID PROPERTIES
        G=W/(PI*(DC/2.0)**2.0)          ! [kg/m^2-s]
        Re=G*DP/flibe_mu(T)
        Pr=flibe_cp(T)*flibe_mu(T)/flibe_k(T)

        ! HEAT TRANSFER CORRELATION
        SELECT CASE(N)
        CASE(1)
            ! Wakao Kaviany's heat transfer handbook Ud=SUPERFICIAL VELOCITY Re:15-8500
            Nu=2.0+1.1*Pr**(1.0/3.0)*Re**0.6 
            IF (Re<=15.0.OR.Re>=8500) THEN
                WRITE(*,*) Re,"Re overflow - Wakao"
            ENDIF
        CASE(2)
            !Kunii and Levenspiel
            Nu=2.0+1.8*Pr**(1.0/3.0)*Re**0.5 
        END SELECT

        HTC_PB=Nu*flibe_k(T)/DP
                
    END FUNCTION HTC_PB
            
            
    !================================================================================
    !  FUNCTION: CALCULATE FRICTION FACTOR IN PEBBLE BED
    !================================================================================
    ! INPUTS
    !   W       :: MASS FLOW RATE [kg/s]
    !   DC      :: CHANNEL DAIMETER [m]
    !   DP      :: PEBBLE DAIMETER [m]
    !   T       :: COOLANT TEMPERATURE [Celcius]
    !   N       :: LOGICAL TRIP
    !       N=1 Ergu's Law 
    !       N=2 
    ! OUTPUT    
    !   FF_PB  ::  friction factor of the pebble bed [none]           
    ! REFERENCE
    !
    !================================================================================
    REAL(8) FUNCTION FF_PB(W,DC,DP,T,N)
        IMPLICIT NONE
        REAL(8), INTENT(IN)     :: W    ! Mass flow rate [kg/s]
        REAL(8), INTENT(IN)     :: DC   ! Channel diameter [m]
        REAL(8), INTENT(IN)     :: DP   ! Pebble diameter [m]
        REAL(8), INTENT(IN)     :: T    ! Temperature [?]
        INTEGER, INTENT(IN)     :: N    ! Correlation Selector
        REAL(8)                 :: E    ! Porosity
        REAL(8)                 :: G    ! Mass flux = mass flow rate / area
        REAL(8)                 :: Re   ! Reynolds number
        REAL(8)                 :: Pr   ! Prandtl number
        REAL(8)                 :: A
        REAL(8)                 :: B

        ! AVERAGE FLUID PROPERTIES
        E=porosity_pb(DC/DP)
        G=W/(PI*(DC/2.0)**2.0)      ! [kg/m^2-s]
        Re=G*DP/flibe_mu(T)
        ! Pr=flibe_cp(T)*flibe_mu(T)/flibe_k(T)

        ! FRICTION FACTOR CORRELATION
        SELECT CASE(N)
        CASE(1)
            ! Ergun's law
            A=180.0
            B=1.8
            FF_PB=(1.0-E)/E**(3.0)*(A*(1-E)/Re+B) 
        CASE(2)
        END SELECT

    END FUNCTION FF_PB
            
END MODULE pebblebed