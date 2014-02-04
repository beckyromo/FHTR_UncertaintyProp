!********************************************************************************
!
!  MODULE:  prismatic
!
!  PURPOSE: Contains functions for prismatic core calculations and subroutines 
!           to do LSSS calculations for a prismatic core
!
!  FUNCTIONS:
!  HTC
!  FF
!
!  SUBROUTINES:
!  prismaticcore        calculates core temperatures for given power, mass flow 
!                       rate, and inlet temperature for either average or hot 
!                       channel, prints core temperature results and outputs the 
!                       outlet temperature, maximum coolant temperature, and 
!                       maximum core/fuel temperature
!  prismaticLSSSloop    calls subroutine pebblebedcore iterating on power to 
!                       print LSSS for given fuel, coolant, and outlet 
!                       temperature limits
!
!******************************************************************************** 
MODULE prismatic
            
    USE global
        
    IMPLICIT NONE

 
CONTAINS
      
    !================================================================================
    !  SUBROUTINE: prismaticcore
    !================================================================================
    ! INPUTS
    !   POWER           :: Power [W]
    !   W_core          :: Mass flow rate [kg/s]
    !   Q_core          :: Volumetric flow rate [m^3/s]
    !   T_in            :: Inlet temperature [Celsius]
    !   channel         :: Selects for average (1) or hot (2) channel calculation
    ! OUTPUTS
    !   output (output_type) -->
    !       T_out           :: Outlet temperature [Celsius]
    !       T_coolant_max   :: Maximum coolant temperature [Celsius]
    !       T_core_max      :: Maximum core/fuel temperature [Celsius]
    ! REFERENCES
    !   Yao's Code, which uses power distribution data, LS-VHTR
    !================================================================================
    SUBROUTINE prismaticcore(POWER,Q_core,T_in,output,channel)
    
        USE global
        USE flibeprop, ONLY: flibe_cp, flibe_k, flibe_rho, flibe_mu, flibe_enthalpy, flibe_temperature
        USE trisoprop, ONLY: k_graphite,k_TRISO, k_SiC, K_densePyC, k_PyC, k_UO2, k_TRISOlayer
        !USE pipes, ONLY: FF_PIPE, dP_PIPE
        
        IMPLICIT NONE

        !================================================================================
        ! Declare Variables
        !--> Declare/Initialize, Declare arrays, initialize array inputs
        !================================================================================
        
        integer,intent(in)  :: channel          ! 1 for average, 2 for hot channel
        
        integer             :: I                ! Loop counter
        integer             :: m                ! Loop counter
        integer             :: N_core           ! Number of nodes(CVs) the core is split into
        integer             :: N_blocks         ! Number of fuel assembly blocks
        integer             :: N_cool           ! Number of coolant channels per block
        integer             :: N_fuel           ! Number of fuel channels per block
        
        real(8)             :: AFP(21)
        
        ! UO2,PyC,DPyC,SiC,DPyC,TRISO+Graphite,Graphite
        real(8)             :: RPF(7)           ! Fuel radii [m]
        real(8)             :: DPF(7)           ! Fuel diameter [m]
        real(8)             :: KPF(7)           ! Fuel thermal conductivity [W/m-C]
        real(8)             :: RHOPF(7)         ! Fuel density [kg/m^3]
        real(8)             :: CPPF(7)          ! Fuel heat capacity [J/kg-C]

        real(8)             :: DD1              ! Core: Face to face distance 1 of hexagonal core
        real(8)             :: DD2              ! Core: Face to face distance 2 of hexagonal core
        real(8)             :: H_core           ! Core: height [m]
        real(8)             :: A_core           ! Core: area of the core [m^2]
        real(8)             :: H_block          ! Height of fuel assemply block [m]
        real(8)             :: D_core           ! Core: hydraulic diameter of the core [m]
        real(8)             :: D_fuel           ! Diameter of fuel channel [m]
        real(8)             :: D_cool           ! Dimaeter of coolant channel [m]       
        real(8)             :: d_cell           ! Unit cell: flat to flat distance (2x short radius) [m] d=sqrt(3)*t
        real(8)             :: t_cell           ! Unit cell: side length of hexagon
        real(8)             :: A_cell           ! Unit cell: area
        real(8)             :: A_fuel           ! Unit cell: fuel area
        real(8)             :: A_cool           ! Unit cell: coolant area
        real(8)             :: A_graphite       ! Unit cell: graphite block area
         
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
        real(8)             :: W_core           ! Mass flow rate in the core [kg/s] based on input
        real(8)             :: W                ! Mass flow rate in the core [kg/s]
        real(8),intent(in)  :: T_in             ! Inlet temperature of the core [Celsius]
        
        real(8)             :: NFC              ! Number of fuel channels (1 for unit cell)
        real(8)             :: NTRISO           ! Number of TRISO particles per fuel channel control volume in unit cell
        real(8)             :: PF               ! TRISO particle packing fraction
        real(8)             :: PPFC             ! Power per fuel channel [W]
        real(8)             :: PPT              ! Power per TRISO particle [W]
        
        real(8)             :: VolPF1           ! Volume of first layer of TRISO (UO2) [m^3]
        real(8)             :: VolPF5           ! Volume of first 5 layers of TRISO [m^3]
        real(8)             :: VolPF6           ! Volume of fuel matrix (TRISOs and graphite) of pebble [m^3]
        real(8)             :: VolPF7           ! Volume of pebble (fuel matrix plus graphite cladding) [m^3]
       
        real(8)             :: FF_core          ! Friction factor in core 
        real(8)             :: HTC_core         ! Heat transfer coefficient in core [W/m^2-C]
        real(8)             :: A_HTC            ! Area of core with convection heat transfer [m^2]
        real(8)             :: ACTUAL_POWER     ! Actual Power [W]
        real(8)             :: dP_core          ! Core pressure drop [Pa] 
        real(8),allocatable :: enthalpy_core(:) ! Enthalpy in each node [J?/kg]
        real(8),allocatable :: T_w_core(:)      ! Pebble surface temperature at each node  [Celsius]
        real(8),allocatable :: T_CL_core(:)     ! Centerline temperature at each node  [Celsius]
        real(8),allocatable :: T_CL_TRISO(:)    ! Centerline temperature of centerline TRISO at each node [C]
        real(8)             :: T_coolant_max    ! Maximum coolant temperature [Celsius]
        real(8)             :: T_core_max       ! Maximum core temperature  [Celsius]
        real(8)             :: T_out            ! Outlet temperature of the core [Celsius]
        type(output_type)   :: output           ! T_out, T_coolant_max, T_out [Celsius]
        real(8)             :: TC               ! Temperature of coolant in core [Celsius]
        real(8)             :: TG               ! Temperature at inner edge of graphite cladding [Celsius]
        real(8)             :: TW               ! Temperature at the pebble surface [Celsius]
        real(8)             :: T_temp           ! Temperature place holder [Celsius]
        real(8)             :: TOL              ! While loop error tolerance
        real(8)             :: U                ! Thermal resistance
        real(8)             :: mu_w             ! Viscosity at the wall
        real(8)             :: mu_w_temp        ! Viscosity at the wall place holder
        real(8)             :: R_cool           ! Radius of coolant in unit cell
        real(8)             :: R_cg             ! Radius of graphite plus coolant in unit cell
        real(8)             :: R_cgf            ! Radius of unit cell (coolant plus graphite plus fuel)
        real(8)             :: k_mod            ! Thermal conductivity factor multiplier for fuel (cyl vs annular)

        
        
!========================================================================================
! Initialize Variables / Inputs
!========================================================================================
        
        ! POWER DISTRIBUTION FROM SINAP REPORT
        DATA AFP/1.3522E-002,4.7941E-002,4.6349E-002,4.7030E-002,4.8748E-002,5.0615E-002,5.2291E-002,&
                 5.3749E-002,5.4762E-002,5.5298E-002,5.5540E-002,5.5100E-002,5.4286E-002,5.3050E-002,&
                 5.1321E-002,4.9228E-002,4.7067E-002,4.4800E-002,4.2871E-002,4.2380E-002,3.4045E-002/ ! SUMS TO 1
        
        !DATA AFP/0.0476,0.0476,0.0476,0.0476,0.0476,0.0476,0.0476,0.0476,0.0476,0.0476,0.0476,0.0476,&
         !       0.0476,0.0476,0.0476,0.0476,0.0476,0.0476,0.0476,0.0476,0.0476/ ! SUMS TO 1
        
        !================================================================================
        ! ENGINERRING HOT CHANNEL FACTORS
        !
        ! See ICONE paper of Yao's for MITR values, see Yao's report for HTC-10 values
        !--------------------------------------------------------------------------------       
        FFDF=0.8_8 
        FH=1.173_8          ! 1.173 is only MITR value, 1.172 uses HTC-10 pebblebed
        FDTW=1.275_8        ! 1.275 is only MITR value, 1.268 uses HTC-10 pebblebed
        FDTF=1.123_8        ! 1.123 is only MITR value, 1.172 uses HTC-10 pebblebed        
        FCORE=0.965_8       ! NOT USED
        FFUEL=0.940_8       ! NOT USED
        FKZ=1.4_8
        FKTRISO=1.0_8       ! NOT USED
        !================================================================================

        ! FUEL: UO2,PyC,DPyC,SiC,DPyC,TRISO+Graphite,Graphite (
        DATA	DPF/0.350E-3,0.550E-3,0.620E-3,0.690E-3,0.770E-3,25.0E-3,30.0E-3/   ! Diameter [m]
        DPF(6)=D_fuel
        RPF=DPF/2.0                                                                 ! Radius [m]
        DATA	KPF/2.5,0.5,4.0,13.9,4.0,30.0,30.0/                                 ! Thermal conductivity [W/m-C] NOTE OVERWRITTEN
        !DATA	RHOPF/10.4E3,1.0E3,1.85E3,3.2E3,1.8E3,1.74E3,1.74E3/                 ! Density [kg/m^3]
        
        ! Core parameters
        DECAY_HEAT=0.0_8        ! [W]
        POWER_NORM=1.0_8        ! [W]
        PF=0.075_8              ! TRISO particle packing fraction
        
        ! Core Geometry
        DD1=134.69E-2_8         ! Face to face distance 1 [m]
        DD2=139.0E-2_8          ! Face to face distance 2 [m]
        H_block=0.793           ! Height of a block [m]
        H_core=H_block*1.0_8   ! Height of the core [m] !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! one block high core
        N_core=21               ! Number of nodes(CVs) in the core
        D_fuel=12.7E-3_8        ! Diameter of fuel channel [m]
        D_cool=9.53E-3_8        ! Diameter of coolant channel [m] or revised value of 1.4 cm
        N_blocks=61             ! Number of fuel assembly blocks
        N_cool=108              ! Number of coolang channels per block
        N_fuel=216              ! Number of fuel channels per block
        
        ! Unit Cell
        !d_cell=0.24_8                          ! Flat to flat distance (2x short radius) [m] d=sqrt(3)*t
        t_cell=18.8E-3_8                        ! Unit cell side length of hexagon [m]
        A_cell=3.0*sqrt(3.0)/2.0*t_cell**2.0    ! sqrt(3)/2.0*d_cell**2.0
        A_fuel=PI*D_fuel**2.0/4.0*6.0/3.0
        A_cool=PI*D_cool**2.0/4.0
        A_graphite=A_cell-A_fuel-A_cool
        ! Convert unit cell to concentric circles
        R_cool=D_cool/2
        R_cg=sqrt((A_cool+A_graphite)/PI)
        R_cgf=sqrt(A_cell/PI)
        ! Ratio of exact annular to cyl soln to account for smallar radius in concentric circle model (Davis and Hawkes 2006)
        k_mod=(2.0*R_cgf**2.0*log(R_cgf/R_cg)-(R_cgf**2.0-R_cg**2.0)) / (D_fuel/2.0)**2.0 
    
 
        ! Calculated inputs
        VolPF1=4.0/3.0*PI*RPF(1)**3.0   ! Fuel kernel volume
        VolPF5=4.0/3.0*PI*RPF(5)**3.0   ! TRISO volume
        VolPF6=A_fuel*H_core/N_core     ! Fuel matrix volume (volume of fuel channel)
        NTRISO=VolPF6*PF/VolPF5         ! Number of TRISO particles per fuel channel CV in unit cell = PF*V_channel/V_TRISO
        
        !A_core=DD1**2 -((DD1*2.0**0.5-DD2)/2.0**0.5)**2*2.0-0.35*0.35
        D_core=9.2 !(4.0*A_core/PI)**0.5                                    ! HYDRAULIC RADIUS OF THE CORE [m]
        A_core=PI*D_core**2.0/4.0                                           ! Used to calculate pressure drop!
        
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
    W_core=Q_core*flibe_rho(T_in)
    IF (channel==1) THEN
        W=W_core
        FH=1.0_8
        FKZ=1.0_8
        FDTW=1.0_8
        FDTF=1.0_8
        FFDF=1.0_8
    ELSE IF (channel==2) THEN
        W=W_core*FFDF
    END IF
            
 
       
!================================================================================       
! Core Loop Calculations
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
            !IF (channel==1) THEN
            !    FF_core=FF(W,D_cool,TC,1)
            !    dP_core=dP_core+0.5*FF_core*(H_core/N_core)*(W/A_core)**2.0/(flibe_rho(TC)*D_cool) &
            !            + flibe_rho(TC)*GRAVITY*H_core/N_core
            !END IF
            
            ! =======================================================
            ! Find wall temperature and iterate with HTC correlations
            ! =======================================================
            TOL=1.0e-6
            TW=TC
            mu_w=flibe_mu(TW)
            mu_w_temp=0
            DO WHILE (abs(mu_w_temp - mu_w) > TOL) 
            ! Calculate heat transfer coefficient in node 
            ! For laminar flow select 1 for Nu=4.364, 2 for Sieder-Tate, and 3 for Rea (entrance length)
            HTC_core=HTC(W/(N_cool*N_blocks),D_cool,H_core,H_core/N_core*I,TC,TW,3) 
            mu_w_temp=flibe_mu(TW)
            ! Find wall temperature
            A_HTC=PI*D_cool*H_core/N_core*N_cool*N_blocks
            TW=ACTUAL_POWER*AFP(I)*FKZ/(HTC_core*A_HTC) + TC
            TW=(TW-TC)*FDTW+TC
            mu_w=flibe_mu(TW)
            END DO
            T_w_core(I)=TW
            ! --------------------------------------------------------
            
            
            
            ! Check if temperature is maximum
            IF (TW>T_coolant_max) THEN
                T_coolant_max=TW
            END IF
            
            ! Calculate power in a fuel channel and centerline TRISO particle
            NFC=2.0_8                                           ! Number of fuel channels: 2 for unit cell
            PPFC=ACTUAL_POWER/(N_blocks*N_fuel)*AFP(I)*FKZ      ! Power per fuel channel 
            PPT=PPFC/NTRISO                                     ! Power per TRISO particle

            ! Calculate thermal conductivties and find temperatures
            TOL=1.0e-5           
                        
            ! GRAPHITE CLADDING
            T_temp=TW
            TG=0.0
            do while (abs(T_temp-TG)>TOL)
                T_temp=TG
                KPF(7)=k_graphite( (TW+T_temp)/2.0 )
                U=log( R_cg/(D_cool/2) ) / (2.0*PI*KPF(7)*H_core/N_core) ! Annular model
                TG=2.0*PPFC*U + TW
            end do
            
            ! FUEL MATRIX (TRISOs IN GRAPHITE MATRIX)
            T_temp=TG
            T_CL_core(I)=0.0
            do while (abs(T_temp-T_CL_core(I))>TOL)
                T_temp=T_CL_core(I)
                KPF(6)=k_TRISO( (TG+T_temp)/2.0 ) * k_mod
                ! Unit cell model decoupled
                !U=1.0/(4.0*PI*KPF(6)*H_core/N_core)                
                ! T&K Eq. 8-97 for annular model
                !U=1.0/(4.0*PI*KPF(6)*H_core/N_core)*(R_cg/R_cgf)**2.0*( (1.0-(R_cgf/R_cg)**2.0) &
                    !- ( (R_cgf/R_cg)**2.0*log((R_cg/R_cgf)**2.0) ) ) 
                ! My Eq. 8-97 (matches T&K)
                U=1.0/(4.0*PI*KPF(6)*H_core/N_core)*(2.0*log(R_cgf/R_cg)-(1.0-(R_cg/R_cgf)**2.0))
                T_CL_core(I)=2.0*PPFC*U + TG
            end do
                           
            ! CENTERLINE TRISO PARTICLE
            T_temp=T_CL_core(I)
            U=0.0
            ! Loop through TRISO coating layers to get their total thermal resistance
            do m=2,5
                U=1.0/(4.0*PI*k_trisolayer(m,T_temp)) * (1.0/RPF(m-1) - 1.0/RPF(m)) + U
            end do
            ! Add fuel kernel to total thermal resistance and calculate centerline temperature of centermost TRISO particle
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
                write(10,*) "  HEIGHT      TC       TW       TG      TCL     TCLT"
            end if
            write(10,'(F9.3,F9.3,F9.3,F9.3,F9.3,F9.3)') I*H_core/N_core, TC, TW, TG, T_CL_core(I), T_CL_TRISO(I)
            ! write to command line        
            ! write(*,*) PPP, PPT, TC, TW, TG, T_CL_core(I), T_CL_TRISO(I)
            
        end do
        
        T_out=flibe_temperature(enthalpy_core(N_core))    
        
        output%T_out=T_out
        output%T_coolant_max=T_coolant_max
        output%T_core_max=T_core_max
        
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
        write(10,*) "   POWER      W      T_IN    T_OUT    TC_MAX   TF_MAX"
        write(10,'(F9.2,F9.3,F9.3,F9.3,F9.3,F9.3)') &
                    POWER/1.0E6,W,T_in,T_out,T_coolant_max,T_core_max
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
        
    END SUBROUTINE prismaticcore
    
    
    
    
    
    
    
    
    
    
    !================================================================================
    !  SUBROUTINE: prismaticLSSSloop  [W/m-K]
    !================================================================================
    ! INPUTS    
    !   input (input_type) -->
    !       POWER               :: not used
    !       W_core              :: not used
    !       Q_core
    !       T_in
    !   limits (limits_type) -->
    !       T_fuel_limit        :: LSSS temoerature limit for fuel (max fuel temp in hot channel) 
    !       T_coolant_limit     :: LSSS temperature limit for coolant (max coolant temp in hot channel)
    !       T_in_limit          :: LSSS temperature limit for T_in (min T_in temp in average channel)
    !       T_out_limit         :: LSSS temperature limit for T_out (max T_out temp in average channel)
    ! OUTPUTS   
    !   output (output_type) -->
    !       T_out               :: Outlet temperature [Celsius]
    !       T_coolant_max       :: Maximum coolant temperature [Celsius]
    !       T_core_max          :: Maximum core/fuel temperature [Celsius]
    ! REFERENCES
    !   none
    !================================================================================
    SUBROUTINE prismaticLSSSloop(input,limits,LSSS)
        
        IMPLICIT NONE
        
        integer                         :: I                    ! Loop counter
        logical                         :: fuelflag             ! Flag if max fuel temperature has been met in hot channel
        logical                         :: coolantflag          ! Flag if max coolant temperature has been met in hot channel
        logical                         :: Toutflag             ! Flag if max coolant outlet temperature has been met in average channel
        type(input_type),intent(in)     :: input    
        real(8)                         :: POWER                ! Reactor power [W]
        real(8)                         :: POWER_HIGH           ! Top power for loop [W]
        real(8)                         :: POWER_LOW            ! Low power for loop [W]
        !real(8)                         :: W_core               ! Mass flow rate in the core [kg/s]
        real(8)                         :: Q_core               ! Volumetric flow rate in the core [m^3/s]        
        real(8)                         :: T_in                 ! Inlet temperature of the core [Celsius]        
        real(8)                         :: T_out_avg            ! Holder for T_out in average channel
        !real(8)                         :: T_out                ! Maximum coolant temperature [Celsius]
        !real(8)                         :: T_coolant_max        ! Maximum core temperature  [Celsius]
        !real(8)                         :: T_core_max           ! Outlet temperature of the core [Celsius]
        real(8)                         :: step                 ! Step size for power do loop
        real(8)                         :: T_in_coolantrange    
        real(8)                         :: POWER_coolantrange   
        real(8)                         :: T_in_fuelrange       
        real(8)                         :: POWER_fuelrange
        type(LSSS_type),intent(out)     :: LSSS
        type(limits_type),intent(in)    :: limits
        real(8)                         :: T_fuel_limit         ! LSSS temoerature limit for fuel (max fuel temp in hot channel) 
        real(8)                         :: T_coolant_limit      ! LSSS temperature limit for coolant (max coolant temp in hot channel)
        real(8)                         :: T_in_limit           ! LSSS temperature limit for T_in (min T_in temp in average channel)
        real(8)                         :: T_out_limit          ! LSSS temperature limit for T_out (max T_out temp in average channel)
                
        write(*,*)
        write(*,*) "Calculating LSSS Data Points..."
        
        
        T_fuel_limit=limits%T_fuel_limit
        T_coolant_limit=limits%T_coolant_limit
        T_in_limit=limits%T_in_limit
        T_out_limit=limits%T_out_limit
        T_in=input%T_in
        Q_core=input%Q_core
        POWER_LOW=1.0E6
        POWER_HIGH=70.0E6
        
        
          
        !write(20,*) "                                               LSSS Results"
        !write(20,*)
        !write(20,'(A54)',ADVANCE='NO') "                           AVG CHANNEL                             "    
        !write(20,'(A54)')              "                           HOT CHANNEL                             "
        !write(20,'(A54)',ADVANCE='NO') " =================================================================="    
        !write(20,'(A54)')              " =================================================================="
        !write(20,'(A54)',ADVANCE='NO') "   POWER      W      T_IN    T_OUT    TC_MAX   TF_MAX"    
        !write(20,'(A54)')              "   POWER      W      T_IN    T_OUT    TC_MAX   TF_MAX"
        !write(20,'(A54)',ADVANCE='NO') " ------------------------------------------------------------------"    
        !write(20,'(A54)')              " ------------------------------------------------------------------"
            
        !================================================================================
        ! Find power level for LSSS limits for given T_in and mass flow rate 
        !================================================================================
        ! Reset flags to .false.
        fuelflag=.false.    
        coolantflag=.false.
        Toutflag=.false.
        !write(*,*)
        !write(*,'(A,F5.2,A,F7.2,A)') "For primary mass flow rate of ", Q_core, " and  inlet temperature of ", T_in, ":"
        !write(*,'(A)') "======================================================================"
        !write(*,*)
        do POWER=POWER_LOW,POWER_HIGH,0.1E6
            !================================================================================
            ! Calculate core temperatures for average channel
            call prismaticcore(POWER,Q_core,T_in,output,1)
            !write(20,'(F9.2,F9.3,F9.3,F9.3,F9.3,F9.3)', ADVANCE='NO') POWER/1.0E6,Q_core,T_in,output
            T_out_avg=output%T_out
            ! Check for exceedance of LSSS T_out temperature limit
            if (output%T_out>=T_out_limit .AND. Toutflag==.false.) then
                !write(*,'(A,F7.2,A)') "The maximum outlet average channel temperature exceeds ", T_out_limit, " at:"
                !write(*,'(F5.2,A,F7.2,A,F7.2)') POWER/1.0E6, " MW    T_in = ", T_in, "   T_out = ", output%T_out
                !write(*,*) 
                Toutflag=.true.
                LSSS%INmaxoutPOWER=POWER/1.0E6
                LSSS%INmaxoutTin=T_in
                LSSS%INmaxoutTout=output%T_out
                LSSS%minPOWER=POWER/1.0E6
            end if
            !================================================================================
            ! Calculate core temperatures for hot channel
            call prismaticcore(POWER,Q_core,T_in,output,2)
            !write(20,'(F9.2,F9.3,F9.3,F9.3,F9.3,F9.3)') POWER/1.0E6,Q_core,T_in,output
            ! Check for exceedance of LSSS coolant temperature limit   
            if (output%T_coolant_max>=T_coolant_limit .AND. coolantflag==.false.) then
                !write(*,'(A,F7.2,A)') "The maximum  coolant hot channel   temperature exceeds ", T_coolant_limit, " at:"
                !write(*,'(F5.2,A,F7.2,A,F7.2,A,F7.2)') POWER/1.0E6, " MW    T_in = ", T_in, "   T_out = ", output%T_out, "   TCMAX = ", output%T_coolant_max
                !write(*,'(A,F7.2)') "T_out_avg = ", T_out_avg
                !write(*,*) 
                coolantflag=.true.
                LSSS%INmaxcoolPOWER=POWER/1.0E6
                LSSS%INmaxcoolTin=T_in
                LSSS%INmaxcoolTout=output%T_out
                LSSS%INmaxcoolTmax=output%T_coolant_max
                LSSS%INmaxcoolToutavg=T_out_avg
                if (POWER/1.0E6 < LSSS%minPOWER) then
                    LSSS%minPOWER=POWER/1.0E6
                end if
            
            end if
		    ! Check for exceedance of LSSS fuel temperature limit
            if (output%T_core_max>=T_fuel_limit .AND. fuelflag==.false.) then
                !write(*,'(A,F7.2,A)') "The maximum    fuel hot channel    temperature exceeds ", T_fuel_limit, " at:"
                !write(*,'(F5.2,A,F7.2,A,F7.2,A,F7.2)') POWER/1.0E6, " MW    T_in = ", T_in, "   T_out = ", output%T_out, "   TFMAX = ", output%T_core_max
                !write(*,'(A,F7.2)') "T_out_avg = ", T_out_avg
                !write(*,*) 
                fuelflag=.true.
                LSSS%INmaxfuelPOWER=POWER/1.0E6
                LSSS%INmaxfuelTin=T_in
                LSSS%INmaxfuelTout=output%T_out
                LSSS%INmaxfuelTmax=output%T_core_max
                LSSS%INmaxfuelToutavg=T_out_avg
                if (POWER/1.0E6 < LSSS%minPOWER) then
                    LSSS%minPOWER=POWER/1.0E6
                end if
            end if
            
        end do
   
        
        
        !================================================================================
        ! Find power level for LSSS limits for given T_out and mass flow rate :: Find T_in and POWER range
        !================================================================================
        ! Reset flags to .false.
        fuelflag=.false. 
        coolantflag=.false.
        Toutflag=.false.
        !write(*,*)
        !write(*,'(A,F5.2,A,F7.2,A)') "For primary mass flow rate of ", Q_core, " and outlet temperature of ", T_out_limit, ":"
        !write(*,'(A)') "======================================================================"
        do T_in=470.0_8,T_out_limit,1.0_8 
            ! Check if all LSSS limits met
            if (Toutflag==.true. .AND. coolantflag==.true. .AND. fuelflag==.true.) then
                write(*,*) "LSSS Limits Range Found"
                EXIT
            else          
                ! Loop through power level to find T_out=T_out_limit
                ! Reset T_out flags to .false.
                Toutflag=.false.
                do POWER=POWER_LOW,POWER_HIGH,1.0E6
                    if (Toutflag==.false.) then
                        ! Calculate core temperatures for average channel
                        call prismaticcore(POWER,Q_core,T_in,output,1)
                        T_out_avg=output%T_out
                        ! Check for exceedance of LSSS T_out temperature limit
                        if (output%T_out>=T_out_limit .AND. Toutflag==.false.) then
                            Toutflag=.true.
                            ! Calculate core temperatures for hot channel
                            call prismaticcore(POWER,Q_core,T_in,output,2)
                            if (output%T_coolant_max <= T_coolant_limit .AND. coolantflag==.false.) then
                                !write(*,*) T_in, T_out_avg
                                T_in_coolantrange=T_in-5.0
                                POWER_coolantrange=POWER
                                !write(*,'(A,F7.2,A,F5.2,A)') & 
                                !"For coolant LSSS start T_in at ", T_in_coolantrange, " and the power level at ", POWER/1.0E6, " MW" 
                                coolantflag=.true.
                            end if
                            if (output%T_core_max <= T_fuel_limit .AND. fuelflag==.false.) then
                                !write(*,*) T_in,  T_out_avg
                                T_in_fuelrange=T_in-7.0
                                POWER_fuelrange=POWER
                                !write(*,'(A,F7.2,A,F5.2,A)') & 
                                !"For fuel LSSS start T_in at ", T_in_fuelrange, " and the power level at ", POWER/1.0E6, " MW" 
                                !write(*,*) output%T_core_max
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
        !write(*,'(A,F5.2,A,F7.2,A)') "For primary mass flow rate of ", Q_core, " and outlet temperature of ", T_out_limit, ":"
        do T_in=T_in_coolantrange,T_out_limit,0.01_8 
            ! Check if all LSSS limits met
            if (Toutflag==.true. .AND. coolantflag==.true.) then
                !write(*,*) "LSSS Coolant Limit Found"
                EXIT
            else          
                ! Loop through power level to find T_out=T_out_limit
                ! Reset T_out flags to .false.
                Toutflag=.false.
                do POWER=POWER_coolantrange,POWER_HIGH,0.01E6_8
                    if (Toutflag==.false.) then
                        ! Calculate core temperatures for average channel
                        call prismaticcore(POWER,Q_core,T_in,output,1)
                        T_out_avg=output%T_out
                        ! Check for exceedance of LSSS T_out temperature limit
                        if (output%T_out>=T_out_limit .AND. Toutflag==.false.) then
                            Toutflag=.true.
                            ! Calculate core temperatures for hot channel
                            call prismaticcore(POWER,Q_core,T_in,output,2)
                            if (output%T_coolant_max <= T_coolant_limit .AND. coolantflag==.false.) then
                                !write(*,'(A,F7.2,A)') "The maximum   coolant hot channel  temperature exceeds ", T_coolant_limit, " at:"
                                !write(*,'(F5.2,A,F7.2,A,F7.2,A,F7.2)') POWER/1.0E6, " MW    T_in = ", T_in, "   T_out = ", T_out_avg, "   TCMAX = ", output%T_coolant_max
                                !write(*,*)
                                coolantflag=.true.
                                LSSS%OUTmaxcoolPOWER=POWER/1.0E6
                                LSSS%OUTmaxcoolTin=T_in
                                LSSS%OUTmaxcoolTout=output%T_out
                                LSSS%OUTmaxcoolTmax=output%T_coolant_max
                                LSSS%OUTmaxcoolToutavg=T_out_avg
                                if (POWER/1.0E6 < LSSS%minPOWER) then
                                    LSSS%minPOWER=POWER/1.0E6
                                end if
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
        !write(*,'(A,F5.2,A,F7.2,A)') "For primary mass flow rate of ", Q_core, " and outlet temperature of ", T_out_limit, ":"
        do T_in=T_in_fuelrange,T_out_limit,0.01_8
            ! Check if all LSSS limits met
            if (Toutflag==.true. .AND. fuelflag==.true.) then
                !write(*,*) "LSSS Fuel Limit Found"
                EXIT
            else          
                ! Loop through power level to find T_out=T_out_limit
                ! Reset flags to .false.
                Toutflag=.false.
                do POWER=POWER_fuelrange,POWER_HIGH,0.01E6
                    if (Toutflag==.false.) then
                        ! Calculate core temperatures for average channel
                        call prismaticcore(POWER,Q_core,T_in,output,1)
                        T_out_avg=output%T_out
                        ! Check for exceedance of LSSS T_out temperature limit
                        if (output%T_out>=T_out_limit .AND. Toutflag==.false.) then
                            Toutflag=.true.
                            ! Calculate core temperatures for hot channel
                            call prismaticcore(POWER,Q_core,T_in,output,2)
                            if (output%T_core_max <= T_fuel_limit .AND. fuelflag==.false.) then
                                !write(*,'(A,F7.2,A)') "The maximum    fuel hot channel    temperature exceeds ", T_fuel_limit, " at:"
                                !write(*,'(F5.2,A,F7.2,A,F7.2,A,F7.2)') POWER/1.0E6, " MW    T_in = ", T_in, "   T_out = ", T_out_avg, "   TFMAX = ", output%T_core_max 
                                !write(*,*)
                                fuelflag=.true.
                                LSSS%OUTmaxfuelPOWER=POWER/1.0E6
                                LSSS%OUTmaxfuelTin=T_in
                                LSSS%OUTmaxfuelTout=output%T_out
                                LSSS%OUTmaxfuelTmax=output%T_core_max
                                LSSS%OUTmaxfuelToutavg=T_out_avg
                                if (POWER/1.0E6 < LSSS%minPOWER) then
                                    LSSS%minPOWER=POWER/1.0E6
                                end if
                            end if
                        end if
                    end if
                end do
            end if
        end do
    
        !write(*,*) "The minimum LSSS power is:"
        !write(*,'(F7.2)') LSSS%minPOWER
        !write(*,*)
        
        
        write(*,*) "LSSS Calculation Complete"
        write(*,*)
    END SUBROUTINE prismaticLSSSloop
    
    
    
    
    
            
            
    !================================================================================
    !  FUNCTION: CALCULATE HEAT TRANSFER COEFFICIENT
    !================================================================================
    ! INPUTS
    !   W       :: MASS FLOW RATE [kg/s]
    !   DC      :: CHANNEL DAIMETER [m]
    !   DP      :: PEBBLE DAIMETER [m]
    !   T       :: COOLANT TEMPERATURE [Celsius]
    !   N       :: LOGICAL TRIP
    !       N=1 Wakao Correlation 
    !       N=2 Kunii and Levenspiel Correlation
    !       n=3 Dittus-Boelter
    ! OUTPUT    
    !   HTC  ::  Heat transfer coefficient of the pebble bed [W/m^2-K]            
    ! REFERENCE
    !
    !================================================================================
    REAL(8) FUNCTION HTC(W,DC,L,x,T,TW,N)
    
        USE flibeprop, ONLY: flibe_k, flibe_cp, flibe_mu
        
        IMPLICIT NONE
        REAL(8), INTENT(IN)     :: W        ! Mass flow rate [kg/s]
        REAL(8), INTENT(IN)     :: DC       ! Channel diameter [m]
        REAL(8), INTENT(IN)     :: L        ! Length of entire control volume channel [m]
        REAL(8), INTENT(IN)     :: x        ! Height of control volume channel [m]
        REAL(8), INTENT(IN)     :: T        ! Bulk Temperature [Celsius]
        REAL(8), INTENT(IN)     :: TW       ! Wall Temperature [Celsius]
        INTEGER, INTENT(IN)     :: N        ! Laminar Correlation Selector
        REAL(8)                 :: G        ! Mass flux = mass flow rate / area
        REAL(8)                 :: Re       ! Reynolds number
        REAL(8)                 :: Pr       ! Prandtl number
        REAL(8)                 :: Nu       ! Nusselt number
        REAL(8)                 :: mu_b     ! Viscosity of flibe in bulk
        REAL(8)                 :: mu_w     ! Viscosity of flibe at the wall
        REAL(8)                 :: xplus    ! x plus = 2*(L/D)/(Re*Pr)
                
        ! AVERAGE FLUID PROPERTIES
        G=W/(PI*(DC/2.0)**2.0)          ! [kg/m^2-s]
        Re=G*DC/flibe_mu(T)
        Pr=flibe_cp(T)*flibe_mu(T)/flibe_k(T)
        mu_b=flibe_mu(T)
        mu_w=flibe_mu(TW)
        
        !write(*,*) "Re=", Re
        !write(*,*) "Pr=", Pr
        !write(*,*) mu_b
        
        

        ! HEAT TRANSFER CORRELATION - Molten Salt (Based on HTC lit review in dissertation)

        if (Re<2000) then
            ! Select Laminar Correlation
            if (N==1) then
                Nu=4.364
            else if (N==2) then
                ! Sieder-Tate (laminar)
                Nu = 1.86 * ( Re*Pr*(DC/L) )**(1.0/3.0) * (mu_b/mu_w)**(0.14)          
            else if (N==3) then
                ! Rea 2008 (HMT6829)
                xplus=2.0*(x/DC)/(Re*Pr)
                !write(*,*) xplus
                if (xplus>0.003) then
                    Nu=4.364+0.263*(xplus/2.0)**(-0.506)*exp(-41.0*xplus/2.0)
                else 
                    Nu=1.302*(xplus/2.0)**(-1.0/3.0)
                end if
            end if
        else if (Re<10000) then
            ! Hausen (transitional)
            Nu = 0.037 * (Re**0.75 - 180.0) * Pr**0.42 * (1+(DC/L)**(2/3)) * (mu_b/mu_w)**(0.14)
        else if (Re<120000) then
            ! Gnielinski
            Nu = 0.012 * (Re**0.87 - 280.0) * Pr**0.4 * (1+(DC/L)**(2/3)) * (mu_b/mu_w)**(0.11)
        else 
            write(*,*) "Reynold's number out of range for HTC correlation.", Re
            Nu=100
        end if

        !write(*,*) "Nu=", Nu
        HTC=Nu*flibe_k(T)/DC
                
    END FUNCTION HTC
            
            
    !================================================================================
    !  FUNCTION: CALCULATE FRICTION FACTOR IN PEBBLE BED
    !================================================================================
    ! INPUTS
    !   W       :: MASS FLOW RATE [kg/s]
    !   DC      :: CHANNEL DAIMETER [m]
    !   T       :: COOLANT TEMPERATURE [Celsius]
    !   N       :: LOGICAL TRIP
    !       N=1 pipe flow
    !       N=2 churchill
    ! OUTPUT    
    !   FF  ::  friction factor [none]           
    ! REFERENCE
    !
    !================================================================================
    REAL(8) FUNCTION FF(W,DC,T,N)
    
        USE flibeprop, ONLY: flibe_mu
    
        IMPLICIT NONE
        REAL(8), INTENT(IN)     :: W    ! Mass flow rate [kg/s]
        REAL(8), INTENT(IN)     :: DC   ! Channel diameter [m]
        REAL(8), INTENT(IN)     :: T    ! Temperature [?]
        INTEGER, INTENT(IN)     :: N    ! Correlation Selector
        REAL(8)                 :: G    ! Mass flux = mass flow rate / area
        REAL(8)                 :: Re   ! Reynolds number
        REAL(8)                 :: A
        REAL(8)                 :: B
        REAL(8)                 :: e    ! Surface roughness
        

        ! AVERAGE FLUID PROPERTIES
        G=W/(PI*(DC/2.0)**2.0)      ! [kg/m^2-s]
        Re=G*DC/flibe_mu(T)

        ! FRICTION FACTOR CORRELATION
        SELECT CASE(N)
        CASE(1)
            ! Laminar and turblent pipe flow
            if (Re<=4000) then
                FF=64.0/Re
            else if (Re>4000 .AND. Re<30000) then ! Blasius 
                FF=0.316*Re**(-0.25)
            else if (Re>30000 .AND. Re<1000000) then ! McAdams
                FF=0.184*Re**(-0.2)
            else
                write(*,*) "Reynold's number out of range for friction factor correlation.", Re
            end if             
        CASE(2)
            ! Churchill
            if (Re<2100) then
                FF=64.0/Re
            else
                A=(-2.457*log( (7.0/Re)**0.9 + 0.27*e/DC ) )**16.0
                B=(37530.0/Re)**16.0
                e=30*10**(-6.0) ! [m] surface roughness of HastN from weighted average of Ni and Mo
                FF=8.0*( (8.0/Re)**12.0 + 1.0/((A+B)**1.5) )**(1.0/12.0)
            end if
        END SELECT

    END FUNCTION FF
            
END MODULE prismatic