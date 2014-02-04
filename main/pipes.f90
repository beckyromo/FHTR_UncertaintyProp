    !********************************************************************************
    !
    !  MODULE:  pipes
    !
    !  PURPOSE: Calculates pressure drop of the (4) pipes
    !
    !  SUBROUTINES:
    !
    !  FUNCTIONS:
    !  dP_PIPE          pressure drop in a pipe
    !  FF_PIPE          friction factor for pipe flow
    !
    !******************************************************************************** 
    MODULE pipes
        USE global
        USE flibeprop, ONLY: flibe_mu, flibe_rho, flibe_enthalpy, flibe_temperature
        IMPLICIT NONE

    CONTAINS
    
        !================================================================================
        !  FUNCTION: CALCULATE PRESSURE DROP IN A PIPE
        !================================================================================
        ! INPUTS
        !   W       :: MASS FLOW RATE [kg/s]
        !   D       :: PIPE DAIMETER [m]
        !   T       :: FLUID TEMPERATURE [Celcius]
        !   L       :: TOTAL LENGTH OF PIPE [m]          
        !   L_vert  :: LENGTH OF PIPE IN VERTICAL DIRECTION [m]
        ! OUTPUT    
        !   dP_PIPE ::  pressure drop in a pipe [none]           
        ! REFERENCE
        !
        !================================================================================   
        REAL(8) FUNCTION dP_PIPE(W,D,T,L,L_vert,K)
            IMPLICIT NONE
            REAL(8), INTENT(IN)     :: W        ! Mass flow rate [kg/s]
            REAL(8), INTENT(IN)     :: D        ! Pipe diameter [m]
            REAL(8), INTENT(IN)     :: T        ! Fluid Temperature [Celcius]
            REAL(8), INTENT(IN)     :: L        ! Total length of pipe [m]
            REAL(8), INTENT(IN)     :: L_vert   ! Lenght of pipe in vertical direction [m]
            REAL(8), INTENT(IN)     :: K        ! Hydraulic resistance coefficient
            REAL(8)                 :: dP_fric  ! Pressure drop - friction
            REAL(8)                 :: dP_grav  ! Pressure drop - gravity
            REAL(8)                 :: dP_form  ! Pressure drop - form loss
            REAL(8)                 :: G        ! Mass flux = mass flow rate / area
            REAL(8)                 :: density  ! Density [kg/m^3]
                
            G=W/(PI*(D/2.0)**2.0)               ! [kg/m^2-s]
            density=flibe_rho(T)                ! [kg/m^3]
            dP_fric=0.5*FF_PIPE(W,D,T)*L*G**2.0/(density*D)
            dP_grav=density*gravity*L_vert
            dP_form=0.5*K*G**2.0/density
                
            dP_PIPE=dP_fric+dP_grav+dP_form
            
        END FUNCTION dP_PIPE
            
        !================================================================================
        !  FUNCTION: CALCULATE FRICTION FACTOR IN PIPE
        !================================================================================
        ! INPUTS
        !   W       :: MASS FLOW RATE [kg/s]
        !   D       :: PIPE DAIMETER [m]
        !   T       :: COOLANT TEMPERATURE [Celcius] 
        ! OUTPUT    
        !   FF_PIPE  ::  friction factor of a pipe [none]           
        ! REFERENCE
        !
        !================================================================================
        REAL(8) FUNCTION FF_PIPE(W,D,T)
            IMPLICIT NONE
            REAL(8), INTENT(IN)     :: W    ! Mass flow rate [kg/s]
            REAL(8), INTENT(IN)     :: D    ! Pipe diameter [m]
            REAL(8), INTENT(IN)     :: T    ! Temperature [?]
            REAL(8)                 :: G    ! Mass flux = mass flow rate / area
            REAL(8)                 :: Re   ! Reynolds number
            REAL(8)                 :: Pr   ! Prandtl number                

            ! AVERAGE FLUID PROPERTIES
            G=W/(PI*(D/2.0)**2.0)           ! [kg/m^2-s]
            Re=G*D/flibe_mu(T)
            ! Pr=flibe_cp(T)*flibe_mu(T)/flibe_k(T)

            ! FRICTION FACTOR CORRELATION
            IF (Re<=1800.0) THEN
                ! Laminar Flow
                FF_PIPE=64.0/Re
            ELSEIF (Re<=2500.0) THEN
                FF_PIPE=0.3164/2500.0**0.25
            ELSE
                FF_PIPE=0.3164/Re**0.25
            ENDIF
                
        END FUNCTION FF_PIPE
            
        !================================================================================
        !  SUBROUTINE: CALCULATE dP of IHX and PIPES
        !================================================================================
        ! INPUTS
        !   H_core  :: HEIGHT OF THE CORE [m]
        !   W_core  :: MASS FLOW RATE [kg/s]
        !   T_out   :: TEMPERATURE AT OUTLET OF CORE
        !   POWER   :: CORE POWER [W]
        !================================================================================
        SUBROUTINE pipes_dP(H_core,W_core,T_out,POWER)
        
            IMPLICIT NONE
                       
            !================================================================================
            ! Declare Variables
            !
            ! Note: skip pipes/IHX calculations, moved to pipes module since already coded
            !
            !================================================================================
            
            real(8), intent(in)     :: H_core
            real(8), intent(in)     :: W_core
            real(8), intent(in)     :: T_out

            real(8), intent(in)     :: POWER
            
            integer             :: I                ! Loop index
            real(8)             :: L_pipe1          ! Length of pipe out of reactor (vertical)
            real(8)             :: L_pipe2          ! Length of pipe out of reactor (horizontal) to IHX
            real(8)             :: L_pipe3          ! Length of pipe out of IHX (vertical)
            real(8)             :: L_pipe4          ! Length of pipe out of IHX (horizontal) to reactor
            real(8)             :: D_pipe           ! Diameter of the pipes
            real(8)             :: K_elbow          ! Hydrualic resistance of an elbow
            real(8)             :: dP_pipe1         ! Pressure drop of pipe 1
            real(8)             :: dP_pipe2         ! Pressure drop of pipe 2
            real(8)             :: dP_pipe3         ! Pressure drop of pipe 3
            real(8)             :: dP_pipe4         ! Pressure drop of pipe 4
            real(8)             :: H_IHX            ! Height of IHX [m]
            real(8)             :: D_IHX            ! Diameter of IHX [m]
            integer             :: N_IHX            ! Number of nodes(CVs) IHX is split into
            real(8)             :: D_IHX_tube_out   ! Outside diameter of IHX tube [m]
            real(8)             :: t_IHX            ! Thickness of IHX tube [m]
            real(8)             :: D_IHX_tube_in    ! Inside diameter of IHX tube [m]
            real(8)             :: N_tubes_IHX      ! Number of tubes in IHX
            real(8)             :: A_IHX_P          ! Area on primary side of IHX [m^2]
            real(8)             :: D_IHX_P          ! Total diameter on primary side of IHX [m]
            real(8)             :: W_IHX_S          ! Mass flow rate on secondary side of IHX [kg/s]
            real(8)             :: A_IHX_S          ! Area on secondary side of IHX [m^2]
            real(8)             :: dP_IHX           ! Pressure drop of IHX        
            real(8),allocatable :: enthalpy_IHX(:)  ! Enthalpy of each node(CV) of IHX
            real(8)             :: dP_loop          ! Pressure drop of the loop
            
            real(8)             :: TC               ! Temperature of coolant in primary loop [Celcius]
 

            !================================================================================
            ! Initialize Variables
            !================================================================================
              
            ! PIPE DIMENSIONS -----------------------------------------------------
            L_pipe1=0.2
            L_pipe2=165.0e-2+650.0e-3/2.0
            L_pipe3=H_CORE+L_pipe1
            L_pipe4=L_pipe2
            D_pipe=16.0e-2
            K_elbow=0.8
        
            ! IHX -----------------------------------------------------------------
            H_IHX=1220.75                           ! Height of IHX [m]
            D_IHX=650.0e-3                          ! Diameter of IHX [m]
            N_IHX=10                                ! Number of nodes(CVs) IHX is split into
            D_IHX_tube_out=25.0e-3                  ! Outside diameter of IHX tube [m]
            t_IHX=2.5e-3                            ! Thickness of IHX tube [m]
            D_IHX_tube_in=D_IHX_tube_out-2.0*t_IHX  ! Inside diameter of IHX tube [m]
            N_tubes_IHX=222                         ! Number of tubes in IHX
            ! Primary Side IHX
            A_IHX_P=PI*D_IHX**2.0/4.0 - PI*D_IHX_tube_out**2.0/4.0*N_tubes_IHX*2.0
            D_IHX_P=4.0*A_IHX_P/(N_tubes_IHX*PI*D_IHX_tube_out*2.0+PI*D_IHX)
            ! Secondary Side IHX
            W_IHX_S=42.0
            A_IHX_S=PI*D_IHX_tube_in**2.0*4.0*N_tubes_IHX
        
        
        
        
            !================================================================================       
            ! Calculate pressure drops for pipes 1 and 2 before IHX at T_out of the core
            !================================================================================       
            dP_pipe1=dP_PIPE(W_core,D_pipe,T_out,L_pipe1,L_pipe1,K_elbow)
            dP_pipe2=dP_PIPE(W_core,D_pipe,T_out,L_pipe2,0.0_8,K_elbow)
        
        
            !================================================================================       
            ! Calculate pressure drop for IHX
            !================================================================================
            ! Initialize variables for loop
            dP_IHX=0.0
            ! Allocate arrays
            allocate(enthalpy_IHX(N_IHX))
            ! Loop through each IHX control volume node
            DO I=1,N_IHX
            	    IF (I==1) THEN
                        ! Y(NEIHX+IHXC)=HFLB(TCPP2)-POWER/10.0/WCORR
                        enthalpy_IHX(I)=flibe_enthalpy(T_out) - POWER/W_core
                        TC=flibe_temperature(enthalpy_IHX(I))
                    ELSE
	                    ! Y(NEIHX+IHXC)=Y(NEIHX+IHXC-1)-POWER/10.0/WCORR                    
                        enthalpy_IHX(I)=enthalpy_IHX(I-1) - POWER/W_core
                        ! TC = TFLB(Y(NEIHX+IHXC))
                        TC=flibe_temperature(enthalpy_IHX(I))
	                ENDIF
	                !TC=TFLB(Y(NEIHX+IHXC))

	                ! FRIC=FRICSL(WCORE,DIHXP,AIHXP,TC,1)
	                ! DPFIHX=DPFIHX+0.5D0*FRIC*(HIHX/NDIHX)*GIHXP**2/(DEN*DIHXP) ---> fix area vs diameter issue
	                ! DPGIHX=DPGIHX-DEN*GRAVITY*HIHX/NDIHX => L_vert = -L_vert
	                ! DPLIHX=0.0 -> K=0
                    dP_IHX=dP_IHX+dP_PIPE(W_core,D_IHX_P,TC,H_IHX/N_IHX,-H_IHX/N_IHX,0.0_8)
            ENDDO
            
            !================================================================================       
            ! Calculate pressure drops for pipes 3 and 4 after IHX at T_out of the IHX (TC)
            !================================================================================       
            dP_pipe3=dP_PIPE(W_core,D_pipe,TC,L_pipe3,L_pipe1,K_elbow)
            dP_pipe4=dP_PIPE(W_core,D_pipe,TC,L_pipe4,0.0_8,K_elbow)
                              
            !================================================================================       
            ! Calculate pressure drop for the loop
            !================================================================================
            dP_loop=0.5*(W_core/(PI*D_pipe**2.0/4.0))**2.0/flibe_rho(TC)*19.0   ! ---------need to check temp
            dP_loop=0.1149*W_core**3.0-7.6066*W_core**2.0+746.126*W_core        ! pressure drop for the primary loop except the core
                 
            !================================================================================
            ! Deallocate all arrays
            !================================================================================
            deallocate(enthalpy_IHX)
        
        END SUBROUTINE pipes_dP

        
    END MODULE pipes