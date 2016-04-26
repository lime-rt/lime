      PROGRAM GRID_INTERPOL
c
c                                            Version December 1 2008
c
c                            Simon Bruderer
c                            Institute of Astronomy 
c                            ETH Zurich
c                            8093 Zurich, Switzerland
c        
c                            simonbr-at-astro.phys.ethz.ch
c
c     This Fortran program does the interpolation of the
c     abundance from the grid of chemical models.
c
c     You can find a documentation of the program at the
c     homepage of the grid of chemical models at
c     http://www.astro.phys.ethz.ch/cgi-bin/chemgrid
c
c
c     Although this program has been thoroughly tested, the
c     authors do not claim that it is free of errors and
c     gives correct results in all situations.
c
c     Publications using this program should make a reference
c     to our paper: S.Bruderer, S.D. Doty and A.O. Benz
c
c
c     This code contains the following routines:
c
c       PROGRAM GRID_INTERPOL
c         Main routine, reads the inputfile and writes the output
c
c       REAL FUNCTION Interpol(Param,Time)
c         Does the multidimensional interpolation
c
c       SUBROUTINE Read_File(Molec,Time)
c         Reads the grid file with the molecular data in binary form
c
c
c     UV-Coordinate system (G0,tau) <-> (alpha,beta)
c
c       SUBROUTINE UVGRID_INV(G0,Tau,Alpha,Beta)
c         Calculates alpha and beta from a given G0,tau         
c
c       SUBROUTINE UVGRID(Alpha,Beta,G0,Tau)
c         Calculates G0, tau from alpha, beta
c
c       SUBROUTINE UVGRID_ESTI(G0,Tau,Alpha,Beta)
c         Calculates a first estimatio of alpha and beta
c         for G0 and tau. This is a routine is needed by UVGRID_INV.
c
c
c     X-ray ionization rate 
c
c       REAL FUNCTION XRATE(FLX,NH,TX)
c         Calculates the ionization rate by X-rays
c
c       REAL FUNCTION SIGMA_PH(EN)
c         Photoionization cross-section for a given photon energy
c
c       REAL FUNCTION SIGMA_CO(EN)
c         Compton ionzation cross-section for a given photon energy
c
c       REAL FUNCTION XRATE_INT(FLX,NH,TX)
c         Interpolation of the ionization rate from a lookup-table
c
c
c     The remainder of the code are data-tables:
c
c       BLOCK DATA UVGRIDDATA
c         Data of the (alpha,beta) - coordinate system
c
c       BLOCK DATA UV_ESTIMATE
c         Data used by the subroutine UVGRID_ESTI to obtain
c         a first guess of (G0,tau) --> (alpha,beta)
c
c       BLOCK DATA XIONDATA_1keV
c         Table of the X-ray ionization rate. The integral 
c         for the ionization rate is evaluated between 1 kev
c         and 100 keV
c
c       BLOCK DATA XIONDATA_100eV
c         Table of the X-ray ionization rate. The integral 
c         for the ionization rate is evaluated between 100 ev
c         and 100 keV
c
c
c
c     You'll find frequently find comments starting with 'cdebug'
c     in the code. These are meant for debuging and maintenance of
c     the program.
c
c     This version of the code works with the grid version 22.0,
c     published in the paper:
c
c     Parameters: 5 (Density, Temperature, Alpha, Beta, Ionization rate)
c     Models    : 113850
c     Header    : 'GRIDMOL_N22.0'
c
c





c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     Main program: Reads input file, Writes output files
c       

      IMPLICIT NONE       
  
c     Constants
c
c       Nr_Param: Number of dimensions of the grid
c       Nr_Mod:   Number of modelsn of the grid
c       Input_File: Name of the input-file
c       Output_File: Name of the output-file
      INTEGER Nr_param, Nr_mod
      PARAMETER(Nr_Param=5)
      PARAMETER(Nr_Mod=24480)

      CHARACTER*(*) Input_File
      PARAMETER(Input_File='grid_in.dat')
      CHARACTER*(*) Output_File
      PARAMETER(Output_File='grid_out.dat')
       
c     Common-Block 'GRID' : Stores the grid-data at two
c                           different chemical ages
c       Param_Range: Number of points along each grid dimension
c       Param_Space: Values of the grid-points along each grid dimension
c       Param_Name: Name of the parameter
c       Grid_A: Grid-data at time-step before the requested chemical age
c       Grid_B: Grid-data at time-step after the requested chemical age
c       Time_Low:  Time in s of GRID_A
c       Time_High: Time in s of GRID_B
      INTEGER Param_Range(Nr_Param),MaxPRange
      parameter(MaxPRange=40)
      REAL Param_Space(Nr_Param,MaxPRange)     
      CHARACTER*9 Param_Name(Nr_Param)
      REAL Grid_A(Nr_Mod)
      REAL Grid_B(Nr_Mod)
      REAL Time_Low, Time_High                    
      COMMON /Grid/ Grid_A, Grid_B, Param_Range, Param_Space,
     &  Param_Name, Time_Low, Time_High  
          
c     Functions
      REAL Interpol,XRATE,XRATE_INT     
      EXTERNAL Interpol,XRATE,XRATE_INT
            
c     Variables
c       Params: Parameters for the interpolation Parms(1)=Density,
c               Parms(2)=Temp etc
c       Model_Time: Chemical age in s
c       Model_Dens: Density in cm^-3
c       Model_T: Temperature in K
c       Model_Zeta: Cosmic ray ionization rate s^-1
c       Model_Alpha: Value of Alpha in the
c                   (Alpha,Beta) - coordinate system
c       Model_Beta: Value of Beta in the
c                   (Alpha,Beta) - coordinate system
c       Model_G0: FUV flux ISRF
c       Model_Tau: FUV attenuation A_V
c       Model_Fx: X-ray flux erg s^-1 cm^-2
c       Model_NH: Attenuating column density cm^-2
c       Model_Tx: Plasma temperature
c       Model_XION: X-ray inization rate
c       Check_G0: Check if the conversion
c                 (G0,Tau)->(Alpha,Beta) was accurate
c       Check_TAU: Check if the conversion
c                  (G0,Tau)->(Alpha,Beta) was accurate
c       Val_Interpol: Interpolated fractional abundance
c                     (-1 if something went wrong)
c       dummy: Read a line of the input-file that is later not
c              used in the code
c       Model_Molecule: Name of the Molecule
c       Nr_Points: Number of lines (parameters) in the input file
c       Model_DO: 0 if the grid should not do an interpolation for
c                 a particular line
      REAL Params(Nr_Param)
      REAL Model_Time, Model_Dens, Model_T, Model_Zeta
      REAL Model_Alpha, Model_Beta, Model_G0, Model_Tau
      REAL Model_Fx, Model_NH, Model_Tx, Model_XION
      REAL Check_G0,Check_TAU
      REAL Val_Interpol              
      CHARACTER*200 dummy
      CHARACTER*10 Model_Molecule       
      INTEGER Nr_Points, I, Model_DO
                             
cccccccccccccccccccccc            
c     Code

c     Say hello to the world...

c      WRITE(*,*) 'GRID:'
c      WRITE(*,*) 'GRID: Interpolation of chemical abundances'
c      WRITE(*,*) 'GRID:'
c      WRITE(*,*) 'GRID: S. Bruderer, A.O. Benz, S.D. Doty 2009'
c      WRITE(*,*) 'GRID:                    ApJS (...)'
c      WRITE(*,*) 'GRID:'
c      WRITE(*,*) 'GRID: This code works with datafiles of version'
c      WRITE(*,*) 'GRID:            >>GRIDMOL_R0.01<<'            
c      WRITE(*,*) 'GRID:'      
c      WRITE(*,*) 'GRID: Nov 25 2009: This is the interpolation'
c      WRITE(*,*) 'GRID: method for Ruud Vissers version of the code'            
c      WRITE(*,*) 'GRID:'            
      
c     open input file and read header
                 
      CALL SYSTEM('rm -f '//Output_File)
      OPEN(20,File=Input_File,Status='OLD',Form='FORMATTED')         
      READ(20,'(A)') dummy
      READ(20,'(A)') dummy
c      WRITE(*,*) 'GRID: Source: '//dummy(11:31)
      READ(20,'(A)') Model_Molecule      
c      WRITE(*,*) 'GRID: Species: '//
c     &           Model_Molecule(1:10)
      READ(20,'(A)') dummy
      READ(20,*) Model_Time              
c      WRITE(*,'(A,1PE9.3)') ' GRID: Time (s): ',Model_Time
c      WRITE(*,*) 'GRID:'      
      READ(20,'(A)') dummy
      READ(20,*) Nr_Points
      READ(20,'(A)') dummy
      
c     open output file and write header      
             
      OPEN(30,File=Output_File,Status='NEW',Form='FORMATTED')         
c      WRITE(30,'(A)') '! GRID-OUTPUTFILE'
c      WRITE(30,'(A)') '! MOLECULE'
c      WRITE(30,'(A)') Model_Molecule      
c      WRITE(30,'(A)') '! TIME'
c      WRITE(30,'(1PE12.3)') Model_Time              
c      WRITE(30,'(A)') '! NR OF POINTS'
c      WRITE(30,'(I5)') Nr_Points      
c      WRITE(30,'(A)') '! Abundance    N(H_TOT)      T         G0     '//
c     & '     Tau         Zeta        Fx          NH          Tx     '//
c     & '     XION'
       
c     read molecular data from the binary file       
       
      CALL read_file(Model_Molecule,Model_Time)
                              
c      WRITE(*,*) 'GRID: Number of Points:',Nr_Mod	           

       
c     Do the interpolation for all lines in the input file      
      
      DO 100 I=1,Nr_Points

cdebug If you want to see if the code died...       
c        IF(MOD(I+1,10000).EQ.0) THEN
c	  WRITE(*,*) 'GRID: Interpolated ',I+1,' points...'
c	END IF
                     
		     
c     Read input parameters
		     
        READ(20,*) Model_DO, Model_Dens, Model_T,
     &    Model_G0,Model_Tau,Model_Zeta
c     &    Model_Fx,Model_NH,Model_Tx
         Model_Fx=0.0
         Model_NH=0.0
         Model_Tx=0.0            
         Model_XION=0.0

c     Do not calculate interpolation for these input parameters
          
 	IF(Model_DO.EQ.0) THEN
	  Val_Interpol=-9.0
	  Model_Alpha=0.0
	  Model_Beta=0.0
	  Check_G0=0.0
	  Check_Tau=0.0
	  Model_XION=0.0
   	  GOTO 300	 
	END IF
	
	
c     Calculate the ionization rate by X-rays: Use the lookup table
c     for plasma temperatures above 3e6 K	

csb251109: switched off for ruuds grid
c	IF(MODEL_TX.GT.3.0E6) THEN
c          Model_XION=XRATE_INT(Model_Fx,Model_NH,Model_Tx)
c        ELSE
c          Model_XION=XRATE(Model_Fx,Model_NH,Model_Tx)
c        END IF

cdebug If you dont want to use the lookup table in XRATE_INT, you may
c      uncomment the next line. This makes sense if you assume low p
c      plasma temperatures, the interpolation is however significantly
c      slower, since an integral has to be calculated for each interpolation
c       Model_XION=XRATE(Model_Fx,Model_NH,Model_Tx)



c     Calculate the values of (Alpha,Beta) from the input values of
c     (G0,Tau) FUV-flux/attenuation

	 	 	
csb251109 <---------------------------From here
c           switched off for ruuds grid----------------------to here>								 
c        IF(Model_Tau.GT.10.0) Model_Tau=30.0
c	IF(Model_Tau.LT.0.1) Model_Tau=0.1
c	IF(Model_G0.LE.1.0E-4) Then
c	  Model_G0=1.0E-4
c	END IF
	                 	 
c	CALL UVGRID_INV(ALOG10(Model_G0),ALOG10(Model_Tau),
c     &	  Model_Alpha,Model_Beta)

c	CALL UVGRID(Model_Alpha,Model_Beta,Check_G0,Check_Tau)	
	 	 
c     switch to the model without FUV
c        IF(Model_Alpha.GT.6.1) THEN
c	  Model_Alpha=10.0
c	  Model_Beta=0.0	 
c	END IF
	
c      correct rounding problems at the edge of the grid
c        IF((Model_Alpha.LE.-0.5).AND.(Model_Alpha.GT.-0.501)) 
c     &	  Model_Alpha=-0.5
c        IF((Model_Beta.LE.0.0).AND.(Model_Beta.GT.-0.001)) 
c     &	  Model_Beta=0.0	 
c        IF((Model_Beta.GT.1.0).AND.(Model_Beta.LE.1.001)) 
c     &	  Model_Beta=1.0

cdebug Print values of G0/Tau, Alpha/Beta
c	WRITE(*,*) 'Model_G0,Model_Tau=',Model_G0,Model_Tau
c	WRITE(*,*) 'Model_Alpha,Model_Beta=',Model_Alpha,Model_Beta
c	WRITE(*,*) 'Check_G0,Check_Tau=',10.0**Check_G0,10.0**Check_Tau
csb251109 <---------------------------to here
c           switched off for ruuds grid----------------------to here>

	 	 	        
c     The grid uses N(H2) as input for historical reasons (Steves Code), 
c     thus divide the total density by two				
        Params(1)=Model_Dens  
        Params(2)=Model_T ! T
csb251109: adpted for ruuds grid
	Params(3)=Model_G0
	Params(4)=Model_TAU
	Params(5)=Model_Zeta



c	Params(3)=10.0**Model_Alpha
c	Params(4)=10.0**Model_Beta
c        Params(5)=Model_Zeta+Model_XION ! Zeta

cdebug To check if the grid interpolation works fine
c      Using these lines, the input paramters Model_G0, Model_Tau correspond to (Alpha,Beta)	
c       Params(3)=10.0**Model_G0
c	Params(4)=10.0**Model_TAU
c	Params(5)=Model_Zeta
 	 	  
c        Val_Interpol=1.0   
        Val_Interpol=Interpol(Params,Model_Time)

	 
300     CONTINUE


c     Write out the result along with the physical parameters and 
c     values to check the conversion of the (G0/Tau) --> (Alpha,Beta)
c     coordinate system.
	 	 	

csb251109: adpted for ruuds grid
        WRITE(30,201) Val_Interpol
c        WRITE(30,200) Val_Interpol,Model_Dens,Model_T,
c     &    Model_G0,
c     &    Model_Tau,
c     &    Model_Zeta,Model_Fx,Model_NH,
c     &    Model_Tx,Model_XION

200     FORMAT(1PE12.3,1PE12.3,0PF12.3,1PE12.3,
     &    1PE12.3,1PE12.3,1PE12.3,
     &    1PE12.3,1PE12.3,1PE12.3)
201     FORMAT(1PE12.3)


       
c        WRITE(30,200) Val_Interpol,Model_Dens,Model_T,
c     &    Model_G0,Model_Alpha,10**Check_G0,
c     &    Model_Tau,Model_Beta,10**Check_Tau,
c     &    Model_Zeta,Model_Fx,Model_NH,
c     &    Model_Tx,Model_XION

c200     FORMAT(1PE12.3,1PE12.3,0PF12.3,1PE12.3,0PF12.3,1PE12.3,
c     &    0PF12.3,0PF12.3,0PF12.3,1PE12.3,1PE12.3,
c     &    1PE12.3,1PE12.3,1PE12.3)

       
100   CONTINUE 

c     Bye...
      
c      WRITE(*,*) 'GRID: Number of interpolations done:',Nr_Points
      CLOSE(20)
      ClOSE(30)       
      
c     Say goodbye to the world...      

c      WRITE(*,*) 'GRID:'
c      WRITE(*,*) 'GRID: Finished!'
c      WRITE(*,*) 'GRID:'      
      

                  
      END




c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c REAL FUNCTION Interpol(Param,Time)
c
c     Interpolates the chemical abundances at the the 'time' with 
c     physical parameters (NH, T...) in Param.
c     (e.g. an array REAL Param(4) with 4 elements corresponding to
c     the density, temperature, G0 and Tau)
c
c     The function returns the fractional abundance or -1 if an error
c     occured (e.g. parameters are outside the grid range).
c
      REAL FUNCTION Interpol(Param,Time)
             
      IMPLICIT NONE
             	     
c     Constants
c       Nr_Param: Number of dimensions of the grid
c       Nr_Mod:   Number of modelsn of the grid
      INTEGER Nr_Param, Nr_Mod
      PARAMETER(Nr_Param=5)
      PARAMETER(Nr_Mod=24480)
          	         
c     Common-Block 'GRID' : Stores the grid-data at two
c                           different chemical ages
c       Param_Range: Number of points along each grid dimension
c       Param_Space: Values of the grid-points along each grid dimension
c       Param_Name: Name of the parameter
c       Grid_A: Grid-data at time-step before the requested chemical age
c       Grid_B: Grid-data at time-step after the requested chemical age
c       Time_Low:  Time in s of GRID_A
c       Time_High: Time in s of GRID_B
      INTEGER Param_Range(Nr_Param),MaxPRange
      parameter(MaxPRange=40)
      REAL Param_Space(Nr_Param,MaxPRange)     
      CHARACTER*9 Param_Name(Nr_Param)
      REAL Grid_A(Nr_Mod)
      REAL Grid_B(Nr_Mod)
      REAL Time_Low, Time_High                    
      COMMON /Grid/ Grid_A, Grid_B, Param_Range, Param_Space,
     &  Param_Name, Time_Low, Time_High
     
c     Parameters
c       Param : Values of the parameters for the interpolation
c       Time : Chemical age (must be between Time_Low and Time_high)
      REAL Param(Nr_Param)
      REAL Time
      
c     Variables
c       Weight: Weight of one grid point along one axis
c       Weight_Sum: Total weight of one grid point along all axis
c       LogSum_A: Weighted sum of the logarithm of the abundances
c                 at the time Time_low
c       LogSum_B: Weighted sum of the logarithm of the abundances
c                 at the time Time_high
c       Index_Param: Position in the grid along one axis
c       Param_Edge: Calculate the index for the model to be summed up
c       Calc_Index: Lookup-table to calculate Index_Grid
c       Index_Grid: The index of the model to be summed up
      REAL Weight(Nr_Param,2)       
      REAL Weight_Sum, LogSum_A,LogSum_B
      INTEGER Index_Param(Nr_Param), Param_Edge(Nr_Param)
      INTEGER Calc_Index(Nr_Param), Index_Grid, HasNeg
     
      INTEGER I, J, K, L, M
      REAL A, B, C, D, E
       

cccccccccccccccccccccc       
c     Code

      Interpol=-10.0
                                     	      	  			      
c     Time is not within the Time_Low and Time_High as read
c     from the database: Stop - should however not happen if 
c     this function is called from the standard code

      IF((Time.GT.Time_High).OR.(Time.LT.Time_Low)) Then
        WRITE(*,'(A)') 'GRID-ERR: Time outside interval'
        WRITE(30,'(1PE12.3)') -1.0
        CLOSE(20)
        CLOSE(30)
	STOP
      END IF	   

      
c     Find the indices, ie. the values on the grid between those
c     the given set of parameter lies.

      DO 100 I=1,Nr_Param          	  
        Index_Param(I)=1
csb251109 adapted for ruuds grid
        if(param_range(i).gt.1) then
          Index_Param(I)=-1
  	  DO 200 J=1,Param_Range(I)-1
cdebug Write out the grid-range
c           WRITE(*,*) 'I,J,Param_Space(I,J)=',I,J,Param_Space(I,J)
c           WRITE(*,'(1PE10.3)') Param_Space(I,J)	             
	    IF((Param(I).GE.Param_Space(I,J)).AND.
     &	      (Param(I).LE.Param_Space(I,J+1))) THEN
              Index_Param(I)=J
            END IF     
200         CONTINUE	          
	    IF(Index_Param(I).EQ.-1) THEN
	      WRITE(*,'(A,A,A,1PE10.3,A)') 'GRID-ERR: Parameter ',
     &	        Param_Name(I),':',Param(I),' outside range'
              RETURN
            END IF  	         
         end if
100   CONTINUE

c     Do the interpolation in logarithmic space, give the 
c     weight 0 to the grid-point on the lower side of an
c     interval in the parameter space (... and 1 to the 
c     upper side)

      DO 600 I=1,Nr_Param   

	Weight(I,1)=1.0
        Weight(I,2)=0.0

csb251109 adapted for ruuds grid
        if(param_range(i).gt.1) then  


c     Log_10 of the lower/upper bound in the interval of the
c     parameter(i) and Log_10 of the parameter

          A=ALOG10(Param_Space(I,Index_Param(I)))
          B=ALOG10(Param_Space(I,Index_Param(I)+1))
          D=ALOG10(Param(I))	   
	
c     Find the position within the interval

          C=(D-A)/(B-A)   	   
          Weight(I,1)=1.0-C 
          Weight(I,2)=C

c Option: use linear interpolation
c         A=Param_Space(I,Index_Param(I))
c         B=Param_Space(I,Index_Param(I)+1)
c         D=Param(I)	   
c         C=(D-A)/(B-A)   	   
c         Weight(I,1)=1.0-C 
c         Weight(I,2)=C	   
	   
	   	   
cdebug Write out the neighbouring points along all grid dimensions
c         WRITE(30,'(I4,A,A15,F10.3,1PE10.3,1PE10.3)') 
c     &     I,' ',Param_Name(I),C,
c     &     Param_Space(I,Index_Param(I)),
c     &     Param_Space(I,Index_Param(I)+1)

        end if

600   CONTINUE       


c     The index in Grid_A(1..Nr_Mod) is obtained using the fact, that
c     the last index runs fastest in the database:

      J=Nr_Mod
      DO 800 I=1,Nr_Param
        J=J/Param_Range(I)
        Calc_Index(I)=J       
800   CONTINUE       

c     Calculate the indices for the interpolation routine: This is
c     illustrated at the example of a interpolation in two dimensions:
c
c     A linear 2d-Interpolation V(x,y) can be written down as
c    
c     V(X,Y)=V_00*(1-x)*(1-y)+V_01*(1-x)*y+V_10*x*(1-y)+V_11*x*y  
c
c     where V_ab are the grid-points. The indices can thus
c     easily be obtained by a conversation of a decimal number 
c     0..2**(Nr_Params-1) into the base 2.
c
c     The array Param_Edge(J=1..Nr_Param) gives the representation 
c     in the base 2.

      LogSum_A=0.0
      LogSum_B=0.0       
      K=2**Nr_Param       
       
      HasNeg=0
      DO 400 I=1,K         
        M=I-1        
        DO 500 J=Nr_Param,1,-1    	  
	  Param_Edge(J)=1
	  L=2**(J-1)	    	     
	  IF(M/L.GT.0) Param_Edge(J)=2
	  M=Mod(M,L)

          if((param_range(j).eq.1).and.
     &       (param_edge(j).eq.2)) goto 400
	     	     	   	     
500     CONTINUE
              
        Weight_Sum=1.0
	Index_Grid=1
		  
        DO 700 J=1,Nr_Param
 	  Weight_Sum=Weight_Sum*Weight(J,Param_Edge(J))
          M=Index_Param(J)+Param_Edge(J)-1    	     
	  Index_Grid=Index_Grid+(M-1)*Calc_Index(J)
700     CONTINUE

 
        IF(Grid_A(Index_Grid).LT.0.0.OR.
     &       Grid_B(Index_Grid).LT.0.0) THEN
          WRITE(*,'(A,A)') 'GRID-ERR: Negative abundance ',
     &      'in one of the grid points'
          WRITE(*,*) index_grid,Grid_A(Index_Grid),Grid_B(Index_Grid)
          write(47,*) index_grid
          HasNeg=1
c          RETURN
        END IF  	         

c option: set a minmum for the minimum fractional abundance 	  
        A=MAX(1.0E-30,Grid_A(Index_Grid))
        B=MAX(1.0E-30,Grid_B(Index_Grid))

c        A=Grid_A(Index_Grid)
c        B=Grid_B(Index_Grid)

cdebug Prints out all summed values
c       WRITE(*,*) Param_Edge
c       WRITE(*,'(1PE10.3,I8)') Weight_Sum, Index_Grid
c       WRITE(*,'(1PE15.3,A,1PE15.3)') A, '    ',B
	  	
        LogSum_A=LogSum_A+Weight_Sum*ALOG10(A)
        LogSum_B=LogSum_B+Weight_Sum*ALOG10(B)

c        end if
         
400   CONTINUE       

      IF(HasNeg.EQ.1) RETURN

c     Interpolation in time also in the Log-Log space

      A=ALOG10(Max(Time_Low,1.0E-30))
      B=ALOG10(Time_High)       
      D=ALOG10(Max(Time,1.0E-30))
      C=(D-A)/(B-A)
      E=(1.0-C)*LogSum_A+C*LogSum_B
      Interpol=10**E
              
cdebug Write out a summary of the input/output data of this subroutine
c      DO 900 I=1,Nr_Param
c        WRITE(*,'(A9,A1,1PE10.3)') Param_Name(I),':',Param(I)       
c 900  CONTINUE       
c      WRITE(*,'(A9,A1,1PE10.3)') 'Time     ',':',Time
c      WRITE(*,'(1PE10.3,A,1PE10.3,A,1PE18.6)') 10.0**LogSum_A,
c     &  ' ',10.0**LogSum_B,' ',10**E
      
      RETURN
                    
      END





c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c SUBROUTINE Read_File(Molec,Time)
c
c     Reads the data from a molecular data file of the grid
c     (binary format)
c
c
      SUBROUTINE Read_File(Molec,Time)
       
      IMPLICIT NONE
    
c     Constants
c       Nr_Param: Number of dimensions of the grid
c       Nr_Mod:   Number of modelsn of the grid
      INTEGER Nr_Param, Nr_Mod
      PARAMETER(Nr_Param=5)
      PARAMETER(Nr_Mod=24480)
              	      
c     Common-Block 'GRID' : Stores the grid-data at two
c                           different chemical ages
c       Param_Range: Number of points along each grid dimension
c       Param_Space: Values of the grid-points along each grid dimension
c       Param_Name: Name of the parameter
c       Grid_A: Grid-data at time-step before the requested chemical age
c       Grid_B: Grid-data at time-step after the requested chemical age
c       Time_Low:  Time in s of GRID_A
c       Time_High: Time in s of GRID_B
      INTEGER Param_Range(Nr_Param),MaxPRange
      parameter(MaxPRange=40)
      REAL Param_Space(Nr_Param,MaxPRange)     
      CHARACTER*9 Param_Name(Nr_Param)
      REAL Grid_A(Nr_Mod)
      REAL Grid_B(Nr_Mod)
      REAL Time_Low, Time_High                    
      COMMON /Grid/ Grid_A, Grid_B, Param_Range, Param_Space,
     &  Param_Name, Time_Low, Time_High   

c     Parameters
c       Molec: Name of the molecule
c       Time: Chemical age
      CHARACTER*(*) Molec       
      REAL Time     
     
c     Variables
c       Read_TimeAbu: Array to read in one data-line of the grid
c       IoOK: IO-Status
c       Read_Nr_Param: Number of parameters to be read in
c       Nr_Time: Nr of timesteps in the grid
c       Index_Time: In what index of the timesteps is the requested
c       Check_Mod: Check the number of models
c       Filename: Filename of the binary molecule data
c       Mol: Name of the molecule
c       Header: Read one line of the header
c       Molname: Read the name of the molecule from the file
      REAL Read_TimeAbu(400)
      INTEGER IoOk, Read_Nr_Param
      INTEGER Nr_Time, Index_Time, Check_Mod      
      INTEGER I, J
      CHARACTER*100 Filename, Mol       
      CHARACTER*80 Header
      CHARACTER*9 Molname       
          
cccccccccccccccccccccc              
c     Code

c     Put the filename together and open the binary file

      Mol=molec
      Filename='mol_grid_' // Mol(1:Index(Mol,' ' )-1) // '.bin'
      OPEN(10,File=Filename,Status='OLD',Form='UNFORMATTED',
     &  Iostat=IoOk, Err=300)
      IF(IoOk.NE.0) GOTO 300

c     Read header and check for the correct version of the grid-data

      READ(10,ERR=350) Header  
      IF(Header.NE.'GRIDMOL_R0.01') THEN
        WRITE(*,'(A)')  'GRID-ERR: Wrong grid-version'
        WRITE(*,'(A)')  '.'//Header//'.'
        WRITE(30,'(1PE12.3)') -2.0
        CLOSE(20)
        CLOSE(30)
	STOP
      END IF

c     Read all physical parameters in
       
      READ(10,ERR=350) Molname
      READ(10,ERR=350) Read_Nr_Param
      IF(Read_Nr_Param.NE.Nr_Param) THEN
        WRITE(*,'(A)')  'GRID-ERR: Number of parameters incorrect'
        WRITE(30,'(1PE12.3)') -3.0
        CLOSE(20)
        CLOSE(30)
	STOP
      END IF
      Check_Mod=1
      DO 100 I=1, Nr_Param
        READ(10,ERR=350) Param_Name(I)
	READ(10,ERR=350) Param_Range(I)
        IF(Param_Range(I).gt.MaxPRange) then
          WRITE(*,'(A)') ' GRID-ERR:'
          WRITE(*,'(A)') ' GRID-ERR: Parameter range too large'
          WRITE(*,'(2A)') ' GRID-ERR:    for parameter ',Param_Name(I)
          WRITE(*,'(A)') ' GRID-ERR:'
          WRITE(30,'(1PE12.3)') -4.0
          CLOSE(20)
          CLOSE(30)
          STOP 2241
        ENDIF
        Check_Mod=Check_Mod*Param_Range(I)
	DO 200 J=1, Param_Range(I)
	  READ(10,ERR=350) Param_Space(I,J)
cdebug Write out the parameter space of the grid-file
c         WRITE(*,*) 'I,J,Param_Space(I,J)=',I,J,Param_Space(I,J)	    
200     CONTINUE	            
100   CONTINUE   

c     Check if the number of models agrees to the parameter Nr_Mod,
c     used for the allocation of the arrays.

      IF(Check_Mod.NE.Nr_Mod) THEN
        WRITE(*,'(A)') 'GRID-ERR: Wrong number of models'
        WRITE(*,*) Check_Mod,Nr_Mod
        WRITE(30,'(1PE12.3)') -5.0
        CLOSE(20)
        CLOSE(30)
	STOP
      END IF

c     Read in the time steps
 
      READ(10,ERR=350) Nr_Time                
      DO 500 I=1,Nr_Time
        READ(10,ERR=350) Read_TimeAbu(I)             
500   CONTINUE      



c     Find out what time-steps should be read out of the grid
      Index_Time=-1
      DO 600 I=1,Nr_Time-1
        IF(Time.GE.(Read_TimeAbu(I)) .AND.
     &    (Time.LE.Read_TimeAbu(I+1))) THEN
          Index_Time=I
        END IF	         
600   CONTINUE


c     Stop if it's outside the grid

      IF(Index_Time.EQ.-1) THEN
        WRITE(*,'(A)')  'GRID-ERR: Time outside Grid'
        WRITE(30,'(1PE12.3)') -6.0
        CLOSE(20)
        CLOSE(30)
	STOP
      END IF
      
      Time_Low=Read_TimeAbu(Index_Time)
      Time_High=Read_TimeAbu(Index_Time+1)
      
c     Read in the grid-data
      
      DO 700 I=1,Nr_Mod
        READ(10,ERR=350) (Read_TimeAbu(J),J=1,Nr_Time)
	Grid_A(I)=Read_TimeAbu(Index_Time)	
	Grid_B(I)=Read_TimeAbu(Index_Time+1) 
cdebug Write out the data just obtained from the data file
c	WRITE(*,'(I4,1PE10.3,1PE10.3)') I,Grid_A(I),Grid_B(I)           
700   CONTINUE    

      CLOSE(10)
                                   
      RETURN

300   WRITE(*,'(A)') 'GRID-ERR: opening: '//Filename       
      WRITE(30,'(1PE12.3)') -7.0
      CLOSE(20)
      CLOSE(30)
      STOP
       
350   CLOSE(10)
      WRITE(*,'(A)')  'GRID-ERR: reading: '//Filename
      WRITE(30,'(1PE12.3)') -8.0
      CLOSE(20)
      CLOSE(30)
      STOP
     
      END
       




c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c SUBROUTINE UVGRID_INV(G0,Tau,Alpha,Beta)
c
c     This subroutine does the coordinate transformation from
c     (G0,Tau) --> (Alpha,Beta) as it is used to read out the 
c     correct value of the grid.
c
c     Scheme of the coordinate system (see. Fig. 2 in the paper)
c
c     Beta (Tau)
c     |     ----1.8750-----
c     |     |             |
c     |     8           -3
c     |     |             |
c     |     ---(-1.125)----
c     | 
c     --------------- Alpha (G0)
c
c     This subroutine uses a gradient descent method. The initial values
c     are taken from a table using the subroutine UVGRID_ESTI. Using this
c     accurate initial guess, only very few steps are required for a good
c     solution. 
c
c
c
      SUBROUTINE UVGRID_INV(G0,Tau,Alpha,Beta)
      
      IMPLICIT NONE
        
c     Common-Block 'UV_DATA' : Coordinate system for the coordinate
c                              transformation (G0,Tau) --> (Alpha,Beta)
c       N_Alpha: Number of points in direction of Alpha
c       N_Beta: Nubmer of points in direction of Beta
c       UV_TABLE: Table containing with lines containing
c                 (log(Alpha),log(Beta),log(g0),log(Tau))
c       Alpha_G: Points of the UV_TABLE in log(Alpha)-direction
c       Beta_G: Points of the UV_TABLE in log(Beta)-direction
      INTEGER N_Alpha,N_Beta
      PARAMETER(N_Alpha=27)
      PARAMETER(N_Beta=28)      
      REAL UV_TABLE(N_Alpha*N_Beta,4),Alpha_G(N_Alpha),Beta_G(N_Beta)
      COMMON/UV_DATA/ UV_TABLE,Alpha_G,Beta_G
            
c     Parameters   
c       G0: FUV Flux in ISRF (Input)
c       Tau: FUV Attenuation in A_V (Input)
c       Alpha: Alpha (Output)
c       Beta: Beta (Output)
      REAL G0, Tau, Alpha, Beta

c     Variables
c       A_UV, B_UV: Find the next point of the tabulated coordinate system
c       G0_UV, Tau_UV: Find the next point of the tabulated coordinate system
c       Dist: Distance in the (G0,Tau)-space
c       Min_Dist: Minimum distance to next point of the tabulated coordinate system
c       A_Guess, B_Guess: First guess of (Alpha,Beta) obtained from UVGRID_ESTI
c       Norm: Normalization of the descent-vector for the descent
c       EPS_A, EPS_B: Difference for the calculation of the gradient
c       Delta, Gamma: Direction of the gradient
c       G0_1, G0_2, G0_3: Derivation in direction of G0
c       Tau_1, Tau_2, Tau_3: Derivation in direction of G0
c       G0_Check, Tau_Check: G0 and TAU obtained from the calculated (Alpha,Beta)
c       Ratio_G0, Ratio_Tau: Ratio between G0_Check and input G0 and Tau_Check / Tau
c       STEPS: Number of steps of the gradient descent
      REAL A_UV,B_UV,G0_UV,Tau_UV,Dist,Min_Dist
      REAL A_Guess,B_Guess,Norm
      REAL EPS_A, EPS_B, Delta, Gamma
      REAL G0_1,G0_2,G0_3,Tau_1,Tau_2,Tau_3
      REAL Ratio_G0,Ratio_Tau,G0_Check,Tau_Check
      INTEGER STEPS, I
      

cccccccccccccccccccccc       
c     Code

c     This is certainly outside the data-table

      IF(ABS(G0).GT.40.OR.ABS(Tau).GT.15) THEN
        WRITE(*,'(A)') 'ERROR IN UVGRID_INV: STOP'
	STOP
      END IF

c     Get an initial guess

      A_Guess=-100.0
      B_Guess=-100.0
      CALL UVGRID_ESTI(G0,Tau,A_Guess,B_Guess)
      
c     If the table did not obtain a result (outside the range)
c     we just use the next point of UV_TABLE as initial point    
                  
      IF(A_Guess.LT.-99.0.OR.B_Guess.LT.-99.0) THEN

cdebug
c	WRITE(*,*) 'WARNING: UVGRID_ESTI out of range'

        Min_Dist=1E10
      
        DO 100 I=1,N_Alpha*N_Beta
          A_UV=UV_TABLE(I,1)
          B_UV=UV_TABLE(I,2)
          G0_UV=UV_TABLE(I,3)
          Tau_UV=UV_TABLE(I,4)
          Dist=SQRT((G0-G0_UV)**2+(Tau-Tau_UV)**2)
          IF(Dist.LE.Min_Dist) THEN
	    A_Guess=A_UV
	    B_Guess=B_UV
	    Min_Dist=Dist
	  END IF      
100     CONTINUE      
      END IF

C     Start gradient descent

cdebug
c     WRITE(*,*) 'A_Guess,B_Guess=',A_Guess,B_Guess

      STEPS=0

200   CONTINUE

c     Calculate derivatives in direction of G0 and Tau

      EPS_A=1E-2
      EPS_B=1E-2
      CALL UVGRID(A_Guess,B_Guess,G0_1,Tau_1)
      CALL UVGRID(A_Guess+EPS_A,B_Guess,G0_2,Tau_2)
      CALL UVGRID(A_Guess,B_Guess+EPS_B,G0_3,Tau_3)

c     Stop code, if we land outside the range of the table
      IF(G0_1.LT.-1E4.OR.G0_2.LT.-1E4.OR.G0_3.LT.-1E4)
     &  GOTO 300
      
c     Calculate the gradient

      Norm=SQRT((G0_2-G0_1)**2+(Tau_2-Tau_1)**2)
      G0_2=(G0_2-G0_1)/Norm
      Tau_2=(Tau_2-Tau_1)/Norm      

      Norm=SQRT((G0_3-G0_1)**2+(Tau_3-Tau_1)**2)
      G0_3=(G0_3-G0_1)/Norm
      Tau_3=(Tau_3-Tau_1)/Norm      

      Gamma=-(-G0_2*Tau+G0_2*Tau_1+Tau_2*G0-Tau_2*G0_1)/
     &  (G0_2*Tau_3-Tau_2*G0_3)     

      Delta=(G0*Tau_3-G0_1*Tau_3-G0_3*Tau+G0_3*Tau_1)/
     &  (G0_2*Tau_3-Tau_2*G0_3)
 
cdebug          
c     WRITE(*,*) 'Delta,Gamma,Norm=',Delta,Gamma,Norm

c     Next step. Use a somewhat "damped" descent to ensure 
c     convergence. Due to the usually very good initial guesses
c     from the table, the algorithm still needs only a few steps
c     until it converges.

      A_Guess=A_Guess+Delta*0.1
      B_Guess=B_Guess+Gamma*0.1
 
cdebug      
c     WRITE(*,*) 'A_Guess,B_Guess=',A_Guess,B_Guess

      STEPS=STEPS+1

c     Did it converge? Check for the ratio between the input values
c     of (G0,Tau) and 

      CALL UVGRID(A_Guess,B_Guess,G0_Check,Tau_Check)
      Ratio_G0=10.0**G0/10.0**G0_Check
      Ratio_Tau=10.0**Tau/10.0**Tau_Check
      Ratio_G0=MAX(Ratio_G0,1.0/Ratio_G0)
      Ratio_Tau=MAX(Ratio_Tau,1.0/Ratio_Tau)

      IF(Ratio_G0.GE.1.001.OR.Ratio_Tau.GE.1.001.AND.STEPS.LE.10000) 
     &  GOTO 200

300   CONTINUE
 
      IF(STEPS.GT.9990) THEN
        WRITE(*,'(A)') 'UV_INV did not converge!'         
        STOP
      END IF
 
cdebug
c      WRITE(*,*) 'I=',I,RATIO_G0,RATIO_Tau,A_Guess,B_Guess
            
      A_Guess=MIN(A_Guess,Alpha_G(N_Alpha))
      A_Guess=MAX(A_Guess,Alpha_G(1))
      B_Guess=MIN(B_Guess,Beta_G(N_Beta))
      B_Guess=MAX(B_Guess,Beta_G(1))

      Alpha=A_Guess
      Beta=B_Guess
      
      RETURN

      END





c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c SUBROUTINE UVGRID(Alpha,Beta,G0,Tau)
c
c
c     Coordinate transformation (Alpha,Beta) --> (G0,TAU): Since the
c     table in 'UV_DATA' contains a list for the (Alpha,Beta) values
c     this can be done by a simple interpolation in that table.
c
c     Returns (G0,TAU) = (-1.0e5,-1.0e5) if (Alpha,Beta) is
c     outside to grid.
c
      SUBROUTINE UVGRID(Alpha,Beta,G0,Tau)

      IMPLICIT NONE

c     Common-Block 'UV_DATA' : Coordinate system for the coordinate
c                              transformation (G0,Tau) --> (Alpha,Beta)
c       N_Alpha: Number of points in direction of Alpha
c       N_Beta: Nubmer of points in direction of Beta
c       UV_TABLE: Table containing with lines containing
c                 (log(Alpha),log(Beta),log(g0),log(Tau))
c       Alpha_G: Points of the UV_TABLE in log(Alpha)-direction
c       Beta_G: Points of the UV_TABLE in log(Beta)-direction
      INTEGER N_Alpha,N_Beta
      PARAMETER(N_Alpha=27)
      PARAMETER(N_Beta=28)
      REAL UV_TABLE(N_Alpha*N_Beta,4),Alpha_G(N_Alpha),Beta_G(N_Beta)
      COMMON/UV_DATA/ UV_TABLE,Alpha_G,Beta_G
 
c     Parameters   
c       Alpha: Alpha (Input)
c       Beta: Beta (Input)
c       G0: FUV Flux in ISRF (Output)
c       Tau: FUV Attenuation in A_V (Output)
      REAL Alpha, Beta, G0, Tau

c     Variables
c       G0_DD,G0_DU,G0_UD,G0_UU: G0 values at the 4 surrounding points
c                                of the (Alpha/Beta)-table
c       Tau_DD,Tau_DU,Tau_UD,Tau_UU: Tau values at the 4 surrounding
c                                    points of the (Alpha/Beta)-table
c       D_Alpha,D_Beta: Position in the intervals of the tables
c       A_PT,B_PT: Index in the intervals of the tables
      REAL G0_DD,G0_DU,G0_UD,G0_UU
      REAL Tau_DD,Tau_DU,Tau_UD,Tau_UU
      REAL D_Alpha,D_Beta
      INTEGER A_PT,B_PT,I

cccccccccccccccccccccc       
c     Code


c     Find the index within the table for both directions Alpha and Beta

      A_PT=-100
      B_PT=-100
      DO 100 I=1,N_Alpha-1
        IF(Alpha.GE.Alpha_G(I).AND.Alpha.LE.Alpha_G(I+1)) A_PT=I
100   CONTINUE
      DO 200 I=1,N_Beta-1
        IF(Beta.GE.Beta_G(I).AND.Beta.LE.Beta_G(I+1)) B_PT=I
200   CONTINUE

      IF(A_PT.EQ.-100.OR.B_PT.EQ.-100) THEN 
cdebug
c       WRITE(*,*) 'OUTSIDE UV-GRID'
        G0=-1.E5
        Tau=-1.E5
        RETURN
      END IF
      
c     Read the values of G0 and Tau from the data: Have to calculate the 
c     index.

      G0_DD=UV_TABLE((A_PT-1)*N_Beta+B_PT,3)
      G0_DU=UV_TABLE((A_PT-1)*N_Beta+(B_PT+1),3)
      G0_UD=UV_TABLE(A_PT*N_Beta+B_PT,3)
      G0_UU=UV_TABLE(A_PT*N_Beta+(B_PT+1),3)
      Tau_DD=UV_TABLE((A_PT-1)*N_Beta+B_PT,4)
      Tau_DU=UV_TABLE((A_PT-1)*N_Beta+(B_PT+1),4)
      Tau_UD=UV_TABLE(A_PT*N_Beta+B_PT,4)
      Tau_UU=UV_TABLE(A_PT*N_Beta+(B_PT+1),4)

c     Find the position within the interval normalized to [0,1]

      D_Alpha=(Alpha-Alpha_G(A_PT))/(Alpha_G(A_PT+1)-Alpha_G(A_PT))
      D_Beta=(Beta-Beta_G(B_PT))/(Beta_G(B_PT+1)-Beta_G(B_PT))

c     Do a 2 dimensional interpolation of the values of G0 and Tau

      G0=(1.0-D_Alpha)*(1.0-D_Beta)*G0_DD +
     &  (    D_Alpha)*(1.0-D_Beta)*G0_UD +
     &  (1.0-D_Alpha)*(    D_Beta)*G0_DU +
     &  (    D_Alpha)*(    D_Beta)*G0_UU
      Tau=(1.0-D_Alpha)*(1.0-D_Beta)*Tau_DD +
     &  (    D_Alpha)*(1.0-D_Beta)*Tau_UD +
     &  (1.0-D_Alpha)*(    D_Beta)*Tau_DU +
     &  (    D_Alpha)*(    D_Beta)*Tau_UU

      RETURN

      END





c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c SUBROUTINE UVGRID_ESTI(G0,Tau,Alpha,Beta)
c
c     Calculates a first estimatio of Alpha and Beta for G0 and Tau.
c     This is simply done by a pre-calculated look-up table.
c
c     Returns (Alpha,Beta)=(-100,-100) if (G0,Tau) is not in the range
c     of the look-up table.
c
      SUBROUTINE UVGRID_ESTI(G0,Tau,Alpha,Beta)

      IMPLICIT NONE

c     Common-Block 'UV_ESTIMAT' : Estimation of the values of
c                                 (Alpha,Beta) from (G0,Tau)
c       UV_TAU: Gridpoints of UV_A_EST/UV_B_EST in direction of Tau
c       UV_G0: Gridpoints of UV_A_EST/UV_B_EST in direction of G0
c       UV_A_EST: Estimated values of Alpha
c       UV_B_EST: Estimated values of Beta
      REAL UV_Tau(60),UV_G0(60),UV_A_EST(60,60),UV_B_EST(60,60)
      COMMON/UV_ESTIMAT/ UV_Tau,UV_G0,UV_A_EST,UV_B_EST

c     Parameters   
c       G0: FUV Flux in ISRF (Input)
c       Tau: FUV Attenuation in A_V (Input)
c       Alpha: Alpha (Output)
c       Beta: Beta (Output)
      REAL G0,Tau,Alpha,Beta

c     Variables
c       D_X,D_Y: Position in the intervals of the tables
c                UV_A_EST/UV_B_EST
c       X_Pt,Y_Pt: Index in the intervals of the tables
c                  UV_A_EST/UV_B_EST
      REAL D_X,D_Y
      INTEGER X_PT,Y_PT,I

cccccccccccccccccccccc       
c     Code


c     Find the index within the table for both directions G0 and Tau

      X_PT=-100
      Y_PT=-100
      DO 100 I=1,59
         IF(G0.GE.UV_G0(I).AND.G0.LE.UV_G0(I+1)) X_PT=I
 100  CONTINUE
      DO 200 I=1,59
        IF(Tau.GE.UV_Tau(I).AND.Tau.LE.UV_Tau(I+1)) Y_PT=I
 200  CONTINUE

c     Outside the tabulated range

      IF((X_PT.LT.-10).OR.(Y_PT.LT.-10)) THEN 
        Alpha=-100.0
        Beta=-100.0
	RETURN
      END IF

c     Find the position within the interval normalized to [0,1]

      D_X=(G0-UV_G0(X_PT))/(UV_G0(X_PT+1)-UV_G0(X_PT))
      D_Y=(Tau-UV_Tau(Y_PT))/(UV_Tau(Y_PT+1)-UV_Tau(Y_PT))

c     Do a 2 dimensional interpolation of the values of Alpha and Beta

      Alpha=(1.0-D_X)*(1.0-D_Y)*UV_A_EST(X_PT,Y_PT)+
     &  (    D_X)*(1.0-D_Y)*UV_A_EST(X_PT+1,Y_PT)+
     &  (1.0-D_X)*(    D_Y)*UV_A_EST(X_PT,Y_PT+1)+
     &  (    D_X)*(    D_Y)*UV_A_EST(X_PT+1,Y_PT+1)

      Beta= (1.0-D_X)*(1.0-D_Y)*UV_B_EST(X_PT,Y_PT)+
     &  (    D_X)*(1.0-D_Y)*UV_B_EST(X_PT+1,Y_PT)+
     &  (1.0-D_X)*(    D_Y)*UV_B_EST(X_PT,Y_PT+1)+
     &  (    D_X)*(    D_Y)*UV_B_EST(X_PT+1,Y_PT+1)

      RETURN      
      
      END





c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c REAL FUNCTION XRate(FLX,NH,TX)
c
c     Calculates the ionization rate by X-rays following the integral in
c     Appendix A of the grid-paper.
c
c
      REAL FUNCTION XRate(FLX,NH,TX)

      IMPLICIT NONE
      
c     Constants      
c       Boltz: Boltzmann constant erg K^-1
c       Boltzev: Boltzmann constant in units of eV
c       Xe_Min: Minimum photon energy (lower bound of the integral)
c       Xe_Max: Maximum photon energy (upper bound of the integral)
c       I_Steps: Number of energies for the integration
      
      REAL Boltz,BoltzEv,Xe_Min,Xe_Max
      INTEGER I_Steps
      PARAMETER(Boltz=1.38054E-16)
      PARAMETER(BoltzEv=8.61664041E-05)
      PARAMETER(Xe_Min=1.0E3)
      PARAMETER(Xe_Max=1.0E5)
      PARAMETER(I_Steps=2000)
     
c     Functions
      REAL SIGMA_PH,SIGMA_CO
      EXTERNAL SIGMA_PH,SIGMA_CO
     
c     Parameters
c       FLX: X-ray flux erg s^-1 cm^-2
c       NH: Attenuating column density cm^-2
c       TX: Plasma temmperature K
      REAL FLX,NH,TX 

c     Variables
c       Fako: Normalization for a thermal X-ray emission (See paper)
c       F_FX: Flux at a given photon energy
c       F_E: Ionization rate at a given energy
c       Rate: Summed up ionization rate
c       W_Mean: Mean energy per ion-pair
c       Energy,Energy2: Photon energy for the summation
c       S_PH: Photo ionization cross-section
c       S_CO: Compoton ionization cross-section
c       Lx_Min: Log10 of the minimum photon energy Xe_Min
c       Lx_Max: Log10 of the maximum photon energy Xe_Max
c       Lx_Step: Step width in the log space
      REAL Fako,F_FX,F_E,Rate,W_Mean
      REAL Energy,Energy2,S_PH,S_CO
      REAL Lx_Min,Lx_Max,Lx_Step
      INTEGER I
              
cccccccccccccccccccccc       
c     Code

      Lx_Min=ALOG10(Xe_Min)
      Lx_Max=ALOG10(Xe_Max)
      Lx_Step=(Lx_Max-Lx_Min)/I_Steps

c     Calculate the normalization of the thermal spectrum

      Fako=FLX/(Boltz*TX)
     1  /(EXP(-Xe_Min/(BoltzEv*TX))-EXP(-Xe_Max/(BoltzEv*TX)))

      Rate=0.0
      
c     Replace the integral Eq (5) in the paper by a summation

      DO 100 I=0,I_Steps-1
      
        Energy=10.0**(Lx_Min+I*Lx_Step)
        Energy2=10.0**(Lx_Min+(I+1)*Lx_Step)
        S_PH=SIGMA_PH(Energy)
        S_CO=SIGMA_CO(Energy)
	
c     Mean photon energy per ion-pair (Appendix A)

	IF(Energy.LT.1.e3) THEN
	  W_Mean=23.65-(ALOG10(Energy)-2.0)*2.7
	ELSE
	  W_Mean=20.95	   
	END IF
	
c     Normalized X-ray spectra at the photon energy Energy

        F_FX=Fako*EXP(-Energy/(BoltzEv*TX))*EXP(-S_PH*NH)*1.602E-12
	
c     Ionization rate at the photon energy

        F_E=F_FX*(S_PH+S_CO)*1.0/(2.0*W_Mean*0.5)*6.2415E11

        Rate=Rate+(Energy2-Energy)*F_E
	
100   CONTINUE

      XRate=Rate
      
      RETURN

      END





c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c REAL FUNCTION SIGMA_PH(EN)
c
c     Interpolation of the photo ionization from table A.1 (paper)
c
c
      REAL FUNCTION SIGMA_PH(EN)
      IMPLICIT NONE
      
c     Parameters
c       En: Photon energy in eV
      REAL En

c     Variables
c       E_Min: Energy at the lower bound of the interval that contains
c              the energy
c       E_Max: Energy at the upper bound of the interval that contains
c              the energy
c       Sig_Min: Photoionization cross-section at E_min cm^2
c       Sig_Max: Photoionization cross-section at E_min cm^2
c       Alpha: Position in the interval E_min to E_max for the
c              interpolation
      REAL E_Min,E_Max,Sig_Min,Sig_Max,Alpha
      REAL Sig_Tab(32,2)
      INTEGER I      
            
c     Data:
c       Sig_Tab: Table containing photo ionization cross-sections
c                (A.1 in the paper). The first column gives the photon
c                energy in eV, the second column the cross-section in
c                cm^2.
      DATA (Sig_Tab(1,I),I=1,2)/ 100.0     ,6.020e-20 /
      DATA (Sig_Tab(2,I),I=1,2)/ 170.0     ,1.330e-20 /
      DATA (Sig_Tab(3,I),I=1,2)/ 291.0     ,2.710e-21 /
      DATA (Sig_Tab(4,I),I=1,2)/ 291.1     ,3.030e-21 /
      DATA (Sig_Tab(5,I),I=1,2)/ 404.7     ,1.150e-21 /
      DATA (Sig_Tab(6,I),I=1,2)/ 404.8     ,1.220e-21 /
      DATA (Sig_Tab(7,I),I=1,2)/ 537.9     ,5.290e-22 /
      DATA (Sig_Tab(8,I),I=1,2)/ 538.0     ,7.970e-22 /
      DATA (Sig_Tab(9,I),I=1,2)/ 723.9     ,3.590e-22 /
      DATA (Sig_Tab(10,I),I=1,2)/ 724.0    ,3.990e-22 /
      DATA (Sig_Tab(11,I),I=1,2)/ 856.9    ,2.540e-22 /
      DATA (Sig_Tab(12,I),I=1,2)/ 857.0    ,2.590e-22 /
      DATA (Sig_Tab(13,I),I=1,2)/ 870.0    ,2.480e-22 /
      DATA (Sig_Tab(14,I),I=1,2)/ 870.1    ,2.890e-22 /
      DATA (Sig_Tab(15,I),I=1,2)/ 1311.0   ,9.750e-23 /
      DATA (Sig_Tab(16,I),I=1,2)/ 1311.1   ,1.060e-22 /
      DATA (Sig_Tab(17,I),I=1,2)/ 1846.0   ,4.130e-23 /
      DATA (Sig_Tab(18,I),I=1,2)/ 1846.1   ,4.610e-23 /
      DATA (Sig_Tab(19,I),I=1,2)/ 2477.0   ,2.030e-23 /
      DATA (Sig_Tab(20,I),I=1,2)/ 2477.1   ,2.230e-23 /
      DATA (Sig_Tab(21,I),I=1,2)/ 3203.0   ,1.090e-23 /
      DATA (Sig_Tab(22,I),I=1,2)/ 3203.1   ,1.110e-23 /
      DATA (Sig_Tab(23,I),I=1,2)/ 4043.0   ,5.760e-24 /
      DATA (Sig_Tab(24,I),I=1,2)/ 4043.1   ,5.890e-24 /
      DATA (Sig_Tab(25,I),I=1,2)/ 5995.9   ,1.890e-24 /
      DATA (Sig_Tab(26,I),I=1,2)/ 5996.0   ,1.910e-24 /
      DATA (Sig_Tab(27,I),I=1,2)/ 7123.9   ,1.150e-24 /
      DATA (Sig_Tab(28,I),I=1,2)/ 7124.0   ,2.240e-24 /
      DATA (Sig_Tab(29,I),I=1,2)/ 8347.9   ,1.450e-24 /
      DATA (Sig_Tab(30,I),I=1,2)/ 8348.0   ,1.500e-24 /
      DATA (Sig_Tab(31,I),I=1,2)/ 28900.0  ,4.240e-26 /
      DATA (Sig_Tab(32,I),I=1,2)/ 100000.0 ,1.040e-27 /

cccccccccccccccccccccc       
c     Code


c     Find the index of the interval in the table for the given
c     photon energy.

      IF(En.LT.Sig_Tab(1,1)) THEN
        SIGMA_PH=Sig_Tab(1,2)
        RETURN
      END IF
      IF(En.GE.Sig_Tab(32,1)) THEN
        SIGMA_PH=Sig_Tab(32,2)
        RETURN
      END IF

      DO 100 I=1,31
        IF((Sig_Tab(I,1).LE.En).AND.(En.LT.Sig_Tab(I+1,1))) THEN
          E_Min=Sig_Tab(I,1)
          E_Max=Sig_Tab(I+1,1)
          Sig_Min=Sig_Tab(I,2)
          Sig_Max=Sig_Tab(I+1,2)
        END IF
100   CONTINUE

c     Do an interpolation in the log-log space

      Alpha=(ALOG10(En)-ALOG10(E_Min))/(ALOG10(E_Max)-ALOG10(E_Min))
      SIGMA_PH=10.0**(ALOG10(Sig_Min)*(1.0-Alpha)
     &  +ALOG10(Sig_Max)*Alpha)

      RETURN
      
      END





c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c REAL FUNCTION SIGMA_PH(EN)
c
c     Approximation for the compton-cross section
c
c
      REAL FUNCTION SIGMA_CO(EN)

      IMPLICIT NONE
      
c     Parameters
c       En: Photon energy in eV
      REAL EN

c     Variables
c       Enlog: Log of the energy
      REAL Enlog

cccccccccccccccccccccc       
c     Code

      Enlog=ALOG10(En)
      
      IF(Enlog.LT.3.0) Enlog=3.0
      IF(Enlog.GT.5.0) Enlog=5.0

c     Use equation A.3 from the grid-paper ot calculate the compton
c     cross-section.

      IF(Enlog.LT.4.0) THEN
        SIGMA_CO=(2.869674e-23)-(2.6364914e-23)*Enlog
     &       +(7.931175e-24)*Enlog**2-(7.74014e-25)*Enlog**3
      ELSE
        SIGMA_CO=(-2.374894e-24)+(1.423853e-24)*Enlog
     &       -(1.70095e-25)*Enlog**2
      END IF

      RETURN      
      
      END





c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c REAL FUNCTION XRATE_INT(FLX,NH,TX)
c
c     Interpolation of the X-ray ionization rate from a table of pre-
c     calculated ionization rates. The X-ray flux FLX, the attenuating
c     column density NH and the plasma temperature TX enter this
c     function.
c
c     The integration used for the x-ray ionization rate takes a
c     significant amount of calculation time. instead of a full
c     calculation, we interpolate on a grid of precalculated ionization
c     rates and using the fact that the ionization rate is proportional
c     to the flux. we thus need only a 2 dimensional interpolation.
c     The routine is proofed for tx > 1e6 to be accurate within 5 % and
c     much better for tx > 1e7.
c
c
      REAL FUNCTION XRATE_INT(FLX,NH,TX)

      IMPLICIT NONE

c     Common-Block 'X_IONDAT_100eV' or 'X_IONDAT_1keV':
c       X-ray ionization rate at an X-ray flux of 1 erg s^-1 cm^-2
c       for different values of the attenuating column density and
c       the plasma temperature.
c       Select one of the two data-sets 'X_IONDAT_100eV'/'X_IONDAT_1keV'
c       to chose the lower bound of the integral for the ionizaton rate
c       X_SEC: Table of ionization rates s^-1
c       X_NH: Griding along the axis for the attenuating
c             column density cm^-2
c       X_TX: Griding along the axis for the plasma temperature K
      REAL X_SEC(40,40), X_NH(40), X_TX(40)
c      COMMON/X_IONDAT_100eV/ X_SEC, X_NH, X_TX
      COMMON/X_IONDAT_1keV/ X_SEC, X_NH, X_TX 

c     Parameters
c       FLX: X-ray flux erg s^-1 cm^-2
c       NH: Attenuating column density cm^-2
c       TX: Plasma temmperature K
      REAL FLX,NH,TX

c     Variables
c       D_X,D_Y: Position in the interval between two table-points,
c                normalized to 0 - 1.
c       Xion_Log: Interpolated fractional abundance
c       X_Pt, Y_Pt: Index in the table X_SEC
      INTEGER X_Pt,Y_Pt
      REAL D_X,D_Y,Xion_Log
      INTEGER I
            
cccccccccccccccccccccc                  
c     Code


cdebug: only if you want to see warnings...
c      IF(TX.LT.1.0E6) THEN
c         WRITE(*,*) 'GRID-WARNING: INTERPOLATION OF XRATE'
c         WRITE(*,*) 'GRID-WARNING: ONLY ACCURATE FOR TX > 1E6!'
c         WRITE(*,*) 'GRID-WARNING: USE XRATE INSTEAD OF XRATE_INT!'
c      END IF

c     Find the index of the interval containing the plasma temperature
c     and the attenuating column density

      X_Pt=-100
      Y_Pt=-100
      DO 100 I=1,39
        IF(TX.GE.X_TX(I).AND.TX.LE.X_TX(I+1)) X_Pt=I
 100  CONTINUE
      DO 200 I=1,39
        IF(NH.GE.X_NH(I).AND.TX.LE.X_NH(I+1)) Y_Pt=I
 200  CONTINUE

      IF(X_Pt.EQ.-100.OR.Y_Pt.EQ.-100) THEN 
        WRITE(*,'(A)') 'OUTSIDE XION-GRID'
        STOP
      END IF

c     Do an interpolation in the log-log space

      D_X=(ALOG10(TX)-ALOG10(X_TX(X_Pt)))/
     & (ALOG10(X_TX(X_Pt+1))-ALOG10(X_TX(X_Pt)))
      D_Y=(ALOG10(NH)-ALOG10(X_NH(Y_Pt)))/
     & (ALOG10(X_NH(Y_Pt+1))-ALOG10(X_NH(Y_Pt)))

      Xion_Log=(1.0-D_X)*(1.0-D_Y)*ALOG10(X_SEC(X_Pt,Y_Pt))+
     &  (    D_X)*(1.0-D_Y)*ALOG10(X_SEC(X_Pt+1,Y_Pt))+
     &  (1.0-D_X)*(    D_Y)*ALOG10(X_SEC(X_Pt,Y_Pt+1))+
     &  (    D_X)*(    D_Y)*ALOG10(X_SEC(X_Pt+1,Y_Pt+1))

c     The rate is given for a flux of 1 erg s^-1 cm^-2 in the table,
c     thus it has to be scaled it with FLX. The look-up table stores
c     the log10 of the ionization rates.

      XRATE_INT=FLX*10.0**Xion_Log

      RETURN
      
      END





c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
c
c Block Data: The remainder of this code consists of data for the
c             UV-Coordinate system and the X-ray ionization rate.
c
c             I have chosen this rather nasty BLOCK DATA statements
c             instead of external files in order to avoid confusions
c             with different versions of the grid.





c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c BLOCK DATA UVGRIDDATA
c
c     coordinate system for uv-grid
c
c
      BLOCK DATA UVGRIDDATA

      IMPLICIT NONE

c     Common-Block 'UV_DATA' : Coordinate system for the coordinate
c                              transformation (G0,Tau) --> (Alpha,Beta)
c       N_Alpha: Number of points in direction of Alpha
c       N_Beta: Nubmer of points in direction of Beta
c       UV_TABLE: Table containing with lines containing
c                 (log(Alpha),log(Beta),log(g0),log(Tau))
c       Alpha_G: Points of the UV_TABLE in log(Alpha)-direction
c       Beta_G: Points of the UV_TABLE in log(Beta)-direction
      INTEGER N_ALPHA,N_BETA
      PARAMETER(N_ALPHA=27)
      PARAMETER(N_BETA=28)
      REAL UV_TABLE(N_ALPHA*N_BETA,4),ALPHA_G(N_ALPHA),BETA_G(N_BETA)
      COMMON/UV_DATA/ UV_TABLE,ALPHA_G,BETA_G
      
      INTEGER I

      DATA (ALPHA_G(I),I=1,N_ALPHA) /-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,
     &    0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,4.8,5.05,5.2,
     &    5.32,5.45,5.7,6.0,6.5,7.0,7.5,8.0/
      DATA (BETA_G(I),I=1,N_BETA) /-1.125,-1.0,-0.875,-0.75,-0.625,
     &   -0.5,-0.375,-0.25,-0.125,0.0,0.125,0.25,0.3,
     &    0.375,0.5,0.6,0.625,0.75,0.8,0.875,1.0,1.125,
     &    1.25,1.375,1.5,1.625,1.75,1.875/
      DATA(UV_TABLE(1,I),I=1,4)/-3.0000,-1.1250,8.8139,-1.6531/
      DATA(UV_TABLE(2,I),I=1,4)/-3.0000,-1.0000,8.8437,-1.5805/
      DATA(UV_TABLE(3,I),I=1,4)/-3.0000,-0.8750,8.8735,-1.5078/
      DATA(UV_TABLE(4,I),I=1,4)/-3.0000,-0.7500,8.9034,-1.4352/
      DATA(UV_TABLE(5,I),I=1,4)/-3.0000,-0.6250,8.9332,-1.3625/
      DATA(UV_TABLE(6,I),I=1,4)/-3.0000,-0.5000,8.9630,-1.2899/
      DATA(UV_TABLE(7,I),I=1,4)/-3.0000,-0.3750,8.9929,-1.2172/
      DATA(UV_TABLE(8,I),I=1,4)/-3.0000,-0.2500,9.0227,-1.1446/
      DATA(UV_TABLE(9,I),I=1,4)/-3.0000,-0.1250,9.0562,-1.0722/
      DATA(UV_TABLE(10,I),I=1,4)/-3.0000,0.0000,9.0963,-1.0003/
      DATA(UV_TABLE(11,I),I=1,4)/-3.0000,0.1250,9.1401,-0.9288/
      DATA(UV_TABLE(12,I),I=1,4)/-3.0000,0.2500,9.1837,-0.8573/
      DATA(UV_TABLE(13,I),I=1,4)/-3.0000,0.3000,9.2011,-0.8287/
      DATA(UV_TABLE(14,I),I=1,4)/-3.0000,0.3750,9.2272,-0.7858/
      DATA(UV_TABLE(15,I),I=1,4)/-3.0000,0.5000,9.2708,-0.7142/
      DATA(UV_TABLE(16,I),I=1,4)/-3.0000,0.6000,9.3056,-0.6570/
      DATA(UV_TABLE(17,I),I=1,4)/-3.0000,0.6250,9.3143,-0.6427/
      DATA(UV_TABLE(18,I),I=1,4)/-3.0000,0.7500,9.3579,-0.5712/
      DATA(UV_TABLE(19,I),I=1,4)/-3.0000,0.8000,9.3942,-0.5456/
      DATA(UV_TABLE(20,I),I=1,4)/-3.0000,0.8750,9.4486,-0.5071/
      DATA(UV_TABLE(21,I),I=1,4)/-3.0000,1.0000,9.5392,-0.4430/
      DATA(UV_TABLE(22,I),I=1,4)/-3.0000,1.1250,9.6299,-0.3789/
      DATA(UV_TABLE(23,I),I=1,4)/-3.0000,1.2500,9.7420,-0.3206/
      DATA(UV_TABLE(24,I),I=1,4)/-3.0000,1.3750,9.8541,-0.2622/
      DATA(UV_TABLE(25,I),I=1,4)/-3.0000,1.5000,9.9662,-0.2038/
      DATA(UV_TABLE(26,I),I=1,4)/-3.0000,1.6250,10.0851,-0.1481/
      DATA(UV_TABLE(27,I),I=1,4)/-3.0000,1.7500,10.2338,-0.1048/
      DATA(UV_TABLE(28,I),I=1,4)/-3.0000,1.8750,10.3826,-0.0614/
      DATA(UV_TABLE(29,I),I=1,4)/-2.5000,-1.1250,8.3007,-1.6539/
      DATA(UV_TABLE(30,I),I=1,4)/-2.5000,-1.0000,8.3305,-1.5812/
      DATA(UV_TABLE(31,I),I=1,4)/-2.5000,-0.8750,8.3604,-1.5086/
      DATA(UV_TABLE(32,I),I=1,4)/-2.5000,-0.7500,8.3902,-1.4359/
      DATA(UV_TABLE(33,I),I=1,4)/-2.5000,-0.6250,8.4200,-1.3633/
      DATA(UV_TABLE(34,I),I=1,4)/-2.5000,-0.5000,8.4499,-1.2906/
      DATA(UV_TABLE(35,I),I=1,4)/-2.5000,-0.3750,8.4797,-1.2180/
      DATA(UV_TABLE(36,I),I=1,4)/-2.5000,-0.2500,8.5095,-1.1453/
      DATA(UV_TABLE(37,I),I=1,4)/-2.5000,-0.1250,8.5395,-1.0727/
      DATA(UV_TABLE(38,I),I=1,4)/-2.5000,0.0000,8.5777,-1.0006/
      DATA(UV_TABLE(39,I),I=1,4)/-2.5000,0.1250,8.6191,-0.9289/
      DATA(UV_TABLE(40,I),I=1,4)/-2.5000,0.2500,8.6629,-0.8574/
      DATA(UV_TABLE(41,I),I=1,4)/-2.5000,0.3000,8.6803,-0.8288/
      DATA(UV_TABLE(42,I),I=1,4)/-2.5000,0.3750,8.7065,-0.7859/
      DATA(UV_TABLE(43,I),I=1,4)/-2.5000,0.5000,8.7500,-0.7143/
      DATA(UV_TABLE(44,I),I=1,4)/-2.5000,0.6000,8.7849,-0.6571/
      DATA(UV_TABLE(45,I),I=1,4)/-2.5000,0.6250,8.7936,-0.6428/
      DATA(UV_TABLE(46,I),I=1,4)/-2.5000,0.7500,8.8395,-0.5716/
      DATA(UV_TABLE(47,I),I=1,4)/-2.5000,0.8000,8.8757,-0.5460/
      DATA(UV_TABLE(48,I),I=1,4)/-2.5000,0.8750,8.9301,-0.5076/
      DATA(UV_TABLE(49,I),I=1,4)/-2.5000,1.0000,9.0208,-0.4435/
      DATA(UV_TABLE(50,I),I=1,4)/-2.5000,1.1250,9.1115,-0.3794/
      DATA(UV_TABLE(51,I),I=1,4)/-2.5000,1.2500,9.2207,-0.3203/
      DATA(UV_TABLE(52,I),I=1,4)/-2.5000,1.3750,9.3328,-0.2619/
      DATA(UV_TABLE(53,I),I=1,4)/-2.5000,1.5000,9.4449,-0.2035/
      DATA(UV_TABLE(54,I),I=1,4)/-2.5000,1.6250,9.5637,-0.1478/
      DATA(UV_TABLE(55,I),I=1,4)/-2.5000,1.7500,9.7125,-0.1045/
      DATA(UV_TABLE(56,I),I=1,4)/-2.5000,1.8750,9.8612,-0.0611/
      DATA(UV_TABLE(57,I),I=1,4)/-2.0000,-1.1250,7.7888,-1.6537/
      DATA(UV_TABLE(58,I),I=1,4)/-2.0000,-1.0000,7.8186,-1.5810/
      DATA(UV_TABLE(59,I),I=1,4)/-2.0000,-0.8750,7.8484,-1.5084/
      DATA(UV_TABLE(60,I),I=1,4)/-2.0000,-0.7500,7.8783,-1.4357/
      DATA(UV_TABLE(61,I),I=1,4)/-2.0000,-0.6250,7.9081,-1.3631/
      DATA(UV_TABLE(62,I),I=1,4)/-2.0000,-0.5000,7.9379,-1.2904/
      DATA(UV_TABLE(63,I),I=1,4)/-2.0000,-0.3750,7.9678,-1.2178/
      DATA(UV_TABLE(64,I),I=1,4)/-2.0000,-0.2500,7.9976,-1.1451/
      DATA(UV_TABLE(65,I),I=1,4)/-2.0000,-0.1250,8.0274,-1.0725/
      DATA(UV_TABLE(66,I),I=1,4)/-2.0000,0.0000,8.0610,-1.0001/
      DATA(UV_TABLE(67,I),I=1,4)/-2.0000,0.1250,8.1001,-0.9281/
      DATA(UV_TABLE(68,I),I=1,4)/-2.0000,0.2500,8.1431,-0.8566/
      DATA(UV_TABLE(69,I),I=1,4)/-2.0000,0.3000,8.1605,-0.8280/
      DATA(UV_TABLE(70,I),I=1,4)/-2.0000,0.3750,8.1867,-0.7850/
      DATA(UV_TABLE(71,I),I=1,4)/-2.0000,0.5000,8.2303,-0.7135/
      DATA(UV_TABLE(72,I),I=1,4)/-2.0000,0.6000,8.2651,-0.6563/
      DATA(UV_TABLE(73,I),I=1,4)/-2.0000,0.6250,8.2738,-0.6420/
      DATA(UV_TABLE(74,I),I=1,4)/-2.0000,0.7500,8.3197,-0.5708/
      DATA(UV_TABLE(75,I),I=1,4)/-2.0000,0.8000,8.3560,-0.5452/
      DATA(UV_TABLE(76,I),I=1,4)/-2.0000,0.8750,8.4104,-0.5067/
      DATA(UV_TABLE(77,I),I=1,4)/-2.0000,1.0000,8.5011,-0.4427/
      DATA(UV_TABLE(78,I),I=1,4)/-2.0000,1.1250,8.5917,-0.3786/
      DATA(UV_TABLE(79,I),I=1,4)/-2.0000,1.2500,8.7009,-0.3194/
      DATA(UV_TABLE(80,I),I=1,4)/-2.0000,1.3750,8.8130,-0.2610/
      DATA(UV_TABLE(81,I),I=1,4)/-2.0000,1.5000,8.9251,-0.2026/
      DATA(UV_TABLE(82,I),I=1,4)/-2.0000,1.6250,9.0440,-0.1470/
      DATA(UV_TABLE(83,I),I=1,4)/-2.0000,1.7500,9.1927,-0.1036/
      DATA(UV_TABLE(84,I),I=1,4)/-2.0000,1.8750,9.3415,-0.0603/
      DATA(UV_TABLE(85,I),I=1,4)/-1.5000,-1.1250,7.2767,-1.6545/
      DATA(UV_TABLE(86,I),I=1,4)/-1.5000,-1.0000,7.3065,-1.5818/
      DATA(UV_TABLE(87,I),I=1,4)/-1.5000,-0.8750,7.3363,-1.5092/
      DATA(UV_TABLE(88,I),I=1,4)/-1.5000,-0.7500,7.3662,-1.4365/
      DATA(UV_TABLE(89,I),I=1,4)/-1.5000,-0.6250,7.3960,-1.3639/
      DATA(UV_TABLE(90,I),I=1,4)/-1.5000,-0.5000,7.4258,-1.2912/
      DATA(UV_TABLE(91,I),I=1,4)/-1.5000,-0.3750,7.4557,-1.2186/
      DATA(UV_TABLE(92,I),I=1,4)/-1.5000,-0.2500,7.4855,-1.1459/
      DATA(UV_TABLE(93,I),I=1,4)/-1.5000,-0.1250,7.5153,-1.0733/
      DATA(UV_TABLE(94,I),I=1,4)/-1.5000,0.0000,7.5453,-1.0006/
      DATA(UV_TABLE(95,I),I=1,4)/-1.5000,0.1250,7.5818,-0.9285/
      DATA(UV_TABLE(96,I),I=1,4)/-1.5000,0.2500,7.6221,-0.8567/
      DATA(UV_TABLE(97,I),I=1,4)/-1.5000,0.3000,7.6396,-0.8281/
      DATA(UV_TABLE(98,I),I=1,4)/-1.5000,0.3750,7.6659,-0.7851/
      DATA(UV_TABLE(99,I),I=1,4)/-1.5000,0.5000,7.7095,-0.7136/
      DATA(UV_TABLE(100,I),I=1,4)/-1.5000,0.6000,7.7443,-0.6564/
      DATA(UV_TABLE(101,I),I=1,4)/-1.5000,0.6250,7.7530,-0.6421/
      DATA(UV_TABLE(102,I),I=1,4)/-1.5000,0.7500,7.7989,-0.5709/
      DATA(UV_TABLE(103,I),I=1,4)/-1.5000,0.8000,7.8352,-0.5453/
      DATA(UV_TABLE(104,I),I=1,4)/-1.5000,0.8750,7.8896,-0.5068/
      DATA(UV_TABLE(105,I),I=1,4)/-1.5000,1.0000,7.9803,-0.4428/
      DATA(UV_TABLE(106,I),I=1,4)/-1.5000,1.1250,8.0709,-0.3787/
      DATA(UV_TABLE(107,I),I=1,4)/-1.5000,1.2500,8.1801,-0.3195/
      DATA(UV_TABLE(108,I),I=1,4)/-1.5000,1.3750,8.2922,-0.2611/
      DATA(UV_TABLE(109,I),I=1,4)/-1.5000,1.5000,8.4043,-0.2027/
      DATA(UV_TABLE(110,I),I=1,4)/-1.5000,1.6250,8.5263,-0.1482/
      DATA(UV_TABLE(111,I),I=1,4)/-1.5000,1.7500,8.6750,-0.1048/
      DATA(UV_TABLE(112,I),I=1,4)/-1.5000,1.8750,8.8238,-0.0614/
      DATA(UV_TABLE(113,I),I=1,4)/-1.0000,-1.1250,6.7644,-1.6543/
      DATA(UV_TABLE(114,I),I=1,4)/-1.0000,-1.0000,6.7942,-1.5816/
      DATA(UV_TABLE(115,I),I=1,4)/-1.0000,-0.8750,6.8241,-1.5090/
      DATA(UV_TABLE(116,I),I=1,4)/-1.0000,-0.7500,6.8539,-1.4363/
      DATA(UV_TABLE(117,I),I=1,4)/-1.0000,-0.6250,6.8837,-1.3637/
      DATA(UV_TABLE(118,I),I=1,4)/-1.0000,-0.5000,6.9136,-1.2910/
      DATA(UV_TABLE(119,I),I=1,4)/-1.0000,-0.3750,6.9434,-1.2184/
      DATA(UV_TABLE(120,I),I=1,4)/-1.0000,-0.2500,6.9732,-1.1457/
      DATA(UV_TABLE(121,I),I=1,4)/-1.0000,-0.1250,7.0031,-1.0731/
      DATA(UV_TABLE(122,I),I=1,4)/-1.0000,0.0000,7.0329,-1.0004/
      DATA(UV_TABLE(123,I),I=1,4)/-1.0000,0.1250,7.0659,-0.9280/
      DATA(UV_TABLE(124,I),I=1,4)/-1.0000,0.2500,7.1046,-0.8560/
      DATA(UV_TABLE(125,I),I=1,4)/-1.0000,0.3000,7.1214,-0.8273/
      DATA(UV_TABLE(126,I),I=1,4)/-1.0000,0.3750,7.1466,-0.7843/
      DATA(UV_TABLE(127,I),I=1,4)/-1.0000,0.5000,7.1904,-0.7128/
      DATA(UV_TABLE(128,I),I=1,4)/-1.0000,0.6000,7.2253,-0.6556/
      DATA(UV_TABLE(129,I),I=1,4)/-1.0000,0.6250,7.2340,-0.6413/
      DATA(UV_TABLE(130,I),I=1,4)/-1.0000,0.7500,7.2952,-0.5724/
      DATA(UV_TABLE(131,I),I=1,4)/-1.0000,0.8000,7.3294,-0.5463/
      DATA(UV_TABLE(132,I),I=1,4)/-1.0000,0.8750,7.3807,-0.5072/
      DATA(UV_TABLE(133,I),I=1,4)/-1.0000,1.0000,7.4714,-0.4431/
      DATA(UV_TABLE(134,I),I=1,4)/-1.0000,1.1250,7.5621,-0.3790/
      DATA(UV_TABLE(135,I),I=1,4)/-1.0000,1.2500,7.6684,-0.3191/
      DATA(UV_TABLE(136,I),I=1,4)/-1.0000,1.3750,7.7804,-0.2607/
      DATA(UV_TABLE(137,I),I=1,4)/-1.0000,1.5000,7.8925,-0.2023/
      DATA(UV_TABLE(138,I),I=1,4)/-1.0000,1.6250,8.0152,-0.1475/
      DATA(UV_TABLE(139,I),I=1,4)/-1.0000,1.7500,8.1607,-0.1025/
      DATA(UV_TABLE(140,I),I=1,4)/-1.0000,1.8750,8.3094,-0.0591/
      DATA(UV_TABLE(141,I),I=1,4)/-0.5000,-1.1250,6.1478,-1.6453/
      DATA(UV_TABLE(142,I),I=1,4)/-0.5000,-1.0000,6.1978,-1.5744/
      DATA(UV_TABLE(143,I),I=1,4)/-0.5000,-0.8750,6.2478,-1.5036/
      DATA(UV_TABLE(144,I),I=1,4)/-0.5000,-0.7500,6.2978,-1.4327/
      DATA(UV_TABLE(145,I),I=1,4)/-0.5000,-0.6250,6.3362,-1.3607/
      DATA(UV_TABLE(146,I),I=1,4)/-0.5000,-0.5000,6.3740,-1.2886/
      DATA(UV_TABLE(147,I),I=1,4)/-0.5000,-0.3750,6.4118,-1.2166/
      DATA(UV_TABLE(148,I),I=1,4)/-0.5000,-0.2500,6.4496,-1.1445/
      DATA(UV_TABLE(149,I),I=1,4)/-0.5000,-0.1250,6.4874,-1.0725/
      DATA(UV_TABLE(150,I),I=1,4)/-0.5000,0.0000,6.5252,-1.0004/
      DATA(UV_TABLE(151,I),I=1,4)/-0.5000,0.1250,6.5732,-0.9295/
      DATA(UV_TABLE(152,I),I=1,4)/-0.5000,0.2500,6.6328,-0.8598/
      DATA(UV_TABLE(153,I),I=1,4)/-0.5000,0.3000,6.6578,-0.8321/
      DATA(UV_TABLE(154,I),I=1,4)/-0.5000,0.3750,6.6953,-0.7906/
      DATA(UV_TABLE(155,I),I=1,4)/-0.5000,0.5000,6.7595,-0.7216/
      DATA(UV_TABLE(156,I),I=1,4)/-0.5000,0.6000,6.8113,-0.6665/
      DATA(UV_TABLE(157,I),I=1,4)/-0.5000,0.6250,6.8243,-0.6527/
      DATA(UV_TABLE(158,I),I=1,4)/-0.5000,0.7500,6.8982,-0.5853/
      DATA(UV_TABLE(159,I),I=1,4)/-0.5000,0.8000,6.9298,-0.5587/
      DATA(UV_TABLE(160,I),I=1,4)/-0.5000,0.8750,6.9771,-0.5188/
      DATA(UV_TABLE(161,I),I=1,4)/-0.5000,1.0000,7.0616,-0.4533/
      DATA(UV_TABLE(162,I),I=1,4)/-0.5000,1.1250,7.1523,-0.3893/
      DATA(UV_TABLE(163,I),I=1,4)/-0.5000,1.2500,7.2570,-0.3288/
      DATA(UV_TABLE(164,I),I=1,4)/-0.5000,1.3750,7.3691,-0.2705/
      DATA(UV_TABLE(165,I),I=1,4)/-0.5000,1.5000,7.4817,-0.2122/
      DATA(UV_TABLE(166,I),I=1,4)/-0.5000,1.6250,7.6053,-0.1577/
      DATA(UV_TABLE(167,I),I=1,4)/-0.5000,1.7500,7.7425,-0.1088/
      DATA(UV_TABLE(168,I),I=1,4)/-0.5000,1.8750,7.8912,-0.0654/
      DATA(UV_TABLE(169,I),I=1,4)/0.0000,-1.1250,5.2834,-1.9828/
      DATA(UV_TABLE(170,I),I=1,4)/0.0000,-1.0000,5.3635,-1.8736/
      DATA(UV_TABLE(171,I),I=1,4)/0.0000,-0.8750,5.4435,-1.7645/
      DATA(UV_TABLE(172,I),I=1,4)/0.0000,-0.7500,5.5236,-1.6553/
      DATA(UV_TABLE(173,I),I=1,4)/0.0000,-0.6250,5.6036,-1.5462/
      DATA(UV_TABLE(174,I),I=1,4)/0.0000,-0.5000,5.6836,-1.4370/
      DATA(UV_TABLE(175,I),I=1,4)/0.0000,-0.3750,5.7637,-1.3279/
      DATA(UV_TABLE(176,I),I=1,4)/0.0000,-0.2500,5.8445,-1.2188/
      DATA(UV_TABLE(177,I),I=1,4)/0.0000,-0.1250,5.9227,-1.1095/
      DATA(UV_TABLE(178,I),I=1,4)/0.0000,0.0000,6.0000,-1.0000/
      DATA(UV_TABLE(179,I),I=1,4)/0.0000,0.1250,6.0806,-0.8909/
      DATA(UV_TABLE(180,I),I=1,4)/0.0000,0.2500,6.1696,-0.7829/
      DATA(UV_TABLE(181,I),I=1,4)/0.0000,0.3000,6.2084,-0.7401/
      DATA(UV_TABLE(182,I),I=1,4)/0.0000,0.3750,6.2665,-0.6760/
      DATA(UV_TABLE(183,I),I=1,4)/0.0000,0.5000,6.3776,-0.5713/
      DATA(UV_TABLE(184,I),I=1,4)/0.0000,0.6000,6.4799,-0.4900/
      DATA(UV_TABLE(185,I),I=1,4)/0.0000,0.6250,6.5055,-0.4697/
      DATA(UV_TABLE(186,I),I=1,4)/0.0000,0.7500,6.6507,-0.3720/
      DATA(UV_TABLE(187,I),I=1,4)/0.0000,0.8000,6.7159,-0.3347/
      DATA(UV_TABLE(188,I),I=1,4)/0.0000,0.8750,6.8138,-0.2788/
      DATA(UV_TABLE(189,I),I=1,4)/0.0000,1.0000,6.9958,-0.1915/
      DATA(UV_TABLE(190,I),I=1,4)/0.0000,1.1250,7.1982,-0.1118/
      DATA(UV_TABLE(191,I),I=1,4)/0.0000,1.2500,7.4232,-0.0424/
      DATA(UV_TABLE(192,I),I=1,4)/0.0000,1.3750,7.6552,0.0232/
      DATA(UV_TABLE(193,I),I=1,4)/0.0000,1.5000,7.9034,0.0787/
      DATA(UV_TABLE(194,I),I=1,4)/0.0000,1.6250,8.1462,0.1379/
      DATA(UV_TABLE(195,I),I=1,4)/0.0000,1.7500,8.3886,0.1974/
      DATA(UV_TABLE(196,I),I=1,4)/0.0000,1.8750,8.6309,0.2570/
      DATA(UV_TABLE(197,I),I=1,4)/0.5000,-1.1250,4.3435,-2.2343/
      DATA(UV_TABLE(198,I),I=1,4)/0.5000,-1.0000,4.4683,-2.0974/
      DATA(UV_TABLE(199,I),I=1,4)/0.5000,-0.8750,4.5931,-1.9605/
      DATA(UV_TABLE(200,I),I=1,4)/0.5000,-0.7500,4.7193,-1.8238/
      DATA(UV_TABLE(201,I),I=1,4)/0.5000,-0.6250,4.8462,-1.6871/
      DATA(UV_TABLE(202,I),I=1,4)/0.5000,-0.5000,4.9725,-1.5504/
      DATA(UV_TABLE(203,I),I=1,4)/0.5000,-0.3750,5.0979,-1.4136/
      DATA(UV_TABLE(204,I),I=1,4)/0.5000,-0.2500,5.2232,-1.2767/
      DATA(UV_TABLE(205,I),I=1,4)/0.5000,-0.1250,5.3412,-1.1388/
      DATA(UV_TABLE(206,I),I=1,4)/0.5000,0.0000,5.4545,-1.0003/
      DATA(UV_TABLE(207,I),I=1,4)/0.5000,0.1250,5.5722,-0.8624/
      DATA(UV_TABLE(208,I),I=1,4)/0.5000,0.2500,5.7039,-0.7265/
      DATA(UV_TABLE(209,I),I=1,4)/0.5000,0.3000,5.7637,-0.6734/
      DATA(UV_TABLE(210,I),I=1,4)/0.5000,0.3750,5.8533,-0.5937/
      DATA(UV_TABLE(211,I),I=1,4)/0.5000,0.5000,6.0291,-0.4661/
      DATA(UV_TABLE(212,I),I=1,4)/0.5000,0.6000,6.1936,-0.3699/
      DATA(UV_TABLE(213,I),I=1,4)/0.5000,0.6250,6.2347,-0.3458/
      DATA(UV_TABLE(214,I),I=1,4)/0.5000,0.7500,6.4620,-0.2320/
      DATA(UV_TABLE(215,I),I=1,4)/0.5000,0.8000,6.5650,-0.1908/
      DATA(UV_TABLE(216,I),I=1,4)/0.5000,0.8750,6.7196,-0.1291/
      DATA(UV_TABLE(217,I),I=1,4)/0.5000,1.0000,6.9995,-0.0359/
      DATA(UV_TABLE(218,I),I=1,4)/0.5000,1.1250,7.3014,0.0454/
      DATA(UV_TABLE(219,I),I=1,4)/0.5000,1.2500,7.6175,0.1178/
      DATA(UV_TABLE(220,I),I=1,4)/0.5000,1.3750,7.9280,0.1941/
      DATA(UV_TABLE(221,I),I=1,4)/0.5000,1.5000,8.2398,0.2693/
      DATA(UV_TABLE(222,I),I=1,4)/0.5000,1.6250,8.5831,0.3181/
      DATA(UV_TABLE(223,I),I=1,4)/0.5000,1.7500,8.9265,0.3669/
      DATA(UV_TABLE(224,I),I=1,4)/0.5000,1.8750,9.2698,0.4157/
      DATA(UV_TABLE(225,I),I=1,4)/1.0000,-1.1250,3.4793,-2.5005/
      DATA(UV_TABLE(226,I),I=1,4)/1.0000,-1.0000,3.6529,-2.3366/
      DATA(UV_TABLE(227,I),I=1,4)/1.0000,-0.8750,3.8066,-2.1696/
      DATA(UV_TABLE(228,I),I=1,4)/1.0000,-0.7500,3.9604,-2.0026/
      DATA(UV_TABLE(229,I),I=1,4)/1.0000,-0.6250,4.1141,-1.8355/
      DATA(UV_TABLE(230,I),I=1,4)/1.0000,-0.5000,4.2668,-1.6683/
      DATA(UV_TABLE(231,I),I=1,4)/1.0000,-0.3750,4.4192,-1.5011/
      DATA(UV_TABLE(232,I),I=1,4)/1.0000,-0.2500,4.5717,-1.3339/
      DATA(UV_TABLE(233,I),I=1,4)/1.0000,-0.1250,4.7260,-1.1669/
      DATA(UV_TABLE(234,I),I=1,4)/1.0000,0.0000,4.8809,-1.0000/
      DATA(UV_TABLE(235,I),I=1,4)/1.0000,0.1250,5.0385,-0.8336/
      DATA(UV_TABLE(236,I),I=1,4)/1.0000,0.2500,5.2126,-0.6698/
      DATA(UV_TABLE(237,I),I=1,4)/1.0000,0.3000,5.2961,-0.6070/
      DATA(UV_TABLE(238,I),I=1,4)/1.0000,0.3750,5.4213,-0.5127/
      DATA(UV_TABLE(239,I),I=1,4)/1.0000,0.5000,5.6653,-0.3639/
      DATA(UV_TABLE(240,I),I=1,4)/1.0000,0.6000,5.8907,-0.2539/
      DATA(UV_TABLE(241,I),I=1,4)/1.0000,0.6250,5.9470,-0.2264/
      DATA(UV_TABLE(242,I),I=1,4)/1.0000,0.7500,6.2653,-0.1022/
      DATA(UV_TABLE(243,I),I=1,4)/1.0000,0.8000,6.4062,-0.0587/
      DATA(UV_TABLE(244,I),I=1,4)/1.0000,0.8750,6.6175,0.0065/
      DATA(UV_TABLE(245,I),I=1,4)/1.0000,1.0000,6.9997,0.0974/
      DATA(UV_TABLE(246,I),I=1,4)/1.0000,1.1250,7.3819,0.1885/
      DATA(UV_TABLE(247,I),I=1,4)/1.0000,1.2500,7.7714,0.2739/
      DATA(UV_TABLE(248,I),I=1,4)/1.0000,1.3750,8.1907,0.3335/
      DATA(UV_TABLE(249,I),I=1,4)/1.0000,1.5000,8.6100,0.3931/
      DATA(UV_TABLE(250,I),I=1,4)/1.0000,1.6250,9.0299,0.4520/
      DATA(UV_TABLE(251,I),I=1,4)/1.0000,1.7500,9.4614,0.4955/
      DATA(UV_TABLE(252,I),I=1,4)/1.0000,1.8750,9.8929,0.5390/
      DATA(UV_TABLE(253,I),I=1,4)/1.5000,-1.1250,2.6183,-2.8105/
      DATA(UV_TABLE(254,I),I=1,4)/1.5000,-1.0000,2.7925,-2.6074/
      DATA(UV_TABLE(255,I),I=1,4)/1.5000,-0.8750,2.9621,-2.4037/
      DATA(UV_TABLE(256,I),I=1,4)/1.5000,-0.7500,3.1373,-2.2008/
      DATA(UV_TABLE(257,I),I=1,4)/1.5000,-0.6250,3.3467,-2.0030/
      DATA(UV_TABLE(258,I),I=1,4)/1.5000,-0.5000,3.5561,-1.8053/
      DATA(UV_TABLE(259,I),I=1,4)/1.5000,-0.3750,3.7499,-1.6051/
      DATA(UV_TABLE(260,I),I=1,4)/1.5000,-0.2500,3.9354,-1.4035/
      DATA(UV_TABLE(261,I),I=1,4)/1.5000,-0.1250,4.1208,-1.2020/
      DATA(UV_TABLE(262,I),I=1,4)/1.5000,0.0000,4.3050,-1.0003/
      DATA(UV_TABLE(263,I),I=1,4)/1.5000,0.1250,4.4971,-0.7998/
      DATA(UV_TABLE(264,I),I=1,4)/1.5000,0.2500,4.7105,-0.6028/
      DATA(UV_TABLE(265,I),I=1,4)/1.5000,0.3000,4.8153,-0.5278/
      DATA(UV_TABLE(266,I),I=1,4)/1.5000,0.3750,4.9724,-0.4154/
      DATA(UV_TABLE(267,I),I=1,4)/1.5000,0.5000,5.2820,-0.2401/
      DATA(UV_TABLE(268,I),I=1,4)/1.5000,0.6000,5.5715,-0.1134/
      DATA(UV_TABLE(269,I),I=1,4)/1.5000,0.6250,5.6439,-0.0817/
      DATA(UV_TABLE(270,I),I=1,4)/1.5000,0.7500,6.0538,0.0568/
      DATA(UV_TABLE(271,I),I=1,4)/1.5000,0.8000,6.2386,0.1004/
      DATA(UV_TABLE(272,I),I=1,4)/1.5000,0.8750,6.5159,0.1657/
      DATA(UV_TABLE(273,I),I=1,4)/1.5000,1.0000,6.9915,0.2651/
      DATA(UV_TABLE(274,I),I=1,4)/1.5000,1.1250,7.4902,0.3439/
      DATA(UV_TABLE(275,I),I=1,4)/1.5000,1.2500,7.9960,0.4158/
      DATA(UV_TABLE(276,I),I=1,4)/1.5000,1.3750,8.5077,0.4800/
      DATA(UV_TABLE(277,I),I=1,4)/1.5000,1.5000,9.0283,0.5325/
      DATA(UV_TABLE(278,I),I=1,4)/1.5000,1.6250,9.5489,0.5850/
      DATA(UV_TABLE(279,I),I=1,4)/1.5000,1.7500,10.0695,0.6375/
      DATA(UV_TABLE(280,I),I=1,4)/1.5000,1.8750,10.5901,0.6900/
      DATA(UV_TABLE(281,I),I=1,4)/2.0000,-1.1250,1.7685,-3.0993/
      DATA(UV_TABLE(282,I),I=1,4)/2.0000,-1.0000,1.9879,-2.8663/
      DATA(UV_TABLE(283,I),I=1,4)/2.0000,-0.8750,2.2019,-2.6325/
      DATA(UV_TABLE(284,I),I=1,4)/2.0000,-0.7500,2.4159,-2.3986/
      DATA(UV_TABLE(285,I),I=1,4)/2.0000,-0.6250,2.6300,-2.1648/
      DATA(UV_TABLE(286,I),I=1,4)/2.0000,-0.5000,2.8280,-1.9288/
      DATA(UV_TABLE(287,I),I=1,4)/2.0000,-0.3750,3.0247,-1.6925/
      DATA(UV_TABLE(288,I),I=1,4)/2.0000,-0.2500,3.2547,-1.4613/
      DATA(UV_TABLE(289,I),I=1,4)/2.0000,-0.1250,3.4975,-1.2320/
      DATA(UV_TABLE(290,I),I=1,4)/2.0000,0.0000,3.7247,-1.0003/
      DATA(UV_TABLE(291,I),I=1,4)/2.0000,0.1250,3.9509,-0.7684/
      DATA(UV_TABLE(292,I),I=1,4)/2.0000,0.2500,4.2060,-0.5414/
      DATA(UV_TABLE(293,I),I=1,4)/2.0000,0.3000,4.3324,-0.4557/
      DATA(UV_TABLE(294,I),I=1,4)/2.0000,0.3750,4.5220,-0.3271/
      DATA(UV_TABLE(295,I),I=1,4)/2.0000,0.5000,4.9039,-0.1307/
      DATA(UV_TABLE(296,I),I=1,4)/2.0000,0.6000,5.2634,0.0069/
      DATA(UV_TABLE(297,I),I=1,4)/2.0000,0.6250,5.3533,0.0413/
      DATA(UV_TABLE(298,I),I=1,4)/2.0000,0.7500,5.8763,0.1754/
      DATA(UV_TABLE(299,I),I=1,4)/2.0000,0.8000,6.0963,0.2220/
      DATA(UV_TABLE(300,I),I=1,4)/2.0000,0.8750,6.4262,0.2919/
      DATA(UV_TABLE(301,I),I=1,4)/2.0000,1.0000,6.9989,0.3885/
      DATA(UV_TABLE(302,I),I=1,4)/2.0000,1.1250,7.5839,0.4736/
      DATA(UV_TABLE(303,I),I=1,4)/2.0000,1.2500,8.1844,0.5385/
      DATA(UV_TABLE(304,I),I=1,4)/2.0000,1.3750,8.7880,0.5994/
      DATA(UV_TABLE(305,I),I=1,4)/2.0000,1.5000,9.3916,0.6602/
      DATA(UV_TABLE(306,I),I=1,4)/2.0000,1.6250,9.9947,0.7219/
      DATA(UV_TABLE(307,I),I=1,4)/2.0000,1.7500,10.5940,0.7893/
      DATA(UV_TABLE(308,I),I=1,4)/2.0000,1.8750,11.1932,0.8568/
      DATA(UV_TABLE(309,I),I=1,4)/2.5000,-1.1250,0.7653,-3.3438/
      DATA(UV_TABLE(310,I),I=1,4)/2.5000,-1.0000,1.0759,-3.0922/
      DATA(UV_TABLE(311,I),I=1,4)/2.5000,-0.8750,1.3802,-2.8394/
      DATA(UV_TABLE(312,I),I=1,4)/2.5000,-0.7500,1.6765,-2.5851/
      DATA(UV_TABLE(313,I),I=1,4)/2.5000,-0.6250,1.9249,-2.3227/
      DATA(UV_TABLE(314,I),I=1,4)/2.5000,-0.5000,2.1692,-2.0596/
      DATA(UV_TABLE(315,I),I=1,4)/2.5000,-0.3750,2.4105,-1.7961/
      DATA(UV_TABLE(316,I),I=1,4)/2.5000,-0.2500,2.6489,-1.5322/
      DATA(UV_TABLE(317,I),I=1,4)/2.5000,-0.1250,2.8706,-1.2660/
      DATA(UV_TABLE(318,I),I=1,4)/2.5000,0.0000,3.0997,-1.0009/
      DATA(UV_TABLE(319,I),I=1,4)/2.5000,0.1250,3.3779,-0.7434/
      DATA(UV_TABLE(320,I),I=1,4)/2.5000,0.2500,3.6870,-0.4918/
      DATA(UV_TABLE(321,I),I=1,4)/2.5000,0.3000,3.8394,-0.3978/
      DATA(UV_TABLE(322,I),I=1,4)/2.5000,0.3750,4.0679,-0.2567/
      DATA(UV_TABLE(323,I),I=1,4)/2.5000,0.5000,4.5318,-0.0467/
      DATA(UV_TABLE(324,I),I=1,4)/2.5000,0.6000,4.9826,0.0855/
      DATA(UV_TABLE(325,I),I=1,4)/2.5000,0.6250,5.0953,0.1186/
      DATA(UV_TABLE(326,I),I=1,4)/2.5000,0.7500,5.7046,0.2563/
      DATA(UV_TABLE(327,I),I=1,4)/2.5000,0.8000,5.9585,0.3036/
      DATA(UV_TABLE(328,I),I=1,4)/2.5000,0.8750,6.3393,0.3746/
      DATA(UV_TABLE(329,I),I=1,4)/2.5000,1.0000,7.0010,0.4668/
      DATA(UV_TABLE(330,I),I=1,4)/2.5000,1.1250,7.6725,0.5470/
      DATA(UV_TABLE(331,I),I=1,4)/2.5000,1.2500,8.3528,0.6156/
      DATA(UV_TABLE(332,I),I=1,4)/2.5000,1.3750,9.0330,0.6842/
      DATA(UV_TABLE(333,I),I=1,4)/2.5000,1.5000,9.7101,0.7576/
      DATA(UV_TABLE(334,I),I=1,4)/2.5000,1.6250,10.3854,0.8336/
      DATA(UV_TABLE(335,I),I=1,4)/2.5000,1.7500,11.0607,0.9097/
      DATA(UV_TABLE(336,I),I=1,4)/2.5000,1.8750,11.7361,0.9857/
      DATA(UV_TABLE(337,I),I=1,4)/3.0000,-1.1250,-0.4383,-3.5372/
      DATA(UV_TABLE(338,I),I=1,4)/3.0000,-1.0000,-0.0810,-3.2611/
      DATA(UV_TABLE(339,I),I=1,4)/3.0000,-0.8750,0.2786,-2.9855/
      DATA(UV_TABLE(340,I),I=1,4)/3.0000,-0.7500,0.6481,-2.7120/
      DATA(UV_TABLE(341,I),I=1,4)/3.0000,-0.6250,1.0051,-2.4360/
      DATA(UV_TABLE(342,I),I=1,4)/3.0000,-0.5000,1.3422,-2.1558/
      DATA(UV_TABLE(343,I),I=1,4)/3.0000,-0.3750,1.6733,-1.8747/
      DATA(UV_TABLE(344,I),I=1,4)/3.0000,-0.2500,1.9485,-1.5840/
      DATA(UV_TABLE(345,I),I=1,4)/3.0000,-0.1250,2.2206,-1.2928/
      DATA(UV_TABLE(346,I),I=1,4)/3.0000,0.0000,2.4878,-1.0009/
      DATA(UV_TABLE(347,I),I=1,4)/3.0000,0.1250,2.7749,-0.7120/
      DATA(UV_TABLE(348,I),I=1,4)/3.0000,0.2500,3.1302,-0.4359/
      DATA(UV_TABLE(349,I),I=1,4)/3.0000,0.3000,3.3175,-0.3369/
      DATA(UV_TABLE(350,I),I=1,4)/3.0000,0.3750,3.5984,-0.1885/
      DATA(UV_TABLE(351,I),I=1,4)/3.0000,0.5000,4.1801,0.0168/
      DATA(UV_TABLE(352,I),I=1,4)/3.0000,0.6000,4.7145,0.1436/
      DATA(UV_TABLE(353,I),I=1,4)/3.0000,0.6250,4.8481,0.1753/
      DATA(UV_TABLE(354,I),I=1,4)/3.0000,0.7500,5.5346,0.3207/
      DATA(UV_TABLE(355,I),I=1,4)/3.0000,0.8000,5.8233,0.3661/
      DATA(UV_TABLE(356,I),I=1,4)/3.0000,0.8750,6.2564,0.4341/
      DATA(UV_TABLE(357,I),I=1,4)/3.0000,1.0000,6.9992,0.5244/
      DATA(UV_TABLE(358,I),I=1,4)/3.0000,1.1250,7.7515,0.6022/
      DATA(UV_TABLE(359,I),I=1,4)/3.0000,1.2500,8.5051,0.6782/
      DATA(UV_TABLE(360,I),I=1,4)/3.0000,1.3750,9.2552,0.7595/
      DATA(UV_TABLE(361,I),I=1,4)/3.0000,1.5000,10.0033,0.8438/
      DATA(UV_TABLE(362,I),I=1,4)/3.0000,1.6250,10.7514,0.9280/
      DATA(UV_TABLE(363,I),I=1,4)/3.0000,1.7500,11.4995,1.0122/
      DATA(UV_TABLE(364,I),I=1,4)/3.0000,1.8750,12.2406,1.1056/
      DATA(UV_TABLE(365,I),I=1,4)/3.5000,-1.1250,-1.6399,-3.7841/
      DATA(UV_TABLE(366,I),I=1,4)/3.5000,-1.0000,-1.2417,-3.4764/
      DATA(UV_TABLE(367,I),I=1,4)/3.5000,-0.8750,-0.8434,-3.1687/
      DATA(UV_TABLE(368,I),I=1,4)/3.5000,-0.7500,-0.4452,-2.8610/
      DATA(UV_TABLE(369,I),I=1,4)/3.5000,-0.6250,-0.0469,-2.5533/
      DATA(UV_TABLE(370,I),I=1,4)/3.5000,-0.5000,0.3580,-2.2470/
      DATA(UV_TABLE(371,I),I=1,4)/3.5000,-0.3750,0.7698,-1.9421/
      DATA(UV_TABLE(372,I),I=1,4)/3.5000,-0.2500,1.1577,-1.6324/
      DATA(UV_TABLE(373,I),I=1,4)/3.5000,-0.1250,1.5334,-1.3202/
      DATA(UV_TABLE(374,I),I=1,4)/3.5000,0.0000,1.8633,-1.0001/
      DATA(UV_TABLE(375,I),I=1,4)/3.5000,0.1250,2.1815,-0.6779/
      DATA(UV_TABLE(376,I),I=1,4)/3.5000,0.2500,2.5740,-0.3695/
      DATA(UV_TABLE(377,I),I=1,4)/3.5000,0.3000,2.7859,-0.2603/
      DATA(UV_TABLE(378,I),I=1,4)/3.5000,0.3750,3.1037,-0.0964/
      DATA(UV_TABLE(379,I),I=1,4)/3.5000,0.5000,3.7783,0.1192/
      DATA(UV_TABLE(380,I),I=1,4)/3.5000,0.6000,4.3875,0.2514/
      DATA(UV_TABLE(381,I),I=1,4)/3.5000,0.6250,4.5398,0.2844/
      DATA(UV_TABLE(382,I),I=1,4)/3.5000,0.7500,5.3275,0.4270/
      DATA(UV_TABLE(383,I),I=1,4)/3.5000,0.8000,5.6581,0.4683/
      DATA(UV_TABLE(384,I),I=1,4)/3.5000,0.8750,6.1541,0.5303/
      DATA(UV_TABLE(385,I),I=1,4)/3.5000,1.0000,6.9964,0.6107/
      DATA(UV_TABLE(386,I),I=1,4)/3.5000,1.1250,7.8370,0.6942/
      DATA(UV_TABLE(387,I),I=1,4)/3.5000,1.2500,8.6710,0.7878/
      DATA(UV_TABLE(388,I),I=1,4)/3.5000,1.3750,9.5048,0.8817/
      DATA(UV_TABLE(389,I),I=1,4)/3.5000,1.5000,10.3387,0.9756/
      DATA(UV_TABLE(390,I),I=1,4)/3.5000,1.6250,11.1619,1.0835/
      DATA(UV_TABLE(391,I),I=1,4)/3.5000,1.7500,11.9847,1.1918/
      DATA(UV_TABLE(392,I),I=1,4)/3.5000,1.8750,12.8076,1.3000/
      DATA(UV_TABLE(393,I),I=1,4)/4.0000,-1.1250,-3.5191,-4.1607/
      DATA(UV_TABLE(394,I),I=1,4)/4.0000,-1.0000,-3.0642,-3.8091/
      DATA(UV_TABLE(395,I),I=1,4)/4.0000,-0.8750,-2.6092,-3.4575/
      DATA(UV_TABLE(396,I),I=1,4)/4.0000,-0.7500,-2.1542,-3.1060/
      DATA(UV_TABLE(397,I),I=1,4)/4.0000,-0.6250,-1.6992,-2.7544/
      DATA(UV_TABLE(398,I),I=1,4)/4.0000,-0.5000,-1.2442,-2.4029/
      DATA(UV_TABLE(399,I),I=1,4)/4.0000,-0.3750,-0.7892,-2.0513/
      DATA(UV_TABLE(400,I),I=1,4)/4.0000,-0.2500,-0.3342,-1.6997/
      DATA(UV_TABLE(401,I),I=1,4)/4.0000,-0.1250,0.1224,-1.3485/
      DATA(UV_TABLE(402,I),I=1,4)/4.0000,0.0000,0.5929,-1.0002/
      DATA(UV_TABLE(403,I),I=1,4)/4.0000,0.1250,1.1150,-0.6643/
      DATA(UV_TABLE(404,I),I=1,4)/4.0000,0.2500,1.7497,-0.3610/
      DATA(UV_TABLE(405,I),I=1,4)/4.0000,0.3000,2.0385,-0.2529/
      DATA(UV_TABLE(406,I),I=1,4)/4.0000,0.3750,2.4718,-0.0907/
      DATA(UV_TABLE(407,I),I=1,4)/4.0000,0.5000,3.2888,0.1312/
      DATA(UV_TABLE(408,I),I=1,4)/4.0000,0.6000,3.9879,0.2794/
      DATA(UV_TABLE(409,I),I=1,4)/4.0000,0.6250,4.1627,0.3164/
      DATA(UV_TABLE(410,I),I=1,4)/4.0000,0.7500,5.0877,0.4557/
      DATA(UV_TABLE(411,I),I=1,4)/4.0000,0.8000,5.4662,0.5019/
      DATA(UV_TABLE(412,I),I=1,4)/4.0000,0.8750,6.0339,0.5712/
      DATA(UV_TABLE(413,I),I=1,4)/4.0000,1.0000,6.9981,0.6604/
      DATA(UV_TABLE(414,I),I=1,4)/4.0000,1.1250,7.9628,0.7469/
      DATA(UV_TABLE(415,I),I=1,4)/4.0000,1.2500,8.9154,0.8542/
      DATA(UV_TABLE(416,I),I=1,4)/4.0000,1.3750,9.8669,0.9629/
      DATA(UV_TABLE(417,I),I=1,4)/4.0000,1.5000,10.8071,1.0866/
      DATA(UV_TABLE(418,I),I=1,4)/4.0000,1.6250,11.7472,1.2103/
      DATA(UV_TABLE(419,I),I=1,4)/4.0000,1.7500,12.6873,1.3340/
      DATA(UV_TABLE(420,I),I=1,4)/4.0000,1.8750,13.6275,1.4577/
      DATA(UV_TABLE(421,I),I=1,4)/4.5000,-1.1250,-5.3705,-4.6510/
      DATA(UV_TABLE(422,I),I=1,4)/4.5000,-1.0000,-4.8456,-4.2454/
      DATA(UV_TABLE(423,I),I=1,4)/4.5000,-0.8750,-4.3206,-3.8398/
      DATA(UV_TABLE(424,I),I=1,4)/4.5000,-0.7500,-3.7957,-3.4341/
      DATA(UV_TABLE(425,I),I=1,4)/4.5000,-0.6250,-3.2707,-3.0285/
      DATA(UV_TABLE(426,I),I=1,4)/4.5000,-0.5000,-2.7457,-2.6229/
      DATA(UV_TABLE(427,I),I=1,4)/4.5000,-0.3750,-2.2208,-2.2173/
      DATA(UV_TABLE(428,I),I=1,4)/4.5000,-0.2500,-1.6958,-1.8117/
      DATA(UV_TABLE(429,I),I=1,4)/4.5000,-0.1250,-1.1709,-1.4060/
      DATA(UV_TABLE(430,I),I=1,4)/4.5000,0.0000,-0.6459,-1.0004/
      DATA(UV_TABLE(431,I),I=1,4)/4.5000,0.1250,-0.0495,-0.6118/
      DATA(UV_TABLE(432,I),I=1,4)/4.5000,0.2500,0.7269,-0.2778/
      DATA(UV_TABLE(433,I),I=1,4)/4.5000,0.3000,1.0909,-0.1679/
      DATA(UV_TABLE(434,I),I=1,4)/4.5000,0.3750,1.6370,-0.0030/
      DATA(UV_TABLE(435,I),I=1,4)/4.5000,0.5000,2.6653,0.1922/
      DATA(UV_TABLE(436,I),I=1,4)/4.5000,0.6000,3.5001,0.3392/
      DATA(UV_TABLE(437,I),I=1,4)/4.5000,0.6250,3.7088,0.3759/
      DATA(UV_TABLE(438,I),I=1,4)/4.5000,0.7500,4.7870,0.5259/
      DATA(UV_TABLE(439,I),I=1,4)/4.5000,0.8000,5.2264,0.5753/
      DATA(UV_TABLE(440,I),I=1,4)/4.5000,0.8750,5.8856,0.6493/
      DATA(UV_TABLE(441,I),I=1,4)/4.5000,1.0000,6.9999,0.7481/
      DATA(UV_TABLE(442,I),I=1,4)/4.5000,1.1250,8.1002,0.8700/
      DATA(UV_TABLE(443,I),I=1,4)/4.5000,1.2500,9.1868,1.0103/
      DATA(UV_TABLE(444,I),I=1,4)/4.5000,1.3750,10.2715,1.1530/
      DATA(UV_TABLE(445,I),I=1,4)/4.5000,1.5000,11.3562,1.2957/
      DATA(UV_TABLE(446,I),I=1,4)/4.5000,1.6250,12.4409,1.4384/
      DATA(UV_TABLE(447,I),I=1,4)/4.5000,1.7500,13.5256,1.5811/
      DATA(UV_TABLE(448,I),I=1,4)/4.5000,1.8750,14.6103,1.7238/
      DATA(UV_TABLE(449,I),I=1,4)/4.8000,-1.1250,-7.5873,-5.2520/
      DATA(UV_TABLE(450,I),I=1,4)/4.8000,-1.0000,-6.9752,-4.7790/
      DATA(UV_TABLE(451,I),I=1,4)/4.8000,-0.8750,-6.3629,-4.3060/
      DATA(UV_TABLE(452,I),I=1,4)/4.8000,-0.7500,-5.7508,-3.8329/
      DATA(UV_TABLE(453,I),I=1,4)/4.8000,-0.6250,-5.1386,-3.3599/
      DATA(UV_TABLE(454,I),I=1,4)/4.8000,-0.5000,-4.5264,-2.8869/
      DATA(UV_TABLE(455,I),I=1,4)/4.8000,-0.3750,-3.9142,-2.4139/
      DATA(UV_TABLE(456,I),I=1,4)/4.8000,-0.2500,-3.3020,-1.9409/
      DATA(UV_TABLE(457,I),I=1,4)/4.8000,-0.1250,-2.6852,-1.4688/
      DATA(UV_TABLE(458,I),I=1,4)/4.8000,0.0000,-2.0506,-1.0006/
      DATA(UV_TABLE(459,I),I=1,4)/4.8000,0.1250,-1.3874,-0.5392/
      DATA(UV_TABLE(460,I),I=1,4)/4.8000,0.2500,-0.4631,-0.1631/
      DATA(UV_TABLE(461,I),I=1,4)/4.8000,0.3000,-0.0184,-0.0486/
      DATA(UV_TABLE(462,I),I=1,4)/4.8000,0.3750,0.6487,0.1230/
      DATA(UV_TABLE(463,I),I=1,4)/4.8000,0.5000,1.8945,0.3052/
      DATA(UV_TABLE(464,I),I=1,4)/4.8000,0.6000,2.9041,0.4349/
      DATA(UV_TABLE(465,I),I=1,4)/4.8000,0.6250,3.1565,0.4673/
      DATA(UV_TABLE(466,I),I=1,4)/4.8000,0.7500,4.4274,0.6221/
      DATA(UV_TABLE(467,I),I=1,4)/4.8000,0.8000,4.9363,0.6837/
      DATA(UV_TABLE(468,I),I=1,4)/4.8000,0.8750,5.6997,0.7760/
      DATA(UV_TABLE(469,I),I=1,4)/4.8000,1.0000,6.9935,0.8968/
      DATA(UV_TABLE(470,I),I=1,4)/4.8000,1.1250,8.2647,1.0549/
      DATA(UV_TABLE(471,I),I=1,4)/4.8000,1.2500,9.5304,1.2203/
      DATA(UV_TABLE(472,I),I=1,4)/4.8000,1.3750,10.7953,1.3868/
      DATA(UV_TABLE(473,I),I=1,4)/4.8000,1.5000,12.0602,1.5532/
      DATA(UV_TABLE(474,I),I=1,4)/4.8000,1.6250,13.3252,1.7196/
      DATA(UV_TABLE(475,I),I=1,4)/4.8000,1.7500,14.5902,1.8860/
      DATA(UV_TABLE(476,I),I=1,4)/4.8000,1.8750,15.8551,2.0525/
      DATA(UV_TABLE(477,I),I=1,4)/5.0500,-1.1250,-9.3350,-5.7313/
      DATA(UV_TABLE(478,I),I=1,4)/5.0500,-1.0000,-8.6532,-5.2045/
      DATA(UV_TABLE(479,I),I=1,4)/5.0500,-0.8750,-7.9714,-4.6777/
      DATA(UV_TABLE(480,I),I=1,4)/5.0500,-0.7500,-7.2896,-4.1509/
      DATA(UV_TABLE(481,I),I=1,4)/5.0500,-0.6250,-6.6078,-3.6241/
      DATA(UV_TABLE(482,I),I=1,4)/5.0500,-0.5000,-5.9260,-3.0973/
      DATA(UV_TABLE(483,I),I=1,4)/5.0500,-0.3750,-5.2442,-2.5705/
      DATA(UV_TABLE(484,I),I=1,4)/5.0500,-0.2500,-4.5621,-2.0438/
      DATA(UV_TABLE(485,I),I=1,4)/5.0500,-0.1250,-3.8690,-1.5194/
      DATA(UV_TABLE(486,I),I=1,4)/5.0500,0.0000,-3.1492,-1.0007/
      DATA(UV_TABLE(487,I),I=1,4)/5.0500,0.1250,-2.4295,-0.4821/
      DATA(UV_TABLE(488,I),I=1,4)/5.0500,0.2500,-1.3814,-0.0735/
      DATA(UV_TABLE(489,I),I=1,4)/5.0500,0.3000,-0.8786,0.0488/
      DATA(UV_TABLE(490,I),I=1,4)/5.0500,0.3750,-0.1244,0.2322/
      DATA(UV_TABLE(491,I),I=1,4)/5.0500,0.5000,1.2865,0.4125/
      DATA(UV_TABLE(492,I),I=1,4)/5.0500,0.6000,2.4307,0.5344/
      DATA(UV_TABLE(493,I),I=1,4)/5.0500,0.6250,2.7167,0.5649/
      DATA(UV_TABLE(494,I),I=1,4)/5.0500,0.7500,4.1411,0.7239/
      DATA(UV_TABLE(495,I),I=1,4)/5.0500,0.8000,4.7059,0.7948/
      DATA(UV_TABLE(496,I),I=1,4)/5.0500,0.8750,5.5625,0.8609/
      DATA(UV_TABLE(497,I),I=1,4)/5.0500,1.0000,6.9903,0.9711/
      DATA(UV_TABLE(498,I),I=1,4)/5.0500,1.1250,8.3989,1.2249/
      DATA(UV_TABLE(499,I),I=1,4)/5.0500,1.2500,9.8077,1.4103/
      DATA(UV_TABLE(500,I),I=1,4)/5.0500,1.3750,11.2164,1.5957/
      DATA(UV_TABLE(501,I),I=1,4)/5.0500,1.5000,12.6252,1.7810/
      DATA(UV_TABLE(502,I),I=1,4)/5.0500,1.6250,14.0340,1.9664/
      DATA(UV_TABLE(503,I),I=1,4)/5.0500,1.7500,15.4428,2.1517/
      DATA(UV_TABLE(504,I),I=1,4)/5.0500,1.8750,16.8516,2.3371/
      DATA(UV_TABLE(505,I),I=1,4)/5.2000,-1.1250,-10.1443,-5.9672/
      DATA(UV_TABLE(506,I),I=1,4)/5.2000,-1.0000,-9.4282,-5.4139/
      DATA(UV_TABLE(507,I),I=1,4)/5.2000,-0.8750,-8.7120,-4.8605/
      DATA(UV_TABLE(508,I),I=1,4)/5.2000,-0.7500,-7.9959,-4.3072/
      DATA(UV_TABLE(509,I),I=1,4)/5.2000,-0.6250,-7.2797,-3.7539/
      DATA(UV_TABLE(510,I),I=1,4)/5.2000,-0.5000,-6.5636,-3.2005/
      DATA(UV_TABLE(511,I),I=1,4)/5.2000,-0.3750,-5.8474,-2.6472/
      DATA(UV_TABLE(512,I),I=1,4)/5.2000,-0.2500,-5.1300,-2.0941/
      DATA(UV_TABLE(513,I),I=1,4)/5.2000,-0.1250,-4.3917,-1.5455/
      DATA(UV_TABLE(514,I),I=1,4)/5.2000,0.0000,-3.6357,-1.0006/
      DATA(UV_TABLE(515,I),I=1,4)/5.2000,0.1250,-2.8797,-0.4559/
      DATA(UV_TABLE(516,I),I=1,4)/5.2000,0.2500,-1.7562,-0.0344/
      DATA(UV_TABLE(517,I),I=1,4)/5.2000,0.3000,-1.2406,0.1026/
      DATA(UV_TABLE(518,I),I=1,4)/5.2000,0.3750,-0.5394,0.3238/
      DATA(UV_TABLE(519,I),I=1,4)/5.2000,0.5000,0.7776,0.5852/
      DATA(UV_TABLE(520,I),I=1,4)/5.2000,0.6000,2.0197,0.7278/
      DATA(UV_TABLE(521,I),I=1,4)/5.2000,0.6250,2.3416,0.7508/
      DATA(UV_TABLE(522,I),I=1,4)/5.2000,0.7500,3.9512,0.8655/
      DATA(UV_TABLE(523,I),I=1,4)/5.2000,0.8000,4.5951,0.9114/
      DATA(UV_TABLE(524,I),I=1,4)/5.2000,0.8750,5.4932,0.9595/
      DATA(UV_TABLE(525,I),I=1,4)/5.2000,1.0000,6.9901,1.0396/
      DATA(UV_TABLE(526,I),I=1,4)/5.2000,1.1250,8.4727,1.3655/
      DATA(UV_TABLE(527,I),I=1,4)/5.2000,1.2500,9.9524,1.5601/
      DATA(UV_TABLE(528,I),I=1,4)/5.2000,1.3750,11.4322,1.7549/
      DATA(UV_TABLE(529,I),I=1,4)/5.2000,1.5000,12.9120,1.9496/
      DATA(UV_TABLE(530,I),I=1,4)/5.2000,1.6250,14.3917,2.1443/
      DATA(UV_TABLE(531,I),I=1,4)/5.2000,1.7500,15.8715,2.3389/
      DATA(UV_TABLE(532,I),I=1,4)/5.2000,1.8750,17.3512,2.5337/
      DATA(UV_TABLE(533,I),I=1,4)/5.3200,-1.1250,-10.7918,-6.1559/
      DATA(UV_TABLE(534,I),I=1,4)/5.3200,-1.0000,-10.0482,-5.5814/
      DATA(UV_TABLE(535,I),I=1,4)/5.3200,-0.8750,-9.3045,-5.0068/
      DATA(UV_TABLE(536,I),I=1,4)/5.3200,-0.7500,-8.5610,-4.4322/
      DATA(UV_TABLE(537,I),I=1,4)/5.3200,-0.6250,-7.8173,-3.8576/
      DATA(UV_TABLE(538,I),I=1,4)/5.3200,-0.5000,-7.0737,-3.2830/
      DATA(UV_TABLE(539,I),I=1,4)/5.3200,-0.3750,-6.3300,-2.7085/
      DATA(UV_TABLE(540,I),I=1,4)/5.3200,-0.2500,-5.5844,-2.1343/
      DATA(UV_TABLE(541,I),I=1,4)/5.3200,-0.1250,-4.8100,-1.5663/
      DATA(UV_TABLE(542,I),I=1,4)/5.3200,0.0000,-4.0249,-1.0006/
      DATA(UV_TABLE(543,I),I=1,4)/5.3200,0.1250,-3.2399,-0.4349/
      DATA(UV_TABLE(544,I),I=1,4)/5.3200,0.2500,-2.0560,-0.0030/
      DATA(UV_TABLE(545,I),I=1,4)/5.3200,0.3000,-1.5303,0.1456/
      DATA(UV_TABLE(546,I),I=1,4)/5.3200,0.3750,-0.9000,0.4034/
      DATA(UV_TABLE(547,I),I=1,4)/5.3200,0.5000,0.5325,0.6684/
      DATA(UV_TABLE(548,I),I=1,4)/5.3200,0.6000,1.8218,0.8209/
      DATA(UV_TABLE(549,I),I=1,4)/5.3200,0.6250,2.1574,0.8439/
      DATA(UV_TABLE(550,I),I=1,4)/5.3200,0.7500,3.8353,0.9587/
      DATA(UV_TABLE(551,I),I=1,4)/5.3200,0.8000,4.5065,1.0047/
      DATA(UV_TABLE(552,I),I=1,4)/5.3200,0.8750,5.4389,1.0670/
      DATA(UV_TABLE(553,I),I=1,4)/5.3200,1.0000,6.9930,1.1708/
      DATA(UV_TABLE(554,I),I=1,4)/5.3200,1.1250,8.5317,1.4779/
      DATA(UV_TABLE(555,I),I=1,4)/5.3200,1.2500,10.0683,1.6800/
      DATA(UV_TABLE(556,I),I=1,4)/5.3200,1.3750,11.6048,1.8822/
      DATA(UV_TABLE(557,I),I=1,4)/5.3200,1.5000,13.1414,2.0844/
      DATA(UV_TABLE(558,I),I=1,4)/5.3200,1.6250,14.6779,2.2866/
      DATA(UV_TABLE(559,I),I=1,4)/5.3200,1.7500,16.2145,2.4887/
      DATA(UV_TABLE(560,I),I=1,4)/5.3200,1.8750,17.7510,2.6909/
      DATA(UV_TABLE(561,I),I=1,4)/5.4500,-1.1250,-11.4932,-6.3604/
      DATA(UV_TABLE(562,I),I=1,4)/5.4500,-1.0000,-10.7199,-5.7628/
      DATA(UV_TABLE(563,I),I=1,4)/5.4500,-0.8750,-9.9464,-5.1652/
      DATA(UV_TABLE(564,I),I=1,4)/5.4500,-0.7500,-9.1731,-4.5676/
      DATA(UV_TABLE(565,I),I=1,4)/5.4500,-0.6250,-8.3996,-3.9701/
      DATA(UV_TABLE(566,I),I=1,4)/5.4500,-0.5000,-7.6263,-3.3724/
      DATA(UV_TABLE(567,I),I=1,4)/5.4500,-0.3750,-6.8528,-2.7749/
      DATA(UV_TABLE(568,I),I=1,4)/5.4500,-0.2500,-6.0766,-2.1779/
      DATA(UV_TABLE(569,I),I=1,4)/5.4500,-0.1250,-5.2630,-1.5889/
      DATA(UV_TABLE(570,I),I=1,4)/5.4500,0.0000,-4.4465,-1.0005/
      DATA(UV_TABLE(571,I),I=1,4)/5.4500,0.1250,-3.6301,-0.4122/
      DATA(UV_TABLE(572,I),I=1,4)/5.4500,0.2500,-2.3808,0.0309/
      DATA(UV_TABLE(573,I),I=1,4)/5.4500,0.3000,-1.8440,0.1922/
      DATA(UV_TABLE(574,I),I=1,4)/5.4500,0.3750,-1.1532,0.4593/
      DATA(UV_TABLE(575,I),I=1,4)/5.4500,0.5000,0.3057,0.7453/
      DATA(UV_TABLE(576,I),I=1,4)/5.4500,0.6000,1.6388,0.9058/
      DATA(UV_TABLE(577,I),I=1,4)/5.4500,0.6250,1.9853,0.9308/
      DATA(UV_TABLE(578,I),I=1,4)/5.4500,0.7500,3.7176,1.0558/
      DATA(UV_TABLE(579,I),I=1,4)/5.4500,0.8000,4.4105,1.1058/
      DATA(UV_TABLE(580,I),I=1,4)/5.4500,0.8750,5.3798,1.1696/
      DATA(UV_TABLE(581,I),I=1,4)/5.4500,1.0000,6.9952,1.2758/
      DATA(UV_TABLE(582,I),I=1,4)/5.4500,1.1250,8.5957,1.5997/
      DATA(UV_TABLE(583,I),I=1,4)/5.4500,1.2500,10.1937,1.8099/
      DATA(UV_TABLE(584,I),I=1,4)/5.4500,1.3750,11.7918,2.0202/
      DATA(UV_TABLE(585,I),I=1,4)/5.4500,1.5000,13.3899,2.2305/
      DATA(UV_TABLE(586,I),I=1,4)/5.4500,1.6250,14.9880,2.4407/
      DATA(UV_TABLE(587,I),I=1,4)/5.4500,1.7500,16.5860,2.6510/
      DATA(UV_TABLE(588,I),I=1,4)/5.4500,1.8750,18.1840,2.8613/
      DATA(UV_TABLE(589,I),I=1,4)/5.7000,-1.1250,-12.8381,-6.7519/
      DATA(UV_TABLE(590,I),I=1,4)/5.7000,-1.0000,-12.0077,-6.1103/
      DATA(UV_TABLE(591,I),I=1,4)/5.7000,-0.8750,-11.1773,-5.4686/
      DATA(UV_TABLE(592,I),I=1,4)/5.7000,-0.7500,-10.3469,-4.8270/
      DATA(UV_TABLE(593,I),I=1,4)/5.7000,-0.6250,-9.5164,-4.1853/
      DATA(UV_TABLE(594,I),I=1,4)/5.7000,-0.5000,-8.6860,-3.5437/
      DATA(UV_TABLE(595,I),I=1,4)/5.7000,-0.3750,-7.8555,-2.9020/
      DATA(UV_TABLE(596,I),I=1,4)/5.7000,-0.2500,-7.0087,-2.2640/
      DATA(UV_TABLE(597,I),I=1,4)/5.7000,-0.1250,-6.1319,-1.6322/
      DATA(UV_TABLE(598,I),I=1,4)/5.7000,0.0000,-5.2553,-1.0005/
      DATA(UV_TABLE(599,I),I=1,4)/5.7000,0.1250,-4.3780,-0.3690/
      DATA(UV_TABLE(600,I),I=1,4)/5.7000,0.2500,-3.0062,0.0967/
      DATA(UV_TABLE(601,I),I=1,4)/5.7000,0.3000,-2.4388,0.2752/
      DATA(UV_TABLE(602,I),I=1,4)/5.7000,0.3750,-1.5878,0.5431/
      DATA(UV_TABLE(603,I),I=1,4)/5.7000,0.5000,0.0594,0.8286/
      DATA(UV_TABLE(604,I),I=1,4)/5.7000,0.6000,1.4404,0.9971/
      DATA(UV_TABLE(605,I),I=1,4)/5.7000,0.6250,1.7857,1.0392/
      DATA(UV_TABLE(606,I),I=1,4)/5.7000,0.7500,3.5314,1.2181/
      DATA(UV_TABLE(607,I),I=1,4)/5.7000,0.8000,4.2202,1.3031/
      DATA(UV_TABLE(608,I),I=1,4)/5.7000,0.8750,5.2534,1.4306/
      DATA(UV_TABLE(609,I),I=1,4)/5.7000,1.0000,6.9975,1.6088/
      DATA(UV_TABLE(610,I),I=1,4)/5.7000,1.1250,8.7134,1.8346/
      DATA(UV_TABLE(611,I),I=1,4)/5.7000,1.2500,10.4292,2.0603/
      DATA(UV_TABLE(612,I),I=1,4)/5.7000,1.3750,12.1451,2.2861/
      DATA(UV_TABLE(613,I),I=1,4)/5.7000,1.5000,13.8610,2.5119/
      DATA(UV_TABLE(614,I),I=1,4)/5.7000,1.6250,15.5769,2.7376/
      DATA(UV_TABLE(615,I),I=1,4)/5.7000,1.7500,17.2927,2.9634/
      DATA(UV_TABLE(616,I),I=1,4)/5.7000,1.8750,19.0086,3.1892/
      DATA(UV_TABLE(617,I),I=1,4)/6.0000,-1.1250,-14.4507,-7.2213/
      DATA(UV_TABLE(618,I),I=1,4)/6.0000,-1.0000,-13.5519,-6.5268/
      DATA(UV_TABLE(619,I),I=1,4)/6.0000,-0.8750,-12.6532,-5.8324/
      DATA(UV_TABLE(620,I),I=1,4)/6.0000,-0.7500,-11.7544,-5.1379/
      DATA(UV_TABLE(621,I),I=1,4)/6.0000,-0.6250,-10.8556,-4.4434/
      DATA(UV_TABLE(622,I),I=1,4)/6.0000,-0.5000,-9.9568,-3.7490/
      DATA(UV_TABLE(623,I),I=1,4)/6.0000,-0.3750,-9.0580,-3.0545/
      DATA(UV_TABLE(624,I),I=1,4)/6.0000,-0.2500,-8.1228,-2.3679/
      DATA(UV_TABLE(625,I),I=1,4)/6.0000,-0.1250,-7.1739,-1.6842/
      DATA(UV_TABLE(626,I),I=1,4)/6.0000,0.0000,-6.2251,-1.0004/
      DATA(UV_TABLE(627,I),I=1,4)/6.0000,0.1250,-5.2748,-0.3172/
      DATA(UV_TABLE(628,I),I=1,4)/6.0000,0.2500,-3.7570,0.1757/
      DATA(UV_TABLE(629,I),I=1,4)/6.0000,0.3000,-3.1500,0.3729/
      DATA(UV_TABLE(630,I),I=1,4)/6.0000,0.3750,-2.2396,0.6688/
      DATA(UV_TABLE(631,I),I=1,4)/6.0000,0.5000,-0.5090,1.0210/
      DATA(UV_TABLE(632,I),I=1,4)/6.0000,0.6000,0.9825,1.2077/
      DATA(UV_TABLE(633,I),I=1,4)/6.0000,0.6250,1.3554,1.2544/
      DATA(UV_TABLE(634,I),I=1,4)/6.0000,0.7500,3.2439,1.4499/
      DATA(UV_TABLE(635,I),I=1,4)/6.0000,0.8000,3.9902,1.5409/
      DATA(UV_TABLE(636,I),I=1,4)/6.0000,0.8750,5.1096,1.6773/
      DATA(UV_TABLE(637,I),I=1,4)/6.0000,1.0000,6.9958,1.8723/
      DATA(UV_TABLE(638,I),I=1,4)/6.0000,1.1250,8.8530,2.1167/
      DATA(UV_TABLE(639,I),I=1,4)/6.0000,1.2500,10.7101,2.3610/
      DATA(UV_TABLE(640,I),I=1,4)/6.0000,1.3750,12.5672,2.6054/
      DATA(UV_TABLE(641,I),I=1,4)/6.0000,1.5000,14.4243,2.8497/
      DATA(UV_TABLE(642,I),I=1,4)/6.0000,1.6250,16.2814,3.0941/
      DATA(UV_TABLE(643,I),I=1,4)/6.0000,1.7500,18.1385,3.3384/
      DATA(UV_TABLE(644,I),I=1,4)/6.0000,1.8750,19.9956,3.5828/
      DATA(UV_TABLE(645,I),I=1,4)/6.5000,-1.1250,-17.1409,-8.0050/
      DATA(UV_TABLE(646,I),I=1,4)/6.5000,-1.0000,-16.1280,-7.2224/
      DATA(UV_TABLE(647,I),I=1,4)/6.5000,-0.8750,-15.1152,-6.4398/
      DATA(UV_TABLE(648,I),I=1,4)/6.5000,-0.7500,-14.1023,-5.6572/
      DATA(UV_TABLE(649,I),I=1,4)/6.5000,-0.6250,-13.0895,-4.8746/
      DATA(UV_TABLE(650,I),I=1,4)/6.5000,-0.5000,-12.0766,-4.0920/
      DATA(UV_TABLE(651,I),I=1,4)/6.5000,-0.3750,-11.0510,-3.3121/
      DATA(UV_TABLE(652,I),I=1,4)/6.5000,-0.2500,-9.9817,-2.5417/
      DATA(UV_TABLE(653,I),I=1,4)/6.5000,-0.1250,-8.9124,-1.7712/
      DATA(UV_TABLE(654,I),I=1,4)/6.5000,0.0000,-7.8431,-1.0007/
      DATA(UV_TABLE(655,I),I=1,4)/6.5000,0.1250,-6.7220,-0.2475/
      DATA(UV_TABLE(656,I),I=1,4)/6.5000,0.2500,-5.0117,0.3079/
      DATA(UV_TABLE(657,I),I=1,4)/6.5000,0.3000,-4.3275,0.5301/
      DATA(UV_TABLE(658,I),I=1,4)/6.5000,0.3750,-3.3013,0.8634/
      DATA(UV_TABLE(659,I),I=1,4)/6.5000,0.5000,-1.4502,1.3316/
      DATA(UV_TABLE(660,I),I=1,4)/6.5000,0.6000,0.2232,1.5509/
      DATA(UV_TABLE(661,I),I=1,4)/6.5000,0.6250,0.6416,1.6057/
      DATA(UV_TABLE(662,I),I=1,4)/6.5000,0.7500,2.7682,1.8295/
      DATA(UV_TABLE(663,I),I=1,4)/6.5000,0.8000,3.6083,1.9332/
      DATA(UV_TABLE(664,I),I=1,4)/6.5000,0.8750,4.8685,2.0888/
      DATA(UV_TABLE(665,I),I=1,4)/6.5000,1.0000,6.9937,2.3096/
      DATA(UV_TABLE(666,I),I=1,4)/6.5000,1.1250,9.0865,2.5850/
      DATA(UV_TABLE(667,I),I=1,4)/6.5000,1.2500,11.1793,2.8603/
      DATA(UV_TABLE(668,I),I=1,4)/6.5000,1.3750,13.2721,3.1357/
      DATA(UV_TABLE(669,I),I=1,4)/6.5000,1.5000,15.3649,3.4110/
      DATA(UV_TABLE(670,I),I=1,4)/6.5000,1.6250,17.4577,3.6864/
      DATA(UV_TABLE(671,I),I=1,4)/6.5000,1.7500,19.5505,3.9618/
      DATA(UV_TABLE(672,I),I=1,4)/6.5000,1.8750,21.6434,4.2371/
      DATA(UV_TABLE(673,I),I=1,4)/7.0000,-1.1250,-19.8441,-8.7915/
      DATA(UV_TABLE(674,I),I=1,4)/7.0000,-1.0000,-18.7168,-7.9204/
      DATA(UV_TABLE(675,I),I=1,4)/7.0000,-0.8750,-17.5894,-7.0494/
      DATA(UV_TABLE(676,I),I=1,4)/7.0000,-0.7500,-16.4621,-6.1783/
      DATA(UV_TABLE(677,I),I=1,4)/7.0000,-0.6250,-15.3347,-5.3072/
      DATA(UV_TABLE(678,I),I=1,4)/7.0000,-0.5000,-14.2074,-4.4362/
      DATA(UV_TABLE(679,I),I=1,4)/7.0000,-0.3750,-13.0403,-3.5736/
      DATA(UV_TABLE(680,I),I=1,4)/7.0000,-0.2500,-11.8502,-2.7160/
      DATA(UV_TABLE(681,I),I=1,4)/7.0000,-0.1250,-10.6600,-1.8584/
      DATA(UV_TABLE(682,I),I=1,4)/7.0000,0.0000,-9.4699,-1.0008/
      DATA(UV_TABLE(683,I),I=1,4)/7.0000,0.1250,-8.1772,-0.1776/
      DATA(UV_TABLE(684,I),I=1,4)/7.0000,0.2500,-6.2735,0.4407/
      DATA(UV_TABLE(685,I),I=1,4)/7.0000,0.3000,-5.5120,0.6880/
      DATA(UV_TABLE(686,I),I=1,4)/7.0000,0.3750,-4.3698,1.0590/
      DATA(UV_TABLE(687,I),I=1,4)/7.0000,0.5000,-2.3864,1.6317/
      DATA(UV_TABLE(688,I),I=1,4)/7.0000,0.6000,-0.5396,1.8934/
      DATA(UV_TABLE(689,I),I=1,4)/7.0000,0.6250,-0.0779,1.9588/
      DATA(UV_TABLE(690,I),I=1,4)/7.0000,0.7500,2.2905,2.2064/
      DATA(UV_TABLE(691,I),I=1,4)/7.0000,0.8000,3.2247,2.3232/
      DATA(UV_TABLE(692,I),I=1,4)/7.0000,0.8750,4.6259,2.4985/
      DATA(UV_TABLE(693,I),I=1,4)/7.0000,1.0000,6.9911,2.7448/
      DATA(UV_TABLE(694,I),I=1,4)/7.0000,1.1250,9.3205,3.0513/
      DATA(UV_TABLE(695,I),I=1,4)/7.0000,1.2500,11.6499,3.3578/
      DATA(UV_TABLE(696,I),I=1,4)/7.0000,1.3750,13.9793,3.6643/
      DATA(UV_TABLE(697,I),I=1,4)/7.0000,1.5000,16.3087,3.9707/
      DATA(UV_TABLE(698,I),I=1,4)/7.0000,1.6250,18.6381,4.2772/
      DATA(UV_TABLE(699,I),I=1,4)/7.0000,1.7500,20.9675,4.5837/
      DATA(UV_TABLE(700,I),I=1,4)/7.0000,1.8750,23.2969,4.8902/
      DATA(UV_TABLE(701,I),I=1,4)/7.5000,-1.1250,-22.4814,-9.5584/
      DATA(UV_TABLE(702,I),I=1,4)/7.5000,-1.0000,-21.2422,-8.6010/
      DATA(UV_TABLE(703,I),I=1,4)/7.5000,-0.8750,-20.0031,-7.6436/
      DATA(UV_TABLE(704,I),I=1,4)/7.5000,-0.7500,-18.7640,-6.6861/
      DATA(UV_TABLE(705,I),I=1,4)/7.5000,-0.6250,-17.5249,-5.7287/
      DATA(UV_TABLE(706,I),I=1,4)/7.5000,-0.5000,-16.2857,-4.7713/
      DATA(UV_TABLE(707,I),I=1,4)/7.5000,-0.3750,-14.9807,-3.8280/
      DATA(UV_TABLE(708,I),I=1,4)/7.5000,-0.2500,-13.6725,-2.8853/
      DATA(UV_TABLE(709,I),I=1,4)/7.5000,-0.1250,-12.3644,-1.9427/
      DATA(UV_TABLE(710,I),I=1,4)/7.5000,0.0000,-11.0562,-1.0001/
      DATA(UV_TABLE(711,I),I=1,4)/7.5000,0.1250,-9.5894,-0.1107/
      DATA(UV_TABLE(712,I),I=1,4)/7.5000,0.2500,-7.4969,0.5689/
      DATA(UV_TABLE(713,I),I=1,4)/7.5000,0.3000,-6.6599,0.8407/
      DATA(UV_TABLE(714,I),I=1,4)/7.5000,0.3750,-5.4044,1.2485/
      DATA(UV_TABLE(715,I),I=1,4)/7.5000,0.5000,-3.2661,1.9026/
      DATA(UV_TABLE(716,I),I=1,4)/7.5000,0.6000,-1.2702,2.2217/
      DATA(UV_TABLE(717,I),I=1,4)/7.5000,0.6250,-0.7712,2.3015/
      DATA(UV_TABLE(718,I),I=1,4)/7.5000,0.7500,1.8238,2.5866/
      DATA(UV_TABLE(719,I),I=1,4)/7.5000,0.8000,2.8528,2.7117/
      DATA(UV_TABLE(720,I),I=1,4)/7.5000,0.8750,4.3962,2.8994/
      DATA(UV_TABLE(721,I),I=1,4)/7.5000,1.0000,6.9968,3.1691/
      DATA(UV_TABLE(722,I),I=1,4)/7.5000,1.1250,9.5572,3.5060/
      DATA(UV_TABLE(723,I),I=1,4)/7.5000,1.2500,12.1175,3.8429/
      DATA(UV_TABLE(724,I),I=1,4)/7.5000,1.3750,14.6779,4.1797/
      DATA(UV_TABLE(725,I),I=1,4)/7.5000,1.5000,17.2382,4.5166/
      DATA(UV_TABLE(726,I),I=1,4)/7.5000,1.6250,19.7986,4.8535/
      DATA(UV_TABLE(727,I),I=1,4)/7.5000,1.7500,22.3589,5.1904/
      DATA(UV_TABLE(728,I),I=1,4)/7.5000,1.8750,24.9193,5.5273/
      DATA(UV_TABLE(729,I),I=1,4)/8.0000,-1.1250,-25.1328,-10.3316/
      DATA(UV_TABLE(730,I),I=1,4)/8.0000,-1.0000,-23.7812,-9.2872/
      DATA(UV_TABLE(731,I),I=1,4)/8.0000,-0.8750,-22.4295,-8.2429/
      DATA(UV_TABLE(732,I),I=1,4)/8.0000,-0.7500,-21.0779,-7.1985/
      DATA(UV_TABLE(733,I),I=1,4)/8.0000,-0.6250,-19.7263,-6.1542/
      DATA(UV_TABLE(734,I),I=1,4)/8.0000,-0.5000,-18.3580,-5.1134/
      DATA(UV_TABLE(735,I),I=1,4)/8.0000,-0.3750,-16.9311,-4.0852/
      DATA(UV_TABLE(736,I),I=1,4)/8.0000,-0.2500,-15.5042,-3.0570/
      DATA(UV_TABLE(737,I),I=1,4)/8.0000,-0.1250,-14.0773,-2.0288/
      DATA(UV_TABLE(738,I),I=1,4)/8.0000,0.0000,-12.6504,-1.0006/
      DATA(UV_TABLE(739,I),I=1,4)/8.0000,0.1250,-11.0154,-0.0421/
      DATA(UV_TABLE(740,I),I=1,4)/8.0000,0.2500,-8.7330,0.6991/
      DATA(UV_TABLE(741,I),I=1,4)/8.0000,0.3000,-7.8200,0.9956/
      DATA(UV_TABLE(742,I),I=1,4)/8.0000,0.3750,-6.4506,1.4404/
      DATA(UV_TABLE(743,I),I=1,4)/8.0000,0.5000,-4.1473,2.1705/
      DATA(UV_TABLE(744,I),I=1,4)/8.0000,0.6000,-2.0084,2.5524/
      DATA(UV_TABLE(745,I),I=1,4)/8.0000,0.6250,-1.4737,2.6479/
      DATA(UV_TABLE(746,I),I=1,4)/8.0000,0.7500,1.3517,2.9675/
      DATA(UV_TABLE(747,I),I=1,4)/8.0000,0.8000,2.4774,3.0993/
      DATA(UV_TABLE(748,I),I=1,4)/8.0000,0.8750,4.1659,3.2969/
      DATA(UV_TABLE(749,I),I=1,4)/8.0000,1.0000,6.9982,3.5984/
      DATA(UV_TABLE(750,I),I=1,4)/8.0000,1.1250,9.7910,3.9659/
      DATA(UV_TABLE(751,I),I=1,4)/8.0000,1.2500,12.5837,4.3334/
      DATA(UV_TABLE(752,I),I=1,4)/8.0000,1.3750,15.3765,4.7008/
      DATA(UV_TABLE(753,I),I=1,4)/8.0000,1.5000,18.1693,5.0683/
      DATA(UV_TABLE(754,I),I=1,4)/8.0000,1.6250,20.9621,5.4358/
      DATA(UV_TABLE(755,I),I=1,4)/8.0000,1.7500,23.7549,5.8032/
      DATA(UV_TABLE(756,I),I=1,4)/8.0000,1.8750,26.5477,6.1707/
      END





c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c BLOCK DATA UV_ESTIMATE
c
c     inversion of the coordinate system for the uv-grid
c
c
      BLOCK DATA UV_ESTIMATE
      
      IMPLICIT NONE
          
c     Common-Block 'UV_ESTIMAT' : Estimations of the values of
c                                 (Alpha,Beta) from (G0,Tau)
c       UV_TAU: Gridpoints of UV_A_EST/UV_B_EST in direction of Tau
c       UV_G0: Gridpoints of UV_A_EST/UV_B_EST in direction of G0
c       UV_A_EST: Estimated values of Alpha
c       UV_B_EST: Estimated values of Beta
      REAL UV_TAU(60),UV_G0(60),UV_A_EST(60,60),UV_B_EST(60,60)
      COMMON/UV_ESTIMAT/ UV_TAU,UV_G0,UV_A_EST,UV_B_EST
      
      INTEGER I
            
      DATA (UV_TAU(I),I=1,60) /     -1.000E+00,-9.580E-01,-9.160E-01,
     &                              -8.740E-01,-8.321E-01,-7.901E-01,
     &                              -7.481E-01,-7.061E-01,-6.641E-01,
     &                              -6.221E-01,-5.801E-01,-5.382E-01,
     &                              -4.962E-01,-4.542E-01,-4.122E-01,
     &                              -3.702E-01,-3.282E-01,-2.863E-01,
     &                              -2.443E-01,-2.023E-01,-1.603E-01,
     &                              -1.183E-01,-7.633E-02,-3.434E-02,
     &                               7.643E-03, 4.963E-02, 9.161E-02,
     &                               1.336E-01, 1.756E-01, 2.176E-01,
     &                               2.596E-01, 3.015E-01, 3.435E-01,
     &                               3.855E-01, 4.275E-01, 4.695E-01,
     &                               5.115E-01, 5.534E-01, 5.954E-01,
     &                               6.374E-01, 6.794E-01, 7.214E-01,
     &                               7.634E-01, 8.054E-01, 8.473E-01,
     &                               8.893E-01, 9.313E-01, 9.733E-01,
     &                               1.015E+00, 1.057E+00, 1.099E+00,
     &                               1.141E+00, 1.183E+00, 1.225E+00,
     &                               1.267E+00, 1.309E+00, 1.351E+00,
     &                               1.393E+00, 1.435E+00, 1.477E+00/
      DATA (UV_G0(I),I=1,60)    /   -3.000E+00,-2.831E+00,-2.661E+00,
     &                              -2.492E+00,-2.322E+00,-2.153E+00,
     &                              -1.983E+00,-1.814E+00,-1.644E+00,
     &                              -1.475E+00,-1.305E+00,-1.136E+00,
     &                              -9.661E-01,-7.966E-01,-6.271E-01,
     &                              -4.576E-01,-2.881E-01,-1.186E-01,
     &                               5.085E-02, 2.203E-01, 3.898E-01,
     &                               5.593E-01, 7.288E-01, 8.983E-01,
     &                               1.068E+00, 1.237E+00, 1.407E+00,
     &                               1.576E+00, 1.746E+00, 1.915E+00,
     &                               2.085E+00, 2.254E+00, 2.424E+00,
     &                               2.593E+00, 2.763E+00, 2.932E+00,
     &                               3.102E+00, 3.271E+00, 3.441E+00,
     &                               3.610E+00, 3.780E+00, 3.949E+00,
     &                               4.119E+00, 4.288E+00, 4.458E+00,
     &                               4.627E+00, 4.797E+00, 4.966E+00,
     &                               5.136E+00, 5.305E+00, 5.475E+00,
     &                               5.644E+00, 5.814E+00, 5.983E+00,
     &                               6.153E+00, 6.322E+00, 6.492E+00,
     &                               6.661E+00, 6.831E+00, 7.000E+00/
      DATA (UV_A_EST(1,I),I=1,60) /  5.016E+00, 5.030E+00, 5.043E+00,
     &                               5.058E+00, 5.076E+00, 5.094E+00,
     &                               5.112E+00, 5.130E+00, 5.148E+00,
     &                               5.166E+00, 5.184E+00, 5.202E+00,
     &                               5.220E+00, 5.238E+00, 5.269E+00,
     &                               5.302E+00, 5.335E+00, 5.369E+00,
     &                               5.404E+00, 5.439E+00, 5.474E+00,
     &                               5.510E+00, 5.546E+00, 5.582E+00,
     &                               5.619E+00, 5.656E+00, 5.694E+00,
     &                               5.733E+00, 5.773E+00, 5.812E+00,
     &                               5.851E+00, 5.890E+00, 5.929E+00,
     &                               5.967E+00, 6.005E+00, 6.044E+00,
     &                               6.083E+00, 6.122E+00, 6.161E+00,
     &                               6.200E+00, 6.239E+00, 6.278E+00,
     &                               6.317E+00, 6.356E+00, 6.401E+00,
     &                               6.447E+00, 6.491E+00, 6.535E+00,
     &                               6.578E+00, 6.621E+00, 6.663E+00,
     &                               6.704E+00, 6.745E+00, 6.786E+00,
     &                               6.826E+00, 6.865E+00, 6.905E+00,
     &                               6.944E+00, 6.983E+00, 7.022E+00/
      DATA (UV_A_EST(2,I),I=1,60) /  4.978E+00, 4.991E+00, 5.004E+00,
     &                               5.018E+00, 5.031E+00, 5.044E+00,
     &                               5.060E+00, 5.078E+00, 5.096E+00,
     &                               5.114E+00, 5.132E+00, 5.150E+00,
     &                               5.168E+00, 5.187E+00, 5.219E+00,
     &                               5.252E+00, 5.286E+00, 5.319E+00,
     &                               5.353E+00, 5.388E+00, 5.423E+00,
     &                               5.459E+00, 5.495E+00, 5.531E+00,
     &                               5.568E+00, 5.605E+00, 5.643E+00,
     &                               5.683E+00, 5.722E+00, 5.762E+00,
     &                               5.801E+00, 5.840E+00, 5.878E+00,
     &                               5.917E+00, 5.956E+00, 5.994E+00,
     &                               6.032E+00, 6.071E+00, 6.110E+00,
     &                               6.149E+00, 6.188E+00, 6.226E+00,
     &                               6.266E+00, 6.310E+00, 6.357E+00,
     &                               6.402E+00, 6.447E+00, 6.491E+00,
     &                               6.534E+00, 6.576E+00, 6.618E+00,
     &                               6.660E+00, 6.701E+00, 6.741E+00,
     &                               6.782E+00, 6.822E+00, 6.861E+00,
     &                               6.900E+00, 6.938E+00, 6.976E+00/
      DATA (UV_A_EST(3,I),I=1,60) /  4.939E+00, 4.953E+00, 4.966E+00,
     &                               4.979E+00, 4.992E+00, 5.006E+00,
     &                               5.019E+00, 5.032E+00, 5.045E+00,
     &                               5.061E+00, 5.079E+00, 5.098E+00,
     &                               5.115E+00, 5.138E+00, 5.170E+00,
     &                               5.203E+00, 5.236E+00, 5.269E+00,
     &                               5.303E+00, 5.337E+00, 5.373E+00,
     &                               5.408E+00, 5.443E+00, 5.479E+00,
     &                               5.516E+00, 5.552E+00, 5.592E+00,
     &                               5.633E+00, 5.672E+00, 5.711E+00,
     &                               5.751E+00, 5.790E+00, 5.828E+00,
     &                               5.867E+00, 5.905E+00, 5.944E+00,
     &                               5.981E+00, 6.020E+00, 6.059E+00,
     &                               6.098E+00, 6.137E+00, 6.175E+00,
     &                               6.218E+00, 6.266E+00, 6.312E+00,
     &                               6.358E+00, 6.403E+00, 6.446E+00,
     &                               6.490E+00, 6.532E+00, 6.575E+00,
     &                               6.616E+00, 6.657E+00, 6.698E+00,
     &                               6.738E+00, 6.778E+00, 6.817E+00,
     &                               6.856E+00, 6.894E+00, 6.933E+00/
      DATA (UV_A_EST(4,I),I=1,60) /  4.900E+00, 4.914E+00, 4.928E+00,
     &                               4.941E+00, 4.954E+00, 4.968E+00,
     &                               4.981E+00, 4.994E+00, 5.007E+00,
     &                               5.020E+00, 5.033E+00, 5.046E+00,
     &                               5.063E+00, 5.089E+00, 5.121E+00,
     &                               5.153E+00, 5.186E+00, 5.219E+00,
     &                               5.253E+00, 5.287E+00, 5.321E+00,
     &                               5.356E+00, 5.392E+00, 5.428E+00,
     &                               5.464E+00, 5.501E+00, 5.542E+00,
     &                               5.582E+00, 5.622E+00, 5.662E+00,
     &                               5.701E+00, 5.740E+00, 5.779E+00,
     &                               5.817E+00, 5.856E+00, 5.894E+00,
     &                               5.932E+00, 5.970E+00, 6.008E+00,
     &                               6.047E+00, 6.085E+00, 6.126E+00,
     &                               6.175E+00, 6.222E+00, 6.269E+00,
     &                               6.315E+00, 6.359E+00, 6.403E+00,
     &                               6.446E+00, 6.489E+00, 6.531E+00,
     &                               6.573E+00, 6.614E+00, 6.655E+00,
     &                               6.695E+00, 6.734E+00, 6.774E+00,
     &                               6.813E+00, 6.851E+00, 6.889E+00/
      DATA (UV_A_EST(5,I),I=1,60) /  4.862E+00, 4.875E+00, 4.889E+00,
     &                               4.902E+00, 4.916E+00, 4.929E+00,
     &                               4.943E+00, 4.956E+00, 4.969E+00,
     &                               4.982E+00, 4.996E+00, 5.009E+00,
     &                               5.022E+00, 5.042E+00, 5.071E+00,
     &                               5.103E+00, 5.135E+00, 5.168E+00,
     &                               5.202E+00, 5.236E+00, 5.270E+00,
     &                               5.305E+00, 5.340E+00, 5.376E+00,
     &                               5.412E+00, 5.451E+00, 5.492E+00,
     &                               5.532E+00, 5.573E+00, 5.612E+00,
     &                               5.652E+00, 5.690E+00, 5.729E+00,
     &                               5.768E+00, 5.806E+00, 5.845E+00,
     &                               5.883E+00, 5.921E+00, 5.959E+00,
     &                               5.996E+00, 6.034E+00, 6.083E+00,
     &                               6.132E+00, 6.180E+00, 6.226E+00,
     &                               6.272E+00, 6.317E+00, 6.361E+00,
     &                               6.404E+00, 6.446E+00, 6.488E+00,
     &                               6.530E+00, 6.571E+00, 6.612E+00,
     &                               6.652E+00, 6.692E+00, 6.731E+00,
     &                               6.770E+00, 6.808E+00, 6.846E+00/
      DATA (UV_A_EST(6,I),I=1,60) /  4.823E+00, 4.837E+00, 4.851E+00,
     &                               4.864E+00, 4.878E+00, 4.891E+00,
     &                               4.904E+00, 4.918E+00, 4.931E+00,
     &                               4.944E+00, 4.958E+00, 4.971E+00,
     &                               4.984E+00, 5.007E+00, 5.029E+00,
     &                               5.053E+00, 5.085E+00, 5.118E+00,
     &                               5.151E+00, 5.185E+00, 5.219E+00,
     &                               5.253E+00, 5.288E+00, 5.323E+00,
     &                               5.359E+00, 5.401E+00, 5.442E+00,
     &                               5.483E+00, 5.523E+00, 5.563E+00,
     &                               5.601E+00, 5.639E+00, 5.679E+00,
     &                               5.719E+00, 5.757E+00, 5.796E+00,
     &                               5.834E+00, 5.872E+00, 5.910E+00,
     &                               5.947E+00, 5.992E+00, 6.041E+00,
     &                               6.090E+00, 6.138E+00, 6.185E+00,
     &                               6.230E+00, 6.275E+00, 6.319E+00,
     &                               6.362E+00, 6.405E+00, 6.446E+00,
     &                               6.487E+00, 6.529E+00, 6.569E+00,
     &                               6.610E+00, 6.649E+00, 6.688E+00,
     &                               6.727E+00, 6.766E+00, 6.804E+00/
      DATA (UV_A_EST(7,I),I=1,60) /  4.786E+00, 4.799E+00, 4.812E+00,
     &                               4.826E+00, 4.839E+00, 4.853E+00,
     &                               4.866E+00, 4.880E+00, 4.893E+00,
     &                               4.906E+00, 4.920E+00, 4.933E+00,
     &                               4.949E+00, 4.971E+00, 4.993E+00,
     &                               5.016E+00, 5.039E+00, 5.067E+00,
     &                               5.100E+00, 5.133E+00, 5.167E+00,
     &                               5.201E+00, 5.236E+00, 5.271E+00,
     &                               5.310E+00, 5.352E+00, 5.394E+00,
     &                               5.434E+00, 5.474E+00, 5.513E+00,
     &                               5.549E+00, 5.586E+00, 5.626E+00,
     &                               5.666E+00, 5.708E+00, 5.746E+00,
     &                               5.785E+00, 5.823E+00, 5.861E+00,
     &                               5.904E+00, 5.953E+00, 6.000E+00,
     &                               6.049E+00, 6.097E+00, 6.143E+00,
     &                               6.189E+00, 6.234E+00, 6.278E+00,
     &                               6.321E+00, 6.363E+00, 6.405E+00,
     &                               6.446E+00, 6.487E+00, 6.527E+00,
     &                               6.567E+00, 6.607E+00, 6.647E+00,
     &                               6.685E+00, 6.724E+00, 6.762E+00/
      DATA (UV_A_EST(8,I),I=1,60) /  4.750E+00, 4.763E+00, 4.776E+00,
     &                               4.788E+00, 4.801E+00, 4.815E+00,
     &                               4.828E+00, 4.842E+00, 4.855E+00,
     &                               4.869E+00, 4.882E+00, 4.895E+00,
     &                               4.913E+00, 4.935E+00, 4.958E+00,
     &                               4.980E+00, 5.003E+00, 5.026E+00,
     &                               5.049E+00, 5.081E+00, 5.115E+00,
     &                               5.149E+00, 5.183E+00, 5.218E+00,
     &                               5.261E+00, 5.304E+00, 5.346E+00,
     &                               5.387E+00, 5.426E+00, 5.461E+00,
     &                               5.496E+00, 5.533E+00, 5.571E+00,
     &                               5.611E+00, 5.653E+00, 5.697E+00,
     &                               5.736E+00, 5.774E+00, 5.814E+00,
     &                               5.865E+00, 5.914E+00, 5.962E+00,
     &                               6.009E+00, 6.057E+00, 6.103E+00,
     &                               6.149E+00, 6.194E+00, 6.237E+00,
     &                               6.281E+00, 6.323E+00, 6.365E+00,
     &                               6.406E+00, 6.447E+00, 6.486E+00,
     &                               6.526E+00, 6.566E+00, 6.605E+00,
     &                               6.644E+00, 6.682E+00, 6.728E+00/
      DATA (UV_A_EST(9,I),I=1,60) /  4.713E+00, 4.726E+00, 4.739E+00,
     &                               4.752E+00, 4.765E+00, 4.778E+00,
     &                               4.791E+00, 4.804E+00, 4.817E+00,
     &                               4.831E+00, 4.844E+00, 4.857E+00,
     &                               4.877E+00, 4.899E+00, 4.921E+00,
     &                               4.944E+00, 4.967E+00, 4.989E+00,
     &                               5.012E+00, 5.036E+00, 5.063E+00,
     &                               5.096E+00, 5.130E+00, 5.168E+00,
     &                               5.213E+00, 5.256E+00, 5.299E+00,
     &                               5.339E+00, 5.377E+00, 5.409E+00,
     &                               5.442E+00, 5.477E+00, 5.515E+00,
     &                               5.555E+00, 5.596E+00, 5.640E+00,
     &                               5.685E+00, 5.725E+00, 5.777E+00,
     &                               5.827E+00, 5.876E+00, 5.924E+00,
     &                               5.971E+00, 6.017E+00, 6.064E+00,
     &                               6.109E+00, 6.154E+00, 6.198E+00,
     &                               6.241E+00, 6.284E+00, 6.325E+00,
     &                               6.366E+00, 6.407E+00, 6.447E+00,
     &                               6.486E+00, 6.525E+00, 6.564E+00,
     &                               6.603E+00, 6.652E+00, 6.701E+00/
      DATA (UV_A_EST(10,I),I=1,60) / 4.677E+00, 4.690E+00, 4.704E+00,
     &                               4.717E+00, 4.729E+00, 4.742E+00,
     &                               4.755E+00, 4.768E+00, 4.781E+00,
     &                               4.793E+00, 4.806E+00, 4.820E+00,
     &                               4.841E+00, 4.863E+00, 4.885E+00,
     &                               4.908E+00, 4.930E+00, 4.953E+00,
     &                               4.976E+00, 4.999E+00, 5.022E+00,
     &                               5.045E+00, 5.077E+00, 5.120E+00,
     &                               5.165E+00, 5.210E+00, 5.252E+00,
     &                               5.293E+00, 5.325E+00, 5.357E+00,
     &                               5.390E+00, 5.423E+00, 5.457E+00,
     &                               5.497E+00, 5.537E+00, 5.580E+00,
     &                               5.624E+00, 5.685E+00, 5.740E+00,
     &                               5.790E+00, 5.839E+00, 5.887E+00,
     &                               5.933E+00, 5.979E+00, 6.025E+00,
     &                               6.071E+00, 6.115E+00, 6.159E+00,
     &                               6.202E+00, 6.245E+00, 6.286E+00,
     &                               6.327E+00, 6.367E+00, 6.407E+00,
     &                               6.447E+00, 6.485E+00, 6.527E+00,
     &                               6.576E+00, 6.625E+00, 6.674E+00/
      DATA (UV_A_EST(11,I),I=1,60) / 4.641E+00, 4.654E+00, 4.668E+00,
     &                               4.681E+00, 4.694E+00, 4.707E+00,
     &                               4.720E+00, 4.732E+00, 4.745E+00,
     &                               4.758E+00, 4.771E+00, 4.784E+00,
     &                               4.805E+00, 4.827E+00, 4.849E+00,
     &                               4.871E+00, 4.894E+00, 4.917E+00,
     &                               4.939E+00, 4.962E+00, 4.985E+00,
     &                               5.009E+00, 5.033E+00, 5.073E+00,
     &                               5.119E+00, 5.163E+00, 5.206E+00,
     &                               5.244E+00, 5.277E+00, 5.307E+00,
     &                               5.337E+00, 5.370E+00, 5.404E+00,
     &                               5.438E+00, 5.476E+00, 5.517E+00,
     &                               5.572E+00, 5.638E+00, 5.704E+00,
     &                               5.754E+00, 5.803E+00, 5.850E+00,
     &                               5.897E+00, 5.943E+00, 5.988E+00,
     &                               6.033E+00, 6.077E+00, 6.121E+00,
     &                               6.164E+00, 6.206E+00, 6.248E+00,
     &                               6.289E+00, 6.329E+00, 6.369E+00,
     &                               6.408E+00, 6.452E+00, 6.501E+00,
     &                               6.550E+00, 6.599E+00, 6.648E+00/
      DATA (UV_A_EST(12,I),I=1,60) / 4.605E+00, 4.618E+00, 4.631E+00,
     &                               4.645E+00, 4.658E+00, 4.671E+00,
     &                               4.684E+00, 4.697E+00, 4.710E+00,
     &                               4.723E+00, 4.735E+00, 4.751E+00,
     &                               4.771E+00, 4.792E+00, 4.813E+00,
     &                               4.835E+00, 4.857E+00, 4.880E+00,
     &                               4.903E+00, 4.926E+00, 4.949E+00,
     &                               4.972E+00, 5.000E+00, 5.033E+00,
     &                               5.073E+00, 5.118E+00, 5.161E+00,
     &                               5.198E+00, 5.232E+00, 5.262E+00,
     &                               5.291E+00, 5.317E+00, 5.349E+00,
     &                               5.383E+00, 5.418E+00, 5.460E+00,
     &                               5.523E+00, 5.589E+00, 5.657E+00,
     &                               5.719E+00, 5.767E+00, 5.815E+00,
     &                               5.861E+00, 5.907E+00, 5.952E+00,
     &                               5.995E+00, 6.040E+00, 6.084E+00,
     &                               6.127E+00, 6.169E+00, 6.210E+00,
     &                               6.251E+00, 6.291E+00, 6.331E+00,
     &                               6.379E+00, 6.427E+00, 6.475E+00,
     &                               6.523E+00, 6.573E+00, 6.621E+00/
      DATA (UV_A_EST(13,I),I=1,60) / 4.568E+00, 4.582E+00, 4.595E+00,
     &                               4.609E+00, 4.622E+00, 4.635E+00,
     &                               4.648E+00, 4.661E+00, 4.674E+00,
     &                               4.687E+00, 4.700E+00, 4.717E+00,
     &                               4.737E+00, 4.758E+00, 4.778E+00,
     &                               4.799E+00, 4.821E+00, 4.843E+00,
     &                               4.866E+00, 4.889E+00, 4.912E+00,
     &                               4.935E+00, 4.967E+00, 5.000E+00,
     &                               5.034E+00, 5.073E+00, 5.115E+00,
     &                               5.154E+00, 5.189E+00, 5.221E+00,
     &                               5.250E+00, 5.277E+00, 5.302E+00,
     &                               5.327E+00, 5.366E+00, 5.418E+00,
     &                               5.473E+00, 5.539E+00, 5.607E+00,
     &                               5.677E+00, 5.733E+00, 5.780E+00,
     &                               5.826E+00, 5.872E+00, 5.916E+00,
     &                               5.960E+00, 6.003E+00, 6.047E+00,
     &                               6.089E+00, 6.131E+00, 6.173E+00,
     &                               6.213E+00, 6.256E+00, 6.305E+00,
     &                               6.353E+00, 6.402E+00, 6.450E+00,
     &                               6.498E+00, 6.546E+00, 6.595E+00/
      DATA (UV_A_EST(14,I),I=1,60) / 4.532E+00, 4.546E+00, 4.559E+00,
     &                               4.573E+00, 4.586E+00, 4.599E+00,
     &                               4.613E+00, 4.626E+00, 4.639E+00,
     &                               4.652E+00, 4.665E+00, 4.683E+00,
     &                               4.704E+00, 4.724E+00, 4.744E+00,
     &                               4.765E+00, 4.785E+00, 4.807E+00,
     &                               4.829E+00, 4.852E+00, 4.875E+00,
     &                               4.900E+00, 4.933E+00, 4.967E+00,
     &                               5.000E+00, 5.035E+00, 5.073E+00,
     &                               5.111E+00, 5.147E+00, 5.181E+00,
     &                               5.212E+00, 5.239E+00, 5.265E+00,
     &                               5.290E+00, 5.326E+00, 5.378E+00,
     &                               5.429E+00, 5.488E+00, 5.556E+00,
     &                               5.626E+00, 5.698E+00, 5.746E+00,
     &                               5.792E+00, 5.838E+00, 5.882E+00,
     &                               5.926E+00, 5.969E+00, 6.011E+00,
     &                               6.053E+00, 6.095E+00, 6.136E+00,
     &                               6.182E+00, 6.231E+00, 6.279E+00,
     &                               6.328E+00, 6.376E+00, 6.424E+00,
     &                               6.472E+00, 6.521E+00, 6.570E+00/
      DATA (UV_A_EST(15,I),I=1,60) / 4.493E+00, 4.510E+00, 4.523E+00,
     &                               4.537E+00, 4.550E+00, 4.564E+00,
     &                               4.577E+00, 4.590E+00, 4.603E+00,
     &                               4.617E+00, 4.630E+00, 4.650E+00,
     &                               4.670E+00, 4.690E+00, 4.710E+00,
     &                               4.731E+00, 4.751E+00, 4.772E+00,
     &                               4.793E+00, 4.815E+00, 4.837E+00,
     &                               4.867E+00, 4.900E+00, 4.933E+00,
     &                               4.967E+00, 5.001E+00, 5.035E+00,
     &                               5.072E+00, 5.108E+00, 5.142E+00,
     &                               5.174E+00, 5.203E+00, 5.229E+00,
     &                               5.263E+00, 5.298E+00, 5.339E+00,
     &                               5.390E+00, 5.439E+00, 5.503E+00,
     &                               5.574E+00, 5.647E+00, 5.713E+00,
     &                               5.759E+00, 5.804E+00, 5.849E+00,
     &                               5.892E+00, 5.935E+00, 5.977E+00,
     &                               6.018E+00, 6.060E+00, 6.108E+00,
     &                               6.157E+00, 6.206E+00, 6.254E+00,
     &                               6.302E+00, 6.351E+00, 6.399E+00,
     &                               6.447E+00, 6.495E+00, 6.544E+00/
      DATA (UV_A_EST(16,I),I=1,60) / 4.424E+00, 4.450E+00, 4.476E+00,
     &                               4.501E+00, 4.515E+00, 4.528E+00,
     &                               4.542E+00, 4.555E+00, 4.568E+00,
     &                               4.581E+00, 4.596E+00, 4.616E+00,
     &                               4.636E+00, 4.656E+00, 4.676E+00,
     &                               4.697E+00, 4.717E+00, 4.738E+00,
     &                               4.759E+00, 4.779E+00, 4.801E+00,
     &                               4.833E+00, 4.866E+00, 4.899E+00,
     &                               4.933E+00, 4.967E+00, 5.001E+00,
     &                               5.036E+00, 5.071E+00, 5.105E+00,
     &                               5.137E+00, 5.168E+00, 5.203E+00,
     &                               5.237E+00, 5.271E+00, 5.306E+00,
     &                               5.351E+00, 5.401E+00, 5.450E+00,
     &                               5.520E+00, 5.593E+00, 5.668E+00,
     &                               5.727E+00, 5.772E+00, 5.816E+00,
     &                               5.859E+00, 5.901E+00, 5.943E+00,
     &                               5.986E+00, 6.035E+00, 6.083E+00,
     &                               6.132E+00, 6.181E+00, 6.229E+00,
     &                               6.277E+00, 6.326E+00, 6.374E+00,
     &                               6.422E+00, 6.470E+00, 6.518E+00/
      DATA (UV_A_EST(17,I),I=1,60) / 4.356E+00, 4.382E+00, 4.408E+00,
     &                               4.434E+00, 4.460E+00, 4.486E+00,
     &                               4.506E+00, 4.520E+00, 4.533E+00,
     &                               4.546E+00, 4.562E+00, 4.582E+00,
     &                               4.602E+00, 4.622E+00, 4.642E+00,
     &                               4.662E+00, 4.683E+00, 4.703E+00,
     &                               4.724E+00, 4.745E+00, 4.770E+00,
     &                               4.800E+00, 4.832E+00, 4.865E+00,
     &                               4.899E+00, 4.933E+00, 4.967E+00,
     &                               5.001E+00, 5.036E+00, 5.070E+00,
     &                               5.103E+00, 5.140E+00, 5.176E+00,
     &                               5.210E+00, 5.244E+00, 5.280E+00,
     &                               5.315E+00, 5.363E+00, 5.412E+00,
     &                               5.465E+00, 5.538E+00, 5.613E+00,
     &                               5.691E+00, 5.740E+00, 5.783E+00,
     &                               5.826E+00, 5.869E+00, 5.915E+00,
     &                               5.962E+00, 6.010E+00, 6.059E+00,
     &                               6.107E+00, 6.156E+00, 6.204E+00,
     &                               6.252E+00, 6.300E+00, 6.349E+00,
     &                               6.396E+00, 6.444E+00, 6.493E+00/
      DATA (UV_A_EST(18,I),I=1,60) / 4.288E+00, 4.314E+00, 4.340E+00,
     &                               4.366E+00, 4.392E+00, 4.418E+00,
     &                               4.444E+00, 4.470E+00, 4.495E+00,
     &                               4.511E+00, 4.528E+00, 4.548E+00,
     &                               4.568E+00, 4.588E+00, 4.608E+00,
     &                               4.628E+00, 4.648E+00, 4.669E+00,
     &                               4.690E+00, 4.710E+00, 4.738E+00,
     &                               4.768E+00, 4.799E+00, 4.831E+00,
     &                               4.865E+00, 4.898E+00, 4.933E+00,
     &                               4.967E+00, 5.002E+00, 5.037E+00,
     &                               5.077E+00, 5.115E+00, 5.151E+00,
     &                               5.184E+00, 5.217E+00, 5.252E+00,
     &                               5.288E+00, 5.326E+00, 5.375E+00,
     &                               5.423E+00, 5.481E+00, 5.557E+00,
     &                               5.635E+00, 5.708E+00, 5.751E+00,
     &                               5.794E+00, 5.842E+00, 5.890E+00,
     &                               5.938E+00, 5.986E+00, 6.034E+00,
     &                               6.082E+00, 6.131E+00, 6.179E+00,
     &                               6.228E+00, 6.276E+00, 6.324E+00,
     &                               6.372E+00, 6.420E+00, 6.467E+00/
      DATA (UV_A_EST(19,I),I=1,60) / 4.219E+00, 4.245E+00, 4.271E+00,
     &                               4.298E+00, 4.324E+00, 4.350E+00,
     &                               4.376E+00, 4.402E+00, 4.428E+00,
     &                               4.454E+00, 4.490E+00, 4.514E+00,
     &                               4.534E+00, 4.554E+00, 4.574E+00,
     &                               4.594E+00, 4.614E+00, 4.634E+00,
     &                               4.655E+00, 4.678E+00, 4.707E+00,
     &                               4.736E+00, 4.766E+00, 4.797E+00,
     &                               4.830E+00, 4.864E+00, 4.898E+00,
     &                               4.932E+00, 4.967E+00, 5.005E+00,
     &                               5.055E+00, 5.093E+00, 5.128E+00,
     &                               5.160E+00, 5.191E+00, 5.224E+00,
     &                               5.260E+00, 5.297E+00, 5.339E+00,
     &                               5.386E+00, 5.433E+00, 5.498E+00,
     &                               5.577E+00, 5.658E+00, 5.722E+00,
     &                               5.770E+00, 5.818E+00, 5.866E+00,
     &                               5.914E+00, 5.962E+00, 6.009E+00,
     &                               6.058E+00, 6.106E+00, 6.155E+00,
     &                               6.203E+00, 6.251E+00, 6.299E+00,
     &                               6.347E+00, 6.395E+00, 6.442E+00/
      DATA (UV_A_EST(20,I),I=1,60) / 4.150E+00, 4.177E+00, 4.203E+00,
     &                               4.229E+00, 4.255E+00, 4.281E+00,
     &                               4.307E+00, 4.333E+00, 4.359E+00,
     &                               4.386E+00, 4.423E+00, 4.461E+00,
     &                               4.500E+00, 4.519E+00, 4.539E+00,
     &                               4.559E+00, 4.580E+00, 4.600E+00,
     &                               4.620E+00, 4.646E+00, 4.675E+00,
     &                               4.704E+00, 4.734E+00, 4.765E+00,
     &                               4.796E+00, 4.829E+00, 4.863E+00,
     &                               4.897E+00, 4.932E+00, 4.979E+00,
     &                               5.030E+00, 5.072E+00, 5.106E+00,
     &                               5.139E+00, 5.169E+00, 5.197E+00,
     &                               5.232E+00, 5.269E+00, 5.306E+00,
     &                               5.350E+00, 5.397E+00, 5.442E+00,
     &                               5.515E+00, 5.605E+00, 5.698E+00,
     &                               5.747E+00, 5.795E+00, 5.842E+00,
     &                               5.890E+00, 5.938E+00, 5.985E+00,
     &                               6.033E+00, 6.082E+00, 6.130E+00,
     &                               6.178E+00, 6.226E+00, 6.274E+00,
     &                               6.322E+00, 6.370E+00, 6.417E+00/
      DATA (UV_A_EST(21,I),I=1,60) / 4.082E+00, 4.108E+00, 4.134E+00,
     &                               4.161E+00, 4.187E+00, 4.213E+00,
     &                               4.239E+00, 4.265E+00, 4.291E+00,
     &                               4.320E+00, 4.357E+00, 4.394E+00,
     &                               4.433E+00, 4.471E+00, 4.505E+00,
     &                               4.525E+00, 4.545E+00, 4.565E+00,
     &                               4.586E+00, 4.613E+00, 4.642E+00,
     &                               4.671E+00, 4.701E+00, 4.732E+00,
     &                               4.763E+00, 4.795E+00, 4.828E+00,
     &                               4.862E+00, 4.902E+00, 4.952E+00,
     &                               5.002E+00, 5.053E+00, 5.087E+00,
     &                               5.118E+00, 5.148E+00, 5.176E+00,
     &                               5.203E+00, 5.240E+00, 5.278E+00,
     &                               5.316E+00, 5.361E+00, 5.408E+00,
     &                               5.467E+00, 5.560E+00, 5.652E+00,
     &                               5.723E+00, 5.771E+00, 5.819E+00,
     &                               5.866E+00, 5.914E+00, 5.961E+00,
     &                               6.009E+00, 6.057E+00, 6.106E+00,
     &                               6.154E+00, 6.201E+00, 6.249E+00,
     &                               6.297E+00, 6.345E+00, 6.392E+00/
      DATA (UV_A_EST(22,I),I=1,60) / 4.014E+00, 4.040E+00, 4.066E+00,
     &                               4.092E+00, 4.119E+00, 4.145E+00,
     &                               4.171E+00, 4.197E+00, 4.223E+00,
     &                               4.253E+00, 4.290E+00, 4.327E+00,
     &                               4.365E+00, 4.403E+00, 4.442E+00,
     &                               4.481E+00, 4.510E+00, 4.530E+00,
     &                               4.553E+00, 4.581E+00, 4.609E+00,
     &                               4.638E+00, 4.668E+00, 4.699E+00,
     &                               4.730E+00, 4.761E+00, 4.793E+00,
     &                               4.827E+00, 4.875E+00, 4.924E+00,
     &                               4.975E+00, 5.027E+00, 5.068E+00,
     &                               5.099E+00, 5.129E+00, 5.156E+00,
     &                               5.182E+00, 5.210E+00, 5.248E+00,
     &                               5.287E+00, 5.330E+00, 5.382E+00,
     &                               5.434E+00, 5.515E+00, 5.607E+00,
     &                               5.699E+00, 5.747E+00, 5.795E+00,
     &                               5.842E+00, 5.890E+00, 5.937E+00,
     &                               5.985E+00, 6.032E+00, 6.081E+00,
     &                               6.129E+00, 6.177E+00, 6.225E+00,
     &                               6.272E+00, 6.320E+00, 6.367E+00/
      DATA (UV_A_EST(23,I),I=1,60) / 3.947E+00, 3.972E+00, 3.998E+00,
     &                               4.024E+00, 4.050E+00, 4.077E+00,
     &                               4.103E+00, 4.129E+00, 4.155E+00,
     &                               4.187E+00, 4.223E+00, 4.260E+00,
     &                               4.297E+00, 4.335E+00, 4.373E+00,
     &                               4.412E+00, 4.452E+00, 4.491E+00,
     &                               4.521E+00, 4.548E+00, 4.576E+00,
     &                               4.605E+00, 4.635E+00, 4.665E+00,
     &                               4.696E+00, 4.727E+00, 4.759E+00,
     &                               4.799E+00, 4.847E+00, 4.896E+00,
     &                               4.947E+00, 4.999E+00, 5.052E+00,
     &                               5.082E+00, 5.111E+00, 5.138E+00,
     &                               5.163E+00, 5.188E+00, 5.218E+00,
     &                               5.262E+00, 5.307E+00, 5.357E+00,
     &                               5.409E+00, 5.470E+00, 5.562E+00,
     &                               5.654E+00, 5.724E+00, 5.772E+00,
     &                               5.819E+00, 5.867E+00, 5.914E+00,
     &                               5.961E+00, 6.008E+00, 6.056E+00,
     &                               6.104E+00, 6.152E+00, 6.200E+00,
     &                               6.248E+00, 6.295E+00, 6.343E+00/
      DATA (UV_A_EST(24,I),I=1,60) / 3.880E+00, 3.904E+00, 3.929E+00,
     &                               3.955E+00, 3.981E+00, 4.008E+00,
     &                               4.035E+00, 4.061E+00, 4.087E+00,
     &                               4.120E+00, 4.156E+00, 4.192E+00,
     &                               4.229E+00, 4.267E+00, 4.305E+00,
     &                               4.343E+00, 4.382E+00, 4.423E+00,
     &                               4.477E+00, 4.515E+00, 4.543E+00,
     &                               4.572E+00, 4.601E+00, 4.631E+00,
     &                               4.661E+00, 4.693E+00, 4.725E+00,
     &                               4.771E+00, 4.818E+00, 4.868E+00,
     &                               4.919E+00, 4.971E+00, 5.025E+00,
     &                               5.066E+00, 5.094E+00, 5.120E+00,
     &                               5.146E+00, 5.171E+00, 5.198E+00,
     &                               5.241E+00, 5.286E+00, 5.332E+00,
     &                               5.384E+00, 5.436E+00, 5.517E+00,
     &                               5.609E+00, 5.701E+00, 5.748E+00,
     &                               5.796E+00, 5.843E+00, 5.890E+00,
     &                               5.938E+00, 5.985E+00, 6.032E+00,
     &                               6.080E+00, 6.127E+00, 6.175E+00,
     &                               6.223E+00, 6.270E+00, 6.317E+00/
      DATA (UV_A_EST(25,I),I=1,60) / 3.813E+00, 3.836E+00, 3.860E+00,
     &                               3.884E+00, 3.910E+00, 3.936E+00,
     &                               3.964E+00, 3.992E+00, 4.019E+00,
     &                               4.053E+00, 4.089E+00, 4.125E+00,
     &                               4.161E+00, 4.198E+00, 4.236E+00,
     &                               4.274E+00, 4.312E+00, 4.357E+00,
     &                               4.410E+00, 4.464E+00, 4.510E+00,
     &                               4.538E+00, 4.567E+00, 4.596E+00,
     &                               4.627E+00, 4.658E+00, 4.696E+00,
     &                               4.742E+00, 4.790E+00, 4.840E+00,
     &                               4.890E+00, 4.943E+00, 4.996E+00,
     &                               5.051E+00, 5.078E+00, 5.104E+00,
     &                               5.131E+00, 5.158E+00, 5.185E+00,
     &                               5.220E+00, 5.265E+00, 5.309E+00,
     &                               5.359E+00, 5.411E+00, 5.473E+00,
     &                               5.565E+00, 5.656E+00, 5.725E+00,
     &                               5.772E+00, 5.819E+00, 5.867E+00,
     &                               5.914E+00, 5.961E+00, 6.008E+00,
     &                               6.055E+00, 6.103E+00, 6.151E+00,
     &                               6.198E+00, 6.247E+00, 6.296E+00/
      DATA (UV_A_EST(26,I),I=1,60) / 3.747E+00, 3.768E+00, 3.790E+00,
     &                               3.813E+00, 3.838E+00, 3.863E+00,
     &                               3.889E+00, 3.916E+00, 3.944E+00,
     &                               3.984E+00, 4.021E+00, 4.057E+00,
     &                               4.093E+00, 4.129E+00, 4.166E+00,
     &                               4.204E+00, 4.243E+00, 4.290E+00,
     &                               4.342E+00, 4.395E+00, 4.451E+00,
     &                               4.504E+00, 4.532E+00, 4.562E+00,
     &                               4.592E+00, 4.623E+00, 4.667E+00,
     &                               4.713E+00, 4.761E+00, 4.810E+00,
     &                               4.861E+00, 4.914E+00, 4.967E+00,
     &                               5.022E+00, 5.063E+00, 5.091E+00,
     &                               5.118E+00, 5.146E+00, 5.173E+00,
     &                               5.200E+00, 5.244E+00, 5.288E+00,
     &                               5.335E+00, 5.386E+00, 5.438E+00,
     &                               5.520E+00, 5.611E+00, 5.701E+00,
     &                               5.748E+00, 5.795E+00, 5.843E+00,
     &                               5.890E+00, 5.937E+00, 5.984E+00,
     &                               6.031E+00, 6.079E+00, 6.128E+00,
     &                               6.177E+00, 6.226E+00, 6.276E+00/
      DATA (UV_A_EST(27,I),I=1,60) / 3.680E+00, 3.700E+00, 3.721E+00,
     &                               3.743E+00, 3.765E+00, 3.789E+00,
     &                               3.814E+00, 3.839E+00, 3.867E+00,
     &                               3.904E+00, 3.944E+00, 3.986E+00,
     &                               4.024E+00, 4.060E+00, 4.097E+00,
     &                               4.134E+00, 4.173E+00, 4.222E+00,
     &                               4.273E+00, 4.326E+00, 4.380E+00,
     &                               4.437E+00, 4.494E+00, 4.526E+00,
     &                               4.556E+00, 4.593E+00, 4.637E+00,
     &                               4.682E+00, 4.730E+00, 4.780E+00,
     &                               4.832E+00, 4.884E+00, 4.938E+00,
     &                               4.993E+00, 5.051E+00, 5.079E+00,
     &                               5.106E+00, 5.133E+00, 5.160E+00,
     &                               5.187E+00, 5.223E+00, 5.268E+00,
     &                               5.311E+00, 5.361E+00, 5.413E+00,
     &                               5.475E+00, 5.567E+00, 5.658E+00,
     &                               5.725E+00, 5.772E+00, 5.819E+00,
     &                               5.867E+00, 5.914E+00, 5.961E+00,
     &                               6.009E+00, 6.058E+00, 6.107E+00,
     &                               6.156E+00, 6.206E+00, 6.255E+00/
      DATA (UV_A_EST(28,I),I=1,60) / 3.613E+00, 3.632E+00, 3.651E+00,
     &                               3.672E+00, 3.693E+00, 3.715E+00,
     &                               3.738E+00, 3.763E+00, 3.789E+00,
     &                               3.824E+00, 3.861E+00, 3.900E+00,
     &                               3.943E+00, 3.988E+00, 4.027E+00,
     &                               4.064E+00, 4.105E+00, 4.154E+00,
     &                               4.203E+00, 4.256E+00, 4.309E+00,
     &                               4.364E+00, 4.422E+00, 4.481E+00,
     &                               4.521E+00, 4.563E+00, 4.606E+00,
     &                               4.652E+00, 4.699E+00, 4.749E+00,
     &                               4.801E+00, 4.854E+00, 4.908E+00,
     &                               4.965E+00, 5.026E+00, 5.067E+00,
     &                               5.094E+00, 5.121E+00, 5.148E+00,
     &                               5.175E+00, 5.203E+00, 5.247E+00,
     &                               5.291E+00, 5.337E+00, 5.388E+00,
     &                               5.439E+00, 5.521E+00, 5.615E+00,
     &                               5.702E+00, 5.749E+00, 5.796E+00,
     &                               5.843E+00, 5.891E+00, 5.940E+00,
     &                               5.989E+00, 6.037E+00, 6.086E+00,
     &                               6.136E+00, 6.185E+00, 6.234E+00/
      DATA (UV_A_EST(29,I),I=1,60) / 3.546E+00, 3.564E+00, 3.582E+00,
     &                               3.601E+00, 3.621E+00, 3.641E+00,
     &                               3.663E+00, 3.685E+00, 3.711E+00,
     &                               3.743E+00, 3.777E+00, 3.814E+00,
     &                               3.854E+00, 3.896E+00, 3.942E+00,
     &                               3.991E+00, 4.037E+00, 4.084E+00,
     &                               4.133E+00, 4.184E+00, 4.237E+00,
     &                               4.291E+00, 4.348E+00, 4.407E+00,
     &                               4.480E+00, 4.532E+00, 4.575E+00,
     &                               4.620E+00, 4.668E+00, 4.718E+00,
     &                               4.770E+00, 4.823E+00, 4.878E+00,
     &                               4.937E+00, 4.999E+00, 5.055E+00,
     &                               5.082E+00, 5.109E+00, 5.136E+00,
     &                               5.163E+00, 5.189E+00, 5.226E+00,
     &                               5.270E+00, 5.313E+00, 5.363E+00,
     &                               5.418E+00, 5.489E+00, 5.579E+00,
     &                               5.664E+00, 5.725E+00, 5.774E+00,
     &                               5.822E+00, 5.871E+00, 5.919E+00,
     &                               5.968E+00, 6.017E+00, 6.066E+00,
     &                               6.115E+00, 6.164E+00, 6.213E+00/
      DATA (UV_A_EST(30,I),I=1,60) / 3.459E+00, 3.492E+00, 3.513E+00,
     &                               3.530E+00, 3.548E+00, 3.568E+00,
     &                               3.588E+00, 3.608E+00, 3.632E+00,
     &                               3.662E+00, 3.694E+00, 3.728E+00,
     &                               3.764E+00, 3.804E+00, 3.846E+00,
     &                               3.891E+00, 3.953E+00, 4.014E+00,
     &                               4.062E+00, 4.112E+00, 4.163E+00,
     &                               4.217E+00, 4.272E+00, 4.330E+00,
     &                               4.411E+00, 4.500E+00, 4.543E+00,
     &                               4.588E+00, 4.636E+00, 4.685E+00,
     &                               4.737E+00, 4.791E+00, 4.849E+00,
     &                               4.908E+00, 4.970E+00, 5.034E+00,
     &                               5.070E+00, 5.098E+00, 5.124E+00,
     &                               5.151E+00, 5.177E+00, 5.206E+00,
     &                               5.250E+00, 5.296E+00, 5.346E+00,
     &                               5.402E+00, 5.461E+00, 5.546E+00,
     &                               5.629E+00, 5.705E+00, 5.754E+00,
     &                               5.802E+00, 5.851E+00, 5.899E+00,
     &                               5.948E+00, 5.996E+00, 6.045E+00,
     &                               6.094E+00, 6.144E+00, 6.193E+00/
      DATA (UV_A_EST(31,I),I=1,60) / 3.323E+00, 3.357E+00, 3.390E+00,
     &                               3.423E+00, 3.456E+00, 3.489E+00,
     &                               3.512E+00, 3.531E+00, 3.553E+00,
     &                               3.580E+00, 3.610E+00, 3.641E+00,
     &                               3.674E+00, 3.711E+00, 3.749E+00,
     &                               3.791E+00, 3.850E+00, 3.915E+00,
     &                               3.985E+00, 4.039E+00, 4.089E+00,
     &                               4.142E+00, 4.196E+00, 4.260E+00,
     &                               4.341E+00, 4.428E+00, 4.510E+00,
     &                               4.555E+00, 4.602E+00, 4.652E+00,
     &                               4.703E+00, 4.760E+00, 4.819E+00,
     &                               4.879E+00, 4.941E+00, 5.006E+00,
     &                               5.059E+00, 5.086E+00, 5.112E+00,
     &                               5.139E+00, 5.165E+00, 5.193E+00,
     &                               5.234E+00, 5.282E+00, 5.331E+00,
     &                               5.386E+00, 5.442E+00, 5.519E+00,
     &                               5.599E+00, 5.677E+00, 5.734E+00,
     &                               5.782E+00, 5.831E+00, 5.879E+00,
     &                               5.927E+00, 5.976E+00, 6.025E+00,
     &                               6.074E+00, 6.123E+00, 6.172E+00/
      DATA (UV_A_EST(32,I),I=1,60) / 3.187E+00, 3.220E+00, 3.254E+00,
     &                               3.288E+00, 3.321E+00, 3.354E+00,
     &                               3.387E+00, 3.420E+00, 3.457E+00,
     &                               3.499E+00, 3.525E+00, 3.554E+00,
     &                               3.584E+00, 3.617E+00, 3.652E+00,
     &                               3.690E+00, 3.746E+00, 3.807E+00,
     &                               3.873E+00, 3.945E+00, 4.014E+00,
     &                               4.065E+00, 4.119E+00, 4.189E+00,
     &                               4.268E+00, 4.353E+00, 4.446E+00,
     &                               4.521E+00, 4.568E+00, 4.617E+00,
     &                               4.670E+00, 4.727E+00, 4.788E+00,
     &                               4.849E+00, 4.912E+00, 4.976E+00,
     &                               5.044E+00, 5.074E+00, 5.101E+00,
     &                               5.128E+00, 5.156E+00, 5.184E+00,
     &                               5.221E+00, 5.268E+00, 5.317E+00,
     &                               5.370E+00, 5.425E+00, 5.493E+00,
     &                               5.571E+00, 5.647E+00, 5.714E+00,
     &                               5.762E+00, 5.810E+00, 5.859E+00,
     &                               5.907E+00, 5.956E+00, 6.004E+00,
     &                               6.053E+00, 6.102E+00, 6.151E+00/
      DATA (UV_A_EST(33,I),I=1,60) / 3.052E+00, 3.085E+00, 3.118E+00,
     &                               3.152E+00, 3.185E+00, 3.218E+00,
     &                               3.252E+00, 3.285E+00, 3.324E+00,
     &                               3.366E+00, 3.408E+00, 3.450E+00,
     &                               3.491E+00, 3.523E+00, 3.555E+00,
     &                               3.589E+00, 3.642E+00, 3.699E+00,
     &                               3.761E+00, 3.829E+00, 3.901E+00,
     &                               3.981E+00, 4.044E+00, 4.115E+00,
     &                               4.193E+00, 4.276E+00, 4.368E+00,
     &                               4.468E+00, 4.533E+00, 4.582E+00,
     &                               4.635E+00, 4.693E+00, 4.755E+00,
     &                               4.818E+00, 4.881E+00, 4.946E+00,
     &                               5.014E+00, 5.062E+00, 5.090E+00,
     &                               5.117E+00, 5.146E+00, 5.175E+00,
     &                               5.207E+00, 5.255E+00, 5.303E+00,
     &                               5.355E+00, 5.409E+00, 5.469E+00,
     &                               5.545E+00, 5.619E+00, 5.690E+00,
     &                               5.742E+00, 5.790E+00, 5.839E+00,
     &                               5.887E+00, 5.935E+00, 5.984E+00,
     &                               6.033E+00, 6.082E+00, 6.131E+00/
      DATA (UV_A_EST(34,I),I=1,60) / 2.915E+00, 2.949E+00, 2.983E+00,
     &                               3.016E+00, 3.049E+00, 3.083E+00,
     &                               3.116E+00, 3.149E+00, 3.190E+00,
     &                               3.233E+00, 3.275E+00, 3.317E+00,
     &                               3.359E+00, 3.401E+00, 3.443E+00,
     &                               3.485E+00, 3.537E+00, 3.591E+00,
     &                               3.649E+00, 3.712E+00, 3.780E+00,
     &                               3.854E+00, 3.947E+00, 4.039E+00,
     &                               4.115E+00, 4.197E+00, 4.287E+00,
     &                               4.385E+00, 4.492E+00, 4.546E+00,
     &                               4.600E+00, 4.658E+00, 4.719E+00,
     &                               4.785E+00, 4.849E+00, 4.915E+00,
     &                               4.983E+00, 5.051E+00, 5.079E+00,
     &                               5.107E+00, 5.136E+00, 5.165E+00,
     &                               5.196E+00, 5.241E+00, 5.289E+00,
     &                               5.340E+00, 5.394E+00, 5.448E+00,
     &                               5.520E+00, 5.592E+00, 5.662E+00,
     &                               5.722E+00, 5.770E+00, 5.819E+00,
     &                               5.867E+00, 5.916E+00, 5.964E+00,
     &                               6.013E+00, 6.062E+00, 6.111E+00/
      DATA (UV_A_EST(35,I),I=1,60) / 2.776E+00, 2.812E+00, 2.846E+00,
     &                               2.881E+00, 2.915E+00, 2.948E+00,
     &                               2.982E+00, 3.015E+00, 3.058E+00,
     &                               3.100E+00, 3.142E+00, 3.184E+00,
     &                               3.226E+00, 3.269E+00, 3.311E+00,
     &                               3.361E+00, 3.421E+00, 3.480E+00,
     &                               3.536E+00, 3.594E+00, 3.657E+00,
     &                               3.726E+00, 3.815E+00, 3.933E+00,
     &                               4.035E+00, 4.115E+00, 4.202E+00,
     &                               4.298E+00, 4.403E+00, 4.510E+00,
     &                               4.564E+00, 4.622E+00, 4.683E+00,
     &                               4.749E+00, 4.816E+00, 4.883E+00,
     &                               4.951E+00, 5.022E+00, 5.067E+00,
     &                               5.096E+00, 5.125E+00, 5.156E+00,
     &                               5.187E+00, 5.228E+00, 5.276E+00,
     &                               5.325E+00, 5.379E+00, 5.432E+00,
     &                               5.498E+00, 5.567E+00, 5.636E+00,
     &                               5.703E+00, 5.751E+00, 5.799E+00,
     &                               5.848E+00, 5.896E+00, 5.944E+00,
     &                               5.992E+00, 6.041E+00, 6.090E+00/
      DATA (UV_A_EST(36,I),I=1,60) / 2.637E+00, 2.673E+00, 2.709E+00,
     &                               2.744E+00, 2.779E+00, 2.814E+00,
     &                               2.848E+00, 2.884E+00, 2.926E+00,
     &                               2.968E+00, 3.009E+00, 3.052E+00,
     &                               3.094E+00, 3.136E+00, 3.178E+00,
     &                               3.237E+00, 3.296E+00, 3.356E+00,
     &                               3.416E+00, 3.475E+00, 3.534E+00,
     &                               3.597E+00, 3.683E+00, 3.796E+00,
     &                               3.918E+00, 4.029E+00, 4.114E+00,
     &                               4.208E+00, 4.313E+00, 4.435E+00,
     &                               4.527E+00, 4.584E+00, 4.646E+00,
     &                               4.711E+00, 4.781E+00, 4.849E+00,
     &                               4.918E+00, 4.989E+00, 5.054E+00,
     &                               5.084E+00, 5.114E+00, 5.145E+00,
     &                               5.177E+00, 5.214E+00, 5.262E+00,
     &                               5.311E+00, 5.364E+00, 5.417E+00,
     &                               5.476E+00, 5.545E+00, 5.611E+00,
     &                               5.677E+00, 5.731E+00, 5.779E+00,
     &                               5.828E+00, 5.876E+00, 5.924E+00,
     &                               5.972E+00, 6.020E+00, 6.069E+00/
      DATA (UV_A_EST(37,I),I=1,60) / 2.499E+00, 2.536E+00, 2.572E+00,
     &                               2.608E+00, 2.643E+00, 2.678E+00,
     &                               2.713E+00, 2.752E+00, 2.793E+00,
     &                               2.835E+00, 2.877E+00, 2.919E+00,
     &                               2.962E+00, 3.004E+00, 3.053E+00,
     &                               3.112E+00, 3.171E+00, 3.230E+00,
     &                               3.290E+00, 3.350E+00, 3.410E+00,
     &                               3.470E+00, 3.550E+00, 3.659E+00,
     &                               3.776E+00, 3.901E+00, 4.022E+00,
     &                               4.113E+00, 4.221E+00, 4.342E+00,
     &                               4.470E+00, 4.545E+00, 4.606E+00,
     &                               4.672E+00, 4.742E+00, 4.815E+00,
     &                               4.884E+00, 4.956E+00, 5.028E+00,
     &                               5.071E+00, 5.102E+00, 5.134E+00,
     &                               5.167E+00, 5.200E+00, 5.249E+00,
     &                               5.297E+00, 5.349E+00, 5.402E+00,
     &                               5.457E+00, 5.523E+00, 5.588E+00,
     &                               5.651E+00, 5.711E+00, 5.759E+00,
     &                               5.807E+00, 5.856E+00, 5.904E+00,
     &                               5.952E+00, 6.000E+00, 6.048E+00/
      DATA (UV_A_EST(38,I),I=1,60) / 2.364E+00, 2.399E+00, 2.436E+00,
     &                               2.472E+00, 2.508E+00, 2.544E+00,
     &                               2.580E+00, 2.619E+00, 2.660E+00,
     &                               2.701E+00, 2.743E+00, 2.785E+00,
     &                               2.827E+00, 2.870E+00, 2.927E+00,
     &                               2.987E+00, 3.046E+00, 3.105E+00,
     &                               3.164E+00, 3.224E+00, 3.284E+00,
     &                               3.349E+00, 3.433E+00, 3.522E+00,
     &                               3.634E+00, 3.754E+00, 3.885E+00,
     &                               4.016E+00, 4.128E+00, 4.247E+00,
     &                               4.374E+00, 4.504E+00, 4.565E+00,
     &                               4.631E+00, 4.701E+00, 4.778E+00,
     &                               4.850E+00, 4.922E+00, 4.995E+00,
     &                               5.058E+00, 5.089E+00, 5.122E+00,
     &                               5.156E+00, 5.190E+00, 5.235E+00,
     &                               5.284E+00, 5.335E+00, 5.387E+00,
     &                               5.440E+00, 5.502E+00, 5.565E+00,
     &                               5.628E+00, 5.689E+00, 5.739E+00,
     &                               5.787E+00, 5.836E+00, 5.884E+00,
     &                               5.932E+00, 5.980E+00, 6.027E+00/
      DATA (UV_A_EST(39,I),I=1,60) / 2.227E+00, 2.263E+00, 2.298E+00,
     &                               2.334E+00, 2.370E+00, 2.407E+00,
     &                               2.444E+00, 2.486E+00, 2.527E+00,
     &                               2.568E+00, 2.609E+00, 2.650E+00,
     &                               2.692E+00, 2.739E+00, 2.796E+00,
     &                               2.855E+00, 2.916E+00, 2.979E+00,
     &                               3.038E+00, 3.098E+00, 3.157E+00,
     &                               3.238E+00, 3.321E+00, 3.406E+00,
     &                               3.493E+00, 3.606E+00, 3.731E+00,
     &                               3.872E+00, 4.033E+00, 4.151E+00,
     &                               4.277E+00, 4.410E+00, 4.522E+00,
     &                               4.587E+00, 4.659E+00, 4.737E+00,
     &                               4.816E+00, 4.888E+00, 4.961E+00,
     &                               5.035E+00, 5.076E+00, 5.110E+00,
     &                               5.144E+00, 5.180E+00, 5.222E+00,
     &                               5.271E+00, 5.320E+00, 5.372E+00,
     &                               5.425E+00, 5.482E+00, 5.544E+00,
     &                               5.605E+00, 5.665E+00, 5.718E+00,
     &                               5.767E+00, 5.814E+00, 5.862E+00,
     &                               5.909E+00, 5.957E+00, 6.003E+00/
      DATA (UV_A_EST(40,I),I=1,60) / 2.092E+00, 2.125E+00, 2.160E+00,
     &                               2.194E+00, 2.230E+00, 2.266E+00,
     &                               2.303E+00, 2.345E+00, 2.387E+00,
     &                               2.430E+00, 2.473E+00, 2.515E+00,
     &                               2.557E+00, 2.607E+00, 2.663E+00,
     &                               2.720E+00, 2.780E+00, 2.841E+00,
     &                               2.904E+00, 2.969E+00, 3.045E+00,
     &                               3.126E+00, 3.208E+00, 3.292E+00,
     &                               3.378E+00, 3.467E+00, 3.575E+00,
     &                               3.715E+00, 3.896E+00, 4.052E+00,
     &                               4.176E+00, 4.308E+00, 4.448E+00,
     &                               4.542E+00, 4.618E+00, 4.697E+00,
     &                               4.778E+00, 4.853E+00, 4.926E+00,
     &                               5.001E+00, 5.062E+00, 5.097E+00,
     &                               5.132E+00, 5.168E+00, 5.208E+00,
     &                               5.257E+00, 5.306E+00, 5.358E+00,
     &                               5.410E+00, 5.463E+00, 5.523E+00,
     &                               5.583E+00, 5.641E+00, 5.697E+00,
     &                               5.744E+00, 5.791E+00, 5.838E+00,
     &                               5.885E+00, 5.933E+00, 5.980E+00/
      DATA (UV_A_EST(41,I),I=1,60) / 1.953E+00, 1.989E+00, 2.022E+00,
     &                               2.056E+00, 2.090E+00, 2.125E+00,
     &                               2.161E+00, 2.202E+00, 2.244E+00,
     &                               2.286E+00, 2.329E+00, 2.372E+00,
     &                               2.417E+00, 2.474E+00, 2.529E+00,
     &                               2.584E+00, 2.642E+00, 2.701E+00,
     &                               2.762E+00, 2.833E+00, 2.921E+00,
     &                               3.012E+00, 3.094E+00, 3.178E+00,
     &                               3.264E+00, 3.351E+00, 3.440E+00,
     &                               3.559E+00, 3.739E+00, 3.922E+00,
     &                               4.072E+00, 4.203E+00, 4.342E+00,
     &                               4.500E+00, 4.576E+00, 4.655E+00,
     &                               4.737E+00, 4.818E+00, 4.891E+00,
     &                               4.966E+00, 5.043E+00, 5.082E+00,
     &                               5.119E+00, 5.156E+00, 5.195E+00,
     &                               5.243E+00, 5.293E+00, 5.342E+00,
     &                               5.395E+00, 5.447E+00, 5.503E+00,
     &                               5.560E+00, 5.615E+00, 5.670E+00,
     &                               5.720E+00, 5.767E+00, 5.815E+00,
     &                               5.863E+00, 5.910E+00, 5.958E+00/
      DATA (UV_A_EST(42,I),I=1,60) / 1.807E+00, 1.842E+00, 1.878E+00,
     &                               1.913E+00, 1.948E+00, 1.984E+00,
     &                               2.020E+00, 2.060E+00, 2.100E+00,
     &                               2.141E+00, 2.183E+00, 2.226E+00,
     &                               2.275E+00, 2.330E+00, 2.387E+00,
     &                               2.445E+00, 2.502E+00, 2.560E+00,
     &                               2.619E+00, 2.698E+00, 2.783E+00,
     &                               2.874E+00, 2.972E+00, 3.061E+00,
     &                               3.146E+00, 3.233E+00, 3.329E+00,
     &                               3.440E+00, 3.583E+00, 3.765E+00,
     &                               3.949E+00, 4.094E+00, 4.232E+00,
     &                               4.401E+00, 4.533E+00, 4.613E+00,
     &                               4.695E+00, 4.779E+00, 4.856E+00,
     &                               4.932E+00, 5.008E+00, 5.067E+00,
     &                               5.105E+00, 5.143E+00, 5.183E+00,
     &                               5.228E+00, 5.278E+00, 5.328E+00,
     &                               5.379E+00, 5.431E+00, 5.485E+00,
     &                               5.539E+00, 5.593E+00, 5.646E+00,
     &                               5.698E+00, 5.745E+00, 5.792E+00,
     &                               5.840E+00, 5.887E+00, 5.934E+00/
      DATA (UV_A_EST(43,I),I=1,60) / 1.660E+00, 1.695E+00, 1.730E+00,
     &                               1.766E+00, 1.801E+00, 1.837E+00,
     &                               1.876E+00, 1.916E+00, 1.957E+00,
     &                               1.997E+00, 2.037E+00, 2.079E+00,
     &                               2.131E+00, 2.185E+00, 2.241E+00,
     &                               2.298E+00, 2.356E+00, 2.414E+00,
     &                               2.481E+00, 2.558E+00, 2.641E+00,
     &                               2.730E+00, 2.825E+00, 2.926E+00,
     &                               3.025E+00, 3.120E+00, 3.229E+00,
     &                               3.339E+00, 3.453E+00, 3.606E+00,
     &                               3.790E+00, 3.976E+00, 4.133E+00,
     &                               4.302E+00, 4.477E+00, 4.569E+00,
     &                               4.651E+00, 4.736E+00, 4.820E+00,
     &                               4.896E+00, 4.972E+00, 5.050E+00,
     &                               5.088E+00, 5.128E+00, 5.169E+00,
     &                               5.214E+00, 5.264E+00, 5.314E+00,
     &                               5.365E+00, 5.416E+00, 5.469E+00,
     &                               5.522E+00, 5.573E+00, 5.624E+00,
     &                               5.673E+00, 5.720E+00, 5.768E+00,
     &                               5.815E+00, 5.863E+00, 5.911E+00/
      DATA (UV_A_EST(44,I),I=1,60) / 1.514E+00, 1.549E+00, 1.584E+00,
     &                               1.619E+00, 1.654E+00, 1.689E+00,
     &                               1.729E+00, 1.769E+00, 1.810E+00,
     &                               1.850E+00, 1.891E+00, 1.935E+00,
     &                               1.987E+00, 2.040E+00, 2.094E+00,
     &                               2.150E+00, 2.206E+00, 2.265E+00,
     &                               2.339E+00, 2.416E+00, 2.494E+00,
     &                               2.579E+00, 2.671E+00, 2.769E+00,
     &                               2.881E+00, 3.018E+00, 3.127E+00,
     &                               3.239E+00, 3.352E+00, 3.467E+00,
     &                               3.629E+00, 3.814E+00, 4.032E+00,
     &                               4.202E+00, 4.376E+00, 4.525E+00,
     &                               4.608E+00, 4.693E+00, 4.780E+00,
     &                               4.858E+00, 4.935E+00, 5.012E+00,
     &                               5.070E+00, 5.112E+00, 5.155E+00,
     &                               5.199E+00, 5.250E+00, 5.301E+00,
     &                               5.352E+00, 5.402E+00, 5.452E+00,
     &                               5.503E+00, 5.554E+00, 5.603E+00,
     &                               5.651E+00, 5.698E+00, 5.744E+00,
     &                               5.792E+00, 5.839E+00, 5.886E+00/
      DATA (UV_A_EST(45,I),I=1,60) / 1.368E+00, 1.403E+00, 1.438E+00,
     &                               1.473E+00, 1.507E+00, 1.543E+00,
     &                               1.582E+00, 1.622E+00, 1.662E+00,
     &                               1.702E+00, 1.743E+00, 1.792E+00,
     &                               1.844E+00, 1.896E+00, 1.948E+00,
     &                               2.000E+00, 2.056E+00, 2.123E+00,
     &                               2.195E+00, 2.270E+00, 2.348E+00,
     &                               2.427E+00, 2.508E+00, 2.602E+00,
     &                               2.736E+00, 2.879E+00, 3.023E+00,
     &                               3.135E+00, 3.249E+00, 3.364E+00,
     &                               3.481E+00, 3.661E+00, 3.895E+00,
     &                               4.100E+00, 4.275E+00, 4.455E+00,
     &                               4.563E+00, 4.648E+00, 4.736E+00,
     &                               4.821E+00, 4.897E+00, 4.973E+00,
     &                               5.049E+00, 5.094E+00, 5.140E+00,
     &                               5.186E+00, 5.236E+00, 5.287E+00,
     &                               5.336E+00, 5.388E+00, 5.439E+00,
     &                               5.487E+00, 5.535E+00, 5.583E+00,
     &                               5.629E+00, 5.675E+00, 5.721E+00,
     &                               5.768E+00, 5.816E+00, 5.863E+00/
      DATA (UV_A_EST(46,I),I=1,60) / 1.220E+00, 1.256E+00, 1.290E+00,
     &                               1.325E+00, 1.360E+00, 1.397E+00,
     &                               1.436E+00, 1.476E+00, 1.514E+00,
     &                               1.554E+00, 1.597E+00, 1.648E+00,
     &                               1.699E+00, 1.751E+00, 1.803E+00,
     &                               1.855E+00, 1.914E+00, 1.980E+00,
     &                               2.049E+00, 2.122E+00, 2.199E+00,
     &                               2.277E+00, 2.357E+00, 2.459E+00,
     &                               2.584E+00, 2.727E+00, 2.880E+00,
     &                               3.028E+00, 3.142E+00, 3.259E+00,
     &                               3.377E+00, 3.506E+00, 3.736E+00,
     &                               3.994E+00, 4.171E+00, 4.351E+00,
     &                               4.515E+00, 4.601E+00, 4.691E+00,
     &                               4.782E+00, 4.860E+00, 4.936E+00,
     &                               5.012E+00, 5.072E+00, 5.121E+00,
     &                               5.171E+00, 5.222E+00, 5.272E+00,
     &                               5.323E+00, 5.373E+00, 5.424E+00,
     &                               5.473E+00, 5.519E+00, 5.565E+00,
     &                               5.610E+00, 5.654E+00, 5.697E+00,
     &                               5.745E+00, 5.792E+00, 5.839E+00/
      DATA (UV_A_EST(47,I),I=1,60) / 1.073E+00, 1.107E+00, 1.142E+00,
     &                               1.177E+00, 1.212E+00, 1.250E+00,
     &                               1.289E+00, 1.329E+00, 1.368E+00,
     &                               1.407E+00, 1.455E+00, 1.503E+00,
     &                               1.554E+00, 1.605E+00, 1.656E+00,
     &                               1.711E+00, 1.775E+00, 1.841E+00,
     &                               1.907E+00, 1.974E+00, 2.045E+00,
     &                               2.122E+00, 2.215E+00, 2.321E+00,
     &                               2.435E+00, 2.564E+00, 2.716E+00,
     &                               2.879E+00, 3.033E+00, 3.153E+00,
     &                               3.278E+00, 3.407E+00, 3.570E+00,
     &                               3.825E+00, 4.066E+00, 4.246E+00,
     &                               4.430E+00, 4.557E+00, 4.648E+00,
     &                               4.738E+00, 4.822E+00, 4.899E+00,
     &                               4.975E+00, 5.054E+00, 5.106E+00,
     &                               5.158E+00, 5.210E+00, 5.260E+00,
     &                               5.309E+00, 5.360E+00, 5.411E+00,
     &                               5.460E+00, 5.505E+00, 5.549E+00,
     &                               5.592E+00, 5.634E+00, 5.677E+00,
     &                               5.722E+00, 5.769E+00, 5.817E+00/
      DATA (UV_A_EST(48,I),I=1,60) / 9.261E-01, 9.605E-01, 9.952E-01,
     &                               1.029E+00, 1.064E+00, 1.103E+00,
     &                               1.141E+00, 1.181E+00, 1.220E+00,
     &                               1.263E+00, 1.312E+00, 1.361E+00,
     &                               1.410E+00, 1.460E+00, 1.509E+00,
     &                               1.571E+00, 1.634E+00, 1.699E+00,
     &                               1.765E+00, 1.832E+00, 1.899E+00,
     &                               1.978E+00, 2.071E+00, 2.175E+00,
     &                               2.288E+00, 2.407E+00, 2.539E+00,
     &                               2.701E+00, 2.885E+00, 3.050E+00,
     &                               3.177E+00, 3.307E+00, 3.442E+00,
     &                               3.647E+00, 3.929E+00, 4.138E+00,
     &                               4.332E+00, 4.516E+00, 4.607E+00,
     &                               4.696E+00, 4.784E+00, 4.863E+00,
     &                               4.944E+00, 5.030E+00, 5.091E+00,
     &                               5.145E+00, 5.199E+00, 5.248E+00,
     &                               5.296E+00, 5.346E+00, 5.397E+00,
     &                               5.449E+00, 5.492E+00, 5.534E+00,
     &                               5.576E+00, 5.617E+00, 5.658E+00,
     &                               5.698E+00, 5.746E+00, 5.793E+00/
      DATA (UV_A_EST(49,I),I=1,60) / 7.785E-01, 8.120E-01, 8.459E-01,
     &                               8.803E-01, 9.157E-01, 9.551E-01,
     &                               9.947E-01, 1.033E+00, 1.072E+00,
     &                               1.119E+00, 1.167E+00, 1.216E+00,
     &                               1.265E+00, 1.315E+00, 1.373E+00,
     &                               1.433E+00, 1.493E+00, 1.555E+00,
     &                               1.620E+00, 1.686E+00, 1.764E+00,
     &                               1.847E+00, 1.931E+00, 2.019E+00,
     &                               2.130E+00, 2.249E+00, 2.376E+00,
     &                               2.525E+00, 2.716E+00, 2.911E+00,
     &                               3.071E+00, 3.204E+00, 3.340E+00,
     &                               3.478E+00, 3.737E+00, 4.037E+00,
     &                               4.237E+00, 4.443E+00, 4.564E+00,
     &                               4.653E+00, 4.741E+00, 4.826E+00,
     &                               4.911E+00, 5.000E+00, 5.075E+00,
     &                               5.131E+00, 5.187E+00, 5.238E+00,
     &                               5.285E+00, 5.333E+00, 5.384E+00,
     &                               5.436E+00, 5.479E+00, 5.520E+00,
     &                               5.561E+00, 5.600E+00, 5.640E+00,
     &                               5.679E+00, 5.721E+00, 5.769E+00/
      DATA (UV_A_EST(50,I),I=1,60) / 6.300E-01, 6.625E-01, 6.955E-01,
     &                               7.289E-01, 7.651E-01, 8.036E-01,
     &                               8.426E-01, 8.820E-01, 9.252E-01,
     &                               9.743E-01, 1.022E+00, 1.070E+00,
     &                               1.119E+00, 1.175E+00, 1.233E+00,
     &                               1.294E+00, 1.353E+00, 1.414E+00,
     &                               1.476E+00, 1.551E+00, 1.630E+00,
     &                               1.711E+00, 1.796E+00, 1.881E+00,
     &                               1.969E+00, 2.075E+00, 2.221E+00,
     &                               2.376E+00, 2.543E+00, 2.739E+00,
     &                               2.941E+00, 3.095E+00, 3.232E+00,
     &                               3.373E+00, 3.531E+00, 3.888E+00,
     &                               4.140E+00, 4.348E+00, 4.523E+00,
     &                               4.612E+00, 4.699E+00, 4.785E+00,
     &                               4.874E+00, 4.967E+00, 5.059E+00,
     &                               5.117E+00, 5.175E+00, 5.227E+00,
     &                               5.273E+00, 5.320E+00, 5.371E+00,
     &                               5.422E+00, 5.467E+00, 5.506E+00,
     &                               5.545E+00, 5.585E+00, 5.623E+00,
     &                               5.662E+00, 5.700E+00, 5.747E+00/
      DATA (UV_A_EST(51,I),I=1,60) / 4.822E-01, 5.137E-01, 5.454E-01,
     &                               5.779E-01, 6.132E-01, 6.507E-01,
     &                               6.888E-01, 7.274E-01, 7.736E-01,
     &                               8.213E-01, 8.701E-01, 9.201E-01,
     &                               9.771E-01, 1.035E+00, 1.093E+00,
     &                               1.152E+00, 1.211E+00, 1.276E+00,
     &                               1.347E+00, 1.419E+00, 1.493E+00,
     &                               1.572E+00, 1.655E+00, 1.740E+00,
     &                               1.835E+00, 1.941E+00, 2.069E+00,
     &                               2.225E+00, 2.388E+00, 2.564E+00,
     &                               2.764E+00, 2.969E+00, 3.117E+00,
     &                               3.270E+00, 3.432E+00, 3.709E+00,
     &                               4.041E+00, 4.248E+00, 4.460E+00,
     &                               4.571E+00, 4.658E+00, 4.745E+00,
     &                               4.833E+00, 4.931E+00, 5.034E+00,
     &                               5.101E+00, 5.161E+00, 5.216E+00,
     &                               5.261E+00, 5.307E+00, 5.356E+00,
     &                               5.407E+00, 5.456E+00, 5.495E+00,
     &                               5.533E+00, 5.571E+00, 5.609E+00,
     &                               5.647E+00, 5.684E+00, 5.727E+00/
      DATA (UV_A_EST(52,I),I=1,60) / 3.268E-01, 3.586E-01, 3.909E-01,
     &                               4.237E-01, 4.605E-01, 4.983E-01,
     &                               5.350E-01, 5.745E-01, 6.191E-01,
     &                               6.650E-01, 7.122E-01, 7.651E-01,
     &                               8.239E-01, 8.846E-01, 9.468E-01,
     &                               1.009E+00, 1.074E+00, 1.143E+00,
     &                               1.214E+00, 1.287E+00, 1.361E+00,
     &                               1.435E+00, 1.511E+00, 1.609E+00,
     &                               1.712E+00, 1.820E+00, 1.932E+00,
     &                               2.068E+00, 2.231E+00, 2.399E+00,
     &                               2.581E+00, 2.797E+00, 3.020E+00,
     &                               3.181E+00, 3.345E+00, 3.524E+00,
     &                               3.892E+00, 4.148E+00, 4.361E+00,
     &                               4.532E+00, 4.618E+00, 4.704E+00,
     &                               4.789E+00, 4.888E+00, 4.998E+00,
     &                               5.082E+00, 5.145E+00, 5.206E+00,
     &                               5.250E+00, 5.295E+00, 5.342E+00,
     &                               5.394E+00, 5.445E+00, 5.484E+00,
     &                               5.521E+00, 5.558E+00, 5.595E+00,
     &                               5.632E+00, 5.669E+00, 5.707E+00/
      DATA (UV_A_EST(53,I),I=1,60) / 1.704E-01, 2.009E-01, 2.320E-01,
     &                               2.647E-01, 3.008E-01, 3.376E-01,
     &                               3.752E-01, 4.177E-01, 4.624E-01,
     &                               5.070E-01, 5.534E-01, 6.087E-01,
     &                               6.659E-01, 7.249E-01, 7.870E-01,
     &                               8.582E-01, 9.333E-01, 1.008E+00,
     &                               1.079E+00, 1.151E+00, 1.224E+00,
     &                               1.307E+00, 1.394E+00, 1.480E+00,
     &                               1.581E+00, 1.691E+00, 1.806E+00,
     &                               1.925E+00, 2.063E+00, 2.236E+00,
     &                               2.422E+00, 2.638E+00, 2.878E+00,
     &                               3.081E+00, 3.249E+00, 3.422E+00,
     &                               3.707E+00, 4.048E+00, 4.263E+00,
     &                               4.481E+00, 4.577E+00, 4.664E+00,
     &                               4.752E+00, 4.848E+00, 4.960E+00,
     &                               5.063E+00, 5.128E+00, 5.194E+00,
     &                               5.240E+00, 5.284E+00, 5.328E+00,
     &                               5.380E+00, 5.431E+00, 5.473E+00,
     &                               5.510E+00, 5.546E+00, 5.582E+00,
     &                               5.619E+00, 5.655E+00, 5.691E+00/
      DATA (UV_A_EST(54,I),I=1,60) / 1.516E-02, 4.392E-02, 7.350E-02,
     &                               1.047E-01, 1.391E-01, 1.744E-01,
     &                               2.116E-01, 2.540E-01, 2.968E-01,
     &                               3.408E-01, 3.933E-01, 4.483E-01,
     &                               5.041E-01, 5.613E-01, 6.302E-01,
     &                               7.019E-01, 7.773E-01, 8.538E-01,
     &                               9.326E-01, 1.019E+00, 1.101E+00,
     &                               1.184E+00, 1.271E+00, 1.359E+00,
     &                               1.447E+00, 1.544E+00, 1.662E+00,
     &                               1.792E+00, 1.929E+00, 2.088E+00,
     &                               2.277E+00, 2.471E+00, 2.706E+00,
     &                               2.963E+00, 3.149E+00, 3.330E+00,
     &                               3.524E+00, 3.898E+00, 4.157E+00,
     &                               4.386E+00, 4.544E+00, 4.631E+00,
     &                               4.718E+00, 4.807E+00, 4.921E+00,
     &                               5.038E+00, 5.110E+00, 5.179E+00,
     &                               5.230E+00, 5.273E+00, 5.315E+00,
     &                               5.366E+00, 5.418E+00, 5.463E+00,
     &                               5.499E+00, 5.535E+00, 5.570E+00,
     &                               5.606E+00, 5.641E+00, 5.676E+00/
      DATA (UV_A_EST(55,I),I=1,60) /-1.448E-01,-1.154E-01,-8.621E-02,
     &                              -5.495E-02,-2.196E-02, 1.053E-02,
     &                               4.723E-02, 8.629E-02, 1.267E-01,
     &                               1.736E-01, 2.242E-01, 2.781E-01,
     &                               3.345E-01, 4.006E-01, 4.697E-01,
     &                               5.394E-01, 6.125E-01, 6.932E-01,
     &                               7.807E-01, 8.736E-01, 9.684E-01,
     &                               1.055E+00, 1.142E+00, 1.230E+00,
     &                               1.326E+00, 1.427E+00, 1.535E+00,
     &                               1.670E+00, 1.809E+00, 1.947E+00,
     &                               2.117E+00, 2.314E+00, 2.521E+00,
     &                               2.791E+00, 3.044E+00, 3.228E+00,
     &                               3.417E+00, 3.721E+00, 4.071E+00,
     &                               4.302E+00, 4.512E+00, 4.599E+00,
     &                               4.686E+00, 4.772E+00, 4.880E+00,
     &                               5.000E+00, 5.091E+00, 5.162E+00,
     &                               5.220E+00, 5.262E+00, 5.304E+00,
     &                               5.352E+00, 5.404E+00, 5.454E+00,
     &                               5.489E+00, 5.524E+00, 5.558E+00,
     &                               5.593E+00, 5.628E+00, 5.662E+00/
      DATA (UV_A_EST(56,I),I=1,60) /-3.066E-01,-2.782E-01,-2.488E-01,
     &                              -2.159E-01,-1.828E-01,-1.487E-01,
     &                              -1.129E-01,-7.724E-02,-3.869E-02,
     &                               1.580E-03, 4.884E-02, 1.018E-01,
     &                               1.624E-01, 2.265E-01, 2.953E-01,
     &                               3.685E-01, 4.491E-01, 5.311E-01,
     &                               6.188E-01, 7.126E-01, 8.084E-01,
     &                               9.121E-01, 1.018E+00, 1.116E+00,
     &                               1.217E+00, 1.323E+00, 1.429E+00,
     &                               1.543E+00, 1.682E+00, 1.825E+00,
     &                               1.967E+00, 2.156E+00, 2.363E+00,
     &                               2.592E+00, 2.891E+00, 3.127E+00,
     &                               3.328E+00, 3.566E+00, 3.974E+00,
     &                               4.217E+00, 4.450E+00, 4.568E+00,
     &                               4.654E+00, 4.740E+00, 4.837E+00,
     &                               4.960E+00, 5.070E+00, 5.145E+00,
     &                               5.211E+00, 5.252E+00, 5.293E+00,
     &                               5.338E+00, 5.390E+00, 5.442E+00,
     &                               5.479E+00, 5.513E+00, 5.547E+00,
     &                               5.581E+00, 5.615E+00, 5.650E+00/
      DATA (UV_A_EST(57,I),I=1,60) /-4.681E-01,-4.403E-01,-4.119E-01,
     &                              -3.781E-01,-3.442E-01,-3.088E-01,
     &                              -2.732E-01,-2.339E-01,-1.944E-01,
     &                              -1.544E-01,-1.100E-01,-6.397E-02,
     &                              -1.688E-02, 4.018E-02, 1.126E-01,
     &                               1.892E-01, 2.716E-01, 3.558E-01,
     &                               4.450E-01, 5.429E-01, 6.496E-01,
     &                               7.604E-01, 8.779E-01, 9.959E-01,
     &                               1.099E+00, 1.209E+00, 1.318E+00,
     &                               1.430E+00, 1.551E+00, 1.695E+00,
     &                               1.847E+00, 1.997E+00, 2.213E+00,
     &                               2.441E+00, 2.731E+00, 3.031E+00,
     &                               3.240E+00, 3.452E+00, 3.819E+00,
     &                               4.131E+00, 4.365E+00, 4.537E+00,
     &                               4.623E+00, 4.708E+00, 4.794E+00,
     &                               4.917E+00, 5.047E+00, 5.125E+00,
     &                               5.203E+00, 5.242E+00, 5.283E+00,
     &                               5.325E+00, 5.376E+00, 5.428E+00,
     &                               5.469E+00, 5.503E+00, 5.536E+00,
     &                               5.570E+00, 5.603E+00, 5.637E+00/
      DATA (UV_A_EST(58,I),I=1,60) /-6.335E-01,-6.073E-01,-5.782E-01,
     &                              -5.423E-01,-5.029E-01,-4.677E-01,
     &                              -4.309E-01,-3.929E-01,-3.533E-01,
     &                              -3.111E-01,-2.672E-01,-2.184E-01,
     &                              -1.672E-01,-1.146E-01,-6.083E-02,
     &                              -5.127E-03, 7.830E-02, 1.716E-01,
     &                               2.697E-01, 3.752E-01, 4.813E-01,
     &                               5.928E-01, 7.124E-01, 8.406E-01,
     &                               9.723E-01, 1.089E+00, 1.206E+00,
     &                               1.323E+00, 1.444E+00, 1.577E+00,
     &                               1.731E+00, 1.886E+00, 2.065E+00,
     &                               2.303E+00, 2.559E+00, 2.887E+00,
     &                               3.144E+00, 3.364E+00, 3.661E+00,
     &                               4.044E+00, 4.280E+00, 4.506E+00,
     &                               4.591E+00, 4.677E+00, 4.764E+00,
     &                               4.872E+00, 5.006E+00, 5.103E+00,
     &                               5.188E+00, 5.233E+00, 5.273E+00,
     &                               5.313E+00, 5.362E+00, 5.414E+00,
     &                               5.460E+00, 5.493E+00, 5.526E+00,
     &                               5.558E+00, 5.591E+00, 5.625E+00/
      DATA (UV_A_EST(59,I),I=1,60) /-8.005E-01,-7.778E-01,-7.521E-01,
     &                              -7.214E-01,-6.876E-01,-6.508E-01,
     &                              -6.101E-01,-5.672E-01,-5.206E-01,
     &                              -4.737E-01,-4.266E-01,-3.775E-01,
     &                              -3.243E-01,-2.681E-01,-2.081E-01,
     &                              -1.484E-01,-8.713E-02,-2.351E-02,
     &                               6.890E-02, 1.836E-01, 3.044E-01,
     &                               4.251E-01, 5.497E-01, 6.880E-01,
     &                               8.291E-01, 9.757E-01, 1.097E+00,
     &                               1.217E+00, 1.340E+00, 1.461E+00,
     &                               1.607E+00, 1.770E+00, 1.937E+00,
     &                               2.150E+00, 2.403E+00, 2.713E+00,
     &                               3.039E+00, 3.269E+00, 3.501E+00,
     &                               3.921E+00, 4.193E+00, 4.436E+00,
     &                               4.560E+00, 4.645E+00, 4.732E+00,
     &                               4.824E+00, 4.962E+00, 5.079E+00,
     &                               5.168E+00, 5.224E+00, 5.263E+00,
     &                               5.303E+00, 5.348E+00, 5.400E+00,
     &                               5.452E+00, 5.484E+00, 5.516E+00,
     &                               5.548E+00, 5.580E+00, 5.613E+00/
      DATA (UV_A_EST(60,I),I=1,60) /-9.678E-01,-9.477E-01,-9.264E-01,
     &                              -9.002E-01,-8.720E-01,-8.406E-01,
     &                              -8.064E-01,-7.693E-01,-7.301E-01,
     &                              -6.806E-01,-6.243E-01,-5.581E-01,
     &                              -4.935E-01,-4.323E-01,-3.679E-01,
     &                              -3.015E-01,-2.318E-01,-1.625E-01,
     &                              -9.217E-02,-2.082E-02, 9.483E-02,
     &                               2.312E-01, 3.701E-01, 5.037E-01,
     &                               6.613E-01, 8.221E-01, 9.796E-01,
     &                               1.105E+00, 1.229E+00, 1.356E+00,
     &                               1.481E+00, 1.642E+00, 1.817E+00,
     &                               1.989E+00, 2.246E+00, 2.521E+00,
     &                               2.891E+00, 3.165E+00, 3.412E+00,
     &                               3.771E+00, 4.105E+00, 4.350E+00,
     &                               4.530E+00, 4.614E+00, 4.700E+00,
     &                               4.784E+00, 4.911E+00, 5.050E+00,
     &                               5.143E+00, 5.214E+00, 5.253E+00,
     &                               5.292E+00, 5.333E+00, 5.386E+00,
     &                               5.439E+00, 5.474E+00, 5.506E+00,
     &                               5.537E+00, 5.569E+00, 5.601E+00/
      DATA (UV_B_EST(1,I),I=1,60) /  1.680E-04, 1.038E-02, 2.047E-02,
     &                               3.044E-02, 4.029E-02, 5.002E-02,
     &                               5.964E-02, 6.914E-02, 7.853E-02,
     &                               8.781E-02, 9.699E-02, 1.061E-01,
     &                               1.150E-01, 1.239E-01, 1.343E-01,
     &                               1.447E-01, 1.550E-01, 1.651E-01,
     &                               1.751E-01, 1.849E-01, 1.946E-01,
     &                               2.040E-01, 2.134E-01, 2.225E-01,
     &                               2.315E-01, 2.404E-01, 2.491E-01,
     &                               2.578E-01, 2.663E-01, 2.746E-01,
     &                               2.827E-01, 2.907E-01, 2.984E-01,
     &                               3.060E-01, 3.134E-01, 3.207E-01,
     &                               3.279E-01, 3.350E-01, 3.419E-01,
     &                               3.487E-01, 3.553E-01, 3.619E-01,
     &                               3.683E-01, 3.746E-01, 3.813E-01,
     &                               3.878E-01, 3.941E-01, 4.003E-01,
     &                               4.063E-01, 4.121E-01, 4.178E-01,
     &                               4.234E-01, 4.289E-01, 4.342E-01,
     &                               4.394E-01, 4.445E-01, 4.495E-01,
     &                               4.544E-01, 4.591E-01, 4.639E-01/
      DATA (UV_B_EST(2,I),I=1,60) /  1.671E-04, 1.056E-02, 2.082E-02,
     &                               3.097E-02, 4.099E-02, 5.090E-02,
     &                               6.069E-02, 7.035E-02, 7.989E-02,
     &                               8.933E-02, 9.865E-02, 1.079E-01,
     &                               1.170E-01, 1.262E-01, 1.369E-01,
     &                               1.474E-01, 1.578E-01, 1.681E-01,
     &                               1.781E-01, 1.880E-01, 1.978E-01,
     &                               2.073E-01, 2.167E-01, 2.260E-01,
     &                               2.351E-01, 2.440E-01, 2.529E-01,
     &                               2.617E-01, 2.703E-01, 2.787E-01,
     &                               2.869E-01, 2.949E-01, 3.027E-01,
     &                               3.103E-01, 3.178E-01, 3.251E-01,
     &                               3.323E-01, 3.394E-01, 3.464E-01,
     &                               3.532E-01, 3.599E-01, 3.665E-01,
     &                               3.729E-01, 3.797E-01, 3.864E-01,
     &                               3.930E-01, 3.993E-01, 4.055E-01,
     &                               4.115E-01, 4.173E-01, 4.231E-01,
     &                               4.286E-01, 4.341E-01, 4.394E-01,
     &                               4.446E-01, 4.497E-01, 4.547E-01,
     &                               4.596E-01, 4.644E-01, 4.691E-01/
      DATA (UV_B_EST(3,I),I=1,60) /  1.661E-04, 1.074E-02, 2.119E-02,
     &                               3.150E-02, 4.167E-02, 5.177E-02,
     &                               6.173E-02, 7.157E-02, 8.130E-02,
     &                               9.090E-02, 1.004E-01, 1.097E-01,
     &                               1.190E-01, 1.288E-01, 1.396E-01,
     &                               1.502E-01, 1.607E-01, 1.711E-01,
     &                               1.812E-01, 1.912E-01, 2.010E-01,
     &                               2.107E-01, 2.202E-01, 2.295E-01,
     &                               2.387E-01, 2.477E-01, 2.568E-01,
     &                               2.657E-01, 2.744E-01, 2.829E-01,
     &                               2.911E-01, 2.992E-01, 3.071E-01,
     &                               3.148E-01, 3.223E-01, 3.297E-01,
     &                               3.368E-01, 3.439E-01, 3.510E-01,
     &                               3.578E-01, 3.646E-01, 3.712E-01,
     &                               3.780E-01, 3.850E-01, 3.918E-01,
     &                               3.983E-01, 4.047E-01, 4.108E-01,
     &                               4.168E-01, 4.227E-01, 4.284E-01,
     &                               4.340E-01, 4.395E-01, 4.448E-01,
     &                               4.500E-01, 4.551E-01, 4.601E-01,
     &                               4.650E-01, 4.698E-01, 4.745E-01/
      DATA (UV_B_EST(4,I),I=1,60) /  1.652E-04, 1.093E-02, 2.156E-02,
     &                               3.206E-02, 4.243E-02, 5.268E-02,
     &                               6.280E-02, 7.280E-02, 8.268E-02,
     &                               9.245E-02, 1.021E-01, 1.117E-01,
     &                               1.211E-01, 1.314E-01, 1.423E-01,
     &                               1.531E-01, 1.637E-01, 1.741E-01,
     &                               1.844E-01, 1.945E-01, 2.044E-01,
     &                               2.141E-01, 2.237E-01, 2.332E-01,
     &                               2.424E-01, 2.516E-01, 2.609E-01,
     &                               2.699E-01, 2.787E-01, 2.872E-01,
     &                               2.955E-01, 3.037E-01, 3.116E-01,
     &                               3.194E-01, 3.270E-01, 3.344E-01,
     &                               3.416E-01, 3.487E-01, 3.557E-01,
     &                               3.626E-01, 3.694E-01, 3.761E-01,
     &                               3.834E-01, 3.905E-01, 3.972E-01,
     &                               4.038E-01, 4.102E-01, 4.164E-01,
     &                               4.224E-01, 4.282E-01, 4.339E-01,
     &                               4.395E-01, 4.450E-01, 4.503E-01,
     &                               4.555E-01, 4.606E-01, 4.656E-01,
     &                               4.705E-01, 4.753E-01, 4.800E-01/
      DATA (UV_B_EST(5,I),I=1,60) /  1.642E-04, 1.113E-02, 2.195E-02,
     &                               3.264E-02, 4.319E-02, 5.361E-02,
     &                               6.390E-02, 7.407E-02, 8.412E-02,
     &                               9.405E-02, 1.039E-01, 1.136E-01,
     &                               1.231E-01, 1.341E-01, 1.452E-01,
     &                               1.560E-01, 1.667E-01, 1.773E-01,
     &                               1.876E-01, 1.978E-01, 2.078E-01,
     &                               2.177E-01, 2.274E-01, 2.369E-01,
     &                               2.462E-01, 2.557E-01, 2.651E-01,
     &                               2.742E-01, 2.831E-01, 2.917E-01,
     &                               3.001E-01, 3.083E-01, 3.163E-01,
     &                               3.241E-01, 3.318E-01, 3.392E-01,
     &                               3.465E-01, 3.536E-01, 3.606E-01,
     &                               3.675E-01, 3.743E-01, 3.818E-01,
     &                               3.891E-01, 3.961E-01, 4.029E-01,
     &                               4.095E-01, 4.159E-01, 4.220E-01,
     &                               4.280E-01, 4.339E-01, 4.396E-01,
     &                               4.451E-01, 4.506E-01, 4.559E-01,
     &                               4.611E-01, 4.662E-01, 4.712E-01,
     &                               4.761E-01, 4.809E-01, 4.856E-01/
      DATA (UV_B_EST(6,I),I=1,60) /  1.632E-04, 1.133E-02, 2.236E-02,
     &                               3.323E-02, 4.397E-02, 5.458E-02,
     &                               6.505E-02, 7.539E-02, 8.560E-02,
     &                               9.569E-02, 1.057E-01, 1.155E-01,
     &                               1.253E-01, 1.367E-01, 1.480E-01,
     &                               1.590E-01, 1.699E-01, 1.805E-01,
     &                               1.910E-01, 2.013E-01, 2.114E-01,
     &                               2.213E-01, 2.311E-01, 2.407E-01,
     &                               2.501E-01, 2.600E-01, 2.695E-01,
     &                               2.787E-01, 2.877E-01, 2.964E-01,
     &                               3.048E-01, 3.130E-01, 3.211E-01,
     &                               3.290E-01, 3.367E-01, 3.442E-01,
     &                               3.516E-01, 3.587E-01, 3.658E-01,
     &                               3.726E-01, 3.800E-01, 3.876E-01,
     &                               3.949E-01, 4.019E-01, 4.087E-01,
     &                               4.153E-01, 4.217E-01, 4.279E-01,
     &                               4.339E-01, 4.397E-01, 4.454E-01,
     &                               4.509E-01, 4.564E-01, 4.617E-01,
     &                               4.669E-01, 4.720E-01, 4.770E-01,
     &                               4.819E-01, 4.867E-01, 4.914E-01/
      DATA (UV_B_EST(7,I),I=1,60) /  1.611E-04, 1.152E-02, 2.277E-02,
     &                               3.385E-02, 4.479E-02, 5.558E-02,
     &                               6.623E-02, 7.675E-02, 8.714E-02,
     &                               9.740E-02, 1.075E-01, 1.175E-01,
     &                               1.279E-01, 1.394E-01, 1.508E-01,
     &                               1.620E-01, 1.730E-01, 1.838E-01,
     &                               1.944E-01, 2.048E-01, 2.150E-01,
     &                               2.250E-01, 2.349E-01, 2.446E-01,
     &                               2.545E-01, 2.646E-01, 2.742E-01,
     &                               2.834E-01, 2.924E-01, 3.012E-01,
     &                               3.096E-01, 3.179E-01, 3.261E-01,
     &                               3.340E-01, 3.418E-01, 3.493E-01,
     &                               3.567E-01, 3.640E-01, 3.710E-01,
     &                               3.784E-01, 3.862E-01, 3.936E-01,
     &                               4.009E-01, 4.080E-01, 4.148E-01,
     &                               4.213E-01, 4.277E-01, 4.339E-01,
     &                               4.399E-01, 4.457E-01, 4.514E-01,
     &                               4.569E-01, 4.623E-01, 4.676E-01,
     &                               4.728E-01, 4.779E-01, 4.829E-01,
     &                               4.878E-01, 4.925E-01, 4.973E-01/
      DATA (UV_B_EST(8,I),I=1,60) /  1.576E-04, 1.176E-02, 2.320E-02,
     &                               3.449E-02, 4.554E-02, 5.662E-02,
     &                               6.746E-02, 7.817E-02, 8.873E-02,
     &                               9.917E-02, 1.095E-01, 1.196E-01,
     &                               1.305E-01, 1.422E-01, 1.537E-01,
     &                               1.650E-01, 1.761E-01, 1.871E-01,
     &                               1.979E-01, 2.084E-01, 2.187E-01,
     &                               2.289E-01, 2.388E-01, 2.486E-01,
     &                               2.591E-01, 2.693E-01, 2.791E-01,
     &                               2.884E-01, 2.975E-01, 3.061E-01,
     &                               3.146E-01, 3.230E-01, 3.311E-01,
     &                               3.392E-01, 3.470E-01, 3.546E-01,
     &                               3.621E-01, 3.694E-01, 3.767E-01,
     &                               3.848E-01, 3.925E-01, 3.999E-01,
     &                               4.071E-01, 4.142E-01, 4.210E-01,
     &                               4.275E-01, 4.339E-01, 4.401E-01,
     &                               4.461E-01, 4.519E-01, 4.575E-01,
     &                               4.631E-01, 4.684E-01, 4.737E-01,
     &                               4.789E-01, 4.839E-01, 4.889E-01,
     &                               4.938E-01, 4.985E-01, 5.036E-01/
      DATA (UV_B_EST(9,I),I=1,60) /  1.539E-04, 1.199E-02, 2.365E-02,
     &                               3.518E-02, 4.651E-02, 5.770E-02,
     &                               6.874E-02, 7.964E-02, 9.039E-02,
     &                               1.010E-01, 1.115E-01, 1.218E-01,
     &                               1.332E-01, 1.450E-01, 1.567E-01,
     &                               1.681E-01, 1.793E-01, 1.904E-01,
     &                               2.013E-01, 2.120E-01, 2.225E-01,
     &                               2.328E-01, 2.429E-01, 2.531E-01,
     &                               2.640E-01, 2.743E-01, 2.842E-01,
     &                               2.937E-01, 3.027E-01, 3.113E-01,
     &                               3.197E-01, 3.281E-01, 3.364E-01,
     &                               3.445E-01, 3.524E-01, 3.601E-01,
     &                               3.676E-01, 3.749E-01, 3.833E-01,
     &                               3.913E-01, 3.991E-01, 4.065E-01,
     &                               4.136E-01, 4.206E-01, 4.273E-01,
     &                               4.339E-01, 4.403E-01, 4.464E-01,
     &                               4.524E-01, 4.582E-01, 4.638E-01,
     &                               4.694E-01, 4.747E-01, 4.800E-01,
     &                               4.851E-01, 4.901E-01, 4.951E-01,
     &                               4.999E-01, 5.053E-01, 5.105E-01/
      DATA (UV_B_EST(10,I),I=1,60) / 1.512E-04, 1.222E-02, 2.412E-02,
     &                               3.585E-02, 4.741E-02, 5.881E-02,
     &                               7.005E-02, 8.115E-02, 9.213E-02,
     &                               1.029E-01, 1.136E-01, 1.241E-01,
     &                               1.360E-01, 1.480E-01, 1.597E-01,
     &                               1.713E-01, 1.826E-01, 1.938E-01,
     &                               2.048E-01, 2.157E-01, 2.263E-01,
     &                               2.368E-01, 2.470E-01, 2.581E-01,
     &                               2.691E-01, 2.796E-01, 2.896E-01,
     &                               2.992E-01, 3.081E-01, 3.169E-01,
     &                               3.254E-01, 3.337E-01, 3.418E-01,
     &                               3.499E-01, 3.579E-01, 3.657E-01,
     &                               3.733E-01, 3.817E-01, 3.901E-01,
     &                               3.981E-01, 4.058E-01, 4.132E-01,
     &                               4.203E-01, 4.272E-01, 4.339E-01,
     &                               4.405E-01, 4.468E-01, 4.529E-01,
     &                               4.589E-01, 4.647E-01, 4.703E-01,
     &                               4.758E-01, 4.812E-01, 4.864E-01,
     &                               4.915E-01, 4.965E-01, 5.015E-01,
     &                               5.070E-01, 5.123E-01, 5.174E-01/
      DATA (UV_B_EST(11,I),I=1,60) / 1.438E-04, 1.247E-02, 2.460E-02,
     &                               3.656E-02, 4.835E-02, 5.996E-02,
     &                               7.142E-02, 8.271E-02, 9.385E-02,
     &                               1.048E-01, 1.157E-01, 1.266E-01,
     &                               1.389E-01, 1.510E-01, 1.628E-01,
     &                               1.745E-01, 1.860E-01, 1.973E-01,
     &                               2.084E-01, 2.194E-01, 2.301E-01,
     &                               2.407E-01, 2.514E-01, 2.633E-01,
     &                               2.745E-01, 2.851E-01, 2.952E-01,
     &                               3.048E-01, 3.139E-01, 3.227E-01,
     &                               3.312E-01, 3.396E-01, 3.478E-01,
     &                               3.557E-01, 3.636E-01, 3.715E-01,
     &                               3.799E-01, 3.886E-01, 3.971E-01,
     &                               4.051E-01, 4.128E-01, 4.201E-01,
     &                               4.272E-01, 4.341E-01, 4.407E-01,
     &                               4.472E-01, 4.535E-01, 4.596E-01,
     &                               4.656E-01, 4.713E-01, 4.770E-01,
     &                               4.824E-01, 4.878E-01, 4.930E-01,
     &                               4.980E-01, 5.033E-01, 5.087E-01,
     &                               5.141E-01, 5.194E-01, 5.245E-01/
      DATA (UV_B_EST(12,I),I=1,60) / 1.344E-04, 1.272E-02, 2.511E-02,
     &                               3.731E-02, 4.932E-02, 6.116E-02,
     &                               7.284E-02, 8.434E-02, 9.568E-02,
     &                               1.069E-01, 1.179E-01, 1.294E-01,
     &                               1.418E-01, 1.540E-01, 1.660E-01,
     &                               1.779E-01, 1.895E-01, 2.009E-01,
     &                               2.121E-01, 2.232E-01, 2.341E-01,
     &                               2.448E-01, 2.563E-01, 2.685E-01,
     &                               2.801E-01, 2.909E-01, 3.011E-01,
     &                               3.108E-01, 3.200E-01, 3.289E-01,
     &                               3.375E-01, 3.459E-01, 3.541E-01,
     &                               3.621E-01, 3.699E-01, 3.780E-01,
     &                               3.870E-01, 3.958E-01, 4.042E-01,
     &                               4.123E-01, 4.199E-01, 4.273E-01,
     &                               4.343E-01, 4.412E-01, 4.477E-01,
     &                               4.541E-01, 4.604E-01, 4.665E-01,
     &                               4.724E-01, 4.782E-01, 4.838E-01,
     &                               4.892E-01, 4.945E-01, 4.997E-01,
     &                               5.053E-01, 5.107E-01, 5.161E-01,
     &                               5.213E-01, 5.265E-01, 5.316E-01/
      DATA (UV_B_EST(13,I),I=1,60) / 1.374E-04, 1.299E-02, 2.563E-02,
     &                               3.808E-02, 5.034E-02, 6.241E-02,
     &                               7.431E-02, 8.603E-02, 9.759E-02,
     &                               1.090E-01, 1.202E-01, 1.322E-01,
     &                               1.448E-01, 1.572E-01, 1.693E-01,
     &                               1.813E-01, 1.930E-01, 2.046E-01,
     &                               2.159E-01, 2.271E-01, 2.381E-01,
     &                               2.489E-01, 2.613E-01, 2.736E-01,
     &                               2.857E-01, 2.970E-01, 3.072E-01,
     &                               3.170E-01, 3.264E-01, 3.354E-01,
     &                               3.441E-01, 3.526E-01, 3.608E-01,
     &                               3.689E-01, 3.770E-01, 3.857E-01,
     &                               3.943E-01, 4.031E-01, 4.115E-01,
     &                               4.196E-01, 4.273E-01, 4.346E-01,
     &                               4.416E-01, 4.484E-01, 4.550E-01,
     &                               4.613E-01, 4.675E-01, 4.735E-01,
     &                               4.794E-01, 4.851E-01, 4.907E-01,
     &                               4.961E-01, 5.016E-01, 5.073E-01,
     &                               5.128E-01, 5.182E-01, 5.235E-01,
     &                               5.287E-01, 5.338E-01, 5.388E-01/
      DATA (UV_B_EST(14,I),I=1,60) / 1.329E-04, 1.326E-02, 2.618E-02,
     &                               3.889E-02, 5.140E-02, 6.372E-02,
     &                               7.584E-02, 8.779E-02, 9.956E-02,
     &                               1.112E-01, 1.226E-01, 1.352E-01,
     &                               1.479E-01, 1.604E-01, 1.727E-01,
     &                               1.848E-01, 1.967E-01, 2.084E-01,
     &                               2.199E-01, 2.311E-01, 2.422E-01,
     &                               2.537E-01, 2.664E-01, 2.788E-01,
     &                               2.910E-01, 3.028E-01, 3.137E-01,
     &                               3.236E-01, 3.331E-01, 3.422E-01,
     &                               3.510E-01, 3.595E-01, 3.678E-01,
     &                               3.761E-01, 3.851E-01, 3.937E-01,
     &                               4.021E-01, 4.106E-01, 4.190E-01,
     &                               4.271E-01, 4.349E-01, 4.421E-01,
     &                               4.491E-01, 4.559E-01, 4.624E-01,
     &                               4.687E-01, 4.748E-01, 4.808E-01,
     &                               4.866E-01, 4.923E-01, 4.978E-01,
     &                               5.036E-01, 5.093E-01, 5.150E-01,
     &                               5.205E-01, 5.258E-01, 5.311E-01,
     &                               5.362E-01, 5.412E-01, 5.462E-01/
      DATA (UV_B_EST(15,I),I=1,60) / 1.280E-04, 1.352E-02, 2.674E-02,
     &                               3.973E-02, 5.250E-02, 6.507E-02,
     &                               7.744E-02, 8.962E-02, 1.016E-01,
     &                               1.134E-01, 1.251E-01, 1.382E-01,
     &                               1.510E-01, 1.637E-01, 1.761E-01,
     &                               1.884E-01, 2.004E-01, 2.122E-01,
     &                               2.239E-01, 2.353E-01, 2.465E-01,
     &                               2.589E-01, 2.717E-01, 2.842E-01,
     &                               2.964E-01, 3.083E-01, 3.200E-01,
     &                               3.306E-01, 3.401E-01, 3.493E-01,
     &                               3.581E-01, 3.667E-01, 3.751E-01,
     &                               3.843E-01, 3.933E-01, 4.020E-01,
     &                               4.103E-01, 4.184E-01, 4.266E-01,
     &                               4.347E-01, 4.425E-01, 4.499E-01,
     &                               4.568E-01, 4.635E-01, 4.700E-01,
     &                               4.763E-01, 4.823E-01, 4.882E-01,
     &                               4.940E-01, 4.996E-01, 5.056E-01,
     &                               5.115E-01, 5.172E-01, 5.228E-01,
     &                               5.282E-01, 5.335E-01, 5.387E-01,
     &                               5.437E-01, 5.487E-01, 5.536E-01/
      DATA (UV_B_EST(16,I),I=1,60) / 1.214E-04, 1.381E-02, 2.731E-02,
     &                               4.057E-02, 5.360E-02, 6.646E-02,
     &                               7.910E-02, 9.153E-02, 1.038E-01,
     &                               1.158E-01, 1.280E-01, 1.413E-01,
     &                               1.543E-01, 1.671E-01, 1.797E-01,
     &                               1.921E-01, 2.042E-01, 2.162E-01,
     &                               2.280E-01, 2.395E-01, 2.511E-01,
     &                               2.642E-01, 2.770E-01, 2.896E-01,
     &                               3.019E-01, 3.139E-01, 3.257E-01,
     &                               3.371E-01, 3.474E-01, 3.566E-01,
     &                               3.655E-01, 3.741E-01, 3.835E-01,
     &                               3.928E-01, 4.017E-01, 4.104E-01,
     &                               4.188E-01, 4.268E-01, 4.346E-01,
     &                               4.426E-01, 4.503E-01, 4.577E-01,
     &                               4.647E-01, 4.714E-01, 4.778E-01,
     &                               4.840E-01, 4.900E-01, 4.959E-01,
     &                               5.017E-01, 5.078E-01, 5.137E-01,
     &                               5.195E-01, 5.252E-01, 5.307E-01,
     &                               5.360E-01, 5.413E-01, 5.464E-01,
     &                               5.514E-01, 5.563E-01, 5.611E-01/
      DATA (UV_B_EST(17,I),I=1,60) / 1.146E-04, 1.407E-02, 2.782E-02,
     &                               4.138E-02, 5.474E-02, 6.792E-02,
     &                               8.086E-02, 9.353E-02, 1.060E-01,
     &                               1.183E-01, 1.310E-01, 1.445E-01,
     &                               1.576E-01, 1.706E-01, 1.833E-01,
     &                               1.959E-01, 2.082E-01, 2.203E-01,
     &                               2.322E-01, 2.439E-01, 2.563E-01,
     &                               2.696E-01, 2.825E-01, 2.952E-01,
     &                               3.075E-01, 3.196E-01, 3.314E-01,
     &                               3.430E-01, 3.543E-01, 3.642E-01,
     &                               3.732E-01, 3.829E-01, 3.923E-01,
     &                               4.014E-01, 4.103E-01, 4.190E-01,
     &                               4.274E-01, 4.354E-01, 4.431E-01,
     &                               4.507E-01, 4.584E-01, 4.658E-01,
     &                               4.728E-01, 4.794E-01, 4.858E-01,
     &                               4.919E-01, 4.979E-01, 5.040E-01,
     &                               5.101E-01, 5.160E-01, 5.219E-01,
     &                               5.276E-01, 5.332E-01, 5.387E-01,
     &                               5.440E-01, 5.492E-01, 5.542E-01,
     &                               5.591E-01, 5.640E-01, 5.687E-01/
      DATA (UV_B_EST(18,I),I=1,60) / 1.075E-04, 1.433E-02, 2.835E-02,
     &                               4.216E-02, 5.577E-02, 6.919E-02,
     &                               8.241E-02, 9.545E-02, 1.084E-01,
     &                               1.208E-01, 1.342E-01, 1.477E-01,
     &                               1.611E-01, 1.742E-01, 1.871E-01,
     &                               1.998E-01, 2.122E-01, 2.245E-01,
     &                               2.365E-01, 2.483E-01, 2.616E-01,
     &                               2.750E-01, 2.881E-01, 3.009E-01,
     &                               3.133E-01, 3.255E-01, 3.374E-01,
     &                               3.490E-01, 3.603E-01, 3.714E-01,
     &                               3.821E-01, 3.920E-01, 4.014E-01,
     &                               4.104E-01, 4.192E-01, 4.278E-01,
     &                               4.361E-01, 4.442E-01, 4.519E-01,
     &                               4.594E-01, 4.667E-01, 4.740E-01,
     &                               4.810E-01, 4.877E-01, 4.940E-01,
     &                               5.001E-01, 5.064E-01, 5.126E-01,
     &                               5.186E-01, 5.244E-01, 5.302E-01,
     &                               5.359E-01, 5.414E-01, 5.468E-01,
     &                               5.520E-01, 5.571E-01, 5.621E-01,
     &                               5.670E-01, 5.718E-01, 5.764E-01/
      DATA (UV_B_EST(19,I),I=1,60) / 1.001E-04, 1.461E-02, 2.890E-02,
     &                               4.298E-02, 5.685E-02, 7.051E-02,
     &                               8.398E-02, 9.725E-02, 1.103E-01,
     &                               1.232E-01, 1.373E-01, 1.511E-01,
     &                               1.646E-01, 1.779E-01, 1.910E-01,
     &                               2.038E-01, 2.164E-01, 2.288E-01,
     &                               2.410E-01, 2.534E-01, 2.670E-01,
     &                               2.805E-01, 2.937E-01, 3.067E-01,
     &                               3.192E-01, 3.315E-01, 3.434E-01,
     &                               3.551E-01, 3.665E-01, 3.784E-01,
     &                               3.917E-01, 4.015E-01, 4.108E-01,
     &                               4.198E-01, 4.283E-01, 4.368E-01,
     &                               4.451E-01, 4.532E-01, 4.609E-01,
     &                               4.683E-01, 4.755E-01, 4.826E-01,
     &                               4.895E-01, 4.961E-01, 5.025E-01,
     &                               5.089E-01, 5.152E-01, 5.213E-01,
     &                               5.272E-01, 5.330E-01, 5.386E-01,
     &                               5.442E-01, 5.497E-01, 5.550E-01,
     &                               5.601E-01, 5.652E-01, 5.701E-01,
     &                               5.749E-01, 5.796E-01, 5.842E-01/
      DATA (UV_B_EST(20,I),I=1,60) / 9.247E-05, 1.490E-02, 2.948E-02,
     &                               4.383E-02, 5.797E-02, 7.189E-02,
     &                               8.561E-02, 9.912E-02, 1.124E-01,
     &                               1.256E-01, 1.401E-01, 1.543E-01,
     &                               1.683E-01, 1.818E-01, 1.950E-01,
     &                               2.080E-01, 2.207E-01, 2.333E-01,
     &                               2.456E-01, 2.588E-01, 2.726E-01,
     &                               2.861E-01, 2.994E-01, 3.124E-01,
     &                               3.253E-01, 3.376E-01, 3.497E-01,
     &                               3.614E-01, 3.729E-01, 3.864E-01,
     &                               4.002E-01, 4.113E-01, 4.205E-01,
     &                               4.294E-01, 4.379E-01, 4.461E-01,
     &                               4.543E-01, 4.623E-01, 4.700E-01,
     &                               4.774E-01, 4.846E-01, 4.915E-01,
     &                               4.983E-01, 5.050E-01, 5.115E-01,
     &                               5.179E-01, 5.241E-01, 5.301E-01,
     &                               5.359E-01, 5.416E-01, 5.472E-01,
     &                               5.527E-01, 5.581E-01, 5.633E-01,
     &                               5.684E-01, 5.734E-01, 5.782E-01,
     &                               5.830E-01, 5.876E-01, 5.922E-01/
      DATA (UV_B_EST(21,I),I=1,60) / 8.448E-05, 1.520E-02, 3.007E-02,
     &                               4.471E-02, 5.913E-02, 7.332E-02,
     &                               8.730E-02, 1.011E-01, 1.146E-01,
     &                               1.283E-01, 1.429E-01, 1.573E-01,
     &                               1.715E-01, 1.855E-01, 1.991E-01,
     &                               2.122E-01, 2.252E-01, 2.378E-01,
     &                               2.504E-01, 2.644E-01, 2.782E-01,
     &                               2.918E-01, 3.051E-01, 3.183E-01,
     &                               3.312E-01, 3.439E-01, 3.560E-01,
     &                               3.678E-01, 3.805E-01, 3.945E-01,
     &                               4.083E-01, 4.214E-01, 4.305E-01,
     &                               4.392E-01, 4.477E-01, 4.559E-01,
     &                               4.638E-01, 4.717E-01, 4.794E-01,
     &                               4.868E-01, 4.939E-01, 5.008E-01,
     &                               5.076E-01, 5.142E-01, 5.207E-01,
     &                               5.270E-01, 5.331E-01, 5.390E-01,
     &                               5.448E-01, 5.504E-01, 5.559E-01,
     &                               5.613E-01, 5.665E-01, 5.717E-01,
     &                               5.767E-01, 5.817E-01, 5.864E-01,
     &                               5.911E-01, 5.957E-01, 6.002E-01/
      DATA (UV_B_EST(22,I),I=1,60) / 7.615E-05, 1.551E-02, 3.069E-02,
     &                               4.564E-02, 6.034E-02, 7.482E-02,
     &                               8.907E-02, 1.031E-01, 1.169E-01,
     &                               1.311E-01, 1.459E-01, 1.605E-01,
     &                               1.748E-01, 1.889E-01, 2.028E-01,
     &                               2.165E-01, 2.297E-01, 2.426E-01,
     &                               2.559E-01, 2.700E-01, 2.839E-01,
     &                               2.976E-01, 3.110E-01, 3.242E-01,
     &                               3.372E-01, 3.499E-01, 3.625E-01,
     &                               3.744E-01, 3.887E-01, 4.028E-01,
     &                               4.165E-01, 4.299E-01, 4.407E-01,
     &                               4.493E-01, 4.577E-01, 4.658E-01,
     &                               4.737E-01, 4.814E-01, 4.890E-01,
     &                               4.963E-01, 5.034E-01, 5.103E-01,
     &                               5.171E-01, 5.236E-01, 5.300E-01,
     &                               5.362E-01, 5.422E-01, 5.480E-01,
     &                               5.537E-01, 5.593E-01, 5.647E-01,
     &                               5.700E-01, 5.752E-01, 5.802E-01,
     &                               5.852E-01, 5.900E-01, 5.948E-01,
     &                               5.994E-01, 6.039E-01, 6.083E-01/
      DATA (UV_B_EST(23,I),I=1,60) / 7.078E-05, 1.573E-02, 3.131E-02,
     &                               4.660E-02, 6.160E-02, 7.637E-02,
     &                               9.090E-02, 1.052E-01, 1.193E-01,
     &                               1.339E-01, 1.489E-01, 1.637E-01,
     &                               1.782E-01, 1.925E-01, 2.066E-01,
     &                               2.204E-01, 2.340E-01, 2.474E-01,
     &                               2.616E-01, 2.758E-01, 2.897E-01,
     &                               3.035E-01, 3.170E-01, 3.303E-01,
     &                               3.433E-01, 3.561E-01, 3.687E-01,
     &                               3.826E-01, 3.971E-01, 4.111E-01,
     &                               4.248E-01, 4.382E-01, 4.511E-01,
     &                               4.596E-01, 4.679E-01, 4.760E-01,
     &                               4.838E-01, 4.914E-01, 4.989E-01,
     &                               5.062E-01, 5.132E-01, 5.201E-01,
     &                               5.267E-01, 5.331E-01, 5.394E-01,
     &                               5.455E-01, 5.514E-01, 5.572E-01,
     &                               5.628E-01, 5.683E-01, 5.736E-01,
     &                               5.788E-01, 5.839E-01, 5.889E-01,
     &                               5.938E-01, 5.985E-01, 6.032E-01,
     &                               6.077E-01, 6.121E-01, 6.165E-01/
      DATA (UV_B_EST(24,I),I=1,60) / 6.614E-05, 1.582E-02, 3.150E-02,
     &                               4.712E-02, 6.267E-02, 7.799E-02,
     &                               9.281E-02, 1.074E-01, 1.218E-01,
     &                               1.369E-01, 1.521E-01, 1.670E-01,
     &                               1.817E-01, 1.962E-01, 2.104E-01,
     &                               2.244E-01, 2.382E-01, 2.520E-01,
     &                               2.671E-01, 2.816E-01, 2.957E-01,
     &                               3.095E-01, 3.231E-01, 3.364E-01,
     &                               3.496E-01, 3.625E-01, 3.752E-01,
     &                               3.906E-01, 4.055E-01, 4.196E-01,
     &                               4.333E-01, 4.467E-01, 4.597E-01,
     &                               4.701E-01, 4.783E-01, 4.863E-01,
     &                               4.941E-01, 5.017E-01, 5.091E-01,
     &                               5.163E-01, 5.232E-01, 5.299E-01,
     &                               5.364E-01, 5.428E-01, 5.490E-01,
     &                               5.550E-01, 5.608E-01, 5.665E-01,
     &                               5.720E-01, 5.774E-01, 5.826E-01,
     &                               5.877E-01, 5.927E-01, 5.976E-01,
     &                               6.024E-01, 6.071E-01, 6.117E-01,
     &                               6.162E-01, 6.205E-01, 6.248E-01/
      DATA (UV_B_EST(25,I),I=1,60) / 6.153E-05, 1.590E-02, 3.168E-02,
     &                               4.739E-02, 6.303E-02, 7.860E-02,
     &                               9.410E-02, 1.095E-01, 1.243E-01,
     &                               1.399E-01, 1.553E-01, 1.704E-01,
     &                               1.853E-01, 2.000E-01, 2.144E-01,
     &                               2.286E-01, 2.426E-01, 2.570E-01,
     &                               2.722E-01, 2.872E-01, 3.017E-01,
     &                               3.156E-01, 3.293E-01, 3.427E-01,
     &                               3.560E-01, 3.689E-01, 3.832E-01,
     &                               3.987E-01, 4.139E-01, 4.282E-01,
     &                               4.419E-01, 4.552E-01, 4.682E-01,
     &                               4.808E-01, 4.889E-01, 4.968E-01,
     &                               5.046E-01, 5.122E-01, 5.195E-01,
     &                               5.265E-01, 5.333E-01, 5.399E-01,
     &                               5.463E-01, 5.526E-01, 5.586E-01,
     &                               5.646E-01, 5.703E-01, 5.759E-01,
     &                               5.813E-01, 5.866E-01, 5.918E-01,
     &                               5.968E-01, 6.017E-01, 6.065E-01,
     &                               6.112E-01, 6.158E-01, 6.203E-01,
     &                               6.247E-01, 6.291E-01, 6.334E-01/
      DATA (UV_B_EST(26,I),I=1,60) / 5.680E-05, 1.598E-02, 3.185E-02,
     &                               4.766E-02, 6.340E-02, 7.907E-02,
     &                               9.467E-02, 1.102E-01, 1.257E-01,
     &                               1.425E-01, 1.586E-01, 1.740E-01,
     &                               1.891E-01, 2.039E-01, 2.185E-01,
     &                               2.329E-01, 2.470E-01, 2.621E-01,
     &                               2.774E-01, 2.925E-01, 3.073E-01,
     &                               3.219E-01, 3.357E-01, 3.492E-01,
     &                               3.625E-01, 3.756E-01, 3.913E-01,
     &                               4.068E-01, 4.220E-01, 4.369E-01,
     &                               4.506E-01, 4.639E-01, 4.769E-01,
     &                               4.895E-01, 4.997E-01, 5.077E-01,
     &                               5.154E-01, 5.228E-01, 5.300E-01,
     &                               5.369E-01, 5.436E-01, 5.501E-01,
     &                               5.564E-01, 5.625E-01, 5.685E-01,
     &                               5.743E-01, 5.799E-01, 5.854E-01,
     &                               5.908E-01, 5.960E-01, 6.010E-01,
     &                               6.060E-01, 6.108E-01, 6.155E-01,
     &                               6.201E-01, 6.246E-01, 6.292E-01,
     &                               6.336E-01, 6.379E-01, 6.421E-01/
      DATA (UV_B_EST(27,I),I=1,60) / 5.143E-05, 1.607E-02, 3.203E-02,
     &                               4.793E-02, 6.377E-02, 7.955E-02,
     &                               9.525E-02, 1.109E-01, 1.266E-01,
     &                               1.434E-01, 1.602E-01, 1.771E-01,
     &                               1.929E-01, 2.080E-01, 2.228E-01,
     &                               2.373E-01, 2.518E-01, 2.673E-01,
     &                               2.827E-01, 2.978E-01, 3.128E-01,
     &                               3.275E-01, 3.421E-01, 3.558E-01,
     &                               3.691E-01, 3.838E-01, 3.995E-01,
     &                               4.150E-01, 4.303E-01, 4.453E-01,
     &                               4.594E-01, 4.727E-01, 4.857E-01,
     &                               4.983E-01, 5.109E-01, 5.188E-01,
     &                               5.264E-01, 5.336E-01, 5.407E-01,
     &                               5.475E-01, 5.541E-01, 5.604E-01,
     &                               5.666E-01, 5.726E-01, 5.785E-01,
     &                               5.841E-01, 5.897E-01, 5.951E-01,
     &                               6.003E-01, 6.055E-01, 6.104E-01,
     &                               6.153E-01, 6.200E-01, 6.246E-01,
     &                               6.292E-01, 6.338E-01, 6.382E-01,
     &                               6.425E-01, 6.467E-01, 6.509E-01/
      DATA (UV_B_EST(28,I),I=1,60) / 4.712E-05, 1.616E-02, 3.221E-02,
     &                               4.821E-02, 6.415E-02, 8.003E-02,
     &                               9.585E-02, 1.116E-01, 1.274E-01,
     &                               1.442E-01, 1.611E-01, 1.779E-01,
     &                               1.948E-01, 2.117E-01, 2.271E-01,
     &                               2.419E-01, 2.570E-01, 2.726E-01,
     &                               2.881E-01, 3.033E-01, 3.184E-01,
     &                               3.332E-01, 3.478E-01, 3.622E-01,
     &                               3.761E-01, 3.920E-01, 4.078E-01,
     &                               4.233E-01, 4.385E-01, 4.536E-01,
     &                               4.684E-01, 4.817E-01, 4.946E-01,
     &                               5.076E-01, 5.205E-01, 5.301E-01,
     &                               5.375E-01, 5.447E-01, 5.516E-01,
     &                               5.582E-01, 5.647E-01, 5.709E-01,
     &                               5.770E-01, 5.828E-01, 5.886E-01,
     &                               5.941E-01, 5.996E-01, 6.050E-01,
     &                               6.100E-01, 6.151E-01, 6.199E-01,
     &                               6.247E-01, 6.294E-01, 6.341E-01,
     &                               6.386E-01, 6.430E-01, 6.473E-01,
     &                               6.516E-01, 6.557E-01, 6.597E-01/
      DATA (UV_B_EST(29,I),I=1,60) / 4.223E-05, 1.624E-02, 3.239E-02,
     &                               4.849E-02, 6.454E-02, 8.052E-02,
     &                               9.644E-02, 1.123E-01, 1.283E-01,
     &                               1.451E-01, 1.619E-01, 1.787E-01,
     &                               1.956E-01, 2.125E-01, 2.293E-01,
     &                               2.463E-01, 2.623E-01, 2.780E-01,
     &                               2.935E-01, 3.089E-01, 3.240E-01,
     &                               3.390E-01, 3.537E-01, 3.682E-01,
     &                               3.841E-01, 4.004E-01, 4.161E-01,
     &                               4.316E-01, 4.469E-01, 4.620E-01,
     &                               4.768E-01, 4.908E-01, 5.039E-01,
     &                               5.171E-01, 5.300E-01, 5.416E-01,
     &                               5.488E-01, 5.558E-01, 5.626E-01,
     &                               5.691E-01, 5.755E-01, 5.816E-01,
     &                               5.875E-01, 5.932E-01, 5.988E-01,
     &                               6.045E-01, 6.100E-01, 6.151E-01,
     &                               6.200E-01, 6.248E-01, 6.297E-01,
     &                               6.344E-01, 6.391E-01, 6.436E-01,
     &                               6.480E-01, 6.523E-01, 6.566E-01,
     &                               6.607E-01, 6.647E-01, 6.687E-01/
      DATA (UV_B_EST(30,I),I=1,60) / 6.510E-05, 1.636E-02, 3.255E-02,
     &                               4.877E-02, 6.492E-02, 8.102E-02,
     &                               9.705E-02, 1.130E-01, 1.291E-01,
     &                               1.459E-01, 1.627E-01, 1.795E-01,
     &                               1.964E-01, 2.132E-01, 2.301E-01,
     &                               2.470E-01, 2.655E-01, 2.835E-01,
     &                               2.991E-01, 3.145E-01, 3.298E-01,
     &                               3.448E-01, 3.597E-01, 3.743E-01,
     &                               3.914E-01, 4.088E-01, 4.245E-01,
     &                               4.401E-01, 4.554E-01, 4.704E-01,
     &                               4.852E-01, 4.998E-01, 5.135E-01,
     &                               5.267E-01, 5.395E-01, 5.520E-01,
     &                               5.603E-01, 5.672E-01, 5.738E-01,
     &                               5.802E-01, 5.864E-01, 5.924E-01,
     &                               5.982E-01, 6.040E-01, 6.097E-01,
     &                               6.153E-01, 6.206E-01, 6.255E-01,
     &                               6.302E-01, 6.348E-01, 6.396E-01,
     &                               6.443E-01, 6.488E-01, 6.532E-01,
     &                               6.575E-01, 6.617E-01, 6.659E-01,
     &                               6.699E-01, 6.739E-01, 6.777E-01/
      DATA (UV_B_EST(31,I),I=1,60) / 1.543E-04, 1.692E-02, 3.345E-02,
     &                               4.975E-02, 6.582E-02, 8.168E-02,
     &                               9.768E-02, 1.138E-01, 1.300E-01,
     &                               1.468E-01, 1.635E-01, 1.803E-01,
     &                               1.971E-01, 2.140E-01, 2.308E-01,
     &                               2.477E-01, 2.663E-01, 2.852E-01,
     &                               3.041E-01, 3.203E-01, 3.357E-01,
     &                               3.508E-01, 3.658E-01, 3.814E-01,
     &                               3.987E-01, 4.160E-01, 4.330E-01,
     &                               4.486E-01, 4.639E-01, 4.790E-01,
     &                               4.938E-01, 5.088E-01, 5.233E-01,
     &                               5.364E-01, 5.492E-01, 5.616E-01,
     &                               5.720E-01, 5.787E-01, 5.852E-01,
     &                               5.915E-01, 5.975E-01, 6.035E-01,
     &                               6.094E-01, 6.151E-01, 6.207E-01,
     &                               6.261E-01, 6.314E-01, 6.361E-01,
     &                               6.407E-01, 6.451E-01, 6.497E-01,
     &                               6.542E-01, 6.587E-01, 6.630E-01,
     &                               6.672E-01, 6.713E-01, 6.753E-01,
     &                               6.793E-01, 6.831E-01, 6.869E-01/
      DATA (UV_B_EST(32,I),I=1,60) / 2.489E-04, 1.751E-02, 3.452E-02,
     &                               5.129E-02, 6.782E-02, 8.412E-02,
     &                               1.002E-01, 1.160E-01, 1.318E-01,
     &                               1.476E-01, 1.644E-01, 1.811E-01,
     &                               1.979E-01, 2.147E-01, 2.315E-01,
     &                               2.484E-01, 2.671E-01, 2.859E-01,
     &                               3.048E-01, 3.238E-01, 3.417E-01,
     &                               3.569E-01, 3.720E-01, 3.887E-01,
     &                               4.059E-01, 4.232E-01, 4.406E-01,
     &                               4.572E-01, 4.726E-01, 4.877E-01,
     &                               5.026E-01, 5.178E-01, 5.328E-01,
     &                               5.462E-01, 5.590E-01, 5.714E-01,
     &                               5.834E-01, 5.905E-01, 5.968E-01,
     &                               6.029E-01, 6.090E-01, 6.150E-01,
     &                               6.208E-01, 6.264E-01, 6.318E-01,
     &                               6.371E-01, 6.422E-01, 6.469E-01,
     &                               6.513E-01, 6.556E-01, 6.598E-01,
     &                               6.643E-01, 6.686E-01, 6.728E-01,
     &                               6.769E-01, 6.809E-01, 6.848E-01,
     &                               6.887E-01, 6.924E-01, 6.961E-01/
      DATA (UV_B_EST(33,I),I=1,60) / 3.493E-04, 1.814E-02, 3.566E-02,
     &                               5.293E-02, 6.994E-02, 8.670E-02,
     &                               1.032E-01, 1.195E-01, 1.359E-01,
     &                               1.521E-01, 1.680E-01, 1.837E-01,
     &                               1.991E-01, 2.155E-01, 2.323E-01,
     &                               2.491E-01, 2.678E-01, 2.867E-01,
     &                               3.056E-01, 3.245E-01, 3.434E-01,
     &                               3.623E-01, 3.788E-01, 3.959E-01,
     &                               4.131E-01, 4.304E-01, 4.478E-01,
     &                               4.653E-01, 4.813E-01, 4.964E-01,
     &                               5.117E-01, 5.269E-01, 5.419E-01,
     &                               5.562E-01, 5.689E-01, 5.812E-01,
     &                               5.932E-01, 6.024E-01, 6.086E-01,
     &                               6.147E-01, 6.207E-01, 6.266E-01,
     &                               6.323E-01, 6.377E-01, 6.430E-01,
     &                               6.481E-01, 6.531E-01, 6.578E-01,
     &                               6.620E-01, 6.662E-01, 6.702E-01,
     &                               6.745E-01, 6.787E-01, 6.828E-01,
     &                               6.868E-01, 6.907E-01, 6.945E-01,
     &                               6.982E-01, 7.019E-01, 7.054E-01/
      DATA (UV_B_EST(34,I),I=1,60) / 3.968E-04, 1.876E-02, 3.686E-02,
     &                               5.467E-02, 7.219E-02, 8.945E-02,
     &                               1.065E-01, 1.232E-01, 1.401E-01,
     &                               1.568E-01, 1.731E-01, 1.892E-01,
     &                               2.049E-01, 2.204E-01, 2.356E-01,
     &                               2.506E-01, 2.686E-01, 2.875E-01,
     &                               3.063E-01, 3.252E-01, 3.441E-01,
     &                               3.630E-01, 3.835E-01, 4.031E-01,
     &                               4.203E-01, 4.375E-01, 4.549E-01,
     &                               4.724E-01, 4.900E-01, 5.055E-01,
     &                               5.209E-01, 5.361E-01, 5.511E-01,
     &                               5.659E-01, 5.790E-01, 5.912E-01,
     &                               6.031E-01, 6.145E-01, 6.206E-01,
     &                               6.266E-01, 6.325E-01, 6.382E-01,
     &                               6.438E-01, 6.492E-01, 6.543E-01,
     &                               6.593E-01, 6.641E-01, 6.687E-01,
     &                               6.729E-01, 6.769E-01, 6.808E-01,
     &                               6.847E-01, 6.888E-01, 6.928E-01,
     &                               6.967E-01, 7.005E-01, 7.042E-01,
     &                               7.078E-01, 7.114E-01, 7.149E-01/
      DATA (UV_B_EST(35,I),I=1,60) / 4.093E-04, 1.935E-02, 3.799E-02,
     &                               5.635E-02, 7.444E-02, 9.226E-02,
     &                               1.098E-01, 1.272E-01, 1.446E-01,
     &                               1.617E-01, 1.785E-01, 1.950E-01,
     &                               2.112E-01, 2.270E-01, 2.426E-01,
     &                               2.583E-01, 2.741E-01, 2.895E-01,
     &                               3.071E-01, 3.260E-01, 3.448E-01,
     &                               3.637E-01, 3.844E-01, 4.073E-01,
     &                               4.274E-01, 4.446E-01, 4.620E-01,
     &                               4.794E-01, 4.970E-01, 5.147E-01,
     &                               5.301E-01, 5.453E-01, 5.603E-01,
     &                               5.751E-01, 5.892E-01, 6.014E-01,
     &                               6.132E-01, 6.247E-01, 6.327E-01,
     &                               6.385E-01, 6.443E-01, 6.499E-01,
     &                               6.555E-01, 6.607E-01, 6.657E-01,
     &                               6.705E-01, 6.752E-01, 6.797E-01,
     &                               6.838E-01, 6.877E-01, 6.914E-01,
     &                               6.951E-01, 6.991E-01, 7.030E-01,
     &                               7.068E-01, 7.105E-01, 7.141E-01,
     &                               7.176E-01, 7.210E-01, 7.244E-01/
      DATA (UV_B_EST(36,I),I=1,60) / 4.227E-04, 1.997E-02, 3.920E-02,
     &                               5.812E-02, 7.674E-02, 9.508E-02,
     &                               1.131E-01, 1.311E-01, 1.491E-01,
     &                               1.669E-01, 1.843E-01, 2.012E-01,
     &                               2.178E-01, 2.340E-01, 2.500E-01,
     &                               2.665E-01, 2.825E-01, 2.982E-01,
     &                               3.135E-01, 3.284E-01, 3.456E-01,
     &                               3.644E-01, 3.853E-01, 4.084E-01,
     &                               4.312E-01, 4.517E-01, 4.690E-01,
     &                               4.864E-01, 5.042E-01, 5.227E-01,
     &                               5.394E-01, 5.546E-01, 5.697E-01,
     &                               5.845E-01, 5.990E-01, 6.117E-01,
     &                               6.235E-01, 6.348E-01, 6.448E-01,
     &                               6.506E-01, 6.562E-01, 6.618E-01,
     &                               6.672E-01, 6.723E-01, 6.772E-01,
     &                               6.819E-01, 6.864E-01, 6.907E-01,
     &                               6.948E-01, 6.986E-01, 7.022E-01,
     &                               7.058E-01, 7.095E-01, 7.133E-01,
     &                               7.170E-01, 7.205E-01, 7.240E-01,
     &                               7.274E-01, 7.308E-01, 7.340E-01/
      DATA (UV_B_EST(37,I),I=1,60) / 4.125E-04, 2.064E-02, 4.048E-02,
     &                               6.000E-02, 7.919E-02, 9.807E-02,
     &                               1.167E-01, 1.352E-01, 1.536E-01,
     &                               1.718E-01, 1.896E-01, 2.072E-01,
     &                               2.245E-01, 2.415E-01, 2.583E-01,
     &                               2.750E-01, 2.914E-01, 3.074E-01,
     &                               3.230E-01, 3.382E-01, 3.530E-01,
     &                               3.675E-01, 3.863E-01, 4.096E-01,
     &                               4.326E-01, 4.551E-01, 4.760E-01,
     &                               4.934E-01, 5.117E-01, 5.303E-01,
     &                               5.483E-01, 5.640E-01, 5.791E-01,
     &                               5.938E-01, 6.084E-01, 6.222E-01,
     &                               6.339E-01, 6.451E-01, 6.557E-01,
     &                               6.628E-01, 6.683E-01, 6.737E-01,
     &                               6.790E-01, 6.841E-01, 6.888E-01,
     &                               6.933E-01, 6.977E-01, 7.019E-01,
     &                               7.059E-01, 7.095E-01, 7.131E-01,
     &                               7.166E-01, 7.200E-01, 7.237E-01,
     &                               7.272E-01, 7.307E-01, 7.341E-01,
     &                               7.374E-01, 7.406E-01, 7.437E-01/
      DATA (UV_B_EST(38,I),I=1,60) / 3.673E-04, 2.118E-02, 4.170E-02,
     &                               6.191E-02, 8.179E-02, 1.012E-01,
     &                               1.204E-01, 1.395E-01, 1.584E-01,
     &                               1.769E-01, 1.952E-01, 2.132E-01,
     &                               2.309E-01, 2.483E-01, 2.662E-01,
     &                               2.839E-01, 3.008E-01, 3.171E-01,
     &                               3.329E-01, 3.484E-01, 3.635E-01,
     &                               3.785E-01, 3.940E-01, 4.108E-01,
     &                               4.340E-01, 4.568E-01, 4.791E-01,
     &                               5.003E-01, 5.194E-01, 5.380E-01,
     &                               5.561E-01, 5.735E-01, 5.886E-01,
     &                               6.033E-01, 6.179E-01, 6.323E-01,
     &                               6.446E-01, 6.556E-01, 6.661E-01,
     &                               6.751E-01, 6.805E-01, 6.857E-01,
     &                               6.908E-01, 6.958E-01, 7.004E-01,
     &                               7.048E-01, 7.091E-01, 7.131E-01,
     &                               7.170E-01, 7.206E-01, 7.241E-01,
     &                               7.274E-01, 7.307E-01, 7.342E-01,
     &                               7.376E-01, 7.410E-01, 7.442E-01,
     &                               7.474E-01, 7.505E-01, 7.535E-01/
      DATA (UV_B_EST(39,I),I=1,60) / 2.942E-04, 2.170E-02, 4.280E-02,
     &                               6.359E-02, 8.407E-02, 1.042E-01,
     &                               1.241E-01, 1.439E-01, 1.633E-01,
     &                               1.824E-01, 2.011E-01, 2.195E-01,
     &                               2.376E-01, 2.557E-01, 2.740E-01,
     &                               2.920E-01, 3.097E-01, 3.270E-01,
     &                               3.435E-01, 3.592E-01, 3.746E-01,
     &                               3.906E-01, 4.062E-01, 4.214E-01,
     &                               4.362E-01, 4.585E-01, 4.811E-01,
     &                               5.038E-01, 5.273E-01, 5.460E-01,
     &                               5.641E-01, 5.818E-01, 5.981E-01,
     &                               6.129E-01, 6.276E-01, 6.422E-01,
     &                               6.555E-01, 6.664E-01, 6.768E-01,
     &                               6.866E-01, 6.928E-01, 6.979E-01,
     &                               7.029E-01, 7.077E-01, 7.122E-01,
     &                               7.165E-01, 7.206E-01, 7.244E-01,
     &                               7.282E-01, 7.318E-01, 7.351E-01,
     &                               7.384E-01, 7.416E-01, 7.448E-01,
     &                               7.481E-01, 7.514E-01, 7.545E-01,
     &                               7.575E-01, 7.605E-01, 7.634E-01/
      DATA (UV_B_EST(40,I),I=1,60) / 2.168E-04, 2.226E-02, 4.397E-02,
     &                               6.537E-02, 8.644E-02, 1.072E-01,
     &                               1.277E-01, 1.481E-01, 1.682E-01,
     &                               1.879E-01, 2.072E-01, 2.262E-01,
     &                               2.447E-01, 2.634E-01, 2.821E-01,
     &                               3.004E-01, 3.183E-01, 3.359E-01,
     &                               3.532E-01, 3.702E-01, 3.870E-01,
     &                               4.032E-01, 4.190E-01, 4.343E-01,
     &                               4.491E-01, 4.635E-01, 4.830E-01,
     &                               5.066E-01, 5.323E-01, 5.541E-01,
     &                               5.724E-01, 5.901E-01, 6.072E-01,
     &                               6.227E-01, 6.379E-01, 6.524E-01,
     &                               6.661E-01, 6.775E-01, 6.877E-01,
     &                               6.974E-01, 7.052E-01, 7.102E-01,
     &                               7.150E-01, 7.197E-01, 7.241E-01,
     &                               7.282E-01, 7.321E-01, 7.359E-01,
     &                               7.395E-01, 7.430E-01, 7.462E-01,
     &                               7.494E-01, 7.525E-01, 7.555E-01,
     &                               7.587E-01, 7.618E-01, 7.648E-01,
     &                               7.677E-01, 7.706E-01, 7.734E-01/
      DATA (UV_B_EST(41,I),I=1,60) / 1.638E-04, 2.286E-02, 4.521E-02,
     &                               6.726E-02, 8.897E-02, 1.104E-01,
     &                               1.315E-01, 1.525E-01, 1.732E-01,
     &                               1.934E-01, 2.133E-01, 2.329E-01,
     &                               2.521E-01, 2.715E-01, 2.905E-01,
     &                               3.091E-01, 3.274E-01, 3.453E-01,
     &                               3.629E-01, 3.805E-01, 3.987E-01,
     &                               4.163E-01, 4.322E-01, 4.476E-01,
     &                               4.626E-01, 4.771E-01, 4.911E-01,
     &                               5.097E-01, 5.362E-01, 5.606E-01,
     &                               5.808E-01, 5.986E-01, 6.158E-01,
     &                               6.331E-01, 6.484E-01, 6.629E-01,
     &                               6.765E-01, 6.888E-01, 6.989E-01,
     &                               7.084E-01, 7.174E-01, 7.226E-01,
     &                               7.273E-01, 7.317E-01, 7.361E-01,
     &                               7.400E-01, 7.438E-01, 7.474E-01,
     &                               7.509E-01, 7.542E-01, 7.574E-01,
     &                               7.604E-01, 7.634E-01, 7.664E-01,
     &                               7.694E-01, 7.723E-01, 7.752E-01,
     &                               7.781E-01, 7.808E-01, 7.835E-01/
      DATA (UV_B_EST(42,I),I=1,60) / 1.707E-04, 2.381E-02, 4.699E-02,
     &                               6.971E-02, 9.199E-02, 1.138E-01,
     &                               1.356E-01, 1.572E-01, 1.785E-01,
     &                               1.994E-01, 2.198E-01, 2.399E-01,
     &                               2.600E-01, 2.799E-01, 2.993E-01,
     &                               3.183E-01, 3.366E-01, 3.551E-01,
     &                               3.730E-01, 3.915E-01, 4.098E-01,
     &                               4.278E-01, 4.455E-01, 4.616E-01,
     &                               4.766E-01, 4.912E-01, 5.057E-01,
     &                               5.204E-01, 5.403E-01, 5.654E-01,
     &                               5.884E-01, 6.073E-01, 6.246E-01,
     &                               6.430E-01, 6.593E-01, 6.737E-01,
     &                               6.872E-01, 6.999E-01, 7.103E-01,
     &                               7.197E-01, 7.285E-01, 7.352E-01,
     &                               7.397E-01, 7.440E-01, 7.481E-01,
     &                               7.519E-01, 7.556E-01, 7.590E-01,
     &                               7.623E-01, 7.656E-01, 7.686E-01,
     &                               7.716E-01, 7.745E-01, 7.774E-01,
     &                               7.802E-01, 7.830E-01, 7.858E-01,
     &                               7.885E-01, 7.912E-01, 7.937E-01/
      DATA (UV_B_EST(43,I),I=1,60) / 1.781E-04, 2.484E-02, 4.900E-02,
     &                               7.266E-02, 9.585E-02, 1.186E-01,
     &                               1.410E-01, 1.631E-01, 1.846E-01,
     &                               2.059E-01, 2.268E-01, 2.475E-01,
     &                               2.683E-01, 2.887E-01, 3.087E-01,
     &                               3.281E-01, 3.471E-01, 3.656E-01,
     &                               3.840E-01, 4.027E-01, 4.212E-01,
     &                               4.394E-01, 4.572E-01, 4.747E-01,
     &                               4.913E-01, 5.064E-01, 5.216E-01,
     &                               5.362E-01, 5.502E-01, 5.707E-01,
     &                               5.943E-01, 6.160E-01, 6.349E-01,
     &                               6.534E-01, 6.704E-01, 6.849E-01,
     &                               6.983E-01, 7.108E-01, 7.221E-01,
     &                               7.312E-01, 7.399E-01, 7.481E-01,
     &                               7.523E-01, 7.563E-01, 7.603E-01,
     &                               7.640E-01, 7.674E-01, 7.707E-01,
     &                               7.739E-01, 7.770E-01, 7.799E-01,
     &                               7.828E-01, 7.856E-01, 7.884E-01,
     &                               7.911E-01, 7.938E-01, 7.965E-01,
     &                               7.991E-01, 8.016E-01, 8.041E-01/
      DATA (UV_B_EST(44,I),I=1,60) / 1.861E-04, 2.596E-02, 5.118E-02,
     &                               7.588E-02, 1.001E-01, 1.237E-01,
     &                               1.471E-01, 1.700E-01, 1.924E-01,
     &                               2.143E-01, 2.356E-01, 2.567E-01,
     &                               2.774E-01, 2.982E-01, 3.186E-01,
     &                               3.386E-01, 3.581E-01, 3.771E-01,
     &                               3.963E-01, 4.149E-01, 4.329E-01,
     &                               4.512E-01, 4.693E-01, 4.869E-01,
     &                               5.047E-01, 5.229E-01, 5.380E-01,
     &                               5.526E-01, 5.665E-01, 5.799E-01,
     &                               6.007E-01, 6.229E-01, 6.457E-01,
     &                               6.642E-01, 6.812E-01, 6.964E-01,
     &                               7.097E-01, 7.221E-01, 7.337E-01,
     &                               7.431E-01, 7.515E-01, 7.593E-01,
     &                               7.650E-01, 7.689E-01, 7.726E-01,
     &                               7.761E-01, 7.794E-01, 7.825E-01,
     &                               7.856E-01, 7.885E-01, 7.913E-01,
     &                               7.941E-01, 7.968E-01, 7.995E-01,
     &                               8.022E-01, 8.048E-01, 8.073E-01,
     &                               8.098E-01, 8.122E-01, 8.146E-01/
      DATA (UV_B_EST(45,I),I=1,60) / 1.440E-04, 2.723E-02, 5.366E-02,
     &                               7.945E-02, 1.047E-01, 1.294E-01,
     &                               1.538E-01, 1.776E-01, 2.009E-01,
     &                               2.236E-01, 2.458E-01, 2.677E-01,
     &                               2.891E-01, 3.099E-01, 3.301E-01,
     &                               3.498E-01, 3.698E-01, 3.897E-01,
     &                               4.092E-01, 4.282E-01, 4.466E-01,
     &                               4.644E-01, 4.817E-01, 4.996E-01,
     &                               5.188E-01, 5.375E-01, 5.551E-01,
     &                               5.696E-01, 5.834E-01, 5.967E-01,
     &                               6.094E-01, 6.310E-01, 6.552E-01,
     &                               6.755E-01, 6.925E-01, 7.081E-01,
     &                               7.215E-01, 7.338E-01, 7.452E-01,
     &                               7.553E-01, 7.634E-01, 7.709E-01,
     &                               7.779E-01, 7.816E-01, 7.851E-01,
     &                               7.884E-01, 7.915E-01, 7.944E-01,
     &                               7.973E-01, 8.001E-01, 8.030E-01,
     &                               8.058E-01, 8.085E-01, 8.111E-01,
     &                               8.136E-01, 8.160E-01, 8.183E-01,
     &                               8.207E-01, 8.230E-01, 8.252E-01/
      DATA (UV_B_EST(46,I),I=1,60) / 9.000E-05, 2.866E-02, 5.649E-02,
     &                               8.363E-02, 1.101E-01, 1.360E-01,
     &                               1.613E-01, 1.860E-01, 2.101E-01,
     &                               2.338E-01, 2.569E-01, 2.797E-01,
     &                               3.018E-01, 3.233E-01, 3.442E-01,
     &                               3.644E-01, 3.842E-01, 4.035E-01,
     &                               4.230E-01, 4.423E-01, 4.610E-01,
     &                               4.792E-01, 4.967E-01, 5.146E-01,
     &                               5.332E-01, 5.521E-01, 5.704E-01,
     &                               5.873E-01, 6.010E-01, 6.142E-01,
     &                               6.268E-01, 6.396E-01, 6.640E-01,
     &                               6.873E-01, 7.043E-01, 7.198E-01,
     &                               7.338E-01, 7.458E-01, 7.572E-01,
     &                               7.676E-01, 7.757E-01, 7.830E-01,
     &                               7.897E-01, 7.946E-01, 7.978E-01,
     &                               8.009E-01, 8.040E-01, 8.070E-01,
     &                               8.099E-01, 8.127E-01, 8.154E-01,
     &                               8.180E-01, 8.205E-01, 8.229E-01,
     &                               8.251E-01, 8.273E-01, 8.294E-01,
     &                               8.316E-01, 8.338E-01, 8.359E-01/
      DATA (UV_B_EST(47,I),I=1,60) / 3.204E-05, 3.025E-02, 5.966E-02,
     &                               8.828E-02, 1.162E-01, 1.435E-01,
     &                               1.701E-01, 1.959E-01, 2.211E-01,
     &                               2.455E-01, 2.695E-01, 2.928E-01,
     &                               3.157E-01, 3.379E-01, 3.594E-01,
     &                               3.803E-01, 4.008E-01, 4.205E-01,
     &                               4.395E-01, 4.579E-01, 4.764E-01,
     &                               4.949E-01, 5.134E-01, 5.316E-01,
     &                               5.492E-01, 5.671E-01, 5.855E-01,
     &                               6.033E-01, 6.194E-01, 6.325E-01,
     &                               6.453E-01, 6.575E-01, 6.730E-01,
     &                               6.965E-01, 7.167E-01, 7.321E-01,
     &                               7.461E-01, 7.587E-01, 7.698E-01,
     &                               7.798E-01, 7.885E-01, 7.954E-01,
     &                               8.020E-01, 8.082E-01, 8.114E-01,
     &                               8.145E-01, 8.174E-01, 8.202E-01,
     &                               8.228E-01, 8.254E-01, 8.280E-01,
     &                               8.304E-01, 8.326E-01, 8.348E-01,
     &                               8.369E-01, 8.389E-01, 8.408E-01,
     &                               8.428E-01, 8.448E-01, 8.468E-01/
      DATA (UV_B_EST(48,I),I=1,60) / 3.433E-05, 3.198E-02, 6.351E-02,
     &                               9.351E-02, 1.229E-01, 1.518E-01,
     &                               1.798E-01, 2.070E-01, 2.333E-01,
     &                               2.590E-01, 2.841E-01, 3.082E-01,
     &                               3.316E-01, 3.541E-01, 3.760E-01,
     &                               3.978E-01, 4.187E-01, 4.389E-01,
     &                               4.584E-01, 4.771E-01, 4.951E-01,
     &                               5.126E-01, 5.310E-01, 5.493E-01,
     &                               5.670E-01, 5.841E-01, 6.011E-01,
     &                               6.189E-01, 6.365E-01, 6.518E-01,
     &                               6.644E-01, 6.764E-01, 6.879E-01,
     &                               7.060E-01, 7.285E-01, 7.450E-01,
     &                               7.594E-01, 7.724E-01, 7.830E-01,
     &                               7.926E-01, 8.014E-01, 8.084E-01,
     &                               8.150E-01, 8.214E-01, 8.251E-01,
     &                               8.281E-01, 8.310E-01, 8.335E-01,
     &                               8.359E-01, 8.382E-01, 8.406E-01,
     &                               8.429E-01, 8.449E-01, 8.469E-01,
     &                               8.488E-01, 8.506E-01, 8.524E-01,
     &                               8.541E-01, 8.560E-01, 8.578E-01/
      DATA (UV_B_EST(49,I),I=1,60) / 1.161E-04, 3.381E-02, 6.668E-02,
     &                               9.874E-02, 1.300E-01, 1.608E-01,
     &                               1.909E-01, 2.193E-01, 2.470E-01,
     &                               2.741E-01, 3.003E-01, 3.254E-01,
     &                               3.497E-01, 3.731E-01, 3.958E-01,
     &                               4.175E-01, 4.385E-01, 4.590E-01,
     &                               4.789E-01, 4.980E-01, 5.165E-01,
     &                               5.342E-01, 5.512E-01, 5.678E-01,
     &                               5.856E-01, 6.028E-01, 6.193E-01,
     &                               6.359E-01, 6.537E-01, 6.703E-01,
     &                               6.841E-01, 6.959E-01, 7.072E-01,
     &                               7.179E-01, 7.384E-01, 7.592E-01,
     &                               7.734E-01, 7.863E-01, 7.969E-01,
     &                               8.060E-01, 8.144E-01, 8.217E-01,
     &                               8.281E-01, 8.344E-01, 8.389E-01,
     &                               8.418E-01, 8.445E-01, 8.469E-01,
     &                               8.490E-01, 8.512E-01, 8.534E-01,
     &                               8.555E-01, 8.574E-01, 8.591E-01,
     &                               8.608E-01, 8.625E-01, 8.641E-01,
     &                               8.656E-01, 8.672E-01, 8.690E-01/
      DATA (UV_B_EST(50,I),I=1,60) / 1.422E-04, 3.582E-02, 7.056E-02,
     &                               1.044E-01, 1.375E-01, 1.699E-01,
     &                               2.013E-01, 2.318E-01, 2.616E-01,
     &                               2.906E-01, 3.183E-01, 3.446E-01,
     &                               3.699E-01, 3.943E-01, 4.176E-01,
     &                               4.400E-01, 4.614E-01, 4.819E-01,
     &                               5.015E-01, 5.208E-01, 5.396E-01,
     &                               5.574E-01, 5.745E-01, 5.907E-01,
     &                               6.062E-01, 6.225E-01, 6.398E-01,
     &                               6.561E-01, 6.720E-01, 6.885E-01,
     &                               7.038E-01, 7.162E-01, 7.272E-01,
     &                               7.376E-01, 7.488E-01, 7.723E-01,
     &                               7.882E-01, 8.006E-01, 8.114E-01,
     &                               8.201E-01, 8.279E-01, 8.351E-01,
     &                               8.414E-01, 8.475E-01, 8.529E-01,
     &                               8.556E-01, 8.582E-01, 8.604E-01,
     &                               8.623E-01, 8.642E-01, 8.662E-01,
     &                               8.682E-01, 8.699E-01, 8.715E-01,
     &                               8.730E-01, 8.745E-01, 8.759E-01,
     &                               8.773E-01, 8.787E-01, 8.803E-01/
      DATA (UV_B_EST(51,I),I=1,60) / 2.628E-04, 3.777E-02, 7.477E-02,
     &                               1.108E-01, 1.459E-01, 1.801E-01,
     &                               2.133E-01, 2.454E-01, 2.769E-01,
     &                               3.073E-01, 3.366E-01, 3.649E-01,
     &                               3.922E-01, 4.178E-01, 4.419E-01,
     &                               4.649E-01, 4.869E-01, 5.079E-01,
     &                               5.278E-01, 5.467E-01, 5.647E-01,
     &                               5.826E-01, 5.997E-01, 6.160E-01,
     &                               6.315E-01, 6.463E-01, 6.617E-01,
     &                               6.779E-01, 6.930E-01, 7.078E-01,
     &                               7.230E-01, 7.370E-01, 7.480E-01,
     &                               7.584E-01, 7.682E-01, 7.852E-01,
     &                               8.037E-01, 8.156E-01, 8.264E-01,
     &                               8.348E-01, 8.421E-01, 8.488E-01,
     &                               8.550E-01, 8.608E-01, 8.664E-01,
     &                               8.694E-01, 8.719E-01, 8.741E-01,
     &                               8.758E-01, 8.774E-01, 8.792E-01,
     &                               8.809E-01, 8.826E-01, 8.840E-01,
     &                               8.853E-01, 8.866E-01, 8.879E-01,
     &                               8.891E-01, 8.903E-01, 8.917E-01/
      DATA (UV_B_EST(52,I),I=1,60) / 2.669E-04, 4.071E-02, 8.001E-02,
     &                               1.182E-01, 1.555E-01, 1.920E-01,
     &                               2.268E-01, 2.609E-01, 2.940E-01,
     &                               3.260E-01, 3.568E-01, 3.864E-01,
     &                               4.149E-01, 4.421E-01, 4.682E-01,
     &                               4.928E-01, 5.153E-01, 5.366E-01,
     &                               5.567E-01, 5.758E-01, 5.938E-01,
     &                               6.109E-01, 6.272E-01, 6.434E-01,
     &                               6.588E-01, 6.733E-01, 6.870E-01,
     &                               7.011E-01, 7.159E-01, 7.298E-01,
     &                               7.434E-01, 7.574E-01, 7.703E-01,
     &                               7.803E-01, 7.896E-01, 7.994E-01,
     &                               8.186E-01, 8.314E-01, 8.416E-01,
     &                               8.501E-01, 8.569E-01, 8.631E-01,
     &                               8.688E-01, 8.743E-01, 8.797E-01,
     &                               8.834E-01, 8.857E-01, 8.878E-01,
     &                               8.893E-01, 8.907E-01, 8.922E-01,
     &                               8.938E-01, 8.953E-01, 8.965E-01,
     &                               8.977E-01, 8.989E-01, 9.000E-01,
     &                               9.011E-01, 9.021E-01, 9.032E-01/
      DATA (UV_B_EST(53,I),I=1,60) /-1.332E-04, 4.344E-02, 8.579E-02,
     &                               1.269E-01, 1.668E-01, 2.054E-01,
     &                               2.428E-01, 2.789E-01, 3.137E-01,
     &                               3.473E-01, 3.796E-01, 4.108E-01,
     &                               4.405E-01, 4.688E-01, 4.958E-01,
     &                               5.215E-01, 5.457E-01, 5.684E-01,
     &                               5.888E-01, 6.079E-01, 6.260E-01,
     &                               6.428E-01, 6.586E-01, 6.734E-01,
     &                               6.882E-01, 7.023E-01, 7.155E-01,
     &                               7.279E-01, 7.405E-01, 7.540E-01,
     &                               7.666E-01, 7.795E-01, 7.923E-01,
     &                               8.029E-01, 8.118E-01, 8.201E-01,
     &                               8.336E-01, 8.481E-01, 8.576E-01,
     &                               8.661E-01, 8.724E-01, 8.781E-01,
     &                               8.833E-01, 8.883E-01, 8.933E-01,
     &                               8.974E-01, 8.995E-01, 9.016E-01,
     &                               9.029E-01, 9.041E-01, 9.054E-01,
     &                               9.068E-01, 9.081E-01, 9.092E-01,
     &                               9.102E-01, 9.112E-01, 9.122E-01,
     &                               9.131E-01, 9.140E-01, 9.149E-01/
      DATA (UV_B_EST(54,I),I=1,60) / 1.031E-05, 4.656E-02, 9.223E-02,
     &                               1.369E-01, 1.799E-01, 2.214E-01,
     &                               2.615E-01, 3.002E-01, 3.373E-01,
     &                               3.729E-01, 4.068E-01, 4.390E-01,
     &                               4.697E-01, 4.993E-01, 5.272E-01,
     &                               5.535E-01, 5.783E-01, 6.016E-01,
     &                               6.235E-01, 6.436E-01, 6.612E-01,
     &                               6.776E-01, 6.930E-01, 7.073E-01,
     &                               7.208E-01, 7.338E-01, 7.464E-01,
     &                               7.581E-01, 7.688E-01, 7.802E-01,
     &                               7.922E-01, 8.032E-01, 8.150E-01,
     &                               8.262E-01, 8.347E-01, 8.424E-01,
     &                               8.501E-01, 8.647E-01, 8.745E-01,
     &                               8.824E-01, 8.887E-01, 8.938E-01,
     &                               8.984E-01, 9.027E-01, 9.072E-01,
     &                               9.114E-01, 9.135E-01, 9.154E-01,
     &                               9.166E-01, 9.177E-01, 9.187E-01,
     &                               9.198E-01, 9.210E-01, 9.220E-01,
     &                               9.229E-01, 9.237E-01, 9.245E-01,
     &                               9.253E-01, 9.261E-01, 9.268E-01/
      DATA (UV_B_EST(55,I),I=1,60) / 6.499E-04, 5.295E-02, 1.026E-01,
     &                               1.504E-01, 1.967E-01, 2.405E-01,
     &                               2.835E-01, 3.252E-01, 3.650E-01,
     &                               4.030E-01, 4.390E-01, 4.731E-01,
     &                               5.053E-01, 5.353E-01, 5.634E-01,
     &                               5.901E-01, 6.155E-01, 6.390E-01,
     &                               6.608E-01, 6.810E-01, 6.998E-01,
     &                               7.162E-01, 7.310E-01, 7.447E-01,
     &                               7.573E-01, 7.689E-01, 7.798E-01,
     &                               7.905E-01, 8.004E-01, 8.094E-01,
     &                               8.195E-01, 8.297E-01, 8.392E-01,
     &                               8.496E-01, 8.587E-01, 8.656E-01,
     &                               8.721E-01, 8.823E-01, 8.928E-01,
     &                               8.997E-01, 9.057E-01, 9.100E-01,
     &                               9.140E-01, 9.177E-01, 9.215E-01,
     &                               9.252E-01, 9.276E-01, 9.293E-01,
     &                               9.305E-01, 9.313E-01, 9.321E-01,
     &                               9.330E-01, 9.340E-01, 9.349E-01,
     &                               9.356E-01, 9.362E-01, 9.369E-01,
     &                               9.375E-01, 9.382E-01, 9.388E-01/
      DATA (UV_B_EST(56,I),I=1,60) /-2.680E-04, 5.948E-02, 1.173E-01,
     &                               1.713E-01, 2.223E-01, 2.706E-01,
     &                               3.165E-01, 3.598E-01, 4.005E-01,
     &                               4.393E-01, 4.775E-01, 5.137E-01,
     &                               5.474E-01, 5.787E-01, 6.078E-01,
     &                               6.345E-01, 6.589E-01, 6.817E-01,
     &                               7.034E-01, 7.233E-01, 7.416E-01,
     &                               7.581E-01, 7.727E-01, 7.852E-01,
     &                               7.965E-01, 8.070E-01, 8.165E-01,
     &                               8.256E-01, 8.343E-01, 8.423E-01,
     &                               8.496E-01, 8.582E-01, 8.665E-01,
     &                               8.747E-01, 8.833E-01, 8.900E-01,
     &                               8.956E-01, 9.018E-01, 9.118E-01,
     &                               9.179E-01, 9.230E-01, 9.268E-01,
     &                               9.301E-01, 9.332E-01, 9.362E-01,
     &                               9.394E-01, 9.419E-01, 9.433E-01,
     &                               9.444E-01, 9.450E-01, 9.456E-01,
     &                               9.463E-01, 9.471E-01, 9.478E-01,
     &                               9.484E-01, 9.489E-01, 9.494E-01,
     &                               9.499E-01, 9.504E-01, 9.509E-01/
      DATA (UV_B_EST(57,I),I=1,60) / 6.014E-04, 6.945E-02, 1.357E-01,
     &                               1.972E-01, 2.552E-01, 3.099E-01,
     &                               3.608E-01, 4.087E-01, 4.529E-01,
     &                               4.943E-01, 5.323E-01, 5.675E-01,
     &                               6.000E-01, 6.318E-01, 6.612E-01,
     &                               6.881E-01, 7.124E-01, 7.345E-01,
     &                               7.542E-01, 7.720E-01, 7.888E-01,
     &                               8.039E-01, 8.175E-01, 8.298E-01,
     &                               8.397E-01, 8.487E-01, 8.568E-01,
     &                               8.641E-01, 8.711E-01, 8.777E-01,
     &                               8.836E-01, 8.891E-01, 8.960E-01,
     &                               9.023E-01, 9.091E-01, 9.154E-01,
     &                               9.200E-01, 9.242E-01, 9.313E-01,
     &                               9.369E-01, 9.409E-01, 9.442E-01,
     &                               9.468E-01, 9.492E-01, 9.514E-01,
     &                               9.539E-01, 9.563E-01, 9.574E-01,
     &                               9.584E-01, 9.588E-01, 9.593E-01,
     &                               9.597E-01, 9.603E-01, 9.608E-01,
     &                               9.613E-01, 9.616E-01, 9.620E-01,
     &                               9.624E-01, 9.627E-01, 9.631E-01/
      DATA (UV_B_EST(58,I),I=1,60) / 1.423E-03, 7.511E-02, 1.493E-01,
     &                               2.247E-01, 3.003E-01, 3.631E-01,
     &                               4.194E-01, 4.721E-01, 5.210E-01,
     &                               5.652E-01, 6.058E-01, 6.424E-01,
     &                               6.750E-01, 7.046E-01, 7.315E-01,
     &                               7.557E-01, 7.783E-01, 7.983E-01,
     &                               8.157E-01, 8.310E-01, 8.443E-01,
     &                               8.572E-01, 8.688E-01, 8.788E-01,
     &                               8.875E-01, 8.944E-01, 9.004E-01,
     &                               9.059E-01, 9.109E-01, 9.155E-01,
     &                               9.198E-01, 9.237E-01, 9.277E-01,
     &                               9.324E-01, 9.368E-01, 9.416E-01,
     &                               9.454E-01, 9.485E-01, 9.524E-01,
     &                               9.569E-01, 9.597E-01, 9.621E-01,
     &                               9.640E-01, 9.657E-01, 9.672E-01,
     &                               9.689E-01, 9.706E-01, 9.717E-01,
     &                               9.724E-01, 9.727E-01, 9.730E-01,
     &                               9.733E-01, 9.736E-01, 9.739E-01,
     &                               9.742E-01, 9.745E-01, 9.747E-01,
     &                               9.749E-01, 9.752E-01, 9.754E-01/
      DATA (UV_B_EST(59,I),I=1,60) /-1.153E-05, 7.311E-02, 1.481E-01,
     &                               2.226E-01, 2.976E-01, 3.733E-01,
     &                               4.497E-01, 5.266E-01, 6.042E-01,
     &                               6.622E-01, 7.041E-01, 7.409E-01,
     &                               7.727E-01, 8.003E-01, 8.246E-01,
     &                               8.460E-01, 8.649E-01, 8.812E-01,
     &                               8.941E-01, 9.045E-01, 9.132E-01,
     &                               9.205E-01, 9.272E-01, 9.337E-01,
     &                               9.392E-01, 9.439E-01, 9.477E-01,
     &                               9.511E-01, 9.541E-01, 9.568E-01,
     &                               9.591E-01, 9.610E-01, 9.626E-01,
     &                               9.648E-01, 9.670E-01, 9.696E-01,
     &                               9.720E-01, 9.738E-01, 9.754E-01,
     &                               9.779E-01, 9.794E-01, 9.807E-01,
     &                               9.817E-01, 9.827E-01, 9.836E-01,
     &                               9.844E-01, 9.854E-01, 9.861E-01,
     &                               9.865E-01, 9.867E-01, 9.868E-01,
     &                               9.869E-01, 9.870E-01, 9.871E-01,
     &                               9.873E-01, 9.874E-01, 9.875E-01,
     &                               9.876E-01, 9.877E-01, 9.878E-01/
      DATA (UV_B_EST(60,I),I=1,60) / 6.452E-04, 7.262E-02, 1.455E-01,
     &                               2.191E-01, 2.931E-01, 3.676E-01,
     &                               4.427E-01, 5.199E-01, 5.960E-01,
     &                               6.744E-01, 7.544E-01, 8.366E-01,
     &                               9.112E-01, 9.329E-01, 9.501E-01,
     &                               9.642E-01, 9.758E-01, 9.856E-01,
     &                               9.939E-01, 1.001E+00, 1.002E+00,
     &                               1.001E+00, 1.001E+00, 1.000E+00,
     &                               1.000E+00, 1.000E+00, 1.000E+00,
     &                               1.001E+00, 1.001E+00, 1.002E+00,
     &                               1.002E+00, 1.002E+00, 1.001E+00,
     &                               1.000E+00, 1.000E+00, 9.998E-01,
     &                               1.000E+00, 1.000E+00, 1.000E+00,
     &                               1.000E+00, 1.000E+00, 1.000E+00,
     &                               1.000E+00, 1.000E+00, 1.000E+00,
     &                               1.001E+00, 1.001E+00, 1.001E+00,
     &                               1.001E+00, 1.001E+00, 1.001E+00,
     &                               1.001E+00, 1.001E+00, 1.000E+00,
     &                               1.000E+00, 1.000E+00, 1.000E+00,
     &                               1.000E+00, 1.000E+00, 1.000E+00/
      END





c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c BLOCK DATA XIONDATA_1keV
c
c     Table of the X-ray ionization rate. The integral 
c     for the ionization rate is evaluated between 1 keV
c     and 100 keV
c
      BLOCK DATA XIONDATA_1keV
            
      IMPLICIT NONE
      
c     Common-Block 'X_IONDAT_100eV':
c       X-ray ionization rate at an X-ray flux of 1 erg s^-1 cm^-2
c       for different values of the attenuating column density and
c       the plasma temperature.
c       X_SEC: Table of ionization rates s^-1
c       X_NH: Griding along the axis for the attenuating
c             column density cm^-2
c       X_TX: Griding along the axis for the plasma temperature K
      REAL X_SEC(40,40), X_NH(40), X_TX(40)
      COMMON/X_IONDAT_1keV/ X_SEC, X_NH, X_TX
      
      INTEGER I
      
      DATA (X_SEC(1,I),I=1,40) /  5.606E-12, 5.601E-12, 5.595E-12,
     &                            5.586E-12, 5.573E-12, 5.554E-12,
     &                            5.528E-12, 5.492E-12, 5.440E-12,
     &                            5.366E-12, 5.264E-12, 5.121E-12,
     &                            4.924E-12, 4.656E-12, 4.299E-12,
     &                            3.838E-12, 3.265E-12, 2.594E-12,
     &                            1.869E-12, 1.172E-12, 6.042E-13,
     &                            2.358E-13, 6.218E-14, 9.484E-15,
     &                            6.823E-16, 1.828E-17, 1.537E-19,
     &                            4.331E-22, 4.251E-25, 1.668E-28,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(2,I),I=1,40) /  5.516E-12, 5.511E-12, 5.505E-12,
     &                            5.496E-12, 5.484E-12, 5.466E-12,
     &                            5.441E-12, 5.405E-12, 5.355E-12,
     &                            5.283E-12, 5.184E-12, 5.044E-12,
     &                            4.853E-12, 4.592E-12, 4.245E-12,
     &                            3.795E-12, 3.235E-12, 2.578E-12,
     &                            1.866E-12, 1.179E-12, 6.141E-13,
     &                            2.436E-13, 6.600E-14, 1.055E-14,
     &                            8.301E-16, 2.679E-17, 3.303E-19,
     &                            1.677E-21, 3.606E-24, 5.438E-27,
     &                            5.378E-30, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(3,I),I=1,40) /  5.414E-12, 5.409E-12, 5.403E-12,
     &                            5.395E-12, 5.383E-12, 5.366E-12,
     &                            5.341E-12, 5.307E-12, 5.258E-12,
     &                            5.189E-12, 5.092E-12, 4.957E-12,
     &                            4.772E-12, 4.519E-12, 4.183E-12,
     &                            3.746E-12, 3.202E-12, 2.561E-12,
     &                            1.864E-12, 1.187E-12, 6.259E-13,
     &                            2.533E-13, 7.090E-14, 1.201E-14,
     &                            1.056E-15, 4.255E-17, 7.806E-19,
     &                            6.923E-21, 3.601E-23, 1.498E-25,
     &                            2.991E-28, 5.771E-32, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(4,I),I=1,40) /  5.297E-12, 5.293E-12, 5.287E-12,
     &                            5.279E-12, 5.267E-12, 5.251E-12,
     &                            5.227E-12, 5.194E-12, 5.147E-12,
     &                            5.081E-12, 4.988E-12, 4.858E-12,
     &                            4.679E-12, 4.436E-12, 4.111E-12,
     &                            3.689E-12, 3.162E-12, 2.540E-12,
     &                            1.860E-12, 1.196E-12, 6.398E-13,
     &                            2.650E-13, 7.714E-14, 1.401E-14,
     &                            1.405E-15, 7.231E-17, 1.964E-18,
     &                            3.102E-20, 3.664E-22, 3.009E-24,
     &                            1.112E-26, 2.307E-29, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(5,I),I=1,40) /  5.164E-12, 5.160E-12, 5.155E-12,
     &                            5.147E-12, 5.136E-12, 5.120E-12,
     &                            5.098E-12, 5.066E-12, 5.021E-12,
     &                            4.957E-12, 4.868E-12, 4.744E-12,
     &                            4.573E-12, 4.340E-12, 4.028E-12,
     &                            3.623E-12, 3.116E-12, 2.514E-12,
     &                            1.855E-12, 1.205E-12, 6.557E-13,
     &                            2.790E-13, 8.501E-14, 1.674E-14,
     &                            1.945E-15, 1.293E-16, 5.204E-18,
     &                            1.461E-19, 3.240E-21, 4.449E-23,
     &                            3.293E-25, 1.828E-27, 5.541E-30,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(6,I),I=1,40) /  5.015E-12, 5.011E-12, 5.006E-12,
     &                            4.998E-12, 4.988E-12, 4.973E-12,
     &                            4.951E-12, 4.921E-12, 4.878E-12,
     &                            4.818E-12, 4.734E-12, 4.616E-12,
     &                            4.453E-12, 4.231E-12, 3.934E-12,
     &                            3.547E-12, 3.061E-12, 2.484E-12,
     &                            1.847E-12, 1.215E-12, 6.735E-13,
     &                            2.955E-13, 9.479E-14, 2.043E-14,
     &                            2.782E-15, 2.403E-16, 1.429E-17,
     &                            6.646E-19, 2.355E-20, 5.224E-22,
     &                            7.896E-24, 8.746E-26, 5.323E-28,
     &                            2.268E-30, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(7,I),I=1,40) /  4.847E-12, 4.844E-12, 4.839E-12,
     &                            4.832E-12, 4.822E-12, 4.808E-12,
     &                            4.788E-12, 4.759E-12, 4.719E-12,
     &                            4.662E-12, 4.582E-12, 4.471E-12,
     &                            4.317E-12, 4.107E-12, 3.826E-12,
     &                            3.459E-12, 2.997E-12, 2.446E-12,
     &                            1.835E-12, 1.223E-12, 6.927E-13,
     &                            3.144E-13, 1.068E-13, 2.538E-14,
     &                            4.080E-15, 4.577E-16, 3.920E-17,
     &                            2.741E-18, 1.427E-19, 5.119E-21,
     &                            1.410E-22, 2.683E-24, 3.394E-26,
     &                            3.527E-28, 2.548E-30, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(8,I),I=1,40) /  4.662E-12, 4.659E-12, 4.654E-12,
     &                            4.648E-12, 4.639E-12, 4.625E-12,
     &                            4.607E-12, 4.580E-12, 4.542E-12,
     &                            4.489E-12, 4.414E-12, 4.310E-12,
     &                            4.166E-12, 3.969E-12, 3.705E-12,
     &                            3.359E-12, 2.924E-12, 2.401E-12,
     &                            1.818E-12, 1.230E-12, 7.127E-13,
     &                            3.356E-13, 1.212E-13, 3.193E-14,
     &                            6.079E-15, 8.766E-16, 1.033E-16,
     &                            1.003E-17, 7.397E-19, 4.153E-20,
     &                            1.835E-21, 5.783E-23, 1.431E-24,
     &                            2.878E-26, 4.199E-28, 4.139E-30,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(9,I),I=1,40) /  4.459E-12, 4.456E-12, 4.452E-12,
     &                            4.446E-12, 4.437E-12, 4.425E-12,
     &                            4.408E-12, 4.383E-12, 4.348E-12,
     &                            4.299E-12, 4.229E-12, 4.132E-12,
     &                            3.998E-12, 3.815E-12, 3.569E-12,
     &                            3.246E-12, 2.838E-12, 2.347E-12,
     &                            1.796E-12, 1.234E-12, 7.324E-13,
     &                            3.587E-13, 1.381E-13, 4.047E-14,
     &                            9.094E-15, 1.651E-15, 2.547E-16,
     &                            3.267E-17, 3.321E-18, 2.739E-19,
     &                            1.796E-20, 9.149E-22, 3.933E-23,
     &                            1.360E-24, 3.535E-26, 6.652E-28,
     &                            9.225E-30, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(10,I),I=1,40) / 4.239E-12, 4.236E-12, 4.232E-12,
     &                            4.227E-12, 4.219E-12, 4.208E-12,
     &                            4.192E-12, 4.169E-12, 4.137E-12,
     &                            4.091E-12, 4.027E-12, 3.938E-12,
     &                            3.814E-12, 3.645E-12, 3.418E-12,
     &                            3.119E-12, 2.740E-12, 2.283E-12,
     &                            1.765E-12, 1.233E-12, 7.508E-13,
     &                            3.830E-13, 1.575E-13, 5.131E-14,
     &                            1.349E-14, 3.005E-15, 5.816E-16,
     &                            9.519E-17, 1.291E-17, 1.471E-18,
     &                            1.372E-19, 1.071E-20, 7.198E-22,
     &                            3.934E-23, 1.678E-24, 5.543E-26,
     &                            1.435E-27, 2.733E-29, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(11,I),I=1,40) / 4.003E-12, 4.000E-12, 3.997E-12,
     &                            3.992E-12, 3.984E-12, 3.974E-12,
     &                            3.960E-12, 3.939E-12, 3.909E-12,
     &                            3.868E-12, 3.810E-12, 3.728E-12,
     &                            3.615E-12, 3.461E-12, 3.252E-12,
     &                            2.978E-12, 2.630E-12, 2.207E-12,
     &                            1.725E-12, 1.226E-12, 7.663E-13,
     &                            4.075E-13, 1.791E-13, 6.460E-14,
     &                            1.962E-14, 5.218E-15, 1.227E-15,
     &                            2.491E-16, 4.347E-17, 6.529E-18,
     &                            8.352E-19, 9.350E-20, 9.093E-21,
     &                            7.347E-22, 4.825E-23, 2.586E-24,
     &                            1.112E-25, 3.557E-27, 7.483E-29,
     &                            7.032E-31/
      DATA (X_SEC(12,I),I=1,40) / 3.753E-12, 3.750E-12, 3.747E-12,
     &                            3.743E-12, 3.736E-12, 3.727E-12,
     &                            3.714E-12, 3.695E-12, 3.668E-12,
     &                            3.631E-12, 3.578E-12, 3.504E-12,
     &                            3.402E-12, 3.262E-12, 3.073E-12,
     &                            2.824E-12, 2.506E-12, 2.119E-12,
     &                            1.676E-12, 1.211E-12, 7.773E-13,
     &                            4.309E-13, 2.020E-13, 8.023E-14,
     &                            2.774E-14, 8.596E-15, 2.391E-15,
     &                            5.874E-16, 1.274E-16, 2.434E-17,
     &                            4.105E-18, 6.207E-19, 8.240E-20,
     &                            9.323E-21, 8.907E-22, 7.172E-23,
     &                            4.671E-24, 2.245E-25, 6.745E-27,
     &                            9.283E-29/
      DATA (X_SEC(13,I),I=1,40) / 3.491E-12, 3.489E-12, 3.487E-12,
     &                            3.482E-12, 3.477E-12, 3.468E-12,
     &                            3.456E-12, 3.440E-12, 3.416E-12,
     &                            3.382E-12, 3.335E-12, 3.269E-12,
     &                            3.177E-12, 3.052E-12, 2.882E-12,
     &                            2.658E-12, 2.371E-12, 2.019E-12,
     &                            1.615E-12, 1.187E-12, 7.824E-13,
     &                            4.516E-13, 2.252E-13, 9.775E-14,
     &                            3.789E-14, 1.340E-14, 4.316E-15,
     &                            1.252E-15, 3.277E-16, 7.726E-17,
     &                            1.651E-17, 3.211E-18, 5.558E-19,
     &                            8.413E-20, 1.107E-20, 1.246E-21,
     &                            1.128E-22, 7.313E-24, 2.785E-25,
     &                            4.513E-27/
      DATA (X_SEC(14,I),I=1,40) / 3.223E-12, 3.221E-12, 3.218E-12,
     &                            3.215E-12, 3.209E-12, 3.202E-12,
     &                            3.192E-12, 3.177E-12, 3.156E-12,
     &                            3.126E-12, 3.084E-12, 3.025E-12,
     &                            2.944E-12, 2.832E-12, 2.681E-12,
     &                            2.481E-12, 2.224E-12, 1.909E-12,
     &                            1.544E-12, 1.154E-12, 7.802E-13,
     &                            4.682E-13, 2.474E-13, 1.163E-13,
     &                            4.985E-14, 1.977E-14, 7.232E-15,
     &                            2.425E-15, 7.466E-16, 2.112E-16,
     &                            5.518E-17, 1.327E-17, 2.887E-18,
     &                            5.610E-19, 9.649E-20, 1.420E-20,
     &                            1.649E-21, 1.320E-22, 5.901E-24,
     &                            1.091E-25/
      DATA (X_SEC(15,I),I=1,40) / 2.950E-12, 2.949E-12, 2.946E-12,
     &                            2.943E-12, 2.939E-12, 2.932E-12,
     &                            2.923E-12, 2.910E-12, 2.891E-12,
     &                            2.865E-12, 2.828E-12, 2.777E-12,
     &                            2.705E-12, 2.607E-12, 2.473E-12,
     &                            2.296E-12, 2.069E-12, 1.789E-12,
     &                            1.462E-12, 1.111E-12, 7.696E-13,
     &                            4.792E-13, 2.670E-13, 1.349E-13,
     &                            6.309E-14, 2.763E-14, 1.129E-14,
     &                            4.293E-15, 1.521E-15, 5.032E-16,
     &                            1.559E-16, 4.497E-17, 1.189E-17,
     &                            2.856E-18, 6.110E-19, 1.107E-19,
     &                            1.544E-20, 1.438E-21, 7.340E-23,
     &                            1.705E-24/
      DATA (X_SEC(16,I),I=1,40) / 2.679E-12, 2.677E-12, 2.675E-12,
     &                            2.673E-12, 2.669E-12, 2.663E-12,
     &                            2.655E-12, 2.643E-12, 2.627E-12,
     &                            2.604E-12, 2.572E-12, 2.527E-12,
     &                            2.465E-12, 2.379E-12, 2.262E-12,
     &                            2.107E-12, 1.908E-12, 1.661E-12,
     &                            1.372E-12, 1.058E-12, 7.504E-13,
     &                            4.833E-13, 2.827E-13, 1.522E-13,
     &                            7.680E-14, 3.665E-14, 1.649E-14,
     &                            6.990E-15, 2.796E-15, 1.058E-15,
     &                            3.792E-16, 1.277E-16, 3.992E-17,
     &                            1.143E-17, 2.910E-18, 6.182E-19,
     &                            9.907E-20, 1.048E-20, 6.380E-22,
     &                            2.506E-23/
      DATA (X_SEC(17,I),I=1,40) / 2.412E-12, 2.411E-12, 2.409E-12,
     &                            2.406E-12, 2.403E-12, 2.398E-12,
     &                            2.391E-12, 2.381E-12, 2.367E-12,
     &                            2.348E-12, 2.320E-12, 2.281E-12,
     &                            2.227E-12, 2.152E-12, 2.051E-12,
     &                            1.917E-12, 1.743E-12, 1.527E-12,
     &                            1.274E-12, 9.972E-13, 7.225E-13,
     &                            4.801E-13, 2.933E-13, 1.669E-13,
     &                            9.000E-14, 4.627E-14, 2.264E-14,
     &                            1.054E-14, 4.683E-15, 1.989E-15,
     &                            8.072E-16, 3.104E-16, 1.118E-16,
     &                            3.697E-17, 1.080E-17, 2.603E-18,
     &                            4.697E-19, 5.750E-20, 4.764E-21,
     &                            3.855E-22/
      DATA (X_SEC(18,I),I=1,40) / 2.154E-12, 2.153E-12, 2.151E-12,
     &                            2.149E-12, 2.146E-12, 2.142E-12,
     &                            2.136E-12, 2.128E-12, 2.116E-12,
     &                            2.099E-12, 2.075E-12, 2.041E-12,
     &                            1.995E-12, 1.931E-12, 1.844E-12,
     &                            1.728E-12, 1.578E-12, 1.392E-12,
     &                            1.171E-12, 9.293E-13, 6.867E-13,
     &                            4.694E-13, 2.980E-13, 1.782E-13,
     &                            1.017E-13, 5.578E-14, 2.935E-14,
     &                            1.483E-14, 7.214E-15, 3.384E-15,
     &                            1.528E-15, 6.583E-16, 2.669E-16,
     &                            9.937E-17, 3.250E-17, 8.723E-18,
     &                            1.781E-18, 2.673E-19, 3.398E-20,
     &                            4.938E-21/
      DATA (X_SEC(19,I),I=1,40) / 1.908E-12, 1.907E-12, 1.906E-12,
     &                            1.904E-12, 1.901E-12, 1.898E-12,
     &                            1.893E-12, 1.885E-12, 1.875E-12,
     &                            1.861E-12, 1.840E-12, 1.812E-12,
     &                            1.772E-12, 1.718E-12, 1.643E-12,
     &                            1.544E-12, 1.416E-12, 1.256E-12,
     &                            1.066E-12, 8.563E-13, 6.442E-13,
     &                            4.518E-13, 2.969E-13, 1.852E-13,
     &                            1.112E-13, 6.448E-14, 3.614E-14,
     &                            1.961E-14, 1.032E-14, 5.276E-15,
     &                            2.612E-15, 1.240E-15, 5.553E-16,
     &                            2.282E-16, 8.223E-17, 2.452E-17,
     &                            5.793E-18, 1.128E-18, 2.180E-19,
     &                            4.658E-20/
      DATA (X_SEC(20,I),I=1,40) / 1.677E-12, 1.676E-12, 1.675E-12,
     &                            1.673E-12, 1.671E-12, 1.668E-12,
     &                            1.664E-12, 1.658E-12, 1.649E-12,
     &                            1.637E-12, 1.620E-12, 1.596E-12,
     &                            1.562E-12, 1.516E-12, 1.453E-12,
     &                            1.369E-12, 1.260E-12, 1.123E-12,
     &                            9.608E-13, 7.805E-13, 5.968E-13,
     &                            4.282E-13, 2.900E-13, 1.879E-13,
     &                            1.178E-13, 7.174E-14, 4.250E-14,
     &                            2.454E-14, 1.383E-14, 7.620E-15,
     &                            4.084E-15, 2.107E-15, 1.027E-15,
     &                            4.596E-16, 1.811E-16, 6.032E-17,
     &                            1.689E-17, 4.328E-18, 1.161E-18,
     &                            3.229E-19/
      DATA (X_SEC(21,I),I=1,40) / 1.463E-12, 1.462E-12, 1.461E-12,
     &                            1.460E-12, 1.458E-12, 1.456E-12,
     &                            1.452E-12, 1.447E-12, 1.440E-12,
     &                            1.430E-12, 1.415E-12, 1.395E-12,
     &                            1.366E-12, 1.327E-12, 1.274E-12,
     &                            1.203E-12, 1.111E-12, 9.955E-13,
     &                            8.577E-13, 7.039E-13, 5.462E-13,
     &                            4.000E-13, 2.783E-13, 1.863E-13,
     &                            1.213E-13, 7.717E-14, 4.800E-14,
     &                            2.926E-14, 1.752E-14, 1.030E-14,
     &                            5.914E-15, 3.277E-15, 1.719E-15,
     &                            8.304E-16, 3.571E-16, 1.338E-16,
     &                            4.476E-17, 1.465E-17, 4.999E-18,
     &                            1.698E-18/
      DATA (X_SEC(22,I),I=1,40) / 1.268E-12, 1.267E-12, 1.266E-12,
     &                            1.265E-12, 1.264E-12, 1.262E-12,
     &                            1.259E-12, 1.254E-12, 1.248E-12,
     &                            1.240E-12, 1.228E-12, 1.211E-12,
     &                            1.187E-12, 1.154E-12, 1.110E-12,
     &                            1.050E-12, 9.725E-13, 8.752E-13,
     &                            7.590E-13, 6.287E-13, 4.944E-13,
     &                            3.687E-13, 2.627E-13, 1.810E-13,
     &                            1.219E-13, 8.057E-14, 5.232E-14,
     &                            3.347E-14, 2.113E-14, 1.315E-14,
     &                            8.014E-15, 4.725E-15, 2.645E-15,
     &                            1.372E-15, 6.431E-16, 2.715E-16,
     &                            1.076E-16, 4.297E-17, 1.757E-17,
     &                            7.017E-18/
      DATA (X_SEC(23,I),I=1,40) / 1.092E-12, 1.091E-12, 1.091E-12,
     &                            1.090E-12, 1.089E-12, 1.087E-12,
     &                            1.084E-12, 1.081E-12, 1.076E-12,
     &                            1.068E-12, 1.058E-12, 1.044E-12,
     &                            1.024E-12, 9.970E-13, 9.598E-13,
     &                            9.100E-13, 8.452E-13, 7.638E-13,
     &                            6.662E-13, 5.565E-13, 4.429E-13,
     &                            3.357E-13, 2.443E-13, 1.728E-13,
     &                            1.199E-13, 8.195E-14, 5.530E-14,
     &                            3.692E-14, 2.442E-14, 1.597E-14,
     &                            1.025E-14, 6.386E-15, 3.790E-15,
     &                            2.101E-15, 1.072E-15, 5.083E-16,
     &                            2.339E-16, 1.093E-16, 5.153E-17,
     &                            2.356E-17/
      DATA (X_SEC(24,I),I=1,40) / 9.349E-13, 9.345E-13, 9.340E-13,
     &                            9.333E-13, 9.322E-13, 9.308E-13,
     &                            9.287E-13, 9.257E-13, 9.215E-13,
     &                            9.155E-13, 9.071E-13, 8.953E-13,
     &                            8.789E-13, 8.562E-13, 8.252E-13,
     &                            7.839E-13, 7.299E-13, 6.620E-13,
     &                            5.805E-13, 4.887E-13, 3.931E-13,
     &                            3.024E-13, 2.243E-13, 1.623E-13,
     &                            1.157E-13, 8.149E-14, 5.690E-14,
     &                            3.946E-14, 2.719E-14, 1.858E-14,
     &                            1.249E-14, 8.165E-15, 5.112E-15,
     &                            3.019E-15, 1.671E-15, 8.808E-16,
     &                            4.598E-16, 2.437E-16, 1.289E-16,
     &                            6.613E-17/
      DATA (X_SEC(25,I),I=1,40) / 7.966E-13, 7.963E-13, 7.959E-13,
     &                            7.953E-13, 7.944E-13, 7.932E-13,
     &                            7.914E-13, 7.890E-13, 7.855E-13,
     &                            7.806E-13, 7.736E-13, 7.638E-13,
     &                            7.502E-13, 7.314E-13, 7.058E-13,
     &                            6.715E-13, 6.268E-13, 5.704E-13,
     &                            5.026E-13, 4.260E-13, 3.461E-13,
     &                            2.698E-13, 2.035E-13, 1.504E-13,
     &                            1.097E-13, 7.946E-14, 5.720E-14,
     &                            4.103E-14, 2.932E-14, 2.082E-14,
     &                            1.458E-14, 9.958E-15, 6.550E-15,
     &                            4.104E-15, 2.449E-15, 1.419E-15,
     &                            8.230E-16, 4.828E-16, 2.813E-16,
     &                            1.592E-16/
      DATA (X_SEC(26,I),I=1,40) / 6.758E-13, 6.755E-13, 6.752E-13,
     &                            6.747E-13, 6.740E-13, 6.730E-13,
     &                            6.715E-13, 6.695E-13, 6.666E-13,
     &                            6.626E-13, 6.568E-13, 6.488E-13,
     &                            6.375E-13, 6.220E-13, 6.008E-13,
     &                            5.725E-13, 5.355E-13, 4.889E-13,
     &                            4.327E-13, 3.691E-13, 3.025E-13,
     &                            2.387E-13, 1.829E-13, 1.377E-13,
     &                            1.027E-13, 7.619E-14, 5.637E-14,
     &                            4.166E-14, 3.076E-14, 2.260E-14,
     &                            1.641E-14, 1.167E-14, 8.031E-15,
     &                            5.317E-15, 3.399E-15, 2.139E-15,
     &                            1.353E-15, 8.628E-16, 5.453E-16,
     &                            3.357E-16/
      DATA (X_SEC(27,I),I=1,40) / 5.712E-13, 5.710E-13, 5.707E-13,
     &                            5.703E-13, 5.697E-13, 5.689E-13,
     &                            5.677E-13, 5.660E-13, 5.636E-13,
     &                            5.603E-13, 5.556E-13, 5.489E-13,
     &                            5.397E-13, 5.269E-13, 5.095E-13,
     &                            4.861E-13, 4.556E-13, 4.171E-13,
     &                            3.707E-13, 3.181E-13, 2.629E-13,
     &                            2.097E-13, 1.630E-13, 1.248E-13,
     &                            9.490E-14, 7.201E-14, 5.462E-14,
     &                            4.147E-14, 3.152E-14, 2.388E-14,
     &                            1.792E-14, 1.321E-14, 9.483E-15,
     &                            6.606E-15, 4.492E-15, 3.032E-15,
     &                            2.063E-15, 1.410E-15, 9.550E-16,
     &                            6.320E-16/
      DATA (X_SEC(28,I),I=1,40) / 4.813E-13, 4.811E-13, 4.809E-13,
     &                            4.805E-13, 4.800E-13, 4.794E-13,
     &                            4.784E-13, 4.770E-13, 4.751E-13,
     &                            4.723E-13, 4.684E-13, 4.630E-13,
     &                            4.554E-13, 4.449E-13, 4.305E-13,
     &                            4.113E-13, 3.863E-13, 3.546E-13,
     &                            3.163E-13, 2.730E-13, 2.273E-13,
     &                            1.832E-13, 1.442E-13, 1.121E-13,
     &                            8.685E-14, 6.725E-14, 5.217E-14,
     &                            4.060E-14, 3.167E-14, 2.467E-14,
     &                            1.908E-14, 1.454E-14, 1.085E-14,
     &                            7.912E-15, 5.681E-15, 4.074E-15,
     &                            2.944E-15, 2.135E-15, 1.533E-15,
     &                            1.079E-15/
      DATA (X_SEC(29,I),I=1,40) / 4.045E-13, 4.044E-13, 4.042E-13,
     &                            4.039E-13, 4.035E-13, 4.029E-13,
     &                            4.021E-13, 4.010E-13, 3.994E-13,
     &                            3.972E-13, 3.940E-13, 3.895E-13,
     &                            3.833E-13, 3.746E-13, 3.629E-13,
     &                            3.471E-13, 3.265E-13, 3.005E-13,
     &                            2.691E-13, 2.334E-13, 1.957E-13,
     &                            1.592E-13, 1.269E-13, 1.001E-13,
     &                            7.885E-14, 6.221E-14, 4.927E-14,
     &                            3.921E-14, 3.132E-14, 2.503E-14,
     &                            1.989E-14, 1.562E-14, 1.207E-14,
     &                            9.180E-15, 6.915E-15, 5.219E-15,
     &                            3.969E-15, 3.024E-15, 2.283E-15,
     &                            1.695E-15/
      DATA (X_SEC(30,I),I=1,40) / 3.394E-13, 3.392E-13, 3.391E-13,
     &                            3.389E-13, 3.385E-13, 3.381E-13,
     &                            3.374E-13, 3.365E-13, 3.352E-13,
     &                            3.334E-13, 3.307E-13, 3.271E-13,
     &                            3.220E-13, 3.149E-13, 3.053E-13,
     &                            2.924E-13, 2.755E-13, 2.541E-13,
     &                            2.283E-13, 1.990E-13, 1.680E-13,
     &                            1.380E-13, 1.112E-13, 8.894E-14,
     &                            7.115E-14, 5.712E-14, 4.612E-14,
     &                            3.746E-14, 3.059E-14, 2.501E-14,
     &                            2.038E-14, 1.646E-14, 1.314E-14,
     &                            1.037E-14, 8.139E-15, 6.416E-15,
     &                            5.093E-15, 4.046E-15, 3.189E-15,
     &                            2.477E-15/
      DATA (X_SEC(31,I),I=1,40) / 2.845E-13, 2.844E-13, 2.842E-13,
     &                            2.841E-13, 2.838E-13, 2.834E-13,
     &                            2.829E-13, 2.821E-13, 2.811E-13,
     &                            2.796E-13, 2.774E-13, 2.744E-13,
     &                            2.702E-13, 2.645E-13, 2.566E-13,
     &                            2.460E-13, 2.322E-13, 2.147E-13,
     &                            1.935E-13, 1.695E-13, 1.440E-13,
     &                            1.193E-13, 9.721E-14, 7.877E-14,
     &                            6.394E-14, 5.219E-14, 4.289E-14,
     &                            3.551E-14, 2.958E-14, 2.471E-14,
     &                            2.060E-14, 1.707E-14, 1.402E-14,
     &                            1.143E-14, 9.304E-15, 7.609E-15,
     &                            6.261E-15, 5.154E-15, 4.212E-15,
     &                            3.400E-15/
      DATA (X_SEC(32,I),I=1,40) / 2.386E-13, 2.385E-13, 2.384E-13,
     &                            2.383E-13, 2.381E-13, 2.378E-13,
     &                            2.373E-13, 2.367E-13, 2.358E-13,
     &                            2.346E-13, 2.328E-13, 2.304E-13,
     &                            2.270E-13, 2.222E-13, 2.158E-13,
     &                            2.071E-13, 1.958E-13, 1.815E-13,
     &                            1.641E-13, 1.444E-13, 1.235E-13,
     &                            1.032E-13, 8.494E-14, 6.969E-14,
     &                            5.738E-14, 4.756E-14, 3.975E-14,
     &                            3.349E-14, 2.842E-14, 2.421E-14,
     &                            2.061E-14, 1.747E-14, 1.472E-14,
     &                            1.235E-14, 1.037E-14, 8.746E-15,
     &                            7.417E-15, 6.289E-15, 5.298E-15,
     &                            4.417E-15/
      DATA (X_SEC(33,I),I=1,40) / 2.008E-13, 2.007E-13, 2.006E-13,
     &                            2.005E-13, 2.003E-13, 2.001E-13,
     &                            1.997E-13, 1.992E-13, 1.985E-13,
     &                            1.975E-13, 1.960E-13, 1.940E-13,
     &                            1.912E-13, 1.873E-13, 1.820E-13,
     &                            1.749E-13, 1.656E-13, 1.538E-13,
     &                            1.396E-13, 1.234E-13, 1.062E-13,
     &                            8.943E-14, 7.440E-14, 6.178E-14,
     &                            5.155E-14, 4.336E-14, 3.681E-14,
     &                            3.153E-14, 2.721E-14, 2.358E-14,
     &                            2.046E-14, 1.770E-14, 1.526E-14,
     &                            1.312E-14, 1.130E-14, 9.782E-15,
     &                            8.504E-15, 7.390E-15, 6.385E-15,
     &                            5.466E-15/
      DATA (X_SEC(34,I),I=1,40) / 1.700E-13, 1.699E-13, 1.699E-13,
     &                            1.698E-13, 1.696E-13, 1.694E-13,
     &                            1.691E-13, 1.687E-13, 1.681E-13,
     &                            1.673E-13, 1.661E-13, 1.644E-13,
     &                            1.621E-13, 1.589E-13, 1.545E-13,
     &                            1.487E-13, 1.410E-13, 1.313E-13,
     &                            1.195E-13, 1.061E-13, 9.187E-14,
     &                            7.800E-14, 6.554E-14, 5.505E-14,
     &                            4.652E-14, 3.967E-14, 3.416E-14,
     &                            2.970E-14, 2.602E-14, 2.292E-14,
     &                            2.021E-14, 1.780E-14, 1.565E-14,
     &                            1.374E-14, 1.209E-14, 1.069E-14,
     &                            9.481E-15, 8.406E-15, 7.414E-15,
     &                            6.484E-15/
      DATA (X_SEC(35,I),I=1,40) / 1.454E-13, 1.453E-13, 1.453E-13,
     &                            1.452E-13, 1.451E-13, 1.449E-13,
     &                            1.447E-13, 1.443E-13, 1.438E-13,
     &                            1.431E-13, 1.421E-13, 1.407E-13,
     &                            1.388E-13, 1.361E-13, 1.325E-13,
     &                            1.276E-13, 1.212E-13, 1.131E-13,
     &                            1.033E-13, 9.212E-14, 8.026E-14,
     &                            6.867E-14, 5.825E-14, 4.946E-14,
     &                            4.230E-14, 3.653E-14, 3.187E-14,
     &                            2.808E-14, 2.494E-14, 2.226E-14,
     &                            1.992E-14, 1.782E-14, 1.592E-14,
     &                            1.422E-14, 1.273E-14, 1.145E-14,
     &                            1.032E-14, 9.300E-15, 8.338E-15,
     &                            7.420E-15/
      DATA (X_SEC(36,I),I=1,40) / 1.260E-13, 1.259E-13, 1.259E-13,
     &                            1.258E-13, 1.257E-13, 1.256E-13,
     &                            1.254E-13, 1.251E-13, 1.246E-13,
     &                            1.241E-13, 1.232E-13, 1.221E-13,
     &                            1.204E-13, 1.182E-13, 1.151E-13,
     &                            1.110E-13, 1.056E-13, 9.874E-14,
     &                            9.046E-14, 8.102E-14, 7.099E-14,
     &                            6.119E-14, 5.237E-14, 4.492E-14,
     &                            3.884E-14, 3.392E-14, 2.994E-14,
     &                            2.669E-14, 2.398E-14, 2.166E-14,
     &                            1.962E-14, 1.777E-14, 1.609E-14,
     &                            1.458E-14, 1.324E-14, 1.206E-14,
     &                            1.102E-14, 1.006E-14, 9.136E-15,
     &                            8.241E-15/
      DATA (X_SEC(37,I),I=1,40) / 1.108E-13, 1.108E-13, 1.108E-13,
     &                            1.107E-13, 1.106E-13, 1.105E-13,
     &                            1.103E-13, 1.101E-13, 1.097E-13,
     &                            1.092E-13, 1.085E-13, 1.075E-13,
     &                            1.061E-13, 1.042E-13, 1.015E-13,
     &                            9.798E-14, 9.335E-14, 8.749E-14,
     &                            8.039E-14, 7.229E-14, 6.369E-14,
     &                            5.528E-14, 4.770E-14, 4.129E-14,
     &                            3.605E-14, 3.180E-14, 2.835E-14,
     &                            2.553E-14, 2.317E-14, 2.114E-14,
     &                            1.934E-14, 1.770E-14, 1.620E-14,
     &                            1.484E-14, 1.363E-14, 1.256E-14,
     &                            1.159E-14, 1.068E-14, 9.804E-15,
     &                            8.938E-15/
      DATA (X_SEC(38,I),I=1,40) / 9.910E-14, 9.907E-14, 9.903E-14,
     &                            9.898E-14, 9.890E-14, 9.879E-14,
     &                            9.864E-14, 9.842E-14, 9.811E-14,
     &                            9.767E-14, 9.704E-14, 9.617E-14,
     &                            9.495E-14, 9.327E-14, 9.097E-14,
     &                            8.789E-14, 8.385E-14, 7.874E-14,
     &                            7.255E-14, 6.549E-14, 5.798E-14,
     &                            5.064E-14, 4.402E-14, 3.841E-14,
     &                            3.382E-14, 3.010E-14, 2.707E-14,
     &                            2.458E-14, 2.249E-14, 2.069E-14,
     &                            1.909E-14, 1.762E-14, 1.627E-14,
     &                            1.504E-14, 1.393E-14, 1.294E-14,
     &                            1.204E-14, 1.119E-14, 1.035E-14,
     &                            9.515E-15/
      DATA (X_SEC(39,I),I=1,40) / 9.001E-14, 8.999E-14, 8.995E-14,
     &                            8.991E-14, 8.984E-14, 8.974E-14,
     &                            8.960E-14, 8.941E-14, 8.913E-14,
     &                            8.874E-14, 8.819E-14, 8.742E-14,
     &                            8.634E-14, 8.484E-14, 8.281E-14,
     &                            8.007E-14, 7.649E-14, 7.196E-14,
     &                            6.647E-14, 6.020E-14, 5.354E-14,
     &                            4.702E-14, 4.114E-14, 3.615E-14,
     &                            3.207E-14, 2.875E-14, 2.604E-14,
     &                            2.381E-14, 2.194E-14, 2.032E-14,
     &                            1.887E-14, 1.754E-14, 1.631E-14,
     &                            1.518E-14, 1.416E-14, 1.324E-14,
     &                            1.240E-14, 1.159E-14, 1.079E-14,
     &                            9.985E-15/
      DATA (X_SEC(40,I),I=1,40) / 8.297E-14, 8.295E-14, 8.292E-14,
     &                            8.288E-14, 8.281E-14, 8.273E-14,
     &                            8.260E-14, 8.243E-14, 8.218E-14,
     &                            8.183E-14, 8.133E-14, 8.063E-14,
     &                            7.966E-14, 7.831E-14, 7.647E-14,
     &                            7.401E-14, 7.078E-14, 6.669E-14,
     &                            6.174E-14, 5.609E-14, 5.008E-14,
     &                            4.419E-14, 3.888E-14, 3.438E-14,
     &                            3.068E-14, 2.767E-14, 2.522E-14,
     &                            2.320E-14, 2.149E-14, 2.001E-14,
     &                            1.868E-14, 1.746E-14, 1.633E-14,
     &                            1.528E-14, 1.433E-14, 1.347E-14,
     &                            1.268E-14, 1.191E-14, 1.114E-14,
     &                            1.036E-14/
      DATA (X_NH(I),I=1,40) /     1.000E+19, 1.425E+19, 2.031E+19,
     &                            2.894E+19, 4.125E+19, 5.878E+19,
     &                            8.377E+19, 1.194E+20, 1.701E+20,
     &                            2.424E+20, 3.455E+20, 4.924E+20,
     &                            7.017E+20, 1.000E+21, 1.425E+21,
     &                            2.031E+21, 2.894E+21, 4.125E+21,
     &                            5.878E+21, 8.377E+21, 1.194E+22,
     &                            1.701E+22, 2.424E+22, 3.455E+22,
     &                            4.924E+22, 7.017E+22, 1.000E+23,
     &                            1.425E+23, 2.031E+23, 2.894E+23,
     &                            4.125E+23, 5.878E+23, 8.377E+23,
     &                            1.194E+24, 1.701E+24, 2.424E+24,
     &                            3.455E+24, 4.924E+24, 7.017E+24,
     &                            1.000E+25/
      DATA (X_TX(I),I=1,40) /     3.162E+05, 3.888E+05, 4.781E+05,
     &                            5.878E+05, 7.227E+05, 8.886E+05,
     &                            1.093E+06, 1.343E+06, 1.652E+06,
     &                            2.031E+06, 2.497E+06, 3.070E+06,
     &                            3.775E+06, 4.642E+06, 5.707E+06,
     &                            7.017E+06, 8.628E+06, 1.061E+07,
     &                            1.304E+07, 1.604E+07, 1.972E+07,
     &                            2.424E+07, 2.981E+07, 3.665E+07,
     &                            4.507E+07, 5.541E+07, 6.813E+07,
     &                            8.377E+07, 1.030E+08, 1.266E+08,
     &                            1.557E+08, 1.914E+08, 2.354E+08,
     &                            2.894E+08, 3.559E+08, 4.375E+08,
     &                            5.380E+08, 6.615E+08, 8.133E+08,
     &                            1.000E+09/
      END





c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c BLOCK DATA XIONDATA_100eV
c
c     Table of the X-ray ionization rate. The integral 
c     for the ionization rate is evaluated between 100 ev
c     and 100 keV
c
      BLOCK DATA XIONDATA_100eV
            
      IMPLICIT NONE
      
c     Common-Block 'X_IONDAT_100eV':
c       X-ray ionization rate at an X-ray flux of 1 erg s^-1 cm^-2
c       for different values of the attenuating column density and
c       the plasma temperature.
c       X_SEC: Table of ionization rates s^-1
c       X_NH: Griding along the axis for the attenuating
c             column density cm^-2
c       X_TX: Griding along the axis for the plasma temperature K
      REAL X_SEC(40,40), X_NH(40), X_TX(40)
      COMMON/X_IONDAT_100ev/ X_SEC, X_NH, X_TX
      
      INTEGER I      
            
      DATA (X_SEC(1,I),I=1,40) /  6.349E-10, 5.359E-10, 4.232E-10,
     &                            3.059E-10, 1.974E-10, 1.112E-10,
     &                            5.421E-11, 2.320E-11, 8.991E-12,
     &                            3.228E-12, 1.078E-12, 3.330E-13,
     &                            9.403E-14, 2.371E-14, 5.158E-15,
     &                            9.379E-16, 1.434E-16, 1.914E-17,
     &                            2.148E-18, 1.858E-19, 1.123E-20,
     &                            3.779E-22, 5.179E-24, 2.864E-26,
     &                            1.302E-28, 3.801E-31, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(2,I),I=1,40) /  5.961E-10, 5.071E-10, 4.052E-10,
     &                            2.982E-10, 1.978E-10, 1.162E-10,
     &                            6.026E-11, 2.806E-11, 1.206E-11,
     &                            4.864E-12, 1.845E-12, 6.542E-13,
     &                            2.138E-13, 6.302E-14, 1.628E-14,
     &                            3.626E-15, 7.062E-16, 1.216E-16,
     &                            1.751E-17, 1.935E-18, 1.447E-19,
     &                            5.940E-21, 1.168E-22, 1.498E-24,
     &                            1.559E-26, 1.079E-28, 6.728E-31,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(3,I),I=1,40) /  5.546E-10, 4.754E-10, 3.843E-10,
     &                            2.877E-10, 1.959E-10, 1.197E-10,
     &                            6.568E-11, 3.296E-11, 1.549E-11,
     &                            6.906E-12, 2.921E-12, 1.164E-12,
     &                            4.313E-13, 1.455E-13, 4.379E-14,
     &                            1.167E-14, 2.789E-15, 5.907E-16,
     &                            1.041E-16, 1.396E-17, 1.255E-18,
     &                            6.575E-20, 2.129E-21, 5.526E-23,
     &                            1.070E-24, 1.734E-26, 2.569E-28,
     &                            2.229E-30, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(4,I),I=1,40) /  5.112E-10, 4.414E-10, 3.607E-10,
     &                            2.745E-10, 1.916E-10, 1.214E-10,
     &                            7.008E-11, 3.755E-11, 1.907E-11,
     &                            9.263E-12, 4.301E-12, 1.895E-12,
     &                            7.817E-13, 2.968E-13, 1.022E-13,
     &                            3.184E-14, 9.023E-15, 2.264E-15,
     &                            4.709E-16, 7.445E-17, 8.070E-18,
     &                            5.745E-19, 3.078E-20, 1.347E-21,
     &                            4.676E-23, 1.505E-24, 3.781E-26,
     &                            5.947E-28, 6.985E-30, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(5,I),I=1,40) /  4.667E-10, 4.059E-10, 3.351E-10,
     &                            2.590E-10, 1.848E-10, 1.210E-10,
     &                            7.314E-11, 4.151E-11, 2.253E-11,
     &                            1.178E-11, 5.926E-12, 2.845E-12,
     &                            1.288E-12, 5.422E-13, 2.101E-13,
     &                            7.491E-14, 2.446E-14, 7.063E-15,
     &                            1.692E-15, 3.118E-16, 4.153E-17,
     &                            4.109E-18, 3.355E-19, 2.265E-20,
     &                            1.332E-21, 7.053E-23, 2.756E-24,
     &                            7.727E-26, 1.713E-27, 2.345E-29,
     &                            1.885E-32, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(6,I),I=1,40) /  4.222E-10, 3.695E-10, 3.080E-10,
     &                            2.414E-10, 1.759E-10, 1.187E-10,
     &                            7.466E-11, 4.455E-11, 2.561E-11,
     &                            1.427E-11, 7.687E-12, 3.974E-12,
     &                            1.951E-12, 8.994E-13, 3.864E-13,
     &                            1.545E-13, 5.681E-14, 1.849E-14,
     &                            5.024E-15, 1.077E-15, 1.781E-16,
     &                            2.394E-17, 2.757E-18, 2.734E-19,
     &                            2.467E-20, 1.910E-21, 1.123E-22,
     &                            5.221E-24, 1.849E-25, 4.100E-27,
     &                            6.197E-29, 6.584E-31, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(7,I),I=1,40) /  3.783E-10, 3.331E-10, 2.801E-10,
     &                            2.225E-10, 1.652E-10, 1.145E-10,
     &                            7.460E-11, 4.649E-11, 2.809E-11,
     &                            1.653E-11, 9.448E-12, 5.209E-12,
     &                            2.745E-12, 1.370E-12, 6.437E-13,
     &                            2.838E-13, 1.154E-13, 4.170E-14,
     &                            1.274E-14, 3.173E-15, 6.471E-16,
     &                            1.134E-16, 1.748E-17, 2.414E-18,
     &                            3.047E-19, 3.242E-20, 2.764E-21,
     &                            1.932E-22, 1.005E-23, 3.607E-25,
     &                            9.944E-27, 1.892E-28, 2.361E-30,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(8,I),I=1,40) /  3.360E-10, 2.975E-10, 2.522E-10,
     &                            2.027E-10, 1.531E-10, 1.086E-10,
     &                            7.304E-11, 4.728E-11, 2.982E-11,
     &                            1.840E-11, 1.107E-11, 6.455E-12,
     &                            3.619E-12, 1.936E-12, 9.831E-13,
     &                            4.712E-13, 2.091E-13, 8.290E-14,
     &                            2.828E-14, 8.136E-15, 2.009E-15,
     &                            4.412E-16, 8.737E-17, 1.589E-17,
     &                            2.626E-18, 3.689E-19, 4.339E-20,
     &                            4.214E-21, 3.105E-22, 1.744E-23,
     &                            7.710E-25, 2.431E-26, 6.014E-28,
     &                            1.199E-29, 0.000E+00, 0.000E+00,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(9,I),I=1,40) /  2.959E-10, 2.633E-10, 2.249E-10,
     &                            1.826E-10, 1.401E-10, 1.015E-10,
     &                            7.017E-11, 4.694E-11, 3.074E-11,
     &                            1.975E-11, 1.243E-11, 7.611E-12,
     &                            4.504E-12, 2.560E-12, 1.391E-12,
     &                            7.163E-13, 3.430E-13, 1.480E-13,
     &                            5.594E-14, 1.840E-14, 5.378E-15,
     &                            1.432E-15, 3.509E-16, 8.003E-17,
     &                            1.654E-17, 2.964E-18, 4.568E-19,
     &                            5.858E-20, 5.955E-21, 4.913E-22,
     &                            3.223E-23, 1.642E-24, 7.057E-26,
     &                            2.440E-27, 6.334E-29, 1.018E-30,
     &                            0.000E+00, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(10,I),I=1,40) / 2.585E-10, 2.311E-10, 1.987E-10,
     &                            1.629E-10, 1.267E-10, 9.356E-11,
     &                            6.624E-11, 4.561E-11, 3.085E-11,
     &                            2.055E-11, 1.344E-11, 8.590E-12,
     &                            5.331E-12, 3.195E-12, 1.840E-12,
     &                            1.008E-12, 5.162E-13, 2.404E-13,
     &                            9.991E-14, 3.708E-14, 1.254E-14,
     &                            3.935E-15, 1.155E-15, 3.173E-16,
     &                            7.935E-17, 1.755E-17, 3.396E-18,
     &                            5.558E-19, 7.540E-20, 8.592E-21,
     &                            8.017E-22, 6.255E-23, 4.205E-24,
     &                            2.298E-25, 9.801E-27, 3.237E-28,
     &                            8.245E-30, 0.000E+00, 0.000E+00,
     &                            0.000E+00/
      DATA (X_SEC(11,I),I=1,40) / 2.242E-10, 2.013E-10, 1.741E-10,
     &                            1.440E-10, 1.134E-10, 8.510E-11,
     &                            6.154E-11, 4.345E-11, 3.023E-11,
     &                            2.077E-11, 1.406E-11, 9.328E-12,
     &                            6.035E-12, 3.787E-12, 2.293E-12,
     &                            1.327E-12, 7.208E-13, 3.597E-13,
     &                            1.629E-13, 6.726E-14, 2.580E-14,
     &                            9.307E-15, 3.176E-15, 1.020E-15,
     &                            3.004E-16, 7.961E-17, 1.871E-17,
     &                            3.800E-18, 6.631E-19, 9.963E-20,
     &                            1.275E-20, 1.427E-21, 1.387E-22,
     &                            1.121E-23, 7.361E-25, 3.946E-26,
     &                            1.696E-27, 5.417E-29, 1.222E-30,
     &                            0.000E+00/
      DATA (X_SEC(12,I),I=1,40) / 1.931E-10, 1.740E-10, 1.513E-10,
     &                            1.262E-10, 1.004E-10, 7.648E-11,
     &                            5.634E-11, 4.067E-11, 2.900E-11,
     &                            2.048E-11, 1.428E-11, 9.789E-12,
     &                            6.567E-12, 4.289E-12, 2.713E-12,
     &                            1.646E-12, 9.420E-13, 5.001E-13,
     &                            2.444E-13, 1.108E-13, 4.734E-14,
     &                            1.923E-14, 7.450E-15, 2.730E-15,
     &                            9.258E-16, 2.863E-16, 7.963E-17,
     &                            1.957E-17, 4.243E-18, 8.110E-19,
     &                            1.368E-19, 2.068E-20, 2.745E-21,
     &                            3.106E-22, 2.967E-23, 2.389E-24,
     &                            1.556E-25, 7.476E-27, 2.240E-28,
     &                            3.019E-30/
      DATA (X_SEC(13,I),I=1,40) / 1.653E-10, 1.495E-10, 1.306E-10,
     &                            1.096E-10, 8.810E-11, 6.799E-11,
     &                            5.092E-11, 3.746E-11, 2.730E-11,
     &                            1.974E-11, 1.413E-11, 9.964E-12,
     &                            6.898E-12, 4.664E-12, 3.064E-12,
     &                            1.936E-12, 1.161E-12, 6.513E-13,
     &                            3.405E-13, 1.674E-13, 7.840E-14,
     &                            3.522E-14, 1.519E-14, 6.228E-15,
     &                            2.385E-15, 8.425E-16, 2.713E-16,
     &                            7.872E-17, 2.060E-17, 4.857E-18,
     &                            1.038E-18, 2.018E-19, 3.494E-20,
     &                            5.289E-21, 6.961E-22, 7.835E-23,
     &                            7.091E-24, 4.596E-25, 1.749E-26,
     &                            2.829E-28/
      DATA (X_SEC(14,I),I=1,40) / 1.408E-10, 1.276E-10, 1.120E-10,
     &                            9.460E-11, 7.667E-11, 5.985E-11,
     &                            4.547E-11, 3.402E-11, 2.527E-11,
     &                            1.865E-11, 1.365E-11, 9.871E-12,
     &                            7.022E-12, 4.892E-12, 3.320E-12,
     &                            2.175E-12, 1.359E-12, 8.000E-13,
     &                            4.437E-13, 2.340E-13, 1.185E-13,
     &                            5.801E-14, 2.740E-14, 1.237E-14,
     &                            5.257E-15, 2.083E-15, 7.621E-16,
     &                            2.556E-16, 7.868E-17, 2.226E-17,
     &                            5.815E-18, 1.399E-18, 3.042E-19,
     &                            5.912E-20, 1.017E-20, 1.496E-21,
     &                            1.737E-22, 1.390E-23, 6.211E-25,
     &                            1.147E-26/
      DATA (X_SEC(15,I),I=1,40) / 1.192E-10, 1.084E-10, 9.550E-11,
     &                            8.109E-11, 6.622E-11, 5.221E-11,
     &                            4.017E-11, 3.050E-11, 2.302E-11,
     &                            1.730E-11, 1.292E-11, 9.543E-12,
     &                            6.951E-12, 4.970E-12, 3.470E-12,
     &                            2.345E-12, 1.518E-12, 9.329E-13,
     &                            5.448E-13, 3.050E-13, 1.653E-13,
     &                            8.701E-14, 4.441E-14, 2.177E-14,
     &                            1.012E-14, 4.431E-15, 1.811E-15,
     &                            6.886E-16, 2.440E-16, 8.072E-17,
     &                            2.501E-17, 7.213E-18, 1.908E-18,
     &                            4.581E-19, 9.800E-20, 1.775E-20,
     &                            2.476E-21, 2.305E-22, 1.176E-23,
     &                            2.729E-25/
      DATA (X_SEC(16,I),I=1,40) / 1.005E-10, 9.159E-11, 8.097E-11,
     &                            6.908E-11, 5.679E-11, 4.518E-11,
     &                            3.514E-11, 2.702E-11, 2.069E-11,
     &                            1.579E-11, 1.199E-11, 9.028E-12,
     &                            6.711E-12, 4.907E-12, 3.510E-12,
     &                            2.438E-12, 1.628E-12, 1.039E-12,
     &                            6.342E-13, 3.740E-13, 2.146E-13,
     &                            1.202E-13, 6.556E-14, 3.449E-14,
     &                            1.734E-14, 8.272E-15, 3.723E-15,
     &                            1.578E-15, 6.313E-16, 2.389E-16,
     &                            8.560E-17, 2.883E-17, 9.011E-18,
     &                            2.580E-18, 6.568E-19, 1.395E-19,
     &                            2.235E-20, 2.363E-21, 1.438E-22,
     &                            5.651E-24/
      DATA (X_SEC(17,I),I=1,40) / 8.435E-11, 7.704E-11, 6.831E-11,
     &                            5.853E-11, 4.840E-11, 3.880E-11,
     &                            3.047E-11, 2.369E-11, 1.837E-11,
     &                            1.421E-11, 1.095E-11, 8.376E-12,
     &                            6.337E-12, 4.723E-12, 3.451E-12,
     &                            2.454E-12, 1.684E-12, 1.110E-12,
     &                            7.044E-13, 4.342E-13, 2.618E-13,
     &                            1.546E-13, 8.923E-14, 4.990E-14,
     &                            2.683E-14, 1.379E-14, 6.746E-15,
     &                            3.142E-15, 1.396E-15, 5.928E-16,
     &                            2.406E-16, 9.252E-17, 3.331E-17,
     &                            1.102E-17, 3.220E-18, 7.755E-19,
     &                            1.399E-19, 1.712E-20, 1.418E-21,
     &                            1.149E-22/
      DATA (X_SEC(18,I),I=1,40) / 7.053E-11, 6.454E-11, 5.738E-11,
     &                            4.935E-11, 4.102E-11, 3.310E-11,
     &                            2.621E-11, 2.058E-11, 1.613E-11,
     &                            1.262E-11, 9.850E-12, 7.639E-12,
     &                            5.867E-12, 4.446E-12, 3.308E-12,
     &                            2.401E-12, 1.687E-12, 1.144E-12,
     &                            7.507E-13, 4.807E-13, 3.022E-13,
     &                            1.868E-13, 1.132E-13, 6.672E-14,
     &                            3.801E-14, 2.084E-14, 1.097E-14,
     &                            5.541E-15, 2.695E-15, 1.264E-15,
     &                            5.709E-16, 2.460E-16, 9.971E-17,
     &                            3.713E-17, 1.214E-17, 3.258E-18,
     &                            6.652E-19, 9.979E-20, 1.269E-20,
     &                            1.845E-21/
      DATA (X_SEC(19,I),I=1,40) / 5.877E-11, 5.386E-11, 4.800E-11,
     &                            4.142E-11, 3.458E-11, 2.807E-11,
     &                            2.238E-11, 1.772E-11, 1.402E-11,
     &                            1.108E-11, 8.743E-12, 6.861E-12,
     &                            5.339E-12, 4.104E-12, 3.102E-12,
     &                            2.291E-12, 1.644E-12, 1.142E-12,
     &                            7.713E-13, 5.104E-13, 3.327E-13,
     &                            2.138E-13, 1.350E-13, 8.331E-14,
     &                            4.992E-14, 2.895E-14, 1.623E-14,
     &                            8.803E-15, 4.634E-15, 2.369E-15,
     &                            1.173E-15, 5.567E-16, 2.493E-16,
     &                            1.025E-16, 3.692E-17, 1.101E-17,
     &                            2.600E-18, 5.065E-19, 9.789E-20,
     &                            2.092E-20/
      DATA (X_SEC(20,I),I=1,40) / 4.881E-11, 4.480E-11, 4.001E-11,
     &                            3.462E-11, 2.901E-11, 2.367E-11,
     &                            1.900E-11, 1.515E-11, 1.208E-11,
     &                            9.630E-12, 7.669E-12, 6.080E-12,
     &                            4.784E-12, 3.722E-12, 2.852E-12,
     &                            2.139E-12, 1.562E-12, 1.108E-12,
     &                            7.676E-13, 5.224E-13, 3.512E-13,
     &                            2.334E-13, 1.529E-13, 9.809E-14,
     &                            6.141E-14, 3.741E-14, 2.216E-14,
     &                            1.279E-14, 7.213E-15, 3.973E-15,
     &                            2.130E-15, 1.098E-15, 5.355E-16,
     &                            2.396E-16, 9.443E-17, 3.144E-17,
     &                            8.802E-18, 2.256E-18, 6.053E-19,
     &                            1.684E-19/
      DATA (X_SEC(21,I),I=1,40) / 4.043E-11, 3.716E-11, 3.324E-11,
     &                            2.883E-11, 2.425E-11, 1.987E-11,
     &                            1.603E-11, 1.286E-11, 1.032E-11,
     &                            8.293E-12, 6.657E-12, 5.324E-12,
     &                            4.229E-12, 3.326E-12, 2.578E-12,
     &                            1.959E-12, 1.453E-12, 1.050E-12,
     &                            7.430E-13, 5.181E-13, 3.578E-13,
     &                            2.447E-13, 1.653E-13, 1.098E-13,
     &                            7.144E-14, 4.544E-14, 2.826E-14,
     &                            1.723E-14, 1.032E-14, 6.066E-15,
     &                            3.483E-15, 1.929E-15, 1.012E-15,
     &                            4.890E-16, 2.103E-16, 7.874E-17,
     &                            2.635E-17, 8.623E-18, 2.943E-18,
     &                            9.998E-19/
      DATA (X_SEC(22,I),I=1,40) / 3.341E-11, 3.073E-11, 2.753E-11,
     &                            2.393E-11, 2.019E-11, 1.660E-11,
     &                            1.346E-11, 1.085E-11, 8.762E-12,
     &                            7.084E-12, 5.726E-12, 4.614E-12,
     &                            3.696E-12, 2.932E-12, 2.296E-12,
     &                            1.765E-12, 1.327E-12, 9.742E-13,
     &                            7.021E-13, 5.000E-13, 3.533E-13,
     &                            2.478E-13, 1.721E-13, 1.178E-13,
     &                            7.926E-14, 5.237E-14, 3.401E-14,
     &                            2.176E-14, 1.374E-14, 8.548E-15,
     &                            5.210E-15, 3.071E-15, 1.719E-15,
     &                            8.915E-16, 4.180E-16, 1.765E-16,
     &                            6.995E-17, 2.793E-17, 1.142E-17,
     &                            4.562E-18/
      DATA (X_SEC(23,I),I=1,40) / 2.754E-11, 2.536E-11, 2.275E-11,
     &                            1.981E-11, 1.675E-11, 1.382E-11,
     &                            1.125E-11, 9.113E-12, 7.394E-12,
     &                            6.010E-12, 4.887E-12, 3.963E-12,
     &                            3.197E-12, 2.556E-12, 2.019E-12,
     &                            1.567E-12, 1.192E-12, 8.874E-13,
     &                            6.500E-13, 4.714E-13, 3.398E-13,
     &                            2.436E-13, 1.733E-13, 1.218E-13,
     &                            8.446E-14, 5.773E-14, 3.896E-14,
     &                            2.601E-14, 1.720E-14, 1.125E-14,
     &                            7.224E-15, 4.499E-15, 2.670E-15,
     &                            1.480E-15, 7.553E-16, 3.580E-16,
     &                            1.647E-16, 7.700E-17, 3.630E-17,
     &                            1.660E-17/
      DATA (X_SEC(24,I),I=1,40) / 2.266E-11, 2.089E-11, 1.876E-11,
     &                            1.636E-11, 1.386E-11, 1.147E-11,
     &                            9.364E-12, 7.617E-12, 6.206E-12,
     &                            5.069E-12, 4.142E-12, 3.378E-12,
     &                            2.741E-12, 2.206E-12, 1.755E-12,
     &                            1.375E-12, 1.056E-12, 7.958E-13,
     &                            5.911E-13, 4.356E-13, 3.196E-13,
     &                            2.336E-13, 1.697E-13, 1.221E-13,
     &                            8.699E-14, 6.129E-14, 4.280E-14,
     &                            2.968E-14, 2.045E-14, 1.397E-14,
     &                            9.393E-15, 6.141E-15, 3.845E-15,
     &                            2.270E-15, 1.256E-15, 6.624E-16,
     &                            3.458E-16, 1.833E-16, 9.698E-17,
     &                            4.974E-17/
      DATA (X_SEC(25,I),I=1,40) / 1.862E-11, 1.717E-11, 1.543E-11,
     &                            1.348E-11, 1.144E-11, 9.491E-12,
     &                            7.770E-12, 6.342E-12, 5.186E-12,
     &                            4.252E-12, 3.490E-12, 2.859E-12,
     &                            2.332E-12, 1.888E-12, 1.512E-12,
     &                            1.193E-12, 9.248E-13, 7.043E-13,
     &                            5.296E-13, 3.957E-13, 2.949E-13,
     &                            2.192E-13, 1.623E-13, 1.193E-13,
     &                            8.705E-14, 6.303E-14, 4.538E-14,
     &                            3.254E-14, 2.326E-14, 1.651E-14,
     &                            1.156E-14, 7.899E-15, 5.195E-15,
     &                            3.255E-15, 1.943E-15, 1.126E-15,
     &                            6.528E-16, 3.829E-16, 2.231E-16,
     &                            1.263E-16/
      DATA (X_SEC(26,I),I=1,40) / 1.527E-11, 1.409E-11, 1.268E-11,
     &                            1.109E-11, 9.425E-12, 7.834E-12,
     &                            6.430E-12, 5.262E-12, 4.317E-12,
     &                            3.552E-12, 2.926E-12, 2.407E-12,
     &                            1.972E-12, 1.605E-12, 1.293E-12,
     &                            1.027E-12, 8.018E-13, 6.163E-13,
     &                            4.685E-13, 3.544E-13, 2.678E-13,
     &                            2.021E-13, 1.522E-13, 1.141E-13,
     &                            8.504E-14, 6.310E-14, 4.669E-14,
     &                            3.451E-14, 2.548E-14, 1.872E-14,
     &                            1.360E-14, 9.663E-15, 6.652E-15,
     &                            4.404E-15, 2.816E-15, 1.771E-15,
     &                            1.121E-15, 7.147E-16, 4.516E-16,
     &                            2.781E-16/
      DATA (X_SEC(27,I),I=1,40) / 1.251E-11, 1.155E-11, 1.040E-11,
     &                            9.103E-12, 7.749E-12, 6.452E-12,
     &                            5.307E-12, 4.354E-12, 3.582E-12,
     &                            2.955E-12, 2.443E-12, 2.017E-12,
     &                            1.659E-12, 1.356E-12, 1.097E-12,
     &                            8.766E-13, 6.894E-13, 5.343E-13,
     &                            4.101E-13, 3.137E-13, 2.399E-13,
     &                            1.836E-13, 1.404E-13, 1.071E-13,
     &                            8.142E-14, 6.178E-14, 4.686E-14,
     &                            3.558E-14, 2.704E-14, 2.049E-14,
     &                            1.538E-14, 1.133E-14, 8.136E-15,
     &                            5.667E-15, 3.854E-15, 2.601E-15,
     &                            1.770E-15, 1.210E-15, 8.193E-16,
     &                            5.423E-16/
      DATA (X_SEC(28,I),I=1,40) / 1.024E-11, 9.458E-12, 8.520E-12,
     &                            7.464E-12, 6.362E-12, 5.305E-12,
     &                            4.371E-12, 3.594E-12, 2.963E-12,
     &                            2.452E-12, 2.032E-12, 1.683E-12,
     &                            1.389E-12, 1.139E-12, 9.263E-13,
     &                            7.439E-13, 5.886E-13, 4.595E-13,
     &                            3.558E-13, 2.749E-13, 2.127E-13,
     &                            1.648E-13, 1.278E-13, 9.905E-14,
     &                            7.668E-14, 5.937E-14, 4.606E-14,
     &                            3.585E-14, 2.796E-14, 2.178E-14,
     &                            1.685E-14, 1.284E-14, 9.576E-15,
     &                            6.985E-15, 5.016E-15, 3.596E-15,
     &                            2.599E-15, 1.885E-15, 1.353E-15,
     &                            9.530E-16/
      DATA (X_SEC(29,I),I=1,40) / 8.374E-12, 7.737E-12, 6.974E-12,
     &                            6.114E-12, 5.216E-12, 4.355E-12,
     &                            3.595E-12, 2.961E-12, 2.446E-12,
     &                            2.028E-12, 1.685E-12, 1.399E-12,
     &                            1.158E-12, 9.536E-13, 7.784E-13,
     &                            6.280E-13, 4.997E-13, 3.927E-13,
     &                            3.065E-13, 2.390E-13, 1.868E-13,
     &                            1.465E-13, 1.151E-13, 9.051E-14,
     &                            7.125E-14, 5.622E-14, 4.452E-14,
     &                            3.543E-14, 2.830E-14, 2.262E-14,
     &                            1.797E-14, 1.412E-14, 1.091E-14,
     &                            8.296E-15, 6.249E-15, 4.716E-15,
     &                            3.586E-15, 2.732E-15, 2.063E-15,
     &                            1.532E-15/
      DATA (X_SEC(30,I),I=1,40) / 6.844E-12, 6.325E-12, 5.704E-12,
     &                            5.004E-12, 4.273E-12, 3.572E-12,
     &                            2.952E-12, 2.435E-12, 2.016E-12,
     &                            1.675E-12, 1.394E-12, 1.161E-12,
     &                            9.635E-13, 7.956E-13, 6.517E-13,
     &                            5.280E-13, 4.223E-13, 3.340E-13,
     &                            2.626E-13, 2.065E-13, 1.630E-13,
     &                            1.292E-13, 1.028E-13, 8.194E-14,
     &                            6.552E-14, 5.261E-14, 4.247E-14,
     &                            3.450E-14, 2.817E-14, 2.303E-14,
     &                            1.877E-14, 1.516E-14, 1.210E-14,
     &                            9.545E-15, 7.496E-15, 5.909E-15,
     &                            4.690E-15, 3.726E-15, 2.937E-15,
     &                            2.282E-15/
      DATA (X_SEC(31,I),I=1,40) / 5.592E-12, 5.169E-12, 4.663E-12,
     &                            4.093E-12, 3.498E-12, 2.927E-12,
     &                            2.422E-12, 2.001E-12, 1.659E-12,
     &                            1.381E-12, 1.152E-12, 9.610E-13,
     &                            7.997E-13, 6.623E-13, 5.443E-13,
     &                            4.427E-13, 3.558E-13, 2.830E-13,
     &                            2.241E-13, 1.776E-13, 1.415E-13,
     &                            1.133E-13, 9.119E-14, 7.369E-14,
     &                            5.980E-14, 4.881E-14, 4.011E-14,
     &                            3.321E-14, 2.766E-14, 2.311E-14,
     &                            1.927E-14, 1.596E-14, 1.311E-14,
     &                            1.069E-14, 8.701E-15, 7.116E-15,
     &                            5.855E-15, 4.820E-15, 3.939E-15,
     &                            3.180E-15/
      DATA (X_SEC(32,I),I=1,40) / 4.572E-12, 4.228E-12, 3.816E-12,
     &                            3.351E-12, 2.866E-12, 2.400E-12,
     &                            1.988E-12, 1.645E-12, 1.366E-12,
     &                            1.139E-12, 9.515E-13, 7.954E-13,
     &                            6.634E-13, 5.509E-13, 4.542E-13,
     &                            3.708E-13, 2.993E-13, 2.395E-13,
     &                            1.909E-13, 1.525E-13, 1.225E-13,
     &                            9.911E-14, 8.067E-14, 6.601E-14,
     &                            5.433E-14, 4.504E-14, 3.764E-14,
     &                            3.172E-14, 2.691E-14, 2.292E-14,
     &                            1.951E-14, 1.654E-14, 1.394E-14,
     &                            1.170E-14, 9.817E-15, 8.281E-15,
     &                            7.023E-15, 5.955E-15, 5.017E-15,
     &                            4.182E-15/
      DATA (X_SEC(33,I),I=1,40) / 3.750E-12, 3.468E-12, 3.131E-12,
     &                            2.751E-12, 2.354E-12, 1.973E-12,
     &                            1.636E-12, 1.355E-12, 1.127E-12,
     &                            9.408E-13, 7.876E-13, 6.596E-13,
     &                            5.514E-13, 4.590E-13, 3.796E-13,
     &                            3.110E-13, 2.522E-13, 2.028E-13,
     &                            1.627E-13, 1.310E-13, 1.062E-13,
     &                            8.672E-14, 7.135E-14, 5.910E-14,
     &                            4.930E-14, 4.147E-14, 3.520E-14,
     &                            3.015E-14, 2.602E-14, 2.256E-14,
     &                            1.957E-14, 1.693E-14, 1.459E-14,
     &                            1.255E-14, 1.081E-14, 9.355E-15,
     &                            8.133E-15, 7.068E-15, 6.107E-15,
     &                            5.228E-15/
      DATA (X_SEC(34,I),I=1,40) / 3.093E-12, 2.861E-12, 2.584E-12,
     &                            2.271E-12, 1.945E-12, 1.632E-12,
     &                            1.354E-12, 1.123E-12, 9.348E-13,
     &                            7.816E-13, 6.553E-13, 5.499E-13,
     &                            4.606E-13, 3.844E-13, 3.188E-13,
     &                            2.621E-13, 2.135E-13, 1.726E-13,
     &                            1.394E-13, 1.130E-13, 9.238E-14,
     &                            7.618E-14, 6.334E-14, 5.308E-14,
     &                            4.485E-14, 3.824E-14, 3.293E-14,
     &                            2.863E-14, 2.509E-14, 2.209E-14,
     &                            1.948E-14, 1.716E-14, 1.508E-14,
     &                            1.324E-14, 1.165E-14, 1.030E-14,
     &                            9.139E-15, 8.103E-15, 7.147E-15,
     &                            6.251E-15/
      DATA (X_SEC(35,I),I=1,40) / 2.576E-12, 2.383E-12, 2.153E-12,
     &                            1.893E-12, 1.622E-12, 1.362E-12,
     &                            1.131E-12, 9.391E-13, 7.827E-13,
     &                            6.553E-13, 5.503E-13, 4.626E-13,
     &                            3.883E-13, 3.248E-13, 2.701E-13,
     &                            2.229E-13, 1.823E-13, 1.482E-13,
     &                            1.204E-13, 9.834E-14, 8.104E-14,
     &                            6.743E-14, 5.664E-14, 4.799E-14,
     &                            4.103E-14, 3.543E-14, 3.091E-14,
     &                            2.724E-14, 2.419E-14, 2.160E-14,
     &                            1.932E-14, 1.728E-14, 1.544E-14,
     &                            1.379E-14, 1.235E-14, 1.110E-14,
     &                            1.001E-14, 9.021E-15, 8.088E-15,
     &                            7.197E-15/
      DATA (X_SEC(36,I),I=1,40) / 2.173E-12, 2.011E-12, 1.817E-12,
     &                            1.599E-12, 1.370E-12, 1.151E-12,
     &                            9.573E-13, 7.955E-13, 6.637E-13,
     &                            5.564E-13, 4.680E-13, 3.940E-13,
     &                            3.314E-13, 2.779E-13, 2.317E-13,
     &                            1.919E-13, 1.576E-13, 1.287E-13,
     &                            1.052E-13, 8.656E-14, 7.189E-14,
     &                            6.034E-14, 5.116E-14, 4.379E-14,
     &                            3.786E-14, 3.306E-14, 2.918E-14,
     &                            2.601E-14, 2.337E-14, 2.112E-14,
     &                            1.912E-14, 1.732E-14, 1.568E-14,
     &                            1.421E-14, 1.290E-14, 1.176E-14,
     &                            1.074E-14, 9.803E-15, 8.905E-15,
     &                            8.033E-15/
      DATA (X_SEC(37,I),I=1,40) / 1.862E-12, 1.724E-12, 1.558E-12,
     &                            1.371E-12, 1.176E-12, 9.886E-13,
     &                            8.228E-13, 6.843E-13, 5.716E-13,
     &                            4.798E-13, 4.041E-13, 3.408E-13,
     &                            2.872E-13, 2.413E-13, 2.018E-13,
     &                            1.676E-13, 1.382E-13, 1.135E-13,
     &                            9.328E-14, 7.724E-14, 6.462E-14,
     &                            5.468E-14, 4.676E-14, 4.040E-14,
     &                            3.527E-14, 3.111E-14, 2.774E-14,
     &                            2.498E-14, 2.267E-14, 2.068E-14,
     &                            1.892E-14, 1.732E-14, 1.585E-14,
     &                            1.452E-14, 1.333E-14, 1.228E-14,
     &                            1.134E-14, 1.045E-14, 9.591E-15,
     &                            8.744E-15/
      DATA (X_SEC(38,I),I=1,40) / 1.624E-12, 1.503E-12, 1.359E-12,
     &                            1.196E-12, 1.026E-12, 8.634E-13,
     &                            7.191E-13, 5.987E-13, 5.006E-13,
     &                            4.207E-13, 3.548E-13, 2.997E-13,
     &                            2.530E-13, 2.130E-13, 1.786E-13,
     &                            1.488E-13, 1.232E-13, 1.016E-13,
     &                            8.395E-14, 6.994E-14, 5.891E-14,
     &                            5.021E-14, 4.328E-14, 3.770E-14,
     &                            3.319E-14, 2.953E-14, 2.656E-14,
     &                            2.412E-14, 2.207E-14, 2.030E-14,
     &                            1.873E-14, 1.729E-14, 1.596E-14,
     &                            1.475E-14, 1.367E-14, 1.270E-14,
     &                            1.181E-14, 1.097E-14, 1.015E-14,
     &                            9.335E-15/
      DATA (X_SEC(39,I),I=1,40) / 1.440E-12, 1.333E-12, 1.206E-12,
     &                            1.062E-12, 9.115E-13, 7.671E-13,
     &                            6.394E-13, 5.327E-13, 4.459E-13,
     &                            3.751E-13, 3.168E-13, 2.679E-13,
     &                            2.266E-13, 1.912E-13, 1.606E-13,
     &                            1.342E-13, 1.115E-13, 9.232E-14,
     &                            7.668E-14, 6.423E-14, 5.444E-14,
     &                            4.670E-14, 4.053E-14, 3.556E-14,
     &                            3.153E-14, 2.827E-14, 2.561E-14,
     &                            2.341E-14, 2.157E-14, 1.998E-14,
     &                            1.855E-14, 1.724E-14, 1.604E-14,
     &                            1.493E-14, 1.392E-14, 1.302E-14,
     &                            1.219E-14, 1.139E-14, 1.061E-14,
     &                            9.818E-15/
      DATA (X_SEC(40,I),I=1,40) / 1.299E-12, 1.203E-12, 1.088E-12,
     &                            9.583E-13, 8.229E-13, 6.929E-13,
     &                            5.779E-13, 4.819E-13, 4.037E-13,
     &                            3.400E-13, 2.874E-13, 2.434E-13,
     &                            2.061E-13, 1.742E-13, 1.467E-13,
     &                            1.229E-13, 1.024E-13, 8.513E-14,
     &                            7.102E-14, 5.979E-14, 5.094E-14,
     &                            4.394E-14, 3.836E-14, 3.387E-14,
     &                            3.022E-14, 2.726E-14, 2.484E-14,
     &                            2.285E-14, 2.117E-14, 1.971E-14,
     &                            1.840E-14, 1.720E-14, 1.608E-14,
     &                            1.505E-14, 1.412E-14, 1.327E-14,
     &                            1.249E-14, 1.173E-14, 1.097E-14,
     &                            1.021E-14/
      DATA (X_NH(I),I=1,40) /     1.000E+19, 1.425E+19, 2.031E+19,
     &                            2.894E+19, 4.125E+19, 5.878E+19,
     &                            8.377E+19, 1.194E+20, 1.701E+20,
     &                            2.424E+20, 3.455E+20, 4.924E+20,
     &                            7.017E+20, 1.000E+21, 1.425E+21,
     &                            2.031E+21, 2.894E+21, 4.125E+21,
     &                            5.878E+21, 8.377E+21, 1.194E+22,
     &                            1.701E+22, 2.424E+22, 3.455E+22,
     &                            4.924E+22, 7.017E+22, 1.000E+23,
     &                            1.425E+23, 2.031E+23, 2.894E+23,
     &                            4.125E+23, 5.878E+23, 8.377E+23,
     &                            1.194E+24, 1.701E+24, 2.424E+24,
     &                            3.455E+24, 4.924E+24, 7.017E+24,
     &                            1.000E+25/
      DATA (X_TX(I),I=1,40) /     3.162E+05, 3.888E+05, 4.781E+05,
     &                            5.878E+05, 7.227E+05, 8.886E+05,
     &                            1.093E+06, 1.343E+06, 1.652E+06,
     &                            2.031E+06, 2.497E+06, 3.070E+06,
     &                            3.775E+06, 4.642E+06, 5.707E+06,
     &                            7.017E+06, 8.628E+06, 1.061E+07,
     &                            1.304E+07, 1.604E+07, 1.972E+07,
     &                            2.424E+07, 2.981E+07, 3.665E+07,
     &                            4.507E+07, 5.541E+07, 6.813E+07,
     &                            8.377E+07, 1.030E+08, 1.266E+08,
     &                            1.557E+08, 1.914E+08, 2.354E+08,
     &                            2.894E+08, 3.559E+08, 4.375E+08,
     &                            5.380E+08, 6.615E+08, 8.133E+08,
     &                            1.000E+09/
      END
