      PROGRAM CONVERT_GRID
c
c                                            Version November 26 2008
c
c                            Simon Bruderer
c                            Institute of Astronomy 
c                            ETH Zurich
c                            8093 Zurich, Switzerland
c        
c                            simonbr-at-astro.phys.ethz.ch
c
c
c     This program reads a ASCII grid file and converts them into binary form
c       
c     The format of the binary file is:
c
c     1 * CHARACTER*80 --> Header
c         CHARACTER*9  --> Name of species
c         INTEGER      --> n, Number of Parameters (NH2, T, FX, G0,...)
c
c     n times:    
c         CHARACTER*9  --> Name of the Parameter
c         INTEGER      --> m, Number of different values for this parameter
c         m * REAL     --> Parameter values
c
c         INTEGER      --> q, Number of timesteps
c         q* REAL      --> Time of gridpoints
c
c     For each model (The first parameter is the slowest running index)
c         q* REAL      --> Abundance of the molecule
c
c
c     The start-parameters for this script are:
c     ./convert [inputfile] [outputfile] [version of the grid]
c
       
      IMPLICIT NONE
       
      CHARACTER*150 File_In, File_Out
      CHARACTER*80 Oneline
      CHARACTER*9 Name_Molec, Name_Param
      REAL Read_Real, Abund(40),Dummies(15)       
      INTEGER Read_Int, Nr_Time, Nr_Param, Nr_Val, I, J, K, Nr_Mod
      
      LOGICAL PARAM
      
c change this parameter to .FALSE. if you want to convert files without 
c using the shell script convert.sh
            
      PARAMETER(PARAM=.TRUE.)
      
      IF(PARAM) THEN
        CALL GETARG(1, FILE_IN)
        CALL GETARG(2, FILE_OUT)               
      ELSE
        FILE_IN="mol_full_HCO+.dat"
        FILE_OUT="mol_full_HCO+.bin"	    
      END IF
            	   	          
       
      WRITE(*,*) 'Convert: ',FILE_IN,' --> ',FILE_OUT
	   
      OPEN(20,File=File_In,Status='OLD')   
       
      CALL System('rm -f '// File_Out)
               
      OPEN(30,File=File_Out,Status='NEW',FORM='UNFORMATTED')
       
      READ(20,*) Oneline ! Header
      READ(20,*) Name_Molec ! Name of Species
      READ(20,*) Oneline ! Header
      READ(20,*) Nr_Param
      READ(20,*) Oneline ! Header     
                  
      IF(PARAM) THEN
        CALL GETARG(3,ONELINE)
      ELSE
        ONELINE='GRIDMOL_R0.01'      
      END IF
		  
      WRITE(30) Oneline    

      WRITE(30) Name_Molec
      WRITE(30) Nr_Param    
       
      Nr_Mod=1
                    
      DO 100 I=1,Nr_Param
        READ(20,*) Name_Param
	READ(20,*) Nr_Val
	WRITE(30) Name_Param
	WRITE(30) Nr_Val 	  
	DO 200 J=1,Nr_Val
	  READ(20,*) Read_Real	     	      
c	  IF(I.EQ.3) Read_Real=10.0**Read_Real	      
c         IF(I.EQ.4) Read_Real=10.0**Read_Real	      	      	      	      
	  WRITE(30) Read_Real	      
200     CONTINUE        

c     Calculate the number of models by Nr_Param1 * Nr_Param2 ...

      Nr_Mod=Nr_Mod*Nr_Val          
100   CONTINUE
              
      READ(20,*) Oneline 
      READ(20,*) Nr_Time
       
      WRITE(30) Nr_Time
              
      DO 300 I=1,Nr_Time
        READ(20,*) Read_Real
	WRITE(30) Read_Real	             
300   CONTINUE       

      READ(20,*) Oneline 

      DO 400 I=1,Nr_Mod
        READ(20,*) (Dummies(J),J=1,Nr_Param), (Abund(K),K=1,Nr_Time)
        WRITE(30) (Abund(K),K=1,Nr_Time)             
400   CONTINUE

      WRITE(*,*) 'Converted! Nr_Models=',Nr_Mod

      CLOSE(20)
      CLOSE(30)       
      END
