C
C     FULLY-IMPLICIT DECOUPLING METHOD OF  
C     THE INCOMPRESSIBLE NAVIER-STOKES EQUATIONS
C
C     DNS/LES OF TURBULENT BOUNDARY LAYER 
C     WITH INFLOW GENERATION ROUTINE FOR SPATIALLY-DEVELOPING
C     BOUNDARY LAYER SIMULATION
C
C        
C                                                   KYOUNGYOUN KIM
C                                          FLOW CONTROL LABORATORY
C                             DEPARTMENT OF MECHANICAL ENGINEERING
C                                                            KAIST 


c     In order to syncronize the phase, ntime_ctd = ntime+ntime_tr is introduced.
c     Subroutine WRITEUP2 is added, in which instantaneous flow field is written at moni_ph.
c     time_avg_routine is added. in 'param.h', MT should be assined as 0        04/06/15  
c     OMP is implemented, however, TRIDP routine is excluded due to unexpected error. 04/08/15

      PROGRAM MAIN

      INCLUDE 'PARAM.H'
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/FINOUT/IDTOPT,NWRITE,NREAD,IAVG,IAVG_TRIG,NPRN,INSF,NINS
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/LES_OPT/ILES
      COMMON/CALC/ICONT
      common/phaseavg/nphavg(0:mt),nphavg_les(0:mt),ntime_ctd

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL EV(0:M1,0:M2,0:M3)   
      REAL EVI(0:M1,0:M2,0:M3)   
      REAL EVJ(0:M1,0:M2,0:M3)   
      REAL EVK(0:M1,0:M2,0:M3)   

      REAL*4 T1_check,T1_start
      

      nphavg(:) = 0
      nphavg_les(:) = 0

      CALL SETUP

      CALL INIUP(Q1,Q2,Q3,P,EV)
      IF(NREAD.EQ.1) CALL READUP(Q1,Q2,Q3,P)

      ntst_tr=0
c      ntst_tr=2000
      time=ntst_tr*dt  ! for the same dt 
    
      CALL DIVCHECK(Q1,Q2,Q3,DIVMAX)
      CALL CFL(Q1,Q2,Q3,CFLM)
      WRITE(*,*)
      WRITE (*,100) DIVMAX,CFLM*DT

      IMORE=0                ! INDEX OF WRITTEN FILE
      NAVG=0                
      DTR=DT

C---------------------------------------------------------
C     TO CHECK CPU TIME DURING THE OVERALL COMPUTATION
      T1_start = 0.0
      T1_check=SECNDS(T1_start)
C---------------------------------------------------------

      DO 10 NTIME=1,NTST

         ntime_ctd=ntime+ntst_tr

      IF (ICONT.EQ.0.AND.NTIME.LE.20) THEN
c pass for main
          DT=DTR/20.   ! FOR STABLE START WITH RANDOM NOISE FIELD
      ELSE  
          CALL CFL(Q1,Q2,Q3,CFLM)
          IF (CFLM*DTR .GE.CFLMAX.AND.IDTOPT.EQ.1) THEN
                 DT=CFLMAX/CFLM
          ELSE
                 DT=DTR
          ENDIF

      ENDIF

      TIME=TIME+DT

           TRANSIENT_TIME=0.
      IF(TIME.LT.TRANSIENT_TIME) IAVG_TRIG=0
      IF(IAVG.EQ.1.AND.TIME.GE.TRANSIENT_TIME) IAVG_TRIG=1

      WRITE(*,*)
      WRITE(*,99) 
      WRITE(*,101) NTIME,TIME
      WRITE(*,99) 
      WRITE(*,*)

C      IF (ILES.EQ.1) CALL EDDY_VISCOS_SM(EV,Q1,Q2,Q3)
      IF (ILES.EQ.1) CALL EDDY_VISCOS_DM(EV,Q1,Q2,Q3)
      CALL EV_INTERPOLATION(EV,EVI,EVJ,EVK)

      CALL GETUP(Q1,Q2,Q3,P,TIME,EV,EVI,EVJ,EVK)
      CALL THIST(TIME,Q1,Q2,Q3,P)
      CALL DIVCHECK(Q1,Q2,Q3,DIVMAX)
      CALL CFL(Q1,Q2,Q3,CFLM)

C      CALL SAVE_INFLOW(Q1,Q2,Q3)

      WRITE(*,*)
      WRITE (*,100) DIVMAX,CFLM*DT
      WRITE(*,*)

c      call write_time_signal(q1,q2,q3,p)

      IF(IAVG_TRIG.EQ.1.AND.MOD(NTIME,1).EQ.0) THEN       ! AVERAGING PER 1 TIME STEP
            if (MT.eq.0) then 
                 call time_avg(q1,q2,q3,p)
            else 
                 call phase_avg(q1,q2,q3,p)
            endif
      ENDIF

      IF(MOD(NTIME,NPRN).EQ.0.AND.NWRITE.EQ.1) THEN
      CALL WRITEUP(Q1,Q2,Q3,P,IMORE)
      ENDIF

c      num_moni_ph=16
c      iphase = mod(ntime_ctd,mt+1)
c      if( mod(iphase,int((mt+1)/num_moni_ph)).eq.0
c     &  ) then
c         CALL WRITEUP2(Q1,Q2,Q3,P)
c         write(*,*) iphase, q2(203,1,n3m/2)
c      endif


      IF(INSF.EQ.1.AND.MOD(NTIME,NINS).EQ.0) CALL INSFIELD(Q1,Q2,Q3,P)
      IF(INSF.EQ.1.AND.MOD(NTIME,NINS).EQ.0) CALL INS_2D(Q1,Q2,Q3,P)
 
 10   CONTINUE

C---------------------------------------------------------
C     TO CHECK CPU TIME DURING THE OVERALL COMPUTATION
      T1_start=SECNDS(T1_check)
      WRITE(*,*) 'OVERALL CPU TIME : ',T1_start
C---------------------------------------------------------

99    FORMAT (80('='))
100   FORMAT('MAXIMUM DIVERGENCE=',E12.5,X,'MAXIMUM CFL NUMBER=',E12.5)
101   FORMAT( I6,'TH STEP   TIME : ',E12.5) 

      STOP
      END


c--------------- write_time_signal ------------------------
      subroutine write_time_signal(q1,q2,q3,p)
      INCLUDE 'PARAM.H'

      parameter(m_probe=99)
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/LES_OPT/ILES

      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)

      common/moni_spec/mij(m_probe,2),len_signal  ! (number of probe, ij index at which probe is located.)

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      real xy_moni(m_probe,2)    ! e.g. xy_moni(1,2)-> 1st probe , y location

      real q1c(0:m3)
      real q2c(0:m3)
      real q3c(0:m3)

      character*30 file1

      if (ntime.eq.1) then

c-------------------------------
         len_signal=2000

         num_probe=42

         do jj=1,6
            xy_moni(jj   ,1)=10.
            xy_moni(jj+6 ,1)=50.
            xy_moni(jj+12,1)=109.
            xy_moni(jj+18,1)=114.
            xy_moni(jj+24,1)=123.
            xy_moni(jj+30,1)=133.
            xy_moni(jj+36,1)=160.
         enddo


         do ii=1,7
            utau_ii=-0.737321E-5*xy_moni((ii-1)*6+1,1)+0.044131
            xy_moni((ii-1)*6+1,2)=0.5*(y(1)+y(2))  !wall
            xy_moni((ii-1)*6+2,2)=15. /1410./utau_ii
            xy_moni((ii-1)*6+3,2)=40. /1410./utau_ii
            xy_moni((ii-1)*6+4,2)=60. /1410./utau_ii
            xy_moni((ii-1)*6+5,2)=150./1410./utau_ii
            xy_moni((ii-1)*6+6,2)=300./1410./utau_ii
         enddo

cc 1st probe
c         xy_moni(1,1)= 50.0
c         xy_moni(1,2)=  1.0
cc 2nd probe
c         xy_moni(2,1)= 50.0
c         xy_moni(2,2)=  3.0
cc 3rd probe
c         xy_moni(3,1)= 50.0
c         xy_moni(3,2)=  0.5*(y(1)+y(2))  ! wall
cc 4th probe
c         xy_moni(4,1)=109.0
c         xy_moni(4,2)=  1.0
cc 5th probe
c         xy_moni(5,1)=109.0
c         xy_moni(5,2)=  3.0
cc 6th probe
c         xy_moni(6,1)=109.0
c         xy_moni(6,2)=  0.5*(y(1)+y(2))  ! wall
c-------------------------------


      do 1 mm=1,num_probe

        do i=1,n1m-1
         temp=(xy_moni(mm,1)-x(i  ))
     &       *(xy_moni(mm,1)-x(i+1))
         if(temp.le.0.0) then
           mij(mm,1)=i
           goto 2
         endif
        enddo

2       continue

        do j=1,n2m-1
         temp=(xy_moni(mm,2)-y(j  ))
     &       *(xy_moni(mm,2)-y(j+1))
         if(temp.le.0.0) then
           mij(mm,2)=j
           goto 1
         endif
        enddo

1     continue

      do k=1,num_probe
       write(*,*) k
     &           ,xy_moni(k,1)
     &           ,mij(k,1)
     &          ,0.5*(x(mij(k,1))+x(mij(k,1)+1))
     &           ,xy_moni(k,2)
     &           ,mij(k,2)
     &          ,0.5*(y(mij(k,2))+y(mij(k,2)+1))
      enddo

      endif  ! (ntime.eq.1)
c---------------------------------------------



      do mm=1,num_probe

        file1='probe'
        nn=index(file1,'e')
        write(unit=file1(nn+1:),fmt='(bn,i2.2,a1)') mm,'.'

        nn=index(file1,'.')
        write(unit=file1(nn+1:),fmt='(bn,i3.3)')
     &       (ntime-1)/len_signal

        open(100,file=file1,status='unknown',action='write'
     &          ,form='unformatted',position='append')


c--- write signals at cell center

      i=mij(mm,1)
      j=mij(mm,2)

c--- in order to check post code for spectrum
c      do k=1,n3m
c       pi=acos(-1.0)
c       q1(i,j,k)=sin(3. *real(ntime)*0.2)
c     &          *sin(10.*real(k-0.5)/dx3)
c      enddo

      do k=1,n3m
        q1c(k)=0.5*(q1(i,j,k)+q1(i+1,j,k))
        q2c(k)=0.5*(q2(i,j,k)+q2(i,j+1,k))
        q3c(k)=0.5*(q3(i,j,k)+q3(i,j,kpa(k)))
      enddo

      write(100) (q1c(k),k=1,n3m)
      write(100) (q2c(k),k=1,n3m)
      write(100) (q3c(k),k=1,n3m)
      write(100) (p(i,j,k),k=1,n3m)

      close(100)

      enddo  ! mm


      return
      end


C***************** SETUP ***********************     
      SUBROUTINE SETUP
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/VPERIN/VPER
      COMMON/SIZE/ALX,ALY,ALZ
      COMMON/FINOUT/IDTOPT,NWRITE,NREAD,IAVG,IAVG_TRIG,NPRN,INSF,NINS
      COMMON/FILENAME/FILEINI,FILEOUT,FILEAVG
      COMMON/CALC/ICONT
      COMMON/LES_OPT/ILES
C---
      COMMON/MLV/MLEV,IWMG,MLW,MLCMG,MHCMG,EPSM,RATIO,NSOR
      COMMON/POIS/IK0,IPH0
      COMMON/MAXIT/NLMAXI,NHMAXI,OMEG
C---
      COMMON/GRD_FILE/FILE_X_GRD,FILE_Y_GRD

      CHARACTER*20 FILE_X_GRD,FILE_Y_GRD
      CHARACTER*20 FILEINI,FILEOUT,FILEAVG

      CHARACTER*65 DUMMY

      OPEN(1,FILE='tbl.set',STATUS='OLD')
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,301) DUMMY,ILES
      WRITE(*,301) DUMMY,ILES
      READ (1,301) DUMMY,N1
      WRITE(*,301) DUMMY,N1
      READ (1,301) DUMMY,N2
      WRITE(*,301) DUMMY,N2
      READ (1,301) DUMMY,N3
      WRITE(*,301) DUMMY,N3
      READ (1,302) DUMMY,RE
      WRITE(*,302) DUMMY,RE
      READ (1,302) DUMMY,ALX
      WRITE(*,302) DUMMY,ALX
      READ (1,302) DUMMY,ALY
      WRITE(*,302) DUMMY,ALY
      READ (1,302) DUMMY,ALZ
      WRITE(*,302) DUMMY,ALZ
      READ (1,301) DUMMY,NTST
      WRITE(*,301) DUMMY,NTST
      READ (1,302) DUMMY,VPER
      WRITE(*,302) DUMMY,VPER
      READ (1,302) DUMMY,DT
      WRITE(*,302) DUMMY,DT
      READ (1,301) DUMMY,IDTOPT
      WRITE(*,301) DUMMY,IDTOPT
      READ (1,301) DUMMY,ICONT
      WRITE(*,301) DUMMY,ICONT
      READ (1,302) DUMMY,CFLMAX
      WRITE(*,302) DUMMY,CFLMAX
      READ (1,301) DUMMY,NWRITE
      WRITE(*,301) DUMMY,NWRITE
      READ (1,301) DUMMY,NREAD
      WRITE(*,301) DUMMY,NREAD
      READ (1,301) DUMMY,IAVG
      WRITE(*,301) DUMMY,IAVG
      READ (1,301) DUMMY,NPRN
      WRITE(*,301) DUMMY,NPRN
      READ (1,301) DUMMY,INSF
      WRITE(*,301) DUMMY,INSF
      READ (1,301) DUMMY,NINS
      WRITE(*,301) DUMMY,NINS
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,301) DUMMY,MLEV
      WRITE(*,301) DUMMY,MLEV
      READ (1,301) DUMMY,IWMG
      WRITE(*,301) DUMMY,IWMG
      READ (1,301) DUMMY,MLW
      WRITE(*,301) DUMMY,MLW
      READ (1,301) DUMMY,MHCMG
      WRITE(*,301) DUMMY,MHCMG
      READ (1,301) DUMMY,MLCMG
      WRITE(*,301) DUMMY,MLCMG
      READ (1,301) DUMMY,NHMAXI
      WRITE(*,301) DUMMY,NHMAXI
      READ (1,301) DUMMY,NLMAXI
      WRITE(*,301) DUMMY,NLMAXI
      READ (1,302) DUMMY,OMEG
      WRITE(*,302) DUMMY,OMEG
      READ (1,302) DUMMY,EPSM
      WRITE(*,302) DUMMY,EPSM
      READ (1,302) DUMMY,RATIO
      WRITE(*,302) DUMMY,RATIO
      READ (1,301) DUMMY,NSOR
      WRITE(*,301) DUMMY,NSOR
      READ (1,301) DUMMY,IK0
      WRITE(*,301) DUMMY,IK0
      READ (1,301) DUMMY,IPH0
      WRITE(*,301) DUMMY,IPH0
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,303) DUMMY,FILEINI
      WRITE(*,303) DUMMY,FILEINI
      READ (1,303) DUMMY,FILE_X_GRD
      WRITE(*,303) DUMMY,FILE_X_GRD
      READ (1,303) DUMMY,FILE_Y_GRD
      WRITE(*,303) DUMMY,FILE_Y_GRD
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY

      CLOSE(1)

  300 FORMAT(A65)
  301 FORMAT(A45,I15)
  302 FORMAT(A45,E15.7)
  303 FORMAT(A45,A20)

      N1M = N1 - 1
      N2M = N2 - 1
      N3M = N3 - 1

      CALL MESH
      CALL INDICES 
      CALL INIPOISSON

      RETURN
      END

C***************** MESH ***********************     
      SUBROUTINE MESH
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/SIZE/ALX,ALY,ALZ
      COMMON/FILENAME/FILEINI,FILEOUT,FILEAVG
      COMMON/GRD_FILE/FILE_X_GRD,FILE_Y_GRD

      CHARACTER*20 FILE_X_GRD,FILE_Y_GRD
      CHARACTER*20 FILEINI,FILEOUT,FILEAVG
 
      PI=ACOS(-1.0)

C     CREATE THE GRID IN X1 DIRECTION      
         X(0)=0.0 
      DO I=1,N1 
         X(I)=REAL(I-1)*ALX/REAL(N1M)
C         X(I)=0.5*ALX*(1.0-COS(PI*REAL(I-1)/REAL(N1M)))
C         X(I)=ALX*(1.0-COS(0.5*PI*REAL(I-1)/REAL(N1M)))
      ENDDO

c         X(0)=0.0 
C      OPEN(10,FILE='GRID_X200.GRD',STATUS='UNKNOWN')
c      OPEN(10,FILE=FILE_X_GRD,STATUS='UNKNOWN')
c      DO I=1,N1
c      READ(10,*) DUMMY,X(I),DUMMY
c      ENDDO
c      CLOSE(10)

         DX(0)=0.0
      DO I=1,N1M
         DX(I)=X(I+1)-X(I)
      ENDDO 
         DX(N1)=0.0

cc---  for consistent with uniform grid at inlet and outlet
c         DX(0 ) =DX(1)
c         DX(N1) =DX(1)


      DO I=1,N1 
         G(I)=0.5*(DX(I)+DX(I-1))
      ENDDO

C     CREATE THE GRID IN X2 DIRECTION      

C         Y(0)=0.0 
C      DO J=1,N2 
C         Y(J)=REAL(J-1)*ALY/REAL(N2M)
C         Y(J)=0.5*ALY*(1.0-COS(PI*REAL(J-1)/REAL(N2M)))
C         Y(J)=ALY*(1.0-COS(0.5*PI*REAL(J-1)/REAL(N2M)))
C         Y(J)=ALY*(PI*REAL(J-1)/REAL(N2M)-SIN(PI*REAL(J-1)/REAL(N2M)))
C      ENDDO


C     CREATE THE TANGENT HYPERBOLIC GRID IN X2 DIRECTION
c      CS=0.048
c      CS=0.05
      CS=0.0550
      A=TANH(CS*REAL(N2M))

         Y(0)=0.0
      DO J=1,N2
         Y(J)=ALY*(A-TANH(CS*REAL(N2-J)))/A
      ENDDO


c         Y(0)=0.0
c      OPEN(10,FILE=FILE_Y_GRD,STATUS='UNKNOWN')
c      DO J=1,N2
c      READ(10,*) DUMMY,Y(J),DUMMY
c      ENDDO
c      CLOSE(10)
c      ALY=Y(N2)


         DY(0)=0.0
      DO J=1,N2M
         DY(J)=Y(J+1)-Y(J)
      ENDDO 
         DY(N2)=0.0

      DO J=1,N2 
         H(J)=0.5*(DY(J)+DY(J-1))
      ENDDO

       
      DX3=DBLE(N3M)/ALZ
      DX3Q=DX3**2.0


      OPEN(10,FILE='GRID.PLT',STATUS='UNKNOWN')
      DO J=1,N2
      WRITE(10,*) J,Y(J),DY(J)
      ENDDO
      CLOSE(10)

      OPEN(10,FILE='GRID_X.PLT',STATUS='UNKNOWN')
      DO I=1,N1
      WRITE(10,*) I,X(I),DX(I)
      ENDDO
      CLOSE(10)


      RETURN
      END

C***************** INDICES ***********************     
      SUBROUTINE INDICES
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      DO 10 IC=1,N1M
      IPA(IC)=IC+1
      IMU(IC)=IC-1
 10   IMV(IC)=IC-1
      IPA(N1M)=N1M
      IMU(2)=2
      IMV(1)=1

      DO 20 KC=1,N3M
      KPA(KC)=KC+1
 20   KMA(KC)=KC-1
      KPA(N3M)=1
      KMA(1)=N3M

      DO 30 JC=1,N2M
      JPA(JC)=JC+1
      JMU(JC)=JC-1
 30   JMV(JC)=JC-1
      JPA(N2M)=N2M
      JMU(1)=1
      JMV(2)=2

      RETURN
      END


C*************** INIUP **********************
      SUBROUTINE INIUP(Q1,Q2,Q3,P,EV)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M

      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/VPERIN/VPER
      COMMON/PARA/RE

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)
      REAL EV(0:M1,0:M2,0:M3)   

      REAL UM(0:M1,0:M2,3)

      DO 10 K=0,N3
      DO 10 J=0,N2
      DO 10 I=0,N1
        Q1(I,J,K)=0.0
        Q2(I,J,K)=0.0
        Q3(I,J,K)=0.0
         P(I,J,K)=0.0
        EV(I,J,K)=0.0
10    CONTINUE

c      THETA=1.0    ! FOR MAIN SIMULATION
c      THETA=270./RE ! FOR INFLOW, WE SHOULD SPECIFY THE INLET RE    
      THETA=620./RE ! FOR INFLOW, WE SHOULD SPECIFY THE INLET RE    

C---- CALCULATE U MEAN 

      DO 20 I=1,N1
      IF (I.GT.1) THETA=THETA+(1./RAMDA**2)*DX(MIN(I,N1M))   ! FROM EQ.(35) OF LUND ET AL.
      RE_TH=RE*THETA
      CALL GET_RAMDA(RE_TH,RAMDA)
C      WRITE(*,*) RE_TH,RAMDA
      CALL GET_UM(UM,RAMDA,I)
20    CONTINUE

      DO I=1,N1
       UM(I,0,1)=0.0
       UM(I,N2,1)=UM(I,N2M,1)
      ENDDO

C---- CALCULATE V MEAN 

      DO 40 I=1,N1M
      UM(I,1,2)=0.0   ! AT THE WALL
      DO 41 J=2,N2M
      UM(I,J,2)=UM(I,J-1,2)-DY(J-1)/DX(I)
     >         *(UM(I+1,J-1,1)-UM(I,J-1,1)) 
41    CONTINUE
      UM(I,N2,2)=UM(I,N2M,2)
40    CONTINUE

      DO 42 J=1,N2
      UM(0,J,2)=UM(1,J,2)
42    UM(N1,J,2)=UM(N1M,J,2)

C---- IMPOSE W MEAN AS ZERO

      DO 50 J=0,N2
      DO 50 I=0,N1
50    UM(I,J,3)=0.0


C---  RANDOM FLUCTUATIONS

c      CALL RNSET(1234)

      DO K=1,N3M
      DO J=1,N2M*3/4
      DO I=0,N1
c-- absoft
c       Q1(I,J,K)=VPER*DRNUNF()
c       Q2(I,J,K)=VPER*DRNUNF()
c       Q3(I,J,K)=VPER*DRNUNF()
c--- hp
       Q1(I,J,K)=VPER*rand()
       Q2(I,J,K)=VPER*rand()
       Q3(I,J,K)=VPER*rand()

      ENDDO
      ENDDO
      ENDDO

C---  FOR CORRECT U_TAU BASED ON THE SPANWISE AVERAGED MEAN VELOCITY
C---  REMOVE THE SPANMEAN OF RANDOM FLUCTUATIONS

      DO 71 J=1,N2
      DO 71 I=0,N1

         U_Z=0.0
         V_Z=0.0
         W_Z=0.0

      DO K=1,N3M
         U_Z=U_Z+Q1(I,J,K)/REAL(N3M)    
         V_Z=V_Z+Q2(I,J,K)/REAL(N3M)    
         W_Z=W_Z+Q3(I,J,K)/REAL(N3M)    
      ENDDO

      DO K=1,N3M
         Q1(I,J,K)=Q1(I,J,K)-U_Z  
         Q2(I,J,K)=Q2(I,J,K)-V_Z  
         Q3(I,J,K)=Q3(I,J,K)-W_Z  
      ENDDO

71    CONTINUE

      DO 140 K=1,N3M
      DO 140 J=1,N2
      DO 140 I=0,N1
        Q1(I,J,K)=UM(I,J,1)+Q1(I,J,K)
        Q2(I,J,K)=UM(I,J,2)+Q2(I,J,K)
        Q3(I,J,K)=UM(I,J,3)+Q3(I,J,K)
140   CONTINUE


      OPEN(20,FILE='INIUP.PLT',STATUS='UNKNOWN')
      WRITE(20,*) 'VARIABLES=X,Y,U,V,W'
      WRITE(20,*) 'ZONE I=',N1,', J=',N2,', F=POINT'      
      DO J=1,N2 
      DO I=1,N1
      WRITE(20,100) X(I),Y(J),Q1(I,J,N3M),Q2(I,J,N3M),Q3(I,J,N3M) 
100   FORMAT(5(E12.5,X))
      ENDDO
      ENDDO
      CLOSE(20)

      RETURN
      END

C*************** GET_RAMDA ******************
      SUBROUTINE GET_RAMDA(RE_TH,RAMDA)
      COMMON/TBL_PARA/C_KAPA,B,PII

      R1=18.0    !INITIAL GUESS FOR RAMDA
      R2=19.0

1     CONTINUE

      F_1=RE_THETA(R1)
      F_2=RE_THETA(R2)

      ERR=MAX(ABS(F_1-RE_TH),ABS(F_2-RE_TH))

      IF (ERR.GE.1E-10) THEN
      TEMP=(R1-R2)/(F_1-F_2)*(RE_TH-F_1)+R1
      R2=R1
      R1=TEMP
      GOTO 1
      ENDIF

      RAMDA=R1

      RETURN
      END 


C****************** RE_THETA(RAMDA) ***********************
      REAL FUNCTION RE_THETA(RAMDA)
      COMMON/TBL_PARA/C_KAPA,B,PII

      C_KAPA=0.41
      B=5.0
      PII=0.5 ! PARAMETER THAT DEPENDS ON THE PRESSURE GRADIENT

      RE_D=RAMDA*EXP(C_KAPA*(RAMDA-B)-2.*PII) 
      RE_THETA=((1+PII)/C_KAPA/RAMDA-1./C_KAPA**2/RAMDA**2*
     >        (2.+2.*PII*1.58951+
     >             1.5*PII**2))*RE_D
      END



C*************** GET_UM ******************************** 
      SUBROUTINE GET_UM(UM,RAMDA,IC)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/PARA/RE

      COMMON/TBL_PARA/C_KAPA,B,PII

      REAL UM(0:M1,0:M2,3)

      PI=ACOS(-1.0)
      U_T=1.0/RAMDA 

      TEMP=0.0 

      DO 10 J=1,N2M

C     INNER PROFILE
      YP=RE/RAMDA*(Y(J)+Y(J+1))*0.5
      U1=TEMP
      U2=RAMDA
100   CONTINUE

      Y1=SPALDING(U1)
      Y2=SPALDING(U2)

      ERR=ABS(Y1-SPALDING(U2))

      IF (ERR.GE.1E-8) THEN
      TEMP=(U2-U1)/(Y2-Y1)*(YP-Y1)+U1 
      U1=U2
      U2=TEMP
      GOTO 100
      ENDIF

      UM_INNER=TEMP/RAMDA

C     OUTER PROFILE

C---  COLE'S LAW OF THE WAKE 
C     'VISCOUS FLUID FLOW,' WHITE,  EQ.(6-47)

      DELTA=RAMDA*EXP(C_KAPA*(RAMDA-B)-2.*PII)/RE
      ETHA=(Y(J)+Y(J+1))*0.5/DELTA
      YP=RE/RAMDA*(Y(J)+Y(J+1))*0.5
      UM_OUTER=(1.0/C_KAPA*LOG(YP)+B+2.*PII/C_KAPA*(SIN(PI/2.*ETHA))**2)
     >        /RAMDA                                               
      IF (ETHA.GE.1.) UM_OUTER=1.0

C     WEIGHTED AVERAGE THE INNER AND OUTER PROFILES
      ALPHA=4.0
      BETA=0.2
      W=0.5*(1.+TANH(ALPHA*(ETHA-BETA)/((1.-2*BETA)*ETHA+BETA))
     >       /TANH(ALPHA))

      IF (ETHA.GE.1.) W=1.0

      UM(IC,J,1)=UM_INNER*(1.-W)+UM_OUTER*W

10    CONTINUE

      RETURN
      END 



C****************** SPALDING(UPLUS) ***********************
      REAL FUNCTION SPALDING(UPLUS)
C---                VISCOUS FLUID FLOW, WHITE EQ.(6-41)

      CK=0.40 
      BB=5.5
      SPALDING=UPLUS
     >        +EXP(-CK*BB)* 
     >         (EXP(CK*UPLUS)
     >         - 1.0  
     >         - CK*UPLUS
     >         -(CK*UPLUS)**2/2.0
     >         -(CK*UPLUS)**3/6.0)

      END



C  ************************  READUP **********************
C     READ FLOW FIELD AND BOUNDARY CONDITIONS
C     AND MEAN PRESSURE GRADIENT

      SUBROUTINE READUP(Q1,Q2,Q3,P)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/FILENAME/FILEINI,FILEOUT,FILEAVG
      CHARACTER*20 FILEINI,FILEOUT,FILEAVG

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      OPEN(3,FILE=FILEINI,STATUS='OLD',FORM='UNFORMATTED')
      READ(3) (((Q1(I,J,K),Q2(I,J,K),Q3(I,J,K)
     >            ,K=1,N3M),J=0,N2),I=0,N1)
      READ(3) (((P(I,J,K),K=1,N3M),J=1,N2M),I=1,N1M)
      CLOSE(3)

C---  CHECK THE PRSSSURE CONTRAINT
C---  (1) MULTI GRID SOLVER
c         TEMP=0.0
c      DO K=1,N3M
c         TEMP=TEMP+P(1,N2M,K)/REAL(N3M)
c      ENDDO 

C---  (2) COS TRANSFORM SOLVER
         TEMP=0.0
       DO K=1,N3M
       DO I=1,N1M
          TEMP=TEMP+P(I,N2M,K)/REAL(N3M)/REAL(N1M)
       ENDDO 
       ENDDO 
      
	DO K=1,N3M
	DO J=1,N2M
	DO I=1,N1M
         P(I,J,K)=P(I,J,K)-TEMP
	ENDDO
	ENDDO
	ENDDO

      RETURN
      END


C  ************************  WRITEUP **********************
C     WRITE FLOW FIELD AND BOUNDARY CONDITIONS
C     AND MEAN PRESSURE GRADIENT

      SUBROUTINE WRITEUP(Q1,Q2,Q3,P,IMORE)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/LES_OPT/ILES

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      NFILE=40+IMORE
      OPEN(NFILE,STATUS='UNKNOWN',FORM='UNFORMATTED')
      WRITE(NFILE) (((Q1(I,J,K),Q2(I,J,K),Q3(I,J,K)
     >            ,K=1,N3M),J=0,N2),I=0,N1)
      WRITE(NFILE) (((P(I,J,K),K=1,N3M),J=1,N2M),I=1,N1M)
      IMORE=IMORE+1
      CLOSE(NFILE)

      RETURN
      END

C  ************************  WRITEUP2 **********************
C     WRITE FLOW FIELD AND BOUNDARY CONDITIONS
C     AND MEAN PRESSURE GRADIENT
c     time step is specifeid in the written file name.

      SUBROUTINE WRITEUP2(Q1,Q2,Q3,P)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/LES_OPT/ILES
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      common/phaseavg/nphavg(0:mt),nphavg_les(0:mt),ntime_ctd

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      CHARACTER*20 FILENAME

      iphase = mod(ntime_ctd,mt+1)

      FILENAME='ntime.'
      NN=INDEX(FILENAME,'.')
      WRITE(UNIT=FILENAME(NN+1:),FMT='(BN,I5.5,a1,I3.3)') NTIME
     &                        ,'_',iphase

      OPEN(10,file=filename,STATUS='UNKNOWN',FORM='UNFORMATTED'
     &       ,action='write')
      WRITE(10) (((Q1(I,J,K),Q2(I,J,K),Q3(I,J,K)
     >            ,K=1,N3M),J=0,N2),I=0,N1)
      WRITE(10) (((P(I,J,K),K=1,N3M),J=1,N2M),I=1,N1M)
      CLOSE(10)

      RETURN
      END

C  ************************  SAVE_INFLOW **********************
C     WRITE INFLOW DATA
      SUBROUTINE SAVE_INFLOW(Q1,Q2,Q3)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)

      CHARACTER*10 FILENAME

c      I_SAVE=65         ! LES RE=1410
      I_SAVE=57       
      LEN_INFLOW = 2000       ! NUMER OF TIME STEPS PER ONE INFLOW DATA FILE

      if(mod(ntime-1,len_inflow).eq.0) then

        if(ntime-1.ne.0) then
           close(80)
           write(*,*) 'close file'
        endif

        filename='INFLOW.'

        nn=index(filename,'.')
        write(unit=filename(nn+1:),fmt='(bn,i3.3)')
     &       (ntime-1)/len_inflow

        open(80,file=filename,status='new',form='unformatted')

        write(*,*) filename

      endif

      WRITE(80)   ((Q1(I_SAVE  ,J,K),J=1,N2M),K=1,N3M)
     >           ,((Q2(I_SAVE-1,J,K),J=2,N2M),K=1,N3M)
     >           ,((Q3(I_SAVE-1,J,K),J=1,N2M),K=1,N3M)


      WRITE(*,20) mod(ntime-1,len_inflow)+1,'th data of'
     &          ,(ntime-1)/len_inflow+1,'th file is being written.'
20    FORMAT (I4,A10,I4,A25)

      RETURN
      END


C***************** GETUP ***********************     
      SUBROUTINE GETUP(Q1,Q2,Q3,P,TIME,EV,EVI,EVJ,EVK)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL QH1(0:M1,0:M2,0:M3)
      REAL QH2(0:M1,0:M2,0:M3)
      REAL QH3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)
      REAL DP(0:M1,0:M2,0:M3)

      REAL EV(0:M1,0:M2,0:M3)
      REAL EVI(0:M1,0:M2,0:M3)   
      REAL EVJ(0:M1,0:M2,0:M3)   
      REAL EVK(0:M1,0:M2,0:M3)   

C     INITIALIZE THE INTERMEDIATE VELOCITY AND PRESSURE
      DO 10 K=1,N3M
      DO 10 J=0,N2
      DO 10 I=0,N1
      QH1(I,J,K)=0.0
      QH2(I,J,K)=0.0
10    QH3(I,J,K)=0.0
      DO 20 K=1,N3M
      DO 20 J=1,N2M
      DO 20 I=1,N1M
 20   DP(I,J,K)=0.0

      CALL BCOND(Q1,Q2,Q3,TIME)

C     CALCULATE THE INTERMEDIATE VELOCITY
      CALL UHCALC(Q1,Q2,Q3,QH1,QH2,QH3,P,EV,EVI,EVJ,EVK)
C      CALL UHCALC_I(Q1,Q2,Q3,QH1,QH2,QH3,P,EV,EVI,EVJ,EVK)

C     CALCULATE DP
      CALL DPCALC(QH1,QH2,QH3,DP)

C     UPDATE THE N+1 TIME STEP VELOCITY AND PRESSURE
      CALL UPCALC(Q1,Q2,Q3,P,QH1,QH2,QH3,DP)

      RETURN
      END

C***************** BCOND ***********************     
      SUBROUTINE BCOND(Q1,Q2,Q3,TIME)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/BCON/UBC1(M1,M3,3),UBC2(M1,M3,3)            
     >           ,UBC3(M2,M3,3),UBC4(M2,M3,3)
      COMMON/SIZE/ALX,ALY,ALZ
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/PARA/RE

      CHARACTER*80 filename


      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)


C     LOWER AND UPPER BOUNDARY CONDITIONS  : UBC1,UBC2
C----------------------------------------------------------------------
      DO 10 K=1,N3M
      DO 10 I=1,N1 

      UBC1(I,K,1)=0.0  ! Q1(I,0,K)  LOWER WALL : NO-SLIP CONDITION
      UBC1(I,K,2)=0.0  ! Q2(I,1,K) 
      UBC1(I,K,3)=0.0  ! Q3(I,0,K)

C---  LOCAL FORCING
      PI=ACOS(-1.0)
c--Re=300
c      IS=98
c      IE=106
c      A0=0.0267  ! Re=300
c--Re=1410
c      IS=199
c      IE=203
c      A0=0.02266  ! Re=1410
c
c      PERIOD=real(mt+1)*DT
c      IF (I.GE.IS.AND.I.LT.IE) 
c     >    UBC1(I,K,2)=A0*(1.-cos(2.0*PI*TIME/PERIOD))
C---

c      UBC2(I,K,1)=Q1(I,N2M,K)  ! Q1(I,N2,K)  UPPER WALL 
c      UBC2(I,K,2)=Q2(I,N2M,K)  ! Q2(I,N2,K) 
c      UBC2(I,K,3)=Q3(I,N2M,K)  ! Q3(I,N2,K)


C     FOR INFLOW GENERATION
      UBC2(I,K,1)=Q1(I,N2M,K)  ! Q1(I,N2,K)
C      UBC2(I,K,2)=Q2(I,N2M,K) ! Q2(I,N2,K)  -> SEE GET_INFLOW
      UBC2(I,K,3)=Q3(I,N2M,K)  ! Q3(I,N2,K)

C     FOR MAIN SIMULATION
c      UBC2(I,K,1)=1.0  ! Q1(I,N2,K)
c      UBC2(I,K,2)=Q2(I,N2M,K)  ! Q2(I,N2,K)
c      UBC2(I,K,3)=Q3(I,N2M,K)  ! Q3(I,N2,K)

 10   CONTINUE


C     SIDE WALL BOUNDARY CONDITIONS  : UBC3,UBC4
C----------------------------------------------------------------------
C
C      IF (NTIME.EQ.1) CALL BLASIUS(Q1,Q2)
C      GOTO 88
C
C      DO 20 K=1,N3M
C      DO 20 J=1,N2M
C
C      UBC3(J,K,1)=1.0
C      UBC3(J,K,2)=0.0
C      UBC3(J,K,3)=0.0
C
C      UBC4(J,K,1)=Q1(N1M,J,K)
C      UBC4(J,K,2)=Q2(N1M,J,K)
C      UBC4(J,K,3)=Q3(N1M,J,K)
C
C20    CONTINUE
C
C      GOTO 89
C     INLET BOUNDARY CONDITIONS   : UBC3
C----------------------------------------------------------------------

      CALL RESCALING(Q1,Q2,Q3,TIME)              ! INFLOW SIMULATION
      GOTO 88

cc--- read time averaged inflow data u+ from the file re300_2d.in
c
c      if (ntime.eq.1) then
c    
c      open(80,file='re1410_2d.in',status='old',action='read'
c     &       ,form='formatted')
c 
c      READ (80,*) ((UBC3(J,K,1),J=1,N2M),K=1,N3M) ! U(1,J,K,1) INLET BOUNDARY
c     >         ,((UBC3(J,K,2),J=2,N2M),K=1,N3M) ! U(0,J,K,2) INLET BOUNDARY
c
c
c      UBC3(:,:,3)=0.0 ! U(0,J,K,3) INLET BOUNDARY
c
c      close(80)
c      endif
c      GOTO 88
cc---

C     READ THE INFLOW DATA FROM THE SAVED FILES

      LEN_INFLOW = 2000      ! INFLOW DATA LENGTH

      IF (MOD(NTIME-1,LEN_INFLOW).EQ.0) THEN

         IF (NTIME-1.NE.0)  THEN
            CLOSE(80)
            WRITE(*,*) 'THE PREVIOUS INFLOW DATA FILE IS CLOSED.'
         ENDIF

        filename=
     & '/hdd1/kkim/inflow1410_les/02_save/INFLOW.'
c     & '/home/kkim/pb_tbl/re1410/les/inflow6/02_save/INFLOW.'
c     & '/home/kkim/EXP_SHPARK/inflow60/INFLOW.'
c     & '/home/kkim/EXP_SHPARK/01_inflow_z225/02_save/INFLOW.'

        nn=index(filename,'.')
        write(unit=filename(nn+1:),fmt='(bn,i3.3)')
     &       (ntime-1)/len_inflow
c     &       +2

        open(80,file=filename,status='old',form='unformatted'
     &         ,action='read')

        write(*,*) filename

      ENDIF

      READ (80) ((UBC3(J,K,1),J=1,N2M),K=1,N3M) ! U(1,J,K,1) INLET BOUNDARY
     >         ,((UBC3(J,K,2),J=2,N2M),K=1,N3M) ! U(0,J,K,2) INLET BOUNDARY
     >         ,((UBC3(J,K,3),J=1,N2M),K=1,N3M) ! U(0,J,K,3) INLET BOUNDARY

88    CONTINUE 

C     CONVECTIVE BOUNDARY CONDITION FOR EXIT IN X1 DIRECTION  : UBC4
C----------------------------------------------------------------------
      UC=0.0     
      DO 31 K=1,N3M
      DO 31 J=1,N2M
      UC=UC+Q1(N1,J,K)*DY(J)/DX3
31    CONTINUE
      UC=UC/ALZ/ALY            ! BULK VELOCITY
      DO 32 K=1,N3M
      DO 32 J=1,N2M
      UBC4(J,K,1)=Q1(N1,J,K)-DT/DX(N1M)*UC*(Q1(N1,J,K)-Q1(N1M,J,K))  ! Q1(N1,J,K,1)
      UBC4(J,K,2)=Q2(N1,J,K)-DT/G(N1)*UC*(Q2(N1,J,K)-Q2(N1M,J,K))  ! Q2(N1,J,K,2)
32    UBC4(J,K,3)=Q3(N1,J,K)-DT/G(N1)*UC*(Q3(N1,J,K)-Q3(N1M,J,K))  ! Q3(N1,J,K,3)



C     THE INTERGAL OF MASS FLUX MUST BE ZERO!!!!
      Q_IN=0.
      Q_UP=0.
      Q_DOWN=0.
      Q_EX=0.
      DO 33 K=1,N3M
      DO 33 J=1,N2M
      Q_IN=Q_IN+UBC3(J,K,1)*DY(J)/DX3
33    Q_EX=Q_EX+UBC4(J,K,1)*DY(J)/DX3

      DO 35 K=1,N3M
      DO 35 I=1,N1M
      Q_UP=Q_UP+UBC2(I,K,2)*DX(I)/DX3
35    Q_DOWN=Q_DOWN+UBC1(I,K,2)*DX(I)/DX3
     
      WRITE(*,*) 
      WRITE(*,100) Q_IN,Q_UP,Q_EX
100   FORMAT('MASS FLUX : Q_IN=',E10.4,X,'Q_UP=',E10.4,X,'Q_EX=',E10.4)

      OPEN(22,FILE='FLOW.PLT',STATUS='UNKNOWN',POSITION='APPEND') 
      WRITE(22,99) TIME,Q_IN,Q_UP,Q_EX
99    FORMAT(4(E15.7,X))
      CLOSE(22)

      RATE=(Q_IN+Q_DOWN-Q_UP)/Q_EX
 
      DO 34 K=1,N3M
      DO 34 J=1,N2M
      UBC4(J,K,1)=RATE*UBC4(J,K,1)
      UBC4(J,K,2)=RATE*UBC4(J,K,2)
34    UBC4(J,K,3)=RATE*UBC4(J,K,3)

89    CONTINUE

      RETURN
      END

C*************** BLASIUS ***********************     
      SUBROUTINE BLASIUS(Q1,Q2)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/BCON/UBC1(M1,M3,3),UBC2(M1,M3,3)            
     >           ,UBC3(M2,M3,3),UBC4(M2,M3,3)
      COMMON/SIZE/ALX,ALY,ALZ
      COMMON/FINOUT/IDTOPT,NWRITE,NREAD,IAVG,IAVG_TRIG,NPRN,INSF,NINS

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)

      REAL F1(M2),F2(M2),F3(M2)  ! F1=F", F2=F', F3=F

      F2(1)=0.0
      F3(1)=0.0

      F1_KM=0.1
	F1(1)=F1_KM

      DO J=1,N2M
          D=DY(J)*F1(1)
CX0          D=DY(J)/SQRT(2.)*SQRT(RE)

          F1(J+1)=F1(J)*(1.-0.125*D**3*F1(J)-0.5*D**2*F2(J)-0.5*D*F3(J))
     >           /(1.+0.125*D**3*F1(J)+0.5*D*F3(J))

          F2(J+1)=(1.*D*F1(J)+1.*F2(J)-0.125*D**3*F1(J)*F2(J)
     >	         +0.5*D*F2(J)*F3(J))
     >           /(1.+0.125*D**3*F1(J)+0.5*D*F3(J))

          F3(J+1)=(0.5*D**2*F1(J)+1.*D*F2(J)+1.*F3(J)
     >	       +0.125*D**3*F1(J)*F3(J)+0.5*D**2*F2(J)*F3(J)
     >		   +0.5*D*F3(J)**2)
     >		   /(1.+0.125*D**3*F1(J)+0.5*D*F3(J))
      ENDDO

      F2_KM=F2(N2)


      F1_KC=0.5

10    CONTINUE
	F1(1)=F1_KC

      DO J=1,N2M
	    D=DY(J)*F1(1)
CX0          D=DY(J)/SQRT(2.)*SQRT(RE)

          F1(J+1)=F1(J)*(1.-0.125*D**3*F1(J)-0.5*D**2*F2(J)-0.5*D*F3(J))
     >           /(1.+0.125*D**3*F1(J)+0.5*D*F3(J))

          F2(J+1)=(1.*D*F1(J)+1.*F2(J)-0.125*D**3*F1(J)*F2(J)
     >	         +0.5*D*F2(J)*F3(J))
     >           /(1.+0.125*D**3*F1(J)+0.5*D*F3(J))

          F3(J+1)=(0.5*D**2*F1(J)+1.*D*F2(J)+1.*F3(J)
     >	       +0.125*D**3*F1(J)*F3(J)+0.5*D**2*F2(J)*F3(J)
     >		   +0.5*D*F3(J)**2)
     >		   /(1.+0.125*D**3*F1(J)+0.5*D*F3(J))
      ENDDO

      F2_KC=F2(N2)

      WRITE(*,*) F1_KC,F2_KC

      IF (ABS(F2_KC-1.).LE.1.E-7) GOTO 20

      TEMP=F1_KC+(F1_KC-F1_KM)/(F2_KC-F2_KM)*(1.-F2_KC)

      F1_KM=F1_KC
      F2_KM=F2_KC
	F1_KC=TEMP
      
	GOTO 10

20    CONTINUE


C---  IMPOSE THE BLASIUS SOLUTION TO INLET BOUNDARY CONDITION

      DO K=1,N3M
      DO J=1,N2M
         UBC3(J,K,1)=0.5*(F2(J)+F2(J+1))
         ETHA=Y(J)*F1(1)
         UBC3(J,K,2)=F1(1)/RE*(ETHA*F2(J)-F3(J))
CX0         ETHA=Y(J)/SQRT(2.)*SQRT(RE)
CX0         UBC3(J,K,2)=1./SQRT(2.*RE)*(ETHA*F2(J)-F3(J))
         UBC3(J,K,3)=0.0
      ENDDO
	ENDDO

C---  IF NREAD=0, IMPOSE THE INITIAL FIELD AS BLASIUS SOLUTION ROUGHLY
C     IN ORDER TO PREVENT THE MASS FLOW RATE OVER UPPER PLANE FROM BEING NEGATIVE

      IF (NREAD.EQ.1) RETURN

      DO K=1,N3M
      DO I=0,N1
         DO J=1,N2M
            Q1(I,J,K)=UBC3(J,K,1)
            Q2(I,J,K)=UBC3(J,K,2)
         ENDDO
            Q1(I,N2,K)=UBC3(N2M,K,1)
            Q2(I,N2,K)=UBC3(N2M,K,2)
      ENDDO
      ENDDO


      OPEN(10,FILE='INLET.PLT',STATUS='UNKNOWN')
      DO J=1,N2M-1
	   WRITE(10,101) 0.5*(Y(J)+Y(J+1))
     > 	            ,UBC3(J,1,1)
     >                ,0.5*(UBC3(J,1,2)+UBC3(J+1,1,2))
101      FORMAT(3(E12.5,X))
      ENDDO
      CLOSE(10)


      OPEN(10,FILE='BLASIUS.PLT',STATUS='UNKNOWN')
      DO J=1,N2
	   WRITE(10,100) F1(1)*Y(J), F1(J),F2(J),F3(J)
100       FORMAT(4(E12.5,X))
      ENDDO
      CLOSE(10)

      RETURN
      END

C***************** UHCALC ***********************     
      SUBROUTINE UHCALC(Q1,Q2,Q3,QH1,QH2,QH3,P,EV,EVI,EVJ,EVK)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL QH1(0:M1,0:M2,0:M3)
      REAL QH2(0:M1,0:M2,0:M3)
      REAL QH3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL RDUH1(0:M1,0:M2,0:M3)
      REAL RDUH2(0:M1,0:M2,0:M3)
      REAL RDUH3(0:M1,0:M2,0:M3)

      REAL EV(0:M1,0:M2,0:M3)
      REAL EVI(0:M1,0:M2,0:M3)
      REAL EVJ(0:M1,0:M2,0:M3)
      REAL EVK(0:M1,0:M2,0:M3)


      CALL RHS1(Q1,Q2,Q3,P,RDUH1,EV,EVI,EVJ,EVK)
      CALL RHS2(Q1,Q2,Q3,P,RDUH2,EV,EVI,EVJ,EVK)
      CALL RHS3(Q1,Q2,Q3,P,RDUH3,EV,EVI,EVJ,EVK)

      QH1(:,:,:)=0.0
      QH2(:,:,:)=0.0
      QH3(:,:,:)=0.0

      WRITE(*,*) 'GET INTERMEDIATE VELOCITY BY DECOUPLING'

      CALL GETDUH1(Q1,Q2,Q3,QH1,QH2,QH3,RDUH1,EV,EVI,EVJ,EVK)
      CALL GETDUH2(Q1,Q2,Q3,QH1,QH2,QH3,RDUH2,EV,EVI,EVJ,EVK)
      CALL GETDUH3(Q1,Q2,Q3,QH1,QH2,QH3,RDUH3,EV,EVI,EVJ,EVK)

C     INTERMEDIATE VELOCITY, U*=U^N+DU*

!$omp parallel do shared(QH1)
      DO 300 K=1,N3M
      DO 300 J=1,N2M
      DO 300 I=2,N1M
300   QH1(I,J,K)=Q1(I,J,K)+QH1(I,J,K)

!$omp parallel do shared(QH2)
      DO 310 K=1,N3M
      DO 310 J=2,N2M
      DO 310 I=1,N1M
310   QH2(I,J,K)=Q2(I,J,K)+QH2(I,J,K)

!$omp parallel do shared(QH3)
      DO 320 K=1,N3M
      DO 320 J=1,N2M
      DO 320 I=1,N1M
320   QH3(I,J,K)=Q3(I,J,K)+QH3(I,J,K)

 
      RETURN
      END
      
C***************** UHCALC_I ***********************     
      SUBROUTINE UHCALC_I(Q1,Q2,Q3,QH1,QH2,QH3,P,EV,EVI,EVJ,EVK)
      INCLUDE 'PARAM.H'
      INCLUDE 'PARAM_ITER.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL QH1(0:M1,0:M2,0:M3)
      REAL QH2(0:M1,0:M2,0:M3)
      REAL QH3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL RDUH1(0:M1,0:M2,0:M3)
      REAL RDUH2(0:M1,0:M2,0:M3)
      REAL RDUH3(0:M1,0:M2,0:M3)

      REAL EV(0:M1,0:M2,0:M3)
      REAL EVI(0:M1,0:M2,0:M3)   
      REAL EVJ(0:M1,0:M2,0:M3)   
      REAL EVK(0:M1,0:M2,0:M3)   


      CALL RHS1(Q1,Q2,Q3,P,RDUH1,EV,EVI,EVJ,EVK)
      CALL RHS2(Q1,Q2,Q3,P,RDUH2,EV,EVI,EVJ,EVK)
      CALL RHS3(Q1,Q2,Q3,P,RDUH3,EV,EVI,EVJ,EVK)


      QH1(:,:,:)=0.0
      QH2(:,:,:)=0.0
      QH3(:,:,:)=0.0


      WRITE(*,*) 'FIRST GUESS INTERMEDIATE VELOCITY BY DECOUPLING'

      CALL GETDUH1(Q1,Q2,Q3,QH1,QH2,QH3,RDUH1,EV,EVI,EVJ,EVK)
      CALL GETDUH2(Q1,Q2,Q3,QH1,QH2,QH3,RDUH2,EV,EVI,EVJ,EVK)
      CALL GETDUH3(Q1,Q2,Q3,QH1,QH2,QH3,RDUH3,EV,EVI,EVJ,EVK)

      WRITE(*,*)
      WRITE(*,*) 'GET INTERMEDIATE VELOCITY BY ITERATION'
      WRITE(*,*) '  ITER    RES1             '
      WRITE(*,*) '---------------------------'

      NITER=0
      RES1=0.0 
 
1     CONTINUE

      NITER=NITER+1

      CALL GETDUH1_I(Q1,Q2,Q3,QH1,QH2,QH3,RDUH1,EV,EVI,EVJ,EVK)
      CALL GETDUH2_I(Q1,Q2,Q3,QH1,QH2,QH3,RDUH2,EV,EVI,EVJ,EVK)
      CALL GETDUH3_I(Q1,Q2,Q3,QH1,QH2,QH3,RDUH3,EV,EVI,EVJ,EVK)

C--   CHECK RES1 IS SUFFICIENT. NOTE THE RESEARCH NOTE
      CALL CHECK_RES1_AF(Q1,Q2,Q3,QH1,QH2,QH3,RDUH1,EV,EVI,EVJ,EVK,RES1)
     
      WRITE(*,100) NITER,RES1
100   FORMAT (I5,X,1(E12.5,X))

      IF(NITER.GE.NITERMAX) THEN
         WRITE(*,*) 'INTERMEDIATE VELOCITIES ARE NOT CONVERGED!'
         GOTO 2
      ENDIF

      IF(RES1.GE.ERRMAX) GOTO 1

     
2     CONTINUE


C     INTERMEDIATE VELOCITY, U*=U^N+DU*

      DO 300 K=1,N3M
      DO 300 J=1,N2M
      DO 300 I=2,N1M
300   QH1(I,J,K)=Q1(I,J,K)+QH1(I,J,K)

      DO 310 K=1,N3M
      DO 310 J=2,N2M
      DO 310 I=1,N1M
310   QH2(I,J,K)=Q2(I,J,K)+QH2(I,J,K)

      DO 320 K=1,N3M
      DO 320 J=1,N2M
      DO 320 I=1,N1M
320   QH3(I,J,K)=Q3(I,J,K)+QH3(I,J,K)

 
      RETURN
      END      

C***************** DPCALC ***********************     
      SUBROUTINE DPCALC(QH1,QH2,QH3,DP)
      INCLUDE 'PARAM.H'
      
      REAL QH1(0:M1,0:M2,0:M3)
      REAL QH2(0:M1,0:M2,0:M3)
      REAL QH3(0:M1,0:M2,0:M3)
      REAL DP(0:M1,0:M2,0:M3)
      REAL RDP(0:M1,0:M2,0:M3)

      CALL RHSDP(RDP,QH1,QH2,QH3)
      CALL TAKEDP(DP,RDP)

      RETURN
      END

C***************** UPCALC ***********************     
      SUBROUTINE UPCALC(Q1,Q2,Q3,P,QH1,QH2,QH3,DP)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/BCON/UBC1(M1,M3,3),UBC2(M1,M3,3)            
     >           ,UBC3(M2,M3,3),UBC4(M2,M3,3)
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/SIZE/ALX,ALY,ALZ

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL QH1(0:M1,0:M2,0:M3)
      REAL QH2(0:M1,0:M2,0:M3)
      REAL QH3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)
      REAL DP(0:M1,0:M2,0:M3)
      
C     U1 VELOCITY UPDATE
     
!$omp parallel do shared(Q1)
      DO 12 K=1,N3M
      DO I=1,N1
      Q1(I,0,K)=UBC1(I,K,1)
      ENDDO
      DO J=1,N2M
         Q1(1,J,K)=UBC3(J,K,1)
      DO I=2,N1M
         Q1(I,J,K)=QH1(I,J,K)
     >          -DT*(DP(I,J,K)-DP(I-1,J,K))/G(I)
      ENDDO
         Q1(N1,J,K)=UBC4(J,K,1)
      ENDDO
      DO I=1,N1
         Q1(I,N2,K)=UBC2(I,K,1)
      ENDDO
       
  12  CONTINUE

C     U2 VELOCITY UPDATE

!$omp parallel do shared(Q2)
      DO 22 K=1,N3M
      DO I=1,N1
         Q2(I,1,K)=UBC1(I,K,2)
      ENDDO
         J=1
         Q2(0,J,K)=UBC3(J,K,2)
         Q2(N1,J,K)=UBC4(J,K,2)
      DO J=2,N2M
         Q2(0,J,K)=UBC3(J,K,2)
      DO I=1,N1M
         Q2(I,J,K)=QH2(I,J,K)
     >          -DT*(DP(I,J,K)-DP(I,J-1,K))/H(J)
      ENDDO
         Q2(N1,J,K)=UBC4(J,K,2)
      ENDDO
      DO I=1,N1
         Q2(I,N2,K)=UBC2(I,K,2)
      ENDDO
       
  22  CONTINUE

C     U3 VELOCITY UPDATE

!$omp parallel do shared(Q3) private(KM)
      DO 32 K=1,N3M
      KM=KMA(K)
      DO I=1,N1M
         Q3(I,0,K)=UBC1(I,K,3)
      ENDDO
      DO J=1,N2M
         Q3(0,J,K)=UBC3(J,K,3)
      DO I=1,N1M
         Q3(I,J,K)=QH3(I,J,K)
     >          -DT*(DP(I,J,K)-DP(I,J,KM))*DX3
      ENDDO
         Q3(N1,J,K)=UBC4(J,K,3)
      ENDDO
      DO I=1,N1M
         Q3(I,N2,K)=UBC2(I,K,3)
      ENDDO
       
  32  CONTINUE

C     PRESSURE UPDATE     
 
!$omp parallel do shared(P)
      DO 40 K=1,N3M
      DO 40 J=1,N2M
      DO 40 I=1,N1M
      P(I,J,K)=P(I,J,K)+DP(I,J,K)
  40  CONTINUE
      
      RETURN
      END
      
C***************** RHS1 ***********************     
      SUBROUTINE RHS1(Q1,Q2,Q3,P,RDUH1,EV,EVI,EVJ,EVK)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/BCON/UBC1(M1,M3,3),UBC2(M1,M3,3)            
     >           ,UBC3(M2,M3,3),UBC4(M2,M3,3)
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)
      REAL RDUH1(0:M1,0:M2,0:M3)

      REAL EV(0:M1,0:M2,0:M3)
      REAL EVI(0:M1,0:M2,0:M3)   
      REAL EVJ(0:M1,0:M2,0:M3)   
      REAL EVK(0:M1,0:M2,0:M3)   

!$omp parallel do private(KP,KM,JP,JM,IP,IM,JUM,JUP,IUM,IUP,
!$omp&            VIS1,VIS2,
!$omp$            VISCOS_11,BC_IN_11,BC_OUT_11,
!$omp$            VISCOS_12,BC_DOWN_12,BC_UP_12,
!$omp$            VISCOS_13,
!$omp$            VISCOS,PRESSG1,
!$omp$            BC_IN,BC_OUT,BC_DOWN,BC_UP,BC,
!$omp$            U1,U2,V1,V2,W1,W2,
!$omp$            RM11U_N,RM12V_N,RM13W_N,
!$omp$            API,ACI,AMI,APJ,ACJ,AMJ,APK,ACK,AMK)
!$omp$         shared(RDUH1)


      DO 10 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 10 J=1,N2M
      JP=J+1
      JM=J-1
      JUM=J-JMU(J)
      JUP=JPA(J)-J
      DO 10 I=2,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMU(I)
      IUP=IPA(I)-I


C---  VISCOS_11

      VIS1 = EV(IM,J,K) + 1./RE
      VIS2 = EV(I ,J,K) + 1./RE

      VISCOS_11 = 1.0/G(I)*       VIS2/DX(I )      *Q1(IP,J,K)
     >           -1.0/G(I)*(VIS2/DX(I)+VIS1/DX(IM))*Q1(I ,J ,K)
     >           +1.0/G(I)*       VIS1/DX(IM)      *Q1(IM,J,K)

      BC_IN_11  =
     >           +1.0/G(I)*       VIS1/DX(IM)      *UBC3(J,K,1)
      BC_OUT_11 = 
     >           +1.0/G(I)*       VIS2/DX(I )      *UBC4(J,K,1)

C---  VISCOS_12

      VIS1 = EVK(I,J ,K) + 1.0/RE
      VIS2 = EVK(I,JP,K) + 1.0/RE


      VISCOS_12 = 0.5/DY(J)*       VIS2/H(JP)     *Q1(I,JP,K)
     >           -0.5/DY(J)*(VIS2/H(JP)+VIS1/H(J))*Q1(I,J ,K)
     >           +0.5/DY(J)*       VIS1/H(J )     *Q1(I,JM,K)
     >        
     >           +0.5/DY(J)/G(I)*( VIS2*(Q2(I,JP,K)-Q2(IM,JP,K))
     >                            -VIS1*(Q2(I,J, K)-Q2(IM,J ,K)))

      BC_DOWN_12 =  
     >            +0.5/DY(J)*       VIS1/H(J )     *UBC1(I,K,1)
     >            +0.5/DY(J)/G(I)*(-VIS1*(UBC1(I,K,2)-UBC1(IM,K,2)))

      BC_UP_12   =  
     >            +0.5/DY(J)*       VIS2/H(JP)     *UBC2(I,K,1)
     >            +0.5/DY(J)/G(I)*( VIS2*(UBC2(I,K,2)-UBC2(IM,K,2)))



C---  VISCOS_13

      VIS1 = EVJ(I,J,K ) + 1.0/RE
      VIS2 = EVJ(I,J,KP) + 1.0/RE

      VISCOS_13 = 0.5*DX3Q*    VIS2   *Q1(I,J,KP)
     >           -0.5*DX3Q*(VIS2+VIS1)*Q1(I,J,K )
     >           +0.5*DX3Q*    VIS1   *Q1(I,J,KM)
     >        
     >           +0.5*DX3/G(I)*(VIS2*(Q3(I,J,KP)-Q3(IM,J,KP))
     >                         -VIS1*(Q3(I,J,K )-Q3(IM,J,K )))


      VISCOS=VISCOS_11+VISCOS_12+VISCOS_13


      V1=1.0/G(I)*(0.5*DX(IM)*Q2(I,J ,K)+0.5*DX(I)*Q2(IM,J ,K))
      V2=1.0/G(I)*(0.5*DX(IM)*Q2(I,JP,K)+0.5*DX(I)*Q2(IM,JP,K))

      U1=1.0/H(J )*(0.5*DY(JM)*Q1(I,J ,K)+0.5*DY(J )*Q1(I,JM,K))
      U2=1.0/H(JP)*(0.5*DY(J )*Q1(I,JP,K)+0.5*DY(JP)*Q1(I,J ,K))

      BC_DOWN=BC_DOWN_12
     >+0.5/DY(J)*V1/H(J )*0.5*DY(J)*UBC1(I,K,1)
     >+0.5/DY(J)*U1/G(I)*(0.5*DX(IM)*UBC1(I,K,2)+0.5*DX(I)*UBC1(IM,K,2)) 

      BC_UP  =BC_UP_12
     >-0.5/DY(J)*V2/H(JP)*0.5*DY(J)*UBC2(I,K,1)
     >-0.5/DY(J)*U2/G(I)*(0.5*DX(IM)*UBC2(I,K,2)+0.5*DX(I)*UBC2(IM,K,2)) 


      U1=0.5*(Q1(I, J,K)+Q1(IM,J,K)) 
      U2=0.5*(Q1(IP,J,K)+Q1(I, J,K))

      BC_IN = 1.0/G(I)*U1*0.5*UBC3(J,K,1)
     >       +BC_IN_11
      BC_OUT=-1.0/G(I)*U2*0.5*UBC4(J,K,1)
     >       +BC_OUT_11

      BC=(1-JUM)*BC_DOWN
     >  +(1-JUP)*BC_UP 
     >  +(1-IUM)*BC_IN 
     >  +(1-IUP)*BC_OUT 


      PRESSG1=(P(I,J,K)-P(IM,J,K))/G(I)

      RDUH1(I,J,K)=1./DT*Q1(I,J,K)
     >           -PRESSG1+VISCOS
     >           +BC


C     R1=R1-AU^N

C--   M12V^N

      VIS1 = EVK(I,J ,K) + 1.0/RE
      VIS2 = EVK(I,JP,K) + 1.0/RE

      V1=1.0/G(I)*(0.5*DX(IM)*Q2(I,J ,K)+0.5*DX(I)*Q2(IM,J ,K))
      V2=1.0/G(I)*(0.5*DX(IM)*Q2(I,JP,K)+0.5*DX(I)*Q2(IM,JP,K))

      U1=1.0/H(J )*(0.5*DY(JM)*Q1(I,J ,K)+0.5*DY(J )*Q1(I,JM,K))
      U2=1.0/H(JP)*(0.5*DY(J )*Q1(I,JP,K)+0.5*DY(JP)*Q1(I,J ,K))

      RM12V_N=JUP*0.5/DY(J)*U2*V2
     >       -JUM*0.5/DY(J)*U1*V1
     >
     >           -0.5/DY(J)/G(I)*( JUP*VIS2*(Q2(I,JP,K)-Q2(IM,JP,K))
     >                            -JUM*VIS1*(Q2(I,J ,K)-Q2(IM,J ,K)))
     
C--   M11U^N_Y      

      APJ=JUP*(
     >      -0.5/DY(J)*VIS2/H(JP)
     >      +0.5/DY(J)*V2/H(JP)*DY(J)/2.0
     >      )

      ACJ= 0.5/DY(J)*(VIS2/H(JP)+VIS1/H(J))
     >    +0.5/DY(J)*(JUP*V2/H(JP)*DY(JP)/2.0
     >               -JUM*V1/H(J )*DY(JM)/2.0)

      AMJ=JUM*(
     >      -0.5/DY(J)*VIS1/H(J )
     >      -0.5/DY(J)*V1/H(J)*DY(J)/2.0
     >      )

C--   M11U^N_X      

      VIS1 = EV(IM,J,K) + 1.0/RE
      VIS2 = EV(I ,J,K) + 1.0/RE

      U2=0.5*(Q1(IP,J,K)+Q1(I, J,K))
      U1=0.5*(Q1(I, J,K)+Q1(IM,J,K)) 

      API= (-1.0/G(I)*       VIS2/DX(I )      
     >      +1.0/G(I)*U2*0.5)*IUP
      ACI=   1.0/G(I)*(VIS2/DX(I)+VIS1/DX(IM))
     >      +1.0/G(I)*(U2*0.5-U1*0.5)
      AMI= (-1.0/G(I)*       VIS1/DX(IM)  
     >      -1.0/G(I)*U1*0.5)*IUM


C--   M13W^N
 
      VIS1 = EVJ(I,J,K ) + 1.0/RE
      VIS2 = EVJ(I,J,KP) + 1.0/RE

      W1=1.0/G(I)*(0.5*DX(IM)*Q3(I,J,K )+0.5*DX(I)*Q3(IM,J,K ))
      W2=1.0/G(I)*(0.5*DX(IM)*Q3(I,J,KP)+0.5*DX(I)*Q3(IM,J,KP))

      U1=0.5*(Q1(I,J,K )+Q1(I,J,KM))
      U2=0.5*(Q1(I,J,KP)+Q1(I,J,K ))

      RM13W_N=0.5*DX3*(U2*W2-U1*W1)
     >
     >           -0.5*DX3/G(I)*(VIS2*(Q3(I,J,KP)-Q3(IM,J,KP))
     >                         -VIS1*(Q3(I,J,K )-Q3(IM,J,K )))

C--   M11U^N_Z      
      APK=   -0.5*DX3Q*    VIS2  
     >       +0.5*DX3*W2*0.5
      ACK=   +0.5*DX3Q*(VIS2+VIS1)
     >       +0.5*DX3*(W2*0.5-W1*0.5)
      AMK=   -0.5*DX3Q*    VIS1  
     >       -0.5*DX3*W1*0.5


      RM11U_N=APJ*Q1(I,JP,K)
     >       +ACJ*Q1(I,J ,K)
     >       +AMJ*Q1(I,JM,K)
     >       +API*Q1(IP,J,K)
     >       +ACI*Q1(I, J,K)
     >       +AMI*Q1(IM,J,K)
     >       +APK*Q1(I,J,KP)
     >       +ACK*Q1(I,J,K )
     >       +AMK*Q1(I,J,KM)


      RDUH1(I,J,K)=RDUH1(I,J,K)
     >           -1./DT*Q1(I,J,K)
     >           -RM11U_N-RM12V_N-RM13W_N

 10   CONTINUE 


      RETURN
      END



C***************** RHS2 ***********************     
      SUBROUTINE RHS2(Q1,Q2,Q3,P,RDUH2,EV,EVI,EVJ,EVK)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/BCON/UBC1(M1,M3,3),UBC2(M1,M3,3)            
     >           ,UBC3(M2,M3,3),UBC4(M2,M3,3)
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)
      REAL RDUH2(0:M1,0:M2,0:M3)

      REAL EV(0:M1,0:M2,0:M3)
      REAL EVI(0:M1,0:M2,0:M3)   
      REAL EVJ(0:M1,0:M2,0:M3)   
      REAL EVK(0:M1,0:M2,0:M3)   

!$omp parallel do private(KP,KM,JP,JM,IP,IM,JUM,JUP,IUM,IUP,
!$omp&            VIS1,VIS2,
!$omp$            VISCOS_21,BC_IN_21,BC_OUT_21,
!$omp$            VISCOS_22,BC_DOWN_22,BC_UP_22,
!$omp$            VISCOS_23,
!$omp$            VISCOS,PRESSG2,
!$omp$            BC_IN,BC_OUT,BC_DOWN,BC_UP,BC,
!$omp$            U1,U2,V1,V2,W1,W2,
!$omp$            RM21U_N,RM22V_N,RM23W_N,
!$omp$            API,ACI,AMI,APJ,ACJ,AMJ,APK,ACK,AMK)
!$omp$         shared(RDUH2)

      DO 10 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      
      DO 10 J=2,N2M
      JP=J+1
      JM=J-1
      JUM=J-JMV(J)
      JUP=JPA(J)-J

      DO 10 I=1,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMV(I)
      IUP=IPA(I)-I
     
C---  VISCOS_21

      VIS1 = EVK(I ,J,K) + 1.0/RE
      VIS2 = EVK(IP,J,K) + 1.0/RE

      VISCOS_21 = 0.5/DX(I)*       VIS2/G(IP)     *Q2(IP,J,K)
     >           -0.5/DX(I)*(VIS2/G(IP)+VIS1/G(I))*Q2(I ,J,K)
     >           +0.5/DX(I)*       VIS1/G(I )     *Q2(IM,J,K)
     >
     >           +0.5/DX(I)/H(J)*( VIS2*(Q1(IP,J,K)-Q1(IP,JM,K))
     >                            -VIS1*(Q1(I, J,K)-Q1(I ,JM,K)))

      BC_IN_21 = +0.5/DX(I)*       VIS1/G(I )     *UBC3(J,K,2)
     >           +0.5/DX(I)/H(J)*(-VIS1*(UBC3(J,K,1)-UBC3(JM,K,1)))
      BC_OUT_21= +0.5/DX(I)*       VIS2/G(IP)     *UBC4(J,K,2)
     >           +0.5/DX(I)/H(J)*( VIS2*(UBC4(J,K,1)-UBC4(JM,K,1)))

C---  VISCOS_22

      VIS1 = EV(I,JM,K) + 1.0/RE
      VIS2 = EV(I,J ,K) + 1.0/RE

      VISCOS_22 = 1.0/H(J)*       VIS2/DY(J )      *Q2(I,JP,K)
     >           -1.0/H(J)*(VIS2/DY(J)+VIS1/DY(JM))*Q2(I,J ,K)
     >           +1.0/H(J)*       VIS1/DY(JM)      *Q2(I,JM,K)

      BC_UP_22  = 1.0/H(J)*       VIS2/DY(J )      *UBC2(I,K,2)
      BC_DOWN_22= 1.0/H(J)*       VIS1/DY(JM)      *UBC1(I,K,2)

C---  VISCOS_23

      VIS1 = EVI(I,J,K ) + 1.0/RE
      VIS2 = EVI(I,J,KP) + 1.0/RE

      VISCOS_23 = 0.5*DX3Q*    VIS2   *Q2(I,J,KP)
     >           -0.5*DX3Q*(VIS2+VIS1)*Q2(I,J,K )
     >           +0.5*DX3Q*    VIS1   *Q2(I,J,KM)
     >        
     >           +0.5*DX3/H(J)*(VIS2*(Q3(I,J,KP)-Q3(I,JM,KP))
     >                         -VIS1*(Q3(I,J,K )-Q3(I,JM,K )))

      VISCOS=VISCOS_21+VISCOS_22+VISCOS_23

      
      PRESSG2=(P(I,J,K)-P(I,JM,K))/H(J) 
      
      V1=0.5*(Q2(I,J ,K)+Q2(I,JM,K))
      V2=0.5*(Q2(I,JP,K)+Q2(I,J ,K))

      BC_DOWN=BC_DOWN_22          
     >       +1./H(J)*V1*0.5*UBC1(I,K,2)
      BC_UP  =BC_UP_22
     >       -1./H(J)*V2*0.5*UBC2(I,K,2)
 
      U1=1.0/H(J)*(0.5*DY(JM)*Q1(I ,J,K)+0.5*DY(J)*Q1(I ,JM,K))
      U2=1.0/H(J)*(0.5*DY(JM)*Q1(IP,J,K)+0.5*DY(J)*Q1(IP,JM,K))

      V1=1.0/G(I )*(0.5*DX(IM)*Q2(I ,J,K)+0.5*DX(I )*Q2(IM,J,K))
      V2=1.0/G(IP)*(0.5*DX(I )*Q2(IP,J,K)+0.5*DX(IP)*Q2(I ,J,K))

      BC_IN =BC_IN_21
     >      +0.5/DX(I)*U1/G(I)*0.5*DX(I)*UBC3(J,K,2)
     >+0.5/DX(I)*V1/H(J)*(0.5*DY(JM)*UBC3(J,K,1)+0.5*DY(J)*UBC3(JM,K,1))
      BC_OUT=BC_OUT_21
     >      -0.5/DX(I)*U2/G(IP)*0.5*DX(I)*UBC4(J,K,2)
     >-0.5/DX(I)*V2/H(J)*(0.5*DY(JM)*UBC4(J,K,1)+0.5*DY(J)*UBC4(JM,K,1))


      BC=(1-JUM)*BC_DOWN
     >  +(1-JUP)*BC_UP 
     >  +(1-IUM)*BC_IN 
     >  +(1-IUP)*BC_OUT 

      RDUH2(I,J,K)=1./DT*Q2(I,J,K)
     >           -PRESSG2+VISCOS
     >           +BC

C     R2=R2-AU^N      

C--   M22V^N_Y

      VIS1 = EV(I,JM,K) + 1.0/RE
      VIS2 = EV(I,J ,K) + 1.0/RE

      V1=0.5*(Q2(I,J ,K)+Q2(I,JM,K))
      V2=0.5*(Q2(I,JP,K)+Q2(I,J ,K))

      APJ=JUP*(
     >      -1.0/H(J)*       VIS2/DY(J ) 
     >      +1.0/H(J)*V2*0.5
     >      )
      ACJ=
     >      +1.0/H(J)*(VIS2/DY(J)+VIS1/DY(JM))
     >      +1.0/H(J)*(V2*0.5-V1*0.5)
      AMJ=JUM*(
     >      -1.0/H(J)*       VIS1/DY(JM) 
     >      -1.0/H(J)*V1*0.5
     >      )


C--   M21U^N

      VIS1 = EVK(I ,J,K) + 1.0/RE
      VIS2 = EVK(IP,J,K) + 1.0/RE

      U1=1.0/H(J)*(0.5*DY(JM)*Q1(I ,J,K)+0.5*DY(J)*Q1(I ,JM,K))
      U2=1.0/H(J)*(0.5*DY(JM)*Q1(IP,J,K)+0.5*DY(J)*Q1(IP,JM,K))

      V1=1.0/G(I )*(0.5*DX(IM)*Q2(I ,J,K)+0.5*DX(I )*Q2(IM,J,K))
      V2=1.0/G(IP)*(0.5*DX(I )*Q2(IP,J,K)+0.5*DX(IP)*Q2(I ,J,K))

      RM21U_N=0.5/DX(I)*(IUP*V2*U2-IUM*V1*U1)
     >
     >       -0.5/DX(I)/H(J)*( IUP*VIS2*(Q1(IP,J,K)-Q1(IP,JM,K))
     >                        -IUM*VIS1*(Q1(I, J,K)-Q1(I ,JM,K)))


C--   M22V^N_X

      API=IUP*(
     >         -0.5/DX(I)*    VIS2/G(IP) 
     >         +0.5/DX(I)*U2/G(IP)*0.5*DX(I)
     >       )
      ACI=  
     >     +0.5/DX(I)*(VIS2/G(IP)+VIS1/G(I))
     >     +0.5/DX(I)*(U2/G(IP)*0.5*DX(IP)-U1/G(I)*0.5*DX(IM))
      AMI=IUM*(
     >         -0.5/DX(I)*    VIS1/G(I) 
     >         -0.5/DX(I)*U1/G(I)*0.5*DX(I)
     >       )

C--   M23W^N

      VIS1 = EVI(I,J,K ) + 1.0/RE
      VIS2 = EVI(I,J,KP) + 1.0/RE

      W1=1.0/H(J)*(DY(J)/2.0*Q3(I,JM,K )+DY(JM)/2.0*Q3(I,J,K ))
      W2=1.0/H(J)*(DY(J)/2.0*Q3(I,JM,KP)+DY(JM)/2.0*Q3(I,J,KP))
    
      V1=0.5*(Q2(I,J,K )+Q2(I,J,KM))
      V2=0.5*(Q2(I,J,KP)+Q2(I,J,K ))

      RM23W_N=0.5*DX3*(V2*W2-V1*W1)
     >        
     >       -0.5*DX3/H(J)*(VIS2*(Q3(I,J,KP)-Q3(I,JM,KP))
     >                     -VIS1*(Q3(I,J,K )-Q3(I,JM,K )))

C--   M22V^N_Z
      APK=
     >      - 0.5*DX3Q*    VIS2  
     >      +0.5*DX3*W2/2.0
      ACK=
     >      +0.5*DX3Q*(VIS2+VIS1)
     >      +0.5*DX3*(W2/2.0-W1/2.0)
      AMK=
     >      - 0.5*DX3Q*    VIS1  
     >      -0.5*DX3*W1/2.0

      RM22V_N=APJ*Q2(I,JP,K)
     >       +ACJ*Q2(I,J ,K)
     >       +AMJ*Q2(I,JM,K)
     >       +API*Q2(IP,J,K)
     >       +ACI*Q2(I ,J,K)
     >       +AMI*Q2(IM,J,K)
     >       +APK*Q2(I,J,KP)
     >       +ACK*Q2(I,J,K )
     >       +AMK*Q2(I,J,KM)

      RDUH2(I,J,K)=RDUH2(I,J,K)
     >           -1./DT*Q2(I,J,K)
     >           -RM21U_N-RM22V_N-RM23W_N

 10   CONTINUE 

      RETURN
      END

C***************** RHS3 ***********************     
      SUBROUTINE RHS3(Q1,Q2,Q3,P,RDUH3,EV,EVI,EVJ,EVK)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/BCON/UBC1(M1,M3,3),UBC2(M1,M3,3)            
     >           ,UBC3(M2,M3,3),UBC4(M2,M3,3)
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)
      REAL RDUH3(0:M1,0:M2,0:M3)

      REAL EV(0:M1,0:M2,0:M3)
      REAL EVI(0:M1,0:M2,0:M3)   
      REAL EVJ(0:M1,0:M2,0:M3)   
      REAL EVK(0:M1,0:M2,0:M3)   

!$omp parallel do private(KP,KM,JP,JM,IP,IM,JUM,JUP,IUM,IUP,
!$omp&            VIS1,VIS2,
!$omp$            VISCOS_31,BC_IN_31,BC_OUT_31,
!$omp$            VISCOS_32,BC_DOWN_32,BC_UP_32,
!$omp$            VISCOS_33,
!$omp$            VISCOS,PRESSG3,
!$omp$            BC_IN,BC_OUT,BC_DOWN,BC_UP,BC,
!$omp$            U1,U2,V1,V2,W1,W2,
!$omp$            RM31U_N,RM32V_N,RM33W_N,
!$omp$            API,ACI,AMI,APJ,ACJ,AMJ,APK,ACK,AMK)
!$omp$         shared(RDUH3)


      DO 10 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      
      DO 10 J=1,N2M
      JP=J+1
      JM=J-1
      JUM=J-JMU(J)
      JUP=JPA(J)-J

      DO 10 I=1,N1M          
      IP=I+1
      IM=I-1
      IUM=I-IMV(I)
      IUP=IPA(I)-I
      
    
C---  VISCOS_31

      VIS1 = EVJ(I ,J,K) + 1.0/RE
      VIS2 = EVJ(IP,J,K) + 1.0/RE

      VISCOS_31 = 0.5/DX(I)*      VIS2/G(IP)      *Q3(IP,J,K)
     >           -0.5/DX(I)*(VIS2/G(IP)+VIS1/G(I))*Q3(I ,J,K)
     >           +0.5/DX(I)*      VIS1/G(I)       *Q3(IM,J,K)
     >
     >           +0.5/DX(I)*DX3*( VIS2*(Q1(IP,J,K)-Q1(IP,J,KM))
     >                           -VIS1*(Q1(I, J,K)-Q1(I ,J,KM)))

      BC_IN_31 =  0.5/DX(I)*    VIS1/G(I)   *UBC3(J,K,3)
     >           +0.5/DX(I)*DX3*(-VIS1*(UBC3(J,K,1)-UBC3(J,KM,1)))
      BC_OUT_31=  0.5/DX(I)*    VIS2/G(IP)  *UBC4(J,K,3)
     >           +0.5/DX(I)*DX3*( VIS2*(UBC4(J,K,1)-UBC4(J,KM,1)))

C---  VISCOS_32

      VIS1 = EVI(I,J ,K) + 1.0/RE
      VIS2 = EVI(I,JP,K) + 1.0/RE

      VISCOS_32 = 0.5/DY(J)*       VIS2/H(JP)     *Q3(I,JP,K)
     >           -0.5/DY(J)*(VIS2/H(JP)+VIS1/H(J))*Q3(I,J ,K)
     >           +0.5/DY(J)*       VIS1/H(J )     *Q3(I,JM,K)
     >        
     >           +0.5/DY(J)*DX3*( VIS2*(Q2(I,JP,K)-Q2(I,JP,KM))
     >                           -VIS1*(Q2(I,J ,K)-Q2(I,J ,KM)))

      BC_DOWN_32 =  
     >            +0.5/DY(J)*       VIS1/H(J )     *UBC1(I,K,3)
     >            +0.5/DY(J)*DX3*(-VIS1*(UBC1(I,K,2)-UBC1(I,KM,2)))

      BC_UP_32   =  
     >            +0.5/DY(J)*       VIS2/H(JP)     *UBC2(I,K,3)
     >            +0.5/DY(J)*DX3*( VIS2*(UBC2(I,K,2)-UBC2(I,KM,2)))


C---  VISCOS_33

      VIS1 = EV(I,J,KM) + 1.0/RE
      VIS2 = EV(I,J,K ) + 1.0/RE

      VISCOS_33 = DX3Q*    VIS2   *Q3(I,J,KP)
     >           -DX3Q*(VIS2+VIS1)*Q3(I,J,K )
     >           +DX3Q*    VIS1   *Q3(I,J,KM)


      VISCOS=VISCOS_31+VISCOS_32+VISCOS_33

      PRESSG3=DX3*(P(I,J,K)-P(I,J,KM)) 

      V1=0.5*(Q2(I,J ,K)+Q2(I,J ,KM))
      V2=0.5*(Q2(I,JP,K)+Q2(I,JP,KM))
   
      W1=1.0/H(J )*(0.5*DY(JM)*Q3(I,J ,K)+0.5*DY(J )*Q3(I,JM,K))
      W2=1.0/H(JP)*(0.5*DY(J )*Q3(I,JP,K)+0.5*DY(JP)*Q3(I,J ,K))

      BC_DOWN=BC_DOWN_32
     >       +0.5/DY(J)*V1/H(J)*0.5*DY(J)*UBC1(I,K,3)
     >       +0.5/DY(J)*W1*0.5*(UBC1(I,K,2)+UBC1(I,KM,2))

      BC_UP  =BC_UP_32
     >       -0.5/DY(J)*V2/H(J)*0.5*DY(J)*UBC2(I,K,3)
     >       -0.5/DY(J)*W2*0.5*(UBC2(I,K,2)+UBC2(I,KM,2))


      U1=0.5*(Q1(I ,J,K)+Q1(I ,J,KM))
      U2=0.5*(Q1(IP,J,K)+Q1(IP,J,KM))

      W1=1.0/G(I )*(0.5*DX(IM)*Q3(I ,J,K)+0.5*DX(I )*Q3(IM,J,K))
      W2=1.0/G(IP)*(0.5*DX(I )*Q3(IP,J,K)+0.5*DX(IP)*Q3(I ,J,K))

      BC_IN =BC_IN_31
     >      +0.5/DX(I)*U1/G(I)*0.5*DX(I)*UBC3(J,K,3)
     >      +0.5/DX(I)*W1*0.5*(UBC3(J,K,1)+UBC3(J,KM,1))

      BC_OUT=BC_OUT_31
     >      -0.5/DX(I)*U2/G(IP)*0.5*DX(I)*UBC4(J,K,3)
     >      -0.5/DX(I)*W2*0.5*(UBC4(J,K,1)+UBC4(J,KM,1))

      BC=(1-JUM)*BC_DOWN
     >  +(1-JUP)*BC_UP 
     >  +(1-IUM)*BC_IN 
     >  +(1-IUP)*BC_OUT 

      RDUH3(I,J,K)=1./DT*Q3(I,J,K)
     >           -PRESSG3+VISCOS
     >           +BC

C     R3=R3-AU^N      


C--   M32V^N

      VIS1 = EVI(I,J ,K) + 1.0/RE
      VIS2 = EVI(I,JP,K) + 1.0/RE

      V1=0.5*(Q2(I,J ,K)+Q2(I,J ,KM))
      V2=0.5*(Q2(I,JP,K)+Q2(I,JP,KM))
   
      W1=1.0/H(J )*(0.5*DY(JM)*Q3(I,J ,K)+0.5*DY(J )*Q3(I,JM,K))
      W2=1.0/H(JP)*(0.5*DY(J )*Q3(I,JP,K)+0.5*DY(JP)*Q3(I,J ,K))

      RM32V_N=JUP*0.5/DY(J)*W2*V2
     >       -JUM*0.5/DY(J)*W1*V1
     >           -0.5/DY(J)*DX3*( JUP*VIS2*(Q2(I,JP,K)-Q2(I,JP,KM))
     >                           -JUM*VIS1*(Q2(I,J ,K)-Q2(I,J ,KM)))

C--   M33W^N_Y
      APJ=JUP*(
     >     -0.5/DY(J)*       VIS2/H(JP)   
     >      +0.5/DY(J)*V2/H(JP)*DY(J)/2.0 
     >      )

      ACJ= 
     >      +0.5/DY(J)*(VIS2/H(JP)+VIS1/H(J))
     >      +0.5/DY(J)*(JUP*V2/H(JP)*DY(JP)/2.0
     >                 -JUM*V1/H(J)*DY(JM)/2.0) 
      AMJ=JUM*(
     >     -0.5/DY(J)*       VIS1/H(J )     
     >      -0.5/DY(J)*V1/H(J)*DY(J)/2.0 
     >      )

C--   M33W^N_Z
      VIS1 = EV(I,J,KM) + 1.0/RE
      VIS2 = EV(I,J,K ) + 1.0/RE

      W2=0.5*(Q3(I,J,KP)+Q3(I,J,K ))
      W1=0.5*(Q3(I,J,K )+Q3(I,J,KM))

      APK=
     >      -DX3Q*    VIS2 
     >      +DX3*W2/2.0
      ACK=
     >      +DX3Q*(VIS2+VIS1)
     >      +DX3*(W2/2.0-W1/2.0)
      AMK=
     >      -DX3Q*    VIS1 
     >      -DX3*W1/2.0

C---  M31U^N

      VIS1 = EVJ(I ,J,K) + 1.0/RE
      VIS2 = EVJ(IP,J,K) + 1.0/RE

      U1=0.5*(Q1(I ,J,K)+Q1(I ,J,KM))
      U2=0.5*(Q1(IP,J,K)+Q1(IP,J,KM))

      W1=1.0/G(I )*(0.5*DX(IM)*Q3(I ,J,K)+0.5*DX(I )*Q3(IM,J,K))
      W2=1.0/G(IP)*(0.5*DX(I )*Q3(IP,J,K)+0.5*DX(IP)*Q3(I ,J,K))

      RM31U_N=0.5/DX(I)*(IUP*W2*U2-IUM*W1*U1)
     >         -0.5/DX(I)*DX3*( IUP*VIS2*(Q1(IP,J,K)-Q1(IP,J,KM))
     >                         -IUM*VIS1*(Q1(I, J,K)-Q1(I ,J,KM)))

C--   M33W^N_X

      API=IUP*(
     >         -0.5/DX(I)*    VIS2/G(IP)  
     >         +0.5/DX(I)*U2/G(IP)*0.5*DX(I)
     >       )
      ACI=
     >      +0.5/DX(I)*(VIS2/G(IP)+VIS1/G(I))
     >      +0.5/DX(I)*(U2/G(IP)*0.5*DX(IP)-U1/G(I)*0.5*DX(IM))
      AMI=IUM*(
     >         -0.5/DX(I)*    VIS1/G(I)  
     >         -0.5/DX(I)*U1/G(I)*0.5*DX(I)
     >       )

      RM33W_N=APJ*Q3(I,JP,K)
     >       +ACJ*Q3(I,J ,K)
     >       +AMJ*Q3(I,JM,K)
     >       +API*Q3(IP,J,K)
     >       +ACI*Q3(I ,J,K)
     >       +AMI*Q3(IM,J,K)
     >       +APK*Q3(I,J,KP)
     >       +ACK*Q3(I,J,K )
     >       +AMK*Q3(I,J,KM)


      RDUH3(I,J,K)=RDUH3(I,J,K)
     >           -1./DT*Q3(I,J,K)
     >           -RM31U_N-RM32V_N-RM33W_N

 10   CONTINUE 

      RETURN
      END


C***************** RHSDP ***********************     
      SUBROUTINE RHSDP(RDP,QH1,QH2,QH3)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/BCON/UBC1(M1,M3,3),UBC2(M1,M3,3)            
     >           ,UBC3(M2,M3,3),UBC4(M2,M3,3)
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      REAL QH1(0:M1,0:M2,0:M3)
      REAL QH2(0:M1,0:M2,0:M3)
      REAL QH3(0:M1,0:M2,0:M3)
      REAL RDP(0:M1,0:M2,0:M3)

!$omp parallel do private(IP,IM,JP,JM,KP,KM,IUM,IUP,JUM,JUP,
!$omp&            DIVUH,CBC)
!$omp&           shared(RDP)

      DO 10 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
 
      DO 10 J=1,N2M
      JP=J+1
      JM=J-1
      JUM=J-JMU(J)
      JUP=JPA(J)-J

      DO 10 I=1,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMV(I)
      IUP=IPA(I)-I

      DIVUH=(IUP*QH1(IP,J ,K )-IUM*QH1(I,J,K))/DX(I)
     >     +(JUP*QH2(I ,JP,K )-JUM*QH2(I,J,K))/DY(J)
     >     +    (QH3(I ,J ,KP)-    QH3(I,J,K))*DX3

      CBC=(1-JUM)*UBC1(I,K,2)/DY(J)
     >   -(1-JUP)*UBC2(I,K,2)/DY(J)
     >   +(1-IUM)*UBC3(J,K,1)/DX(I)
     >   -(1-IUP)*UBC4(J,K,1)/DX(I)

      RDP(I,J,K)=(DIVUH-CBC)/DT

  10  CONTINUE

      RETURN 
      END


C***************** GETDUH1 ***********************     
      SUBROUTINE GETDUH1(Q1,Q2,Q3,QH1,QH2,QH3,RDUH1,EV,EVI,EVJ,EVK)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL QH1(0:M1,0:M2,0:M3)
      REAL QH2(0:M1,0:M2,0:M3)
      REAL QH3(0:M1,0:M2,0:M3)
      REAL RDUH1(0:M1,0:M2,0:M3)
      REAL EV(0:M1,0:M2,0:M3)
      REAL EVI(0:M1,0:M2,0:M3)   
      REAL EVJ(0:M1,0:M2,0:M3)   
      REAL EVK(0:M1,0:M2,0:M3)   

      REAL API(M1,M2),ACI(M1,M2),AMI(M1,M2)
      REAL APJ(M1,M2),ACJ(M1,M2),AMJ(M1,M2)
      REAL APK(M1,M3),ACK(M1,M3),AMK(M1,M3)
      REAL R1(M1,M2),R2(M1,M2),R3(M1,M3)
      EQUIVALENCE(API,APJ)
      EQUIVALENCE(ACI,ACJ)
      EQUIVALENCE(AMI,AMJ)
      EQUIVALENCE(R1,R2)

!$omp  parallel do default(shared)
!$omp& private(IP,IM,JP,JM,JUM,JUP,
!$omp&         VIS1,VIS2,V1,V2,APJ,ACJ,AMJ,
!$omp&         R2)

      DO 2 K=1,N3M
      DO 20 J=1,N2M
      JP=J+1
      JM=J-1
      JUM=J-JMU(J)
      JUP=JPA(J)-J
      DO 20 I=2,N1M
      IP=I+1
      IM=I-1

      VIS1 = EVK(I,J ,K) + 1.0/RE
      VIS2 = EVK(I,JP,K) + 1.0/RE

      V1=1.0/G(I)*(0.5*DX(IM)*Q2(I,J ,K)+0.5*DX(I)*Q2(IM,J ,K))
      V2=1.0/G(I)*(0.5*DX(IM)*Q2(I,JP,K)+0.5*DX(I)*Q2(IM,JP,K))

      APJ(I,J)=JUP*(
     >      -0.5/DY(J)*VIS2/H(JP)
     >      +0.5/DY(J)*V2/H(JP)*DY(J)/2.0
     >      )*DT

      ACJ(I,J)=( 0.5/DY(J)*(VIS2/H(JP)+VIS1/H(J))
     >    +0.5/DY(J)*(JUP*V2/H(JP)*DY(JP)/2.0
     >               -JUM*V1/H(J )*DY(JM)/2.0))*DT
     >    +1.0

      AMJ(I,J)=JUM*(
     >      -0.5/DY(J)*VIS1/H(J )
     >      -0.5/DY(J)*V1/H(J)*DY(J)/2.0
     >      )*DT

      R2(I,J)=RDUH1(I,J,K)*DT

  20  CONTINUE

      CALL TDMAI(AMJ,ACJ,APJ,R2,R2,1,N2M,2,N1M)
      DO 21 J=1,N2M
      DO 21 I=2,N1M
 21   QH1(I,J,K)=R2(I,J)
  2   CONTINUE


!$omp  parallel do default(shared)
!$omp& private(IP,IM,IUM,IUP,
!$omp&         VIS1,VIS2,U1,U2,API,ACI,AMI,
!$omp&         R1)

      DO 1 K=1,N3M
      DO 10 J=1,N2M
      DO 10 I=2,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMU(I)
      IUP=IPA(I)-I

      VIS1 = EV(IM,J,K) + 1.0/RE
      VIS2 = EV(I ,J,K) + 1.0/RE

      U2=0.5*(Q1(IP,J,K)+Q1(I, J,K))
      U1=0.5*(Q1(I, J,K)+Q1(IM,J,K)) 

      API(I,J)= (-1.0/G(I)*       VIS2/DX(I )      
     >      +1.0/G(I)*U2*0.5)*IUP
     >       *DT
      ACI(I,J)=(  1.0/G(I)*(VIS2/DX(I)+VIS1/DX(IM))
     >      +1.0/G(I)*(U2*0.5-U1*0.5))*DT
     >      +1.0
      AMI(I,J)= (-1.0/G(I)*       VIS1/DX(IM)  
     >      -1.0/G(I)*U1*0.5)*IUM
     >       *DT

      R1(I,J)=QH1(I,J,K) 
  10  CONTINUE

      CALL TDMAJ(AMI,ACI,API,R1,R1,2,N1M,1,N2M)
      DO 11 J=1,N2M
      DO 11 I=2,N1M
 11   QH1(I,J,K)=R1(I,J)
  1   CONTINUE

 
c!$omp  parallel do default(shared)
c!$omp& private(IP,IM,KP,KM,
c!$omp&         VIS1,VIS2,W1,W2,APK,ACK,AMK,
c!$omp&         R3)

      DO 3 J=1,N2M
      DO 30 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 30 I=2,N1M
      IP=I+1
      IM=I-1

      VIS1 = EVJ(I,J,K ) + 1.0/RE
      VIS2 = EVJ(I,J,KP) + 1.0/RE

      W1=1.0/G(I)*(0.5*DX(IM)*Q3(I,J,K )+0.5*DX(I)*Q3(IM,J,K ))
      W2=1.0/G(I)*(0.5*DX(IM)*Q3(I,J,KP)+0.5*DX(I)*Q3(IM,J,KP))

      APK(I,K)=(   -0.5*DX3Q*    VIS2  
     >       +0.5*DX3*W2*0.5)*DT
      ACK(I,K)=(   +0.5*DX3Q*(VIS2+VIS1)
     >       +0.5*DX3*(W2*0.5-W1*0.5))*DT
     >       +1.0
      AMK(I,K)=(   -0.5*DX3Q*    VIS1  
     >       -0.5*DX3*W1*0.5)*DT

   
      R3(I,K)=QH1(I,J,K)
  30  CONTINUE
c      CALL CTDMA3I(AMK,ACK,APK,R3,QH1,J,N3M,2,N1M)  <- why not accepted in OMP
      CALL TRIDP3 (AMK,ACK,APK,R3,1,N3M,2,N1M)
        DO K=1,N3M
        DO I=2,N1M
         QH1(I,J,K)=R3(I,K)
        ENDDO
        ENDDO
  3   CONTINUE

      RETURN
      END



C***************** GETDUH2 ***********************     
      SUBROUTINE GETDUH2(Q1,Q2,Q3,QH1,QH2,QH3,RDUH2,EV,EVI,EVJ,EVK)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL QH1(0:M1,0:M2,0:M3)
      REAL QH2(0:M1,0:M2,0:M3)
      REAL QH3(0:M1,0:M2,0:M3)
      REAL RDUH2(0:M1,0:M2,0:M3)
      REAL EV(0:M1,0:M2,0:M3)
      REAL EVI(0:M1,0:M2,0:M3)   
      REAL EVJ(0:M1,0:M2,0:M3)   
      REAL EVK(0:M1,0:M2,0:M3)   

      REAL API(M1,M2),ACI(M1,M2),AMI(M1,M2)
      REAL APJ(M1,M2),ACJ(M1,M2),AMJ(M1,M2)
      REAL APK(M1,M3),ACK(M1,M3),AMK(M1,M3)
      REAL R1(M1,M2),R2(M1,M2),R3(M1,M3)
      EQUIVALENCE(API,APJ)
      EQUIVALENCE(ACI,ACJ)
      EQUIVALENCE(AMI,AMJ)
      EQUIVALENCE(R1,R2)

!$omp  parallel do default(shared)
!$omp& private(IP,IM,JP,JM,KP,KM,IUM,IUP,JUM,JUP,
!$omp&         VIS1,VIS2,APJ,ACJ,AMJ,V1,V2,
!$omp$         UH1,UH2,RM21UH, 
!$omp&         R2)

      DO 2 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 20 J=2,N2M
      JP=J+1
      JM=J-1
      JUM=J-JMV(J)
      JUP=JPA(J)-J
      DO 20 I=1,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMV(I)
      IUP=IPA(I)-I

      VIS1 = EV(I,JM,K) + 1.0/RE
      VIS2 = EV(I,J ,K) + 1.0/RE

      V2=0.5*(Q2(I,JP,K)+Q2(I,J ,K))
      V1=0.5*(Q2(I,J ,K)+Q2(I,JM,K))

      APJ(I,J)=JUP*(
     >      -1.0/H(J)*       VIS2/DY(J ) 
     >      +1.0/H(J)*V2*0.5
     >      )*DT
      ACJ(I,J)=1.0+(
     >      +1.0/H(J)*(VIS2/DY(J)+VIS1/DY(JM))
     >      +1.0/H(J)*(V2*0.5-V1*0.5)
     >      )*DT
      AMJ(I,J)=JUM*(
     >      -1.0/H(J)*       VIS1/DY(JM) 
     >      -1.0/H(J)*V1*0.5
     >      )*DT

C     M21UH

      VIS1 = EVK(I ,J,K) + 1.0/RE
      VIS2 = EVK(IP,J,K) + 1.0/RE

      UH1=1.0/H(J)*(0.5*DY(JM)*QH1(I ,J,K)+0.5*DY(J)*QH1(I ,JM,K))
      UH2=1.0/H(J)*(0.5*DY(JM)*QH1(IP,J,K)+0.5*DY(J)*QH1(IP,JM,K))

      V1=1.0/G(I )*(0.5*DX(IM)*Q2(I ,J,K)+0.5*DX(I )*Q2(IM,J,K))
      V2=1.0/G(IP)*(0.5*DX(I )*Q2(IP,J,K)+0.5*DX(IP)*Q2(I ,J,K))

      RM21UH=0.5/DX(I)*(IUP*V2*UH2-IUM*V1*UH1)
     >
     >       -0.5/DX(I)/H(J)*( IUP*VIS2*(QH1(IP,J,K)-QH1(IP,JM,K))
     >                        -IUM*VIS1*(QH1(I, J,K)-QH1(I ,JM,K)))

      R2(I,J)=DT*(RDUH2(I,J,K)-RM21UH)
    
  20  CONTINUE

      CALL TDMAI(AMJ,ACJ,APJ,R2,R2,2,N2M,1,N1M)
      DO 21 J=2,N2M
      DO 21 I=1,N1M
 21   QH2(I,J,K)=R2(I,J)
  2   CONTINUE


!$omp  parallel do default(shared)
!$omp& private(IP,IM,JP,JM,KP,KM,IUM,IUP,
!$omp&         VIS1,VIS2,API,ACI,AMI,U1,U2,
!$omp&         R1)

      DO 1 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 10 J=2,N2M
      JP=J+1
      JM=J-1
      DO 10 I=1,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMV(I)
      IUP=IPA(I)-I

      VIS1 = EVK(I ,J,K) + 1.0/RE
      VIS2 = EVK(IP,J,K) + 1.0/RE

      U1=1.0/H(J)*(0.5*DY(JM)*Q1(I ,J,K)+0.5*DY(J)*Q1(I ,JM,K))
      U2=1.0/H(J)*(0.5*DY(JM)*Q1(IP,J,K)+0.5*DY(J)*Q1(IP,JM,K))

      API(I,J)=IUP*(
     >         -0.5/DX(I)*    VIS2/G(IP) 
     >         +0.5/DX(I)*U2/G(IP)*0.5*DX(I)
     >       )*DT
      ACI(I,J)=(  
     >     +0.5/DX(I)*(VIS2/G(IP)+VIS1/G(I))
     >     +0.5/DX(I)*(U2/G(IP)*0.5*DX(IP)-U1/G(I)*0.5*DX(IM)))*DT
     >     +1.0
      AMI(I,J)=IUM*(
     >         -0.5/DX(I)*    VIS1/G(I) 
     >         -0.5/DX(I)*U1/G(I)*0.5*DX(I)
     >       )*DT

      R1(I,J)=QH2(I,J,K)

  10  CONTINUE
      CALL TDMAJ(AMI,ACI,API,R1,R1,1,N1M,2,N2M)
      DO 11 J=2,N2M
      DO 11 I=1,N1M
 11   QH2(I,J,K)=R1(I,J)
  1   CONTINUE


c!$omp  parallel do default(shared)
c!$omp& private(IP,IM,JP,JM,KP,KM,
c!$omp&         VIS1,VIS2,APK,ACK,AMK,W1,W2,
c!$omp&         R3)

      DO 3 J=2,N2M
      JP=J+1
      JM=J-1
      DO 30 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 30 I=1,N1M
      IP=I+1
      IM=I-1

      VIS1 = EVI(I,J,K ) + 1.0/RE
      VIS2 = EVI(I,J,KP) + 1.0/RE

      W2=1.0/H(J)*(DY(J)/2.0*Q3(I,JM,KP)+DY(JM)/2.0*Q3(I,J,KP))
      W1=1.0/H(J)*(DY(J)/2.0*Q3(I,JM,K )+DY(JM)/2.0*Q3(I,J,K ))

      APK(I,K)=(
     >      - 0.5*DX3Q*    VIS2  
     >      +0.5*DX3*W2/2.0
     >       )*DT
      ACK(I,K)=1.+(
     >      +0.5*DX3Q*(VIS2+VIS1)
     >      +0.5*DX3*(W2/2.0-W1/2.0)
     >       )*DT
      AMK(I,K)=(
     >      - 0.5*DX3Q*    VIS1  
     >      -0.5*DX3*W1/2.0
     >       )*DT

      R3(I,K)=QH2(I,J,K)
  30  CONTINUE
c      CALL CTDMA3I(AMK,ACK,APK,R3,QH2,J,N3M,1,N1M)
      CALL TRIDP3 (AMK,ACK,APK,R3,1,N3M,1,N1M)
        DO K=1,N3M
        DO I=1,N1M
          QH2(I,J,K)=R3(I,K)
        ENDDO
        ENDDO

  3   CONTINUE

      RETURN
      END



C***************** GETDUH3 ***********************     
      SUBROUTINE GETDUH3(Q1,Q2,Q3,QH1,QH2,QH3,RDUH3,EV,EVI,EVJ,EVK)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL QH1(0:M1,0:M2,0:M3)
      REAL QH2(0:M1,0:M2,0:M3)
      REAL QH3(0:M1,0:M2,0:M3)
      REAL RDUH3(0:M1,0:M2,0:M3)
      REAL EV(0:M1,0:M2,0:M3)
      REAL EVI(0:M1,0:M2,0:M3)   
      REAL EVJ(0:M1,0:M2,0:M3)   
      REAL EVK(0:M1,0:M2,0:M3)   

      REAL API(M1,M2),ACI(M1,M2),AMI(M1,M2)
      REAL APJ(M1,M2),ACJ(M1,M2),AMJ(M1,M2)
      REAL APK(M1,M3),ACK(M1,M3),AMK(M1,M3)
      REAL R1(M1,M2),R2(M1,M2),R3(M1,M3)
      EQUIVALENCE(API,APJ)
      EQUIVALENCE(ACI,ACJ)
      EQUIVALENCE(AMI,AMJ)
      EQUIVALENCE(R1,R2)

!$omp  parallel do default(shared)
!$omp& private(IP,IM,JP,JM,KP,KM,IUM,IUP,JUM,JUP,
!$omp&         VIS1,VIS2,V1,V2,APJ,ACJ,AMJ,W1,W2,VH1,VH2,RM32VH,
!$omp&         UH1,UH2,RM31UH,R2)

      DO 2 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 20 J=1,N2M
      JP=J+1
      JM=J-1
      JUM=J-JMU(J)
      JUP=JPA(J)-J
      DO 20 I=1,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMV(I)
      IUP=IPA(I)-I

      VIS1 = EVI(I,J ,K) + 1.0/RE
      VIS2 = EVI(I,JP,K) + 1.0/RE

      V2=0.5*(Q2(I,JP,K)+Q2(I,JP,KM))
      V1=0.5*(Q2(I,J ,K)+Q2(I,J ,KM))

      APJ(I,J)=JUP*(
     >      -0.5/DY(J)*       VIS2/H(JP)   
     >      +0.5/DY(J)*V2/H(JP)*DY(J)/2.0 
     >      )*DT

      ACJ(I,J)= (
     >      +0.5/DY(J)*(VIS2/H(JP)+VIS1/H(J))
     >      +0.5/DY(J)*(JUP*V2/H(JP)*DY(JP)/2.0
     >                 -JUM*V1/H(J)*DY(JM)/2.0) )
     >     *DT +1.0
      AMJ(I,J)=JUM*(
     >      -0.5/DY(J)*       VIS1/H(J )     
     >      -0.5/DY(J)*V1/H(J)*DY(J)/2.0 
     >      )*DT

C     M32VH
      W2=1.0/H(JP)*(DY(J )/2.*Q3(I,JP,K)+DY(JP)/2.*Q3(I,J ,K))
      W1=1.0/H(J )*(DY(JM)/2.*Q3(I,J ,K)+DY(J )/2.*Q3(I,JM,K))
      VH2=0.5*(QH2(I,JP,K)+QH2(I,JP,KM))      
      VH1=0.5*(QH2(I,J ,K)+QH2(I,J ,KM))      
      RM32VH=JUP*0.5/DY(J)*W2*VH2
     >      -JUM*0.5/DY(J)*W1*VH1
     >           -0.5/DY(J)*DX3*( JUP*VIS2*(QH2(I,JP,K)-QH2(I,JP,KM))
     >                           -JUM*VIS1*(QH2(I,J, K)-QH2(I,J ,KM)))

C     M31UH

      VIS1 = EVJ(I ,J,K) + 1.0/RE
      VIS2 = EVJ(IP,J,K) + 1.0/RE

      UH1=0.5*(QH1(I ,J,K)+QH1(I ,J,KM))
      UH2=0.5*(QH1(IP,J,K)+QH1(IP,J,KM))

      W1=1.0/G(I )*(0.5*DX(IM)*Q3(I ,J,K)+0.5*DX(I )*Q3(IM,J,K))
      W2=1.0/G(IP)*(0.5*DX(I )*Q3(IP,J,K)+0.5*DX(IP)*Q3(I ,J,K))

      RM31UH=0.5/DX(I)*(IUP*W2*UH2-IUM*W1*UH1)
     >         -0.5/DX(I)*DX3*( IUP*VIS2*(QH1(IP,J,K)-QH1(IP,J,KM))
     >                         -IUM*VIS1*(QH1(I, J,K)-QH1(I ,J,KM)))

      R2(I,J)=DT*(RDUH3(I,J,K)-RM31UH-RM32VH)

  20  CONTINUE
      CALL TDMAI(AMJ,ACJ,APJ,R2,R2,1,N2M,1,N1M)
      DO 21 J=1,N2M
      DO 21 I=1,N1M
 21   QH3(I,J,K)=R2(I,J)
  2   CONTINUE


c!$omp  parallel do default(shared)
c!$omp& private(IP,IM,JP,JM,KP,KM,
c!$omp&         VIS1,VIS2,APK,ACK,AMK,W1,W2,
c!$omp&         R3)

      DO 3 J=1,N2M
      JP=J+1
      JM=J-1
      DO 30 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 30 I=1,N1M
      IP=I+1
      IM=I-1

      VIS1 = EV(I,J,KM) + 1.0/RE
      VIS2 = EV(I,J,K ) + 1.0/RE

      W2=0.5*(Q3(I,J,KP)+Q3(I,J,K ))
      W1=0.5*(Q3(I,J,K )+Q3(I,J,KM))

      APK(I,K)=(
     >      -DX3Q*    VIS2 
     >      +DX3*W2/2.0
     >       )*DT
      ACK(I,K)=1.+(
     >      +DX3Q*(VIS2+VIS1)
     >      +DX3*(W2/2.0-W1/2.0)
     >      )*DT
      AMK(I,K)=(
     >      -DX3Q*    VIS1 
     >      -DX3*W1/2.0
     >       )*DT

      R3(I,K)=QH3(I,J,K)
  30  CONTINUE
c      CALL CTDMA3I(AMK,ACK,APK,R3,QH3,J,N3M,1,N1M) 
      CALL TRIDP3 (AMK,ACK,APK,R3,1,N3M,1,N1M)
        DO K=1,N3M
        DO I=1,N1M
         QH3(I,J,K)=R3(I,K)
        ENDDO
        ENDDO
  3   CONTINUE


!$omp  parallel do default(shared)
!$omp& private(IP,IM,KP,KM,IUM,IUP,
!$omp&         VIS1,VIS2,API,ACI,AMI,U1,U2,
!$omp&         R1)

      DO 1 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 10 J=1,N2M
      DO 10 I=1,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMV(I)
      IUP=IPA(I)-I

      VIS1 = EVJ(I ,J,K) + 1.0/RE
      VIS2 = EVJ(IP,J,K) + 1.0/RE

      U1=0.5*(Q1(I ,J,K)+Q1(I ,J,KM))
      U2=0.5*(Q1(IP,J,K)+Q1(IP,J,KM))

      API(I,J)=IUP*(
     >         -0.5/DX(I)*    VIS2/G(IP)  
     >         +0.5/DX(I)*U2/G(IP)*0.5*DX(I)
     >       )*DT
      ACI(I,J)=(
     >      +0.5/DX(I)*(VIS2/G(IP)+VIS1/G(I))
     >      +0.5/DX(I)*(U2/G(IP)*0.5*DX(IP)-U1/G(I)*0.5*DX(IM)))*DT
     >      +1.0
      AMI(I,J)=IUM*(
     >         -0.5/DX(I)*    VIS1/G(I)  
     >         -0.5/DX(I)*U1/G(I)*0.5*DX(I)
     >       )*DT

      R1(I,J)=QH3(I,J,K)
  10  CONTINUE
      CALL TDMAJ(AMI,ACI,API,R1,R1,1,N1M,1,N2M)
      DO 11 J=1,N2M
      DO 11 I=1,N1M
 11   QH3(I,J,K)=R1(I,J)
  1   CONTINUE



C     DVH UPDATE        

!$omp  parallel do default(shared)
!$omp& private(IP,IM,JP,JM,KP,KM,IUM,IUP,JUM,JUP,
!$omp&         VIS1,VIS2,V1,V2,WH1,WH2,RM23WH)

      DO 110 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)

      DO 110 J=2,N2M
      JP=J+1
      JM=J-1
      JUM=J-JMV(J)
      JUP=JPA(J)-J

      DO 110 I=1,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMV(I)
      IUP=IPA(I)-I

C     M23WH

      VIS1 = EVI(I,J,K ) + 1.0/RE
      VIS2 = EVI(I,J,KP) + 1.0/RE

      V2=0.5*(Q2(I,J,KP)+Q2(I,J,K ))
      V1=0.5*(Q2(I,J,K )+Q2(I,J,KM))
      WH2=1.0/H(J)*(DY(J)/2.*QH3(I,JM,KP)+DY(JM)/2.*QH3(I,J,KP))
      WH1=1.0/H(J)*(DY(J)/2.*QH3(I,JM,K )+DY(JM)/2.*QH3(I,J,K ))
      RM23WH=0.5*DX3*(V2*WH2-V1*WH1)
     >        
     >       -0.5*DX3/H(J)*(VIS2*(QH3(I,J,KP)-QH3(I,JM,KP))
     >                     -VIS1*(QH3(I,J,K )-QH3(I,JM,K )))
      
      QH2(I,J,K)=QH2(I,J,K)-DT*RM23WH 
110   CONTINUE


C     DUH UPDATE

!$omp  parallel do default(shared)
!$omp& private(IP,IM,JP,JM,KP,KM,IUM,IUP,JUM,JUP,
!$omp&         VIS1,VIS2,U1,U2,VH1,VH2,RM12VH,
!$omp&         WH1,WH2,RM13WH )

      DO 210 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)

      DO 210 J=1,N2M
      JP=J+1
      JM=J-1
      JUM=J-JMU(J)
      JUP=JPA(J)-J
      DO 210 I=2,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMU(I)
      IUP=IPA(I)-I

C     M12VH

      VIS1 = EVK(I,J ,K) + 1.0/RE
      VIS2 = EVK(I,JP,K) + 1.0/RE

      VH1=1.0/G(I)*(0.5*DX(IM)*QH2(I,J ,K)+0.5*DX(I)*QH2(IM,J ,K))
      VH2=1.0/G(I)*(0.5*DX(IM)*QH2(I,JP,K)+0.5*DX(I)*QH2(IM,JP,K))

      U1=1.0/H(J )*(0.5*DY(JM)*Q1(I,J ,K)+0.5*DY(J )*Q1(I,JM,K))
      U2=1.0/H(JP)*(0.5*DY(J )*Q1(I,JP,K)+0.5*DY(JP)*Q1(I,J ,K))

      RM12VH=JUP*0.5/DY(J)*U2*VH2
     >       -JUM*0.5/DY(J)*U1*VH1
     >
     >           -0.5/DY(J)/G(I)*( JUP*VIS2*(QH2(I,JP,K)-QH2(IM,JP,K))
     >                            -JUM*VIS1*(QH2(I,J ,K)-QH2(IM,J ,K)))

C     M13WH

      VIS1 = EVJ(I,J,K ) + 1.0/RE
      VIS2 = EVJ(I,J,KP) + 1.0/RE

      WH1=1.0/G(I)*(0.5*DX(IM)*QH3(I,J,K )+0.5*DX(I)*QH3(IM,J,K ))
      WH2=1.0/G(I)*(0.5*DX(IM)*QH3(I,J,KP)+0.5*DX(I)*QH3(IM,J,KP))

      U1=0.5*(Q1(I,J,K )+Q1(I,J,KM))
      U2=0.5*(Q1(I,J,KP)+Q1(I,J,K ))

      RM13WH=0.5*DX3*(U2*WH2-U1*WH1)
     >
     >           -0.5*DX3/G(I)*(VIS2*(QH3(I,J,KP)-QH3(IM,J,KP))
     >                         -VIS1*(QH3(I,J,K )-QH3(IM,J,K )))

      QH1(I,J,K)=QH1(I,J,K)-DT*(RM12VH+RM13WH)
210   CONTINUE

      RETURN
      END


C***************** GETDUH1_I ***********************     
      SUBROUTINE GETDUH1_I(Q1,Q2,Q3,QH1,QH2,QH3,RDUH1,EV,EVI,EVJ,EVK)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL QH1(0:M1,0:M2,0:M3)
      REAL QH2(0:M1,0:M2,0:M3)
      REAL QH3(0:M1,0:M2,0:M3)
      REAL RDUH1(0:M1,0:M2,0:M3)
      REAL EV(0:M1,0:M2,0:M3)
      REAL EVI(0:M1,0:M2,0:M3)   
      REAL EVJ(0:M1,0:M2,0:M3)   
      REAL EVK(0:M1,0:M2,0:M3)   

      REAL API(M1,M2),ACI(M1,M2),AMI(M1,M2)
      REAL APJ(M1,M2),ACJ(M1,M2),AMJ(M1,M2)
      REAL APK(M1,M3),ACK(M1,M3),AMK(M1,M3)
      REAL R1(M1,M2),R2(M1,M2),R3(M1,M3)
      EQUIVALENCE(API,APJ)
      EQUIVALENCE(ACI,ACJ)
      EQUIVALENCE(AMI,AMJ)
      EQUIVALENCE(R1,R2)


      DO 2 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 20 J=1,N2M
      JP=J+1
      JM=J-1
      JUM=J-JMU(J)
      JUP=JPA(J)-J
      DO 20 I=2,N1M
      IP=I+1
      IM=I-1

      VIS1 = EVK(I,J ,K) + 1.0/RE
      VIS2 = EVK(I,JP,K) + 1.0/RE

      V1=1.0/G(I)*(0.5*DX(IM)*Q2(I,J ,K)+0.5*DX(I)*Q2(IM,J ,K))
      V2=1.0/G(I)*(0.5*DX(IM)*Q2(I,JP,K)+0.5*DX(I)*Q2(IM,JP,K))

      APJ(I,J)=JUP*(
     >      -0.5/DY(J)*VIS2/H(JP)
     >      +0.5/DY(J)*V2/H(JP)*DY(J)/2.0
     >      )*DT

      ACJ(I,J)=( 0.5/DY(J)*(VIS2/H(JP)+VIS1/H(J))
     >    +0.5/DY(J)*(JUP*V2/H(JP)*DY(JP)/2.0
     >               -JUM*V1/H(J )*DY(JM)/2.0))*DT
     >    +1.0

      AMJ(I,J)=JUM*(
     >      -0.5/DY(J)*VIS1/H(J )
     >      -0.5/DY(J)*V1/H(J)*DY(J)/2.0
     >      )*DT

C     M12VH

      VIS1 = EVK(I,J ,K) + 1.0/RE
      VIS2 = EVK(I,JP,K) + 1.0/RE

      VH1=1.0/G(I)*(0.5*DX(IM)*QH2(I,J ,K)+0.5*DX(I)*QH2(IM,J ,K))
      VH2=1.0/G(I)*(0.5*DX(IM)*QH2(I,JP,K)+0.5*DX(I)*QH2(IM,JP,K))

      U1=1.0/H(J )*(0.5*DY(JM)*Q1(I,J ,K)+0.5*DY(J )*Q1(I,JM,K))
      U2=1.0/H(JP)*(0.5*DY(J )*Q1(I,JP,K)+0.5*DY(JP)*Q1(I,J ,K))

      RM12VH=JUP*0.5/DY(J)*U2*VH2
     >       -JUM*0.5/DY(J)*U1*VH1
     >
     >           -0.5/DY(J)/G(I)*( JUP*VIS2*(QH2(I,JP,K)-QH2(IM,JP,K))
     >                            -JUM*VIS1*(QH2(I,J ,K)-QH2(IM,J ,K)))


C     M13WH

      VIS1 = EVJ(I,J,K ) + 1.0/RE
      VIS2 = EVJ(I,J,KP) + 1.0/RE

      WH1=1.0/G(I)*(0.5*DX(IM)*QH3(I,J,K )+0.5*DX(I)*QH3(IM,J,K ))
      WH2=1.0/G(I)*(0.5*DX(IM)*QH3(I,J,KP)+0.5*DX(I)*QH3(IM,J,KP))

      U1=0.5*(Q1(I,J,K )+Q1(I,J,KM))
      U2=0.5*(Q1(I,J,KP)+Q1(I,J,K ))

      RM13WH=0.5*DX3*(U2*WH2-U1*WH1)
     >
     >           -0.5*DX3/G(I)*(VIS2*(QH3(I,J,KP)-QH3(IM,J,KP))
     >                         -VIS1*(QH3(I,J,K )-QH3(IM,J,K )))
   
      R2(I,J)=(RDUH1(I,J,K)-RM12VH-RM13WH)*DT
  20  CONTINUE

      CALL TDMAI(AMJ,ACJ,APJ,R2,R2,1,N2M,2,N1M)
      DO 21 J=1,N2M
      DO 21 I=2,N1M
 21   QH1(I,J,K)=R2(I,J)
  2   CONTINUE


      DO 1 K=1,N3M
      DO 10 J=1,N2M
      DO 10 I=2,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMU(I)
      IUP=IPA(I)-I

      VIS1 = EV(IM,J,K) + 1.0/RE
      VIS2 = EV(I ,J,K) + 1.0/RE

      U2=0.5*(Q1(IP,J,K)+Q1(I, J,K))
      U1=0.5*(Q1(I, J,K)+Q1(IM,J,K)) 

      API(I,J)= (-1.0/G(I)*       VIS2/DX(I )      
     >      +1.0/G(I)*U2*0.5)*IUP
     >       *DT
      ACI(I,J)=(  1.0/G(I)*(VIS2/DX(I)+VIS1/DX(IM))
     >      +1.0/G(I)*(U2*0.5-U1*0.5))*DT
     >      +1.0
      AMI(I,J)= (-1.0/G(I)*       VIS1/DX(IM)  
     >      -1.0/G(I)*U1*0.5)*IUM
     >       *DT

      R1(I,J)=QH1(I,J,K) 
  10  CONTINUE

      CALL TDMAJ(AMI,ACI,API,R1,R1,2,N1M,1,N2M)
      DO 11 J=1,N2M
      DO 11 I=2,N1M
 11   QH1(I,J,K)=R1(I,J)
  1   CONTINUE



      DO 3 J=1,N2M
      DO 30 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 30 I=2,N1M
      IP=I+1
      IM=I-1

      VIS1 = EVJ(I,J,K ) + 1.0/RE
      VIS2 = EVJ(I,J,KP) + 1.0/RE

      W1=1.0/G(I)*(0.5*DX(IM)*Q3(I,J,K )+0.5*DX(I)*Q3(IM,J,K ))
      W2=1.0/G(I)*(0.5*DX(IM)*Q3(I,J,KP)+0.5*DX(I)*Q3(IM,J,KP))

      APK(I,K)=(   -0.5*DX3Q*    VIS2  
     >       +0.5*DX3*W2*0.5)*DT
      ACK(I,K)=(   +0.5*DX3Q*(VIS2+VIS1)
     >       +0.5*DX3*(W2*0.5-W1*0.5))*DT
     >       +1.0
      AMK(I,K)=(   -0.5*DX3Q*    VIS1  
     >       -0.5*DX3*W1*0.5)*DT

      R3(I,K)=QH1(I,J,K)
  30  CONTINUE
c      CALL CTDMA3I(AMK,ACK,APK,R3,QH1,J,N3M,2,N1M)
      CALL TRIDP3 (AMK,ACK,APK,R3,1,N3M,2,N1M)
        DO K=1,N3M
        DO I=2,N1M
         QH1(I,J,K)=R3(I,K)
        ENDDO
        ENDDO
  3   CONTINUE

      RETURN
      END

C***************** GETDUH2_I ***********************     
      SUBROUTINE GETDUH2_I(Q1,Q2,Q3,QH1,QH2,QH3,RDUH2,EV,EVI,EVJ,EVK)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL QH1(0:M1,0:M2,0:M3)
      REAL QH2(0:M1,0:M2,0:M3)
      REAL QH3(0:M1,0:M2,0:M3)
      REAL RDUH2(0:M1,0:M2,0:M3)
      REAL EV(0:M1,0:M2,0:M3)
      REAL EVI(0:M1,0:M2,0:M3)   
      REAL EVJ(0:M1,0:M2,0:M3)   
      REAL EVK(0:M1,0:M2,0:M3)   

      REAL API(M1,M2),ACI(M1,M2),AMI(M1,M2)
      REAL APJ(M1,M2),ACJ(M1,M2),AMJ(M1,M2)
      REAL APK(M1,M3),ACK(M1,M3),AMK(M1,M3)
      REAL R1(M1,M2),R2(M1,M2),R3(M1,M3)
      EQUIVALENCE(API,APJ)
      EQUIVALENCE(ACI,ACJ)
      EQUIVALENCE(AMI,AMJ)
      EQUIVALENCE(R1,R2)


      DO 2 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 20 J=2,N2M
      JP=J+1
      JM=J-1
      JUM=J-JMV(J)
      JUP=JPA(J)-J
      DO 20 I=1,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMV(I)
      IUP=IPA(I)-I

      VIS1 = EV(I,JM,K) + 1.0/RE
      VIS2 = EV(I,J ,K) + 1.0/RE

      V2=0.5*(Q2(I,JP,K)+Q2(I,J ,K))
      V1=0.5*(Q2(I,J ,K)+Q2(I,JM,K))

      APJ(I,J)=JUP*(
     >      -1.0/H(J)*       VIS2/DY(J ) 
     >      +1.0/H(J)*V2*0.5
     >      )*DT
      ACJ(I,J)=1.0+(
     >      +1.0/H(J)*(VIS2/DY(J)+VIS1/DY(JM))
     >      +1.0/H(J)*(V2*0.5-V1*0.5)
     >      )*DT
      AMJ(I,J)=JUM*(
     >      -1.0/H(J)*       VIS1/DY(JM) 
     >      -1.0/H(J)*V1*0.5
     >      )*DT

C     M21UH

      VIS1 = EVK(I ,J,K) + 1.0/RE
      VIS2 = EVK(IP,J,K) + 1.0/RE

      UH1=1.0/H(J)*(0.5*DY(JM)*QH1(I ,J,K)+0.5*DY(J)*QH1(I ,JM,K))
      UH2=1.0/H(J)*(0.5*DY(JM)*QH1(IP,J,K)+0.5*DY(J)*QH1(IP,JM,K))

      V1=1.0/G(I )*(0.5*DX(IM)*Q2(I ,J,K)+0.5*DX(I )*Q2(IM,J,K))
      V2=1.0/G(IP)*(0.5*DX(I )*Q2(IP,J,K)+0.5*DX(IP)*Q2(I ,J,K))

      RM21UH=0.5/DX(I)*(IUP*V2*UH2-IUM*V1*UH1)
     >
     >       -0.5/DX(I)/H(J)*( IUP*VIS2*(QH1(IP,J,K)-QH1(IP,JM,K))
     >                        -IUM*VIS1*(QH1(I, J,K)-QH1(I ,JM,K)))

C     M23WH
 
      VIS1 = EVI(I,J,K ) + 1.0/RE
      VIS2 = EVI(I,J,KP) + 1.0/RE

      V2=0.5*(Q2(I,J,KP)+Q2(I,J,K ))
      V1=0.5*(Q2(I,J,K )+Q2(I,J,KM))
      WH2=1.0/H(J)*(DY(J)/2.*QH3(I,JM,KP)+DY(JM)/2.*QH3(I,J,KP))
      WH1=1.0/H(J)*(DY(J)/2.*QH3(I,JM,K )+DY(JM)/2.*QH3(I,J,K ))
      RM23WH=0.5*DX3*(V2*WH2-V1*WH1)
     >        
     >       -0.5*DX3/H(J)*(VIS2*(QH3(I,J,KP)-QH3(I,JM,KP))
     >                     -VIS1*(QH3(I,J,K )-QH3(I,JM,K )))
     
      R2(I,J)=DT*(RDUH2(I,J,K)-RM21UH-RM23WH)
    
  20  CONTINUE

      CALL TDMAI(AMJ,ACJ,APJ,R2,R2,2,N2M,1,N1M)
      DO 21 J=2,N2M
      DO 21 I=1,N1M
 21   QH2(I,J,K)=R2(I,J)
  2   CONTINUE



      DO 1 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 10 J=2,N2M
      JP=J+1
      JM=J-1
      DO 10 I=1,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMV(I)
      IUP=IPA(I)-I

      VIS1 = EVK(I ,J,K) + 1.0/RE
      VIS2 = EVK(IP,J,K) + 1.0/RE

      U1=1.0/H(J)*(0.5*DY(JM)*Q1(I ,J,K)+0.5*DY(J)*Q1(I ,JM,K))
      U2=1.0/H(J)*(0.5*DY(JM)*Q1(IP,J,K)+0.5*DY(J)*Q1(IP,JM,K))

      API(I,J)=IUP*(
     >         -0.5/DX(I)*    VIS2/G(IP) 
     >         +0.5/DX(I)*U2/G(IP)*0.5*DX(I)
     >       )*DT
      ACI(I,J)=(  
     >     +0.5/DX(I)*(VIS2/G(IP)+VIS1/G(I))
     >     +0.5/DX(I)*(U2/G(IP)*0.5*DX(IP)-U1/G(I)*0.5*DX(IM)))*DT
     >     +1.0
      AMI(I,J)=IUM*(
     >         -0.5/DX(I)*    VIS1/G(I) 
     >         -0.5/DX(I)*U1/G(I)*0.5*DX(I)
     >       )*DT

      R1(I,J)=QH2(I,J,K)

  10  CONTINUE
      CALL TDMAJ(AMI,ACI,API,R1,R1,1,N1M,2,N2M)
      DO 11 J=2,N2M
      DO 11 I=1,N1M
 11   QH2(I,J,K)=R1(I,J)
  1   CONTINUE



      DO 3 J=2,N2M
      JP=J+1
      JM=J-1
      DO 30 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 30 I=1,N1M
      IP=I+1
      IM=I-1

      VIS1 = EVI(I,J,K ) + 1.0/RE
      VIS2 = EVI(I,J,KP) + 1.0/RE

      W2=1.0/H(J)*(DY(J)/2.0*Q3(I,JM,KP)+DY(JM)/2.0*Q3(I,J,KP))
      W1=1.0/H(J)*(DY(J)/2.0*Q3(I,JM,K )+DY(JM)/2.0*Q3(I,J,K ))

      APK(I,K)=(
     >      - 0.5*DX3Q*    VIS2  
     >      +0.5*DX3*W2/2.0
     >       )*DT
      ACK(I,K)=1.+(
     >      +0.5*DX3Q*(VIS2+VIS1)
     >      +0.5*DX3*(W2/2.0-W1/2.0)
     >       )*DT
      AMK(I,K)=(
     >      - 0.5*DX3Q*    VIS1  
     >      -0.5*DX3*W1/2.0
     >       )*DT

      R3(I,K)=QH2(I,J,K)
  30  CONTINUE
c      CALL CTDMA3I(AMK,ACK,APK,R3,QH2,J,N3M,1,N1M)
      CALL TRIDP3 (AMK,ACK,APK,R3,1,N3M,1,N1M)
        DO K=1,N3M
        DO I=1,N1M
         QH2(I,J,K)=R3(I,K)
        ENDDO
        ENDDO
  3   CONTINUE

      RETURN
      END

C***************** GETDUH3_I ***********************     
      SUBROUTINE GETDUH3_I(Q1,Q2,Q3,QH1,QH2,QH3,RDUH3,EV,EVI,EVJ,EVK)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL QH1(0:M1,0:M2,0:M3)
      REAL QH2(0:M1,0:M2,0:M3)
      REAL QH3(0:M1,0:M2,0:M3)
      REAL RDUH3(0:M1,0:M2,0:M3)
      REAL EV(0:M1,0:M2,0:M3)
      REAL EVI(0:M1,0:M2,0:M3)   
      REAL EVJ(0:M1,0:M2,0:M3)   
      REAL EVK(0:M1,0:M2,0:M3)   

      REAL API(M1,M2),ACI(M1,M2),AMI(M1,M2)
      REAL APJ(M1,M2),ACJ(M1,M2),AMJ(M1,M2)
      REAL APK(M1,M3),ACK(M1,M3),AMK(M1,M3)
      REAL R1(M1,M2),R2(M1,M2),R3(M1,M3)
      EQUIVALENCE(API,APJ)
      EQUIVALENCE(ACI,ACJ)
      EQUIVALENCE(AMI,AMJ)
      EQUIVALENCE(R1,R2)


      DO 2 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 20 J=1,N2M
      JP=J+1
      JM=J-1
      JUM=J-JMU(J)
      JUP=JPA(J)-J
      DO 20 I=1,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMV(I)
      IUP=IPA(I)-I

      VIS1 = EVI(I,J ,K) + 1.0/RE
      VIS2 = EVI(I,JP,K) + 1.0/RE

      V2=0.5*(Q2(I,JP,K)+Q2(I,JP,KM))
      V1=0.5*(Q2(I,J ,K)+Q2(I,J ,KM))


      APJ(I,J)=JUP*(
     >      -0.5/DY(J)*       VIS2/H(JP)   
     >      +0.5/DY(J)*V2/H(JP)*DY(J)/2.0 
     >      )*DT

      ACJ(I,J)= (
     >      +0.5/DY(J)*(VIS2/H(JP)+VIS1/H(J))
     >      +0.5/DY(J)*(JUP*V2/H(JP)*DY(JP)/2.0
     >                 -JUM*V1/H(J)*DY(JM)/2.0) )
     >     *DT +1.0
      AMJ(I,J)=JUM*(
     >      -0.5/DY(J)*       VIS1/H(J )     
     >      -0.5/DY(J)*V1/H(J)*DY(J)/2.0 
     >      )*DT

C     M32VH
      W2=1.0/H(JP)*(DY(J )/2.*Q3(I,JP,K)+DY(JP)/2.*Q3(I,J ,K))
      W1=1.0/H(J )*(DY(JM)/2.*Q3(I,J ,K)+DY(J )/2.*Q3(I,JM,K))
      VH2=0.5*(QH2(I,JP,K)+QH2(I,JP,KM))      
      VH1=0.5*(QH2(I,J ,K)+QH2(I,J ,KM))      
      RM32VH=JUP*0.5/DY(J)*W2*VH2
     >      -JUM*0.5/DY(J)*W1*VH1
     >           -0.5/DY(J)*DX3*( JUP*VIS2*(QH2(I,JP,K)-QH2(I,JP,KM))
     >                           -JUM*VIS1*(QH2(I,J, K)-QH2(I,J ,KM)))

C     M31UH

      VIS1 = EVJ(I ,J,K) + 1.0/RE
      VIS2 = EVJ(IP,J,K) + 1.0/RE

      UH1=0.5*(QH1(I ,J,K)+QH1(I ,J,KM))
      UH2=0.5*(QH1(IP,J,K)+QH1(IP,J,KM))

      W1=1.0/G(I )*(0.5*DX(IM)*Q3(I ,J,K)+0.5*DX(I )*Q3(IM,J,K))
      W2=1.0/G(IP)*(0.5*DX(I )*Q3(IP,J,K)+0.5*DX(IP)*Q3(I ,J,K))

      RM31UH=0.5/DX(I)*(IUP*W2*UH2-IUM*W1*UH1)
     >         -0.5/DX(I)*DX3*( IUP*VIS2*(QH1(IP,J,K)-QH1(IP,J,KM))
     >                         -IUM*VIS1*(QH1(I, J,K)-QH1(I ,J,KM)))

      R2(I,J)=DT*(RDUH3(I,J,K)-RM31UH-RM32VH)

  20  CONTINUE
      CALL TDMAI(AMJ,ACJ,APJ,R2,R2,1,N2M,1,N1M)
      DO 21 J=1,N2M
      DO 21 I=1,N1M
 21   QH3(I,J,K)=R2(I,J)
  2   CONTINUE



      DO 3 J=1,N2M
      JP=J+1
      JM=J-1
      DO 30 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 30 I=1,N1M
      IP=I+1
      IM=I-1

      VIS1 = EV(I,J,KM) + 1.0/RE
      VIS2 = EV(I,J,K ) + 1.0/RE

      W2=0.5*(Q3(I,J,KP)+Q3(I,J,K ))
      W1=0.5*(Q3(I,J,K )+Q3(I,J,KM))

      APK(I,K)=(
     >      -DX3Q*    VIS2 
     >      +DX3*W2/2.0
     >       )*DT
      ACK(I,K)=1.+(
     >      +DX3Q*(VIS2+VIS1)
     >      +DX3*(W2/2.0-W1/2.0)
     >      )*DT
      AMK(I,K)=(
     >      -DX3Q*    VIS1 
     >      -DX3*W1/2.0
     >       )*DT

      R3(I,K)=QH3(I,J,K)
  30  CONTINUE
c      CALL CTDMA3I(AMK,ACK,APK,R3,QH3,J,N3M,1,N1M) 
      CALL TRIDP3 (AMK,ACK,APK,R3,1,N3M,1,N1M)
        DO K=1,N3M
        DO I=1,N1M
         QH3(I,J,K)=R3(I,K)
        ENDDO
        ENDDO
  3   CONTINUE



      DO 1 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 10 J=1,N2M
      DO 10 I=1,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMV(I)
      IUP=IPA(I)-I


      VIS1 = EVJ(I ,J,K) + 1.0/RE
      VIS2 = EVJ(IP,J,K) + 1.0/RE

      U1=0.5*(Q1(I ,J,K)+Q1(I ,J,KM))
      U2=0.5*(Q1(IP,J,K)+Q1(IP,J,KM))

      API(I,J)=IUP*(
     >         -0.5/DX(I)*    VIS2/G(IP)  
     >         +0.5/DX(I)*U2/G(IP)*0.5*DX(I)
     >       )*DT
      ACI(I,J)=(
     >      +0.5/DX(I)*(VIS2/G(IP)+VIS1/G(I))
     >      +0.5/DX(I)*(U2/G(IP)*0.5*DX(IP)-U1/G(I)*0.5*DX(IM)))*DT
     >      +1.0
      AMI(I,J)=IUM*(
     >         -0.5/DX(I)*    VIS1/G(I)  
     >         -0.5/DX(I)*U1/G(I)*0.5*DX(I)
     >       )*DT

      R1(I,J)=QH3(I,J,K)
  10  CONTINUE
      CALL TDMAJ(AMI,ACI,API,R1,R1,1,N1M,1,N2M)
      DO 11 J=1,N2M
      DO 11 I=1,N1M
 11   QH3(I,J,K)=R1(I,J)
  1   CONTINUE


      RETURN
      END


C*************** CHECK_RES1_AF ***********************     
      SUBROUTINE CHECK_RES1_AF(Q1,Q2,Q3,QH1,QH2,QH3,RDUH1
     >                        ,EV,EVI,EVJ,EVK,RES1)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL QH1(0:M1,0:M2,0:M3)
      REAL QH2(0:M1,0:M2,0:M3)
      REAL QH3(0:M1,0:M2,0:M3)
      REAL RDUH1(0:M1,0:M2,0:M3)
      REAL EV(0:M1,0:M2,0:M3)
      REAL EVI(0:M1,0:M2,0:M3)
      REAL EVJ(0:M1,0:M2,0:M3)
      REAL EVK(0:M1,0:M2,0:M3)

      REAL API(M1,M2),ACI(M1,M2),AMI(M1,M2)
      REAL APJ(M1,M2),ACJ(M1,M2),AMJ(M1,M2)
      REAL APK(M1,M3),ACK(M1,M3),AMK(M1,M3)
      REAL TEMP1(0:M1,0:M2,0:M3)
      REAL TEMP2(0:M1,0:M2,0:M3)
      REAL TEMP3(0:M1,0:M2,0:M3)
      
      EQUIVALENCE(API,APJ)
      EQUIVALENCE(ACI,ACJ)
      EQUIVALENCE(AMI,AMJ)
 
      RES1=0.0

      DO 3 J=1,N2M
      DO 30 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 30 I=2,N1M
      IP=I+1
      IM=I-1

      VIS1 = EVJ(I,J,K ) + 1.0/RE
      VIS2 = EVJ(I,J,KP) + 1.0/RE

      W1=1.0/G(I)*(0.5*DX(IM)*Q3(I,J,K )+0.5*DX(I)*Q3(IM,J,K ))
      W2=1.0/G(I)*(0.5*DX(IM)*Q3(I,J,KP)+0.5*DX(I)*Q3(IM,J,KP))

      APK(I,K)=(   -0.5*DX3Q*    VIS2  
     >       +0.5*DX3*W2*0.5)*DT
      ACK(I,K)=(   +0.5*DX3Q*(VIS2+VIS1)
     >       +0.5*DX3*(W2*0.5-W1*0.5))*DT
     >       +1.0
      AMK(I,K)=(   -0.5*DX3Q*    VIS1  
     >       -0.5*DX3*W1*0.5)*DT

  30  CONTINUE
       
      DO 31 K=1,N3M
      DO 31 I=2,N1M
      KP=KPA(K)
      KM=KMA(K)
      TEMP3(I,J,K)=AMK(I,K)*QH1(I,J,KM)
     >            +ACK(I,K)*QH1(I,J,K )
     >            +APK(I,K)*QH1(I,J,KP)
 31   CONTINUE
  3   CONTINUE


      DO 1 K=1,N3M
      DO 10 J=1,N2M
      DO 10 I=2,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMU(I)
      IUP=IPA(I)-I

      VIS1 = EV(IM,J,K) + 1.0/RE
      VIS2 = EV(I ,J,K) + 1.0/RE

      U2=0.5*(Q1(IP,J,K)+Q1(I, J,K))
      U1=0.5*(Q1(I, J,K)+Q1(IM,J,K)) 

      API(I,J)= (-1.0/G(I)*       VIS2/DX(I )      
     >      +1.0/G(I)*U2*0.5)*IUP
     >       *DT
      ACI(I,J)=(  1.0/G(I)*(VIS2/DX(I)+VIS1/DX(IM))
     >      +1.0/G(I)*(U2*0.5-U1*0.5))*DT
     >      +1.0
      AMI(I,J)= (-1.0/G(I)*       VIS1/DX(IM)  
     >      -1.0/G(I)*U1*0.5)*IUM
     >       *DT

  10  CONTINUE

      DO 11 J=1,N2M
      DO 11 I=2,N1M
      IP=I+1
      IM=I-1
      TEMP1(I,J,K)=AMI(I,J)*TEMP3(IM,J,K)
     >            +ACI(I,J)*TEMP3(I ,J,K)
     >            +API(I,J)*TEMP3(IP,J,K)
 11   CONTINUE
  1   CONTINUE

      DO 2 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 20 J=1,N2M
      JP=J+1
      JM=J-1
      JUM=J-JMU(J)
      JUP=JPA(J)-J
      DO 20 I=2,N1M
      IP=I+1
      IM=I-1

      VIS1 = EVK(I,J ,K) + 1.0/RE
      VIS2 = EVK(I,JP,K) + 1.0/RE

      V1=1.0/G(I)*(0.5*DX(IM)*Q2(I,J ,K)+0.5*DX(I)*Q2(IM,J ,K))
      V2=1.0/G(I)*(0.5*DX(IM)*Q2(I,JP,K)+0.5*DX(I)*Q2(IM,JP,K))

      APJ(I,J)=JUP*(
     >      -0.5/DY(J)*VIS2/H(JP)
     >      +0.5/DY(J)*V2/H(JP)*DY(J)/2.0
     >      )*DT

      ACJ(I,J)=( 0.5/DY(J)*(VIS2/H(JP)+VIS1/H(J))
     >    +0.5/DY(J)*(JUP*V2/H(JP)*DY(JP)/2.0
     >               -JUM*V1/H(J )*DY(JM)/2.0))*DT
     >    +1.0

      AMJ(I,J)=JUM*(
     >      -0.5/DY(J)*VIS1/H(J )
     >      -0.5/DY(J)*V1/H(J)*DY(J)/2.0
     >      )*DT

  20  CONTINUE
      DO 21 J=1,N2M
      JP=J+1
      JM=J-1
      DO 21 I=2,N1M
      TEMP2(I,J,K)=AMJ(I,J)*TEMP1(I,JM,K)
     >            +ACJ(I,J)*TEMP1(I,J ,K)
     >            +APJ(I,J)*TEMP1(I,JP,K)
 21   CONTINUE
  2   CONTINUE


      DO 40 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 40 J=1,N2M
      JP=J+1
      JM=J-1
      JUM=J-JMU(J)
      JUP=JPA(J)-J
      DO 40 I=2,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMU(I)
      IUP=IPA(I)-I

C     M12VH

      VIS1 = EVK(I,J ,K) +1.0/RE
      VIS2 = EVK(I,JP,K) +1.0/RE

      VH1=1.0/G(I)*(0.5*DX(IM)*QH2(I,J ,K)+0.5*DX(I)*QH2(IM,J ,K))
      VH2=1.0/G(I)*(0.5*DX(IM)*QH2(I,JP,K)+0.5*DX(I)*QH2(IM,JP,K))

      U1=1.0/H(J )*(0.5*DY(JM)*Q1(I,J ,K)+0.5*DY(J )*Q1(I,JM,K))
      U2=1.0/H(JP)*(0.5*DY(J )*Q1(I,JP,K)+0.5*DY(JP)*Q1(I,J ,K))

      RM12VH=JUP*0.5/DY(J)*U2*VH2
     >       -JUM*0.5/DY(J)*U1*VH1
     >
     >           -0.5/DY(J)/G(I)*( JUP*VIS2*(QH2(I,JP,K)-QH2(IM,JP,K))
     >                            -JUM*VIS1*(QH2(I,J ,K)-QH2(IM,J ,K)))

C     M13WH

      VIS1 = EVJ(I,J,K ) + 1.0/RE
      VIS2 = EVJ(I,J,KP) + 1.0/RE

      WH1=1.0/G(I)*(0.5*DX(IM)*QH3(I,J,K )+0.5*DX(I)*QH3(IM,J,K ))
      WH2=1.0/G(I)*(0.5*DX(IM)*QH3(I,J,KP)+0.5*DX(I)*QH3(IM,J,KP))

      U1=0.5*(Q1(I,J,K )+Q1(I,J,KM))
      U2=0.5*(Q1(I,J,KP)+Q1(I,J,K ))

      RM13WH=0.5*DX3*(U2*WH2-U1*WH1)
     >
     >           -0.5*DX3/G(I)*(VIS2*(QH3(I,J,KP)-QH3(IM,J,KP))
     >                         -VIS1*(QH3(I,J,K )-QH3(IM,J,K )))

   
      RESID=TEMP2(I,J,K)/DT
     >     +RM12VH+RM13WH
     >     -RDUH1(I,J,K)
      RES1=AMAX1(RES1,ABS(RESID))
c      IF (RESID.EQ.RES1) THEN
c      IMAX=I
c      JMAX=J
c      KMAX=K
c      RMAX=RDUH1(I,J,K)
c      ENDIF

 40   CONTINUE

C      WRITE(*,500) IMAX,JMAX,KMAX,RMAX
500   FORMAT(3(I4,X),E12.5)

      RETURN
      END


C*************** TDMAI **********************
      SUBROUTINE TDMAI(A,B,C,R,X,NS,NF,NIS,NIF)
      INCLUDE 'PARAM.H'
      REAL A(M1,M2),B(M1,M2),C(M1,M2),R(M1,M2),X(M1,M2)
      REAL BET,GAM(M1,0:M2)

      J=NS
      DO I=NIS,NIF
      BET=B(I,J)
      GAM(I,J)=C(I,J)/BET
      X(I,J)=R(I,J)/BET
      ENDDO
      DO 10 J=NS+1,NF
      DO I=NIS,NIF
      BET=B(I,J)-A(I,J)*GAM(I,J-1)
      GAM(I,J)=C(I,J)/BET
      X(I,J)=(R(I,J)-A(I,J)*X(I,J-1))/BET
      ENDDO
   10 CONTINUE
      DO 20 J=NF-1,NS,-1
      DO I=NIS,NIF
      X(I,J)=X(I,J)-GAM(I,J)*X(I,J+1)
      ENDDO
   20 CONTINUE
 
      RETURN
      END

C*************** TDMAJ **********************
      SUBROUTINE TDMAJ(A,B,C,R,X,NS,NF,NJS,NJF)
      INCLUDE 'PARAM.H'
      REAL A(M1,*),B(M1,*),C(M1,*),R(M1,*),X(M1,*)
      REAL BET,GAM(0:M1,M2)

      DO J=NJS,NJF-1,2
         I=NS
         BET1=1./B(I,J)
         BET2=1./B(I,J+1)
         GAM(I,J)=C(I,J)*BET1
         GAM(I,J+1)=C(I,J+1)*BET2
         X(I,J)=R(I,J)*BET1
         X(I,J+1)=R(I,J+1)*BET2
      DO I=NS+1,NF
         BET1=1./(B(I,J)-A(I,J)*GAM(I-1,J))
         BET2=1./(B(I,J+1)-A(I,J+1)*GAM(I-1,J+1))
         GAM(I,J)=C(I,J)*BET1
         GAM(I,J+1)=C(I,J+1)*BET2
         X(I,J)=(R(I,J)-A(I,J)*X(I-1,J))*BET1
         X(I,J+1)=(R(I,J+1)-A(I,J+1)*X(I-1,J+1))*BET2
      ENDDO
      DO I=NF-1,NS,-1
         X(I,J)=X(I,J)-GAM(I,J)*X(I+1,J)
         X(I,J+1)=X(I,J+1)-GAM(I,J+1)*X(I+1,J+1)
      ENDDO
      ENDDO
 
      IF(MOD(NJF-NJS+1,2).EQ.1) THEN 
         J=NJF
         I=NS
         BET1=1./B(I,J)
         GAM(I,J)=C(I,J)*BET1
         X(I,J)=R(I,J)*BET1
      DO I=NS+1,NF
         BET1=1./(B(I,J)-A(I,J)*GAM(I-1,J))
         GAM(I,J)=C(I,J)*BET1
         X(I,J)=(R(I,J)-A(I,J)*X(I-1,J))*BET1
      ENDDO
      DO I=NF-1,NS,-1
         X(I,J)=X(I,J)-GAM(I,J)*X(I+1,J)
      ENDDO
      ENDIF
 
      RETURN
      END

C  ****************************** TRIDP3 **********************
C     INVERT A PERIODIC MATRIX (X-3 TRIDIAGONAL)

      SUBROUTINE TRIDP3 (A,B,C,F,J1,J2,L1,L2)

      INCLUDE 'PARAM.H'
      REAL A(M1,M3),B(M1,M3),C(M1,M3),F(M1,M3),Q(M1,M3),S(M1,M3)
     1    ,QE(M1,M3),FN(M1),P(M1)

      JA=J1+1
      JJ=J1+J2
      DO 10 K=L1,L2
      Q(K,J1)=-C(K,J1)/B(K,J1)
      S(K,J1)=-A(K,J1)/B(K,J1)
      FN(K)=F(K,J2)
      F(K,J1)=F(K,J1)/B(K,J1)
   10 CONTINUE

C     FORWARD ELIMINATION SWEEP

      DO 20 J=JA,J2
      DO 20 K=L1,L2
      P(K)=1./(B(K,J)+A(K,J)*Q(K,J-1))
      Q(K,J)=-C(K,J)*P(K)
      S(K,J)=-A(K,J)*S(K,J-1)*P(K)
      F(K,J)=(F(K,J)-A(K,J)*F(K,J-1))*P(K)
   20 CONTINUE

C     BACKWARD PASS

      DO 30 K=L1,L2
      S(K,J2)=1.
      QE(K,J2)=0.
   30 CONTINUE
      DO 40 I=JA,J2
      J=JJ-I
      DO 40 K=L1,L2
      S(K,J)=S(K,J)+Q(K,J)*S(K,J+1)
      QE(K,J)=F(K,J)+Q(K,J)*QE(K,J+1)
   40 CONTINUE
      DO 50 K=L1,L2
      F(K,J2)=(FN(K)-C(K,J2)*QE(K,J1)-A(K,J2)*QE(K,J2-1))
     &       /(C(K,J2)*S(K,J1)+A(K,J2)*S(K,J2-1)+B(K,J2))
   50 CONTINUE

C     BACKWARD ELIMINATION PASS

      DO 60 I=JA,J2
      J=JJ-I
      DO 60 K=L1,L2
      F(K,J)=F(K,J2)*S(K,J)+QE(K,J)
   60 CONTINUE

      RETURN
      END


C*************** CTDMA3I ***********************    
      SUBROUTINE CTDMA3I(A,B,C,R,X,J,N,NIS,NIF)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      REAL A(M1,M3),B(M1,M3),C(M1,M3),R(M1,M3)
      REAL X(0:M1,0:M2,0:M3)
      REAL GAM(M1,M3)
      REAL P(M1,M3),Q(M1,M3)

      IF (N3M.EQ.1) THEN
         DO I=NIS,NIF
            X(I,J,1)=R(I,1)
         ENDDO   
         RETURN
      ENDIF

      DO I=NIS,NIF
      BET=1./B(I,1)
      GAM(I,1)=C(I,1)*BET
      X(I,J,1)=R(I,1)*BET
      P(I,1)=C(I,N)
      Q(I,1)=A(I,1)*BET
      ENDDO

      DO 10 K=2,N-1
      DO I=NIS,NIF
      BET=1./(B(I,K)-A(I,K)*GAM(I,K-1))
      GAM(I,K)=C(I,K)*BET
      P(I,K)=-P(I,K-1)*GAM(I,K-1)
      Q(I,K)=-A(I,K)*BET*Q(I,K-1)
      X(I,J,K)=(R(I,K)-A(I,K)*X(I,J,K-1))*BET
      ENDDO
 10   CONTINUE

      DO I=NIS,NIF
      BET=1./(B(I,N-1)-A(I,N-1)*GAM(I,N-2))
      P(I,N-1)=A(I,N)-P(I,N-2)*GAM(I,N-2)
      Q(I,N-1)=(C(I,N-1)-A(I,N-1)*Q(I,N-2))*BET
      
      X(I,J,N)=R(I,N)
      P(I,N)=B(I,N)
      ENDDO
      DO 20 K=1,N-1
      DO I=NIS,NIF
      X(I,J,N)=X(I,J,N)-P(I,K)*X(I,J,K)
      P(I,N)=P(I,N)-P(I,K)*Q(I,K)
      ENDDO
 20   CONTINUE
      DO I=NIS,NIF
      X(I,J,N)=X(I,J,N)/P(I,N)

      GAM(I,N-1)=0.0
      ENDDO
      DO 30 K=N-1,1,-1
      DO I=NIS,NIF
      X(I,J,K)=X(I,J,K)-GAM(I,K)*X(I,J,K+1)-Q(I,K)*X(I,J,N)
      ENDDO
 30   CONTINUE
 
      RETURN
      END


C***************** DIVCHECK ***********************    
      SUBROUTINE DIVCHECK(Q1,Q2,Q3,DIVMAX) 
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      
      DIVMAX=0.0

      DO 20 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 20 J=1,N2M
      DO 20 I=1,N1M
      DIV=ABS((Q1(I+1,J,K)-Q1(I,J,K))/DX(I)
     >       +(Q2(I,J+1,K)-Q2(I,J,K))/DY(J)
     >       +(Q3(I,J,KP)-Q3(I,J,K))*DX3)
      IF (DIV.GT.DIVMAX) THEN
      IMAX=I
      JMAX=J
      KMAX=K
      ENDIF
      DIVMAX=AMAX1(DIV,DIVMAX)
  20  CONTINUE
C      WRITE(*,*) IMAX,JMAX,KMAX,DIVMAX 
      RETURN
      END

C*********************** CFL ***********************
C     THIS SUBROUTINE CALCULATE THE MAXIMUM LOCAL CFL NUMBER
C     DEVIDED BY DT
C     AT THE CELL CENTER
C
      SUBROUTINE CFL(Q1,Q2,Q3,CFLM)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)

      CFLM=0.0
      DO 10 K=1,N3M
      KP=KPA(K)
      DO 10 J=1,N2M
      JP=J+1
      DO 10 I=1,N1M
      IP=I+1
      CFLL= ABS(Q1(I,J,K)+Q1(IP,J,K))*0.5/DX(I)
     >     +ABS(Q2(I,J,K)+Q2(I,JP,K))*0.5/DY(J)
     >     +ABS(Q3(I,J,K)+Q3(I,J,KP))*0.5*DX3

      CFLM=AMAX1(CFLM,CFLL)

      IF (CFLL.EQ.CFLM) THEN
         IMAX=I
         JMAX=J
         KMAX=K
      ENDIF

 10   CONTINUE

C      WRITE(*,*) IMAX,JMAX,KMAX
C     >          ,Q1(IMAX,JMAX,KMAX)
C     >          ,Q2(IMAX,JMAX,KMAX)
C     >          ,Q3(IMAX,JMAX,KMAX)

      RETURN
      END


C*************** INSFIELD ***********************
      SUBROUTINE INSFIELD(Q1,Q2,Q3,P)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/SIZE/ALX,ALY,ALZ

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      OPEN (12,FILE='INS_U.PLT',STATUS='UNKNOWN')
      WRITE(12,*) 'ZONE I=',N1,',J=',N2M,',K=',N3M,', F=POINT'
      DO 30 K=1,N3M
      X3=ALZ*REAL(K-0.5)/REAL(N3M)
      DO 30 J=1,N2M
      X2=0.5*(Y(J)+Y(J+1))
      DO 30 I=1,N1
      X1=X(I)
30    WRITE(12,300) X1,X2,X3,Q1(I,J,K)
      CLOSE(12)

      OPEN (13,FILE='INS_V.PLT',STATUS='UNKNOWN')
      WRITE(13,*) 'ZONE I=',N1M,',J=',N2,',K=',N3M,', F=POINT'
      DO 40 K=1,N3M
      X3=ALZ*REAL(K-0.5)/REAL(N3M)
      DO 40 J=1,N2
      X2=Y(J)
      DO 40 I=1,N1M
      X1=0.5*(X(I)+X(I+1))
40    WRITE(13,300) X1,X2,X3,Q2(I,J,K)
      CLOSE(13)

      OPEN (14,FILE='INS_W.PLT',STATUS='UNKNOWN')
      WRITE(14,*) 'ZONE I=',N1M,',J=',N2M,',K=',N3M,', F=POINT'
      DO 50 K=1,N3M
      X3=ALZ*REAL(K-1.0)/REAL(N3M)
      DO 50 J=1,N2M
      X2=0.5*(Y(J)+Y(J+1))
      DO 50 I=1,N1M
      X1=0.5*(X(I)+X(I+1))
50    WRITE(14,300) X1,X2,X3,Q3(I,J,K)
      CLOSE(14)

      OPEN (15,FILE='INS_P.PLT',STATUS='UNKNOWN')
      WRITE(15,*) 'ZONE I=',N1M,',J=',N2M,',K=',N3M,', F=POINT'
      DO 60 K=1,N3M
      X3=ALZ*REAL(K-0.5)/REAL(N3M)
      DO 60 J=1,N2M
      X2=0.5*(Y(J)+Y(J+1))
      DO 60 I=1,N1M
      X1=0.5*(X(I)+X(I+1))
60    WRITE(15,300) X1,X2,X3,P(I,J,K)
      CLOSE(15)

300   FORMAT(4(E12.5,2X))


      RETURN
      END

C*************** INS_2D ***********************
      SUBROUTINE INS_2D(Q1,Q2,Q3,P)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/SIZE/ALX,ALY,ALZ
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL DUDX11(M1,M3),DUDX12(M1,M3),DUDX13(M1,M3)  ! VEL-GRADIENT AT CELL CENTER
      REAL DUDX21(M1,M3),DUDX22(M1,M3),DUDX23(M1,M3)
      REAL DUDX31(M1,M3),DUDX32(M1,M3),DUDX33(M1,M3)

      CHARACTER*11 FILENAME

      FILENAME='INS2D.'
      NN=INDEX(FILENAME,'.')
      WRITE(UNIT=FILENAME(NN+1:),FMT='(BN,I5.5)') NTIME

      OPEN (11,FILE=FILENAME,STATUS='UNKNOWN')
      WRITE(11,*) 'VARIABLES=X,Y,U,V,W,P,OMEGA3'
      WRITE(11,*) 'ZONE I=',N1M,',J=',N2M,', F=POINT'

      K=1
      DO 40 J=1,N2M
      CALL DUDX(Q1,Q2,Q3
     >         ,DUDX11,DUDX12,DUDX13
     >         ,DUDX21,DUDX22,DUDX23
     >         ,DUDX31,DUDX32,DUDX33,J)
      DO 40 I=1,N1M

C     INSTANTANEOUS VORTICITY AT CELL CENTER

      VOR1=DUDX32(I,K)-DUDX23(I,K)
      VOR2=DUDX13(I,K)-DUDX31(I,K)
      VOR3=DUDX21(I,K)-DUDX12(I,K)

      WRITE(11,101) 0.5*(X(I)+X(I+1)),0.5*(Y(J)+Y(J+1)),
     >              0.5*(Q1(I,J,K)+Q1(I+1,J,K)),
     >              0.5*(Q2(I,J,K)+Q2(I,J+1,K)),
     >              0.5*(Q3(I,J,K)+Q3(I,J,K+1)),
     >              P(I,J,K),VOR3
40    CONTINUE
      CLOSE(11)
101   FORMAT(7(E12.5,2X))



      RETURN
      END


C*************** THIST *********************
      SUBROUTINE THIST(TIME,Q1,Q2,Q3,P)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH2/Y(0:M2)

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      OPEN (70,FILE='THIST_U.PLT',STATUS='UNKNOWN',POSITION='APPEND')
      OPEN (71,FILE='THIST_V.PLT',STATUS='UNKNOWN',POSITION='APPEND')
      OPEN (72,FILE='THIST_W.PLT',STATUS='UNKNOWN',POSITION='APPEND')
      OPEN (73,FILE='THIST_P.PLT',STATUS='UNKNOWN',POSITION='APPEND')
    
      IS=98
 
      MONI_X_1=2          ! to check right inflow reading
      MONI_X_2=IS+3       ! to check right boundary condition
      MONI_Y_1=1
      MONI_Y_2=N2M/3
      MONI_Y_3=N2M
      MONI_Y_4=N2
      
      WRITE(70,100) TIME,Q1(MONI_X_1,MONI_Y_1,N3M)
     >                  ,Q1(MONI_X_1,MONI_Y_2,N3M)
     >                  ,Q1(MONI_X_1,MONI_Y_3,N3M)
     >                  ,Q1(MONI_X_1,MONI_Y_4,N3M)
     >                  ,Q1(MONI_X_2,MONI_Y_1,N3M)
     >                  ,Q1(MONI_X_2,MONI_Y_2,N3M)
     >                  ,Q1(MONI_X_2,MONI_Y_3,N3M)
     >                  ,Q1(MONI_X_2,MONI_Y_4,N3M)

      WRITE(71,100) TIME,Q2(MONI_X_1,MONI_Y_1+1,N3M)
     >                  ,Q2(MONI_X_1,MONI_Y_2,N3M)
     >                  ,Q2(MONI_X_1,MONI_Y_3,N3M)
     >                  ,Q2(MONI_X_1,MONI_Y_4,N3M)
     >                  ,Q2(MONI_X_2,MONI_Y_1+1,N3M)
     >                  ,Q2(MONI_X_2,MONI_Y_2,N3M)
     >                  ,Q2(MONI_X_2,MONI_Y_3,N3M)
     >                  ,Q2(MONI_X_2,MONI_Y_4,N3M)

      WRITE(72,100) TIME,Q3(MONI_X_1,MONI_Y_1,N3M)
     >                  ,Q3(MONI_X_1,MONI_Y_2,N3M)
     >                  ,Q3(MONI_X_1,MONI_Y_3,N3M)
     >                  ,Q3(MONI_X_1,MONI_Y_4,N3M)
     >                  ,Q3(MONI_X_2,MONI_Y_1,N3M)
     >                  ,Q3(MONI_X_2,MONI_Y_2,N3M)
     >                  ,Q3(MONI_X_2,MONI_Y_3,N3M)
     >                  ,Q3(MONI_X_2,MONI_Y_4,N3M)


      WRITE(73,100) TIME,P(MONI_X_1,MONI_Y_1,N3M)
     >                  ,P(MONI_X_1,MONI_Y_2,N3M)
     >                  ,P(MONI_X_1,MONI_Y_3,N3M)
     >                  ,P(MONI_X_1,MONI_Y_4,N3M)
     >                  ,P(MONI_X_2,MONI_Y_1,N3M)
     >                  ,P(MONI_X_2,MONI_Y_2,N3M)
     >                  ,P(MONI_X_2,MONI_Y_3,N3M)
     >                  ,P(MONI_X_2,MONI_Y_4,N3M)

100   FORMAT (9(E12.5,2X))
  
      CLOSE(70)
      CLOSE(71)
      CLOSE(72)
      CLOSE(73)

      RETURN
      END

C*************** time_AVG **********************
      SUBROUTINE time_avg(Q1,Q2,Q3,P)

      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/SIZE/ALX,ALY,ALZ
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/FINOUT/IDTOPT,NWRITE,NREAD,IAVG,IAVG_TRIG,NPRN,INSF,NINS
      COMMON/CALC/ICONT           ! ICONT=2 continuous phase_avg with existing files 

      common/phaseavg/nphavg(0:mt),nphavg_les(0:mt),ntime_ctd

      common/t_avg1/VM_t(3,M1,M2),PM_t(M1,M2)
      common/t_avg2/VVM_t(6,M1,M2),PPM_t(M1,M2)
      common/t_avg3/V3M_t(3,M1,M2),V4M_t(3,M1,M2)
      common/t_avg4/DISSP_t(6,M1,M2)
      common/t_avg5/TRANS_t(6,2,M1,M2)
      common/t_avg6/PHI_t(6,M1,M2)
      common/t_avg7/VORM_t(3,M1,M2),VORQM_t(3,M1,M2)    
      common/t_avg8/REDIS_t(6,M1,M2)

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL DUDX11(M1,M3),DUDX12(M1,M3),DUDX13(M1,M3)  ! VEL-GRADIENT AT CELL CENTER
      REAL DUDX21(M1,M3),DUDX22(M1,M3),DUDX23(M1,M3)
      REAL DUDX31(M1,M3),DUDX32(M1,M3),DUDX33(M1,M3)

      REAL NZ

      REAL EPS1(6),EPS2(6),EPS3(6)
      REAL T1(6),T2(6)
      REAL VZ(3),VVZ(6)
      REAL VP(6)
      REAL GP(3)
      REAL V(3),VIP(3),VIM(3)
      REAL VJP(3),VJM(3),VKP(3),VKM(3)
      REAL VOR(3),VORZ(3),VORQZ(3)
      REAL PS(6)
      REAL V3Z(3),V4Z(3)


      iphase = mod(ntime_ctd,mt+1)
    
      if (nphavg(iphase).eq.0) then       ! to avoid NaN for the initial value 

        VM_t(:,:,:)     = 0.0
        VORM_t(:,:,:)   = 0.0
        VORQM_t(:,:,:)  = 0.0
        PM_t(:,:)       = 0.0
        PPM_t(:,:)      = 0.0

        VVM_t(:,:,:)    = 0.0
        V3M_t(:,:,:)    = 0.0
        V4M_t(:,:,:)    = 0.0
        DISSP_t(:,:,:)  = 0.0
        TRANS_t(:,:,:,:)= 0.0
        PHI_t(:,:,:)    = 0.0
        REDIS_t(:,:,:)  = 0.0

      endif

      if(ICONT.eq.2) then
        call read_phaseavg(vm_t,pm_t,vvm_t,ppm_t,dissp_t
     &                    ,v3m_t,v4m_t
     &                    ,trans_t,phi_t,vorm_t,vorqm_t,redis_t)         
        ICONT = -1    ! to read the prevous file only at firt step
      endif

      nphavg(iphase)=nphavg(iphase)+1
      write(*,*) iphase,nphavg(iphase)

      RNAVG=REAL(nphavg(iphase))
      NZ=1.0/DBLE(N3M)

!$omp parallel do 
!$omp$  default(shared)
!$omp&  private(JP,JM,JUM,JUP,IP,IM,IUM,IUP,
!$omp$                DUDX11,DUDX12,DUDX13,
!$omp$                DUDX21,DUDX22,DUDX23,
!$omp$                DUDX31,DUDX32,DUDX33,
!$omp$                VZ,V3Z,V4Z,VORZ,VORQZ,
!$omp$                VVZ,EPS1,EPS2,EPS3,
!$omp$                V,VOR,T1,T2,VP,PS,PZ,PPZ,
!$omp$                KP,KPP,KM,NV1,NV2,P_1,P_2)

      DO 1 J=1,N2M
      JP=J+1
      JM=J-1
      JUM=J-JMU(J)
      JUP=JPA(J)-J

      CALL DUDX(Q1,Q2,Q3
     >         ,DUDX11,DUDX12,DUDX13
     >         ,DUDX21,DUDX22,DUDX23
     >         ,DUDX31,DUDX32,DUDX33,J)

      DO 1 I=1,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMV(I)
      IUP=IPA(I)-I

      DO L=1,3
      VZ(L)=0.0
      V3Z(L)=0.0
      V4Z(L)=0.0
      VORZ(L)=0.0
      VORQZ(L)=0.0
      ENDDO
      DO L=1,6
      VVZ(L)=0.0
      EPS1(L)=0.0 
      EPS2(L)=0.0 
      EPS3(L)=0.0 
      T1(L)=0.0
      T2(L)=0.0
      VP(L)=0.0
      PS(L)=0.0
      ENDDO
      PZ=0.0
      PPZ=0.0

      DO 10 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      KPP=KPA(KP) 

C     INSTANTANEOUS VELOCITY AND PRESSURE AT CELL CENTER AND NEIGHBOR

      V(1)=(Q1(I,J,K)+Q1(IP,J,K))*0.5
      V(2)=(Q2(I,J,K)+Q2(I,JP,K))*0.5
      V(3)=(Q3(I,J,K)+Q3(I,J,KP))*0.5


C     INSTANTANEOUS VORTICITY AT CELL CENTER

      VOR(1)=DUDX32(I,K)-DUDX23(I,K)
      VOR(2)=DUDX13(I,K)-DUDX31(I,K)
      VOR(3)=DUDX21(I,K)-DUDX12(I,K)


C     AVERAGING U AND P IN SPANWISE DIRECTION

      DO L=1,3
      VZ(L)=VZ(L)+V(L)
      V3Z(L)=V3Z(L)+V(L)**3
      V4Z(L)=V4Z(L)+V(L)**4
      VORZ(L)=VORZ(L)+VOR(L)
      VORQZ(L)=VORQZ(L)+VOR(L)**2
      ENDDO
      PZ=PZ+P(I,J,K)
      PPZ=PPZ+P(I,J,K)**2

C     AVERAGING HIGHER ORDER STATISTICS E.G. UU,UUU, ETC.

      DO L=1,6
      IF(L.EQ.1) THEN
      NV1=1
      NV2=1
      ELSEIF (L.EQ.2) THEN
      NV1=2
      NV2=2
      ELSEIF (L.EQ.3) THEN
      NV1=3
      NV2=3
      ELSEIF (L.EQ.4) THEN
      NV1=1
      NV2=2
      ELSEIF (L.EQ.5) THEN
      NV1=1
      NV2=3
      ELSEIF (L.EQ.6) THEN
      NV1=2
      NV2=3
      ENDIF

C     AVERAGING U_I*U_J IN SPANWISE DIRECTION

      VVZ(L)=VVZ(L)+V(NV1)*V(NV2)  

C     AVERAGING (DU_I/DX_K)(DU_J/DX_K) IN SPANWISE DIRECTION

      EPS1(L)=EPS1(L)
     >       + 
     >       ((NV1-2.)*(NV1-3.)/(1.-2.)/(1.-3.)*DUDX11(I,K)
     >       +(NV1-1.)*(NV1-3.)/(2.-1.)/(2.-3.)*DUDX21(I,K)
     >       +(NV1-2.)*(NV1-1.)/(3.-2.)/(3.-1.)*DUDX31(I,K))
     >      *  
     >       ((NV2-2.)*(NV2-3.)/(1.-2.)/(1.-3.)*DUDX11(I,K)
     >       +(NV2-1.)*(NV2-3.)/(2.-1.)/(2.-3.)*DUDX21(I,K)
     >       +(NV2-2.)*(NV2-1.)/(3.-2.)/(3.-1.)*DUDX31(I,K))

      EPS2(L)=EPS2(L)
     >       + 
     >       ((NV1-2.)*(NV1-3.)/(1.-2.)/(1.-3.)*DUDX12(I,K)
     >       +(NV1-1.)*(NV1-3.)/(2.-1.)/(2.-3.)*DUDX22(I,K)
     >       +(NV1-2.)*(NV1-1.)/(3.-2.)/(3.-1.)*DUDX32(I,K))
     >      *  
     >       ((NV2-2.)*(NV2-3.)/(1.-2.)/(1.-3.)*DUDX12(I,K)
     >       +(NV2-1.)*(NV2-3.)/(2.-1.)/(2.-3.)*DUDX22(I,K)
     >       +(NV2-2.)*(NV2-1.)/(3.-2.)/(3.-1.)*DUDX32(I,K))

      EPS3(L)=EPS3(L)
     >       + 
     >       ((NV1-2.)*(NV1-3.)/(1.-2.)/(1.-3.)*DUDX13(I,K)
     >       +(NV1-1.)*(NV1-3.)/(2.-1.)/(2.-3.)*DUDX23(I,K)
     >       +(NV1-2.)*(NV1-1.)/(3.-2.)/(3.-1.)*DUDX33(I,K))
     >      *  
     >       ((NV2-2.)*(NV2-3.)/(1.-2.)/(1.-3.)*DUDX13(I,K)
     >       +(NV2-1.)*(NV2-3.)/(2.-1.)/(2.-3.)*DUDX23(I,K)
     >       +(NV2-2.)*(NV2-1.)/(3.-2.)/(3.-1.)*DUDX33(I,K))


C     AVERAGING U_I*U_J*U_K IN SPANWISE DIRECTION

      T1(L)=T1(L)+V(NV1)*V(NV2)*V(1)
      T2(L)=T2(L)+V(NV1)*V(NV2)*V(2)

C     AVERAGING U_I*DP/DX_J+U_J*DP/DX_I IN SPANWISE DIRECTION 


      P_2=1.0/G(IP)*(0.5*DX(I )*P(IP,J,K)+0.5*DX(IP)*P(I ,J,K))
      P_1=1.0/G(I )*(0.5*DX(IM)*P(I ,J,K)+0.5*DX(I )*P(IM,J,K))
      
      GP(1)=IUP*IUM*(P_2-P_1)/DX(I)
     >     +(1-IUM)*(P_2-P(1,J,K))/(0.5*DX(1))
     >     +(1-IUP)*(P(N1M,J,K)-P_1)/(0.5*DX(N1M))

      P_2=1.0/H(JP)*(0.5*DY(J )*P(I,JP,K)+0.5*DY(JP)*P(I,J ,K))
      P_1=1.0/H(J )*(0.5*DY(JM)*P(I,J ,K)+0.5*DY(J )*P(I,JM,K))
      
      GP(2)=JUP*JUM*(P_2-P_1)/DY(J)
     >     +(1-JUM)*(P_2-P(I,1,K))/(0.5*DY(1))
     >     +(1-JUP)*(P(I,N2M,K)-P_1)/(0.5*DY(N2M))

      GP(3)=(P(I,J,KP)-P(I,J,KM))*DX3/2

      VP(L)=VP(L)+(V(NV1)*GP(NV2)+V(NV2)*GP(NV1))

C     AVERAGING PRESSURE STRAIN CORRELATION TENSOR
C     P*(DU_I/DX_J+DU_J/DX_I)

      PS(L)=PS(L)+P(I,J,K)
     >      *(
     >          (L -2.)*(L -3.)*(L -4.)*(L -5.)*(L -6.)
     >         /(1.-2.)/(1.-3.)/(1.-4.)/(1.-5.)/(1.-6.)
     >         *(DUDX11(I,K)+DUDX11(I,K))
     >        + (L -1.)*(L -3.)*(L -4.)*(L -5.)*(L -6.)
     >         /(2.-1.)/(2.-3.)/(2.-4.)/(2.-5.)/(2.-6.)
     >         *(DUDX22(I,K)+DUDX22(I,K))
     >        + (L -1.)*(L -2.)*(L -4.)*(L -5.)*(L -6.)
     >         /(3.-1.)/(3.-2.)/(3.-4.)/(3.-5.)/(3.-6.)
     >         *(DUDX33(I,K)+DUDX33(I,K))
     >        + (L -1.)*(L -2.)*(L -3.)*(L -5.)*(L -6.)
     >         /(4.-1.)/(4.-2.)/(4.-3.)/(4.-5.)/(4.-6.)
     >         *(DUDX12(I,K)+DUDX21(I,K))
     >        + (L -1.)*(L -2.)*(L -3.)*(L -4.)*(L -6.)
     >         /(5.-1.)/(5.-2.)/(5.-3.)/(5.-4.)/(5.-6.)
     >         *(DUDX13(I,K)+DUDX31(I,K))
     >        + (L -1.)*(L -2.)*(L -3.)*(L -4.)*(L -5.)
     >         /(6.-1.)/(6.-2.)/(6.-3.)/(6.-4.)/(6.-5.)
     >         *(DUDX23(I,K)+DUDX32(I,K))
     >        )

      ENDDO

   10 CONTINUE

C     CALCULATE TIME MEAN FROM SPANWISE MEAN AT THIS TIME STEP

      DO L=1,3
      VM_t(L,I,J)   =(VZ(L)   *NZ+(RNAVG-1.)*VM_t(L,I,J)   )/RNAVG
      V3M_t(L,I,J)  =(V3Z(L)  *NZ+(RNAVG-1.)*V3M_t(L,I,J)  )/RNAVG
      V4M_t(L,I,J)  =(V4Z(L)  *NZ+(RNAVG-1.)*V4M_t(L,I,J)  )/RNAVG
      VORM_t(L,I,J) =(VORZ(L) *NZ+(RNAVG-1.)*VORM_t(L,I,J) )/RNAVG
      VORQM_t(L,I,J)=(VORQZ(L)*NZ+(RNAVG-1.)*VORQM_t(L,I,J))/RNAVG
      ENDDO
      PM_t(I,J) =(PZ *NZ+(RNAVG-1.)*PM_t(I,J) )/RNAVG
      PPM_t(I,J)=(PPZ*NZ+(RNAVG-1.)*PPM_t(I,J))/RNAVG

      DO L=1,6
      VVM_t(L,I,J)  =(VVZ(L)*NZ+(RNAVG-1.)*VVM_t(L,I,J))/RNAVG
      DISSP_t(L,I,J)=((EPS1(L)+EPS2(L)+EPS3(L))*NZ
     >             +(RNAVG-1.)*DISSP_t(L,I,J))/RNAVG
      TRANS_t(L,1,I,J)=(T1(L)*NZ+(RNAVG-1.)*TRANS_t(L,1,I,J))/RNAVG
      TRANS_t(L,2,I,J)=(T2(L)*NZ+(RNAVG-1.)*TRANS_t(L,2,I,J))/RNAVG
      PHI_t(L,I,J)    =(VP(L)*NZ+(RNAVG-1.)*PHI_t(L,I,J)    )/RNAVG
      REDIS_t(L,I,J)  =(PS(L)*NZ+(RNAVG-1.)*REDIS_t(L,I,J)  )/RNAVG
      ENDDO

    1 CONTINUE


      IF(MOD(NTIME,NPRN).EQ.0.AND.NWRITE.EQ.1) THEN
      call write_phaseavg(vm_t,pm_t,vvm_t,ppm_t,dissp_t
     &                   ,v3m_t,v4m_t
     &                   ,trans_t,phi_t,vorm_t,vorqm_t,redis_t)         
      endif

      RETURN
      END

C*************** PHASE_AVG **********************
      SUBROUTINE phase_avg(Q1,Q2,Q3,P)

      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/SIZE/ALX,ALY,ALZ
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/CALC/ICONT           ! ICONT=2 continuous phase_avg with existing files 

      common/phaseavg/nphavg(0:mt),nphavg_les(0:mt),ntime_ctd

      REAL VM(3,M1,M2),PM(M1,M2)
      REAL VVM(6,M1,M2),PPM(M1,M2)
      REAL V3M(3,M1,M2),V4M(3,M1,M2)
      REAL DISSP(6,M1,M2)
      REAL TRANS(6,2,M1,M2)
      REAL PHI(6,M1,M2)
      REAL VORM(3,M1,M2),VORQM(3,M1,M2)    
      REAL REDIS(6,M1,M2)

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL DUDX11(M1,M3),DUDX12(M1,M3),DUDX13(M1,M3)  ! VEL-GRADIENT AT CELL CENTER
      REAL DUDX21(M1,M3),DUDX22(M1,M3),DUDX23(M1,M3)
      REAL DUDX31(M1,M3),DUDX32(M1,M3),DUDX33(M1,M3)

      REAL NZ

      REAL EPS1(6),EPS2(6),EPS3(6)
      REAL T1(6),T2(6)
      REAL VZ(3),VVZ(6)
      REAL VP(6)
      REAL GP(3)
      REAL V(3),VIP(3),VIM(3)
      REAL VJP(3),VJM(3),VKP(3),VKM(3)
      REAL VOR(3),VORZ(3),VORQZ(3)
      REAL PS(6)
      REAL V3Z(3),V4Z(3)


      iphase = mod(ntime_ctd,mt+1)
    
      if (nphavg(iphase).eq.0) then       ! to avoid NaN for the initial value 

        VM(:,:,:)     = 0.0
        VORM(:,:,:)   = 0.0
        VORQM(:,:,:)  = 0.0
        PM(:,:)       = 0.0
        PPM(:,:)      = 0.0

        VVM(:,:,:)    = 0.0
        V3M(:,:,:)    = 0.0
        V4M(:,:,:)    = 0.0
        DISSP(:,:,:)  = 0.0
        TRANS(:,:,:,:)= 0.0
        PHI(:,:,:)    = 0.0
        REDIS(:,:,:)  = 0.0

      endif

      if(ICONT.eq.2 .or. nphavg(iphase).ne.0) then
        call read_phaseavg(vm,pm,vvm,ppm,dissp
     &                    ,v3m,v4m
     &                    ,trans,phi,vorm,vorqm,redis)         
      endif

      nphavg(iphase)=nphavg(iphase)+1 

      write(*,*) iphase,nphavg(iphase)

      RNAVG=REAL(nphavg(iphase))
      NZ=1.0/DBLE(N3M)

!$omp parallel do 
!$omp$  default(shared)
!$omp&  private(JP,JM,JUM,JUP,IP,IM,IUM,IUP,
!$omp$                DUDX11,DUDX12,DUDX13,
!$omp$                DUDX21,DUDX22,DUDX23,
!$omp$                DUDX31,DUDX32,DUDX33,
!$omp$                VZ,V3Z,V4Z,VORZ,VORQZ,
!$omp$                VVZ,EPS1,EPS2,EPS3,
!$omp$                V,VOR,T1,T2,VP,PS,PZ,PPZ,
!$omp$                KP,KPP,KM,NV1,NV2,P_1,P_2)

      DO 1 J=1,N2M
      JP=J+1
      JM=J-1
      JUM=J-JMU(J)
      JUP=JPA(J)-J

      CALL DUDX(Q1,Q2,Q3
     >         ,DUDX11,DUDX12,DUDX13
     >         ,DUDX21,DUDX22,DUDX23
     >         ,DUDX31,DUDX32,DUDX33,J)

      DO 1 I=1,N1M
      IP=I+1
      IM=I-1
      IUM=I-IMV(I)
      IUP=IPA(I)-I

      DO L=1,3
      VZ(L)=0.0
      V3Z(L)=0.0
      V4Z(L)=0.0
      VORZ(L)=0.0
      VORQZ(L)=0.0
      ENDDO
      DO L=1,6
      VVZ(L)=0.0
      EPS1(L)=0.0 
      EPS2(L)=0.0 
      EPS3(L)=0.0 
      T1(L)=0.0
      T2(L)=0.0
      VP(L)=0.0
      PS(L)=0.0
      ENDDO
      PZ=0.0
      PPZ=0.0

      DO 10 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      KPP=KPA(KP) 

C     INSTANTANEOUS VELOCITY AND PRESSURE AT CELL CENTER AND NEIGHBOR

      V(1)=(Q1(I,J,K)+Q1(IP,J,K))*0.5
      V(2)=(Q2(I,J,K)+Q2(I,JP,K))*0.5
      V(3)=(Q3(I,J,K)+Q3(I,J,KP))*0.5


C     INSTANTANEOUS VORTICITY AT CELL CENTER

      VOR(1)=DUDX32(I,K)-DUDX23(I,K)
      VOR(2)=DUDX13(I,K)-DUDX31(I,K)
      VOR(3)=DUDX21(I,K)-DUDX12(I,K)


C     AVERAGING U AND P IN SPANWISE DIRECTION

      DO L=1,3
      VZ(L)=VZ(L)+V(L)
      V3Z(L)=V3Z(L)+V(L)**3
      V4Z(L)=V4Z(L)+V(L)**4
      VORZ(L)=VORZ(L)+VOR(L)
      VORQZ(L)=VORQZ(L)+VOR(L)**2
      ENDDO
      PZ=PZ+P(I,J,K)
      PPZ=PPZ+P(I,J,K)**2

C     AVERAGING HIGHER ORDER STATISTICS E.G. UU,UUU, ETC.

      DO L=1,6
      IF(L.EQ.1) THEN
      NV1=1
      NV2=1
      ELSEIF (L.EQ.2) THEN
      NV1=2
      NV2=2
      ELSEIF (L.EQ.3) THEN
      NV1=3
      NV2=3
      ELSEIF (L.EQ.4) THEN
      NV1=1
      NV2=2
      ELSEIF (L.EQ.5) THEN
      NV1=1
      NV2=3
      ELSEIF (L.EQ.6) THEN
      NV1=2
      NV2=3
      ENDIF

C     AVERAGING U_I*U_J IN SPANWISE DIRECTION

      VVZ(L)=VVZ(L)+V(NV1)*V(NV2)  

C     AVERAGING (DU_I/DX_K)(DU_J/DX_K) IN SPANWISE DIRECTION

      EPS1(L)=EPS1(L)
     >       + 
     >       ((NV1-2.)*(NV1-3.)/(1.-2.)/(1.-3.)*DUDX11(I,K)
     >       +(NV1-1.)*(NV1-3.)/(2.-1.)/(2.-3.)*DUDX21(I,K)
     >       +(NV1-2.)*(NV1-1.)/(3.-2.)/(3.-1.)*DUDX31(I,K))
     >      *  
     >       ((NV2-2.)*(NV2-3.)/(1.-2.)/(1.-3.)*DUDX11(I,K)
     >       +(NV2-1.)*(NV2-3.)/(2.-1.)/(2.-3.)*DUDX21(I,K)
     >       +(NV2-2.)*(NV2-1.)/(3.-2.)/(3.-1.)*DUDX31(I,K))

      EPS2(L)=EPS2(L)
     >       + 
     >       ((NV1-2.)*(NV1-3.)/(1.-2.)/(1.-3.)*DUDX12(I,K)
     >       +(NV1-1.)*(NV1-3.)/(2.-1.)/(2.-3.)*DUDX22(I,K)
     >       +(NV1-2.)*(NV1-1.)/(3.-2.)/(3.-1.)*DUDX32(I,K))
     >      *  
     >       ((NV2-2.)*(NV2-3.)/(1.-2.)/(1.-3.)*DUDX12(I,K)
     >       +(NV2-1.)*(NV2-3.)/(2.-1.)/(2.-3.)*DUDX22(I,K)
     >       +(NV2-2.)*(NV2-1.)/(3.-2.)/(3.-1.)*DUDX32(I,K))

      EPS3(L)=EPS3(L)
     >       + 
     >       ((NV1-2.)*(NV1-3.)/(1.-2.)/(1.-3.)*DUDX13(I,K)
     >       +(NV1-1.)*(NV1-3.)/(2.-1.)/(2.-3.)*DUDX23(I,K)
     >       +(NV1-2.)*(NV1-1.)/(3.-2.)/(3.-1.)*DUDX33(I,K))
     >      *  
     >       ((NV2-2.)*(NV2-3.)/(1.-2.)/(1.-3.)*DUDX13(I,K)
     >       +(NV2-1.)*(NV2-3.)/(2.-1.)/(2.-3.)*DUDX23(I,K)
     >       +(NV2-2.)*(NV2-1.)/(3.-2.)/(3.-1.)*DUDX33(I,K))


C     AVERAGING U_I*U_J*U_K IN SPANWISE DIRECTION

      T1(L)=T1(L)+V(NV1)*V(NV2)*V(1)
      T2(L)=T2(L)+V(NV1)*V(NV2)*V(2)

C     AVERAGING U_I*DP/DX_J+U_J*DP/DX_I IN SPANWISE DIRECTION 


      P_2=1.0/G(IP)*(0.5*DX(I )*P(IP,J,K)+0.5*DX(IP)*P(I ,J,K))
      P_1=1.0/G(I )*(0.5*DX(IM)*P(I ,J,K)+0.5*DX(I )*P(IM,J,K))
      
      GP(1)=IUP*IUM*(P_2-P_1)/DX(I)
     >     +(1-IUM)*(P_2-P(1,J,K))/(0.5*DX(1))
     >     +(1-IUP)*(P(N1M,J,K)-P_1)/(0.5*DX(N1M))

      P_2=1.0/H(JP)*(0.5*DY(J )*P(I,JP,K)+0.5*DY(JP)*P(I,J ,K))
      P_1=1.0/H(J )*(0.5*DY(JM)*P(I,J ,K)+0.5*DY(J )*P(I,JM,K))
      
      GP(2)=JUP*JUM*(P_2-P_1)/DY(J)
     >     +(1-JUM)*(P_2-P(I,1,K))/(0.5*DY(1))
     >     +(1-JUP)*(P(I,N2M,K)-P_1)/(0.5*DY(N2M))

      GP(3)=(P(I,J,KP)-P(I,J,KM))*DX3/2

      VP(L)=VP(L)+(V(NV1)*GP(NV2)+V(NV2)*GP(NV1))

C     AVERAGING PRESSURE STRAIN CORRELATION TENSOR
C     P*(DU_I/DX_J+DU_J/DX_I)

      PS(L)=PS(L)+P(I,J,K)
     >      *(
     >          (L -2.)*(L -3.)*(L -4.)*(L -5.)*(L -6.)
     >         /(1.-2.)/(1.-3.)/(1.-4.)/(1.-5.)/(1.-6.)
     >         *(DUDX11(I,K)+DUDX11(I,K))
     >        + (L -1.)*(L -3.)*(L -4.)*(L -5.)*(L -6.)
     >         /(2.-1.)/(2.-3.)/(2.-4.)/(2.-5.)/(2.-6.)
     >         *(DUDX22(I,K)+DUDX22(I,K))
     >        + (L -1.)*(L -2.)*(L -4.)*(L -5.)*(L -6.)
     >         /(3.-1.)/(3.-2.)/(3.-4.)/(3.-5.)/(3.-6.)
     >         *(DUDX33(I,K)+DUDX33(I,K))
     >        + (L -1.)*(L -2.)*(L -3.)*(L -5.)*(L -6.)
     >         /(4.-1.)/(4.-2.)/(4.-3.)/(4.-5.)/(4.-6.)
     >         *(DUDX12(I,K)+DUDX21(I,K))
     >        + (L -1.)*(L -2.)*(L -3.)*(L -4.)*(L -6.)
     >         /(5.-1.)/(5.-2.)/(5.-3.)/(5.-4.)/(5.-6.)
     >         *(DUDX13(I,K)+DUDX31(I,K))
     >        + (L -1.)*(L -2.)*(L -3.)*(L -4.)*(L -5.)
     >         /(6.-1.)/(6.-2.)/(6.-3.)/(6.-4.)/(6.-5.)
     >         *(DUDX23(I,K)+DUDX32(I,K))
     >        )

      ENDDO

   10 CONTINUE

C     CALCULATE TIME MEAN FROM SPANWISE MEAN AT THIS TIME STEP

      DO L=1,3
      VM(L,I,J)   =(VZ(L)   *NZ+(RNAVG-1.)*VM(L,I,J)   )/RNAVG
      V3M(L,I,J)  =(V3Z(L)  *NZ+(RNAVG-1.)*V3M(L,I,J)  )/RNAVG
      V4M(L,I,J)  =(V4Z(L)  *NZ+(RNAVG-1.)*V4M(L,I,J)  )/RNAVG
      VORM(L,I,J) =(VORZ(L) *NZ+(RNAVG-1.)*VORM(L,I,J) )/RNAVG
      VORQM(L,I,J)=(VORQZ(L)*NZ+(RNAVG-1.)*VORQM(L,I,J))/RNAVG
      ENDDO
      PM(I,J) =(PZ *NZ+(RNAVG-1.)*PM(I,J) )/RNAVG
      PPM(I,J)=(PPZ*NZ+(RNAVG-1.)*PPM(I,J))/RNAVG

      DO L=1,6
      VVM(L,I,J)  =(VVZ(L)*NZ+(RNAVG-1.)*VVM(L,I,J))/RNAVG
      DISSP(L,I,J)=((EPS1(L)+EPS2(L)+EPS3(L))*NZ
     >             +(RNAVG-1.)*DISSP(L,I,J))/RNAVG
      TRANS(L,1,I,J)=(T1(L)*NZ+(RNAVG-1.)*TRANS(L,1,I,J))/RNAVG
      TRANS(L,2,I,J)=(T2(L)*NZ+(RNAVG-1.)*TRANS(L,2,I,J))/RNAVG
      PHI(L,I,J)    =(VP(L)*NZ+(RNAVG-1.)*PHI(L,I,J)    )/RNAVG
      REDIS(L,I,J)  =(PS(L)*NZ+(RNAVG-1.)*REDIS(L,I,J)  )/RNAVG
      ENDDO

    1 CONTINUE

      call write_phaseavg(vm,pm,vvm,ppm,dissp
     &                   ,v3m,v4m
     &                   ,trans,phi,vorm,vorqm,redis)         

      RETURN
      END

C*************** read_phaseavg  **********************
      SUBROUTINE read_phaseavg(vm,pm,vvm,ppm,dissp
     &                        ,v3m,v4m
     &                        ,trans,phi,vorm,vorqm,redis)         
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      common/phaseavg/nphavg(0:mt),nphavg_les(0:mt),ntime_ctd

      real VM(3,M1,M2),PM(M1,M2)
      real VVM(6,M1,M2),PPM(M1,M2)
      real V3M(3,M1,M2),V4M(3,M1,M2)
      real DISSP(6,M1,M2)
      real TRANS(6,2,M1,M2)
      real PHI(6,M1,M2)
      real VORM(3,M1,M2),VORQM(3,M1,M2)    
      real REDIS(6,M1,M2)

      character*30 phasefile 

      iphase = mod(ntime_ctd,mt+1)

      phasefile='phaseavg.'
      nn=index(phasefile,'.')
      write(unit=phasefile(nn+1:),fmt='(bn,i4.4)') iphase

      OPEN(11,FILE=phasefile,STATUS='UNKNOWN',form='unformatted')
      
          read(11) nphavg(iphase) 
      DO I=1,N1M
        DO J=1,N2M
          read(11) (VM(L,I,J),L=1,3),PM(I,J),PPM(I,J)
          read(11) (VVM(L,I,J),L=1,6)
          read(11) (V3M(L,I,J),L=1,3)
          read(11) (V4M(L,I,J),L=1,3)
          read(11) (TRANS(L,1,I,J),L=1,6),
     >             (TRANS(L,2,I,J),L=1,6)
          read(11) (DISSP(L,I,J),L=1,6)
          read(11) (PHI(L,I,J),L=1,6)
          read(11) (VORM(L,I,J),L=1,3),(VORQM(L,I,J),L=1,3)
          read(11) (REDIS(L,I,J),L=1,6)
        ENDDO
      ENDDO

      CLOSE(11)

      RETURN
      END

C*************** write_phaseavg  **********************
      SUBROUTINE write_phaseavg(vm,pm,vvm,ppm,dissp
     &                         ,v3m,v4m
     &                         ,trans,phi,vorm,vorqm,redis)         
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      common/phaseavg/nphavg(0:mt),nphavg_les(0:mt),ntime_ctd

      real VM(3,M1,M2),PM(M1,M2)
      real VVM(6,M1,M2),PPM(M1,M2)
      real V3M(3,M1,M2),V4M(3,M1,M2)
      real DISSP(6,M1,M2)
      real TRANS(6,2,M1,M2)
      real PHI(6,M1,M2)
      real VORM(3,M1,M2),VORQM(3,M1,M2)    
      real REDIS(6,M1,M2)

      character*30 phasefile 

      iphase = mod(ntime_ctd,mt+1)

      phasefile='phaseavg.'
      nn=index(phasefile,'.')
      write(unit=phasefile(nn+1:),fmt='(bn,i4.4)') iphase

      OPEN(11,FILE=phasefile,STATUS='UNKNOWN',form='unformatted')
               
          write(11) nphavg(iphase) 
      DO I=1,N1M
        DO J=1,N2M
          write(11) (VM(L,I,J),L=1,3),PM(I,J),PPM(I,J)
          write(11) (VVM(L,I,J),L=1,6)
          write(11) (V3M(L,I,J),L=1,3)
          write(11) (V4M(L,I,J),L=1,3)
          write(11) (TRANS(L,1,I,J),L=1,6),
     >              (TRANS(L,2,I,J),L=1,6)
          write(11) (DISSP(L,I,J),L=1,6)
          write(11) (PHI(L,I,J),L=1,6)
          write(11) (VORM(L,I,J),L=1,3),(VORQM(L,I,J),L=1,3)
          write(11) (REDIS(L,I,J),L=1,6)
        ENDDO
      ENDDO

      CLOSE(11)

      RETURN
      END

C*************** RESCALING ****************************
C
C    "GENERATION OF TURBULENT INFLOW DATA FOR SPATIALLY-
C     DEVELOPING BOUNDARY LAYER SIMULATION,"
C
C                                     LUND, WU & SQUIRES
C            J. COMPUT. PHYS. VOL. 140, PP.233-258, 1998
C

      SUBROUTINE RESCALING(Q1,Q2,Q3,TIME)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/BCON/UBC1(M1,M3,3),UBC2(M1,M3,3)            
     >           ,UBC3(M2,M3,3),UBC4(M2,M3,3)
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/SIZE/ALX,ALY,ALZ
      COMMON/RECYCLE1/IRC,UT_IN,UT_RC,DELTA_RC
      COMMON/RECYCLE2/THETA_IN,THETA_RC 
      COMMON/RECYCLE3/UMEAN(3,0:M1,0:M2),UM_RC(3,0:M2)

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)

      REAL UF_RC(3,0:M2,0:M3)
      REAL UM_IN(3,2,0:M2)
      REAL UF_IN(3,2,0:M2,0:M3)

C     THE I INDEX OF THE RECYCLE STATION
      IRC=N1M*4/5

c      THETA_IN=270./RE
      THETA_IN=620./RE
     
      CALL GET_UM_RC(Q1,Q2,Q3,UF_RC,TIME)

      UT_RC=(ABS(UM_RC(1,1))/H(1)/RE)**0.5
      
C     FIND THE BOUNDARY LAYER THICKNESS AT THE RECYCLE STATION
C     BY USING LINEAR INTERPOLATION OF MEAN VELOCITY PROFILE

      U_INF=0.99*UM_RC(1,N2)
      DO 10 J=2,N2M
         IF (UM_RC(1,J).GE.U_INF) THEN
             D1=(Y(J-1)+Y(J))*0.5
             D2=(Y(J)+Y(J+1))*0.5
             U1=UM_RC(1,J-1)
             U2=UM_RC(1,J)
             DELTA_RC=(D2-D1)/(U2-U1)*(U_INF-U1)+D1
             GOTO 20
         ENDIF 
10    CONTINUE
20    CONTINUE

C     CALCULATE THE MOMENTUM THICKNESS AT THE RECYCLE STATION
C     BY INTEGRATING THE MEAN VELOCITY PROFILE UP TO DELTA_RC

      THETA_RC=0.0
      DO J=1,N2M    
       IF(Y(J).LE.DELTA_RC)
     & THETA_RC=THETA_RC+UM_RC(1,J)*(1.0-UM_RC(1,J))*DY(J)
      ENDDO

      UT_IN=UT_RC*(THETA_RC/THETA_IN)**(1./2./(5.-1.))

      GAMMA=UT_IN/UT_RC
 
      CALL GET_DELTA_IN(UM_IN,DELTA_IN) 
      CALL GET_UF_IN(UF_RC,UF_IN,DELTA_IN)
      CALL GET_INFLOW(UM_IN,UF_IN)

C---  FOR CHECK THE INLET MOMENTUM THICKNESS
         TH_IN=0.0
      DO 50 J=1,N2M
         U_Z=0.0
         DO K=1,N3M
            U_Z=U_Z+UBC3(J,K,1)/REAL(N3M)
         ENDDO
         IF(Y(J).LE.DELTA_IN)
     &   TH_IN=TH_IN+U_Z*(1.0-U_Z)*DY(J)  
50    CONTINUE

      WRITE(*,100) UT_IN,TH_IN*RE,DELTA_IN,Q1(1,N2,N3M)
100   FORMAT('AT INLET   : U_T=',E12.6,X,'RE_TH=',E12.6,X,'D99=',E12.6
     >       ,X,'U_INF=',E12.6)
      WRITE(*,101) UT_RC,THETA_RC*RE,DELTA_RC,Q1(IRC,N2,N3M)
101   FORMAT('AT RECYCLE : U_T=',E12.6,X,'RE_TH=',E12.6,X,'D99=',E12.6
     >       ,X,'U_INF=',E12.6)

      OPEN(21,FILE='CF.PLT',STATUS='UNKNOWN',POSITION='APPEND')
      WRITE(21,200) TIME,2.*UT_IN**2,DELTA_IN
200   FORMAT(3(E12.5,X))      
      CLOSE(21)

      RETURN
      END

C*************** GET_UM_RC ***********************************
      SUBROUTINE GET_UM_RC(Q1,Q2,Q3,UF_RC,TIME)
      INCLUDE 'PARAM.H'
      COMMON/PARA/RE
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/SIZE/ALX,ALY,ALZ
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/RECYCLE1/IRC,UT_IN,UT_RC,DELTA_RC
      COMMON/RECYCLE2/THETA_IN,THETA_RC 
      COMMON/RECYCLE3/UMEAN(3,0:M1,0:M2),UM_RC(3,0:M2)
      COMMON/FINOUT/IDTOPT,NWRITE,NREAD,IAVG,IAVG_TRIG,NPRN,INSF,NINS
      COMMON/CALC/ICONT

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL UF_RC(3,0:M2,0:M3)

C     FOR CONTINUOUS CALCULATION, OBTAIN UM_RC AND UMEAN FROM SAVED FILES 

      IF (ICONT.EQ.1.AND.NTIME.EQ.1) THEN
          OPEN(10,FILE='RESCALE.DAT',STATUS='OLD',FORM='UNFORMATTED')
          READ(10) ((UM_RC(NV,J),J=1,N2),NV=1,3)
          READ(10) ((UMEAN(1,I,J),J=1,N2M),I=1,N1)
          CLOSE(10)
      ENDIF

      IF (TIME.LE.1000.0) THEN 
          TAVG=100.               ! NOTE THE LENGTH SCALE IS THETA_IN
      ELSE IF (TIME.GE.1000.0 ) THEN 
               TAVG=1000.               
      ENDIF
     
      IF (ICONT.EQ.1) TAVG=1000.

      DO 10 J=1,N2
         TEMP1=0.0
         TEMP2=0.0
         DO 20 K=1,N3M
            TEMP1=TEMP1+Q1(IRC,J,K)/REAL(N3M)         
20          TEMP2=TEMP2+Q2(IRC,J,K)/REAL(N3M)         

         DO 21 K=1,N3M
            UF_RC(1,J,K)=Q1(IRC,J,K)-TEMP1
21          UF_RC(2,J,K)=Q2(IRC,J,K)-TEMP2

         UM_RC(1,J)=DT/TAVG*TEMP1+(1.-DT/TAVG)*UM_RC(1,J) 
         UM_RC(2,J)=DT/TAVG*TEMP2+(1.-DT/TAVG)*UM_RC(2,J) 
         IF (NTIME.EQ.1.AND.ICONT.EQ.0) THEN
             UM_RC(1,J)=TEMP1
             UM_RC(2,J)=TEMP2
         ENDIF

10    CONTINUE

      DO 30 J=1,N2M
      DO 31 K=1,N3M
31       UF_RC(3,J,K)=Q3(IRC,J,K)
30    UM_RC(3,J)=0.0

      DO 40 I=1,N1
      DO 40 J=1,N2M
         TEMP=0.0
         DO 41 K=1,N3M
41          TEMP=TEMP+Q1(I,J,K)/REAL(N3M)         
         UMEAN(1,I,J)=DT/TAVG*TEMP+(1.-DT/TAVG)*UMEAN(1,I,J) 
         IF (NTIME.EQ.1.AND.ICONT.EQ.0) UMEAN(1,I,J)=TEMP 
40    CONTINUE



      IF(MOD(NTIME,NPRN).EQ.0.AND.NWRITE.EQ.1) THEN 
         OPEN(10,FILE='RESCALE.DAT',STATUS='UNKNOWN'
     >          ,FORM='UNFORMATTED')
         WRITE(10) ((UM_RC(NV,J),J=1,N2),NV=1,3)
         WRITE(10) ((UMEAN(1,I,J),J=1,N2M),I=1,N1)
         CLOSE(10)
      ENDIF

      RETURN
      END 

C*************** GET_UM_IN *******************************
      SUBROUTINE GET_UM_IN(UM_IN,DELTA_IN)
      INCLUDE 'PARAM.H'
      COMMON/PARA/RE
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/SIZE/ALX,ALY,ALZ
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/RECYCLE1/IRC,UT_IN,UT_RC,DELTA_RC
      COMMON/RECYCLE2/THETA_IN,THETA_RC 
      COMMON/RECYCLE3/UMEAN(3,0:M1,0:M2),UM_RC(3,0:M2)

      REAL UM_IN(3,2,0:M2)
C          UM_IN(NV,ND,J )
C                NV : VARIABLES
C                ND : NONDIMENSIONALIZED LENGTH SCALE
C                RETURN MEAN PROFILE AS UM_IN(NV,1,J)
C                THROUGH THIS SUBROUTINE
     
      GAMMA=UT_IN/UT_RC

C     CALCULATE THE MEAN VELOCITY PROFILE ; U1

      NV=1

C          JJ=1    FOR CONSISTENCY OF U_T AND UM_IN 
      UM_IN(1,1,1)=UT_IN**2*H(1)*RE    ! README

      DO 1 JJ=2,N2M    
C     LENGTH SCALE : INNER VARIABLES
      YP_IN=RE*UT_IN*0.5*(Y(JJ)+Y(JJ+1))
      YP_RC_MAX=RE*UT_RC*0.5*(Y(N2M)+Y(N2))

      IF (YP_IN.GE.YP_RC_MAX) THEN 
      UM_IN(NV,1,JJ)=1.0     
      GOTO 111
      ENDIF

      DO 10 J=1,N2M    
      YP_RC=RE*UT_RC*0.5*(Y(J)+Y(J+1))
      IF (YP_RC.GT.YP_IN) THEN 
      Y2=YP_RC
      Y1=RE*UT_RC*0.5*(Y(J)+Y(J-1))
      J_RC=J
      GOTO 11
      ENDIF
10    CONTINUE
11    CONTINUE
      U1=UM_RC(NV,J_RC-1)
      U2=UM_RC(NV,J_RC)
      UM_IN(NV,1,JJ)=((U2-U1)/(Y2-Y1)*(YP_IN-Y2)+U2)
     >              *GAMMA
111   CONTINUE

C     LENGTH SCALE : OUTER VARIABLES
      ETHA_IN=0.5*(Y(JJ)+Y(JJ+1))/DELTA_IN
      ETHA_RC_MAX=0.5*(Y(N2M)+Y(N2))/DELTA_RC

      IF (ETHA_IN.GT.ETHA_RC_MAX) THEN 
      UM_IN(NV,2,JJ)=1.0       
      GOTO 222
      ENDIF

      DO 20 J=1,N2M    
      ETHA_RC=0.5*(Y(J)+Y(J+1))/DELTA_RC
      IF (ETHA_RC.GT.ETHA_IN) THEN 
      Y2=ETHA_RC
      Y1=0.5*(Y(J)+Y(J-1))/DELTA_RC
      J_RC=J
      GOTO 21
      ENDIF
20    CONTINUE
21    CONTINUE

      U1=UM_RC(NV,J_RC-1)
      U2=UM_RC(NV,J_RC)
      UM_IN(NV,2,JJ)=((U2-U1)/(Y2-Y1)*(ETHA_IN-Y2)+U2)
     >              *GAMMA+(1.0-GAMMA)
222   CONTINUE   

C     WEIGHTED AVERAGE WITH FUNCTION, W
      W=WEIGHT(ETHA_IN)
      UM_IN(NV,1,JJ)=UM_IN(NV,1,JJ)*(1.-W)+UM_IN(NV,2,JJ)*W
      
1     CONTINUE

C     CALCULATE THE MEAN VELOCITY PROFILE ; U2

      NV=2
      DO 2 JJ=2,N2M   

C     LENGTH SCALE : INNER VARIABLES

      YP_IN=RE*UT_IN*Y(JJ)
      YP_RC_MAX=RE*UT_RC*Y(N2)

      IF (YP_IN.GE.YP_RC_MAX) THEN 
      UM_IN(NV,1,JJ)=UM_IN(NV,1,JJ-1)
      GOTO 333
      ENDIF

      DO 30 J=1,N2    
      YP_RC=RE*UT_RC*Y(J)
      IF (YP_RC.GE.YP_IN) THEN 
      Y2=YP_RC
      Y1=RE*UT_RC*Y(J-1)
      J_RC=J
      GOTO 31
      ENDIF
30    CONTINUE
31    CONTINUE
      U1=UM_RC(NV,J_RC-1)
      U2=UM_RC(NV,J_RC)
      UM_IN(NV,1,JJ)=(U2-U1)/(Y2-Y1)*(YP_IN-Y2)+U2
333   CONTINUE

C     LENGTH SCALE : OUTER VARIABLES

      ETHA_IN=Y(JJ)/DELTA_IN
      ETHA_RC_MAX=Y(N2)/DELTA_RC
      IF (ETHA_IN.GE.ETHA_RC_MAX) THEN 
      UM_IN(NV,2,JJ)=UM_IN(NV,2,JJ-1)
      GOTO 444
      ENDIF

      DO 40 J=1,N2    
      ETHA_RC=Y(J)/DELTA_RC
      IF (ETHA_RC.GE.ETHA_IN) THEN 
      Y2=ETHA_RC
      Y1=Y(J-1)/DELTA_RC
      J_RC=J
      GOTO 41
      ENDIF
40    CONTINUE
41    CONTINUE
      U1=UM_RC(NV,J_RC-1)
      U2=UM_RC(NV,J_RC)
      UM_IN(NV,2,JJ)=(U2-U1)/(Y2-Y1)*(ETHA_IN-Y2)+U2
444   CONTINUE   

C     WEIGHTED AVERAGE WITH FUNCTION, W
      W=WEIGHT(ETHA_IN)
      UM_IN(NV,1,JJ)=UM_IN(NV,1,JJ)*(1.-W)+UM_IN(NV,2,JJ)*W
     
2     CONTINUE
      
C     CALCULATE THE MEAN VELOCITY PROFILE ; U3

      NV=3
      DO 3 JJ=1,N2M   
      UM_IN(NV,1,JJ)=0.0
3     CONTINUE

      RETURN
      END


C****************** WEIGHT ***********************************
      REAL FUNCTION WEIGHT(ETHA)
      ALPHA=4.0
      B=0.2
      WEIGHT=0.5*(1.0+TANH(ALPHA*(ETHA-B)/((1.-2.*B)*ETHA+B))
     >           /TANH(ALPHA))
      IF (ETHA.GE.1.) WEIGHT=1.0 
      
      RETURN
      END

C*************** GET_DELTA_IN ********************************* 
      SUBROUTINE GET_DELTA_IN(UM_IN,DELTA_IN) 
      INCLUDE 'PARAM.H'
      COMMON/PARA/RE
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/SIZE/ALX,ALY,ALZ
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/RECYCLE1/IRC,UT_IN,UT_RC,DELTA_RC
      COMMON/RECYCLE2/THETA_IN,THETA_RC 
      COMMON/RECYCLE3/UMEAN(3,0:M1,0:M2),UM_RC(3,0:M2)

      REAL UM_IN(3,2,0:M2)

C     THIS IS THE NEWTON'S METHOD.
      ITERMAX=100
      D1=10.0  ! INITIAL GUESS FOR NEWTON'S METHOD

      ITER=0 
      WRITE(*,*) 'GET THE DELTA AT THE INLET STATION.' 
      WRITE(*,*) 'ITER    DELTA_IN     RE_THETA' 
      WRITE(*,*) '----------------------------------' 
      
1     CONTINUE
      ITER=ITER+1 
      IF (ITER.GE.ITERMAX) THEN
      OPEN(24,FILE='ALERT.TXT',STATUS='UNKNOWN',POSITION='APPEND') 
      WRITE(24,244) NTIME
244   FORMAT('CANNOT FIND DELTA_IN FOR RE_THETA=300 AT'
     > ,I6,'TH STEP') 
      CLOSE(24)
      GOTO 2
      ENDIF  

      CALL GET_UM_IN(UM_IN,D1)
      THETA=0.0
      DO J=1,N2M 
       IF(Y(J).LE.D1)
     & THETA=THETA+UM_IN(1,1,J)*(1.0-UM_IN(1,1,J))*DY(J)
      ENDDO

      F_D1=THETA-THETA_IN

      IF (ABS(F_D1).GE.1E-5) THEN        !1

      WRITE(*,100) ITER,D1,(F_D1+THETA_IN)*RE
100   FORMAT(I4,2X,2(E12.5,2X)) 

C     TO CALCULATE THE DERIVATIVE AT D1
C     THROUGH CENTRAL DIFFERENCING
      DEL_D=1E-10 
      D_P=D1+DEL_D
      D_M=D1-DEL_D

      CALL GET_UM_IN(UM_IN,D_P)
      THETA=0.0
      DO J=1,N2M      
       IF(Y(J).LE.D_P)
     & THETA=THETA+UM_IN(1,1,J)*(1.0-UM_IN(1,1,J))*DY(J)
      ENDDO
      F_D_P=THETA-THETA_IN

      CALL GET_UM_IN(UM_IN,D_M)
      THETA=0.0
      DO J=1,N2M      
       IF(Y(J).LE.D_M)
     & THETA=THETA+UM_IN(1,1,J)*(1.0-UM_IN(1,1,J))*DY(J)
      ENDDO
      F_D_M=THETA-THETA_IN

      D_F_D1=(F_D_P-F_D_M)/(D_P-D_M) ! GET DERIVATIVE AT D1 


      IF (ABS(D_F_D1).LE.1E-10) THEN
       D_F_D1=1E-10
       OPEN(24,FILE='ALERT.TXT',STATUS='UNKNOWN',POSITION='APPEND') 
       WRITE(24,245) NTIME
245    FORMAT('IN NEWTON`S METHOD, DENOMINATOR IS ZERO AT',I8,'TH')
       CLOSE(24)
      ENDIF

      D1=D1-F_D1/D_F_D1   ! NEWTON'S ITERATION

      GOTO 1
      ENDIF    !1

2     CONTINUE     

      WRITE(*,100) ITER,D1,(F_D1+THETA_IN)*RE
      DELTA_IN=D1
      WRITE(*,*)

      RETURN
      END  

C*************** GET_UF_IN ******************************
      SUBROUTINE GET_UF_IN(UF_RC,UF_IN,DELTA_IN)
      INCLUDE 'PARAM.H'
      COMMON/PARA/RE
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/SIZE/ALX,ALY,ALZ
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/RECYCLE1/IRC,UT_IN,UT_RC,DELTA_RC
      COMMON/RECYCLE2/THETA_IN,THETA_RC 
      COMMON/RECYCLE3/UMEAN(3,0:M1,0:M2),UM_RC(3,0:M2)

      REAL UF_RC(3,0:M2,0:M3)
      REAL UF_IN(3,2,0:M2,0:M3)
C          UF_IN(NV,ND,J ,K )
C                NV : VARIABLES
C                ND : NONDIMENSIONALIZED LENGTH SCALE
C                RETURN RESULTING OUTPUT AS UF_IN(NV,1,J,K)
C                THROUGH THIS SUBROUTINE
    
      GAMMA=UT_IN/UT_RC

C     CALCULATE THE FLUCTUATING COMPONENTS ; U1

      NV=1
      DO 1 K=1,N3M
      DO 1 JJ=1,N2M    

C     LENGTH SCALE : INNER VARIABLES

      YP_IN=RE*UT_IN*0.5*(Y(JJ)+Y(JJ+1))
      YP_RC_MAX=RE*UT_RC*0.5*(Y(N2M)+Y(N2))

      IF (YP_IN.GE.YP_RC_MAX) THEN 
      UF_IN(NV,1,JJ,K)=0.0       
      GOTO 111
      ENDIF

      DO 10 J=1,N2M    
      YP_RC=RE*UT_RC*0.5*(Y(J)+Y(J+1))
      IF (YP_RC.GE.YP_IN) THEN 
      Y2=YP_RC
      Y1=RE*UT_RC*0.5*(Y(J)+Y(J-1))
      J_RC=J
      GOTO 11
      ENDIF
10    CONTINUE
11    CONTINUE
      U1=UF_RC(NV,J_RC-1,K)
      U2=UF_RC(NV,J_RC,K)
      UF_IN(NV,1,JJ,K)=((U2-U1)/(Y2-Y1)*(YP_IN-Y2)+U2)
     >                *GAMMA
111   CONTINUE

C     LENGTH SCALE : OUTER VARIABLES

      ETHA_IN=0.5*(Y(JJ)+Y(JJ+1))/DELTA_IN
      ETHA_RC_MAX=0.5*(Y(N2M)+Y(N2))/DELTA_RC

      IF (ETHA_IN.GE.ETHA_RC_MAX) THEN 
      UF_IN(NV,2,JJ,K)=0.0      
      GOTO 222
      ENDIF

      DO 20 J=1,N2M    
      ETHA_RC=0.5*(Y(J)+Y(J+1))/DELTA_RC
      IF (ETHA_RC.GE.ETHA_IN) THEN 
      Y2=ETHA_RC
      Y1=0.5*(Y(J)+Y(J-1))/DELTA_RC
      J_RC=J
      GOTO 21
      ENDIF
20    CONTINUE
21    CONTINUE
      U1=UF_RC(NV,J_RC-1,K)
      U2=UF_RC(NV,J_RC,K)
      UF_IN(NV,2,JJ,K)=((U2-U1)/(Y2-Y1)*(ETHA_IN-Y2)+U2)
     >                *GAMMA
222   CONTINUE   

C     WEIGHTED AVERAGE WITH FUNCTION, W
      W=WEIGHT(ETHA_IN)
      UF_IN(NV,1,JJ,K)=UF_IN(NV,1,JJ,K)*(1.-W)+UF_IN(NV,2,JJ,K)*W
            
1     CONTINUE
 

C     CALCULATE THE FLUCTUATING COMPONENTS ; U2
      NV=2

      DO 2 K=1,N3M
      DO 2 JJ=2,N2M   

C     LENGTH SCALE : INNER VARIABLES

      YP_IN=RE*UT_IN*Y(JJ)
      YP_RC_MAX=RE*UT_RC*Y(N2)

      IF (YP_IN.GE.YP_RC_MAX) THEN 
      UF_IN(NV,1,JJ,K)=0.0        
      GOTO 333
      ENDIF

      DO 30 J=1,N2    
      YP_RC=RE*UT_RC*Y(J)
      IF (YP_RC.GE.YP_IN) THEN 
      Y2=YP_RC
      Y1=RE*UT_RC*Y(J-1)
      J_RC=J
      GOTO 31
      ENDIF
30    CONTINUE
31    CONTINUE
      U1=UF_RC(NV,J_RC-1,K)
      U2=UF_RC(NV,J_RC,K)
      UF_IN(NV,1,JJ,K)=((U2-U1)/(Y2-Y1)*(YP_IN-Y2)+U2)
     >                *GAMMA
333   CONTINUE

C     LENGTH SCALE : OUTER VARIABLES

      ETHA_IN=Y(JJ)/DELTA_IN
      ETHA_RC_MAX=Y(N2)/DELTA_RC
      IF (ETHA_IN.GE.ETHA_RC_MAX) THEN 
      UF_IN(NV,2,JJ,K)=0.0     
      GOTO 444
      ENDIF

      DO 40 J=1,N2    
      ETHA_RC=Y(J)/DELTA_RC
      IF (ETHA_RC.GE.ETHA_IN) THEN 
      Y2=ETHA_RC
      Y1=Y(J-1)/DELTA_RC
      J_RC=J
      GOTO 41
      ENDIF
40    CONTINUE
41    CONTINUE
      U1=UF_RC(NV,J_RC-1,K)
      U2=UF_RC(NV,J_RC,K)
      UF_IN(NV,2,JJ,K)=((U2-U1)/(Y2-Y1)*(ETHA_IN-Y2)+U2)
     >                *GAMMA
444   CONTINUE   

C     WEIGHTED AVERAGE WITH FUNCTION, W
      W=WEIGHT(ETHA_IN)
      UF_IN(NV,1,JJ,K)=UF_IN(NV,1,JJ,K)*(1.-W)+UF_IN(NV,2,JJ,K)*W
      
2     CONTINUE


C     CALCULATE THE FLUCTUATING COMPONENTS ; U3

      NV=3
      DO 3 K=1,N3M
      DO 3 JJ=1,N2M    

C     LENGTH SCALE : INNER VARIABLES

      YP_IN=RE*UT_IN*0.5*(Y(JJ)+Y(JJ+1))
      YP_RC_MAX=RE*UT_RC*0.5*(Y(N2M)+Y(N2))

      IF (YP_IN.GE.YP_RC_MAX) THEN 
      UF_IN(NV,1,JJ,K)=0.0       
      GOTO 555
      ENDIF

      DO 50 J=1,N2M    
      YP_RC=RE*UT_RC*0.5*(Y(J)+Y(J+1))
      IF (YP_RC.GE.YP_IN) THEN 
      Y2=YP_RC
      Y1=RE*UT_RC*0.5*(Y(J)+Y(J-1))
      J_RC=J
      GOTO 51
      ENDIF
50    CONTINUE
51    CONTINUE
      U1=UF_RC(NV,J_RC-1,K)
      U2=UF_RC(NV,J_RC,K)
      UF_IN(NV,1,JJ,K)=((U2-U1)/(Y2-Y1)*(YP_IN-Y2)+U2)
     >                *GAMMA
555   CONTINUE

C     LENGTH SCALE : OUTER VARIABLES

      ETHA_IN=0.5*(Y(JJ)+Y(JJ+1))/DELTA_IN
      ETHA_RC_MAX=0.5*(Y(N2M)+Y(N2))/DELTA_RC

      IF (ETHA_IN.GE.ETHA_RC_MAX) THEN 
      UF_IN(NV,2,JJ,K)=0.0         
      GOTO 666
      ENDIF

      DO 60 J=1,N2M    
      ETHA_RC=0.5*(Y(J)+Y(J+1))/DELTA_RC
      IF (ETHA_RC.GE.ETHA_IN) THEN 
      Y2=ETHA_RC
      Y1=0.5*(Y(J)+Y(J-1))/DELTA_RC
      J_RC=J
      GOTO 61
      ENDIF
60    CONTINUE
61    CONTINUE
      U1=UF_RC(NV,J_RC-1,K)
      U2=UF_RC(NV,J_RC,K)
      UF_IN(NV,2,JJ,K)=((U2-U1)/(Y2-Y1)*(ETHA_IN-Y2)+U2)
     >                *GAMMA
666   CONTINUE   

C     WEIGHTED AVERAGE WITH FUNCTION, W
      W=WEIGHT(ETHA_IN)
      UF_IN(NV,1,JJ,K)=UF_IN(NV,1,JJ,K)*(1.-W)+UF_IN(NV,2,JJ,K)*W
      
3     CONTINUE

      RETURN
      END

C*************** GET_INFLOW ***********************
      SUBROUTINE GET_INFLOW(UM_IN,UF_IN)
      INCLUDE 'PARAM.H'
      COMMON/PARA/RE
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/SIZE/ALX,ALY,ALZ
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/BCON/UBC1(M1,M3,3),UBC2(M1,M3,3)            
     >           ,UBC3(M2,M3,3),UBC4(M2,M3,3)
      COMMON/RECYCLE1/IRC,UT_IN,UT_RC,DELTA_RC
      COMMON/RECYCLE2/THETA_IN,THETA_RC 
      COMMON/RECYCLE3/UMEAN(3,0:M1,0:M2),UM_RC(3,0:M2)

      REAL UM_IN(3,2,0:M2)
      REAL UF_IN(3,2,0:M2,0:M3)
C          UF_IN(NV,ND,J ,K )
C                NV : VARIABLES
C                ND : NONDIMENSIONALIZED LENGTH SCALE
C                RETURN MEAN PROFILE AS UF_IN(NV,1,J,K)

      REAL DELSTA_LOCAL(M1)

      DO 20 NV=1,3
      DO 20 J=1,N2M
      DO 20 K=1,N3M
20    UBC3(J,K,NV)=UM_IN(NV,1,J)+UF_IN(NV,1,J,K) 

      DO K=1,N3M
      UBC3(1,K,2)=0.0 
      ENDDO

C     IMPOSE  V AT UPPER BOUNDARY
C     CALCULATE THE GROWTH RATE OF DISPLACEMENT THICKNESS
C     AND FIND V_INF

      DELSTA_LOCAL(1)=0.0
      DO 10 J=1,N2M
10    DELSTA_LOCAL(1)=DELSTA_LOCAL(1)+(1.-UM_IN(1,1,J))*DY(J)

      DO 30 I=2,N1M
      DELSTA_LOCAL(I)=0.0
      DO 30 J=1,N2M
30    DELSTA_LOCAL(I)=DELSTA_LOCAL(I)+(1.-UMEAN(1,I,J))*DY(J)

C     LINEAR REGRESSION OF DELSTA_LOCAL(I) AS V_INF*X(I)+B
C     EQ(25) IN LUND'S PAPER

      A11=REAL(N1M)  

      A12=0.0
      DO 32 I=1,N1M
32    A12=A12+X(I) 

      R_1=0.0
      DO 33 I=1,N1M
33    R_1=R_1+DELSTA_LOCAL(I) 

      A21=A12

      A22=0.0
      DO 34 I=1,N1M
34    A22=A22+X(I)**2

      R_2=0.0
      DO 35 I=1,N1M
35    R_2=R_2+X(I)*DELSTA_LOCAL(I) 

      V_INF=(R_1*A21-R_2*A11)/(A21*A12-A22*A11)
      B=(R_1*A22-R_2*A12)/(A22*A11-A12*A21)

      DO 40 I=1,N1M
      DO 40 K=1,N3M
40    UBC2(I,K,2)=V_INF

      WRITE(*,100) V_INF
100   FORMAT('V_INF =',E13.7)
      WRITE(*,*)

      RETURN
      END
   
C*************** EDDY_VISCOS_SM **************************************
      SUBROUTINE EDDY_VISCOS_SM(EV,Q1,Q2,Q3)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/CALC/ICONT
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/FINOUT/IDTOPT,NWRITE,NREAD,IAVG,IAVG_TRIG,NPRN,INSF,NINS

      REAL EV(0:M1,0:M2,0:M3)

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)

      REAL SR11(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR22(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR33(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR12(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR13(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR23(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)

      WRITE(*,*) '- SMAGORINSKY SGS MODEL -'
      WRITE(*,*) 


      CS=0.01

      DO 1 J=1,N2M

      CALL STRAIN_RATE(Q1,Q2,Q3,SR11,SR22,SR33,SR12,SR13,SR23,J)

      DO 1 K=1,N3M
      DO 1 I=1,N1M

         DELTA = (DX(I)*(1.0/DX3)*DY(J))**(1./3.)
         DAMPING = 1.0-EXP(-Y(J)*RE/22.0/25) ! VAN DRIEST DAMPING FUNCTION

         ABS_STRAIN =(2.*SR11(I,K)**2  
     >              + 2.*SR22(I,K)**2  
     >              + 2.*SR33(I,K)**2  
     >              + 4.*SR12(I,K)**2  
     >              + 4.*SR13(I,K)**2  
     >              + 4.*SR23(I,K)**2)**0.5

         EV(I,J,K)= CS*(DELTA*DAMPING)**2*ABS_STRAIN
1     CONTINUE


      DO K=1,N3M
      DO J=1,N2M
         EV(0 ,J,K)=EV(1,J,K)
         EV(N1,J,K)=EV(N1M,J,K)
      ENDDO
      ENDDO

      DO K=1,N3M
      DO I=0,N1
         EV(I,0 ,K)=EV(I,1  ,K)  
         EV(I,N2,K)=EV(I,N2M,K)
      ENDDO
      ENDDO


      IF(IAVG_TRIG.EQ.1.AND.MOD(NTIME,1).EQ.0) THEN       ! AVERAGING PER 1 TIME STEP
         if (MT.eq.0) then 
            call time_avg_les(EV,Q1,Q2,Q3)
         else 
            call phase_avg_les(EV,Q1,Q2,Q3)
         endif
      ENDIF

      RETURN
      END


C*************** EDDY_VISCOS_DM **************************************
      SUBROUTINE EDDY_VISCOS_DM(EV,Q1,Q2,Q3)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/CALC/ICONT
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/FINOUT/IDTOPT,NWRITE,NREAD,IAVG,IAVG_TRIG,NPRN,INSF,NINS

      REAL EV(0:M1,0:M2,0:M3)

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)

      REAL U_TILDE1(M1,M3)
      REAL U_TILDE2(M1,M3)
      REAL U_TILDE3(M1,M3)

      REAL UU_TILDE11(M1,M3)
      REAL UU_TILDE22(M1,M3)
      REAL UU_TILDE33(M1,M3)
      REAL UU_TILDE12(M1,M3)
      REAL UU_TILDE13(M1,M3)
      REAL UU_TILDE23(M1,M3)

      REAL SR11(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR22(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR33(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR12(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR13(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR23(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)

      REAL ALPHA11(M1,M3)  
      REAL ALPHA22(M1,M3)  
      REAL ALPHA33(M1,M3)  
      REAL ALPHA12(M1,M3)  
      REAL ALPHA13(M1,M3)  
      REAL ALPHA23(M1,M3)  

      REAL BETA11(M1,M3)  
      REAL BETA22(M1,M3)  
      REAL BETA33(M1,M3)  
      REAL BETA12(M1,M3)  
      REAL BETA13(M1,M3)  
      REAL BETA23(M1,M3)  

      REAL L11(M1,M3)    ! ANISOTROPIC PART OF RESOLVED SGS STRESS TENSOR
      REAL L22(M1,M3)    ! ANISOTROPIC PART OF RESOLVED SGS STRESS TENSOR
      REAL L33(M1,M3)    ! ANISOTROPIC PART OF RESOLVED SGS STRESS TENSOR
      REAL L12(M1,M3)    ! ANISOTROPIC PART OF RESOLVED SGS STRESS TENSOR
      REAL L13(M1,M3)    ! ANISOTROPIC PART OF RESOLVED SGS STRESS TENSOR
      REAL L23(M1,M3)    ! ANISOTROPIC PART OF RESOLVED SGS STRESS TENSOR

      REAL CD(M1,M2)
      
      WRITE(*,*) '- DYNAMIC SGS MODEL -'
      WRITE(*,*) 

      I_AVG_DEN = 0
      I_CLIP=0

      DO 1000 J=1,N2M
         JP=J+1

C---  CALCULATE L^A_IJ=L_IJ-1/3*DELTA_IJ*L_KK : ANISOTROPIC PART OF RESOLVED SGS STRESS

!$omp parallel do default(shared) private(KP,IP)
      DO 21 K=1,N3M
         KP=KPA(K)
      DO 21 I=1,N1M
         IP=I+1

         U_TILDE1(I,K)=0.5*(Q1(I,J,K)+Q1(IP,J ,K ))
         U_TILDE2(I,K)=0.5*(Q2(I,J,K)+Q2(I ,JP,K ))
         U_TILDE3(I,K)=0.5*(Q3(I,J,K)+Q3(I ,J ,KP))

         UU_TILDE11(I,K)=U_TILDE1(I,K)*U_TILDE1(I,K)
         UU_TILDE22(I,K)=U_TILDE2(I,K)*U_TILDE2(I,K)
         UU_TILDE33(I,K)=U_TILDE3(I,K)*U_TILDE3(I,K)
         UU_TILDE12(I,K)=U_TILDE1(I,K)*U_TILDE2(I,K)
         UU_TILDE13(I,K)=U_TILDE1(I,K)*U_TILDE3(I,K)
         UU_TILDE23(I,K)=U_TILDE2(I,K)*U_TILDE3(I,K)

21    CONTINUE

      CALL TEST_FILTERING(U_TILDE1,U_TILDE1)
      CALL TEST_FILTERING(U_TILDE2,U_TILDE2)
      CALL TEST_FILTERING(U_TILDE3,U_TILDE3)

      CALL TEST_FILTERING(UU_TILDE11,UU_TILDE11)   
      CALL TEST_FILTERING(UU_TILDE22,UU_TILDE22)   
      CALL TEST_FILTERING(UU_TILDE33,UU_TILDE33)   
      CALL TEST_FILTERING(UU_TILDE12,UU_TILDE12)   
      CALL TEST_FILTERING(UU_TILDE13,UU_TILDE13)   
      CALL TEST_FILTERING(UU_TILDE23,UU_TILDE23)   

!$omp parallel do default(shared) private(SUM_L_KK)
      DO 22 K=1,N3M
      DO 22 I=1,N1M

         L11(I,K)=UU_TILDE11(I,K)-U_TILDE1(I,K)*U_TILDE1(I,K)
         L22(I,K)=UU_TILDE22(I,K)-U_TILDE2(I,K)*U_TILDE2(I,K)
         L33(I,K)=UU_TILDE33(I,K)-U_TILDE3(I,K)*U_TILDE3(I,K)
         L12(I,K)=UU_TILDE12(I,K)-U_TILDE1(I,K)*U_TILDE2(I,K)
         L13(I,K)=UU_TILDE13(I,K)-U_TILDE1(I,K)*U_TILDE3(I,K)
         L23(I,K)=UU_TILDE23(I,K)-U_TILDE2(I,K)*U_TILDE3(I,K)

         SUM_L_KK  =L11(I,K)+L22(I,K)+L33(I,K)
         
         L11(I,K)=L11(I,K)-1./3.*SUM_L_KK
         L22(I,K)=L22(I,K)-1./3.*SUM_L_KK
         L33(I,K)=L33(I,K)-1./3.*SUM_L_KK
                  
22    CONTINUE


C---  CALCULATE BETA_IJ

      CALL STRAIN_RATE(Q1,Q2,Q3,SR11,SR22,SR33,SR12,SR13,SR23,J)

      CALL TEST_FILTERING(SR11,SR11)
      CALL TEST_FILTERING(SR22,SR22)
      CALL TEST_FILTERING(SR33,SR33)
      CALL TEST_FILTERING(SR12,SR12)
      CALL TEST_FILTERING(SR13,SR13)
      CALL TEST_FILTERING(SR23,SR23)


!$omp parallel do default(shared) private(TDELTA,ABS_STRAIN)
      DO 12 K=1,N3M
      DO 12 I=1,N1M

         TDELTA = ((2.0*DX(I))*(2.0/DX3)*(1.0*DY(J)))**(1./3.)
      
         ABS_STRAIN =(2.*SR11(I,K)**2  
     >              + 2.*SR22(I,K)**2  
     >              + 2.*SR33(I,K)**2  
     >              + 4.*SR12(I,K)**2  
     >              + 4.*SR13(I,K)**2  
     >              + 4.*SR23(I,K)**2)**0.5

         BETA11(I,K)=TDELTA**2*ABS_STRAIN*SR11(I,K)
         BETA22(I,K)=TDELTA**2*ABS_STRAIN*SR22(I,K)
         BETA33(I,K)=TDELTA**2*ABS_STRAIN*SR33(I,K)
         BETA12(I,K)=TDELTA**2*ABS_STRAIN*SR12(I,K)
         BETA13(I,K)=TDELTA**2*ABS_STRAIN*SR13(I,K)
         BETA23(I,K)=TDELTA**2*ABS_STRAIN*SR23(I,K)
                  
12    CONTINUE

C---  CALCULATE ALPHA_IJ

      CALL STRAIN_RATE(Q1,Q2,Q3,SR11,SR22,SR33,SR12,SR13,SR23,J)

!$omp parallel do default(shared) private(GDELTA,ABS_STRAIN)
      DO 11 K=1,N3M
      DO 11 I=1,N1M
      
         GDELTA = ((1.0*DX(I))*(1.0/DX3)*(1.0*DY(J)))**(1./3.)
         
         ABS_STRAIN =(2.*SR11(I,K)**2  
     >              + 2.*SR22(I,K)**2  
     >              + 2.*SR33(I,K)**2  
     >              + 4.*SR12(I,K)**2  
     >              + 4.*SR13(I,K)**2  
     >              + 4.*SR23(I,K)**2)**0.5

         ALPHA11(I,K)=GDELTA**2*ABS_STRAIN*SR11(I,K)
         ALPHA22(I,K)=GDELTA**2*ABS_STRAIN*SR22(I,K)
         ALPHA33(I,K)=GDELTA**2*ABS_STRAIN*SR33(I,K)
         ALPHA12(I,K)=GDELTA**2*ABS_STRAIN*SR12(I,K)
         ALPHA13(I,K)=GDELTA**2*ABS_STRAIN*SR13(I,K)
         ALPHA23(I,K)=GDELTA**2*ABS_STRAIN*SR23(I,K)
                           
11    CONTINUE


C---  CALCULATE THE SGS MODEL COEFFICIENT

      CALL TEST_FILTERING(ALPHA11,ALPHA11) 
      CALL TEST_FILTERING(ALPHA22,ALPHA22) 
      CALL TEST_FILTERING(ALPHA33,ALPHA33) 
      CALL TEST_FILTERING(ALPHA12,ALPHA12) 
      CALL TEST_FILTERING(ALPHA13,ALPHA13) 
      CALL TEST_FILTERING(ALPHA23,ALPHA23) 


      DO 5 I=1,N1M

         AVG_NUM=0.0
         AVG_DEN=0.0
      DO K=1,N3M
         AVG_NUM = AVG_NUM
     >           + L11(I,K)*(BETA11(I,K)-ALPHA11(I,K))
     >           + L22(I,K)*(BETA22(I,K)-ALPHA22(I,K))
     >           + L33(I,K)*(BETA33(I,K)-ALPHA33(I,K))     
     >           + L12(I,K)*(BETA12(I,K)-ALPHA12(I,K))*2.0
     >           + L13(I,K)*(BETA13(I,K)-ALPHA13(I,K))*2.0
     >           + L23(I,K)*(BETA23(I,K)-ALPHA23(I,K))*2.0     
         AVG_DEN = AVG_DEN
     >           + (BETA11(I,K)-ALPHA11(I,K))**2
     >           + (BETA22(I,K)-ALPHA22(I,K))**2
     >           + (BETA33(I,K)-ALPHA33(I,K))**2
     >           + (BETA12(I,K)-ALPHA12(I,K))**2*2.0
     >           + (BETA13(I,K)-ALPHA13(I,K))**2*2.0
     >           + (BETA23(I,K)-ALPHA23(I,K))**2*2.0
      ENDDO
      
         AVG_NUM=AVG_NUM/REAL(N3M)
         AVG_DEN=AVG_DEN/REAL(N3M)
      
	   IF (ABS(AVG_DEN).LT.1.E-12) THEN
	       AVG_DEN=1.E-12
	       I_AVG_DEN=I_AVG_DEN+1
         ENDIF

         CD(I,J)=-0.5*AVG_NUM/AVG_DEN

5     CONTINUE
 
C---  CALCULATE THE EDDY VISCOSITY AT N TH TIME STEP      

      DO 1 K=1,N3M
      DO 1 I=1,N1M
      
         GDELTA = ((1.0*DX(I))*(1.0/DX3)*(1.0*DY(J)))**(1./3.)

         ABS_STRAIN =(2.*SR11(I,K)**2  
     >              + 2.*SR22(I,K)**2  
     >              + 2.*SR33(I,K)**2  
     >              + 4.*SR12(I,K)**2  
     >              + 4.*SR13(I,K)**2  
     >              + 4.*SR23(I,K)**2)**0.5
     
         EV(I,J,K)= CD(I,J)*GDELTA**2*ABS_STRAIN

         IF (EV(I,J,K).LT.-1.0/RE) I_CLIP=I_CLIP+1
         EV(I,J,K)=AMAX1(EV(I,J,K),-1.0/RE)       ! CLIPPING
           
1     CONTINUE

1000  CONTINUE


	WRITE(*,100) 100.*REAL(I_AVG_DEN)/REAL(N1M)/REAL(N2M)
100   FORMAT('NON-ZERO DENOMINATOR TREATMENT FOR C : ',F10.3,'%')
      WRITE(*,101) 100.*REAL(I_CLIP)/REAL(N1M)/REAL(N2M)/REAL(N3M)
101   FORMAT('CLIPPING OPERATION  : ',F5.3,'%')

!$omp parallel do default(shared)
      DO K=1,N3M
      DO J=1,N2M
         EV(0 ,J,K)=EV(1  ,J,K)
         EV(N1,J,K)=EV(N1M,J,K)
      ENDDO
      ENDDO

!$omp parallel do default(shared)
      DO K=1,N3M
      DO I=0,N1
         EV(I,0 ,K)=EV(I,1  ,K)  
         EV(I,N2,K)=EV(I,N2M,K)
      ENDDO
      ENDDO

      SGS_MAX=0.0

      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M

         IF (EV(I,J,K).GT.SGS_MAX) THEN
             SGS_MAX=EV(I,J,K)
             I_MAX=I  
             J_MAX=J  
             K_MAX=K  
         ENDIF

      ENDDO
      ENDDO
      ENDDO

      WRITE(*,102) SGS_MAX*RE,I_MAX,J_MAX,K_MAX
102   FORMAT('MAXIMUM NUT/NU = ',E12.5,' AT I,J,K=',3(I3,X))
      WRITE(*,*) 

      IF(IAVG_TRIG.EQ.1.AND.MOD(NTIME,1).EQ.0) THEN       ! AVERAGING PER 1 TIME STEP
         if (MT.eq.0) then 
            call time_avg_les(EV,Q1,Q2,Q3)
         else 
            call phase_avg_les(EV,Q1,Q2,Q3)
         endif
      ENDIF

      RETURN
      END

C*************** STRAIN_RATE **************************************
      SUBROUTINE STRAIN_RATE(Q1,Q2,Q3
     >                      ,SR11,SR22,SR33,SR12,SR13,SR23,JJ)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/SIZE/ALX,ALY,ALZ

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)
      REAL SR11(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR22(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR33(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR12(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR13(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR23(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)

c      DO 1 J=1,N2M
           J =JJ
           JP=J+1
           JM=J-1

!$omp parallel do private(IP,IM,KP,KM,U1,U2,V1,V2,W1,W2)
!$omp&            shared(SR11,SR22,SR33,SR12,SR13,SR23)

      DO 1 K=1,N3M
           KP=KPA(K)
           KM=KMA(K)
      DO 1 I=1,N1M
           IP=I+1
           IM=I-1

      SR11(I,K) = (Q1(IP,J,K)-Q1(I ,J,K))/DX(I)
      SR22(I,K) = (Q2(I,JP,K)-Q2(I,J ,K))/DY(J)
      SR33(I,K) = (Q3(I,J,KP)-Q3(I,J,K ))*DX3


C---

      U1=1.0/H(J )*( 
     >    0.5*DY(JM)/DX(I)*(0.5*G(I)*Q1(IP,J ,K)+0.5*G(IP)*Q1(I,J ,K))
     >   +0.5*DY(J )/DX(I)*(0.5*G(I)*Q1(IP,JM,K)+0.5*G(IP)*Q1(I,JM,K))
     >   )

      U2=1.0/H(JP)*( 
     >    0.5*DY(J )/DX(I)*(0.5*G(I)*Q1(IP,JP,K)+0.5*G(IP)*Q1(I,JP,K))
     >   +0.5*DY(JP)/DX(I)*(0.5*G(I)*Q1(IP,J ,K)+0.5*G(IP)*Q1(I,J ,K))
     >   )

      V1=0.5*(
     >   1.0/G(I )*(0.5*DX(IM)*Q2(I ,JP,K)+0.5*DX(I )*Q2(IM,JP,K))
     >  +1.0/G(I )*(0.5*DX(IM)*Q2(I ,J ,K)+0.5*DX(I )*Q2(IM,J ,K))
     >  )

      V2=0.5*(
     >   1.0/G(IP)*(0.5*DX(I )*Q2(IP,JP,K)+0.5*DX(IP)*Q2(I ,JP,K))
     >  +1.0/G(IP)*(0.5*DX(I )*Q2(IP,J ,K)+0.5*DX(IP)*Q2(I ,J ,K))
     >  )

      SR12(I,K) = 0.5*((U2-U1)/DY(J)+(V2-V1)/DX(I))

C---


      U1=0.5*(
     >   1.0/DX(I)*(0.5*G(I)*Q1(IP,J,K )+0.5*G(IP)*Q1(I,J,K ))
     >  +1.0/DX(I)*(0.5*G(I)*Q1(IP,J,KM)+0.5*G(IP)*Q1(I,J,KM))
     >  )

      U2=0.5*(
     >   1.0/DX(I)*(0.5*G(I)*Q1(IP,J,KP)+0.5*G(IP)*Q1(I,J,KP))
     >  +1.0/DX(I)*(0.5*G(I)*Q1(IP,J,K )+0.5*G(IP)*Q1(I,J,K ))
     >  )

      W1=0.5*(
     >   1.0/G(I )*(0.5*DX(IM)*Q3(I ,J,KP)+0.5*DX(I )*Q3(IM,J,KP))
     >  +1.0/G(I )*(0.5*DX(IM)*Q3(I ,J,K )+0.5*DX(I )*Q3(IM,J,K ))
     >  )

      W2=0.5*(
     >   1.0/G(IP)*(0.5*DX(I )*Q3(IP,J,KP)+0.5*DX(IP)*Q3(I ,J,KP))
     >  +1.0/G(IP)*(0.5*DX(I )*Q3(IP,J,K )+0.5*DX(IP)*Q3(I ,J,K ))
     >  )

      SR13(I,K) = 0.5*((U2-U1)*DX3+(W2-W1)/DX(I))

C---
      W1=1.0/H(J )*(0.5*DY(JM)*0.5*(Q3(I,J ,K)+Q3(I,J ,KP))
     >             +0.5*DY(J )*0.5*(Q3(I,JM,K)+Q3(I,JM,KP)))

      W2=1.0/H(JP)*(0.5*DY(J )*0.5*(Q3(I,JP,K)+Q3(I,JP,KP))
     >             +0.5*DY(JP)*0.5*(Q3(I,J ,K)+Q3(I,J ,KP)))

      V1=0.25*(Q2(I,J,K )+Q2(I,JP,K )+Q2(I,JP,KM)+Q2(I,J ,KM))
      V2=0.25*(Q2(I,J,KP)+Q2(I,JP,KP)+Q2(I,JP,K )+Q2(I,J ,K ))
       
      SR23(I,K) = 0.5*((W2-W1)/DY(J)+(V2-V1)*DX3)

1     CONTINUE

      RETURN
      END

C*************** TEST_FILTERING ************************************
      SUBROUTINE TEST_FILTERING(FI,FO)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/SIZE/ALX,ALY,ALZ

      REAL FI(M1,M3)
      REAL FO(M1,M3)
      REAL TEMP(M1)

C---  FILTERING IN X3 DIRECTION
      DO  K=1,N3M
           KP=KPA(K)
           KM=KMA(K)
      DO  I=1,N1M
           FO(I,K)=(FI(I,KP)+4.0*FI(I,K)+FI(I,KM))/6.0 
      ENDDO
      ENDDO

C---  FILTERING IN X1 DIRECTION
      DO K=1,N3M
         DO IC=1,N1M            
            IP=MIN(IC+1,N1M)
            IM=MAX(IC-1,1)
            TEMP(IC)=(FO(IP,K)+4.0*FO(IC,K)+FO(IM,K))/6.0 
         ENDDO
         DO I=1,N1M
            FO(I,K)=TEMP(I)
         ENDDO

      ENDDO

      RETURN
      END
 

C*************** time_avg_les **********************
      SUBROUTINE time_avg_les(EV,Q1,Q2,Q3)

      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/PARA/RE
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/FINOUT/IDTOPT,NWRITE,NREAD,IAVG,IAVG_TRIG,NPRN,INSF,NINS
      COMMON/CALC/ICONT
      common/phaseavg/nphavg(0:mt),nphavg_les(0:mt),ntime_ctd

      common/time_les/TAU_t(M1,M2),EV_ZT_t(M1,M2),EV2_ZT_t(M1,M2)
     >     ,E_SGS_t(M1,M2),E_MOL_t(M1,M2)           ! AKSELVOLL & MOIN TF-63 PP.71
     >     ,F_NEGA_EV_t(M1,M2),F_NEGA_TV_t(M1,M2)   ! AKSELVOLL & MOIN TF-63 PP.70

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)

      REAL EV(0:M1,0:M2,0:M3)
      REAL SR11(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR22(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR33(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR12(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR13(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR23(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)


      iphase = mod(ntime_ctd,mt+1)

      if (nphavg_les(iphase).eq.0) then      ! to avoid NaN for the initial value 

        EV_ZT_t(:,:)    = 0.0
        EV2_ZT_t(:,:)   = 0.0
        TAU_t(:,:)      = 0.0
        E_SGS_t(:,:)    = 0.0
        E_MOL_t(:,:)    = 0.0
        F_NEGA_EV_t(:,:)= 0.0
        F_NEGA_TV_t(:,:)= 0.0

      endif

      if(ICONT.eq.2) then
        call read_phaseavg_les(
     &                  tau_t,ev_zt_t,ev2_zt_t,
     &                  e_sgs_t,e_mol_t,f_nega_ev_t,f_nega_tv_t)
c        ICONT = -1    ! to read the prevous file only at firt step
                       ! this treatment will be done at time_avg.
      endif

      nphavg_les(iphase)=nphavg_les(iphase)+1

      write(*,*) 'les_avg',iphase,nphavg_les(iphase)

      RNAVG_LES=real(NPHAVG_les(IPHASE))


      DO 1 J=1,N2M

      CALL STRAIN_RATE(Q1,Q2,Q3,SR11,SR22,SR33,SR12,SR13,SR23,J)

      DO 1 I=1,N1M

         EV_Z = 0.0
         EV2_Z = 0.0
         TAU_Z = 0.0
         E_SGS_Z = 0.0
         E_MOL_Z = 0.0
         F_NEGA_EV_Z = 0.0
         F_NEGA_TV_Z = 0.0

         DO K=1,N3M

         EV_Z  = EV_Z  + EV(I,J,K)/REAL(N3M)
         EV2_Z  = EV2_Z  + EV(I,J,K)**2/REAL(N3M)
         TAU_Z = TAU_Z - 2.0*EV(I,J,K)*SR12(I,K)/REAL(N3M)

         ABS_STRAIN =(2.*SR11(I,K)**2  
     >              + 2.*SR22(I,K)**2  
     >              + 2.*SR33(I,K)**2  
     >              + 4.*SR12(I,K)**2  
     >              + 4.*SR13(I,K)**2  
     >              + 4.*SR23(I,K)**2)**0.5

         E_SGS_Z = E_SGS_Z + EV(I,J,K)*ABS_STRAIN**2/REAL(N3M)
         E_MOL_Z = E_MOL_Z + ABS_STRAIN**2/RE/REAL(N3M)

         IF (EV(I,J,K).LT.-1.0/RE) 
     >       F_NEGA_TV_Z=F_NEGA_TV_Z+1.0/REAL(N3M)
         IF (EV(I,J,K).LT. 0.0   ) 
     >       F_NEGA_EV_Z=F_NEGA_EV_Z+1.0/REAL(N3M)

         ENDDO

      EV_ZT_t(I,J) =(EV_ZT_t(I,J) *(RNAVG_LES-1.0)+EV_Z   )/RNAVG_LES
      EV2_ZT_t(I,J)=(EV2_ZT_t(I,J)*(RNAVG_LES-1.0)+EV2_Z  )/RNAVG_LES
      TAU_t(I,J)   =(TAU_t(I,J)   *(RNAVG_LES-1.0)+TAU_Z  )/RNAVG_LES
      E_SGS_t(I,J) =(E_SGS_t(I,J) *(RNAVG_LES-1.0)+E_SGS_Z)/RNAVG_LES
      E_MOL_t(I,J) =(E_MOL_t(I,J) *(RNAVG_LES-1.0)+E_MOL_Z)/RNAVG_LES
      F_NEGA_EV_t(I,J)=(F_NEGA_EV_t(I,J)*(RNAVG_LES-1.0)+F_NEGA_EV_Z)
     >               /RNAVG_LES
      F_NEGA_TV_t(I,J)=(F_NEGA_TV_t(I,J)*(RNAVG_LES-1.0)+F_NEGA_TV_Z)
     >               /RNAVG_LES


1     CONTINUE

      IF(MOD(NTIME,NPRN).EQ.0.AND.NWRITE.EQ.1) THEN
        call write_phaseavg_les(
     &                  tau_t,ev_zt_t,ev2_zt_t,
     &                  e_sgs_t,e_mol_t,f_nega_ev_t,f_nega_tv_t)
      endif


      RETURN
      END

C*************** phase_avg_les **********************
      SUBROUTINE phase_avg_les(EV,Q1,Q2,Q3)

      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/PARA/RE
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      COMMON/CALC/ICONT
      common/phaseavg/nphavg(0:mt),nphavg_les(0:mt),ntime_ctd

      real  TAU(M1,M2),EV_ZT(M1,M2),EV2_ZT(M1,M2)
     >     ,E_SGS(M1,M2),E_MOL(M1,M2)           ! AKSELVOLL & MOIN TF-63 PP.71
     >     ,F_NEGA_EV(M1,M2),F_NEGA_TV(M1,M2)   ! AKSELVOLL & MOIN TF-63 PP.70

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)

      REAL EV(0:M1,0:M2,0:M3)
      REAL SR11(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR22(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR33(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR12(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR13(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)
      REAL SR23(M1,M3)   ! STRAIN RATE SR_IJ=1/2*(DU_I/DX_J+DU_J/DX_I)


      iphase = mod(ntime_ctd,mt+1)

      if (nphavg_les(iphase).eq.0) then      ! to avoid NaN for the initial value 

        EV_ZT(:,:)    = 0.0
        EV2_ZT(:,:)   = 0.0
        TAU(:,:)      = 0.0
        E_SGS(:,:)    = 0.0
        E_MOL(:,:)    = 0.0
        F_NEGA_EV(:,:)= 0.0
        F_NEGA_TV(:,:)= 0.0

      endif

      if(ICONT.eq.2 .or. nphavg_les(iphase).ne.0) then
        call read_phaseavg_les(
     &                  tau,ev_zt,ev2_zt,
     &                  e_sgs,e_mol,f_nega_ev,f_nega_tv)
      endif

      nphavg_les(iphase)=nphavg_les(iphase)+1

      write(*,*) 'les_avg',iphase,nphavg_les(iphase)

      RNAVG_LES=real(NPHAVG_les(IPHASE))


      DO 1 J=1,N2M

      CALL STRAIN_RATE(Q1,Q2,Q3,SR11,SR22,SR33,SR12,SR13,SR23,J)

      DO 1 I=1,N1M

         EV_Z = 0.0
         EV2_Z = 0.0
         TAU_Z = 0.0
         E_SGS_Z = 0.0
         E_MOL_Z = 0.0
         F_NEGA_EV_Z = 0.0
         F_NEGA_TV_Z = 0.0

         DO K=1,N3M

         EV_Z  = EV_Z  + EV(I,J,K)/REAL(N3M)
         EV2_Z  = EV2_Z  + EV(I,J,K)**2/REAL(N3M)
         TAU_Z = TAU_Z - 2.0*EV(I,J,K)*SR12(I,K)/REAL(N3M)

         ABS_STRAIN =(2.*SR11(I,K)**2  
     >              + 2.*SR22(I,K)**2  
     >              + 2.*SR33(I,K)**2  
     >              + 4.*SR12(I,K)**2  
     >              + 4.*SR13(I,K)**2  
     >              + 4.*SR23(I,K)**2)**0.5

         E_SGS_Z = E_SGS_Z + EV(I,J,K)*ABS_STRAIN**2/REAL(N3M)
         E_MOL_Z = E_MOL_Z + ABS_STRAIN**2/RE/REAL(N3M)

         IF (EV(I,J,K).LT.-1.0/RE) 
     >       F_NEGA_TV_Z=F_NEGA_TV_Z+1.0/REAL(N3M)
         IF (EV(I,J,K).LT. 0.0   ) 
     >       F_NEGA_EV_Z=F_NEGA_EV_Z+1.0/REAL(N3M)

         ENDDO

      EV_ZT(I,J) =(EV_ZT(I,J) *(RNAVG_LES-1.0)+EV_Z   )/RNAVG_LES
      EV2_ZT(I,J)=(EV2_ZT(I,J)*(RNAVG_LES-1.0)+EV2_Z  )/RNAVG_LES
      TAU(I,J)   =(TAU(I,J)   *(RNAVG_LES-1.0)+TAU_Z  )/RNAVG_LES
      E_SGS(I,J) =(E_SGS(I,J) *(RNAVG_LES-1.0)+E_SGS_Z)/RNAVG_LES
      E_MOL(I,J) =(E_MOL(I,J) *(RNAVG_LES-1.0)+E_MOL_Z)/RNAVG_LES
      F_NEGA_EV(I,J)=(F_NEGA_EV(I,J)*(RNAVG_LES-1.0)+F_NEGA_EV_Z)
     >               /RNAVG_LES
      F_NEGA_TV(I,J)=(F_NEGA_TV(I,J)*(RNAVG_LES-1.0)+F_NEGA_TV_Z)
     >               /RNAVG_LES


1     CONTINUE

        call write_phaseavg_les(
     &                  tau,ev_zt,ev2_zt,
     &                  e_sgs,e_mol,f_nega_ev,f_nega_tv)


      RETURN
      END

C*************** read_phaseavg_les **********************
      SUBROUTINE read_phaseavg_les(
     &                  tau,ev_zt,ev2_zt,
     &                  e_sgs,e_mol,f_nega_ev,f_nega_tv)

      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/PARA/RE
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      common/phaseavg/nphavg(0:mt),nphavg_les(0:mt),ntime_ctd

      real TAU(M1,M2),EV_ZT(M1,M2),EV2_ZT(M1,M2)
     >    ,E_SGS(M1,M2),E_MOL(M1,M2)           ! AKSELVOLL & MOIN TF-63 PP.71
     >    ,F_NEGA_EV(M1,M2),F_NEGA_TV(M1,M2)   ! AKSELVOLL & MOIN TF-63 PP.70

      character*30 phasefile 

      iphase = mod(ntime_ctd,mt+1)

      phasefile='phaseles.'
      nn=index(phasefile,'.')
      write(unit=phasefile(nn+1:),fmt='(bn,i4.4)') iphase

      OPEN(11,FILE=phasefile,STATUS='UNKNOWN',form='unformatted')
      
         read(11) nphavg_les(iphase) 
      DO I=1,N1M
        do J=1,N2M
         read(11) TAU(I,J)
         read(11) EV_ZT(I,J),EV2_ZT(I,J)
         read(11) E_SGS(I,J),E_MOL(I,J)
         read(11) F_NEGA_EV(I,J),F_NEGA_TV(I,J)
        enddo
      enddo

      CLOSE(11)

      RETURN
      END

C*************** write_phaseavg_les **********************
      SUBROUTINE write_phaseavg_les(
     &                  tau,ev_zt,ev2_zt,
     &                  e_sgs,e_mol,f_nega_ev,f_nega_tv)

      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/PARA/RE
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
      common/phaseavg/nphavg(0:mt),nphavg_les(0:mt),ntime_ctd

      real TAU(M1,M2),EV_ZT(M1,M2),EV2_ZT(M1,M2)
     >    ,E_SGS(M1,M2),E_MOL(M1,M2)           ! AKSELVOLL & MOIN TF-63 PP.71
     >    ,F_NEGA_EV(M1,M2),F_NEGA_TV(M1,M2)   ! AKSELVOLL & MOIN TF-63 PP.70

      character*30 phasefile 

      iphase = mod(ntime_ctd,mt+1)

      phasefile='phaseles.'
      nn=index(phasefile,'.')
      write(unit=phasefile(nn+1:),fmt='(bn,i4.4)') iphase

      OPEN(11,FILE=phasefile,STATUS='UNKNOWN',form='unformatted')
      
         write(11) nphavg_les(iphase) 
      DO I=1,N1M
        do J=1,N2M
         WRITE(11) TAU(I,J)
         WRITE(11) EV_ZT(I,J),EV2_ZT(I,J)
         WRITE(11) E_SGS(I,J),E_MOL(I,J)
         WRITE(11) F_NEGA_EV(I,J),F_NEGA_TV(I,J)
        enddo
      enddo

      CLOSE(11)

      RETURN
      END

C*************** EV_INTERPOLATION ********************************
      SUBROUTINE EV_INTERPOLATION(EV,EVI,EVJ,EVK)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      REAL EV(0:M1,0:M2,0:M3)   
      REAL EVI(0:M1,0:M2,0:M3)     ! EV ON I AXIS
      REAL EVJ(0:M1,0:M2,0:M3)     ! EV ON J AXIS 
      REAL EVK(0:M1,0:M2,0:M3)     ! EV ON K AXIS 

C--------------------------------------
      DO K=1,N3M
         KM=KMA(K)
      DO J=1,N2
         JM=J-1
      DO I=1,N1M

      EVI(I,J,K)=1.0/H(J)*(
     >            0.5*(EV(I,J ,K)+EV(I,J ,KM))*0.5*DY(JM)
     >           +0.5*(EV(I,JM,K)+EV(I,JM,KM))*0.5*DY(J )
     >           )
      ENDDO
      ENDDO
      ENDDO

C--------------------------------------
      DO K=1,N3M
         KM=KMA(K)
      DO J=1,N2M
      DO I=1,N1
         IM=I-1

      EVJ(I,J,K)=0.5*(
     >         1.0/G(I)*(0.5*DX(IM)*EV(I,J,K )+0.5*DX(I)*EV(IM,J,K ))
     >        +1.0/G(I)*(0.5*DX(IM)*EV(I,J,KM)+0.5*DX(I)*EV(IM,J,KM))
     >        )

      ENDDO
      ENDDO
      ENDDO

C--------------------------------------
      DO K=1,N3M
      DO J=1,N2
         JM=J-1
      DO I=1,N1
         IM=I-1

      EVK(I,J,K)=1.0/H(J)*(
     >  1.0/G(I)
     >    *(0.5*DX(IM)*EV(I,J ,K)+0.5*DX(I)*EV(IM,J ,K))*0.5*DY(JM)
     > +1.0/G(I)
     >    *(0.5*DX(IM)*EV(I,JM,K)+0.5*DX(I)*EV(IM,JM,K))*0.5*DY(J )
     >   )

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END

C*************** DUDX **************************************
      SUBROUTINE DUDX(Q1,Q2,Q3
     >               ,DUDX11,DUDX12,DUDX13
     >               ,DUDX21,DUDX22,DUDX23
     >               ,DUDX31,DUDX32,DUDX33,JJ)

      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/SIZE/ALX,ALY,ALZ

      REAL Q1(0:M1,0:M2,0:M3)
      REAL Q2(0:M1,0:M2,0:M3)
      REAL Q3(0:M1,0:M2,0:M3)

      REAL DUDX11(M1,M3),DUDX12(M1,M3),DUDX13(M1,M3)  ! VEL-GRADIENT AT CELL CENTER
      REAL DUDX21(M1,M3),DUDX22(M1,M3),DUDX23(M1,M3)
      REAL DUDX31(M1,M3),DUDX32(M1,M3),DUDX33(M1,M3)

C      DO 1 J=1,N2M
           J =JJ
           JP=J+1
           JM=J-1
      DO 1 K=1,N3M
           KP=KPA(K)
           KM=KMA(K)
      DO 1 I=1,N1M
           IP=I+1
           IM=I-1

C---
      DUDX11(I,K) = (Q1(IP,J,K)-Q1(I ,J,K))/DX(I)
      DUDX22(I,K) = (Q2(I,JP,K)-Q2(I,J ,K))/DY(J)
      DUDX33(I,K) = (Q3(I,J,KP)-Q3(I,J,K ))*DX3


C---
      U1=1.0/H(J )*( 
     >    0.5*DY(JM)/DX(I)*(0.5*G(I)*Q1(IP,J ,K)+0.5*G(IP)*Q1(I,J ,K))
     >   +0.5*DY(J )/DX(I)*(0.5*G(I)*Q1(IP,JM,K)+0.5*G(IP)*Q1(I,JM,K))
     >   )

      U2=1.0/H(JP)*( 
     >    0.5*DY(J )/DX(I)*(0.5*G(I)*Q1(IP,JP,K)+0.5*G(IP)*Q1(I,JP,K))
     >   +0.5*DY(JP)/DX(I)*(0.5*G(I)*Q1(IP,J ,K)+0.5*G(IP)*Q1(I,J ,K))
     >   )

      V1=0.5*(
     >   1.0/G(I )*(0.5*DX(IM)*Q2(I ,JP,K)+0.5*DX(I )*Q2(IM,JP,K))
     >  +1.0/G(I )*(0.5*DX(IM)*Q2(I ,J ,K)+0.5*DX(I )*Q2(IM,J ,K))
     >  )

      V2=0.5*(
     >   1.0/G(IP)*(0.5*DX(I )*Q2(IP,JP,K)+0.5*DX(IP)*Q2(I ,JP,K))
     >  +1.0/G(IP)*(0.5*DX(I )*Q2(IP,J ,K)+0.5*DX(IP)*Q2(I ,J ,K))
     >  )

      DUDX12(I,K) = (U2-U1)/DY(J)
      DUDX21(I,K) = (V2-V1)/DX(I)

C---
      U1=0.5*(
     >   1.0/DX(I)*(0.5*G(I)*Q1(IP,J,K )+0.5*G(IP)*Q1(I,J,K ))
     >  +1.0/DX(I)*(0.5*G(I)*Q1(IP,J,KM)+0.5*G(IP)*Q1(I,J,KM))
     >  )

      U2=0.5*(
     >   1.0/DX(I)*(0.5*G(I)*Q1(IP,J,KP)+0.5*G(IP)*Q1(I,J,KP))
     >  +1.0/DX(I)*(0.5*G(I)*Q1(IP,J,K )+0.5*G(IP)*Q1(I,J,K ))
     >  )

      W1=0.5*(
     >   1.0/G(I )*(0.5*DX(IM)*Q3(I ,J,KP)+0.5*DX(I )*Q3(IM,J,KP))
     >  +1.0/G(I )*(0.5*DX(IM)*Q3(I ,J,K )+0.5*DX(I )*Q3(IM,J,K ))
     >  )

      W2=0.5*(
     >   1.0/G(IP)*(0.5*DX(I )*Q3(IP,J,KP)+0.5*DX(IP)*Q3(I ,J,KP))
     >  +1.0/G(IP)*(0.5*DX(I )*Q3(IP,J,K )+0.5*DX(IP)*Q3(I ,J,K ))
     >  )

      DUDX13(I,K) = (U2-U1)*DX3
      DUDX31(I,K) = (W2-W1)/DX(I)

C---
      W1=1.0/H(J )*(0.5*DY(JM)*0.5*(Q3(I,J ,K)+Q3(I,J ,KP))
     >             +0.5*DY(J )*0.5*(Q3(I,JM,K)+Q3(I,JM,KP)))

      W2=1.0/H(JP)*(0.5*DY(J )*0.5*(Q3(I,JP,K)+Q3(I,JP,KP))
     >             +0.5*DY(JP)*0.5*(Q3(I,J ,K)+Q3(I,J ,KP)))

      V1=0.25*(Q2(I,J,K )+Q2(I,JP,K )+Q2(I,JP,KM)+Q2(I,J ,KM))
      V2=0.25*(Q2(I,J,KP)+Q2(I,JP,KP)+Q2(I,JP,K )+Q2(I,J ,K ))
       
      DUDX32(I,K) = (W2-W1)/DY(J)
      DUDX23(I,K) = (V2-V1)*DX3

1     CONTINUE

      RETURN
      END

c      include 'POISSON_MG_CXML.F'
      include 'POISSON_CS_CXML.F'
