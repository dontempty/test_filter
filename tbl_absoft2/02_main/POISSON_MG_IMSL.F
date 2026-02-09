C---------------------------------------------------------------------------
C
C     POISSON SOLVER WITH MULTIGRID
C
C     POISSON EQUATION IS SOLVED BY
C     1. FFT IN X3 DIRECTION
C     2. MULTIGRID IN X1 DIRECTION
C     3. TRIDIAGONAL SOLVER IN X2 DIRECTION
C---------------------------------------------------------------------------
C
C     NON-UNIFORM CODE
C
C     2003.11. 11.
C     KYOUNGYOUN KIM
C     FLOW CONTROL LABORATORY
C     DEPARTMENT OF MECHANICAL ENGINEERING
C     KAIST 
C-----------------------------------------------------------------------
C


C*************** INIPOISSON ***********************     
C
C     GET THE MODIFIED WAVE NUMBER
C     DEFINE BASIC MULTIGRID PARAMTER
C---------------------------------------------------------------------------

      SUBROUTINE INIPOISSON
      INCLUDE 'PARAM.H'
      PARAMETER (M3M=M3-1,M3MH=M3M/2+1)
      PARAMETER (NLEV=11)                    ! MULTIGRID

      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/WAVK3/AK3(M3MH)

      COMMON/CNLEV/N1G(0:NLEV),N12G(0:NLEV)
      COMMON/INLEV/IPLV(M1,NLEV),IMLV(M1,NLEV)
      COMMON/MLV/MLEV,IWMG,MLW,MLCMG,MHCMG,EPSM,RATIO,NSOR

      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      PI=ACOS(-1.0)
 
C     MODIFIED WAVE NUMBER DEFINITION NECESSARY FOR THE X3 DIRECTION
      DO 20 K=1,N3M/2+1
      AK3(K)=2.*(1.-COS((K-1)*2.*PI/REAL(N3M)))*DX3Q
20    CONTINUE

c---
c     N1G(NL) : total grid number at level NL
c     N12G(NL) : sum of grid numer from level 1 to NL
c
      N1G(0)=0
      N12G(0)=0
      DO 1 NL=1,MLEV
      N1G(NL)=N1M*2**FLOAT(-NL+1)
      N12G(NL)=(N1G(NL)+1)*N2M+N12G(NL-1)
      WRITE(*,333)NL,N1G(NL),N2M,N12G(NL)
333   FORMAT(3X,'NL=',I3,3(2X,I7))
c

c---
c     IMLV(IC,NL) : left index of IC at level NL
c     IPLV(IC,NL) : right index of IC at level NL
c
      DO 10 IC=1,N1G(NL)
      IMLV(IC,NL)=IC-1
      IPLV(IC,NL)=IC+1
10    CONTINUE
      IMLV(1,NL)=1
      IPLV(N1G(NL),NL)=N1G(NL)
1     CONTINUE
c
      CALL COPOISSON

      RETURN
      END


C---------------------------------------------------------------------------
C     SUBROUTINE COPOISSON
C
C     DEFINE COEFFICIENT IN POSSION EQUATION
C---------------------------------------------------------------------------
C     COEFFICIENTS FOR THE FIVE POINT STENCIL OF THE POISSON EQ.
C     1=(IP,JC);2=(IM,JC);3=(IC,JP);4=(IC,JM);5=(IC,JC)

      SUBROUTINE COPOISSON

      INCLUDE 'PARAM.H'
      PARAMETER (M1M=M1-1,M2M=M2-1)
      PARAMETER (M3M=M3-1,M3MH=M3M/2+1)
      PARAMETER (NLEV=11,NIJ=M1M*M2M+M1M*M2M/2+M1M*M2M/4
     >      +M1M*M2M/8+M1M*M2M/16+M1M*M2M/32+M1M*M2M/64+M1M*M2M/128
     >      +M1M*M2M/256+M1M*M2M/512+M1M*M2M/1024+M2M*NLEV)

      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/WAVK3/AK3(M3MH)

      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/COMG/CO(6,NIJ)

      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/MESH5/XM(M1,NLEV),DXX(0:M1,NLEV),GX(M1,NLEV)

      COMMON/CNLEV/N1G(0:NLEV),N12G(0:NLEV)
      COMMON/INLEV/IPLV(M1,NLEV),IMLV(M1,NLEV)
      COMMON/MLV/MLEV,IWMG,MLW,MLCMG,MHCMG,EPSM,RATIO,NSOR

      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)


      IGJG(I,J,L)=I+(J-1)*(N1G(L)+1)+N12G(L-1)


C---  DEFINE X GRID VALUE AT EACH MULTI GRID LEVEL

c---
c     XM(IC,NL)  : physical location corresponding of IC at level NL
c     DXX(IC,NL) : DX corresponding of IC at level NL
c     GX(IC,NL)  : GX corresponding of IC at level NL
c
      DO NL=1,MLEV

         DO IC=1,N1G(NL)
            ICJ=(IC-1)*2**(NL-1)+1   ! ICJ=I index at level 1 corresponding of IC at level NL   
            XM(IC,NL)=X(ICJ) 
         ENDDO
            XM(N1G(NL)+1,NL)=X(N1)

            DXX(0,NL)=0.0
         DO IC=1,N1G(NL)
            DXX(IC,NL)=XM(IC+1,NL)-XM(IC,NL)
         ENDDO
            DXX(N1G(NL)+1,NL)=0.0

         DO IC=1,N1G(NL)+1
            GX(IC,NL)=0.5*(DXX(IC,NL)+DXX(IC-1,NL))
         ENDDO

      ENDDO
c
C---  DEFINE THE COEFFICIENT OF POISSON EQ.

      DO 10 NL=1,MLEV

      DO 100 IC=1,N1G(NL)

      IUM=IC-IMLV(IC,NL)
      IUP=IPLV(IC,NL)-IC
      IP=IC+1
      DO 100 JC=1,N2M
      JP=JC+1
      JM=JC-1
      JUM=JC-JMU(JC)
      JUP=JPA(JC)-JC

      IJ=IGJG(IC,JC,NL)


      CO(1,IJ)=IUP*(1.0/DXX(IC,NL)/GX(IP,NL))
      CO(2,IJ)=IUM*(1.0/DXX(IC,NL)/GX(IC,NL))
      CO(3,IJ)=JUP*(1.0/DY(JC)/H(JP))
      CO(4,IJ)=JUM*(1.0/DY(JC)/H(JC))
      CO(5,IJ)=-CO(1,IJ)-CO(2,IJ)-CO(3,IJ)-CO(4,IJ)

100   CONTINUE
10    CONTINUE

      RETURN
      END

C   ********************* TAKEDP ********************************
C     CALCULATE DPH. 
C     PERIODIC DIRECTION ALONG X3 ALLOWS TO USE THE REAL FOURIER TRANSFORM.

      SUBROUTINE TAKEDP(DP,RDP)
C---  FOR CXML
C      INCLUDE 'CXML_INCLUDE.F90'
C---
      INCLUDE 'PARAM.H'
      PARAMETER (M1M=M1-1,M2M=M2-1)
      PARAMETER (M3M=M3-1,M3MH=M3M/2+1,M3MHM=M3M/2)
      PARAMETER (NLEV=11,NIJ=M1M*M2M+M1M*M2M/2+M1M*M2M/4
     >      +M1M*M2M/8+M1M*M2M/16+M1M*M2M/32+M1M*M2M/64+M1M*M2M/128
     >      +M1M*M2M/256+M1M*M2M/512+M1M*M2M/1024+M2M*NLEV)

      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/WAVK3/AK3(M3MH)

      COMMON/N1MEQ1/AJP(M2),AJC(M2,M3MH),AJM(M2)
      COMMON/DAK/N3MH,N3MD,N3MDU,N3MU,NDM

      COMMON/CNLEV/N1G(0:NLEV),N12G(0:NLEV)
      COMMON/POIS/IK0,IPH0

      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      DIMENSION DP(0:M1,0:M2,0:M3),RDP(0:M1,0:M2,0:M3)
      DIMENSION PHM(NIJ),REST(NIJ)

C---   FOR CXML
C      REAL FRDP_R(M1,M2,0:M3MHM)
C      REAL FRDP_I(M1,M2,0:M3MHM)
C      REAL FDP_R(M1,M2,0:M3MHM)
C      REAL FDP_I(M1,M2,0:M3MHM)
C      INTEGER*4 STATUS
C      
C      REAL*8 R2(0:N3M+2)
C---

C---  FOR IMSL
      REAL FRDP_R(M1,M2,0:M3MH)
      REAL FRDP_I(M1,M2,0:M3MH)
      REAL FDP_R(M1,M2,0:M3MH)
      REAL FDP_I(M1,M2,0:M3MH)

      REAL WFFTR(2*M3M+15)
      REAL RFIN(M3M)
      REAL RFEX(M3M)
C---
      
	
      IGJG(I,J,L)=I+(J-1)*(N1G(L)+1)+N12G(L-1)

      WRITE(*,*) 
      WRITE(*,*) 'MULTI-GRID PROCEDURE FOR PRESSURE'
        
      N3MH=N3M/2+1 

C     IN THE CASE N3M=1, WE DON'T HAVE TO DO FFT.
      IF (N3M.EQ.1) THEN
      DO 10 IC=1,N1M
      DO 10 JC=1,N2M
      IJ=IGJG(IC,JC,1)
      PHM(IJ)=DP(IC,JC,1)*(1.-REAL(IPH0))
      REST(IJ)=RDP(IC,JC,1)
10    CONTINUE
      CALL DSOLV(PHM,REST,1)
      DO 11 IC=1,N1M
      DO 11 JC=1,N2M
      IJ=IGJG(IC,JC,1)
      DP(IC,JC,1)=PHM(IJ)
11    CONTINUE
      GO TO 1000
      ENDIF

      PI=ACOS(-1.)

C     FFT OF RDP IN X3 DIRECTION

C---   FOR CXML
C      DO 110 I=1,N1M              
C      DO 110 J=1,N2M
C      DO 111 K=0,N3M-1
C      R2(K)=RDP(I,J,K+1)
C111   CONTINUE
C
C      STATUS=DFFT('R','C','F',R2,R2,N3M,1)
C      
C      DO 112 K=0,N3M/2-1
C      FRDP_R(I,J,K)=R2(K*2)
C      FRDP_I(I,J,K)=-R2(K*2+1)
C112   CONTINUE
C      FRDP_R(I,J,N3M/2)=R2(N3M)
C110   CONTINUE
C---

C---  FOR IMSL
      CALL DFFTRI(N3M,WFFTR)

      DO 110 L=1,N1M
      DO 110 J=1,N2M

      DO K=1,N3M  
         RFIN(K)=RDP(L,J,K)
      END DO

      CALL DF2TRF(N3M,RFIN,RFEX,WFFTR)

         FRDP_R(L,J,0)=RFEX(1)
         FRDP_I(L,J,0)=0.0
      DO M=1,N3M/2-1
         FRDP_R(L,J,M)=RFEX(2*M)
         FRDP_I(L,J,M)=RFEX(2*M+1)
      ENDDO
         FRDP_R(L,J,N3M/2)=RFEX(N3M)

110   CONTINUE

C---


      DO 70 K=1,N3M/2+1
      DO 71 I=1,N1M
      DO 71 J=1,N2M
      IJ=IGJG(I,J,1)
      PHM(IJ)=FDP_R(I,J,K-1)*(1.0-REAL(IPH0))
      REST(IJ)=FRDP_R(I,J,K-1)
71    CONTINUE
      CALL DSOLV(PHM,REST,K)

      DO 72 I=1,N1M
      DO 72 J=1,N2M
      IJ=IGJG(I,J,1)
      FDP_R(I,J,K-1)=PHM(IJ)
72    CONTINUE
70    CONTINUE

      DO 80 K=2,N3M/2
      DO 81 I=1,N1M
      DO 81 J=1,N2M
      IJ=IGJG(I,J,1)
      PHM(IJ)=FDP_I(I,J,K-1)*(1.0-REAL(IPH0))
      REST(IJ)=FRDP_I(I,J,K-1)
81    CONTINUE
      CALL DSOLV(PHM,REST,K)
      DO 82 I=1,N1M
      DO 82 J=1,N2M
      IJ=IGJG(I,J,1)
      FDP_I(I,J,K-1)=PHM(IJ)
82    CONTINUE
80    CONTINUE

C     INVERSE FFT OF FDP IN X3 DIRECTION

C---  FOR CXML
C      DO 200 I=1,N1M
C      DO 200 J=1,N2M
C
C      DO 210 K=0,N3M/2-1
C      R2(2*K)=FDP_R(I,J,K)
C210   CONTINUE
C      DO 220 K=1,N3M/2-1
C      R2(2*K+1)=-FDP_I(I,J,K)
C220   CONTINUE
C      R2(N3M)=FDP_R(I,J,N3M/2)
C      
C      STATUS=DFFT('C','R','B',R2,R2,N3M,1)
C      
C      DO 240 K=1,N3M
C      DP(I,J,K)=R2(K-1)
C240   CONTINUE
C
C200   CONTINUE
C
C---


C---  FOR IMSL
      DO 200 L=1,N1M
      DO 200 J=1,N2M

         RFIN(1)=FDP_R(L,J,0)
      DO M=1,N3M/2-1
         RFIN(2*M)=FDP_R(L,J,M)
         RFIN(2*M+1)=FDP_I(L,J,M)
      ENDDO
         RFIN(N3M)=FDP_R(L,J,N3M/2)

      CALL DF2TRB(N3M,RFIN,RFEX,WFFTR)

      SCALE = 1.0/REAL(N3M)

      DO 23 K=1,N3M
23    DP(L,J,K)=RFEX(K)*SCALE

200   CONTINUE

C---

1000  CONTINUE


      RETURN
      END


C   ********************* DSOLV ******************************

      SUBROUTINE DSOLV(PHM,REST,K)

      INCLUDE 'PARAM.H'
      PARAMETER (M1M=M1-1,M2M=M2-1)
      PARAMETER (M3M=M3-1,M3MH=M3M/2+1,M3MHM=M3M/2)
      PARAMETER (NLEV=11,NIJ=M1M*M2M+M1M*M2M/2+M1M*M2M/4
     >      +M1M*M2M/8+M1M*M2M/16+M1M*M2M/32+M1M*M2M/64+M1M*M2M/128
     >      +M1M*M2M/256+M1M*M2M/512+M1M*M2M/1024+M2M*NLEV)

      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/WAVK3/AK3(M3MH)

      COMMON/CNLEV/N1G(0:NLEV),N12G(0:NLEV)
      COMMON/DAK/N3MH,N3MD,N3MDU,N3MU,NDM
      COMMON/MLV/MLEV,IWMG,MLW,MLCMG,MHCMG,EPSM,RATIO,NSOR

      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      COMMON/POIS/IK0,IPH0

      REAL PHM(NIJ),REST(NIJ)

      NCM=0

      CALL RESMA(1,PHM,REML1,REML2,K,REST)
c      REOOL1=REML1
c      RRM=0.
c      ERRM=1.
c      WRITE(*,337) K,NCM,REML1,REML2,RRM,ERRM
c      IF(REML1.LT.EPSM) GO TO 25


c---  start iteration loop

9     CONTINUE

c     Because high wave number components converges 
c     faster than low wave number ones, 
c     maximum iteration of multigrid is adjusted for wavenumber.
c
      IF(K.LE.NSOR) THEN
         MAXCMG=MLCMG
      ELSE
         MAXCMG=MHCMG                   
      ENDIF
c
      NCM=NCM+1

c      IF (K.GT.NSOR) GO TO 100

      IF(K.GT.NSOR) THEN
c---  case of no multigrid
        CALL ITSOLV(1,PHM,REST,K)
      ELSE
c---  case of multigrid
        CALL CFILA(1,MLEV,PHM,REST,K)

        IF(IWMG.GT.0) THEN ! W-cycle
         DO 30 MW=1,IWMG
           MLWL=MLEV-MW
           IF(IWMG.EQ.1) MLWL=MLW
           CALL CLAFI(MLWL,MLEV,PHM,REST,K)
           CALL CFILA(MLWL,MLEV,PHM,REST,K)
30       CONTINUE
        ENDIF

        CALL CLAFI(1,MLEV,PHM,REST,K)

      ENDIF
c      GO TO 110
c
c100   CALL CFILA(1,1,PHM,REST,K)

110   REMOL1=REML1
c      REMOL1=REML1


c---- check residual

      CALL RESMA(1,PHM,REML1,REML2,K,REST)
      RRM=REML1/REMOL1
      ERRM=ABS(1.-RRM)

c      WRITE(*,337) K,NCM,REML1,REML2,RRM,ERRM

c---  for avoiding waste iteration, 
c     when ABS(RRM)>0.9, escape the iteration loop
      IF (RRM.GT.0.9) GO TO 12
c
      IF (NCM.GT.MAXCMG) GO TO 12
      IF (REML1.GT.EPSM) GO TO 9


c---  end of iteration loop

12    CONTINUE

      if (k.eq.1)
     &WRITE(*,337) K,NCM,REML1,REML2,RRM,ERRM      
337   FORMAT(1X,'K=',I3,1X,'NCM=',I4,1X,'RML1',E12.5,1X,
     >         'RML2',E12.5,1X,'RP',E11.4,' ERRM',E11.4)
25    CONTINUE
      RETURN
      END

C  ****************************** RESMA **********************

C     THIS ROUTINE CALCULATES THE RESIDUAL OF POISSON EQUATION WITH 
C     GIVEN PH VALUES.

c---  REML1 = maximum residual at grid level NL
c---  REML2 = rms of residual at grid level NL

      SUBROUTINE RESMA(NL,PHM,REML1,REML2,K,REST)

      INCLUDE 'PARAM.H'
      PARAMETER (M1M=M1-1,M2M=M2-1)
      PARAMETER (M3M=M3-1,M3MH=M3M/2+1,M3MHM=M3M/2)
      PARAMETER (NLEV=11,NIJ=M1M*M2M+M1M*M2M/2+M1M*M2M/4
     >      +M1M*M2M/8+M1M*M2M/16+M1M*M2M/32+M1M*M2M/64+M1M*M2M/128
     >      +M1M*M2M/256+M1M*M2M/512+M1M*M2M/1024+M2M*NLEV)

      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/WAVK3/AK3(M3MH)
      COMMON/CNLEV/N1G(0:NLEV),N12G(0:NLEV)
      COMMON/INLEV/IPLV(M1,NLEV),IMLV(M1,NLEV)
      COMMON/COMG/CO(6,NIJ)

      COMMON/POIS/IK0,IPH0

      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      REAL PHM(NIJ),REST(NIJ)

      IGJG(I,J,L)=I+(J-1)*(N1G(L)+1)+N12G(L-1)


      REML1=0.
      REML2=0.
      DO 10 IC=1,N1G(NL)
      DO 10 JC=1,N2M
      IJ=IGJG(IC,JC,NL)
      D2PH=+(CO(1,IJ)*PHM(IGJG(IPLV(IC,NL),JC,NL))+
     >       CO(2,IJ)*PHM(IGJG(IMLV(IC,NL),JC,NL))+
     >       CO(3,IJ)*PHM(IGJG(IC,JPA(JC),NL))+
     >       CO(4,IJ)*PHM(IGJG(IC,JMU(JC),NL))+
     >      (CO(5,IJ)-AK3(K))*PHM(IJ))
      RS=(-D2PH+REST(IJ))
      REML1=AMAX1(ABS(RS),REML1)  
      REML2=RS**2+REML2          
10    CONTINUE
   
      REML2=SQRT(REML2/(N1G(NL)*N2M))

      RETURN
      END

C   ********************* CFILA *******************************

C     SOLVE PHM AND CALCULATE RESIDUAL. NEXT, GO TO COARSER GRID
C     AND DO THE INTERPOLATION.

      SUBROUTINE CFILA(MLI,MLEV,PHM,REST,K)
c
c---  input
c     MLI  : finest grid level from which to begin
c     MLEV : coarsest grid level
c     PHM  : dp on finest grid
c     REST : rhs of Poisson eq. on finest grid
c---  ouput
c     PHM  : smoothed dp at level MLEV
c     REST : rhs of Poisson eq. at level MLEV
c

      INCLUDE 'PARAM.H'
      PARAMETER (M1M=M1-1,M2M=M2-1,M3M=M3-1,M3MH=M3M/2+1)
      PARAMETER (NLEV=11,NIJ=M1M*M2M+M1M*M2M/2+M1M*M2M/4
     >      +M1M*M2M/8+M1M*M2M/16+M1M*M2M/32+M1M*M2M/64+M1M*M2M/128
     >      +M1M*M2M/256+M1M*M2M/512+M1M*M2M/1024+M2M*NLEV)

      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/CNLEV/N1G(0:NLEV),N12G(0:NLEV)

      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      REAL PHM(NIJ),RES(NIJ),REST(NIJ)

      IGJG(I,J,L)=I+(J-1)*(N1G(L)+1)+N12G(L-1)

      DO 20 NL=MLI,MLEV

      CALL ITSOLV(NL,PHM,REST,K)
      CALL RESCA(NL,PHM,RES,K,REST)

      DO 19 IC=1,N1G(NL+1)
      DO 19 JC=1,N2M
      IJ=IGJG(IC,JC,NL+1)
      PHM(IJ)=0.
19    CONTINUE

      IF(NL.LT.MLEV) THEN
c---  get interpolated residual REST at level NL+1 (coarser level)
c     from the residual RES in finer mesh
      CALL TFILA(RES,NL,REST)  ! NL->NL+1
c
      ENDIF

20    CONTINUE

      RETURN
      END

C  ****************************** ITSOLV **********************

C     WITH GIVEN PREVIOUS RESULT OF PHIT, THIS ROUTINE GETS NEW PHIT
C     WITH LINE SOR ITERATION. ALSO IF IK0=1, PHIT AT I=1,J=N2M @ K_3=0 
C     IS SET TO ZERO.
C     TRIDIAGONAL SOLVER IN J

      SUBROUTINE ITSOLV(NL,PHIT,RES,K)

      INCLUDE 'PARAM.H'
      PARAMETER (M1M=M1-1,M2M=M2-1)
      PARAMETER (M3M=M3-1,M3MH=M3M/2+1,M3MHM=M3M/2)
      PARAMETER (NLEV=11,NIJ=M1M*M2M+M1M*M2M/2+M1M*M2M/4
     >      +M1M*M2M/8+M1M*M2M/16+M1M*M2M/32+M1M*M2M/64+M1M*M2M/128
     >      +M1M*M2M/256+M1M*M2M/512+M1M*M2M/1024+M2M*NLEV)

      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/WAVK3/AK3(M3MH)
      COMMON/CNLEV/N1G(0:NLEV),N12G(0:NLEV)
      COMMON/INLEV/IPLV(M1,NLEV),IMLV(M1,NLEV)
      COMMON/COMG/CO(6,NIJ)

      COMMON/POIS/IK0,IPH0
      COMMON/MAXIT/NLMAXI,NHMAXI,OMEG

      COMMON/MLV/MLEV,IWMG,MLW,MLCMG,MHCMG,EPSM,RATIO,NSOR

      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      REAL PHIT(NIJ),RES(NIJ)
      REAL AMJ(M2,M1M),ACJ(M2,M1M),APJ(M2,M1M),FJJ(M2,M1M)
      IGJG(I,J,L)=I+(J-1)*(N1G(L)+1)+N12G(L-1)


c     Because high wave number components converges 
c     faster than low wave number ones, 
c     maximum iteration of multigrid is adjusted for wavenumber.
c
      IF(K.LE.NSOR) THEN
         NMAXI=NLMAXI
      ELSE
         NMAXI=NHMAXI
      ENDIF


      DO 1000 NIT=0,NMAXI-1

      DO 10 IC=1,N1G(NL),2
      DO 10 JC=1,N2M
      IJ=IGJG(IC,JC,NL)
      APJ(JC,IC)=CO(3,IJ)
      ACJ(JC,IC)=CO(5,IJ)-AK3(K)
      AMJ(JC,IC)=CO(4,IJ)
      FJJ(JC,IC)
     >  =(-(CO(1,IJ)*PHIT(IGJG(IPLV(IC,NL),JC,NL))
     >     +CO(2,IJ)*PHIT(IGJG(IMLV(IC,NL),JC,NL)))
     >     +RES(IJ))
10    CONTINUE

c---  fixing the pressure at the location (i,j)=(1,n2m)
c     The following way showed very slow convergence as compare the way (**)
c     iko
c
c      IF(IK0.EQ.1.AND.K.EQ.1) THEN
c        AMJ(N2M,1)=0.
c        ACJ(N2M,1)=1.
c        APJ(N2M,1)=0.
c        FJJ(N2M,1)=0.
c      ENDIF

      CALL TRIBI(AMJ,ACJ,APJ,FJJ,N2M,FJJ,N1G(NL),1)

      DO 40 IC=1,N1G(NL),2
      DO 40 JC=1,N2M
      IJ=IGJG(IC,JC,NL)
      PHIT(IJ)=FJJ(JC,IC)

c---  way for fast convergence (**)
      IF (IK0.EQ.1.AND.K.EQ.1) PHIT(IJ)=FJJ(JC,IC)-FJJ(N2M,1)

40    CONTINUE

      DO 110 IC=2,N1G(NL),2
      DO 110 JC=1,N2M
      IJ=IGJG(IC,JC,NL)
      APJ(JC,IC)=CO(3,IJ)
      ACJ(JC,IC)=CO(5,IJ)-AK3(K)
      AMJ(JC,IC)=CO(4,IJ)
      FJJ(JC,IC)
     >  =(-(CO(1,IJ)*PHIT(IGJG(IPLV(IC,NL),JC,NL))
     >     +CO(2,IJ)*PHIT(IGJG(IMLV(IC,NL),JC,NL)))
     >     +RES(IJ))
110   CONTINUE

      CALL TRIBI(AMJ,ACJ,APJ,FJJ,N2M,FJJ,N1G(NL),2)

      DO 140 IC=2,N1G(NL),2
      DO 140 JC=1,N2M
      IJ=IGJG(IC,JC,NL)
      PHIT(IJ)=FJJ(JC,IC)
140   CONTINUE


1000  CONTINUE


      RETURN
      END

C  ****************************** TRIBI  **********************

      SUBROUTINE TRIBI(A,B,C,R,N,U,M,IN)
      INCLUDE 'PARAM.H'
      PARAMETER (M1M=M1-1,M2M=M2-1)

c      DIMENSION GAM(M2,M1M),A(M2,1),B(M2,1),C(M2,1),R(M2,1)
c      DIMENSION BET(M1M),U(M2,1)
      DIMENSION GAM(M2,M1M),A(M2,M1M),B(M2,M1M),C(M2,M1M),R(M2,M1M)
      DIMENSION BET(M1M),U(M2,M1M)
      
      DO 10 I=IN,M,2
      BET(I)=B(1,I)
      U(1,I)=R(1,I)/BET(I)
10    CONTINUE
      DO 11 J=2,N
      DO 21 I=IN,M,2
      GAM(J,I)=C(J-1,I)/BET(I)
      BET(I)=B(J,I)-A(J,I)*GAM(J,I)
      U(J,I)=(R(J,I)-A(J,I)*U(J-1,I))/BET(I)
21    CONTINUE
11    CONTINUE
      DO 12 J=N-1,1,-1
      DO 22 I=IN,M,2
      U(J,I)=U(J,I)-GAM(J+1,I)*U(J+1,I)
22    CONTINUE
12    CONTINUE
      RETURN
      END

C  ****************************** RESCA **********************
C     SAME PROCEDURE AS RESMA.

      SUBROUTINE RESCA(NL,PH,RS,K,RST)
c
c---  input
c     PH  : solution on level NL
c     RST : RHS of Poisson eq
c---  output
c     RS  : residual at level NL
c
      INCLUDE 'PARAM.H'
      PARAMETER (NDD=2,M1M=M1-1,M2M=M2-1)
      PARAMETER (M3M=M3-1,M3MH=M3M/2+1,M3MHM=M3M/2)
      PARAMETER (NLEV=11,NIJ=M1M*M2M+M1M*M2M/2+M1M*M2M/4
     >      +M1M*M2M/8+M1M*M2M/16+M1M*M2M/32+M1M*M2M/64+M1M*M2M/128
     >      +M1M*M2M/256+M1M*M2M/512+M1M*M2M/1024+M2M*NLEV)

      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/CNLEV/N1G(0:NLEV),N12G(0:NLEV)
      COMMON/INLEV/IPLV(M1,NLEV),IMLV(M1,NLEV)
      COMMON/COMG/CO(6,NIJ)
      COMMON/WAVK3/AK3(M3MH)
      COMMON/POIS/IK0,IPH0

      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      REAL PH(NIJ),RS(NIJ),RST(NIJ)

      IGJG(I,J,L)=I+(J-1)*(N1G(L)+1)+N12G(L-1)

      DO 10 IC=1,N1G(NL)
      DO 20 JC=1,N2M
      IJ=IGJG(IC,JC,NL)
      D2PHIJ=+(CO(1,IJ)*PH(IGJG(IPLV(IC,NL),JC,NL))+
     >         CO(2,IJ)*PH(IGJG(IMLV(IC,NL),JC,NL))+
     >         CO(3,IJ)*PH(IGJG(IC,JPA(JC),NL))+
     >         CO(4,IJ)*PH(IGJG(IC,JMU(JC),NL))+
     >        (CO(5,IJ)-AK3(K))*PH(IJ))
      RS(IJ)=(-D2PHIJ+RST(IJ))
20    CONTINUE
10    CONTINUE

c--- following statement is commented due to slow convergence, see ITSOLV
c    iko
c      IF(IK0.EQ.1) THEN
c        IF(K.EQ.1) RS(IGJG(1,N2M,1))=0.
c      ENDIF

      RETURN
      END

C  ****************************** TFILA **********************

C     RESIDUAL VALUES AT COARSER GRID BY INTERPOLATION.

      SUBROUTINE TFILA(RES,NL,REST)
c---  input       
c     RES  : residual at level NL
c---  ouput
c     REST : interpolated residual at level NL+1 (coarser level)

      INCLUDE 'PARAM.H'
      PARAMETER (NDD=2,M1M=M1-1,M2M=M2-1)
      PARAMETER (NLEV=11,NIJ=M1M*M2M+M1M*M2M/2+M1M*M2M/4
     >      +M1M*M2M/8+M1M*M2M/16+M1M*M2M/32+M1M*M2M/64+M1M*M2M/128
     >      +M1M*M2M/256+M1M*M2M/512+M1M*M2M/1024+M2M*NLEV)

      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/CNLEV/N1G(0:NLEV),N12G(0:NLEV)
      COMMON/COMG/CO(6,NIJ)

      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/MESH5/XM(M1,NLEV),DXX(0:M1,NLEV),GX(M1,NLEV)

      REAL RES(NIJ),REST(NIJ)

      IGJG(I,J,L)=I+(J-1)*(N1G(L)+1)+N12G(L-1)

      DO 10 ICL=1,N1G(NL+1)
      ICF=ICL*2-1
      DO 10 JCL=1,N2M
      JCF=JCL
      IJL=IGJG(ICL,JCL,NL+1)
      IJF=IGJG(ICF,JCF,NL)
      IPJF=IJF+1

      XC_L=0.5*(XM(ICL  ,NL+1)+XM(ICL+1,NL+1))
      X1_F=0.5*(XM(ICF  ,NL  )+XM(ICF+1,NL  ))
      X2_F=0.5*(XM(ICF+1,NL  )+XM(ICF+2,NL  ))


      REST(IJL)=(RES(IJF)*(X2_F-XC_L)+RES(IPJF)*(XC_L-X1_F))
     >         /(X2_F-X1_F)


   10 CONTINUE
      RETURN
      END

C   ********************* CLAFI *********************************

C     THE INVERSE DIRECTION WITH CFILA.

      SUBROUTINE CLAFI(MLI,MLEV,PHM,REST,K)

      INCLUDE 'PARAM.H'
      PARAMETER (M1M=M1-1,M2M=M2-1,M3M=M3-1,M3MH=M3M/2+1)
      PARAMETER (NLEV=11,NIJ=M1M*M2M+M1M*M2M/2+M1M*M2M/4
     >      +M1M*M2M/8+M1M*M2M/16+M1M*M2M/32+M1M*M2M/64+M1M*M2M/128
     >      +M1M*M2M/256+M1M*M2M/512+M1M*M2M/1024+M2M*NLEV)

      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

      REAL PHM(NIJ),REST(NIJ)

      CALL ITSOLV(MLEV,PHM,REST,K)
      DO 30 NL=MLEV,MLI+1,-1
      CALL TLAFI(PHM,NL)
      CALL ITSOLV(NL-1,PHM,REST,K)
30    CONTINUE
      RETURN
      END

C  ****************************** TLAFI **********************
C     UPDATE THE PH VALUE AT DENSER GRID BY INTERPOLATING THAT OF
C     PH AT COARSER GRID.

      SUBROUTINE TLAFI(PH,NL)
c---  input
c     NL    : level from which PH is interpolated to NL-1  
c     PH    : solution at level NL
c---  output
c     PH    : updated solution at level NL-1 (finer mesh)
c
c     Note that PH^new[NL-1]=PH^old[NL-1]+PH[NL]_{NL->NL-1}
c
      INCLUDE 'PARAM.H'
      PARAMETER (M1M=M1-1,M2M=M2-1)
      PARAMETER (NLEV=11,NIJ=M1M*M2M+M1M*M2M/2+M1M*M2M/4
     >      +M1M*M2M/8+M1M*M2M/16+M1M*M2M/32+M1M*M2M/64+M1M*M2M/128
     >      +M1M*M2M/256+M1M*M2M/512+M1M*M2M/1024+M2M*NLEV)

      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/CNLEV/N1G(0:NLEV),N12G(0:NLEV)

      COMMON/INDX/IPA(M1),IMU(M1),IMV(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/MESH5/XM(M1,NLEV),DXX(0:M1,NLEV),GX(M1,NLEV)

      REAL PH(NIJ)
      REAL TEMP(2,M2M)

      IGJG(I,J,L)=I+(J-1)*(N1G(L)+1)+N12G(L-1)

      DO 20 IC=2,N1G(NL-1)-2,2
      ICL=IC/2
      DO 25 JCL=1,N2M
      IJ=IGJG(IC,JCL,NL-1)
      IPJ=IJ+1
      IJL=IGJG(ICL,JCL,NL)
      IPJL=IJL+1


      X1_L=0.5*(XM(ICL  ,NL  )+XM(ICL+1,NL  ))
      X2_L=0.5*(XM(ICL+1,NL  )+XM(ICL+2,NL  ))
      X1_F=0.5*(XM(IC   ,NL-1)+XM(IC +1,NL-1))
      X2_F=0.5*(XM(IC +1,NL-1)+XM(IC +2,NL-1))

      TEMP(1,JCL)=PH(IJ)+
     >           ((X2_L-X1_F)*PH(IJL)+(X1_F-X1_L)*PH(IPJL))/(X2_L-X1_L)
      TEMP(2,JCL)=PH(IPJ)+
     >           ((X2_L-X2_F)*PH(IJL)+(X2_F-X1_L)*PH(IPJL))/(X2_L-X1_L)


25    CONTINUE
      DO 26 JCL=1,N2M
      IJ=IGJG(IC,JCL,NL-1)
      IPJ=IJ+1
      PH(IJ)=TEMP(1,JCL)
      PH(IPJ)=TEMP(2,JCL)
26    CONTINUE
20    CONTINUE

      IC1=1
      DO 21 JC=1,N2M
      IJ1=IGJG(IC1,JC,NL-1)
      IJL1=IGJG(IC1,JC,NL)
      TEMP(1,JC)=PH(IJL1)+PH(IJ1)
21    CONTINUE
      DO 22 JC=1,N2M
      IJ1=IGJG(IC1,JC,NL-1)
      PH(IJ1)=TEMP(1,JC)
22    CONTINUE


      IC2=N1G(NL-1)
      ICL2=IC2/2
      DO 31 JC=1,N2M
      IJ2=IGJG(IC2,JC,NL-1)
      IJL2=IGJG(ICL2,JC,NL)
      TEMP(2,JC)=PH(IJL2)+PH(IJ2)
31    CONTINUE
      DO 32 JC=1,N2M
      IJ2=IGJG(IC2,JC,NL-1)
      PH(IJ2)=TEMP(2,JC)
32    CONTINUE
      
      RETURN
      END

