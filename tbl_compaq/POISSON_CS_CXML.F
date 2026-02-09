C---------------------------------------------------------------------------
C
C     POISSON SOLVER WITH COS/FFT TRANSFORM FOR UNIFORM MESHS IN X1, X3 DIRECTIONS
C
C     POISSON EQUATION IS SOLVED BY
C     1. FFT IN X3 DIRECTION
C     2. COS TRANSFORM IN X1 DIRECTION
C     3. TRIDIAGONAL SOLVER IN X2 DIRECTION
C
C     CONTRAINT : PLANE AVERAGE OF PRESSURE OVER XZ PLANE AT J=N2M IS ZERO.
C                 X1,X3 UNIFORM MESH 
C
C     2003. 12. 17.
C     KYOUNGYOUN KIM
C     FLOW CONTROL LABORATORY
C     DEPARTMENT OF MECHANICAL ENGINEERING
C     KAIST 


C*************** INIPOISSON ******************************
      SUBROUTINE INIPOISSON
      INCLUDE 'PARAM.H'
      PARAMETER (M3M=M3-1,M3MH=M3M/2)
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/WAVK13/AK1(0:M1),AK3(0:M3MH)
      COMMON/SIZE/ALX,ALY,ALZ

      PI=ACOS(-1.0)

C     MODIFIED WAVE NUMBER DEFINITION NECESSARY FOR THE X1 DIRECTION
      
	DX1Q = (DBLE(N1M)/ALX)**2
	DO 10 L=0,N1M-1
  10  AK1(L)=2.*(1.-COS(L*PI/N1M))*DX1Q

C     MODIFIED WAVE NUMBER DEFINITION NECESSARY FOR THE X3 DIRECTION
      DO 20 M=0,N3M/2
  20  AK3(M)=2.*(1.-COS(M*2.*PI/N3M))*DX3Q

      CALL METRICPOISSON

      RETURN
      END


C  ************************ METRICPOISSON  **********************
      SUBROUTINE METRICPOISSON

      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/METPOI/PMJ(M2),PCJ(M2),PPJ(M2)

      DO 100 JC=1,N2M
      JP=JC+1
      JM=JC-1
      JUM=JC-JMU(JC)
      JUP=JPA(JC)-JC
      PMJ(JC)=JUM*(1.0/DY(JC)/H(JC))
      PPJ(JC)=JUP*(1.0/DY(JC)/H(JP))
	PCJ(JC)=-(PMJ(JC)+PPJ(JC))
      
  100 CONTINUE

      RETURN
      END

C*************** TAKEDP ********************************
      SUBROUTINE TAKEDP(DP,RDP)
      INCLUDE 'PARAM.H'
      PARAMETER (M1M=M1-1,M3M=M3-1,M3MH=M3M/2)
      PARAMETER (M2M=M2-1)
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/WAVK13/AK1(0:M1),AK3(0:M3MH)
      COMMON/METPOI/PMJ(M2),PCJ(M2),PPJ(M2)

      REAL RDP(0:M1,0:M2,0:M3),DP(0:M1,0:M2,0:M3)
      REAL CRDP(0:M1,M2,M3)
      REAL CDP(0:M1,M2,M3)
      REAL FRDP_R(0:M1,M2,0:M3MH)
      REAL FRDP_I(0:M1,M2,0:M3MH)
      REAL FDP_R(0:M1,M2,0:M3MH)
      REAL FDP_I(0:M1,M2,0:M3MH)
      REAL CMJ(0:M1,M2),CCJ(0:M1,M2),CPJ(0:M1,M2),R(0:M1,M2)

      REAL RCIN(0:M1M-1)
      REAL RCEX(0:M1M-1)
      REAL RFIN(0:M3M+1)
      REAL RFEX(0:M3M+1)


C     COSINE TRANSFORM OF RDP IN X1 DIRECTION

      DO 100 J=1,N2M
      DO 100 K=1,N3M
      DO L=0,N1M-1  
      RCIN(L)=RDP(L+1,J,K)
      END DO

      STATUS=DFCT('F',RCIN,RCEX,N1M,2,1)

      DO  L=0,N1M-1
      CRDP(L,J,K)=RCEX(L)
      END DO
100   CONTINUE


C---  IF N3M=1, WE DON'T HAVE TO FFT IN X3 DIRECTION
      IF(N3M.EQ.1) THEN

        DO J=1,N2M
        DO L=0,N1M-1
          CMJ(L,J)=PMJ(J)
          CCJ(L,J)=PCJ(J)-AK1(L)
          CPJ(L,J)=PPJ(J)
          R(L,J)=CRDP(L,J,1)
        ENDDO
        ENDDO

C     SPECIAL TREATMENT TO REMOVE SINGULARITY WHEN L

         CPJ(0,N2M-1)=0.0
         CALL TDMAI0(CMJ,CCJ,CPJ,R,R,1,N2M-1,0,0)
         CALL TDMAI0(CMJ,CCJ,CPJ,R,R,1,N2M,1,N1M-1)

        DO J=1,N2M
        DO L=0,N1M-1
          CDP(L,J,1)=R(L,J)
        ENDDO
        ENDDO

          CDP(0,N2M,1)=0.0

      GOTO 777
      ENDIF
C---

C     FFT OF RDP IN X3 DIRECTION

      DO 110 L=0,N1M-1
      DO 110 J=1,N2M
        DO K=0,N3M-1  
          RFIN(K)=CRDP(L,J,K+1)
        END DO
        STATUS=DFFT('R','C','F',RFIN,RFEX,N3M,1)   
        DO M=0,N3M/2-1
          FRDP_R(L,J,M)=RFEX(2*M)
          FRDP_I(L,J,M)=-RFEX(2*M+1)
        ENDDO
          FRDP_R(L,J,N3M/2)=RFEX(N3M)
110   CONTINUE


C     SOLUTION OF POISSON EQ. ; REAL PART

!$omp parallel do default(shared) private(CMJ,CCJ,CPJ,R)
      DO 70 M=0,N3M/2
        DO J=1,N2M
        DO L=0,N1M-1
          CMJ(L,J)=PMJ(J)
          CCJ(L,J)=PCJ(J)-AK3(M)-AK1(L)
          CPJ(L,J)=PPJ(J)
          R(L,J)=FRDP_R(L,J,M)
        ENDDO
        ENDDO

C     SPECIAL TREATMENT TO REMOVE SINGULARITY WHEN L=K=0

      IF (M.EQ.0) THEN
         CPJ(0,N2M-1)=0.0
         CALL TDMAI0(CMJ,CCJ,CPJ,R,R,1,N2M-1,0,0)
         CALL TDMAI0(CMJ,CCJ,CPJ,R,R,1,N2M,1,N1M-1)
      ELSE
         CALL TDMAI0(CMJ,CCJ,CPJ,R,R,1,N2M,0,N1M-1)
      ENDIF

        DO J=1,N2M
        DO L=0,N1M-1
          FDP_R(L,J,M)=R(L,J)
        ENDDO
        ENDDO
 70   CONTINUE
          FDP_R(0,N2M,0)=0.0

C     SOLUTION OF POISSON EQ. ; IMAGE PART

!$omp parallel do default(shared) private(CMJ,CCJ,CPJ,R)
      DO 80 M=1,N3M/2-1
        DO J=1,N2M
        DO L=0,N1M-1
          CMJ(L,J)=PMJ(J)
          CCJ(L,J)=PCJ(J)-AK3(M)-AK1(L)
          CPJ(L,J)=PPJ(J)
          R(L,J)=FRDP_I(L,J,M)
        ENDDO
        ENDDO

        CALL TDMAI0(CMJ,CCJ,CPJ,R,R,1,N2M,0,N1M-1)

        DO J=1,N2M
        DO L=0,N1M-1
          FDP_I(L,J,M)=R(L,J)
        ENDDO
        ENDDO
 80   CONTINUE


C     INVERSE FFT OF FDP IN X3 DIRECTION

      DO 200 L=0,N1M-1
      DO 200 J=1,N2M

      DO M=0,N3M/2-1
         RFIN(2*M)=FDP_R(L,J,M)
      ENDDO
      DO M=1,N3M/2-1
         RFIN(2*M+1)=-FDP_I(L,J,M)
      ENDDO
         RFIN(N3M)=FDP_R(L,J,N3M/2)

      STATUS=DFFT('C','R','B',RFIN,RFEX,N3M,1)

      DO 23 K=1,N3M
23    CDP(L,J,K)=RFEX(K-1)

 200  CONTINUE


C      INVERSE COSINE TRANSFORM OF CDP IN X1 DIRECTION
777   CONTINUE

      DO 40 J=1,N2M
      DO 40 K=1,N3M

      DO 41 L=0,N1M-1
41    RCIN(L)=CDP(L,J,K)

      STATUS=DFCT('B',RCIN,RCEX,N1M,2,1)

      DO 42 I=1,N1M
42    DP(I,J,K)=RCEX(I-1)

40    CONTINUE

      RETURN
      END

C*************** TDMAI0 **********************
      SUBROUTINE TDMAI0(A,B,C,R,X,NS,NF,NIS,NIF)
      include 'PARAM.H'
      REAL A(0:M1,M2),B(0:M1,M2),C(0:M1,M2),R(0:M1,M2),X(0:M1,M2)
      REAL BET,GAM(0:M1,0:M2)

c
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
c
c
      RETURN
      END
