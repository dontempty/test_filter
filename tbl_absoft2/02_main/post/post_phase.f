      program post_phase
      include '../PARAM.H'


      call setup
      call set_moni 

c      call t_statistics
c      call o_statistics

      call output_along_x
      call output_moni_x
c      call output_2d

      call output_phase_2d
c      call output_phase_moni_x

      stop
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

      OPEN(1,FILE='../tbl.set',STATUS='OLD')
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
c      CALL INDICES 
c      CALL INIPOISSON

      RETURN
      END


C***************** MESH ***********************     
      SUBROUTINE MESH
      include '../PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/SIZE/ALX,ALY,ALZ
      COMMON/FILENAME/FILEINI,FILEOUT,FILEAVG
      CHARACTER*15 FILEINI,FILEOUT,FILEAVG
 
      PI=ACOS(-1.0)

C     CREATE THE GRID IN X1 DIRECTION      
C         X(0)=0.0 
C      DO I=1,N1 
C         X(I)=REAL(I-1)*ALX/REAL(N1M)
C         X(I)=0.5*ALX*(1.0-COS(PI*REAL(I-1)/REAL(N1M)))
C         X(I)=ALX*(1.0-COS(0.5*PI*REAL(I-1)/REAL(N1M)))
C      ENDDO

         X(0)=0.0 
      OPEN(10,FILE='../GRID_X.PLT',STATUS='UNKNOWN')
      DO I=1,N1
      READ(10,*) DUMMY,X(I),DUMMY
      ENDDO
      CLOSE(10)

         DX(0)=0.0
      DO I=1,N1M
         DX(I)=X(I+1)-X(I)
      ENDDO 
         DX(N1)=0.0

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
c      CS=0.035
c      A=TANH(CS*REAL(N2M))
c
c         Y(0)=0.0
c      DO J=1,N2
c         Y(J)=ALY*(A-TANH(CS*REAL(N2-J)))/A
c      ENDDO

         Y(0)=0.0
      OPEN(10,FILE='../GRID.PLT',STATUS='UNKNOWN')
      DO J=1,N2
      READ(10,*) DUMMY,Y(J),DUMMY
      ENDDO
      CLOSE(10)

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


c      OPEN(10,FILE='GRID.PLT',STATUS='UNKNOWN')
c      DO J=1,N2
c      WRITE(10,*) J,Y(J),DY(J)
c      ENDDO
c      CLOSE(10)
c
c      OPEN(10,FILE='GRID_X.PLT',STATUS='UNKNOWN')
c      DO I=1,N1
c      WRITE(10,*) I,X(I),DX(I)
c      ENDDO
c      CLOSE(10)


      RETURN
      END

c*************** set_moni ***************************************
c     find i-index corresponding on given monitoring x location
      subroutine set_moni 
      include '../PARAM.H'
      COMMON/PARA/RE
      common/dim/n1,n2,n3,n1m,n2m,n3m
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)

      common/moni/moni_i(100),x_moni(100),num_moni

      open(10,file='moni_x.set',status='old')
       read (10,*) num_moni
      do k=1,num_moni
       read (10,*) x_moni(k)
      enddo
      close(10)

                   x_slot=75.0
                   u_tau_in=0.0492

      do 1 k=1,num_moni

        x_moni(k)=x_moni(k)/RE/u_tau_in+x_slot

        do i=1,n1m-1
         temp=(x_moni(k)-0.5*(x(i  )+x(i+1)))
     &       *(x_moni(k)-0.5*(x(i+1)+x(i+2)))
      
         if(temp.le.0.0) then 
           moni_i(k)=i
           goto 1 
         endif
        enddo

1     continue

      do k=1,num_moni
       write(*,*) k,x_moni(k),moni_i(k)
     &       ,0.5*(x(moni_i(k))+x(moni_i(k)+1))
     &       ,(0.5*(x(moni_i(k))+x(moni_i(k)+1))-x_slot)*u_tau_in*RE
      enddo

      return
      end


c*************** t_statistics *********************************
c     t_** : time-averaged quantities

      subroutine t_statistics
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/LES_OPT/ILES

      real t_u1(3,M1,M2),t_p1(M1,M2)        
      real t_u2(6,M1,M2),t_p2(M1,M2)
      real t_u3(3,M1,M2),t_u4(3,M1,M2)
      real t_vor1(3,M1,M2),t_vor2(3,M1,M2)    
      real t_up2(3,M1,M2),t_up3(3,M1,M2),t_up4(3,M1,M2)

      real VM(3,M1,M2),PM(M1,M2)
      real VVM(6,M1,M2),PPM(M1,M2)
      real V3M(3,M1,M2),V4M(3,M1,M2)
      real DISSP(6,M1,M2)
      real TRANS(6,2,M1,M2)
      real PHI(6,M1,M2)
      real VORM(3,M1,M2),VORQM(3,M1,M2)    
      real REDIS(6,M1,M2)

      real t_tau(M1,M2)
      real t_ev1(M1,M2),t_ev2(M1,M2)
      real t_esgs(M1,M2),t_emol(M1,M2)
      real t_f_nega_ev(M1,M2)

      real TAU(M1,M2),EV_ZT(M1,M2),EV2_ZT(M1,M2)
     >    ,E_SGS(M1,M2),E_MOL(M1,M2)           ! AKSELVOLL & MOIN TF-63 PP.71
     >    ,F_NEGA_EV(M1,M2),F_NEGA_TV(M1,M2)   ! AKSELVOLL & MOIN TF-63 PP.70

      t_u1(:,:,:)  =0.0         ! u
      t_u2(:,:,:)  =0.0         ! u"u"
      t_u3(:,:,:)  =0.0         ! u"^3
      t_u4(:,:,:)  =0.0         ! u"^4
      t_vor1(:,:,:)=0.0         ! omega
      t_vor2(:,:,:)=0.0         ! omega"^2
      t_p1(:,:)    =0.0         ! p
      t_p2(:,:)    =0.0         ! p"^2

      t_tau(:,:)      =0.0      ! tau_12
      t_ev1(:,:)      =0.0      ! nu_t
      t_ev2(:,:)      =0.0      ! nu_t"^2
      t_esgs(:,:)     =0.0      ! SGS dissipation
      t_emol(:,:)     =0.0      ! molecular dissipation
      t_f_nega_ev(:,:)=0.0      ! fraction of negative nu_t  

      do 1 iphase=0,mt

        call read_phaseavg(iphase
     &           ,vm,pm,vvm,ppm,v3m,v4m  
     &           ,dissp,trans,phi,vorm,vorqm,redis)
       
        do j=1,n2m
         do i=1,n1m
          do nv=1,3 

           t_u1(nv,i,j)=(t_u1(nv,i,j)*real(iphase)+ vm(nv,i,j))
     &                  /real(iphase+1)

           t_u2(nv,i,j)=(t_u2(nv,i,j)*real(iphase)
     &                     +vvm(nv,i,j)-vm(nv,i,j)**2)
     &                  /real(iphase+1)

           t_u3(nv,i,j)=(t_u3(nv,i,j)*real(iphase)
     &                     +v3m(nv,i,j)+2.0*vm(nv,i,j)**3
     &                     -3.0*vm(nv,i,j)*vvm(nv,i,j))
     &                  /real(iphase+1)

           t_u4(nv,i,j)=(t_u4(nv,i,j)*real(iphase)
     &                     +v4m(nv,i,j)-3.0*vm(nv,i,j)**4
     &                     +6.0*vm(nv,i,j)**2*vvm(nv,i,j) 
     &                     -4.0*vm(nv,i,j)*v3m(nv,i,j))
     &                  /real(iphase+1)


           t_vor1(nv,i,j)=(t_vor1(nv,i,j)*real(iphase)+ vorm(nv,i,j))
     &                  /real(iphase+1)

           t_vor2(nv,i,j)=(t_vor2(nv,i,j)*real(iphase)
     &                     +vorqm(nv,i,j)-vorm(nv,i,j)**2)
     &                  /real(iphase+1)

          enddo  !nv

           t_p1(i,j)=(t_p1(i,j)*real(iphase)+ pm(i,j))
     &                  /real(iphase+1)

           t_p2(i,j)=(t_p2(i,j)*real(iphase)
     &                   + ppm(i,j)-pm(i,j)**2)
     &                  /real(iphase+1)


           t_u2(4,i,j)=(t_u2(4,i,j)*real(iphase)
     &                     +vvm(4,i,j)-vm(1,i,j)*vm(2,i,j))
     &                  /real(iphase+1)

           t_u2(5,i,j)=(t_u2(5,i,j)*real(iphase)
     &                     +vvm(5,i,j)-vm(1,i,j)*vm(3,i,j))
     &                  /real(iphase+1)

           t_u2(6,i,j)=(t_u2(6,i,j)*real(iphase)
     &                     +vvm(6,i,j)-vm(2,i,j)*vm(3,i,j))
     &                  /real(iphase+1)


         enddo
        enddo

1     continue


      t_up2(:,:,:)  =0.0         ! u'u'
      t_up3(:,:,:)  =0.0         ! u'^3
      t_up4(:,:,:)  =0.0         ! u'^4

      do 3 iphase=0,mt

        call read_phaseavg(iphase
     &           ,vm,pm,vvm,ppm,v3m,v4m  
     &           ,dissp,trans,phi,vorm,vorqm,redis)

        do j=1,n2m
         do i=1,n1m
          do nv=1,3 

           t_up2(nv,i,j)=(t_up2(nv,i,j)*real(iphase)
     &                   +vvm(nv,i,j)-vm(nv,i,j)*t_u1(nv,i,j))
     &                   /real(iphase+1)
      
           t_up3(nv,i,j)=(t_up3(nv,i,j)*real(iphase)
     &                   +v3m(nv,i,j)-3.*vvm(nv,i,j)*t_u1(nv,i,j)
     &                   +3.*vm(nv,i,j)*t_u1(nv,i,j)**2-t_u1(nv,i,j)**3)
     &                   /real(iphase+1)

           t_up4(nv,i,j)=(t_up4(nv,i,j)*real(iphase)
     &                   +v4m(nv,i,j)-4.*v3m(nv,i,j)*t_u1(nv,i,j)
     &                   +6.*vvm(nv,i,j)*t_u1(nv,i,j)**2
     &                   -4.*vm(nv,i,j)*t_u1(nv,i,j)**3
     &                   +t_u1(nv,i,j)**4)
     &                   /real(iphase+1)

        enddo
       enddo  
      enddo  

3     continue


      open(10,file='t_staticstics.dat',status='unknown'
     &       ,form='unformatted' 
     &       ,action='write')

      write(10) 
     &  (((t_u1(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u
     & ,(((t_u2(nv,i,j)  ,nv=1,6),i=1,n1m),j=1,n2m)           ! u"u"
     & ,(((t_u3(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^3
     & ,(((t_u4(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^4
     & ,(((t_vor1(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega
     & ,(((t_vor2(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega"^2
     & ,((t_p1(i,j)              ,i=1,n1m),j=1,n2m)           ! p
     & ,((t_p2(i,j)              ,i=1,n1m),j=1,n2m)           ! p"^2
     & ,(((t_up2(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'u' ; u'=u-U, cf. u"=u-<u>
     & ,(((t_up3(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'^3
     & ,(((t_up4(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'^4

      close(10)


c      open(10,file='t_test.plt',status='unknown')
c
c      write(10,*) 'variables=x,y,u,v,p,omega3'
c      write(10,*) 'zone i=',n1m,',j=',n2m,',f=point'
c
c      do j=1,n2m
c       do i=1,n1m
c
c         xc=0.5*(x(i+1)+x(i))
c         yc=0.5*(y(j+1)+y(j))
c
c         write(10,100) xc,yc
c     &                ,(t_u1(nv,i,j),nv=1,2)
c     &                , t_p1(i,j)
c     &                , t_vor1(3,i,j)
c100      format(6(e12.5,x))
c
c       enddo
c      enddo
c
c      close(10)



      if(ILES.eq.1) then

       do 2 iphase=0,mt
        call read_phaseavg_les(iphase,
     &               tau,ev_zt,ev2_zt,
     &               e_sgs,e_mol,f_nega_ev,f_nega_tv)

        do j=1,n2m
         do i=1,n1m

           t_tau(i,j)=(t_tau(i,j)*real(iphase)+ tau(i,j))
     &                  /real(iphase+1)

           t_ev1(i,j)=(t_ev1(i,j)*real(iphase)+ ev_zt(i,j))
     &                  /real(iphase+1)

           t_ev2(i,j)=(t_ev2(i,j)*real(iphase)
     &                   + ev2_zt(i,j)-ev_zt(i,j)**2)
     &                  /real(iphase+1)

           t_esgs(i,j)=(t_esgs(i,j)*real(iphase)+ e_sgs(i,j))
     &                  /real(iphase+1)
           t_emol(i,j)=(t_emol(i,j)*real(iphase)+ e_mol(i,j))
     &                  /real(iphase+1)
           t_f_nega_ev(i,j)
     &           =(t_f_nega_ev(i,j)*real(iphase)+ f_nega_tv(i,j))
     &                  /real(iphase+1)

         enddo
        enddo

2      continue

      open(10,file='t_staticstics_les.dat',status='unknown'
     &       ,form='unformatted' 
     &       ,action='write')
      write(10)
     &  ((t_tau(i,j),i=1,n1m),j=1,n2m)            ! tau_12
     & ,((t_ev1(i,j),i=1,n1m),j=1,n2m)            ! nu_t
     & ,((t_ev2(i,j),i=1,n1m),j=1,n2m)            ! nu_t"^2
     & ,((t_esgs(i,j),i=1,n1m),j=1,n2m)           ! SGS dissipation
     & ,((t_emol(i,j),i=1,n1m),j=1,n2m)           ! molecular dissipation
     & ,((t_f_nega_ev(i,j),i=1,n1m),j=1,n2m)      ! fraction of negative nu_t  

      close(10)

      endif !ILES=1



      return
      end

c*************** o_statistics *********************************
c     t_** : time-averaged quantities
c     o_** : rms of oscillating quantities

      subroutine o_statistics
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)

      real t_u1(3,M1,M2),t_p1(M1,M2)        
      real t_u2(6,M1,M2),t_p2(M1,M2)
      real t_u3(3,M1,M2),t_u4(3,M1,M2)
      real t_vor1(3,M1,M2),t_vor2(3,M1,M2)    
      real t_up2(3,M1,M2),t_up3(3,M1,M2),t_up4(3,M1,M2)

      real o_u1(3,M1,M2),o_p1(M1,M2)
      real o_u2(6,M1,M2),o_p2(M1,M2)
      real o_u3(3,M1,M2),o_u4(3,M1,M2)
      real o_vor1(3,M1,M2),o_vor2(3,M1,M2)    

      real VM(3,M1,M2),PM(M1,M2)
      real VVM(6,M1,M2),PPM(M1,M2)
      real V3M(3,M1,M2),V4M(3,M1,M2)
      real DISSP(6,M1,M2)
      real TRANS(6,2,M1,M2)
      real PHI(6,M1,M2)
      real VORM(3,M1,M2),VORQM(3,M1,M2)    
      real REDIS(6,M1,M2)

      o_u1(:,:,:)  =0.0         ! tilde u
      o_u2(:,:,:)  =0.0         ! tilde u"u"
      o_u3(:,:,:)  =0.0         ! tilde u"^3
      o_u4(:,:,:)  =0.0         ! tilde u"^4
      o_vor1(:,:,:)=0.0         ! tilde omega
      o_vor2(:,:,:)=0.0         ! tilde omega"^2
      o_p1(:,:)    =0.0         ! tilde p
      o_p2(:,:)    =0.0         ! tilde p"^2


      open(10,file='t_staticstics.dat',status='unknown'
     &       ,form='unformatted',action='read') 

      read(10) 
     &  (((t_u1(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u
     & ,(((t_u2(nv,i,j)  ,nv=1,6),i=1,n1m),j=1,n2m)           ! u"u"
     & ,(((t_u3(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^3
     & ,(((t_u4(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^4
     & ,(((t_vor1(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega
     & ,(((t_vor2(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega"^2
     & ,((t_p1(i,j)              ,i=1,n1m),j=1,n2m)           ! p
     & ,((t_p2(i,j)              ,i=1,n1m),j=1,n2m)           ! p"^2
     & ,(((t_up2(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'u' ; u'=u-U, cf. u"=u-<u>
     & ,(((t_up3(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'^3
     & ,(((t_up4(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'^4

      close(10)


      do 2 iphase=0,mt

        call read_phaseavg(iphase
     &           ,vm,pm,vvm,ppm,v3m,v4m  
     &           ,dissp,trans,phi,vorm,vorqm,redis)
       

        do j=1,n2m
         do i=1,n1m
          do nv=1,3 

           o_u1(nv,i,j)=(o_u1(nv,i,j)*real(iphase)
     &                   +(vm(nv,i,j)
     &                    -t_u1(nv,i,j))**2)
     &                  /real(iphase+1)

           o_u2(nv,i,j)=(o_u2(nv,i,j)*real(iphase)
     &                     +(vvm(nv,i,j)-vm(nv,i,j)**2
     &                      -t_u2(nv,i,j))**2)
     &                  /real(iphase+1)

           o_u3(nv,i,j)=(o_u3(nv,i,j)*real(iphase)
     &                     +(v3m(nv,i,j)+2.0*vm(nv,i,j)**3
     &                      -3.0*vm(nv,i,j)*vvm(nv,i,j)
     &                      -t_u3(nv,i,j))**2)
     &                  /real(iphase+1)

           o_u4(nv,i,j)=(o_u4(nv,i,j)*real(iphase)
     &                     +(v4m(nv,i,j)-3.0*vm(nv,i,j)**4
     &                      +6.0*vm(nv,i,j)**2*vvm(nv,i,j) 
     &                      -4.0*vm(nv,i,j)*v3m(nv,i,j)
     &                      -t_u4(nv,i,j))**2)
     &                  /real(iphase+1)


           o_vor1(nv,i,j)=(o_vor1(nv,i,j)*real(iphase)
     &                      +(vorm(nv,i,j)
     &                       -t_vor1(nv,i,j))**2)
     &                  /real(iphase+1)

           o_vor2(nv,i,j)=(o_vor2(nv,i,j)*real(iphase)
     &                     +(vorqm(nv,i,j)-vorm(nv,i,j)**2
     &                      -t_vor2(nv,i,j))**2)
     &                  /real(iphase+1)

          enddo  !nv

           o_p1(i,j)=(o_p1(i,j)*real(iphase)
     &                    +(pm(i,j)
     &                     -t_p1(i,j))**2)
     &                  /real(iphase+1)

           o_p2(i,j)=(o_p2(i,j)*real(iphase)
     &                   + (ppm(i,j)-pm(i,j)**2
     &                     -t_p2(i,j))**2)
     &                  /real(iphase+1)


           o_u2(4,i,j)=(o_u2(4,i,j)*real(iphase)
     &                     +(vvm(4,i,j)-vm(1,i,j)*vm(2,i,j) 
     &                      -t_u2(4,i,j))**2)
     &                  /real(iphase+1)

           o_u2(5,i,j)=(o_u2(5,i,j)*real(iphase)
     &                     +(vvm(5,i,j)-vm(1,i,j)*vm(3,i,j)
     &                      -t_u2(5,i,j))**2)
     &                  /real(iphase+1)

           o_u2(6,i,j)=(o_u2(6,i,j)*real(iphase)
     &                     +(vvm(6,i,j)-vm(2,i,j)*vm(3,i,j)
     &                      -t_u2(6,i,j))**2)
     &                  /real(iphase+1)


         enddo
        enddo

c      if(iphase.eq.0) open(80,file='moni.plt',status='unknown')
c      write(80,200) iphase
c     &        ,vorm(3,175,10),t_vor1(3,175,10),o_vor1(3,175,10) 
c200   format(i5,4(x,e12.5))
c      if(iphase.eq.mt) close(80)

2     continue


      open(10,file='o_staticstics.dat',status='unknown'
     &       ,form='unformatted' 
     &       ,action='write')

      write(10) 
     &  (((o_u1(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u
     & ,(((o_u2(nv,i,j)  ,nv=1,6),i=1,n1m),j=1,n2m)           ! u"u"
     & ,(((o_u3(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^3
     & ,(((o_u4(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^4
     & ,(((o_vor1(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega
     & ,(((o_vor2(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega"^2
     & ,((o_p1(i,j)              ,i=1,n1m),j=1,n2m)           ! p
     & ,((o_p2(i,j)              ,i=1,n1m),j=1,n2m)           ! p"^2

      close(10)

c      open(10,file='o_test.plt',status='unknown')
c
c      write(10,*) 'variables=x,y,u,v,p,omega3'
c      write(10,*) 'zone i=',n1m,',j=',n2m,',f=point'
c
c      do j=1,n2m
c       do i=1,n1m
c
c         xc=0.5*(x(i+1)+x(i))
c         yc=0.5*(y(j+1)+y(j))
c
c         write(10,100) xc,yc
c     &                ,(o_u1(nv,i,j)**0.5,nv=1,2)
c     &                , o_p1(i,j)**0.5
c     &                , o_vor1(3,i,j)**0.5
c100      format(6(e12.5,x))
c
c       enddo
c      enddo
c
c      close(10)

      return
      end

C*************** read_phaseavg  **********************
      SUBROUTINE read_phaseavg(iphase
     &           ,vm,pm,vvm,ppm,v3m,v4m  
     &           ,dissp,trans,phi,vorm,vorqm,redis)

      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m

      real VM(3,M1,M2),PM(M1,M2)
      real VVM(6,M1,M2),PPM(M1,M2)
      real V3M(3,M1,M2),V4M(3,M1,M2)
      real DISSP(6,M1,M2)
      real TRANS(6,2,M1,M2)
      real PHI(6,M1,M2)
      real VORM(3,M1,M2),VORQM(3,M1,M2)    
      real REDIS(6,M1,M2)

      character*30 phasefile 

      phasefile='phaseavg.'
      nn=index(phasefile,'.')
      write(unit=phasefile(nn+1:),fmt='(bn,i4.4)') iphase
      phasefile='../'//phasefile

      open(11,file=phasefile,status='unknown',form='unformatted'
     &       ,action='read')
      
          read(11) nphavg 
      do i=1,n1m
        do j=1,n2m
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
        enddo
      enddo

      close(11)

      write(*,*) 'iphase=',iphase,'nphavg=',nphavg
c      write(*,*) vm(1,100,20),vvm(1,100,20),v3m(1,100,20),v4m(1,100,20)


c      if(iphase.eq.0) open(80,file='moni.plt',status='unknown')
c      write(80,100) iphase,(vm(nv,175,1),nv=1,3),pm(175,1) 
c100   format(i5,4(x,e12.5))
c      if(iphase.eq.mt) close(80)


      return
      end


C*************** read_phaseavg_les **********************
      SUBROUTINE read_phaseavg_les(iphase,
     &                  tau,ev_zt,ev2_zt,
     &                  e_sgs,e_mol,f_nega_ev,f_nega_tv)

      INCLUDE '../PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M

      real TAU(M1,M2),EV_ZT(M1,M2),EV2_ZT(M1,M2)
     >    ,E_SGS(M1,M2),E_MOL(M1,M2)           ! AKSELVOLL & MOIN TF-63 PP.71
     >    ,F_NEGA_EV(M1,M2),F_NEGA_TV(M1,M2)   ! AKSELVOLL & MOIN TF-63 PP.70

      character*30 phasefile


      phasefile='phaseles.'
      nn=index(phasefile,'.')
      write(unit=phasefile(nn+1:),fmt='(bn,i4.4)') iphase
      phasefile='../'//phasefile

      OPEN(11,FILE=phasefile,STATUS='UNKNOWN',form='unformatted'
     &       ,action='read')

         read(11) nphavg_les
      DO I=1,N1M
        do J=1,N2M
         read(11) TAU(I,J)
         read(11) EV_ZT(I,J),EV2_ZT(I,J)
         read(11) E_SGS(I,J),E_MOL(I,J)
         read(11) F_NEGA_EV(I,J),F_NEGA_TV(I,J)
        enddo
      enddo

      CLOSE(11)

      write(*,*) 'iphase=',iphase,'nphavg_les=',nphavg_les

      RETURN
      END

c*************** output_moni_x **************************
      subroutine output_moni_x

      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/PARA/RE
      common/wall/utau(m1),utau_o(m1)
      common/moni/moni_i(100),x_moni(100),num_moni
      COMMON/LES_OPT/ILES

      real t_u1(3,M1,M2),t_p1(M1,M2)        
      real t_u2(6,M1,M2),t_p2(M1,M2)
      real t_u3(3,M1,M2),t_u4(3,M1,M2)
      real t_vor1(3,M1,M2),t_vor2(3,M1,M2)    
      real t_up2(3,M1,M2),t_up3(3,M1,M2),t_up4(3,M1,M2)

      real o_u1(3,M1,M2),o_p1(M1,M2)
      real o_u2(6,M1,M2),o_p2(M1,M2)
      real o_u3(3,M1,M2),o_u4(3,M1,M2)
      real o_vor1(3,M1,M2),o_vor2(3,M1,M2)    

      real t_tau(M1,M2)
      real t_ev1(M1,M2),t_ev2(M1,M2)
      real t_esgs(M1,M2),t_emol(M1,M2)
      real t_f_nega_ev(M1,M2)
 
      character*30 filename

      open(10,file='t_staticstics.dat',status='old'
     &       ,form='unformatted' 
     &       ,action='read')

      read(10) 
     &  (((t_u1(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u
     & ,(((t_u2(nv,i,j)  ,nv=1,6),i=1,n1m),j=1,n2m)           ! u"u"  
     & ,(((t_u3(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^3
     & ,(((t_u4(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^4
     & ,(((t_vor1(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega
     & ,(((t_vor2(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega"^2
     & ,((t_p1(i,j)              ,i=1,n1m),j=1,n2m)           ! p
     & ,((t_p2(i,j)              ,i=1,n1m),j=1,n2m)           ! p"^2
     & ,(((t_up2(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'u' ; u'=u-U, cf. u"=u-<u>
     & ,(((t_up3(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'^3
     & ,(((t_up4(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'^4

      close(10)

      open(10,file='o_staticstics.dat',status='old'
     &       ,form='unformatted' 
     &       ,action='read')

      read(10) 
     &  (((o_u1(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u
     & ,(((o_u2(nv,i,j)  ,nv=1,6),i=1,n1m),j=1,n2m)           ! u"u"
     & ,(((o_u3(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^3
     & ,(((o_u4(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^4
     & ,(((o_vor1(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega
     & ,(((o_vor2(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega"^2
     & ,((o_p1(i,j)              ,i=1,n1m),j=1,n2m)           ! p
     & ,((o_p2(i,j)              ,i=1,n1m),j=1,n2m)           ! p"^2

      close(10)


      t_tau(:,:)      =0.0      ! tau_12
      t_ev1(:,:)      =0.0      ! nu_t
      t_ev2(:,:)      =0.0      ! nu_t"^2
      t_esgs(:,:)     =0.0      ! SGS dissipation
      t_emol(:,:)     =0.0      ! molecular dissipation
      t_f_nega_ev(:,:)=0.0      ! fraction of negative nu_t  

      if(ILES.eq.1) then
      open(10,file='t_staticstics_les.dat',status='old'
     &       ,form='unformatted' 
     &       ,action='read')
      read(10)
     &  ((t_tau(i,j),i=1,n1m),j=1,n2m)            ! tau_12
     & ,((t_ev1(i,j),i=1,n1m),j=1,n2m)            ! nu_t
     & ,((t_ev2(i,j),i=1,n1m),j=1,n2m)            ! nu_t"^2
     & ,((t_esgs(i,j),i=1,n1m),j=1,n2m)           ! SGS dissipation
     & ,((t_emol(i,j),i=1,n1m),j=1,n2m)           ! molecular dissipation
     & ,((t_f_nega_ev(i,j),i=1,n1m),j=1,n2m)      ! fraction of negative nu_t  

      close(10)
      endif

      do k=1,num_moni
         i=moni_i(k) 

c-
      filename='u_x'
      nn=index(filename,'x')
      write(unit=filename(nn+1:),fmt='(bn,i2.2,a4)') k,'.dat'
      open(10,file=filename,status='unknown')
       write(10,102)
102    format('variables=y/`q_i_n,y^+_o,y^+local,u,u^+,u^+local')
       do j=1,n2m
          yc=0.5*(y(j+1)+y(j))
          write(10,200) yc 
     &                 ,yc*utau_o(i)*re
     &                 ,yc*utau(i)*re
     &                 ,t_u1(1,i,j)
     &                 ,t_u1(1,i,j)/utau_o(i)
     &                 ,t_u1(1,i,j)/utau(i)
200       format(6(e12.5,x))
       enddo
      close(10)
c-
c      filename='new_u_x'
c      nn=index(filename,'x')
c      write(unit=filename(nn+1:),fmt='(bn,i2.2,a4)') k,'.dat'
c      open(10,file=filename,status='unknown')
c       write(10,402)
c402    format('variables=y/`q_i_n,y^+_o,y^+local,u,u^+,u^+local,dudy')
c       do j=1,n2m
c          yc=0.5*(y(j+1)+y(j))
c          write(10,400) yc 
c     &                 ,yc*utau_o(i)*re
c     &                 ,yc*utau(i)*re
c     &                 ,t_u1(1,i,j)
c     &                 ,t_u1(1,i,j)/utau_o(i)
c     &                 ,t_u1(1,i,j)/utau(i)
c     &                 ,(t_u1(1,i,j+1)-t_u1(1,i,j))/H(j)
400       format(7(e12.5,x))
c       enddo
c      close(10)

c-
      filename='u_tilde_x'
      nn=index(filename,'x')
      write(unit=filename(nn+1:),fmt='(bn,i2.2,a4)') k,'.dat'
      open(10,file=filename,status='unknown')
       do j=1,n2m
          yc=0.5*(y(j+1)+y(j))
          write(10,203) yc*utau_o(i)*re
     &                 ,o_u1(1,i,j)**0.5/utau_o(i)
     &                 ,o_u1(2,i,j)**0.5/utau_o(i)
203       format(3(e12.5,x))
       enddo
      close(10)

c-
      filename='uf_rms_x'
      nn=index(filename,'x')
      write(unit=filename(nn+1:),fmt='(bn,i2.2,a4)') k,'.dat'
      open(10,file=filename,status='unknown')
       write(10,101)
101    format('variables=y^+_o,u,v,w,uv,uv2')
       do j=1,n2m
          yc=0.5*(y(j+1)+y(j))
          write(10,202) yc*utau_o(i)*re
     &                ,(t_u2(nv,i,j)**0.5/utau_o(i),nv=1,3)
     &                ,-(t_u2(4,i,j)+t_tau(i,j))/utau_o(i)**2
     &                ,-t_u2(4,i,j)/utau_o(i)**2
202       format(5(e12.5,x))
       enddo
      close(10)



cc-
      filename='k_sta_x'
      nn=index(filename,'x')
      write(unit=filename(nn+1:),fmt='(bn,i2.2,a4)') k,'.dat'
      open(10,file=filename,status='unknown')
       do j=1,n2m
          yc=0.5*(y(j+1)+y(j))
          write(10,203) yc,yc*utau_o(i)*re
     &             ,2.0*t_u2(1,i,j)/( t_u2(2,i,j)+t_u2(3,i,j) )
       enddo
      close(10)


c-vorticity fluctuations : omega
      filename='vor_t_x'
      nn=index(filename,'x')
      write(unit=filename(nn+1:),fmt='(bn,i2.2,a4)') k,'.dat'
      open(10,file=filename,status='unknown')
       write(10,301)
301    format('variables=y,y^+_o,`w_x,`w_y,`w_z')
       do j=1,n2m
          yc=0.5*(y(j+1)+y(j))
          write(10,201) yc,yc*utau_o(i)*re
     &                 ,(t_vor2(nv,i,j)**0.5,nv=1,3)
201       format(5(e12.5,x))
       enddo
      close(10)


c      if(ILES.eq.1) then
cc-eddy viscosity : nu_t/nu
c      filename='nu_t_x'
c      nn=index(filename,'x')
c      write(unit=filename(nn+1:),fmt='(bn,i2.2,a4)') k,'.dat'
c      open(10,file=filename,status='unknown')
c       do j=1,n2m
c          yc=0.5*(y(j+1)+y(j))
c          write(10,203) yc*utau_o(i)*re
c     &                 ,t_ev1(i,j)*re
c     &                 ,t_tau(i,j)
c       enddo
c      close(10)
c      endif

      enddo

      return
      end

c*************** output_along_x **************************
      subroutine output_along_x

      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/PARA/RE
      common/wall/utau(m1),utau_o(m1)
      COMMON/LES_OPT/ILES

      real t_u1(3,M1,M2),t_p1(M1,M2)        
      real t_u2(6,M1,M2),t_p2(M1,M2)
      real t_u3(3,M1,M2),t_u4(3,M1,M2)
      real t_vor1(3,M1,M2),t_vor2(3,M1,M2)    
      real t_up2(3,M1,M2),t_up3(3,M1,M2),t_up4(3,M1,M2)

      real t_tau(M1,M2)
      real t_ev1(M1,M2),t_ev2(M1,M2)
      real t_esgs(M1,M2),t_emol(M1,M2)
      real t_f_nega_ev(M1,M2)

      character*12 filename

      open(10,file='t_staticstics.dat',status='old'
     &       ,form='unformatted' 
     &       ,action='read')

      read(10) 
     &  (((t_u1(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u
     & ,(((t_u2(nv,i,j)  ,nv=1,6),i=1,n1m),j=1,n2m)           ! u"u"  
     & ,(((t_u3(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^3
     & ,(((t_u4(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^4
     & ,(((t_vor1(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega
     & ,(((t_vor2(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega"^2
     & ,((t_p1(i,j)              ,i=1,n1m),j=1,n2m)           ! p
     & ,((t_p2(i,j)              ,i=1,n1m),j=1,n2m)           ! p"^2
     & ,(((t_up2(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'u' ; u'=u-U, cf. u"=u-<u>
     & ,(((t_up3(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'^3
     & ,(((t_up4(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'^4

      close(10)

      t_tau(:,:)=0.0

      if(ILES.eq.1) then

      open(10,file='t_staticstics_les.dat',status='old'
     &       ,form='unformatted' 
     &       ,action='read')
      read(10)
     &  ((t_tau(i,j),i=1,n1m),j=1,n2m)            ! tau_12
     & ,((t_ev1(i,j),i=1,n1m),j=1,n2m)            ! nu_t
     & ,((t_ev2(i,j),i=1,n1m),j=1,n2m)            ! nu_t"^2
     & ,((t_esgs(i,j),i=1,n1m),j=1,n2m)           ! SGS dissipation
     & ,((t_emol(i,j),i=1,n1m),j=1,n2m)           ! molecular dissipation
     & ,((t_f_nega_ev(i,j),i=1,n1m),j=1,n2m)      ! fraction of negative nu_t  

      close(10)
 
      endif
c---- get utau
      do i=1,n1m
       utau(i)=sqrt(abs(t_u1(1,i,1)/(0.5*dy(1))/re))
      enddo

      open(10,file='../../../nb/02_main/post3/utau.no'
     &             ,status='unknown')  ! read no forcing u_tau
      do i=1,n1m
c      write(10,400) x(i),utau_o(i)
      read(10,400) dummy,utau_o(i)
      enddo
      close(10)

c----- along x parameter
c-Fig.5

      open(10,file='along_x.dat',status='unknown')
      open(11,file='along_x2.dat',status='unknown')
      open(20,file='../../../nb/02_main/post3/along_no.bin'
     &              ,status='unknown',form='unformatted')

      write(10,*) 
     & 'variables=xx,xp,c_f/2,H,Re_`q,p_w,p_w",`-uv_m_a_x,vor1max'
      write(11,*) 
     & 'variables=xx,xp,rcf,p_w,rprms,ruvmax,rvor1max'

      do i=1,n1m
       xc=0.5*(x(i)+x(i+1))
          u_inf=1.0
          delsta=0.0
          theta =0.0
          uv_max=0.0
       do j=1,n2m
          delsta=delsta + (u_inf-t_u1(1,i,j))*dy(j)
          theta =theta + t_u1(1,i,j)*(u_inf-t_u1(1,i,j))*dy(j) 
          if (t_u1(1,i,j).gt.0.99*u_inf) goto 888
          uv_max=amax1(uv_max,abs(t_u2(4,i,j)+t_tau(i,j)))
       enddo    

888    continue

          vor1_max=0.0
       do j=4,n2m-4
          if (t_vor2(1,i,j).gt.t_vor2(1,i,j-1)
     &      .and. 
     &        t_vor2(1,i,j).gt.t_vor2(1,i,j+1)) then
            vor1_max=amax1(vor1_max,abs(t_vor2(1,i,j)))
          endif
       enddo    


                   x_slot=75.0
                   u_tau_in=0.0492

       write(10,401) (xc-x_slot),(xc-x_slot)*u_tau_in*re
     &            ,utau(i)**2
     &            ,delsta/theta
     &            ,theta*re
     &            ,t_p1(i,1),sqrt(t_p2(i,1)),uv_max
     &            ,sqrt(vor1_max)

       read(20) xi,xi_p
     &          ,cfh_o
     &          ,H_o
     &          ,theta_o
     &          ,pw_o,pwrms_o,uvmax_o
     &          ,vor1max_o

       write(11,402) (xc-x_slot),(xc-x_slot)*u_tau_in*re
     &            ,utau(i)**2/cfh_o
     &            ,t_p1(i,1)
     &            ,sqrt(t_p2(i,1))/pwrms_o
     &            ,uv_max/uvmax_o
     &            ,sqrt(vor1_max)/vor1max_o

402       format(7(e12.5,x))
401       format(9(e12.5,x))
400       format(2(e12.5,x))

      enddo
      close(10)
      close(11)
      close(20)

 
      return
      end

c--------------- output_phase_2d --------------------------------
      subroutine output_phase_2d
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      common/moni/moni_i(100),x_moni(100),num_moni
      COMMON/PARA/RE
      COMMON/LES_OPT/ILES
      common/wall/utau(m1),utau_o(m1)

      real VM(3,M1,M2),PM(M1,M2)
      real VVM(6,M1,M2),PPM(M1,M2)
      real V3M(3,M1,M2),V4M(3,M1,M2)
      real DISSP(6,M1,M2)
      real TRANS(6,2,M1,M2)
      real PHI(6,M1,M2)
      real VORM(3,M1,M2),VORQM(3,M1,M2)    
      real REDIS(6,M1,M2)

      real t_u1_o(3,M1,M2),t_p1_o(M1,M2)        
      real t_u2_o(6,M1,M2),t_p2_o(M1,M2)
      real t_u3_o(3,M1,M2),t_u4_o(3,M1,M2)
      real t_vor1_o(3,M1,M2),t_vor2_o(3,M1,M2)    
      real t_up2_o(3,M1,M2),t_up3_o(3,M1,M2),t_up4_o(3,M1,M2)

      real t_u1(3,M1,M2),t_p1(M1,M2)
      real t_u2(6,M1,M2),t_p2(M1,M2)
      real t_u3(3,M1,M2),t_u4(3,M1,M2)
      real t_vor1(3,M1,M2),t_vor2(3,M1,M2)
      real t_up2(3,M1,M2),t_up3(3,M1,M2),t_up4(3,M1,M2)

      character*30 file_out

      open(10,file='../../../nb/02_main/post3/t_staticstics.dat'
     &       ,status='old'
     &       ,form='unformatted'
     &       ,action='read')

      read(10)
     &  (((t_u1_o(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u
     & ,(((t_u2_o(nv,i,j)  ,nv=1,6),i=1,n1m),j=1,n2m)           ! u"u"
     & ,(((t_u3_o(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^3
     & ,(((t_u4_o(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^4
     & ,(((t_vor1_o(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega
     & ,(((t_vor2_o(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega"^2
     & ,((t_p1_o(i,j)              ,i=1,n1m),j=1,n2m)           ! p
     & ,((t_p2_o(i,j)              ,i=1,n1m),j=1,n2m)           ! p"^2
     & ,(((t_up2_o(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'u' ; u'=u-U, cf. u"=u-<u>
     & ,(((t_up3_o(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'^3
     & ,(((t_up4_o(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'^4

      close(10)


      open(10,file='t_staticstics.dat',status='old'
     &       ,form='unformatted'
     &       ,action='read')

      read(10)
     &  (((t_u1(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u
     & ,(((t_u2(nv,i,j)  ,nv=1,6),i=1,n1m),j=1,n2m)           ! u"u"
     & ,(((t_u3(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^3
     & ,(((t_u4(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^4
     & ,(((t_vor1(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega
     & ,(((t_vor2(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega"^2
     & ,((t_p1(i,j)              ,i=1,n1m),j=1,n2m)           ! p
     & ,((t_p2(i,j)              ,i=1,n1m),j=1,n2m)           ! p"^2
     & ,(((t_up2(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'u' ; u'=u-U, cf. u"=u-<u>
     & ,(((t_up3(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'^3
     & ,(((t_up4(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'^4

      close(10)


      moni_phase = 8

      do 3 k_phase=0,moni_phase-1

        iphase=int((mt+1)*k_phase/moni_phase)
        write(*,*) iphase

        call read_phaseavg(iphase
     &           ,vm,pm,vvm,ppm,v3m,v4m  
     &           ,dissp,trans,phi,vorm,vorqm,redis)


      file_out='2d'
      NN=INDEX(file_out,'d')
      WRITE(UNIT=file_out(NN+1:),FMT='(BN,a6,I2.2,a4)') 
     & '_phase',k_phase,'.dat'
      open(10,file=file_out,status='unknown')
      write(10,301) 
301   format('variables=x,y,u,v,`w_z')

      write(10,102) n1m,n2m
102   format('zone i=',i5,',j=',i5,',f=point')
       do j=1,n2m
        do i=1,n1m
         xc=0.5*(x(i)+x(i+1))
         yc=0.5*(y(j)+y(j+1))

                   x_slot=75.0
                   u_tau_in=0.0492
c                   x_slot=78.90625
c                   u_tau_in=0.0525

         write(10,103)  (xc-x_slot)*u_tau_in*re
     &                , yc*utau_o(i)*re
     &               ,vm(1,i,j)
     &               ,vm(2,i,j)
     &               ,vorm(3,i,j)-t_vor1_o(3,i,j)
c     &               ,vorm(3,i,j)
c     &               ,sqrt(vorqm(1,i,j)-vorm(1,i,j)**2)
       enddo
      enddo

103   format(5(e12.5,x))
      close(10)


c---
      file_out='w1_2d'
      NN=INDEX(file_out,'d')
      WRITE(UNIT=file_out(NN+1:),FMT='(BN,a6,I2.2,a4)') 
     & '_phase',k_phase,'.dat'
      open(10,file=file_out,status='unknown')
      write(10,302) 
302   format('variables=x,y,`w_x_r_m_s')

      write(10,102) n1m,n2m
c102   format('zone i=',i5,',j=',i5,',f=point')
       do j=1,n2m
        do i=1,n1m
         xc=0.5*(x(i)+x(i+1))
         yc=0.5*(y(j)+y(j+1))

                   x_slot=75.0
                   u_tau_in=0.0492
c                   x_slot=78.90625
c                   u_tau_in=0.0525

         write(10,104)  (xc-x_slot)*u_tau_in*re
     &                , yc*utau_o(i)*re
     &               ,sqrt(vorqm(1,i,j)-vorm(1,i,j)**2)
       enddo
      enddo

104   format(3(e12.5,x))
      close(10)

c----
      file_out='2d_vel'
      NN=INDEX(file_out,'l')
      WRITE(UNIT=file_out(NN+1:),FMT='(BN,a6,I2.2,a4)')
     & '_phase',k_phase,'.dat'
      open(10,file=file_out,status='unknown')
      write(10,304)
304   format('variables=x,y,u,v,ut,vt,vor3,vor1')

      write(10,102) n1m-2,n2m-2
c102   format('zone i=',i5,',j=',i5,',f=point')
       do j=2,n2m-1
        do i=2,n1m-1

         xc=0.5*(x(i)+x(i+1))
         yc=0.5*(y(j)+y(j+1))

                   x_slot=75.0
                   u_tau_in=0.0492
c                   x_slot=78.90625
c                   u_tau_in=0.0525

         write(10,106)  (xc-x_slot)*u_tau_in*re
     &                , yc*utau_o(i)*re
     &   ,vm(1,i,j)
     &   ,vm(2,i,j)
     &   ,vm(1,i,j)-t_u1(1,i,j)
     &   ,vm(2,i,j)-t_u1(2,i,j)
     &   ,vorm(3,i,j)
     &   ,sqrt(vorqm(1,i,j)-vorm(1,i,j)**2)

       enddo
      enddo

106   format(7(e12.5,x))
      close(10)

3     continue

      return
      end
      
c--------------- output_phase_moni_x --------------------------------
      subroutine output_phase_moni_x
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/LES_OPT/ILES
      common/moni/moni_i(100),x_moni(100),num_moni

      real VM(3,M1,M2),PM(M1,M2)
      real VVM(6,M1,M2),PPM(M1,M2)
      real V3M(3,M1,M2),V4M(3,M1,M2)
      real DISSP(6,M1,M2)
      real TRANS(6,2,M1,M2)
      real PHI(6,M1,M2)
      real VORM(3,M1,M2),VORQM(3,M1,M2)    
      real REDIS(6,M1,M2)

      character*30 filename

      moni_phase = 8

      do 3 k_phase=0,moni_phase-1

        iphase=int((mt+1)*k_phase/moni_phase)
        write(*,*) iphase

        call read_phaseavg(iphase
     &           ,vm,pm,vvm,ppm,v3m,v4m  
     &           ,dissp,trans,phi,vorm,vorqm,redis)


        k_moni_x=1

      filename='phase_x'
      nn=index(filename,'x')
      write(unit=filename(nn+1:),fmt='(bn,i2.2,a3,i2.2,a4)') 
     &      k_moni_x,'_ph',k_phase,'.dat' 
      open(10,file=filename,status='unknown')
       do j=1,n2m
          i=moni_i(k_moni_x)
          yc=0.5*(y(j+1)+y(j))
          write(10,200) yc
     &                 ,vm(1,i,j)
200       format(2(e12.5,x))
       enddo
      close(10)


3     continue

      return
      end

c*************** output_2d **************************
      subroutine output_2d

      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m
      COMMON/MESH1/DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2)
      COMMON/MESH4/X(0:M1),DX(0:M1),G(M1)
      COMMON/PARA/RE
      COMMON/LES_OPT/ILES
      common/wall/utau(m1),utau_o(m1)
      common/moni/moni_i(100),x_moni(100),num_moni

      real t_u1(3,M1,M2),t_p1(M1,M2)        
      real t_u2(6,M1,M2),t_p2(M1,M2)
      real t_u3(3,M1,M2),t_u4(3,M1,M2)
      real t_vor1(3,M1,M2),t_vor2(3,M1,M2)    
      real t_up2(3,M1,M2),t_up3(3,M1,M2),t_up4(3,M1,M2)

      real t_u1_o(3,M1,M2),t_p1_o(M1,M2)        
      real t_u2_o(6,M1,M2),t_p2_o(M1,M2)
      real t_u3_o(3,M1,M2),t_u4_o(3,M1,M2)
      real t_vor1_o(3,M1,M2),t_vor2_o(3,M1,M2)    
      real t_up2_o(3,M1,M2),t_up3_o(3,M1,M2),t_up4_o(3,M1,M2)

      real o_u1(3,M1,M2),o_p1(M1,M2)
      real o_u2(6,M1,M2),o_p2(M1,M2)
      real o_u3(3,M1,M2),o_u4(3,M1,M2)
      real o_vor1(3,M1,M2),o_vor2(3,M1,M2)    

      real t_tau(M1,M2)
      real t_ev1(M1,M2),t_ev2(M1,M2)
      real t_esgs(M1,M2),t_emol(M1,M2)
      real t_f_nega_ev(M1,M2)
 
      character*12 filename

      open(10,file='t_staticstics.dat',status='old'
     &       ,form='unformatted' 
     &       ,action='read')

      read(10) 
     &  (((t_u1(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u
     & ,(((t_u2(nv,i,j)  ,nv=1,6),i=1,n1m),j=1,n2m)           ! u"u"  
     & ,(((t_u3(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^3
     & ,(((t_u4(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^4
     & ,(((t_vor1(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega
     & ,(((t_vor2(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega"^2
     & ,((t_p1(i,j)              ,i=1,n1m),j=1,n2m)           ! p
     & ,((t_p2(i,j)              ,i=1,n1m),j=1,n2m)           ! p"^2
     & ,(((t_up2(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'u' ; u'=u-U, cf. u"=u-<u>
     & ,(((t_up3(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'^3
     & ,(((t_up4(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'^4

      close(10)

      open(10,file='../../../nb/02_main/post3/t_staticstics.dat'
     &       ,status='old'
     &       ,form='unformatted'
     &       ,action='read')

      read(10)
     &  (((t_u1_o(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u
     & ,(((t_u2_o(nv,i,j)  ,nv=1,6),i=1,n1m),j=1,n2m)           ! u"u"
     & ,(((t_u3_o(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^3
     & ,(((t_u4_o(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^4
     & ,(((t_vor1_o(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega
     & ,(((t_vor2_o(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega"^2
     & ,((t_p1_o(i,j)              ,i=1,n1m),j=1,n2m)           ! p
     & ,((t_p2_o(i,j)              ,i=1,n1m),j=1,n2m)           ! p"^2
     & ,(((t_up2_o(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'u' ; u'=u-U, cf. u"=u-<u>
     & ,(((t_up3_o(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'^3
     & ,(((t_up4_o(nv,i,j) ,nv=1,3),i=1,n1m),j=1,n2m)           ! u'^4

      close(10)

      open(10,file='o_staticstics.dat',status='old'
     &       ,form='unformatted' 
     &       ,action='read')

      read(10) 
     &  (((o_u1(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u
     & ,(((o_u2(nv,i,j)  ,nv=1,6),i=1,n1m),j=1,n2m)           ! u"u"
     & ,(((o_u3(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^3
     & ,(((o_u4(nv,i,j)  ,nv=1,3),i=1,n1m),j=1,n2m)           ! u"^4
     & ,(((o_vor1(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega
     & ,(((o_vor2(nv,i,j),nv=1,3),i=1,n1m),j=1,n2m)           ! omega"^2
     & ,((o_p1(i,j)              ,i=1,n1m),j=1,n2m)           ! p
     & ,((o_p2(i,j)              ,i=1,n1m),j=1,n2m)           ! p"^2

      close(10)


      if(ILES.eq.1) then

      open(10,file='t_staticstics_les.dat',status='old'
     &       ,form='unformatted' 
     &       ,action='read')
      read(10)
     &  ((t_tau(i,j),i=1,n1m),j=1,n2m)            ! tau_12
     & ,((t_ev1(i,j),i=1,n1m),j=1,n2m)            ! nu_t
     & ,((t_ev2(i,j),i=1,n1m),j=1,n2m)            ! nu_t"^2
     & ,((t_esgs(i,j),i=1,n1m),j=1,n2m)           ! SGS dissipation
     & ,((t_emol(i,j),i=1,n1m),j=1,n2m)           ! molecular dissipation
     & ,((t_f_nega_ev(i,j),i=1,n1m),j=1,n2m)      ! fraction of negative nu_t  

      close(10)
 
      endif


c---   
      open(10,file='2d_tilde.dat',status='unknown')
      write(10,101) n1m,n2m
101   format('zone i=',i5,',j=',i5,',f=point')
       do j=1,n2m
        do i=1,n1m
         xx=0.5*(x(i)+x(i+1))
         yy=0.5*(y(j)+y(j+1))
       
         write(10,100) xx,yy
     &                ,o_u1(1,i,j)**0.5
     &                ,o_u1(2,i,j)**0.5/(0.4/sqrt(2.))  ! normalized by A/sqrt(2)
     &                ,o_vor1(3,i,j)**0.5
        enddo
       enddo

100   format(5(e12.5,x))
      close(10)

c---   
      open(10,file='2d_tilde_p.dat',status='unknown')
      write(10,101) n1m,n2m
       do j=1,n2m
        do i=1,n1m
         xx=0.5*(x(i)+x(i+1))
         yy=0.5*(y(j)+y(j+1))
       
         write(10,100) xx,yy
     &                ,o_p1(i,j)**0.5
        enddo
       enddo

      close(10)

c--
      open(10,file='2d_diff.dat',status='unknown')
      write(10,101) n1m,n2m

       do j=1,n2m
        do i=1,n1m
         xc=0.5*(x(i)+x(i+1))
         yc=0.5*(y(j)+y(j+1))
       
                   x_slot=75.0
                   u_tau_in=0.0492
c                   x_slot=78.90625
c                   u_tau_in=0.0525

         write(10,401)  (xc-x_slot)*u_tau_in*re
     &                , yc*utau_o(i)*re
     &                ,(t_u2(1,i,j)**0.5-t_u2_o(1,i,j)**0.5)/utau_o(i)
     &                ,(t_u2(2,i,j)**0.5-t_u2_o(2,i,j)**0.5)/utau_o(i)
     &                ,(t_u2(3,i,j)**0.5-t_u2_o(3,i,j)**0.5)/utau_o(i)
     &                ,(-t_u2(4,i,j)     +t_u2_o(4,i,j) )/utau_o(i)**2    


        enddo
       enddo

401   format(6(e12.5,x))
      close(10)



      return
      end
