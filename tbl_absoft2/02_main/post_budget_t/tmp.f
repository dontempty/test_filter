c
c    The postprocessing program of periodically perturbed t.b.l.
c
c
c
c    2. phase mean budget analysis 
c       : main/phase.***  
c        
c
c                                                   Kyoungyoun KIM
c                                          Flow Control Laboratory
c                             Department of Mechanical Engineering
c                                                            KAIST 



      program budget

      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/moni/ix(100),x_moni(100),num_moni
      common/fileio/ith,k_phase
      common/moni2/moni_ph

      call setup


c      num_moni_ph = 8
c
c      do 3 k_phase=0,num_moni_ph-1
c
c        moni_ph=int((mt+1)*k_phase/num_moni_ph)
c         write(*,*) k_phase,moni_ph


c---- calculate the fluctuating quantities
c      call budget_term 

c---- for time-averaged budget : time averaing of u~_i*u~_j : oscilating term of mean budget
c      call oscil_term 

c---- for phase-averaed budget : time derivatives of <u_i>, <u_i"*u_j"> : acc term of budget
c      call acc_term  



      call read_budget_term_t 
c      call read_budget_term_p 
      call read_u_tau_no
c      goto 1

      call set_moni

      do ith=1,num_moni  ! at every monitoring x position
         call u1budget(ix(ith))
         call u2budget(ix(ith))
         call uu1budget(ix(ith))
         call uu2budget(ix(ith))
         call uu3budget(ix(ith))
         call uu4budget(ix(ith))
c         call k_budget(ix(ith))
         call press_strain(ix(ith))
      enddo

1     continue
c      call two_d_uu1
c      call two_d_uu2
c      call two_d_uu4

c3     continue

      stop
      end



c*************** set_moni ***************************************
c     find i-index corresponding on given monitoring x location
      subroutine set_moni
      include '../PARAM.H'
      COMMON/PARA/RE
      common/dim/n1,n2,n3,n1m,n2m,n3m
      common/fileio/ith,k_phase
      common/mesh1/dx1,dx1q,dx3,dx3q

      common/moni/ix(100),x_moni(100),num_moni


      real x(m1)


      open(10,file='moni_x.set',status='old')
       read (10,*) num_moni
      do k=1,num_moni
       read (10,*) x_moni(k)
      enddo
      close(10)


      do i=1,n1
        x(i)=real(i-1)/dx1
      enddo


                   x_slot=75.0
                   u_tau_in=0.0492

      do 1 k=1,num_moni

        x_moni(k)=x_moni(k)/RE/u_tau_in+x_slot

        do i=1,n1m-1
         temp=(x_moni(k)-0.5*(x(i  )+x(i+1)))
     &       *(x_moni(k)-0.5*(x(i+1)+x(i+2)))

         if(temp.le.0.0) then
           ix(k)=i
           goto 1
         endif
        enddo

1     continue

      do k=1,num_moni
       write(*,*) k,x_moni(k),ix(k)
     &       ,0.5*(x(ix(k))+x(ix(k)+1))
     &       ,(0.5*(x(ix(k))+x(ix(k)+1))-x_slot)*u_tau_in*RE

      enddo

      return
      end


c*************** read_u_tau_no ***********************     
      subroutine read_u_tau_no      

      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/fileio/ith,k_phase
      common/para/re
      common/utau/u_tau_no(m1)

      open(10,file='../../../nb/02_main/post3/utau.no'
     &             ,status='unknown')  ! read no forcing u_tau

      do i=1,n1m
      read(10,400) dummy,u_tau_no(i)
      enddo
      close(10)

400   format(2(e12.5,x))


      return
      end


C***************** SETUP ***********************
      SUBROUTINE SETUP
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME
c      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      COMMON/PARA/RE
      COMMON/VPERIN/VPER
c      COMMON/SIZE/ALX,ALY,ALZ
      common/size/alx,aly,alz,vol
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

      n2mm = n2 - 2

      call mesh
      call indices 

      return
      end



c*************** mesh ***********************     
      subroutine mesh
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/fileio/ith,k_phase
      common/mesh1/dx1,dx1q,dx3,dx3q
      common/mesh2/y(0:m2)
      common/mesh3/dy(0:m2),h(m2),hm(m2),hc(m2),hp(m2)
      common/mesh4/dym(m2),dyc(m2),dyp(m2)
      common/size/alx,aly,alz,vol


c     create the tangent hyperbolic grid in x2 direction      

         Y(0)=0.0
      OPEN(10,FILE='../GRID.PLT',STATUS='UNKNOWN')
      DO J=1,N2
      READ(10,*) DUMMY,Y(J),DUMMY
      ENDDO
      CLOSE(10)

      vol = alx*aly*alz 


      dx1 = dble(n1m)/alx
      dx3 = dble(n3m)/alz

      dx1q = dx1**2.0
      dx3q = dx3**2.0


         dy(1) = y(2)
         h(1)  = 0.5*dy(1)
      do j=2,n2m
         dy(j) = y(j+1)-y(j)
         h(j)  = 0.5*(dy(j)+dy(j-1))
      enddo
         h(n2) = dy(n2m)*0.5


      do j=1,n2m
         hp(j)=2.0/h(j+1)/(h(j)+h(j+1))
         hc(j)=2.0/h(j+1)/h(j)
         hm(j)=2.0/h(j)/(h(j)+h(j+1))
      enddo


      do j=2,n2m
         dyp(j)=2.0/dy(j)/(dy(j-1)+dy(j))
         dyc(j)=2.0/dy(j)/dy(j-1)
         dym(j)=2.0/dy(j-1)/(dy(j-1)+dy(j))
      enddo


      return
      end

c*************** indices ***********************     
      subroutine indices
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/fileio/ith,k_phase
      common/indx1/ipa(m1),imu(m1),imv(m1),kpa(m3),kma(m3)
      common/indx2/jpa(m2),jmu(m2),jmv(m2)

      do i=1,n1m
         ipa(i)   = i + 1
         imu(i)   = i - 1
         imv(i)   = i - 1
      enddo
         ipa(n1m) = n1m
         imu(2)   = 2
         imv(1)   = 1


      do k=1,n3m
         kpa(k)   = k + 1
         kma(k)   = k - 1
      enddo
         kpa(n3m) = 1
         kma(1)   = n3m


      do j=1,n2m
         jpa(j)   = j + 1
         jmu(j)   = j - 1
         jmv(j)   = j - 1
      enddo
         jpa(n2m) = n2m
         jmu(1)   = 1
         jmv(2)   = 2

      
      return
      end
   
c*************** read_budget_term_p *******************************************
      subroutine read_budget_term_p 
     
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/fileio/ith,k_phase
      common/indx1/ipa(m1),imu(m1),imv(m1),kpa(m3),kma(m3)
      common/indx2/jpa(m2),jmu(m2),jmv(m2)
      common/mesh1/dx1,dx1q,dx3,dx3q
      common/mesh2/y(0:m2)
      common/mesh3/dy(0:m2),h(m2),hm(m2),hc(m2),hp(m2)
      common/mesh4/dym(m2),dyc(m2),dyp(m2)

      common/budget1/u(3,0:m1,0:m2),p(0:m1,0:m2),uu(6,m1,m2)
     >              ,uuu(6,2,m1,m2)
      common/budget2/upg(6,m1,m2),diss(6,m1,m2),presstrain(6,m1,m2)
      common/budget3/pp(0:m1,0:m2)
      common/budget4/acc_u1(3,0:m1,0:m2)
      common/budget5/acc_u2(6,m1,m2)
      common/budget6/v1_osciltm(2,m1,m2)      ! overbar{u~_i*u~_j}
     >              ,v2_osciltm_P(6,m1,m2)    ! Production term
     >              ,v2_osciltm_C(6,m1,m2)    ! Convection term
      common/budget7/rss_os(2,m1,m2)


c---- set oscilating term to zeo for phase-averaged budget

      v1_osciltm(:,:,:)=0.0
      v2_osciltm_P(:,:,:)=0.0
      v2_osciltm_C(:,:,:)=0.0
      rss_os(:,:,:)=0.0

c---- read the budget term from data file

      open(10,file='budget.avg',form='unformatted',status='unknown')

      do 1 i=1,n1m
      do 1 j=1,n2m

      read(10) (u(l,i,j),l=1,3)
      read(10) (uu(l,i,j),l=1,6)
      read(10) p(i,j),pp(i,j)
      read(10) ((uuu(l,m,i,j),m=1,2),l=1,6)
      read(10) (upg(l,i,j),l=1,6)
      read(10) (diss(l,i,j),l=1,6)
      read(10) (presstrain(l,i,j),l=1,6)

1     continue

      close(10)

c---- read the derivativesi to 'acc_term.avg'

      open(10,file='acc_term.avg',form='unformatted',status='unknown')

      do 6 i=2,n1m-1
      do 6 j=2,n2m-1
         read(10) (acc_u1(l,i,j),l=1,2)
         read(10) (acc_u2(l,i,j),l=1,6)
6     continue

      close(10)




      return
      end

 
c*************** read_budget_term_t *******************************************
      subroutine read_budget_term_t 
     
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/fileio/ith,k_phase
      common/indx1/ipa(m1),imu(m1),imv(m1),kpa(m3),kma(m3)
      common/indx2/jpa(m2),jmu(m2),jmv(m2)
      common/mesh1/dx1,dx1q,dx3,dx3q
      common/mesh2/y(0:m2)
      common/mesh3/dy(0:m2),h(m2),hm(m2),hc(m2),hp(m2)
      common/mesh4/dym(m2),dyc(m2),dyp(m2)

      common/budget1/u(3,0:m1,0:m2),p(0:m1,0:m2),uu(6,m1,m2)
     >              ,uuu(6,2,m1,m2)
      common/budget2/upg(6,m1,m2),diss(6,m1,m2),presstrain(6,m1,m2)
      common/budget3/pp(0:m1,0:m2)
      common/budget4/acc_u1(3,0:m1,0:m2)
      common/budget5/acc_u2(6,m1,m2)
      common/budget6/v1_osciltm(2,m1,m2)      ! overbar{u~_i*u~_j}
     >              ,v2_osciltm_P(6,m1,m2)    ! Production term
     >              ,v2_osciltm_C(6,m1,m2)    ! Convection term
      common/budget7/rss_os(2,m1,m2)


c---- set acceleration term to zeo for time-averaged budget

      acc_u1(:,:,:)=0.0
      acc_u2(:,:,:)=0.0

c---- read the budget term from data file

      open(10,file='budget.avg',form='unformatted',status='unknown')

      do 1 i=1,n1m
      do 1 j=1,n2m

      read(10) (u(l,i,j),l=1,3)
      read(10) (uu(l,i,j),l=1,6)
      read(10) p(i,j),pp(i,j)
      read(10) ((uuu(l,m,i,j),m=1,2),l=1,6)
      read(10) (upg(l,i,j),l=1,6)
      read(10) (diss(l,i,j),l=1,6)
      read(10) (presstrain(l,i,j),l=1,6)

1     continue

      close(10)



c---  read  the time averaged data 

      open(10,file='oscil_term.avg',form='unformatted',status='unknown') 

      do 3 i=2,n1m-1
      do 3 j=2,n2m-1
         read(10) (v1_osciltm(l,i,j),l=1,2)
         read(10) (v2_osciltm_P(l,i,j),l=1,6)
         read(10) (v2_osciltm_C(l,i,j),l=1,6)
         read(10) (rss_os(l,i,j),l=1,2)
3     continue
      close(10)


      return
      end

 
c*************** u1budget **********************************************
      subroutine u1budget(moni_x)
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/fileio/ith,k_phase
      common/indx1/ipa(m1),imu(m1),imv(m1),kpa(m3),kma(m3)
      common/indx2/jpa(m2),jmu(m2),jmv(m2)
      common/mesh1/dx1,dx1q,dx3,dx3q
      common/mesh2/y(0:m2)
      common/mesh3/dy(0:m2),h(m2),hm(m2),hc(m2),hp(m2)
      common/mesh4/dym(m2),dyc(m2),dyp(m2)
      common/para/re
      common/size/alx,aly,alz,vol

      common/budget1/u(3,0:m1,0:m2),p(0:m1,0:m2),uu(6,m1,m2)
     >              ,uuu(6,2,m1,m2)
      common/budget2/upg(6,m1,m2),diss(6,m1,m2),presstrain(6,m1,m2)
      common/budget3/pp(0:m1,0:m2)
      common/budget4/acc_u1(3,0:m1,0:m2)
      common/budget5/acc_u2(6,m1,m2)
      common/budget6/v1_osciltm(2,m1,m2)      ! overbar{u~_i*u~_j}
     >              ,v2_osciltm_P(6,m1,m2)    ! Production term
     >              ,v2_osciltm_C(6,m1,m2)    ! Convection term

      real convec(m2)
      real viscos(m2)
      real reynol(m2)
      real pressg(m2)

      integer moni_x

      character*30 name

      name='budet_u1x'
      nn=index(name,'x')
      write(unit=name(nn+1:),fmt='(bn,i2.2,a1,i1.1,a4)')
     >           ith,'p',k_phase,'.dat'

      open(10,file=name,status='unknown')
      write(10,*) 
     > 'variables=y/`q_i_n,C,C_o_s,VD,RSS,pg,balance'

      i =moni_x  ! calculating the budget at this x position
      ip=i+1
      im=i-1
     
      do j=2,n2mm
      jp=j+1
      jm=j-1
 

c---- convection term      
      convec(j)=-(
     >   (u(1,ip,j)*u(1,ip,j)-u(1,im,j)*u(1,im,j))
     >    *dx1/2.
     > +
     >    ((h(j)**2)*u(1,i,jp)*u(2,i,jp)
     >   -(h(j)**2-h(jp)**2)*u(1,i,j)*u(2,i,j)
     >   -(h(jp)**2)*u(1,i,jm)*u(2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     > )                                  ! -(u_1*u_k),k

c---- viscous term
      viscos(j)=1./re*(
     >   (u(1,ip,j)-2.*u(1,i,j)+u(1,im,j))*dx1q
     >  +(hp(j)*u(1,i,jp)
     >   -hc(j)*u(1,i,j )
     >   +hm(j)*u(1,i,jm))
     > ) 

c---- reynolds shear stress term
      reynol(j)=-(
     >   (uu(1,ip,j)-uu(1,im,j))*dx1/2.
     > +  (h(j)**2*uu(4,i,jp)
     >   -(h(j)**2-h(jp)**2)*uu(4,i,j)
     >   -h(jp)**2*uu(4,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     > )                                  ! -<u_1u_k>,k

c---- pressuer gradient term
      pressg(j)=-(p(ip,j)-p(im,j))*dx1/2.      

      sum=convec(j)-acc_u1(1,i,j)+viscos(j)+reynol(j)+pressg(j)
     >   +v1_osciltm(1,i,j)      ! overbar{u~_i*u~_j}

      write(10,100) 0.5*(y(j)+y(j+1))
     >             ,convec(j)-acc_u1(1,i,j)+v1_osciltm(1,i,j)
     >             ,v1_osciltm(1,i,j)
     >             ,viscos(j),reynol(j),pressg(j),sum 
100   format(7(e15.7,x))    
      enddo
      close(10)

      return
      end


c*************** u2budget **********************************************
      subroutine u2budget(moni_x)
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/fileio/ith,k_phase
      common/indx1/ipa(m1),imu(m1),imv(m1),kpa(m3),kma(m3)
      common/indx2/jpa(m2),jmu(m2),jmv(m2)
      common/mesh1/dx1,dx1q,dx3,dx3q
      common/mesh2/y(0:m2)
      common/mesh3/dy(0:m2),h(m2),hm(m2),hc(m2),hp(m2)
      common/mesh4/dym(m2),dyc(m2),dyp(m2)
      common/para/re
      common/size/alx,aly,alz,vol

      common/budget1/u(3,0:m1,0:m2),p(0:m1,0:m2),uu(6,m1,m2)
     >              ,uuu(6,2,m1,m2)
      common/budget2/upg(6,m1,m2),diss(6,m1,m2),presstrain(6,m1,m2)
      common/budget3/pp(0:m1,0:m2)
      common/budget4/acc_u1(3,0:m1,0:m2)
      common/budget6/v1_osciltm(2,m1,m2)      ! overbar{u~_i*u~_j}
     >              ,v2_osciltm_P(6,m1,m2)    ! Production term
     >              ,v2_osciltm_C(6,m1,m2)    ! Convection term

      real convec(m2)
      real viscos(m2)
      real reynol(m2)
      real pressg(m2)
      
      integer moni_x

      character*30 name

      name='budet_u2x'
      nn=index(name,'x')
      write(unit=name(nn+1:),fmt='(bn,i2.2,a1,i1.1,a4)')
     >           ith,'p',k_phase,'.dat'

      open(10,file=name,status='unknown')
      write(10,*) 
     > 'variables=y/`q_i_n,C,C_o_s,VD,RSS,pg,balance'

      i =moni_x  ! calculating the budget at this x position
      ip=i+1
      im=i-1
     
      do j=2,n2mm
      jp=j+1
      jm=j-1
 

c---- convection term      
      convec(j)=-(
     >   (u(1,ip,j)*u(2,ip,j)-u(1,im,j)*u(2,im,j))
     >    *dx1/2.
     > +
     >    ((h(j)**2)*u(2,i,jp)*u(2,i,jp)
     >   -(h(j)**2-h(jp)**2)*u(2,i,j)*u(2,i,j)
     >   -(h(jp)**2)*u(2,i,jm)*u(2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     > )                                  ! -(u_1*u_k),k

c---- viscous term
      viscos(j)=1./re*(
     >   (u(2,ip,j)-2.*u(2,i,j)+u(2,im,j))*dx1q
     >  +(hp(j)*u(2,i,jp)
     >   -hc(j)*u(2,i,j )
     >   +hm(j)*u(2,i,jm))
     > ) 

c---- reynolds shear stress term
      reynol(j)=-(
     >   (uu(4,ip,j)-uu(4,im,j))*dx1/2.
     > +  (h(j)**2*uu(2,i,jp)
     >   -(h(j)**2-h(jp)**2)*uu(2,i,j)
     >   -h(jp)**2*uu(2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     > )                                  ! -<u_1u_k>,k

c---- pressuer gradient term
      pressg(j)=-(      
     >    (h(j)**2*p(i,jp)
     >   -(h(j)**2-h(JP)**2)*p(i,j)
     >   -h(jp)**2*p(i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     >  )

      sum=convec(j)-acc_u1(2,i,j)+viscos(j)+reynol(j)+pressg(j)
     >   +v1_osciltm(2,i,j)      ! overbar{u~_i*u~_j}

      write(10,100) 0.5*(y(j)+y(j+1))
     >             ,convec(j)-acc_u1(2,i,j)+v1_osciltm(2,i,j)
     >             ,v1_osciltm(2,i,j)
     >             ,viscos(j),reynol(j),pressg(j),sum 

100   format(7(e15.7,x))    
      enddo
      close(10)

      return
      end


c*************** uu1budget *********************************************
      subroutine uu1budget(moni_x)
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/fileio/ith,k_phase
      common/indx1/ipa(m1),imu(m1),imv(m1),kpa(m3),kma(m3)
      common/indx2/jpa(m2),jmu(m2),jmv(m2)
      common/mesh1/dx1,dx1q,dx3,dx3q
      common/mesh2/y(0:m2)
      common/mesh3/dy(0:m2),h(m2),hm(m2),hc(m2),hp(m2)
      common/mesh4/dym(m2),dyc(m2),dyp(m2)
      common/para/re
      common/size/alx,aly,alz,vol
      common/utau/u_tau_no(m1)

      common/budget1/u(3,0:m1,0:m2),p(0:m1,0:m2),uu(6,m1,m2)
     >              ,uuu(6,2,m1,m2)
      common/budget2/upg(6,m1,m2),diss(6,m1,m2),presstrain(6,m1,m2)
      common/budget3/pp(0:m1,0:m2)
      common/budget5/acc_u2(6,m1,m2)
      common/budget6/v1_osciltm(2,m1,m2)      ! overbar{u~_i*u~_j}
     >              ,v2_osciltm_P(6,m1,m2)    ! Production term
     >              ,v2_osciltm_C(6,m1,m2)    ! Convection term

      real produc(m2)
      real dissip(m2)
      real visdif(m2)
      real convec(m2)
      real transp(m2)
      real vpgrad(m2)
      

      integer moni_x

      character*30 name

      name='budet_uu1x'
      nn=index(name,'x')
      write(unit=name(nn+1:),fmt='(bn,i2.2,a1,i1.1,a4)')
     >           ith,'p',k_phase,'.dat'

      open(10,file=name,status='unknown')
      write(10,11) 
11    format('variables=y,P,`e,D,C,T,`P,P_o_s,C_o_s,sum')

      name='rss1x'
      nn=index(name,'x')
      write(unit=name(nn+1:),fmt='(bn,i2.2,a1,i1.1,a4)')
     >           ith,'p',k_phase,'.dat'

      open(11,file=name,status='unknown')
      write(11,12) 
12    format('variables=y^+_o,P,P_o_s,Px,Py')


      i=moni_x
      ip=i+1
      im=i-1
     
      do j=2,n2mm
      jp=j+1
      jm=j-1

c     production term
      produc(j)=-2.*(
     >   uu(1,i,j)*(u(1,ip,j)-u(1,im,j))*dx1/2.  
     >  +uu(4,i,j)*  
     >    ((h(j)**2)*u(1,i,jp)
     >   -(h(j)**2-h(jp)**2)*u(1,i,j)
     >   -(h(jp)**2)*u(1,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     >  )                                 ! -2<u_1u_k>u1,k

c     dissipation term
      dissip(j)=-2./re*diss(1,i,j)

c     viscous diffusion term
      visdif(j)=1./re*(
     >   (uu(1,ip,j)-2.*uu(1,i,j)+uu(1,im,j))*dx1q
     >  +(hp(j)*uu(1,i,jp)
     >   -hc(j)*uu(1,i,j )
     >   +hm(j)*uu(1,i,jm))
     > ) 
     > 

c     convection term      
      convec(j)=-(
     >    u(1,i,j)*(uu(1,ip,j)-uu(1,im,j))*dx1/2.
     > +  u(2,i,j)*
     >    ((h(j)**2)*uu(1,i,jp)
     >   -(h(j)**2-h(jp)**2)*uu(1,i,j)
     >   -(h(jp)**2)*uu(1,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     >  )                                 ! -u_k*<u_1u_1>,k

c     turbulent transport term
      transp(j)=-(
     >    (uuu(1,1,ip,j)-uuu(1,1,im,j))*dx1/2.
     > +  
     >    ((h(j)**2)*uuu(1,2,i,jp)
     >   -(h(j)**2-h(jp)**2)*uuu(1,2,i,j)
     >   -(h(jp)**2)*uuu(1,2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     > )

c     velocity gradient term
      vpgrad(j)=-upg(1,i,j)


c---
      sum=produc(j)+dissip(j)+visdif(j)+convec(j)
     >   +transp(j)+vpgrad(j)
     >   -acc_u2(1,i,j)
     >   +v2_osciltm_P(1,i,j)
     >   +v2_osciltm_C(1,i,j)

      scale=1./(u_tau_no(i)**4*RE)

      yy=0.5*(y(j)+y(j+1))

      write(10,100) 
     >      yy*u_tau_no(i)*RE
     >     ,(produc(j)+v2_osciltm_P(1,i,j))*scale
     >     ,dissip(j)*scale
     >     ,visdif(j)*scale
     >     ,(convec(j)-acc_u2(1,i,j)+v2_osciltm_C(1,i,j))*scale
     >     ,transp(j)*scale
     >     ,vpgrad(j)*scale
     >     ,v2_osciltm_P(1,i,j)*scale
     >     ,v2_osciltm_C(1,i,j)*scale
     >     ,sum*scale

100   format(10(e12.5,x))    


      dudx = (u(1,ip,j)-u(1,im,j))*dx1/2.
      dudy =( (h(j)**2)*u(1,i,jp)
     >     -(h(j)**2-h(jp)**2)*u(1,i,j)
     >     -(h(jp)**2)*u(1,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))

      write(11,101)
     >      yy*u_tau_no(i)*RE
     >     ,(produc(j)+v2_osciltm_P(1,i,j))*scale
     >     ,v2_osciltm_P(1,i,j)*scale
     >     ,-2.*uu(1,i,j)*dudx*scale
     >     ,-2.*uu(4,i,j)*dudy*scale

101   format(5(e12.5,x))    

c---


      enddo
      close(10)
      close(11)

      return
      end


c*************** uu2budget *********************************************
      subroutine uu2budget(moni_x)
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/fileio/ith,k_phase
      common/indx1/ipa(m1),imu(m1),imv(m1),kpa(m3),kma(m3)
      common/indx2/jpa(m2),jmu(m2),jmv(m2)
      common/mesh1/dx1,dx1q,dx3,dx3q
      common/mesh2/y(0:m2)
      common/mesh3/dy(0:m2),h(m2),hm(m2),hc(m2),hp(m2)
      common/mesh4/dym(m2),dyc(m2),dyp(m2)
      common/para/re
      common/size/alx,aly,alz,vol
      common/utau/u_tau_no(m1)

      common/budget1/u(3,0:m1,0:m2),p(0:m1,0:m2),uu(6,m1,m2)
     >              ,uuu(6,2,m1,m2)
      common/budget2/upg(6,m1,m2),diss(6,m1,m2),presstrain(6,m1,m2)
      common/budget3/pp(0:m1,0:m2)
      common/budget5/acc_u2(6,m1,m2)
      common/budget6/v1_osciltm(2,m1,m2)      ! overbar{u~_i*u~_j}
     >              ,v2_osciltm_P(6,m1,m2)    ! Production term
     >              ,v2_osciltm_C(6,m1,m2)    ! Convection term

      real produc(m2)
      real dissip(m2)
      real visdif(m2)
      real convec(m2)
      real transp(m2)
      real vpgrad(m2)


      integer moni_x

      character*30 name

      name='budet_uu2x'
      nn=index(name,'x')
      write(unit=name(nn+1:),fmt='(bn,i2.2,a1,i1.1,a4)')
     >           ith,'p',k_phase,'.dat'

      open(10,file=name,status='unknown')
      write(10,11) 
11    format('variables=y,P,`e,D,C,T,`P,P_o_s,C_o_s,sum')

      name='rss2x'
      nn=index(name,'x')
      write(unit=name(nn+1:),fmt='(bn,i2.2,a1,i1.1,a4)')
     >           ith,'p',k_phase,'.dat'

      open(11,file=name,status='unknown')
      write(11,12) 
12    format('variables=y^+_o,P,P_o_s,Px,Py')

      i=moni_x
      ip=i+1
      im=i-1
     
      do j=2,n2mm
      jp=j+1
      jm=j-1

c     production term
      produc(j)=-2.*(
     >   uu(4,i,j)*(u(2,ip,j)-u(2,im,j))*dx1/2.  
     >  +uu(2,i,j)*  
     >    ((h(j)**2)*u(2,i,jp)
     >   -(h(j)**2-h(jp)**2)*u(2,i,j)
     >   -(h(jp)**2)*u(2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     >  )                                 ! -2<u_2u_k>u2,k

c     dissipation term
      dissip(j)=-2./re*diss(2,i,j)

c     viscous diffusion term
      visdif(j)=1./re*(
     >   (uu(2,ip,j)-2.*uu(2,i,j)+uu(2,im,j))*dx1q
     >  +(hp(j)*uu(2,i,jp)
     >   -hc(j)*uu(2,i,j )
     >   +hm(j)*uu(2,i,jm))
     > ) 
     > 

c     convection term      
      convec(j)=-(
     >    u(1,i,j)*(uu(2,ip,j)-uu(2,im,j))*dx1/2.
     > +  u(2,i,j)*
     >    ((h(j)**2)*uu(2,i,jp)
     >   -(h(j)**2-h(jp)**2)*uu(2,i,j)
     >   -(h(jp)**2)*uu(2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     >  )                                 ! -u_k*<u_2u_2>,k

c     turbulent transport term
      transp(j)=-(
     >    (uuu(2,1,ip,j)-uuu(2,1,im,j))*dx1/2.
     > +  
     >    ((h(j)**2)*uuu(2,2,i,jp)
     >   -(h(j)**2-h(jp)**2)*uuu(2,2,i,j)
     >   -(h(jp)**2)*uuu(2,2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     > )

c     velocity gradient term
      vpgrad(j)=-upg(2,i,j)

c---
      sum=produc(j)+dissip(j)+visdif(j)+convec(j)
     >   +transp(j)+vpgrad(j)
     >   -acc_u2(2,i,j)
     >   +v2_osciltm_P(2,i,j)
     >   +v2_osciltm_C(2,i,j)

      scale=1./(u_tau_no(i)**4*RE)

      yy=0.5*(y(j)+y(j+1))

      write(10,100) 
     >      yy*u_tau_no(i)*RE
     >     ,(produc(j)+v2_osciltm_P(2,i,j))*scale
     >     ,dissip(j)*scale
     >     ,visdif(j)*scale
     >     ,(convec(j)-acc_u2(2,i,j)+v2_osciltm_C(2,i,j))*scale
     >     ,transp(j)*scale
     >     ,vpgrad(j)*scale
     >     ,v2_osciltm_P(2,i,j)*scale
     >     ,v2_osciltm_C(2,i,j)*scale
     >     ,sum*scale

100   format(10(e12.5,x))    

      dvdx = (u(2,ip,j)-u(2,im,j))*dx1/2.
      dvdy =( (h(j)**2)*u(2,i,jp)
     >     -(h(j)**2-h(jp)**2)*u(2,i,j)
     >     -(h(jp)**2)*u(2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))


      write(11,101)
     >      yy*u_tau_no(i)*RE
     >     ,(produc(j)+v2_osciltm_P(2,i,j))*scale
     >     ,v2_osciltm_P(2,i,j)*scale
     >     ,-2.*uu(4,i,j)*dvdx*scale
     >     ,-2.*uu(2,i,j)*dvdy*scale

101   format(5(e12.5,x))    
c---

      enddo
      close(10)

      return
      end


c*************** uu3budget *********************************************
      subroutine uu3budget(moni_x)
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/fileio/ith,k_phase
      common/indx1/ipa(m1),imu(m1),imv(m1),kpa(m3),kma(m3)
      common/indx2/jpa(m2),jmu(m2),jmv(m2)
      common/mesh1/dx1,dx1q,dx3,dx3q
      common/mesh2/y(0:m2)
      common/mesh3/dy(0:m2),h(m2),hm(m2),hc(m2),hp(m2)
      common/mesh4/dym(m2),dyc(m2),dyp(m2)
      common/para/re
      common/size/alx,aly,alz,vol
      common/utau/u_tau_no(m1)

      common/budget1/u(3,0:m1,0:m2),p(0:m1,0:m2),uu(6,m1,m2)
     >              ,uuu(6,2,m1,m2)
      common/budget2/upg(6,m1,m2),diss(6,m1,m2),presstrain(6,m1,m2)
      common/budget3/pp(0:m1,0:m2)
      common/budget5/acc_u2(6,m1,m2)
      common/budget6/v1_osciltm(2,m1,m2)      ! overbar{u~_i*u~_j}
     >              ,v2_osciltm_P(6,m1,m2)    ! Production term
     >              ,v2_osciltm_C(6,m1,m2)    ! Convection term

      real produc(m2)
      real dissip(m2)
      real visdif(m2)
      real convec(m2)
      real transp(m2)
      real vpgrad(m2)

      integer moni_x

      character*30 name

      name='budet_uu3x'
      nn=index(name,'x')
      write(unit=name(nn+1:),fmt='(bn,i2.2,a1,i1.1,a4)')
     >           ith,'p',k_phase,'.dat'

      open(10,file=name,status='unknown')
      write(10,11) 
11    format('variables=y,P,`e,D,C,T,`P,P_o_s,C_o_s,sum')

      i=moni_x
      ip=i+1
      im=i-1
     
      do j=2,n2mm
      jp=j+1
      jm=j-1

c     production term
      produc(j)=-2.*(
     >   uu(5,i,j)*(u(3,ip,j)-u(3,im,j))*dx1/2.  
     >  +uu(6,i,j)*  
     >    ((h(j)**2)*u(3,i,jp)
     >   -(h(j)**2-h(jp)**2)*u(3,i,j)
     >   -(h(jp)**2)*u(3,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     >  )                                 ! -2<u_3u_k>u3,k


c     dissipation term
      dissip(j)=-2./re*diss(3,i,j)

c     viscous diffusion term
      visdif(j)=1./re*(
     >   (uu(3,ip,j)-2.*uu(3,i,j)+uu(3,im,j))*dx1q
     >  +(hp(j)*uu(3,i,jp)
     >   -hc(j)*uu(3,i,j )
     >   +hm(j)*uu(3,i,jm))
     > ) 
     > 

c     convection term      
      convec(j)=-(
     >    u(1,i,j)*(uu(3,ip,j)-uu(3,im,j))*dx1/2.
     > +  u(2,i,j)*
     >    ((h(j)**2)*uu(3,i,jp)
     >   -(h(j)**2-h(jp)**2)*uu(3,i,j)
     >   -(h(jp)**2)*uu(3,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     >  )                                 ! -u_k*<u_3u_3>,k

c     turbulent transport term
      transp(j)=-(
     >    (uuu(3,1,ip,j)-uuu(3,1,im,j))*dx1/2.
     > +  
     >    ((h(j)**2)*uuu(3,2,i,jp)
     >   -(h(j)**2-h(jp)**2)*uuu(3,2,i,j)
     >   -(h(jp)**2)*uuu(3,2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     > )

c     velocity gradient term
      vpgrad(j)=-upg(3,i,j)

c---
      sum=produc(j)+dissip(j)+visdif(j)+convec(j)
     >   +transp(j)+vpgrad(j)
     >   -acc_u2(3,i,j)
     >   +v2_osciltm_P(3,i,j)
     >   +v2_osciltm_C(3,i,j)

      scale=1./(u_tau_no(i)**4*RE)

      yy=0.5*(y(j)+y(j+1))

      write(10,100) 
     >      yy*u_tau_no(i)*RE
     >     ,(produc(j)+v2_osciltm_P(3,i,j))*scale
     >     ,dissip(j)*scale
     >     ,visdif(j)*scale
     >     ,(convec(j)-acc_u2(3,i,j)+v2_osciltm_C(3,i,j))*scale
     >     ,transp(j)*scale
     >     ,vpgrad(j)*scale
     >     ,v2_osciltm_P(3,i,j)*scale
     >     ,v2_osciltm_C(3,i,j)*scale
     >     ,sum*scale

100   format(10(e12.5,x))    
c---

      enddo
      close(10)

      return
      end


c*************** uu4budget *********************************************
      subroutine uu4budget(moni_x)
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/fileio/ith,k_phase
      common/indx1/ipa(m1),imu(m1),imv(m1),kpa(m3),kma(m3)
      common/indx2/jpa(m2),jmu(m2),jmv(m2)
      common/mesh1/dx1,dx1q,dx3,dx3q
      common/mesh2/y(0:m2)
      common/mesh3/dy(0:m2),h(m2),hm(m2),hc(m2),hp(m2)
      common/mesh4/dym(m2),dyc(m2),dyp(m2)
      common/para/re
      common/size/alx,aly,alz,vol
      common/utau/u_tau_no(m1)

      common/budget1/u(3,0:m1,0:m2),p(0:m1,0:m2),uu(6,m1,m2)
     >              ,uuu(6,2,m1,m2)
      common/budget2/upg(6,m1,m2),diss(6,m1,m2),presstrain(6,m1,m2)
      common/budget3/pp(0:m1,0:m2)
      common/budget5/acc_u2(6,m1,m2)
      common/budget6/v1_osciltm(2,m1,m2)      ! overbar{u~_i*u~_j}
     >              ,v2_osciltm_P(6,m1,m2)    ! Production term
     >              ,v2_osciltm_C(6,m1,m2)    ! Convection term
      common/budget7/rss_os(2,m1,m2)

      real produc(m2)
      real dissip(m2)
      real visdif(m2)
      real convec(m2)
      real transp(m2)
      real vpgrad(m2)

      integer moni_x

      character*30 name

      name='budet_uu4x'
      nn=index(name,'x')
      write(unit=name(nn+1:),fmt='(bn,i2.2,a1,i1.1,a4)')
     >           ith,'p',k_phase,'.dat'

      open(10,file=name,status='unknown')
      write(10,11) 
11    format('variables=y,P,`e,D,C,T,`P,P_o_s,C_o_s,sum')


      name='rss_x'
      nn=index(name,'x')
      write(unit=name(nn+1:),fmt='(bn,i2.2,a1,i1.1,a4)')
     >           ith,'p',k_phase,'.dat'

      open(11,file=name,status='unknown')
      write(11,12) 
12    format('variables=y^+_o,P,Pos1,Pos2,v^2,dUdy,Py,Px')

      name='rss4x'
      nn=index(name,'x')
      write(unit=name(nn+1:),fmt='(bn,i2.2,a1,i1.1,a4)')
     >           ith,'p',k_phase,'.dat'

      open(12,file=name,status='unknown')
      write(12,13) 
13    format('variables=y^+_o,P,P_o_s,Px,Py')

      i=moni_x

      ip=i+1
      im=i-1
     
      do j=2,n2mm
      jp=j+1
      jm=j-1

c     production term
      produc(j)=-(
     >   uu(1,i,j)*(u(2,ip,j)-u(2,im,j))*dx1/2.  
     >  +uu(4,i,j)*  
     >    ((h(j)**2)*u(2,i,jp)
     >   -(h(j)**2-h(jp)**2)*u(2,i,j)
     >   -(h(jp)**2)*u(2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))              !<u_1u_k>u2,k
     >  +uu(4,i,j)*(u(1,ip,j)-u(1,im,j))*dx1/2.  
     >  +uu(2,i,j)*  
     >    ((h(j)**2)*u(1,i,jp)
     >   -(h(j)**2-h(jp)**2)*u(1,i,j)
     >   -(h(jp)**2)*u(1,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))              !<u_2u_k>u1,k
     >  )                                 ! -(<u_1u_k>u2,k+<u_2u_k>u1,k)

c     dissipation term
      dissip(j)=-2./re*diss(4,i,j)

c     viscous diffusion term
      visdif(j)=1./re*(
     >   (uu(4,ip,j)-2.*uu(4,i,j)+uu(4,im,j))*dx1q
     >  +(hp(j)*uu(4,i,jp)
     >   -hc(j)*uu(4,i,j )
     >   +hm(j)*uu(4,i,jm))
     > ) 
     > 

c     convection term      
      convec(j)=-(
     >    u(1,i,j)*(uu(4,ip,j)-uu(4,im,j))*dx1/2.
     > +  u(2,i,j)*
     >    ((h(j)**2)*uu(4,i,jp)
     >   -(h(j)**2-h(jp)**2)*uu(4,i,j)
     >   -(h(jp)**2)*uu(4,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     >  )                                 ! -u_k*<u_1u_2>,k

c     turbulent transport term
      transp(j)=-(
     >    (uuu(4,1,ip,j)-uuu(4,1,im,j))*dx1/2.
     > +  
     >    ((h(j)**2)*uuu(4,2,i,jp)
     >   -(h(j)**2-h(jp)**2)*uuu(4,2,i,j)
     >   -(h(jp)**2)*uuu(4,2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     > )

c     velocity gradient term
      vpgrad(j)=-upg(4,i,j)

c---
      sum=produc(j)+dissip(j)+visdif(j)+convec(j)
     >   +transp(j)+vpgrad(j)
     >   -acc_u2(4,i,j)
     >   +v2_osciltm_P(4,i,j)
     >   +v2_osciltm_C(4,i,j)

      scale=1./(u_tau_no(i)**4*RE)

      yy=0.5*(y(j)+y(j+1))

      write(10,100) 
     >      yy*u_tau_no(i)*RE
     >     ,-(produc(j)+v2_osciltm_P(4,i,j))*scale
     >     ,-dissip(j)*scale
     >     ,-visdif(j)*scale
     >     ,-(convec(j)-acc_u2(4,i,j)+v2_osciltm_C(4,i,j))*scale
     >     ,-transp(j)*scale
     >     ,-vpgrad(j)*scale
     >     ,-v2_osciltm_P(4,i,j)*scale
     >     ,-v2_osciltm_C(4,i,j)*scale
     >     ,-sum*scale

100   format(10(e12.5,x))    
c---

      rss1=uu(2,i,j)
      rss2= 
     >    ( (h(j)**2)*u(1,i,jp)
     >     -(h(j)**2-h(jp)**2)*u(1,i,j)
     >     -(h(jp)**2)*u(1,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))           

      write(11,101) 
     >      yy*u_tau_no(i)*RE
     >     ,-(produc(j)+v2_osciltm_P(4,i,j))*scale
     >     ,-rss_os(1,i,j)*scale
     >     ,-rss_os(2,i,j)*scale
     >     ,rss1/u_tau_no(i)**2
     >     ,rss2/(u_tau_no(i)**2*RE)
     >     ,+rss1*rss2*scale
     >     ,+uu(1,i,j)*(u(2,ip,j)-u(2,im,j))*dx1/2.*scale  

101   format(8(e12.5,x))    

      dvdx = (u(2,ip,j)-u(2,im,j))*dx1/2.
      dudy =( (h(j)**2)*u(1,i,jp)
     >     -(h(j)**2-h(jp)**2)*u(1,i,j)
     >     -(h(jp)**2)*u(1,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))


      write(12,102)
     >      yy*u_tau_no(i)*RE
     >     ,-(produc(j)+v2_osciltm_P(4,i,j))*scale
     >     ,-rss_os(1,i,j)*scale
     >      -rss_os(2,i,j)*scale
     >     , 1.*uu(1,i,j)*dvdx*scale
     >     , 1.*uu(2,i,j)*dudy*scale

102   format(5(e12.5,x))
c---



      enddo
      close(10)
      close(11)
      close(12)

      return
      end


c*************** k_budget **********************************************
      subroutine k_budget(moni_x)
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/fileio/ith,k_phase
      common/indx1/ipa(m1),imu(m1),imv(m1),kpa(m3),kma(m3)
      common/indx2/jpa(m2),jmu(m2),jmv(m2)
      common/mesh1/dx1,dx1q,dx3,dx3q
      common/mesh2/y(0:m2)
      common/mesh3/dy(0:m2),h(m2),hm(m2),hc(m2),hp(m2)
      common/mesh4/dym(m2),dyc(m2),dyp(m2)
      common/para/re
      common/size/alx,aly,alz,vol

      common/budget1/u(3,0:m1,0:m2),p(0:m1,0:m2),uu(6,m1,m2)
     >              ,uuu(6,2,m1,m2)
      common/budget2/upg(6,m1,m2),diss(6,m1,m2),presstrain(6,m1,m2)
      common/budget3/pp(0:m1,0:m2)
      common/budget5/acc_u2(6,m1,m2)
      common/budget6/v1_osciltm(2,m1,m2)      ! overbar{u~_i*u~_j}
     >              ,v2_osciltm_P(6,m1,m2)    ! Production term
     >              ,v2_osciltm_C(6,m1,m2)    ! Convection term

      real produc(m2)
      real dissip(m2)
      real visdif(m2)
      real convec(m2)
      real transp(m2)
      real vpgrad(m2)


      integer moni_x

      character*30 name

      name='budet_kx'
      nn=index(name,'x')
      write(unit=name(nn+1:),fmt='(bn,i2.2,a1,i1.1,a4)')
     >           ith,'p',k_phase,'.dat'

      open(10,file=name,status='unknown')
      write(10,11) 
11    format('variables=y,P,Acc,`e,VD,C,T,`p,sum')

      i=moni_x

      ip=i+1
      im=i-1
     
      do j=2,n2mm
      jp=j+1
      jm=j-1

c     production term
      produc(j)=-(
     >   uu(1,i,j)*(u(1,ip,j)-u(1,im,j))*dx1/2.  
     >  +uu(4,i,j)*  
     >    ((h(j)**2)*u(1,i,jp)
     >   -(h(j)**2-h(jp)**2)*u(1,i,j)
     >   -(h(jp)**2)*u(1,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     >+  uu(4,i,j)*(u(2,ip,j)-u(2,im,j))*dx1/2.  
     >  +uu(2,i,j)*  
     >    ((h(j)**2)*u(2,i,jp)
     >   -(h(j)**2-h(jp)**2)*u(2,i,j)
     >   -(h(jp)**2)*u(2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     >+  uu(5,i,j)*(u(3,ip,j)-u(3,im,j))*dx1/2.  
     >  +uu(6,i,j)*  
     >    ((h(j)**2)*u(3,i,jp)
     >   -(h(j)**2-h(jp)**2)*u(3,i,j)
     >   -(h(jp)**2)*u(3,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     >  )

c     dissipation term
      dissip(j)=-1./re
     >       *(diss(1,i,j)+diss(2,i,j)+diss(3,i,j))

c     viscous diffusion term
      visdif(j)=0.5/re*(
     >   (uu(1,ip,j)-2.*uu(1,i,j)+uu(1,im,j))*dx1q
     >  +(hp(j)*uu(1,i,jp)
     >   -hc(j)*uu(1,i,j )
     >   +hm(j)*uu(1,i,jm))
     >  +(uu(2,ip,j)-2.*uu(2,i,j)+uu(2,im,j))*dx1q
     >  +(hp(j)*uu(2,i,jp)
     >   -hc(j)*uu(2,i,j )
     >   +hm(j)*uu(2,i,jm))
     >  +(uu(3,ip,j)-2.*uu(3,i,j)+uu(3,im,j))*dx1q
     >  +(hp(j)*uu(3,i,jp)
     >   -hc(j)*uu(3,i,j )
     >   +hm(j)*uu(3,i,jm))
     > ) 
     > 

c     convection term      
      convec(j)=-0.5*(
     >    u(1,i,j)*(uu(1,ip,j)-uu(1,im,j))*dx1/2.
     > +  u(1,i,j)*(uu(2,ip,j)-uu(2,im,j))*dx1/2.
     > +  u(1,i,j)*(uu(3,ip,j)-uu(3,im,j))*dx1/2.
     > +  u(2,i,j)*
     >    ((h(j)**2)*uu(1,i,jp)
     >   -(h(j)**2-h(jp)**2)*uu(1,i,j)
     >   -(h(jp)**2)*uu(1,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     > +  u(2,i,j)*
     >    ((h(j)**2)*uu(2,i,jp)
     >   -(h(j)**2-h(jp)**2)*uu(2,i,j)
     >   -(h(jp)**2)*uu(2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     > +  u(2,i,j)*
     >    ((h(j)**2)*uu(3,i,jp)
     >   -(h(j)**2-h(jp)**2)*uu(3,i,j)
     >   -(h(jp)**2)*uu(3,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     >  )

c     turbulent transport term
      transp(j)=-0.5*(
     >    (uuu(1,1,ip,j)-uuu(1,1,im,j))*dx1/2.
     > +  
     >    ((h(j)**2)*uuu(1,2,i,jp)
     >   -(h(j)**2-h(jp)**2)*uuu(1,2,i,j)
     >   -(h(jp)**2)*uuu(1,2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     > +  (uuu(2,1,ip,j)-uuu(2,1,im,j))*dx1/2.
     > +  
     >    ((h(j)**2)*uuu(2,2,i,jp)
     >   -(h(j)**2-h(jp)**2)*uuu(2,2,i,j)
     >   -(h(jp)**2)*uuu(2,2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     > +  (uuu(3,1,ip,j)-uuu(3,1,im,j))*dx1/2.
     > +  
     >    ((h(j)**2)*uuu(3,2,i,jp)
     >   -(h(j)**2-h(jp)**2)*uuu(3,2,i,j)
     >   -(h(jp)**2)*uuu(3,2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     > )

c     velocity gradient term
      vpgrad(j)=-0.5*(upg(1,i,j)+upg(2,i,j)+upg(3,i,j))

      sum=produc(j)+dissip(j)+visdif(j)+convec(j)
     >   +transp(j)+vpgrad(j)
     >   -(acc_u2(1,i,j)+acc_u2(2,i,j)+acc_u2(3,i,j))/2.


      write(10,100) y(j)
     >             ,produc(j)
     >   ,(acc_u2(1,i,j)+acc_u2(2,i,j)+acc_u2(3,i,j))/2.
     >             ,dissip(j)
     >             ,visdif(j)
     >             ,convec(j)
     >   -(acc_u2(1,i,j)+acc_u2(2,i,j)+acc_u2(3,i,j))/2.
     >             ,transp(j)
     >             ,vpgrad(j)
     >             ,sum

100   format(9(e12.5,x))    
      enddo
      close(10)

      return
      end

c*************** press_strain **********************************************
      subroutine press_strain(moni_x)
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/fileio/ith,k_phase
      common/indx1/ipa(m1),imu(m1),imv(m1),kpa(m3),kma(m3)
      common/indx2/jpa(m2),jmu(m2),jmv(m2)
      common/mesh1/dx1,dx1q,dx3,dx3q
      common/mesh2/y(0:m2)
      common/mesh3/dy(0:m2),h(m2),hm(m2),hc(m2),hp(m2)
      common/mesh4/dym(m2),dyc(m2),dyp(m2)
      common/para/re
      common/size/alx,aly,alz,vol
      common/utau/u_tau_no(m1)

      common/budget1/u(3,0:m1,0:m2),p(0:m1,0:m2),uu(6,m1,m2)
     >              ,uuu(6,2,m1,m2)
      common/budget2/upg(6,m1,m2),diss(6,m1,m2),presstrain(6,m1,m2)
      common/budget3/pp(0:m1,0:m2)
      common/budget5/acc_u2(6,m1,m2)


      integer moni_x

      character*30 name

      name='ps_x'
      nn=index(name,'x')
      write(unit=name(nn+1:),fmt='(bn,i2.2,a1,i1.1,a4)')
     >           ith,'p',k_phase,'.dat'

      open(10,file=name,status='unknown')
      write(10,11) 
11    format('variables=y,`F_1_1,`F_2_2,`F_3_3')

      i=moni_x
      do 10 j=1,n2m
      yy=0.5*(y(j)+y(j+1))
      write(10,100) yy*u_tau_no(i)*re
     >             ,presstrain(1,i,j)/(u_tau_no(i)**4*re)
     >             ,presstrain(2,i,j)/(u_tau_no(i)**4*re)
     >             ,presstrain(3,i,j)/(u_tau_no(i)**4*re)
10    continue
100   format(4(e12.5,x))
      close(10)
   
      return
      end


c*************** read_avg **********************************************
      subroutine read_avg(v1m,v2m,v3m,v4m, 
     >                    p1m,p2m,p3m,p4m,
     >                    vor1m,vor2m, 
     >                    trans,dissp,phi,redis,
     >                    iphase,nphavg)

      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm

      real v1m(3,m1,m2),v2m(6,m1,m2),v3m(3,m1,m2),v4m(3,m1,m2)
      real p1m(m1,m2),p2m(m1,m2),p3m(m1,m2),p4m(m1,m2)
      real vor1m(3,m1,m2),vor2m(3,m1,m2)
      real trans(6,2,m1,m2),dissp(6,m1,m2),phi(6,m1,m2),redis(6,m1,m2)

      integer iphase,nphavg

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
          read(11) (V1M(L,I,J),L=1,3),P1M(I,J),P2M(I,J)
          read(11) (V2M(L,I,J),L=1,6)
          read(11) (V3M(L,I,J),L=1,3)
          read(11) (V4M(L,I,J),L=1,3)
          read(11) (TRANS(L,1,I,J),L=1,6),
     >             (TRANS(L,2,I,J),L=1,6)
          read(11) (DISSP(L,I,J),L=1,6)
          read(11) (PHI(L,I,J),L=1,6)
          read(11) (VOR1M(L,I,J),L=1,3),(VOR2M(L,I,J),L=1,3)
          read(11) (REDIS(L,I,J),L=1,6)
        enddo
      enddo

      close(11)

      write(*,*) 'iphase=',iphase,'nphavg=',nphavg

      return
      end

c*************** acc_term ********************************************
      subroutine acc_term 
     
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/mesh1/dx1,dx1q,dx3,dx3q
      common/mesh2/y(0:m2)
      common/mesh3/dy(0:m2),h(m2),hm(m2),hc(m2),hp(m2)
      common/mesh4/dym(m2),dyc(m2),dyp(m2)

      common/moni2/moni_ph
      COMMON/TSTEP/NTST,DT,CFLMAX,NTIME


      double precision um(3,m1,m2)
      double precision pm(m1,m2)
      double precision uum(6,m1,m2)
      double precision uuum(6,2,m1,m2)
      double precision upgm(6,m1,m2),redism(6,m1,m2)
      double precision dissm(6,m1,m2)
      double precision ppm(m1,m2)
      double precision gp(3)

      double precision v3m(3,m1,m2),v4m(3,m1,m2)
      double precision p3m(m1,m2),p4m(m1,m2)
      double precision vor1m(3,m1,m2),vor2m(3,m1,m2)

      real u_np(3,m1,m2),uu_np(6,m1,m2)  ! at n+1 phase
      real u_nm(3,m1,m2),uu_nm(6,m1,m2)  ! at n-1 phase

      real acc_u1(3,m1,m2),acc_u2(6,m1,m2)


c---- read the n+1 phase data file

      iphase = moni_ph + 1
      if (iphase .eq. mt+1) iphase = 0


         call read_avg(um,uum,v3m,v4m, 
     >                 pm,ppm,p3m,p4m,
     >                 vor1m,vor2m, 
     >                 uuum,dissm,upgm,redism,
     >                 iphase,nphavg)


      do 3 i=1,n1m
      do 3 j=1,n2m

         do l=1,3
         u_np(l,i,j)=um(l,i,j)
         enddo

         do l=1,6
         if(l.eq.1) then
            nv1=1
            nv2=1
         endif
         if(l.eq.2) then
            nv1=2
            nv2=2
         endif
         if(l.eq.3) then
            nv1=3
            nv2=3
         endif
         if(l.eq.4) then
            nv1=1
            nv2=2
         endif
         if(l.eq.5) then
            nv1=1
            nv2=3
         endif
         if(l.eq.6) then
            nv1=2
            nv2=3
         endif

         uu_np(l,i,j)=uum(l,i,j)-u_np(nv1,i,j)*u_np(nv2,i,j)

         enddo 
3     continue




c---- read the n-1 phase data file

      iphase = moni_ph - 1
      if (iphase .eq. -1) iphase = mt


         call read_avg(um,uum,v3m,v4m, 
     >                 pm,ppm,p3m,p4m,
     >                 vor1m,vor2m, 
     >                 uuum,dissm,upgm,redism,
     >                 iphase,nphavg)


      do 4 i=1,n1m
      do 4 j=1,n2m

         do l=1,3
         u_nm(l,i,j)=um(l,i,j)
         enddo

         do l=1,6
         if(l.eq.1) then
            nv1=1
            nv2=1
         endif
         if(l.eq.2) then
            nv1=2
            nv2=2
         endif
         if(l.eq.3) then
            nv1=3
            nv2=3
         endif
         if(l.eq.4) then
            nv1=1
            nv2=2
         endif
         if(l.eq.5) then
            nv1=1
            nv2=3
         endif
         if(l.eq.6) then
            nv1=2
            nv2=3
         endif

         uu_nm(l,i,j)=uum(l,i,j)-u_nm(nv1,i,j)*u_nm(nv2,i,j)

         enddo 
4     continue




c---- calculate the derivatives


      do 5 i=1,n1m
      do 5 j=1,n2m

         do l=1,3
            acc_u1(l,i,j)=(u_np(l,i,j)-u_nm(l,i,j))/(2.*dt)
         enddo

         do l=1,6
            acc_u2(l,i,j)=(uu_np(l,i,j)-uu_nm(l,i,j))/(2.*dt)
         enddo

5     continue


c---- write the derivativesi to 'acc_term.avg'

      open(10,file='acc_term.avg',form='unformatted',status='unknown')

      do 6 i=2,n1m-1
      do 6 j=2,n2m-1
         write(10) (acc_u1(l,i,j),l=1,2)
         write(10) (acc_u2(l,i,j),l=1,6)
6     continue

      close(10)


      return
      end



c****************budget_term********************************************
      subroutine budget_term 

      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/indx1/ipa(m1),imu(m1),imv(m1),kpa(m3),kma(m3)
      common/indx2/jpa(m2),jmu(m2),jmv(m2)
      common/mesh1/dx1,dx1q,dx3,dx3q
      common/mesh2/y(0:m2)
      common/mesh3/dy(0:m2),h(m2),hm(m2),hc(m2),hp(m2)
      common/mesh4/dym(m2),dyc(m2),dyp(m2)
      common/moni2/moni_ph

      double precision um(3,m1,m2)
      double precision pm(m1,m2)
      double precision uum(6,m1,m2)
      double precision uuum(6,2,m1,m2)
      double precision upgm(6,m1,m2)
      double precision redism(6,m1,m2)
      double precision dissm(6,m1,m2)
      double precision ppm(m1,m2)
      double precision gp(3)
      double precision strain(3,3)

c---- required dimensions to call read_avg

      double precision v3m(3,m1,m2),v4m(3,m1,m2)
      double precision p3m(m1,m2),p4m(m1,m2)
      double precision vor1m(3,m1,m2),vor2m(3,m1,m2)

      real u(3,0:m1,0:m2),p(0:m1,0:m2),uu(6,m1,m2)
     >    ,uuu(6,2,m1,m2)
      real upg(6,m1,m2),diss(6,m1,m2),presstrain(6,m1,m2)
      real pp(0:m1,0:m2)

      real u_t(3,0:m1,0:m2),p_t(0:m1,0:m2),uu_t(6,m1,m2)
     >    ,uuu_t(6,2,m1,m2)
      real upg_t(6,m1,m2),diss_t(6,m1,m2),presstrain_t(6,m1,m2)
      real pp_t(0:m1,0:m2)

      real nsample

      nsample=0.0

      do 2 iphase=0,mt ! for time mean budget term
                         ! Put command for phase averaged budget

c      iphase=moni_ph

         call read_avg(um,uum,v3m,v4m, 
     >                 pm,ppm,p3m,p4m,
     >                 vor1m,vor2m, 
     >                 uuum,dissm,upgm,redism,
     >                 iphase,nphavg)

         phavg   = real(nphavg)
         nsample = nsample + phavg

      do 3 i=1,n1m
      do 3 j=1,n2m

         do l=1,3
         u(l,i,j)=um(l,i,j)
         enddo
         p(i,j)=pm(i,j)
         pp(i,j)=ppm(i,j)-p(i,j)**2  

         do l=1,6
         if(l.eq.1) then
            nv1=1
            nv2=1
         endif
         if(l.eq.2) then
            nv1=2
            nv2=2
         endif
         if(l.eq.3) then
            nv1=3
            nv2=3
         endif
         if(l.eq.4) then
            nv1=1
            nv2=2
         endif
         if(l.eq.5) then
            nv1=1
            nv2=3
         endif
         if(l.eq.6) then
            nv1=2
            nv2=3
         endif

         uu(l,i,j)=uum(l,i,j)-u(nv1,i,j)*u(nv2,i,j)

         enddo 
3     continue


      do 4 i=1,n1m 
         ip=i+1
         im=i-1
         ium=i-imv(i)
         iup=ipa(i)-i 
      do 4 j=1,n2m 
         jp=j+1
         jm=j-1
         jum=j-jmu(j)
         jup=jpa(j)-j

         do l=1,6
         if(l.eq.1) then
            nv1=1
            nv2=1
         endif
         if(l.eq.2) then
            nv1=2
            nv2=2
         endif
         if(l.eq.3) then
            nv1=3
            nv2=3
         endif
         if(l.eq.4) then
            nv1=1
            nv2=2
         endif
         if(l.eq.5) then
            nv1=1
            nv2=3
         endif
         if(l.eq.6) then
            nv1=2
            nv2=3
         endif

         if (nv2.eq.1) l1=1    ! u'(nv2)*u'(1)*u(nv1) 
         if (nv1.eq.1) l2=1    ! u'(nv1)*u'(1)*u(nv2) 
         if (nv2.eq.2) l1=4 
         if (nv1.eq.2) l2=4 
         if (nv2.eq.3) l1=5 
         if (nv1.eq.3) l2=5 

         uuu(l,1,i,j)=uuum(l,1,i,j)
     >               -u(nv1,i,j)*u(nv2,i,j)*u(1,i,j)
     >               -uu(l,i,j)*u(1,i,j)
     >               -uu(l1,i,j)*u(nv1,i,j)
     >               -uu(l2,i,j)*u(nv2,i,j)  

         if (nv2.eq.1) l1=4    ! u'(nv2)*u'(2)*u(nv2) 
         if (nv1.eq.1) l2=4    ! u'(nv1)*u'(2)*u(nv2) 
         if (nv2.eq.2) l1=2 
         if (nv1.eq.2) l2=2 
         if (nv2.eq.3) l1=6 
         if (nv1.eq.3) l2=6 

         uuu(l,2,i,j)=uuum(l,2,i,j)
     >               -u(nv1,i,j)*u(nv2,i,j)*u(2,i,j)
     >               -uu(l,i,j)*u(2,i,j)
     >               -uu(l1,i,j)*u(nv1,i,j)
     >               -uu(l2,i,j)*u(nv2,i,j)

         eps1=
     >    ium*iup*(u(nv1,ip,j)-u(nv1,im,j))*dx1/2
     >           *(u(nv2,ip,j)-u(nv2,im,j))*dx1/2
     >   +(1-ium)*(u(nv1,ip,j)-u(nv1,i ,j))*dx1
     >           *(u(nv2,ip,j)-u(nv2,i ,j))*dx1
     >   +(1-iup)*(u(nv1,i ,j)-u(nv1,im,j))*dx1
     >           *(u(nv2,i ,j)-u(nv2,im,j))*dx1
         eps2=
     >    ( jum*jup
     >     *( h(j)**2*u(nv1,i,jp)
     >      -(h(j)**2-h(jp)**2)*u(nv1,i,j)
     >      - h(jp)**2*u(nv1,i,jm))
     >     /h(j)/h(jp)/(h(j)+h(jp))
     >     +(1-jum)*(u(nv1,i,jp)-u(nv1,i,j ))/h(jp)
     >     +(1-jup)*(u(nv1,i,j )-u(nv1,i,jm))/h(j)
     >    )
     >   *( jum*jup
     >     *( h(j)**2*u(nv2,i,jp)
     >      -(h(j)**2-h(jp)**2)*u(nv2,i,j)
     >      - h(jp)**2*u(nv2,i,jm))
     >     /h(j)/h(jp)/(h(j)+h(jp))
     >     +(1-jum)*(u(nv2,i,jp)-u(nv2,i,j ))/h(jp)
     >     +(1-jup)*(u(nv2,i,j )-u(nv2,i,jm))/h(j)
     >    )

         diss(l,i,j)=dissm(l,i,j)-eps1-eps2

         gp(1)=ium*iup*(p(ip,j)-p(im,j))*dx1/2
     >        +(1-ium)*(p(ip,j)-p(i ,j))*dx1
     >        +(1-iup)*(p(i ,j)-p(im,j))*dx1
         gp(2)=jum*jup
     >        *(h(j)**2*p(i,jp)
     >        -(h(j)**2-h(jp)**2)*p(i,j)
     >        - h(jp)**2*p(i,jm))
     >        /h(j)/h(jp)/(h(j)+h(jp))
     >      +(1-jum)*(p(i,jp)-p(i,j ))/h(jp)
     >      +(1-jup)*(p(i,j )-p(i,jm))/h(j)
         gp(3)=0.0   

         upg(l,i,j)=upgm(l,i,j)
     >             -u(nv1,i,j)*gp(nv2) 
     >             -u(nv2,i,j)*gp(nv1) 


c        pressure strain correlation

         strain(1,1)=ium*iup*(u(1,ip,j)-u(1,im,j))*dx1/2  ! recall u(1,i,j) is cell center value
     >              +(1-ium)*(u(1,ip,j)-u(1,i ,j))*dx1
     >              +(1-iup)*(u(1,i ,j)-u(1,im,j))*dx1

         strain(2,1)=ium*iup*(u(2,ip,j)-u(2,im,j))*dx1/2
     >              +(1-ium)*(u(2,ip,j)-u(2,i ,j))*dx1
     >              +(1-iup)*(u(2,i ,j)-u(2,im,j))*dx1

         strain(1,2)=jum*jup
     >              *(h(j)**2*u(1,i,jp)
     >              -(h(j)**2-h(jp)**2)*u(1,i,j)
     >              - h(jp)**2*u(1,i,jm))
     >              /h(j)/h(jp)/(h(j)+h(jp))
     >            +(1-jum)*(u(1,i,jp)-u(1,i,j ))/h(jp)
     >            +(1-jup)*(u(1,i,j )-u(1,i,jm))/h(j)

         strain(2,2)=jum*jup
     >              *(h(j)**2*u(2,i,jp)
     >              -(h(j)**2-h(jp)**2)*u(2,i,j)
     >              - h(jp)**2*u(2,i,jm))
     >              /h(j)/h(jp)/(h(j)+h(jp))
     >            +(1-jum)*(u(2,i,jp)-u(2,i,j ))/h(jp)
     >            +(1-jup)*(u(2,i,j )-u(2,i,jm))/h(j)

         strain(3,1)=0.0
         strain(3,2)=0.0
         strain(3,3)=0.0
         strain(1,3)=0.0
         strain(2,3)=0.0

         presstrain(l,i,j) = redism(l,i,j)
     >                     - p(i,j)*(strain(nv1,nv2)+strain(nv2,nv1))


         enddo

4     continue

      do 5 i=1,n1m
      do 5 j=1,n2m

      do l=1,3

         u_t(l,i,j)=(u_t(l,i,j)*(nsample-phavg)  
     >              +u(l,i,j)*phavg)/nsample

      enddo

         p_t(i,j)=(p_t(i,j)*(nsample-phavg)  
     >            +p(i,j)*phavg)/nsample

         pp_t(i,j)=(pp_t(i,j)*(nsample-phavg)  
     >             +pp(i,j)*phavg)/nsample

      do l=1,6

         uu_t(l,i,j)=(uu_t(l,i,j)*(nsample-phavg)  
     >               +uu(l,i,j)*phavg)/nsample

         uuu_t(l,1,i,j)=(uuu_t(l,1,i,j)*(nsample-phavg)  
     >                  +uuu(l,1,i,j)*phavg)/nsample

         uuu_t(l,2,i,j)=(uuu_t(l,2,i,j)*(nsample-phavg)  
     >                  +uuu(l,2,i,j)*phavg)/nsample

         upg_t(l,i,j)=(upg_t(l,i,j)*(nsample-phavg)  
     >                +upg(l,i,j)*phavg)/nsample

         diss_t(l,i,j)=(diss_t(l,i,j)*(nsample-phavg)  
     >                 +diss(l,i,j)*phavg)/nsample

         presstrain_t(l,i,j)=(presstrain_t(l,i,j)*(nsample-phavg)  
     >                       +presstrain(l,i,j)*phavg)/nsample

      enddo
5     continue


2     continue


C     Write the time averaged data

      open(10,file='budget.avg',form='unformatted',status='unknown')

      do 9 i=1,n1m
      do 9 j=1,n2m

      write(10) (u_t(l,i,j),l=1,3)
      write(10) (uu_t(l,i,j),l=1,6)
      write(10) p_t(i,j),pp_t(i,j)
      write(10) ((uuu_t(l,m,i,j),m=1,2),l=1,6) 
      write(10) (upg_t(l,i,j),l=1,6)
      write(10) (diss_t(l,i,j),l=1,6)
      write(10) (presstrain_t(l,i,j),l=1,6)

9     continue

      close(10)

      return
      end

c*************** oscil_term ********************************************
      subroutine oscil_term 
     
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/mesh1/dx1,dx1q,dx3,dx3q
      common/mesh2/y(0:m2)
      common/mesh3/dy(0:m2),h(m2),hm(m2),hc(m2),hp(m2)
      common/mesh4/dym(m2),dyc(m2),dyp(m2)

      real v1tm(3,m1,m2)
      real v2tm(6,m1,m2)
      real v1_osciltm(2,m1,m2)      ! overbar{u~_i*u~_j}
      real v2_osciltm_P(6,m1,m2)    ! Production term
      real v2_osciltm_C(6,m1,m2)    ! Convection term
      real rss_os(2,m1,m2)          ! each comp of OS term of RSS production

      real v1pm(3,m1,m2),v2pm(6,m1,m2),v3pm(3,m1,m2),v4pm(3,m1,m2)
      real p1pm(m1,m2),p2pm(m1,m2),p3pm(m1,m2),p4pm(m1,m2)
      real vor1pm(3,m1,m2),vor2pm(3,m1,m2)
      real trans_p(6,2,m1,m2),dissp_p(6,m1,m2),phi_p(6,m1,m2)
      real redis_p(6,m1,m2)

      real dummy,dummy0(6,2,m1,m2),dummy1(m1,m2),dummy2(6,m1,m2)
 
      integer ll(3,3)

      real nsample

      nsample=0.0

C     read the time averaged data 

      open(10,file='budget.avg',form='unformatted',status='unknown') 
 
      do 2 i=1,n1m
      do 2 j=1,n2m

      read(10) (v1tm(l,i,j),l=1,3)
      read(10) (v2tm(l,i,j),l=1,6)
      read(10) dummy1(i,j),dummy1(i,j)
      read(10) ((dummy0(l,m,i,j),m=1,2),l=1,6)
      read(10) (dummy2(l,i,j),l=1,6)
      read(10) (dummy2(l,i,j),l=1,6)
      read(10) (dummy2(l,i,j),l=1,6)

2     continue
 
      close(10)


      do 1 iphase=0,mt

         call read_avg(v1pm,v2pm,v3pm,v4pm, 
     >                 p1pm,p2pm,p3pm,p4pm,
     >                 vor1pm,vor2pm, 
     >                 trans_p,dissp_p,phi_p,redis_p,
     >                 iphase,nphavg)

         phavg   = real(nphavg)
         nsample = nsample + phavg


c------- calculate <u_i"*u_j"> from overbar{u_i*u_j}
c        and then, put this value to v2pm dimension

         do 5 l=1,6
         do 5 i=1,n1m
         do 5 j=1,n2m


               if(l.eq.1) then
                  nv1=1
                  nv2=1
               endif 
               if(l.eq.2) then
                  nv1=2
                  nv2=2
               endif 
               if(l.eq.3) then
                  nv1=3
                  nv2=3
               endif 
               if(l.eq.4) then
                  nv1=1
                  nv2=2
               endif 
               if(l.eq.5) then
                  nv1=2
                  nv2=3
               endif 
               if(l.eq.6) then
                  nv1=1
                  nv2=3
               endif


               v2pm(l,i,j)=v2pm(l,i,j)
     >                    -v1pm(nv1,i,j)*v1pm(nv2,i,j)


5        continue


         do 10 i=2,n1m-1
         do 10 j=2,n2m-1
            
            ip=i+1
            im=i-1
            jp=j+1
            jm=j-1

c------------- oscil term in U_i budget

            do nv=1,2

c------------- tilde{Conv_i}

               v1_oscilpm= 
     >   -1.0*(
     >         (
     >          (v1pm(nv,ip,j)-v1tm(nv,ip,j))          
     >         *(v1pm(1 ,ip,j)-v1tm(1 ,ip,j))          
     >         -(v1pm(nv,im,j)-v1tm(nv,im,j))          
     >         *(v1pm(1 ,im,j)-v1tm(1 ,im,j))          
     >         )*dx1/2.
     >        +
     >         (
     >                      h(j)**2*(v1pm(nv,i,jp)-v1tm(nv,i,jp))
     >                             *(v1pm(2 ,i,jp)-v1tm(2 ,i,jp))          
     >          -(h(j)**2-h(jp)**2)*(v1pm(nv,i,j )-v1tm(nv,i,j ))
     >                             *(v1pm(2 ,i,j )-v1tm(2 ,i,j ))          
     >          -          h(jp)**2*(v1pm(nv,i,jm)-v1tm(nv,i,jm))
     >                             *(v1pm(2 ,i,jm)-v1tm(2 ,i,jm))          
     >         )
     >         /h(j)/h(jp)/(h(j)+h(jp))
     >        )


               v1_osciltm(nv,i,j)=(v1_osciltm(nv,i,j)*(nsample-phavg)  
     >                            +v1_oscilpm*phavg)/nsample

            enddo



c------------- oscil term in u_i*u_j budget

            do l=1,6

               if(l.eq.1) then
                  nv1=1
                  nv2=1
               endif 
               if(l.eq.2) then
                  nv1=2
                  nv2=2
               endif 
               if(l.eq.3) then
                  nv1=3
                  nv2=3
               endif 
               if(l.eq.4) then
                  nv1=1
                  nv2=2
               endif 
               if(l.eq.5) then
                  nv1=2
                  nv2=3
               endif 
               if(l.eq.6) then
                  nv1=1
                  nv2=3
               endif
 

               ll(1,1)=1
               ll(2,2)=2
               ll(3,3)=3
               ll(1,2)=4
               ll(2,1)=4
               ll(1,3)=5
               ll(3,1)=5
               ll(2,3)=6
               ll(3,2)=6

c------------- tilde{P_ij}

               v2_oscilpm_P1=     
     >   -1.0*(
     >         (v2pm(ll(nv1,1),i,j)-v2tm(ll(nv1,1),i,j))*
     >         ( 
     >            (v1pm(nv2,ip,j)-v1tm(nv2,ip,j))
     >          - (v1pm(nv2,im,j)-v1tm(nv2,im,j))
     >         )*dx1/2.
     >        +
     >         (v2pm(ll(nv2,1),i,j)-v2tm(ll(nv2,1),i,j))*
     >         ( 
     >            (v1pm(nv1,ip,j)-v1tm(nv1,ip,j))
     >          - (v1pm(nv1,im,j)-v1tm(nv1,im,j))
     >         )*dx1/2.
     >        )

               v2_oscilpm_P2=     
     >   -1.0*(
     >         (v2pm(ll(nv1,2),i,j)-v2tm(ll(nv1,2),i,j))*
     >         (
     >                      h(j)**2*(v1pm(nv2,i,jp)-v1tm(nv2,i,jp))
     >          -(h(j)**2-h(jp)**2)*(v1pm(nv2,i,j )-v1tm(nv2,i,j ))
     >          -          h(jp)**2*(v1pm(nv2,i,jm)-v1tm(nv2,i,jm))
     >         )
     >         /h(j)/h(jp)/(h(j)+h(jp))
     >        +
     >         (v2pm(ll(nv2,2),i,j)-v2tm(ll(nv2,2),i,j))*
     >         (
     >                      h(j)**2*(v1pm(nv1,i,jp)-v1tm(nv1,i,jp))
     >          -(h(j)**2-h(jp)**2)*(v1pm(nv1,i,j )-v1tm(nv1,i,j ))
     >          -          h(jp)**2*(v1pm(nv1,i,jm)-v1tm(nv1,i,jm))
     >         )
     >         /h(j)/h(jp)/(h(j)+h(jp))
     >        )
               
               v2_oscilpm_P=v2_oscilpm_P1+v2_oscilpm_P2     

               v2_osciltm_P(l,i,j)=(v2_osciltm_P(l,i,j)*(nsample-phavg)  
     >                             +v2_oscilpm_P*phavg)/nsample

c--- for RSS production
               if (l.eq.4) then
               rss_os(1,i,j)=(rss_os(1,i,j)*(nsample-phavg)
     >                       +v2_oscilpm_P1*phavg)/nsample
               rss_os(2,i,j)=(rss_os(2,i,j)*(nsample-phavg)
     >                       +v2_oscilpm_P2*phavg)/nsample
               endif
c---


c------------- tilde{C_ij}

               v2_oscilpm_C=
     >   -1.0*(
     >         (v1pm(1,i,j)-v1tm(1,i,j))*
     >         ( 
     >            (v2pm(l,ip,j)-v2tm(l,ip,j))
     >          - (v2pm(l,im,j)-v2tm(l,im,j))
     >         )*dx1/2.
     >        + 
     >         (v1pm(2,i,j)-v1tm(2,i,j))*
     >         (
     >                      h(j)**2*(v2pm(l,i,jp)-v2tm(l,i,jp))
     >          -(h(j)**2-h(jp)**2)*(v2pm(l,i,j )-v2tm(l,i,j ))
     >          -          h(jp)**2*(v2pm(l,i,jm)-v2tm(l,i,jm))
     >         )
     >         /h(j)/h(jp)/(h(j)+h(jp))
     >        )

               v2_osciltm_C(l,i,j)=(v2_osciltm_C(l,i,j)*(nsample-phavg)  
     >                             +v2_oscilpm_C*phavg)/nsample

            enddo

10       continue
 
1     continue




C     Write the time averaged data 

      open(10,file='oscil_term.avg',form='unformatted',status='unknown') 
      do 3 i=2,n1m-1
      do 3 j=2,n2m-1

         write(10) (v1_osciltm(l,i,j),l=1,2)
         write(10) (v2_osciltm_P(l,i,j),l=1,6)
         write(10) (v2_osciltm_C(l,i,j),l=1,6)
         write(10) (rss_os(l,i,j),l=1,2)

3     continue
      close(10)

      return
      end

c*************** two_d_uu1 *********************************************
      subroutine two_d_uu1
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/fileio/ith,k_phase
      common/indx1/ipa(m1),imu(m1),imv(m1),kpa(m3),kma(m3)
      common/indx2/jpa(m2),jmu(m2),jmv(m2)
      common/mesh1/dx1,dx1q,dx3,dx3q
      common/mesh2/y(0:m2)
      common/mesh3/dy(0:m2),h(m2),hm(m2),hc(m2),hp(m2)
      common/mesh4/dym(m2),dyc(m2),dyp(m2)
      common/para/re
      common/size/alx,aly,alz,vol
      common/utau/u_tau_no(m1)

      common/budget1/u(3,0:m1,0:m2),p(0:m1,0:m2),uu(6,m1,m2)
     >              ,uuu(6,2,m1,m2)
      common/budget2/upg(6,m1,m2),diss(6,m1,m2),presstrain(6,m1,m2)
      common/budget3/pp(0:m1,0:m2)
      common/budget5/acc_u2(6,m1,m2)
      common/budget6/v1_osciltm(2,m1,m2)      ! overbar{u~_i*u~_j}
     >              ,v2_osciltm_P(6,m1,m2)    ! Production term
     >              ,v2_osciltm_C(6,m1,m2)    ! Convection term

      real produc(m2)
      real dissip(m2)
      real visdif(m2)
      real convec(m2)
      real transp(m2)
      real vpgrad(m2)
      

      integer moni_x

      character*30 name

      name='2d_budet_uu1'
      nn=index(name,'1')
      write(unit=name(nn+1:),fmt='(bn,a2,i1.1,a4)')
     >           '_p',k_phase,'.dat'



      open(10,file=name,status='unknown')
      write(10,11) 
11    format('variables=x,y,P,`e,D,C,T,`P,P_o_s,C_o_s,sum')
      write(10,12) n1m-20, n2m-2
12    format('zone i=',i4,', j=',i4,', f=point')
     
      do j=2,n2m-1
      jp=j+1
      jm=j-1

      do i=11,n1m-10
      ip=i+1
      im=i-1

c     production term
      produc(j)=-2.*(
     >   uu(1,i,j)*(u(1,ip,j)-u(1,im,j))*dx1/2.  
     >  +uu(4,i,j)*  
     >    ((h(j)**2)*u(1,i,jp)
     >   -(h(j)**2-h(jp)**2)*u(1,i,j)
     >   -(h(jp)**2)*u(1,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     >  )                                 ! -2<u_1u_k>u1,k

c     dissipation term
      dissip(j)=-2./re*diss(1,i,j)

c     viscous diffusion term
      visdif(j)=1./re*(
     >   (uu(1,ip,j)-2.*uu(1,i,j)+uu(1,im,j))*dx1q
     >  +(hp(j)*uu(1,i,jp)
     >   -hc(j)*uu(1,i,j )
     >   +hm(j)*uu(1,i,jm))
     > ) 
     > 

c     convection term      
      convec(j)=-(
     >    u(1,i,j)*(uu(1,ip,j)-uu(1,im,j))*dx1/2.
     > +  u(2,i,j)*
     >    ((h(j)**2)*uu(1,i,jp)
     >   -(h(j)**2-h(jp)**2)*uu(1,i,j)
     >   -(h(jp)**2)*uu(1,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     >  )                                 ! -u_k*<u_1u_1>,k

c     turbulent transport term
      transp(j)=-(
     >    (uuu(1,1,ip,j)-uuu(1,1,im,j))*dx1/2.
     > +  
     >    ((h(j)**2)*uuu(1,2,i,jp)
     >   -(h(j)**2-h(jp)**2)*uuu(1,2,i,j)
     >   -(h(jp)**2)*uuu(1,2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     > )

c     velocity gradient term
      vpgrad(j)=-upg(1,i,j)


c---
      sum=produc(j)+dissip(j)+visdif(j)+convec(j)
     >   +transp(j)+vpgrad(j)
     >   -acc_u2(1,i,j)
     >   +v2_osciltm_P(1,i,j)
     >   +v2_osciltm_C(1,i,j)

      scale=1./(u_tau_no(i)**4*RE)

      xc=(i-0.5)/dx1
      yc=0.5*(y(j)+y(j+1))

                   x_slot=75.0
                   u_tau_in=0.0492
c                   x_slot=78.90625
c                   u_tau_in=0.0525

         write(10,100)  (xc-x_slot)*u_tau_in*re
     >                , yc*u_tau_no(i)*re
     >     ,(produc(j)+v2_osciltm_P(1,i,j))*scale
     >     ,dissip(j)*scale
     >     ,visdif(j)*scale
     >     ,(convec(j)-acc_u2(1,i,j)+v2_osciltm_C(1,i,j))*scale
     >     ,transp(j)*scale
     >     ,vpgrad(j)*scale
     >     ,v2_osciltm_P(1,i,j)*scale
     >     ,v2_osciltm_C(1,i,j)*scale
     >     ,sum*scale

100   format(11(e12.5,x))    

c---


      enddo
      enddo
      close(10)

      return
      end

c*************** two_d_uu2 *********************************************
      subroutine two_d_uu2
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/fileio/ith,k_phase
      common/indx1/ipa(m1),imu(m1),imv(m1),kpa(m3),kma(m3)
      common/indx2/jpa(m2),jmu(m2),jmv(m2)
      common/mesh1/dx1,dx1q,dx3,dx3q
      common/mesh2/y(0:m2)
      common/mesh3/dy(0:m2),h(m2),hm(m2),hc(m2),hp(m2)
      common/mesh4/dym(m2),dyc(m2),dyp(m2)
      common/para/re
      common/size/alx,aly,alz,vol
      common/utau/u_tau_no(m1)

      common/budget1/u(3,0:m1,0:m2),p(0:m1,0:m2),uu(6,m1,m2)
     >              ,uuu(6,2,m1,m2)
      common/budget2/upg(6,m1,m2),diss(6,m1,m2),presstrain(6,m1,m2)
      common/budget3/pp(0:m1,0:m2)
      common/budget5/acc_u2(6,m1,m2)
      common/budget6/v1_osciltm(2,m1,m2)      ! overbar{u~_i*u~_j}
     >              ,v2_osciltm_P(6,m1,m2)    ! Production term
     >              ,v2_osciltm_C(6,m1,m2)    ! Convection term

      real produc(m2)
      real dissip(m2)
      real visdif(m2)
      real convec(m2)
      real transp(m2)
      real vpgrad(m2)


      integer moni_x

      character*30 name

      name='2d_budet_v2'
      nn=index(name,'v')
      write(unit=name(nn+2:),fmt='(bn,a2,i1.1,a4)')
     >           '_p',k_phase,'.dat'



      open(10,file=name,status='unknown')
      write(10,11) 
11    format('variables=x,y,P,`e,D,C,T,`P,P_o_s,C_o_s,sum')
      write(10,12) n1m-20, n2m-2
12    format('zone i=',i4,', j=',i4,', f=point')
     
      do j=2,n2m-1
      jp=j+1
      jm=j-1

      do i=11,n1m-10
      ip=i+1
      im=i-1

c     production term
      produc(j)=-2.*(
     >   uu(4,i,j)*(u(2,ip,j)-u(2,im,j))*dx1/2.  
     >  +uu(2,i,j)*  
     >    ((h(j)**2)*u(2,i,jp)
     >   -(h(j)**2-h(jp)**2)*u(2,i,j)
     >   -(h(jp)**2)*u(2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     >  )                                 ! -2<u_2u_k>u2,k

c     dissipation term
      dissip(j)=-2./re*diss(2,i,j)

c     viscous diffusion term
      visdif(j)=1./re*(
     >   (uu(2,ip,j)-2.*uu(2,i,j)+uu(2,im,j))*dx1q
     >  +(hp(j)*uu(2,i,jp)
     >   -hc(j)*uu(2,i,j )
     >   +hm(j)*uu(2,i,jm))
     > ) 
     > 

c     convection term      
      convec(j)=-(
     >    u(1,i,j)*(uu(2,ip,j)-uu(2,im,j))*dx1/2.
     > +  u(2,i,j)*
     >    ((h(j)**2)*uu(2,i,jp)
     >   -(h(j)**2-h(jp)**2)*uu(2,i,j)
     >   -(h(jp)**2)*uu(2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     >  )                                 ! -u_k*<u_2u_2>,k

c     turbulent transport term
      transp(j)=-(
     >    (uuu(2,1,ip,j)-uuu(2,1,im,j))*dx1/2.
     > +  
     >    ((h(j)**2)*uuu(2,2,i,jp)
     >   -(h(j)**2-h(jp)**2)*uuu(2,2,i,j)
     >   -(h(jp)**2)*uuu(2,2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     > )

c     velocity gradient term
      vpgrad(j)=-upg(2,i,j)

c---
      sum=produc(j)+dissip(j)+visdif(j)+convec(j)
     >   +transp(j)+vpgrad(j)
     >   -acc_u2(2,i,j)
     >   +v2_osciltm_P(2,i,j)
     >   +v2_osciltm_C(2,i,j)

      scale=1./(u_tau_no(i)**4*RE)

      xc=(i-0.5)/dx1
      yc=0.5*(y(j)+y(j+1))

                   x_slot=75.0
                   u_tau_in=0.0492
c                   x_slot=78.90625
c                   u_tau_in=0.0525

         write(10,100)  (xc-x_slot)*u_tau_in*re
     >                , yc*u_tau_no(i)*re
     >     ,(produc(j)+v2_osciltm_P(2,i,j))*scale
     >     ,dissip(j)*scale
     >     ,visdif(j)*scale
     >     ,(convec(j)-acc_u2(2,i,j)+v2_osciltm_C(2,i,j))*scale
     >     ,transp(j)*scale
     >     ,vpgrad(j)*scale
     >     ,v2_osciltm_P(2,i,j)*scale
     >     ,v2_osciltm_C(2,i,j)*scale
     >     ,sum*scale

100   format(11(e12.5,x))    
c---

      enddo
      enddo
      close(10)

      return
      end

c*************** two_d_uu4 *********************************************
      subroutine two_d_uu4
      include '../PARAM.H'
      common/dim/n1,n2,n3,n1m,n2m,n3m,n2mm
      common/fileio/ith,k_phase
      common/indx1/ipa(m1),imu(m1),imv(m1),kpa(m3),kma(m3)
      common/indx2/jpa(m2),jmu(m2),jmv(m2)
      common/mesh1/dx1,dx1q,dx3,dx3q
      common/mesh2/y(0:m2)
      common/mesh3/dy(0:m2),h(m2),hm(m2),hc(m2),hp(m2)
      common/mesh4/dym(m2),dyc(m2),dyp(m2)
      common/para/re
      common/size/alx,aly,alz,vol
      common/utau/u_tau_no(m1)

      common/budget1/u(3,0:m1,0:m2),p(0:m1,0:m2),uu(6,m1,m2)
     >              ,uuu(6,2,m1,m2)
      common/budget2/upg(6,m1,m2),diss(6,m1,m2),presstrain(6,m1,m2)
      common/budget3/pp(0:m1,0:m2)
      common/budget5/acc_u2(6,m1,m2)
      common/budget6/v1_osciltm(2,m1,m2)      ! overbar{u~_i*u~_j}
     >              ,v2_osciltm_P(6,m1,m2)    ! Production term
     >              ,v2_osciltm_C(6,m1,m2)    ! Convection term
      common/budget7/rss_os(2,m1,m2)

      real produc(m2)
      real dissip(m2)
      real visdif(m2)
      real convec(m2)
      real transp(m2)
      real vpgrad(m2)

      integer moni_x

      character*30 name

      name='2d_budet_uv'
      nn=index(name,'v')
      write(unit=name(nn+1:),fmt='(bn,a2,i1.1,a4)')
     >           '_p',k_phase,'.dat'



      open(10,file=name,status='unknown')
      write(10,11) 
11    format('variables=x,y,P,`e,D,C,T,`P,P_o_s,C_o_s,sum')
      write(10,12) n1m-20, n2m-2
12    format('zone i=',i4,', j=',i4,', f=point')
     
      do j=2,n2m-1
      jp=j+1
      jm=j-1

      do i=11,n1m-10
      ip=i+1
      im=i-1

c     production term
      produc(j)=-(
     >   uu(1,i,j)*(u(2,ip,j)-u(2,im,j))*dx1/2.  
     >  +uu(4,i,j)*  
     >    ((h(j)**2)*u(2,i,jp)
     >   -(h(j)**2-h(jp)**2)*u(2,i,j)
     >   -(h(jp)**2)*u(2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))              !<u_1u_k>u2,k
     >  +uu(4,i,j)*(u(1,ip,j)-u(1,im,j))*dx1/2.  
     >  +uu(2,i,j)*  
     >    ((h(j)**2)*u(1,i,jp)
     >   -(h(j)**2-h(jp)**2)*u(1,i,j)
     >   -(h(jp)**2)*u(1,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))              !<u_2u_k>u1,k
     >  )                                 ! -(<u_1u_k>u2,k+<u_2u_k>u1,k)

c     dissipation term
      dissip(j)=-2./re*diss(4,i,j)

c     viscous diffusion term
      visdif(j)=1./re*(
     >   (uu(4,ip,j)-2.*uu(4,i,j)+uu(4,im,j))*dx1q
     >  +(hp(j)*uu(4,i,jp)
     >   -hc(j)*uu(4,i,j )
     >   +hm(j)*uu(4,i,jm))
     > ) 
     > 

c     convection term      
      convec(j)=-(
     >    u(1,i,j)*(uu(4,ip,j)-uu(4,im,j))*dx1/2.
     > +  u(2,i,j)*
     >    ((h(j)**2)*uu(4,i,jp)
     >   -(h(j)**2-h(jp)**2)*uu(4,i,j)
     >   -(h(jp)**2)*uu(4,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     >  )                                 ! -u_k*<u_1u_2>,k

c     turbulent transport term
      transp(j)=-(
     >    (uuu(4,1,ip,j)-uuu(4,1,im,j))*dx1/2.
     > +  
     >    ((h(j)**2)*uuu(4,2,i,jp)
     >   -(h(j)**2-h(jp)**2)*uuu(4,2,i,j)
     >   -(h(jp)**2)*uuu(4,2,i,jm))
     >   /h(j)/h(jp)/(h(j)+h(jp))
     > )

c     velocity gradient term
      vpgrad(j)=-upg(4,i,j)

c---
      sum=produc(j)+dissip(j)+visdif(j)+convec(j)
     >   +transp(j)+vpgrad(j)
     >   -acc_u2(4,i,j)
     >   +v2_osciltm_P(4,i,j)
     >   +v2_osciltm_C(4,i,j)

      scale=1./(u_tau_no(i)**4*RE)

      xc=(i-0.5)/dx1
      yc=0.5*(y(j)+y(j+1))

                   x_slot=75.0
                   u_tau_in=0.0492
c                   x_slot=78.90625
c                   u_tau_in=0.0525

         write(10,100)  (xc-x_slot)*u_tau_in*re
     >                , yc*u_tau_no(i)*re
     >     ,-(produc(j)+v2_osciltm_P(4,i,j))*scale
     >     ,-dissip(j)*scale
     >     ,-visdif(j)*scale
     >     ,-(convec(j)-acc_u2(4,i,j)+v2_osciltm_C(4,i,j))*scale
     >     ,-transp(j)*scale
     >     ,-vpgrad(j)*scale
     >     ,-v2_osciltm_P(4,i,j)*scale
     >     ,-v2_osciltm_C(4,i,j)*scale
     >     ,-sum*scale

100   format(11(e12.5,x))    
c---




      enddo
      enddo
      close(10)
      close(11)

      return
      end
