c  **********************  set boundary conditions  *********
c
c      subroutine boucdq(q1,q2,q3,time)          !PARALL
      subroutine boucdq(time,comm)
      include 'param.f'
      include "mpif.h"       !PARALL         

      integer comm,ierr
c      dimension qb2so(m1,m2)          !PARALL
c      dimension cou2s(m1,m2),cou2n(m1,m2)
c      dimension cou3s(m1,m2),cou3n(m1,m2)
c      dimension cou1n(m1,m2)
c      dimension cou3l(m1,m3),cou2l(m1,m3),cou1l(m1,m3)          !PARALL
      pi=2.*asin(1.)
c
c
      x3m(n3)=x3c(n3)
c     cou2m = 0.
      
      do i=1,n1m
c         do k=1,n3m      !LAURA
         do k=1,n3

            dqb3dn(i,k)=uzero*ft*cos(alphax)
c    $           -Ros*wx*(x2m(1)-x20) +Ros*wy*(x1m(i)-x10)-u0z
     $           -qb3dn(i,k)          !Modifica nuova versione NS
            dqb3up(i,k)=uzero*ft*cos(alphax)
c     $        - Ros*wx*(x2m(n2m)-x20) + Ros*wy*(x1m(i)-x10)-u0z
     $           - qb3up(i,k)         !Modifica nuova versione NS
            dqb2dn(i,k)=uzero*ft*sin(alphax)
c     $        +Ros*wx*(x3m(k)-x30)-u0y
     $           -qb2dn(i,k)          !Modifica nuova versione NS
            dqb2up(i,k)=uzero*ft*sin(alphax)
c     $        +Ros*wx*(x3m(k)-x30)-u0y 
     $           -qb2up(i,k)          !Modifica nuova versione NS

         enddo
      enddo

      do i=1,n1m         !condiz outflow coincide condiz fondo Q_2
        dqb2n(i,1)=dqb2dn(i,n3)
      enddo

      j = 1
      do i=1,n1m

         d2d1p1 = 1. + udx3m(n3m)/udx3m(n3m-1)
         dq3dx3 = ( q3(i,j,n3m-1) - d2d1p1*d2d1p1*q3(i,j,n3m)
     $        + qb3n(i,j)*(d2d1p1*d2d1p1 -1.)
     $        )*udx3m(n3m-1)/d2d1p1
         dq3x2 = 
c     $        wx*( qb2n(i,j) + qb2n(i,j+1) )*0.5
c     $        + ( uinfth(i,j)*sin(alphax) 
c     $        +wx*( x3c(n3)-x30 ) - uhy*cos(alphax)  
c     $      )*( (qb3n(i,j+1)+qb3n(i,j))*0.5-qb3dn(i,n3) )*udx2m(j)
     $        ( -dwx*( x2m(j)-x20 ) + wx*uhy*cos(alphax)
     $        + duhy*sin(alphax) )
     $        + ( uinfth(i,j)*cos(alphax) 
     $        -wx*( x2m(j)-x20 ) + uhy*sin(alphax)   
     $        )*dq3dx3
C         dqb3n(i,j)=-dt*(ga*dq3x2+ro*dq3x2o(i,j))
C         dq3x2o(i,j)=dq3x2
         dqb3n(i,j)=-dt*(ga*dq3dx3+ro*dq3x2o(i,j))*uzero    !RV
         dq3x2o(i,j)=dq3dx3    !RV

      enddo

c      do j=1,n2m   !!LAURA
      do j=2,n2m-1
        do i=1,n1m

c    radiative b.c. for radial velocity at the outflow    
C          dq2x2=(qb2n(i,j)-q2(i,j,n3m))*2.*udx3c(n3)   !nuove outflow         
           d2d1p1 = 1. + 2.*udx3c(n3)/udx3c(n3m)
           dq2dx3 = ( q2(i,j,n3m-1) - d2d1p1*d2d1p1*q2(i,j,n3m)
     $          + qb2n(i,j)*(d2d1p1*d2d1p1 -1.)
     $          )*udx3c(n3m)/d2d1p1
           dq2x2 = 
c     $          -wx*( qb3n(i,j) + qb3n(i,j-1) )*0.5 
c     $          + ( uinfth(i,j)*sin(alphax) 
c     $          +wx*( x3c(n3)-x30 ) - uhy*cos(alphax)     
c     $           )*( qb2n(i,j+1) - qb2n(i,j-1) )*0.5*udx2c(j)
     $        ( dwx*( x3c(n3)-x30 ) + wx*uhy*sin(alphax)
     $        - duhy*cos(alphax) )
     $          + ( uinfth(i,j)*cos(alphax) 
     $          -wx*( x2c(j)-x20 ) + uhy*sin(alphax)        
     $           )*( dq2dx3 + wx )
C     $          ( qb2n(i,j) - q2(i,j,n3m) )*2.*udx3c(n3)  !!I ordine

c          dqb2n(i,j)=-dt*(ga*dq2x2+ro*dq2x2o(i,j))*cou
C          dqb2n(i,j)=-dt*(ga*dq2x2+ro*dq2x2o(i,j))
           dqb2n(i,j)=-dt*(ga*dq2dx3+ro*dq2x2o(i,j))*uzero  !RV
         
c          if (i.eq.13.and.j.eq.2)write(94,*) 'boucdq -cost',
c     $         ( uinfth(i,j)*sin(alphax) +
c     $          wx*( x3c(n3)-x30 ) - uhy*cos(alphax) ),
c     $         ( uinfth(i,j)*cos(alphax) -
c     $          wx*( x2c(j)-x20 ) + uhy*sin(alphax) )
c          if (i.eq.13.and.j.eq.2)write(94,*)'boucdq',i,dq2dx3,
c     $         ( qb2n(i,j+1) - qb2n(i,j-1) )*0.5*udx2c(j),
c     $         dq2x2,dq2x2o(i,j),dqb2n(i,j),qb2n(i,j)
c          if (i.eq.14.and.j.eq.2)write(94,*)'boucdq',i,dq2dx3,
c     $         ( qb2n(i,j+1) - qb2n(i,j-1) )*0.5*udx2c(j),
c     $         dq2x2,dq2x2o(i,j),dqb2n(i,j),qb2n(i,j)
c          if (i.eq.15.and.j.eq.2)write(94,*)'boucdq',i,dq2dx3,
c     $         ( qb2n(i,j+1) - qb2n(i,j-1) )*0.5*udx2c(j),
c     $         dq2x2,dq2x2o(i,j),dqb2n(i,j),qb2n(i,j)
c          if (i.eq.16.and.j.eq.2)write(94,*)'boucdq',i,dq2dx3,
c     $         ( qb2n(i,j+1) - qb2n(i,j-1) )*0.5*udx2c(j),
c     $         dq2x2,dq2x2o(i,j),dqb2n(i,j),qb2n(i,j)

c          dq2x2o(i,j)=dq2x2
          dq2x2o(i,j)=dq2dx3  !RV

c    radiative b.c. for axial velocity at the outflow
c          dq3x2=(qb3n(i,j)-q3(i,j,n3m))*udx3c(n3)   !nuove outflow
          d2d1p1 = 1. + udx3m(n3m)/udx3m(n3m-1)
          dq3dx3 = ( q3(i,j,n3m-1) - d2d1p1*d2d1p1*q3(i,j,n3m)
     $        + qb3n(i,j)*(d2d1p1*d2d1p1 -1.)
     $        )*udx3m(n3m-1)/d2d1p1
          dq3x2 = 
c     $         wx*( qb2n(i,j) + qb2n(i,j+1) )*0.5
c     $          + ( uinfth(i,j)*sin(alphax) 
c     $          +wx*( x3c(n3)-x30 ) - uhy*cos(alphax)        
c     $          )*( qb3n(i,j+1) - qb3n(i,j-1) )*0.5*udx2m(j)
     $        ( -dwx*( x2m(j)-x20 ) + wx*uhy*cos(alphax)
     $        + duhy*sin(alphax) )
     $          + ( uinfth(i,j)*cos(alphax) 
     $          -wx*( x2m(j)-x20 ) + uhy*sin(alphax)        
     $          )*dq3dx3
c     $          ( qb3n(i,j) - q3(i,j,n3m) )*udx3m(n3m)  !!I ordine

c          dqb3n(i,j)=-dt*(ga*dq3x2+ro*dq3x2o(i,j))*cou
c          dqb3n(i,j)=-dt*(ga*dq3x2+ro*dq3x2o(i,j))
          dqb3n(i,j)=-dt*(ga*dq3dx3+ro*dq3x2o(i,j))*uzero    !RV

c          if (j.eq.2.and.i.eq.13)write(*,*)'BOUCDQ',I,dq3x2,dqb3n(i,j),
c     $         dq3x2o(i,j)
c          if (j.eq.2.and.i.eq.14)write(*,*)'BOUCDQ',I,dq3x2,dqb3n(i,j),
c     $         dq3x2o(i,j)
c          if (j.eq.2.and.i.eq.15)write(*,*)'BOUCDQ',I,dq3x2,dqb3n(i,j),
c     $         dq3x2o(i,j)
c          if (j.eq.2.and.i.eq.16)write(*,*)'BOUCDQ',I,dq3x2,dqb3n(i,j),
c     $         dq3x2o(i,j)

c          if (i.eq.13.and.j.eq.2)write(98,*)'boucdq',i,dq3dx3,
c     $         ( qb3n(i,j+1) - qb3n(i,j-1) )*0.5*udx2m(j),
c     $         dq3x2,dq3x2o(i,j),dqb3n(i,j),qb3n(i,j)
c          if (i.eq.13.and.j.eq.2)write(98,*)'dq3dx3',i,dq3dx3,
c     $         q3(i,j,n3m-1),q3(i,j,n3m),qb3n(i,j)
c          if (i.eq.14.and.j.eq.2)write(98,*)'boucdq',i,dq3dx3,
c     $         ( qb3n(i,j+1) - qb3n(i,j-1) )*0.5*udx2m(j),
c     $         dq3x2,dq3x2o(i,j),dqb3n(i,j),qb3n(i,j)
c          if (i.eq.14.and.j.eq.2)write(98,*)'dq3dx3',i,dq3dx3,
c     $         q3(i,j,n3m-1),q3(i,j,n3m),qb3n(i,j)
c          if (i.eq.15.and.j.eq.2)write(98,*)'boucdq',i,dq3dx3,
c     $         ( qb3n(i,j+1) - qb3n(i,j-1) )*0.5*udx2m(j),
c     $         dq3x2,dq3x2o(i,j),dqb3n(i,j),qb3n(i,j)
c          if (i.eq.15.and.j.eq.2)write(98,*)'dq3dx3',i,dq3dx3,
c     $         q3(i,j,n3m-1),q3(i,j,n3m),qb3n(i,j)
c          if (i.eq.16.and.j.eq.2)write(98,*)'boucdq',i,dq3dx3,
c     $         ( qb3n(i,j+1) - qb3n(i,j-1) )*0.5*udx2m(j),
c     $         dq3x2,dq3x2o(i,j),dqb3n(i,j),qb3n(i,j)
c          if (i.eq.16.and.j.eq.2)write(98,*)'dq3dx3',i,dq3dx3,
c     $         q3(i,j,n3m-1),q3(i,j,n3m),qb3n(i,j)

c          dq3x2o(i,j)=dq3x2
          dq3x2o(i,j)=dq3dx3  !RV
          
c    imposed condition for radial velocity at the outflow
c         dqb2n(i,j)=uinfth(i,j)*ft*sin(alphax)+
c    %      Ros*wx*(x3m(n3)-x30)-u0y-qb2n(i,j)
c    imposed condition for axial velocity at the outflow
c         dqb3n(i,j)=uinfth(i,j)*ft*cos(alphax)
c    %     -Ros*wx*(x2m(j)-x20)-u0z-qb3n(i,j)
        enddo   
      enddo 

      j = n2m
      do i=1,n1m

           d2d1p1 = 1. + 2.*udx3c(n3)/udx3c(n3m)
           dq2dx3 = ( q2(i,j,n3m-1) - d2d1p1*d2d1p1*q2(i,j,n3m)
     $          + qb2n(i,j)*(d2d1p1*d2d1p1 -1.)
     $          )*udx3c(n3m)/d2d1p1
           dq2x2 = 
c     $          -wx*( qb3n(i,j) + qb3n(i,j-1) )*0.5       
c     $          + ( uinfth(i,j)*sin(alphax) 
c     $          +wx*( x3c(n3)-x30 ) - uhy*cos(alphax)       
c     $          )*( qb2up(i,n3) - qb2n(i,j-1) )*0.5*udx2c(j)
     $        ( dwx*( x3c(n3)-x30 ) + wx*uhy*sin(alphax)
     $        - duhy*cos(alphax) )
     $          + ( uinfth(i,j)*cos(alphax) 
     $          -wx*( x2c(j)-x20 ) + uhy*sin(alphax)     
     $          )*( dq2dx3 + wx )
c     $          ( qb2n(i,j) - q2(i,j,n3m) )*2.*udx3c(n3) !!I ordine
c          dqb2n(i,j)=-dt*(ga*dq2x2+ro*dq2x2o(i,j))
c          dq2x2o(i,j)=dq2x2
          dqb2n(i,j)=-dt*(ga*dq2dx3+ro*dq2x2o(i,j))*uzero !rv
          dq2x2o(i,j)=dq2dx3    !RV

          d2d1p1 = 1. + udx3m(n3m)/udx3m(n3m-1)
          dq3dx3 = ( q3(i,j,n3m-1) - d2d1p1*d2d1p1*q3(i,j,n3m)
     $         + qb3n(i,j)*(d2d1p1*d2d1p1 -1.)
     $         )*udx3m(n3m-1)/d2d1p1
          dq3x2 = 
c     $         wx*( qb2n(i,j) + qb2up(i,n3) )*0.5  
c     $          + ( uinfth(i,j)*sin(alphax)     
c     $          +wx*( x3c(n3)-x30 ) - uhy*cos(alphax)    
c     $           )*( qb3up(i,n3) - 
c     $            ( qb3n(i,j)+qb3n(i,j-1) )*0.5)*udx2m(j)
     $        ( -dwx*( x2m(j)-x20 ) + wx*uhy*cos(alphax)
     $        + duhy*sin(alphax) )
     $          + ( uinfth(i,j)*cos(alphax) 
     $          -wx*( x2m(j)-x20 ) + uhy*sin(alphax)    
     $          )*dq3dx3
c     $          ( qb3n(i,j) - q3(i,j,n3m) )*udx3m(n3m) !!I ordine
c          dqb3n(i,j)=-dt*(ga*dq3x2+ro*dq3x2o(i,j))
c          dq3x2o(i,j)=dq3x2
          dqb3n(i,j)=-dt*(ga*dq3dx3+ro*dq3x2o(i,j))*uzero    !RV
          dq3x2o(i,j)=dq3dx3    !RV

      enddo

      do j=1,n2m
        do i=1,n1m
         uinfth(i,j)=uzero
c    b.c. for radial velocity at the inflow    
          dqb2s(i,j)=+uinfth(i,j)*ft*sin(alphax)
C     %     +Ros*wx*(x3m(1)-x30)-u0y    !Modifica nuova versione NS
     $         -qb2s(i,j)
c    b.c. for axial velocity at the inflow    
            dqb3s(i,j) =uinfth(i,j)*ft*cos(alphax)
C     %      -Ros*wx*(x2m(j)-x20)-u0z   !Modifica nuova versione NS
     $         -qb3s(i,j)

        enddo   
      enddo   

c
c      if(jv3.gt.0) then
c        do j=1,jv3
c         do i=1,n1m
c           dqb3s(i,j) =(uinfth(i,j)*ft*cos(alphax)
c    %      +Ros*wy*(x1m(i)-x10)-Ros*wx*(x2m(j)-x20))-u0z-qb3s(i,j)
c         enddo   
c        enddo   
c       endif
c       do j=jv3+1,n2m
c        do i=1,n1m
c         dqb3s(i,j)=0.
c        enddo   
c       enddo
                
      dq3ct=0.
      area = 0.
      qout=0.
      qinf=0.
      do j=1,n2m
        darea=g2m(j)/dx2/dx1
        do i=1,n1m
          area = area+darea
          qout=qout+dqb3n(i,j)*darea
          qinf=qinf+dqb3s(i,j)*darea
        end do
      end do
c     write(6,*) ' Q_in = ',qinf,' Q_ou = ',qout
      call MPI_ALLREDUCE(qout,qoutall,1,MPI_DOUBLE_PRECISION,
     $     MPI_SUM,comm,ierr)   ! PARALL
      call MPI_ALLREDUCE(qinf,qinfall,1,MPI_DOUBLE_PRECISION,
     $     MPI_SUM,comm,ierr)   ! PARALL
      call MPI_ALLREDUCE(area,areaall,1,MPI_DOUBLE_PRECISION,
     $     MPI_SUM,comm,ierr)   ! PARALL

      cor=(qinfall-qoutall)/areaall    ! PARALL
c      if (myid.eq.0)write(*,*)' boucdq ',QOUT,QINF,AREA
c      if (myid.eq.0)write(*,*)' boucdq ',QOUTall,QINFall,AREAall,cor
c      cor=(qinf-qout)/area
c      write(93,*)time,cor
      qout2=0.
      qoutf=0.
      do j=1,n2m
        darea=g2m(j)/dx2/dx1
        do i=1,n1m
          qcc=cor
          dqb3n(i,j)=dqb3n(i,j)+qcc
          dq3x2o(i,j) = dq3x2o(i,j) - qcc/(dt*ga)   !!!laura
          qout2=qout2+qcc*darea
          qoutf=qoutf+dqb3n(i,j)*darea
        end do
      end do
      if (time.eq.0.)then
         do j=1,n2m
            do i=1,n1m
               dqb2n(i,j)=dqb2s(i,j)
               dq3x2o(i,j) = dq3x2o(i,j) - dqb2s(i,j)/(dt*ga) !!!laura
            enddo
         enddo
      endif
c
c   b.c. for azimuthal velocity from continuity eq.
c
      if(n1m.ne.1) then
         d2d1p1_3 = 1. + 2.*udx3c(n3)/udx3c(n3m)
         do i=1,n1m
            j = 1
            d2d1p1_2 = 1. + udx2c(2)/udx2c(3)
            dq1dx2 = ( - qb1n(i,3) + d2d1p1_2*d2d1p1_2*qb1n(i,2)
     $           - qb1n(i,1)*(d2d1p1_2*d2d1p1_2 -1.)
     $           )*udx2c(3)/d2d1p1_2
            dq1dx3 = ( q1(i,j,n3m-1) - d2d1p1_3*d2d1p1_3*q1(i,j,n3m)
     $           + qb1n(i,j)*(d2d1p1_3*d2d1p1_3 - 1. )
     $           )*udx3c(n3m)/d2d1p1_3
            dq1x2 = 
c     $           ( uinfth(i,j)*sin(alphax) 
c     $          +wx*( x3c(n3)-x30 ) - uhy*cos(alphax)   
c     $           )*dq1dx2 +
C     $          ( qb1n(i,j+1) - qb1n(i,j) )*udx2c(j+1) +   !!I ordine
     $          ( uinfth(i,j)*cos(alphax)
     $           - wx*( x2m(j)-x20 ) + uhy*sin(alphax)    
     $              )*dq1dx3
c     $          ( qb1n(i,j) - q1(i,j,n3m) )*2.*udx3c(n3)  !!I ordine
c            dqb1n(i,j)= -dt*(ga*dq1x2+ro*dq1x2o(i,j))
c            dq1x2o(i,j)=dq1x2
            dqb1n(i,j)= -dt*(ga*dq1dx3+ro*dq1x2o(i,j))*uzero    !RV
            dq1x2o(i,j)=dq1dx3   !rv
            dqb1s(i,j)=0.

c            do j=1,n2m  !!!!laura
            do j=2,n2m-1

c     dq1x2=(qb1n(i,j)-q1(i,j,n3m))*2.*dx3/g3c(n3)   !nuove outflow
               dq1dx3 = ( q1(i,j,n3m-1) - d2d1p1_3*d2d1p1_3*q1(i,j,n3m)
     $              + qb1n(i,j)*(d2d1p1_3*d2d1p1_3 -1.)
     $              )*udx3c(n3m)/d2d1p1_3
               dq1x2 = 
c     $              ( uinfth(i,j)*sin(alphax) 
c     $              +wx*( x3c(n3)-x30 ) - uhy*cos(alphax)    
c     $              )*( qb1n(i,j+1) - qb1n(i,j-1) )*0.5*udx2m(j) +
     $              ( uinfth(i,j)*cos(alphax)
     $              - wx*( x2m(j)-x20 ) + uhy*sin(alphax)     
     $              )*dq1dx3
c     $              ( qb1n(i,j) - q1(i,j,n3m) )*2.*udx3c(n3)  !!I ordine
               
c     dqb1n(i,j)=-dt*(ga*dq1x2+ro*dq1x2o(i,j))*cou
c               dqb1n(i,j)= -dt*(ga*dq1x2+ro*dq1x2o(i,j))
c               dq1x2o(i,j)=dq1x2
               dqb1n(i,j)= -dt*(ga*dq1dx3+ro*dq1x2o(i,j))*uzero    !RV
               dq1x2o(i,j)=dq1dx3    !RV
               
c     Modifica nuova versione NS
c          dqb1s(i,j)=0.-Ros*wy*(x3m(1)-x30)+Ros*wz*(x2m(j)-x20)
c     %    -uinfth(i,j)*ft*sin(alphay)-u0x-qb1s(i,j)
               dqb1s(i,j)=0.    !??????
c     Modifica nuova versione NS
            end do

            j = n2m
            d2d1p1_2 = 1.+udx2c(n2m)/udx2c(n2m-1)
            dq1dx2 = ( qb1n(i,j-2) - d2d1p1_2*d2d1p1_2*qb1n(i,j-1)
     $           + qb1n(i,j)*(d2d1p1_2*d2d1p1_2 -1.)
     $           )*udx2c(n2m-1)/d2d1p1_2
            dq1dx3 = ( q1(i,j,n3m-1) - d2d1p1_3*d2d1p1_3*q1(i,j,n3m)
     $           + qb1n(i,j)*(d2d1p1_3*d2d1p1_3 - 1. )
     $           )*udx3c(n3m)/d2d1p1_3
            dq1x2 = 
c     $           ( uinfth(i,j)*sin(alphax) 
c     $          +wx*( x3c(n3)-x30 ) - uhy*cos(alphax) 
c     $            )*dq1dx2 +
c     $          ( qb1n(i,j) - qb1n(i,j-1) )*udx2c(j) +  !!I ordine
     $          ( uinfth(i,j)*cos(alphax)
     $           - wx*( x2m(j)-x20 ) + uhy*sin(alphax)    
     $               )*dq1dx3
c     $          ( qb1n(i,j) - q1(i,j,n3m) )*2.*udx3c(n3)  !!I ordine
c            dqb1n(i,j)= -dt*(ga*dq1x2+ro*dq1x2o(i,j))
c            dq1x2o(i,j)=dq1x2
            dqb1n(i,j)= -dt*(ga*dq1dx3+ro*dq1x2o(i,j))*uzero   !RV
            dq1x2o(i,j)=dq1dx3   !RV
            dqb1s(i,j)=0.

         end do
      end if
      
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
c      do i=1,n1m
c         do j=1,n2m
c            dqb3n(i,j) = dqb3s(i,j)
c            dqb2n(i,j) = dqb2s(i,j)
c            dqb1n(i,j) = dqb1s(i,j)
c         enddo
c      enddo
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

c     ddq2nm = -10.
c     ddq2sm = -10.
c     ddq3nm = -10.
c     ddq3sm = -10.
c     do j=1,n2m
c       do i=1,n1m
c         ddq2nm = max(ddq2nm,abs(dqb2n(i,j)))
c         ddq2sm = max(ddq2sm,abs(dqb2s(i,j)))
c         ddq3nm = max(ddq3nm,abs(dqb3n(i,j)))
c         ddq3sm = max(ddq3sm,abs(dqb3s(i,j)))
c       end do
c     end do
c     write(6,*) ' In boucon '
c     write(6,*) ddq2nm,ddq2sm,ddq3nm,ddq3sm
c
c     write(6,*)'BOUCDQ'
c     write(6,*) ' ===>>> T = ',time,' Q_3 = ',qb3s(1,n2m/2)
c     write(6,*) ' ===>>> T = ',time,' DQ_3 = ',dqb3s(1,n2m/2)
c     write(6,*) ' ===>>> T = ',time,' Q_1 = ',qb1s(1,n2m/2)
c     write(6,*) ' ===>>> T = ',time,' DQ_1 = ',dqb1s(1,n2m/2)

      alln3 = x3c(n3)
      alll=5.*alln3/6.
      do k=1,n3
c      if(x3c(k).le.alll) then
        visc(k) = 1.
c      else
c       visc(k) = (10.-1.)/(alln3-alll)**2*(x3c(k)-alll)**2+1.
c      endif
      end do
      return
      end
c
c  **********************  set boundary conditions  *********
c
      subroutine boucqt(time)
c      subroutine boucqt(q1,q2,q3,time)          ! PARALL
      include 'param.f'

      common/qiou/qinf,qout     !sembra inutile PARALL
c
      qout=0.
      qinf=0.
      do j=1,n2m
        darea=g2m(j)/dx2/dx1
        do i=1,n1m
          qb1s(i,j) = qb1s(i,j)+dqb1s(i,j)
          qb2s(i,j) = qb2s(i,j)+dqb2s(i,j)
c         qb3s(i,j) = ft*uinfth(i,j)
          qb3s(i,j) = qb3s(i,j)+dqb3s(i,j)
          qb1n(i,j) = qb1n(i,j)+dqb1n(i,j)
          qb2n(i,j) = qb2n(i,j)+dqb2n(i,j)
          qb3n(i,j) = qb3n(i,j)+dqb3n(i,j)
          qout=qout+qb3n(i,j)*darea
          qinf=qinf+qb3s(i,j)*darea
C      if (i.eq.13.and.j.eq.2)write(94,*)'boucqt',i,qb2n(i,j),dqb2n(i,j)          ! PARALL
C      if (i.eq.13.and.j.eq.2)write(98,*)'boucqt',i,qb3n(i,j),dqb3n(i,j)
C      if (i.eq.14.and.j.eq.2)write(94,*)'boucqt',i,qb2n(i,j),dqb2n(i,j)
C      if (i.eq.14.and.j.eq.2)write(98,*)'boucqt',i,qb3n(i,j),dqb3n(i,j)
C      if (i.eq.15.and.j.eq.2)write(94,*)'boucqt',i,qb2n(i,j),dqb2n(i,j)
C      if (i.eq.15.and.j.eq.2)write(98,*)'boucqt',i,qb3n(i,j),dqb3n(i,j)
C      if (i.eq.16.and.j.eq.2)write(94,*)'boucqt',i,qb2n(i,j),dqb2n(i,j)
C      if (i.eq.16.and.j.eq.2)write(98,*)'boucqt',i,qb3n(i,j),dqb3n(i,j)
       enddo
      enddo

c      do k=1,n3m !LAURA
      do k=1,n3
        do i=1,n1m
          qb2dn(i,k) = qb2dn(i,k)+dqb2dn(i,k)
          qb3dn(i,k) = qb3dn(i,k)+dqb3dn(i,k)
          qb2up(i,k) = qb2up(i,k)+dqb2up(i,k)
          qb3up(i,k) = qb3up(i,k)+dqb3up(i,k)
        enddo
      enddo  

      do i = 1, n1m    !PROVO IMPORRE QB2N SUL FONDO
         qb2n(i,1) = qb2dn(i,n3)
      enddo
      
c      write(6,*)'BOUCQT'
c      write(6,*) ' ===>>> T = ',time,' Q_3! = ',qb3s(5,50)
c      write(6,*) ' ===>>> T = ',time,' DQ_3! = ',dqb3s(5,50)
c      write(6,*) ' ===>>> T = ',time,' Q_3_dn = ',qb3dn(1,1)
c      write(6,*) ' ===>>> T = ',time,' DQ_3_dn = ',dqb3dn(1,1)
c      write(6,*) ' ===>>> T = ',time,' Q_3_up = ',qb3up(1,1)
c      write(6,*) ' ===>>> T = ',time,' DQ_3_up = ',dqb3up(1,1)
c      write(6,*) ' ===>>> T = ',time,' Q_1 = ',qb1s(1,1)
c      write(6,*) ' ===>>> T = ',time,' DQ_1 = ',dqb1s(1,1)
c      write(6,*) ' ===>>> T = ',time,' Q_2!! = ',qb2s(5,50),q2(5,50,1)
c      write(6,*) ' ===>>> T = ',time,' DQ_2! = ',dqb2s(5,50)
c      write(6,*) ' ===>>> T = ',time,' Q_2_dn = ',qb2dn(1,1)
c      write(6,*) ' ===>>> T = ',time,' DQ_2_dn = ',dqb2dn(1,1)
c      write(6,*) ' ===>>> T = ',time,' Q_2_up = ',qb2up(1,1)
c      write(6,*) ' ===>>> T = ',time,' DQ_2_up = ',dqb2up(1,1)
      
      return
      end
