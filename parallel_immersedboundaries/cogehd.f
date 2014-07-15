c***********************************************************************
c     particular care is necessary for the metric quantities in the
c     non-linear terms
c***********************************************************************

c      subroutine hdnl1(q1,q2,q3,dq)         !PARALL
      subroutine hdnl1
      include 'param.f'
c
c     h term for the q1 momentum equation at i,j+1/2,k+1/2
c
      udx1 = dx1*0.25
      do 10 kc=1,n3m
      kmm=kmv(kc)
      kpp=kpv(kc)
      kp=kc+1
      do 10 jc=1,n2m
      jmm=jmv(jc)
      jpp=jpv(jc)
      jp=jc+1
      jm=jc-1
c     udx2=dx2/g2m(jc)*0.25 
      do 10 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
c
c    q1 q1 term
c
c
c             1   d  q_t q_t
c            --- -----------
c             r   d   t
c
      h11=( (q1(ip,jc,kc)+q1(ic,jc,kc))*(q1(ip,jc,kc)+q1(ic,jc,kc))
     %     -(q1(im,jc,kc)+q1(ic,jc,kc))*(q1(im,jc,kc)+q1(ic,jc,kc))
     %    )*udx1
C
c
c
c    q1 q2 term
c
c             1   d  q_t q_r     q_t q_r
c            --- -----------  +  --------    jc > 2
c             r   d   r            r^2
c
c
      q1nn=(q1(ic,jpp,kc)+q1(ic,jc,kc))
      q1ss=(q1(ic,jc,kc)+q1(ic,jmm,kc))
      h12=( (q2(ic,jp,kc)+q2(im,jp,kc))*q1nn
     %      -(q2(ic,jc,kc)+q2(im,jc,kc))*q1ss
     %     )*udx2m(jc)*0.25
c
      hq1=h11+h12

c     Non-inertial term
C     Modifica nuova versione NS
c      coriol=du0x+1/Ros*(dwy*(x3m(kc)-x30)-dwz*(x2m(jc)-x20))+
c     %2/Ros*(wy*(q3(ic,jc,kc)+q3(im,jc,kc)+q3(ic,jc,kp)+q3(im,jc,kp))
c     %*.25
c     %-wz*(q2(ic,jc,kc)+q2(im,jc,kc)+q2(ic,jp,kc)+q2(im,jp,kc))*.25)+
c     %1/Ros**2*(wx*wy*(x2m(jc)-x20)-wy**2*(x1c(ic)-x10)+
c     %wx*wz*(x3m(kc)-x30)-wz**2*(x1c(ic)-x10))
      coriol = 0.
C     Modifica nuova versione NS
c
c
      dq(ic,jc,kc)=-hq1-coriol
   10 continue
c
c   add the q1 q3 term
c
c
c                 d  q_t q_x 
c                -----------
c                 d   x      
c
      do 12 kc=2,n3m-1
c     udx3 = dx3/g3m(kc)*0.25
      kmm=kmv(kc)
      kpp=kpv(kc)
      kp=kc+1
      do 12 jc=1,n2m
      do 12 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
      jmm=jmv(jc)
      jpp=jpv(jc)
      h13=((q3(ic,jc,kp)+q3(im,jc,kp))*(q1(ic,jc,kpp)+q1(ic,jc,kc))
     %    -(q3(ic,jc,kc)+q3(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jc,kmm))
     %    )*udx3m(kc)*0.25

      coriol =  
     $   + ( - uhy*cos(alphax) + wx*( x3m(kc)-x30 ) )* 
     $     ( q1(ic,jpp,kc) - q1(ic,jmm,kc) )*0.5*udx2m(jc)
     $   + ( uhy*sin(alphax) - wx*( x2m(jc)-x20 ) )* 
     $     ( q1(ic,jc,kp) - q1(ic,jc,kmm) )*udx3m(kc)*0.5

      dq(ic,jc,kc)=dq(ic,jc,kc)-h13-coriol

   12 continue
c
c  OUTFLOW
c
      kc=n3m
c     ddx3 = dx3*0.5/g3m(kc)
      km=kc-1
      do jc=1,n2m
      do ic=1,n1m
      im=imv(ic)
      jmm=jmv(jc)
      jpp=jpv(jc)
      h13=(qb1n(ic,jc)*(qb3n(ic,jc)+qb3n(im,jc)) 
     1    -(q3(ic,jc,kc)+q3(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jc,km))
     1    *0.5)*udx3m(kc)*0.5

      coriol =  
     $   + ( - uhy*cos(alphax) + wx*( x3m(kc)-x30 ) )* 
     $     ( q1(ic,jpp,kc) - q1(ic,jmm,kc) )*0.5*udx2m(jc)
     $   + ( uhy*sin(alphax) - wx*( x2m(jc)-x20 ) )* 
     $     ( qb1n(ic,jc) - ( q1(ic,jc,kc) + q1(ic,jc,km))*0.5 )
     $     *udx3m(kc)*0.5

      dq(ic,jc,kc)=dq(ic,jc,kc)-h13-coriol
      enddo
      enddo
c
c INFLOW
c
      kc=1
      kp=kc+1
c     ddx3 = dx3*0.5/g3m(kc)
      do jc=1,n2m
      do ic=1,n1m
      im=imv(ic)
      jmm=jmv(jc)
      jpp=jpv(jc)
      h13=-(qb1s(ic,jc)*(qb3s(ic,jc)+qb3s(im,jc)) 
     1    -(q3(ic,jc,kp)+q3(im,jc,kp))*(q1(ic,jc,kc)+q1(ic,jc,kp))
     1    *0.5)*udx3m(kc)      *0.5

      coriol =  
     $   + ( - uhy*cos(alphax) + wx*( x3m(kc)-x30 ) )* 
     $     ( q1(ic,jpp,kc) - q1(ic,jmm,kc) )*0.5*udx2m(jc)
     $   + ( uhy*sin(alphax) - wx*( x2m(jc)-x20 ) )* 
     $     ( (q1(ic,jc,kp) + q1(ic,jc,kc))*0.5 - qb1s(ic,jc) )
     $     *udx3m(kc)*0.5

      dq(ic,jc,kc)=dq(ic,jc,kc)-h13-coriol
      enddo
      enddo
c
c  LES TERMS
c
c    COMMENTATI X PARALL
c
c      if(iles.gt.0) then
c      do kc=1,n3m
c        km = kmv(kc)
c        kp = kpv(kc)
c        do jc=1,n2m
c          jm = jmv(jc)
c          jp = jpv(jc)
c          do ic=1,n1m
c            im = imv(ic)
c            ip = ipv(ic)
c            term1 = (
c     %                   (st(ic,jc,kc,6)+st(ic,jp,kc,6)
c     %                   +st(im,jc,kc,6)+st(im,jp,kc,6))
c     %                  -(st(ic,jc,kc,6)+st(ic,jm,kc,6)
c     %                   +st(im,jc,kc,6)+st(im,jm,kc,6))
c     %                  )*udx2m(jc)*0.25
c            term2 = (
c     %                    st(ic,jc,kc,1)
c     %                 -  st(im,jc,kc,1)
c     %                  )*dx1
c            term3 = (
c     %                   (st(ic,jc,kc,5)+st(ic,jc,kp,5)
c     %                   +st(im,jc,kc,5)+st(im,jc,kp,5))
c     %                 - (st(ic,jc,kc,5)+st(ic,jc,km,5)
c     %                   +st(im,jc,kc,5)+st(im,jc,km,5))
c     %                  )*udx3m(kc)*0.25
c            dq(ic,jc,kc)=dq(ic,jc,kc)+(term1+term2+term3)
c          end do
c        end do
c      end do
c      end if

      return
      end
c
c***********************************************************************
c
c      subroutine hdnl2(q1,q2,q3)   ! PARALL
      subroutine hdnl2
      include 'param.f'
c
c     h term for the q2 momentum equation at i+1/2,j,k+1/2
c
      udx1 = dx1*0.25

      do 20 kc=1,n3m
      kmm=kmv(kc)
      kpp=kpv(kc)
      kp=kc+1
      do 20 jc=2,n2m
      jm=jc-1
      jp=jc+1
c     udx2 = 1./g2c(jc)*dx2*0.25
      do 20 ic=1,n1m
      imm=imv(ic)
      ipp=ipv(ic)
c
c     q2 q1 term
c
c
c             1   d  q_t q_r
c            --- -----------
c             r   d   t
c
      h21=((q1(ipp,jc,kc)+q1(ipp,jm,kc))
     %     *(q2(ipp,jc,kc)+q2(ic,jc,kc))
     %     -(q1(ic,jc,kc)+q1(ic,jm,kc))
     %     *(q2(ic,jc,kc)+q2(imm,jc,kc)))*udx1
C
c
c
c     q2 q2 term
c
c
c                 d  q_r q_r / r
c                ---------------
c                 d   r      
c
      h22=( (q2(ic,jp,kc)+q2(ic,jc,kc))
     %     *(q2(ic,jp,kc)+q2(ic,jc,kc))
     %     -(q2(ic,jc,kc)+q2(ic,jm,kc))
     %     *(q2(ic,jc,kc)+q2(ic,jm,kc))
     %    )*udx2c(jc)*0.25
c
c     non-inertial term
c
C     Modifica nuova versione NS
c      coriol=du0y+1./Ros*(-dwx*(x3m(kc)-x30)+dwz*(x1m(ic)-x10))+
c     %2./Ros*(-wx*(q3(ic,jc,kc)+q3(ic,jm,kc)+q3(ic,jc,kp)+q3(ic,jm,kp))
c     %*.25
c     %+wz*(q1(ic,jc,kc)+q1(ic,jm,kc)+q1(ipp,jc,kc)+q1(ipp,jm,kc))*.25)
c     %+1./Ros**2*(-wx**2.*(x2c(jc)-x20)+wx*wy*(x1m(ic)-x10)+
c     %wy*wz*(x3m(kc)-x30)
c     %-wz**2.*(x2c(jc)-x20))
C     Modifica nuova versione NS
      dph(ic,jc,kc)=-(h21+h22)
C     -coriol      Modifica nuova versione NS
c
   20 continue
c
c
c     add q2 q3 term
c
c
c                 d  q_x q_r 
c                -----------
c                 d   x      
c
      do 22 kc=2,n3m-1
c     udx3 = dx3/g3m(kc)*0.25
      kmm=kmv(kc)
      kpp=kpv(kc)
      kp=kc+1
      do 22 jc=2,n2m
      jm = jmv(jc)
      jp = jpv(jc)
      do 22 ic=1,n1m
      imm=imv(ic)
      ipp=ipv(ic)
      h23=((q3(ic,jc,kp)+q3(ic,jm,kp))*(q2(ic,jc,kpp)+q2(ic,jc,kc))
     %    -(q3(ic,jc,kc)+q3(ic,jm,kc))*(q2(ic,jc,kc)+q2(ic,jc,kmm))
     %    )*udx3m(kc)*0.25
c
C     Modifica nuova versione NS
      coriol = 
c     wx*uhy*sin(alphax) - wx**2.*(x2c(jc)-x20)
     $   - ( q3(ic,jc,kc)+q3(ic,jm,kc)+q3(ic,jc,kp)+q3(ic,jm,kp) )
     $     *wx*0.25
     $   + ( - uhy*cos(alphax) + wx*( x3m(kc)-x30 ) )* 
     $     ( q2(ic,jp,kc) - q2(ic,jm,kc) )*0.5*udx2c(jc)
     $   + ( uhy*sin(alphax) - wx*( x2c(jc)-x20 ) )* 
     $     ( q2(ic,jc,kp) - q2(ic,jc,kmm) )*udx3m(kc)*0.5
C     Modifica nuova versione NS
      
      dph(ic,jc,kc)= dph(ic,jc,kc) -h23-coriol

   22 continue

c
c OUTFLOW
c
      kc=n3m
      km=kc-1
c     ddx3 = dx3*0.5/g3m(kc)
      do jc=2,n2m
      do ic=1,n1m
      jm=jmv(jc)
      jp = jpv(jc)
      h23=( qb2n(ic,jc)*(qb3n(ic,jc)+qb3n(ic,jm))
     1    -(q3(ic,jc,kc)+q3(ic,jm,kc))*(q2(ic,jc,kc)+q2(ic,jc,km))
     1    *0.5 )*udx3m(kc)*0.5

C     Modifica nuova versione NS
      coriol = 
c     wx*uhy*sin(alphax) - wx**2.*(x2c(jc)-x20)
     $   - ( q3(ic,jc,kc)+q3(ic,jm,kc)+qb3n(ic,jc)+qb3n(ic,jm) )
     $     *wx*0.25
     $   + ( - uhy*cos(alphax) + wx*( x3m(kc)-x30 ) )* 
     $     ( q2(ic,jp,kc) - q2(ic,jm,kc) )*0.5*udx2c(jc)
     $   + ( uhy*sin(alphax) - wx*( x2c(jc)-x20 ) )* 
     $     ( qb2n(ic,jc) - ( q2(ic,jc,km)+q2(ic,jc,kc) )*0.5 )
     $     *udx3m(kc)
C     Modifica nuova versione NS

      dph(ic,jc,kc)=dph(ic,jc,kc)-h23-coriol

      enddo
      enddo
c
c INFLOW
c       
      kc=1
      kp=kc+1
c     ddx3 = dx3*0.5/g3m(kc)
      do jc=2,n2m
      jm=jmv(jc)
      jp = jpv(jc)
      do ic=1,n1m
      h23=
     $    ( (q3(ic,jc,kp)+q3(ic,jm,kp))*(q2(ic,jc,kc)+q2(ic,jc,kp))*0.5
     $     -qb2s(ic,jc)*(qb3s(ic,jc)+qb3s(ic,jm)))*udx3m(kc)*0.5

C     Modifica nuova versione NS
      coriol = 
c     wx*uhy*sin(alphax) - wx**2.*(x2c(jc)-x20)
     $   - ( qb3s(ic,jc)+qb3s(ic,jm)+q3(ic,jc,kp)+q3(ic,jm,kp) )
     $     *wx*0.25
     $   + ( - uhy*cos(alphax) + wx*( x3m(kc)-x30 ) )* 
     $     ( q2(ic,jp,kc) - q2(ic,jm,kc) )*0.5*udx2c(jc)
     $   + ( uhy*sin(alphax) - wx*( x2c(jc)-x20 ) )* 
     $     ( ( q2(ic,jc,kp) + q2(ic,jc,kc) )*0.5 - qb2s(ic,jc) )
     $     *udx3m(kc)
C     Modifica nuova versione NS

      dph(ic,jc,kc)=dph(ic,jc,kc)-h23-coriol

      enddo    
      enddo    
c
c  LES TERMS
c
c    COMMENTATI X PARALL
c
c      if(iles.gt.0) then
c      udx1 = dx1*0.25
c      do kc=1,n3m
c        km = kmv(kc)
c        kp = kpv(kc)
cc       udx3 = dx3/g3m(kc)
c        do jc=2,n2m
cc         udx2 = dx2/g2c(jc)
c          jm = jmv(jc)
c          jp = jpv(jc)
c          do ic=1,n1m
c            im = imv(ic)
c            ip = ipv(ic)
c            term1 = (
c     %                st(ic,jc,kc,2)
c     %               -st(ic,jm,kc,2)
c     %                  )*udx2c(jc)
c            term2 = (
c     %                   (st(ic,jc,kc,6)+st(ip,jc,kc,6)
c     %                   +st(ic,jm,kc,6)+st(ip,jm,kc,6))
c     %                 - (st(ic,jc,kc,6)+st(im,jc,kc,6)
c     %                   +st(ic,jm,kc,6)+st(im,jm,kc,6))
c     %                  )*udx1
c            term3 = (
c     %                   (st(ic,jc,kc,4)+st(ic,jc,kp,4)
c     %                   +st(ic,jm,kc,4)+st(ic,jm,kp,4))
c     %                 - (st(ic,jc,kc,4)+st(ic,jc,km,4)
c     %                   +st(ic,jm,kc,4)+st(ic,jm,km,4))
c     %                  )*udx3m(kc)*0.25
c            dph(ic,jc,kc)=dph(ic,jc,kc)+
c     %                   (term1+term2+term3)
c          end do
c        end do
c      end do
c      end if
c     
      return
      end
c
c***********************************************************************
c
c      subroutine hdnl3(q1,q2,q3)            ! PARALL
      subroutine hdnl3
      include 'param.f'
c
c     h term for the q3 momentum equation at i+1/2,j+1/2,k
c
      udx1=dx1*0.25

      do 30 kc=2,n3m
c     udx3 = dx3/g3c(kc)*0.25
      km=kmv(kc)
      kp=kc+1
      kmm=kc-1
      do 30 jc=1,n2m
      jp=jc+1
      jmm=jmv(jc)
      jpp=jpv(jc)
c      if (jc.eq.n2m)jpp=jc !!!! prova c.c. up
c     udx2=dx2/g2m(jc)*0.25
      do 30 ic=1,n1m
      imm=imv(ic)
      ipp=ipv(ic)
c
c    q3 q1 term
c
c
c            1    d  q_x q_t
c           ---  -----------
c            r    d   t
c
      h31=((q1(ipp,jc,kc)+q1(ipp,jc,km))
     %        *(q3(ipp,jc,kc)+q3(ic,jc,kc))
     %        -(q1(ic,jc,kc)+q1(ic,jc,km))
     %        *(q3(ic,jc,kc)+q3(imm,jc,kc)))*udx1

C
c
c
c    q3 q2 term
c
c
c            1    d  q_x q_r 
c           ---  -----------
c            r    d   r      
c
      h32=(((q2(ic,jp,kc)+q2(ic,jp,km))
     %     *(q3(ic,jpp,kc)+q3(ic,jc,kc)))
     %    -((q2(ic,jc,kc)+q2(ic,jc,km))
     %     *(q3(ic,jc,kc)+q3(ic,jmm,kc))))*udx2m(jc)*0.25
c      if (jc.eq.n2m)h32=( ((q2(ic,jp,kc)+q2(ic,jp,km))
c     %     *qb3up(ic,kc) )*0.5
c     %    -((q2(ic,jc,kc)+q2(ic,jc,km))
c     %     *(q3(ic,jc,kc)+q3(ic,jmm,kc)))*0.25 )*udx2m(jc)
c      if (jc.eq.1)h32=( ((q2(ic,jp,kc)+q2(ic,jp,km))
c     %     *(qb3up(ic,kc)+q3(ic,jc,kc)))*0.25
c     %    -((q2(ic,jc,kc)+q2(ic,jc,km))
c     %     *qb3dn(ic,kc))*0.5 )*udx2m(jc)        prova c.c. up & dn
c
c    q3 q3 term
c
c
c                 d  q_x q_x 
c                -----------
c                 d   x      
c
      h33=((q3(ic,jc,kp)+q3(ic,jc,kc))*(q3(ic,jc,kp)+q3(ic,jc,kc))
     %    -(q3(ic,jc,kc)+q3(ic,jc,kmm))*(q3(ic,jc,kc)+q3(ic,jc,kmm))
     %    )*udx3c(kc)*0.25
c
c
c     inertial term

C     Modifica nuova versione NS
C      coriol=du0z+1./Ros*(dwx*(x2m(jc)-x20)-dwy*(x1m(ic)-x10))+2./Ros*
C     %(wx*(q2(ic,jc,kc)+q2(ic,jc,kmm)+q2(ic,jpp,kc)+q2(ic,jpp,kmm))
C     %*.25
C     %-wy*(q1(ic,jc,kc)+q1(ic,jc,kmm)+q1(ipp,jc,kc)+q1(ipp,jc,kmm))
C     %*.25)
C     %+1./Ros**2*(-wx**2*(x3c(kc)-x30)+wx*wz*(x1m(ic)-x10)-
C     %wy**2*(x3c(kc)-x30)+wy*wz*(x2m(jc)-x20))
C     Modifica nuova versione NS
      coriol = 
c     wx*uhy*sin(alphax) - wx**2.*(x2c(jc)-x20)
     $   + ( q2(ic,jc,kc)+q2(ic,jp,kc)+q2(ic,jc,km)+q2(ic,jp,km) )
     $     *wx*0.25
     $   + ( - uhy*cos(alphax) + wx*( x3c(kc)-x30 ) )* 
     $     ( q3(ic,jp,kc) - q3(ic,jmm,kc) )*0.5*udx2m(jc)
     $   + ( uhy*sin(alphax) - wx*( x2m(jc)-x20 ) )* 
     $     ( q3(ic,jc,kp) - q3(ic,jc,km) )*udx3c(kc)*0.5
C     Modifica nuova versione NS

      qcap(ic,jc,kc)=-(h31+h32+h33)-coriol

   30 continue

c
c  LES TERMS
c
c  L.E.S. commentata PER PARALL
c
c      if(iles.gt.0) then
c      udx1=dx1*0.25
c      do kc=2,n3m
c        km = kmv(kc)
c        kp = kpv(kc)
c        do jc=1,n2m
c          jm = jmv(jc)
c          jp = jpv(jc)
c          do ic=1,n1m
c            im = imv(ic)
c            ip = ipv(ic)
c            term1 = (
c     %                (st(ic,jc,kc,4)+st(ic,jp,kc,4)
c     %                +st(ic,jc,km,4)+st(ic,jp,km,4))
c     %               -(st(ic,jc,kc,4)+st(ic,jm,kc,4)
c     %                +st(ic,jc,km,4)+st(ic,jm,km,4))
c     %                  )*udx2m(jc)*0.25
c            term2 = (
c     %                   (st(ic,jc,kc,5)+st(ip,jc,kc,5)
c     %                   +st(ic,jc,km,5)+st(ip,jc,km,5))
c     %                 - (st(ic,jc,kc,5)+st(im,jc,kc,5)
c     %                   +st(ic,jc,km,5)+st(im,jc,km,5))
c     %                  )*udx1
c            term3 = (
c     %                    st(ic,jc,kc,3)
c     %                 -  st(ic,jc,km,3)
c     %                  )*udx3c(kc)
c            qcap(ic,jc,kc)=qcap(ic,jc,kc)+(term1+term2+term3)
c          end do
c        end do
c      end do
c      end if
c     
      return
      end
c
c*************************************************************************
c
      subroutine rotation (tin)
      include 'param.f'

      pi=2.*asin(1.)
c     Sigma e' l' inverso del numero di strouhal      
      sigma=1.

      amp=amplit            !pari a c*/A*
      
      ttilde = pi/(amp*8.)
c      if (tin.lt.ttilde)then
c         uzero = uzero_max*( 1. - cos(8.*amp*tin) )*0.5
c      else 
      uzero = uzero_max
c      endif

      do k = 1,n3
         do j = 1,n2
           do i = 1,n1
             uinfth(i,j) = uzero
c               uinfth(i,j) = uzero*( 1. +
c     $          cos(tin)*(cos( pi*x2c(j)/2.) + 1. )/4. 
c     $                              )
            enddo
         enddo
      enddo

      alphaymax=0.
      alphaxmax=thetamax*pi/180
            
      hy = x20 + sin(amp*tin)/amp
      uhy = cos(amp*tin)        !h_punto
      duhy = -amp*sin(amp*tin)  !h_2punto
c      hy = 0.
c      uhy = 0.
c      duhy = 0.
      uhx=0.
      duhx=0.

c    velocita' di inflow 
c   ATTENZIONE E' NEGATIVA POICHE' E' LA VELOCITA' DEL CORPO
      uinflow=0.
      
      alphax=alphaxmax*cos(amp*tin)
      wx=-alphaxmax*amp*sin(amp*tin)
      dwx=-alphaxmax*amp**2*cos(amp*tin)

      alphay=-alphaymax*cos(sigma*tin)
      wy=alphaymax*sigma*sin(sigma*tin)
      dwy=alphaymax*sigma**2.*cos(sigma*tin)

      wz=0.
      dwz=0.

      u0x=uhx*cos(alphay)
C      u0x=uhx*cos(alphay)-uinflow*sin(alphay)
c      du0x=duhx*cos(alphay)-wy*uhx*sin(alphay)-wy*uinflow*cos(alphay)

      u0y=+uhy*cos(alphax)
c      du0y=+duhy*cos(alphax)

      u0z=uhx*sin(alphay) - uhy*sin(alphax)
C      u0z=uhx*sin(alphay)+uinflow*cos(alphay)-uhy*sin(alphax)
c      du0z=duhx*sin(alphay)-uinflow*wy*sin(alphay)
c     %     -duhy*sin(alphax)
c     Centro di rotazione
      x10=0
      x20=0
      x30=-0.16667
      if (myid.eq.0)write(67,1001)tin,hy,uhy,duhy,alphax,wx,dwx,uzero 
 1001 format(8(2x,e12.5))

c     write(*,*)' Posizione centro di rotazione: ',x10,x20,x30
      return
      end

