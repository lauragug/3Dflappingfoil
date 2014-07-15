
************************************************************************
c
c     in this subroutine the indices ip,im,jp,jm,kp,km are calculated.
c     the indices jpc,jup,jmc,jum,kpc,kup,kmc,kum are indices 
c     for the walls
c
      subroutine indic
      include'param.f'
c
c
c   azimuthal direction
c
      do 1 ic=1,n1m
        ip=ic+1
        imv(ic)=ic-1
        if(ic.eq.1) imv(ic)=n1m
        ipv(ic)=ic+1
        if(ic.eq.n1m) ipv(ic)=1
    1 continue
c
c    symmetrical index for the azimuthal location
c    i + pi
c
      do i=1,n1m
       isym(i) = i + n1m/2
       if(isym(i).gt.n1m) isym(i) = isym(i) - n1m
      enddo
c
c   streamwise direction
c
      do 4 kc=1,n3m
        kmv(kc)=kc-1
        kpv(kc)=kc+1
        if(kc.eq.1) kmv(kc)=kc
        if(kc.eq.n3m) kpv(kc)=kc
    4 continue
c
c     direction normal to non-slip walls
c
      do 150 kc=1,n3m
        kpc(kc)=kpv(kc)-kc
        kmc(kc)=kc-kmv(kc)
        kup(kc)=1-kpc(kc)
        kum(kc)=1-kmc(kc)
  150 continue
c
c   indices for radial direction
c
      do 3 jc=1,n2
        jmv(jc)=jc-1
        jpv(jc)=jc+1
        if(jc.eq.1) jmv(jc)=jc
        if(jc.eq.n2m) jpv(jc)=jc
        if(jc.eq.n2) jpv(jc)=jc
    3 continue
c
c   indices for radial direction
c
      do 15 jc=1,n2m
        jpc(jc)=jpv(jc)-jc
        jmc(jc)=jc-jmv(jc)
        jup(jc)=1-jpc(jc)
        jum(jc)=1-jmc(jc)
   15 continue
c
      return
      end
c
************************************************************************
c
c     this subroutine calculates divg(q).
c
      subroutine divg(q1,q2,q3)
      include'param.f'
      REAL q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)
c   
      usdtal = 1./(dt*al)
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,ipv)
!$OMP$ SHARED(q1,q2,q3,qcap,dx1,udx2m,udx3m,usdtal)
!$OMP$ PRIVATE(kp,jp,ip,dqcap)

      do 11 kc=1,n3m
        kp=kc+1
        do 11 jc=1,n2m
          jp=jc+1
            do 11 ic=1,n1m
              ip=ipv(ic)
              dqcap= (q1(ip,jc,kc)-q1(ic,jc,kc))*dx1
     %              +(q2(ic,jp,kc)-q2(ic,jc,kc))*udx2m(jc)
     %              +(q3(ic,jc,kp)-q3(ic,jc,kc))*udx3m(kc)
              qcap(ic,jc,kc)=dqcap*usdtal
   11 continue
!$OMP  END PARALLEL DO
      return
      end
c
c
c************************************************************************
c
c     this subroutine calculates the solenoidal vel field.
c       q(n+1)=qhat-grad(dph)*dt ,  pr=dph
c    third order runge-kutta is used.
c
      subroutine updvp(q1,q2,q3)
      include'param.f'
      REAL q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)
c
c  ***********  compute the q1 velocity component
c               v1dgf=component 1 of grad(dph)
c   
c
      usurm1 = al*dt*dx1
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,imv,udx2c)
!$OMP$ SHARED(q1,q2,q3,dph,usurm1,al,dt,n2)
!$OMP$ PRIVATE(im,jm,km,usurm2)
      do kc=1,n3m
        do 1 jc=1,n2m
          do 1 ic=1,n1m
            im=imv(ic)
            q1(ic,jc,kc)=q1(ic,jc,kc)-
     %      (dph(ic,jc,kc)-dph(im,jc,kc))*usurm1
    1 continue
c
c
c  ***********  compute the q2 velocity component
c
        do 2 jc=2,n2m
          jm=jc-1
          usurm2 = al*dt*udx2c(jc)
          do 2 ic=1,n1m
            q2(ic,jc,kc)=q2(ic,jc,kc)-
     %      (dph(ic,jc,kc)-dph(ic,jm,kc))*usurm2
    2 continue
        do 21 ic=1,n1m
c         q2(ic,1,kc)=0.
c         q2(ic,n2,kc)=0.
          q2(ic,1,kc)=qb2dn(ic,kc)
          q2(ic,n2,kc)=qb2up(ic,kc)
   21 continue
      end do
!$OMP END PARALLEL DO
c
c  ***********  compute the q3 velocity component
c               q3 is the cartesian component
c               v3dgf=component 3 of grad(dph)
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,udx3c)
!$OMP$ SHARED(q3,dph,al,dt)
!$OMP$ PRIVATE(km,usurm,usurm2)
      do  kc=2,n3m
        usurm = al*dt*udx3c(kc)
        km=kc-1
        do 5 jc=1,n2m
          do 5 ic=1,n1m
            q3(ic,jc,kc)=q3(ic,jc,kc)-
     %      (dph(ic,jc,kc)-dph(ic,jc,km))*usurm
    5 continue
      do 7 ic=1,n1m
c         q3(ic,1,kc)=qb3dn(ic,kc) !LAURA nuova versione NS
         q3(ic,n2,kc)=qb3up(ic,kc)
    7 continue
      enddo
      do jc=1,n2m
       do ic=1,n1m
          q3(ic,jc,1)=qb3s(ic,jc)
          q3(ic,jc,n3)=qb3n(ic,jc)
       enddo
      enddo
!$OMP  END PARALLEL  DO
      return
      end
c
c
c
c
************************************************************************
c   this subroutine performs the calculation of the pressure.
c   this depends on the fractional step
c
      subroutine prcalc(pr)
      include'param.f'
      REAL   pr(m1,m2,m3)
c
c    the pressure is evaluated at the center of the box.
c
c     p^{n+1} = p^{n} + phi^{n+1} - b * Nabla^2 phi^{n+1}
c
      be=al*beta
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,kpv,kmv,jmv,jpv,imv,ipv)
!$OMP$ SHARED(pr,dph,be,dx1q,apphk,acphkk,amphk,apphj,acphj,amphj,visc)
!$OMP$ PRIVATE(kp,km,jm,jp,im,ip)

      do 1 kc=1,n3m
        kp=kpv(kc)
        km=kmv(kc)
        do 1 jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
            do 1 ic=1,n1m
              ip=ipv(ic)
              im=imv(ic)
              pr(ic,jc,kc)=pr(ic,jc,kc)+dph(ic,jc,kc)
     %   -be*(
     %        (dph(ip,jc,kc)-2.*dph(ic,jc,kc)+dph(im,jc,kc))*dx1q + 
     %        (dph(ic,jc,kp)*apphk(kc)
     %        +dph(ic,jc,kc)*acphkk(1,kc)
     %        +dph(ic,jc,km)*amphk(kc)) +
     %        (dph(ic,jp,kc)*apphj(jc)
     %        +dph(ic,jc,kc)*acphj(jc)
     %        +dph(ic,jm,kc)*amphj(jc))           )*visc(kc)
    1 continue
!$OMP  END PARALLEL DO
c
      return
      end
c
c
************************************************************************
c
c   this subroutine performs the inversion of the q1 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q1
c   are treated implicitly
c         direction x1
c
      subroutine solq1i(q1)
      include'param.f'
      REAL q1(m1,m2,m3)
      dimension amil(m2,m1),acil(m2,m1),apil(m2,m1),fil(m2,m1)
c

!$OMP  PARALLEL DO
!$OMP$ SHARED(n3,n2,n1,forclo)
      do k=1,n3
        do j=1,n2
          do i=1,n1
            forclo(i,j,k)=1.
          end do
        end do
      end do
!$OMP  END PARALLEL DO

      if(infig.ne.-1) then
!$OMP  PARALLEL DO
!$OMP$ SHARED(indgeo,npunt,forclo)
       do n=1,npunt
          i=indgeo(1,n,1)
          j=indgeo(1,n,2)
          k=indgeo(1,n,3)
          forclo(i,j,k)= 0.
       end do
!$OMP  END PARALLEL DO
      end if

c
      betadx=beta*dx1q*al

!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,betadx,forclo,rhs)
!$OMP$ PRIVATE(amil,apil,acil,fil,forcl,iadd)

      do k=1,n3m
        do j=1,n2m
          do i=1,n1m
            forcl = betadx*forclo(i,j,k)
            amil(j,i)=-forcl
            acil(j,i)=1.+2.*forcl
            apil(j,i)=-forcl
            iadd = i+(j-1)*n1m+(k-1)*n1m*n2m
            fil(j,i)=rhs(iadd)
          end do
        end do
c
      call trvpki(amil,acil,apil,fil,1,n1m,1,n2m,m1,m2)
c
        do j=1,n2m
          do i=1,n1m
            iadd = i+(j-1)*n1m+(k-1)*n1m*n2m
            rhs(iadd)=fil(j,i)
          end do
        end do
      end do
!$OMP  END PARALLEL DO


c
c
      return
      end
c
************************************************************************
c
c   this subroutine performs the inversion of the q1 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q1
c   are treated implicitly
c
      subroutine solq1j
      include'param.f'
      dimension amjl(m1,m2),apjl(m1,m2),acjl(m1,m2),fjl(m1,m2)
c
c  ************ compute dq1  sweeping along the x2 direction
c               wall boundaries interior points
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3,n2,n1,forclo)
      do k=1,n3
        do j=1,n2
          do i=1,n1
            forclo(i,j,k)=1.
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
       if(infig.ne.-1) then
!$OMP  PARALLEL DO
!$OMP$ SHARED(indgeo,npunt,forclo)
      do n=1,npunt
          i=indgeo(1,n,1)
          j=indgeo(1,n,2)
          k=indgeo(1,n,3)
          forclo(i,j,k)= 0.
      end do
!$OMP  END PARALLEL DO
      end if
c
      betadx=beta*al

!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,betadx,forclo,rhs)
!$OMP$ PRIVATE(amjl,apjl,acjl,fjl,fclo,iadd)

       do k=1,n3m
        do j=1,n2m
          do i=1,n1m
            fclo=forclo(i,j,k)*betadx
            amjl(i,j)=-am1j(j)*fclo
            acjl(i,j)=1.-ac1j(j)*fclo
            apjl(i,j)=-ap1j(j)*fclo
            iadd=i+(j-1)*n1m+(k-1)*n1m*n2m
            fjl(i,j)=rhs(iadd)
          end do
        end do

c
      call trvb(amjl,acjl,apjl,fjl,m1,m2,n2m,n1m)
c
      do j=1,n2m
        do i=1,n1m
          iadd=i+(j-1)*n1m+(k-1)*n1m*n2m
          rhs(iadd) = fjl(i,j)
        end do
       end do
      end do
!$OMP  END PARALLEL DO

      return
      end
c
************************************************************************
c   this subroutine performs the inversion of the q2 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q2
c   are treated implicitly
c   in the first part the rhs is calculated
c        direction x1
c
      subroutine solq2i(q1)
      include'param.f'
      REAL q1(m1,m2,m3)
      dimension amil(m2,m1),acil(m2,m1),apil(m2,m1),fil(m2,m1)
c  ************ compute dq2** sweeping along the x1 direction
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3,n2,n1,forclo)
      do k=1,n3
        do j=1,n2
          do i=1,n1
            forclo(i,j,k)=1.
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
       if(infig.ne.-1) then
!$OMP  PARALLEL DO
!$OMP$ SHARED(indgeo,npunr,forclo)
      do n=1,npunr
        i=indgeo(2,n,1)
        j=indgeo(2,n,2)
        k=indgeo(2,n,3)
        forclo(i,j,k)= 0.
      end do
!$OMP  END PARALLEL DO
      end if
c
      betadx=beta*dx1q*al
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n2,n1m,betadx,forclo,rhs)
!$OMP$ PRIVATE(amil,apil,acil,fil,forcl,iadd)
       do k=1,n3m
        do j=1,n2m
          do i=1,n1m
            forcl = forclo(i,j,k)*betadx
            amil(j,i)=-forcl
            acil(j,i)=1.+2.*forcl
            apil(j,i)=-forcl
            iadd = i+(j-1)*n1m+(k-1)*n1m*n2
            fil(j,i)=rhs(iadd)
          end do
        end do
c
      call trvpki(amil,acil,apil,fil,1,n1m,1,n2m,m1,m2)
c
        do j=1,n2m
          do i=1,n1m
            iadd = i+(j-1)*n1m+(k-1)*n1m*n2
            rhs(iadd)=fil(j,i)
          end do
        end do
      end do
!$OMP  END PARALLEL DO

c
      return
      end
c
************************************************************************
c   this subroutine performs the inversion of the q2 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q2
c   are treated implicitly
c   in the first part the rhs is calculated
c         direction x2
c
      subroutine solq2j
      include'param.f'
      dimension amjl(m1,m2),apjl(m1,m2),acjl(m1,m2),fjl(m1,m2)
c  ********* compute the dq2* sweeping in the x2 direction
c            wall boundaries direction
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3,n2,n1,forclo)
      do k=1,n3
        do j=1,n2
          do i=1,n1
            forclo(i,j,k)=1.
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
       if(infig.ne.-1) then
!$OMP  PARALLEL DO
!$OMP$ SHARED(indgeo,npunr,forclo)
      do n=1,npunr
        i=indgeo(2,n,1)
        j=indgeo(2,n,2)
        k=indgeo(2,n,3)
        forclo(i,j,k)= 0.
      end do
!$OMP  END PARALLEL DO
      end if
c
      betadx=beta*al

!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n2,n1m,betadx,forclo,rhs)
!$OMP$ PRIVATE(amjl,apjl,acjl,fjl,fclo,iadd)
      do k=1,n3m
        do j=2,n2m
          do i=1,n1m
            fclo=forclo(i,j,k)*betadx
            amjl(i,j)=-am2j(j)*fclo
            acjl(i,j)=1.-ac2j(j)*fclo
            apjl(i,j)=-ap2j(j)*fclo
            iadd=i+(j-1)*n1m+(k-1)*n1m*n2
            fjl(i,j)=rhs(iadd)
          end do
        end do
c j=1 and j=n2
        do i=1,n1m
         apjl(i,1)=0.
         acjl(i,1)=1.
         amjl(i,1)=0.
         fjl(i,1)=dqb2dn(i,k)
         apjl(i,n2)=0.
         acjl(i,n2)=1.
         amjl(i,n2)=0.
         fjl(i,n2)=dqb2up(i,k)
        end do

c
      call trvb(amjl,acjl,apjl,fjl,m1,m2,n2,n1m)
c
      do j=1,n2
        do i=1,n1m
          iadd=i+(j-1)*n1m+(k-1)*n1m*n2
          rhs(iadd) = fjl(i,j)
          end do
        end do
      end do
!$OMP  END PARALLEL DO

      return
      end

c
************************************************************************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c       direction x1
c
      subroutine solq3i(q1)
      include'param.f'
      REAL q1(m1,m2,m3)
      dimension amil(m2,m1),acil(m2,m1),apil(m2,m1),fil(m2,m1)
c
c  ************ compute dq3** sweeping along the x1 direction
c               inflow outflow
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3,n2,n1,forclo)
      do k=1,n3
        do j=1,n2
          do i=1,n1
            forclo(i,j,k)=1.
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
       if(infig.ne.-1) then
!$OMP  PARALLEL DO
!$OMP$ SHARED(indgeo,npunz,forclo)
      do n=1,npunz
        i=indgeo(3,n,1)
        j=indgeo(3,n,2)
        k=indgeo(3,n,3)
        forclo(i,j,k)= 0.
      end do
!$OMP  END PARALLEL DO
      end if
c
      betadx=beta*dx1q*al
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,betadx,forclo,rhs)
!$OMP$ PRIVATE(amil,apil,acil,fil,forcl,iadd)
      do k=1,n3m
        do j=1,n2m
          do i=1,n1m
            forcl = betadx*forclo(i,j,k)
            amil(j,i)=-forcl
            acil(j,i)=1.+2.*forcl
            apil(j,i)=-forcl
            iadd = i+(j-1)*n1m+(k-1)*n1m*n2m
            fil(j,i)=rhs(iadd)
          end do
        end do
c
      call trvpki(amil,acil,apil,fil,1,n1m,1,n2m,m1,m2)
c
        do j=1,n2m
          do i=1,n1m
            iadd = i+(j-1)*n1m+(k-1)*n1m*n2m
            rhs(iadd)=fil(j,i)
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
      return
      end
c
************************************************************************
c                       SUBROUTINE INVTR1
c   This subroutine performs the computation of Q~~ for the q1 momentum 
c   equation (azimuthal direction) by a factored implicit scheme.
c   Viscous terms are treated implicitly, nonlinear terms explicitly.
c   
c        ~~     n
c   dQ = Q  -  Q    
c
c         alp dt   
c   bet = -------
c          2 Re
c
c   The first equation in the introduction of THSCHEM then becomes
c   
c                
c              2      [            n         n       n-1   alp       2  n ]
c  (1-bet*nabla ) dQ =[-alp*grad (P ) + gam*H + rho*H   + ----- nabla (Q )]*dt
c                     [                                    Re             ]
c  
c                                                                      3 
c   The left hand side of this equation is then factored at an order dt 
c   as follows
c
c           2           2           2
c  (1 -bet*d_x)*(1-bet*d_r)*(1-bet*d_th) dQ = RHS
c
c           2     2       2
c  Where   d_x,  d_r and d_th  are the discrete differential operators of
c  the viscous terms in the axial, radial nad azimuthal directions
c  respectively. RHS is the righ hand side of the momentum equation as
c  is written above.            
c
      subroutine invtr1(q1,pr,ru1,dq,aldto)
      include'param.f'
      REAL   pr(m1,m2,m3)
      REAL ru1(m1,m2,m3)
      REAL q1(m1,m2,m3)
      REAL dq(m1,m2,m3)
      real etime,t(2)
      timinv=etime(t)
      alre=al/ren
c     rhs_1m = -10.
c
c     do kc=1,n3
c       do jc=1,n2
c         do ic=1,n1 
c           iadd=ic+(jc-1)*n1+(kc-1)*n1*n2
c           rhs(iadd) = 0.
c         end do
c       end do 
c     end do 
c
c
c
c  compute the rhs of the factored equation everything at i,j+1/2,k+1/2
c
      ugkk=4./3.
      udx1=dx1*al
      udx1q = dx1q
      n3mm=n3m-1
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3mm,kmv,kpv,n2m,jpv,jmv,n1m,ipv,imv)
!$OMP$ SHARED(q1,udx1q,ap1j,ac1j,am1j,ap1k,ac1k,am1k,pr,udx1)
!$OMP$ SHARED(ga,ro,dq,ru1,alre,visc,dt)
!$OMP$ PRIVATE(km,kp,jp,jm,ip,im,d11q1,d22q1,d33q1,dcq1,dpx11,iadd)

        do 25 kc=2,n3mm
          km=kmv(kc)
          kp=kpv(kc)
          do 25 jc=1,n2m
             jp=jpv(jc)
             jm=jmv(jc)
             do 25 ic=1,n1m
                ip=ipv(ic)
                im=imv(ic)
c
c   11 second deriv. of q1(n)
c
                d11q1=(q1(ip,jc,kc)-q1(ic,jc,kc)*2.+q1(im,jc,kc))*
     %            udx1q
c
c   22 second deriv. of q1(n)
c
               d22q1=q1(ic,jp,kc)*ap1j(jc)
     %              +q1(ic,jc,kc)*ac1j(jc)
     %              +q1(ic,jm,kc)*am1j(jc)
c
c   33 second deriv. of q1(n)
c
               d33q1=q1(ic,jc,kp)*ap1k(kc)
     %              +q1(ic,jc,kc)*ac1k(kc)
     %              +q1(ic,jc,km)*am1k(kc)
c 
c   all viscous terms
c
                dcq1 = d11q1 + d22q1 + d33q1
c
c   grad(pr) along 1
c
                dpx11=(pr(ic,jc,kc)-pr(im,jc,kc))*udx1
c
c
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c
c   evaluation of RHS
c
            rhs(iadd)=+(ga*dq(ic,jc,kc)+ro*ru1(ic,jc,kc)
     %                 +alre*dcq1*visc(kc)-dpx11)*dt
c
c   updating the nonlinear terms
c
                ru1(ic,jc,kc)=dq(ic,jc,kc)
   25 continue
!$OMP  END PARALLEL DO
c
c
c  including near horizonthal boundaries points.
c
c
c      LOWER WALL
c
        kc=1
        kp=kc+1
        codea = dx3q/(g3c(kp)*g3m(kc))
        codeb = dx3q/(g3c(kc)*g3m(kc))
        do 101 jc=1,n2m
           jp=jpv(jc)
           jm=jmv(jc)
           do 101 ic=1,n1m
              ip=ipv(ic)
              im=imv(ic)
c
c   11 second deriv. of q1(n)
c
              d11q1=(q1(ip,jc,kc)-q1(ic,jc,kc)*2.+q1(im,jc,kc))*
     %           udx1q
c
c   22 second deriv. of q1(n)
c
               d22q1=q1(ic,jp,kc)*ap1j(jc)
     %              +q1(ic,jc,kc)*ac1j(jc)
     %              +q1(ic,jm,kc)*am1j(jc)
c
c   33 second deriv. of q1(n)
c
             d33q1=((q1(ic,jc,kp)-q1(ic,jc,kc))*codea   - 
     %              (q1(ic,jc,kc)-qb1s(ic,jc) )*codeb*2. )
c
            dcq1 = d11q1 + d22q1 + d33q1
c
c   grad(pr) along 1
c
            dpx11=(pr(ic,jc,kc)-pr(im,jc,kc))*udx1
c
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c
c   evaluation of RHS
c
            rhs(iadd)=+(ga*dq(ic,jc,kc)+ro*ru1(ic,jc,kc)
     %                 +alre*dcq1*visc(kc)-dpx11)*dt
c
c
            ru1(ic,jc,kc)=dq(ic,jc,kc)
  101 continue
c
c     UPPER WALL
c
        kc=n3m
        kp=kc+1
        km=kc-1
        codea = dx3q/(g3c(kp)*g3m(kc))
        codeb = dx3q/(g3c(kc)*g3m(kc))
        do 102 jc=1,n2m
           jp=jpv(jc)
           jm=jmv(jc)
           do 102 ic=1,n1m
              ip=ipv(ic)
              im=imv(ic)
c
c   11 second deriv. of q1(n)
c
              d11q1=(q1(ip,jc,kc)-q1(ic,jc,kc)*2.+q1(im,jc,kc))*
     %           udx1q
c
c   22 second deriv. of q1(n)
c
               d22q1=q1(ic,jp,kc)*ap1j(jc)
     %              +q1(ic,jc,kc)*ac1j(jc)
     %              +q1(ic,jm,kc)*am1j(jc)
c
c   33 second deriv. of q1(n)
c
             d33q1=((qb1n(ic,jc)-q1(ic,jc,kc))*codea*2. -
     %              (q1(ic,jc,kc)-q1(ic,jc,km))*codeb ) 
c
            dcq1 = d11q1 + d22q1 + d33q1
c
c   grad(pr) along 1
c
            dpx11=(pr(ic,jc,kc)-pr(im,jc,kc))*udx1
c
c
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c
c   evaluation of RHS
c
            rhs(iadd)=+(ga*dq(ic,jc,kc)+ro*ru1(ic,jc,kc)
     %                 +alre*dcq1*visc(kc)-dpx11)*dt
c
c
            ru1(ic,jc,kc)=dq(ic,jc,kc)
  102 continue
c
       if(infig.ne.-1) then
         usaldto = 1./aldto

!$OMP PARALLEL DO
!$OMP$ SHARED(npunt,indgeo,indgeoe,rhs,distb)
!$OMP$ SHARED(q1,q1bo,al,dt,usaldto,aldto)
!$OMP$ PRIVATE(i,j,k,ie,je,ke,q1e,iadd)
         do n=1,npunt
           i=indgeo(1,n,1)
           j=indgeo(1,n,2)
           k=indgeo(1,n,3)
           ie=indgeoe(1,n,1)
           je=indgeoe(1,n,2)
           ke=indgeoe(1,n,3)
           q1e=((al*dt+aldto)*q1(ie,je,ke)-al*dt*q1bo(n))*usaldto
c
c      punti sul contorno esterno
c
           iadd=i+(j-1)*n1m+(k-1)*n1m*n2m
           rhs(iadd) = -q1(i,j,k) + q1e*distb(1,n)
           q1bo(n)= q1(ie,je,ke)
         end do
!$OMP END PARALLEL DO

       end if
c
      timsol=etime(t)
      call solq1i(q1)
      timsol2=etime(t)
      call solq1j
      timsol3=etime(t)
      call solq1k(q1)
      timsol4=etime(t)
C     print *,'        in invtr1',-timinv+timsol4,timsol2-timsol,
C    c                 timsol3-timsol2,timsol4-timsol3

      return
      end
************************************************************************
c                       SUBROUTINE INVTR2
c   This subroutine performs the computation of Q~~ for the q2 momentum 
c   equation (radial direction) by a factored implicit scheme.
c   For details see the introduction of INVTR1
c   
c
      subroutine invtr2(q2,pr,ru2,q1,aldto)
      include'param.f'
      REAL   pr(m1,m2,m3)
      REAL q2(m1,m2,m3)
      REAL ru2(m1,m2,m3)
      REAL q1(m1,m2,m3)
      real etime, t(2)
      timinv=etime(t)
      alre=al/ren
c
c
c  compute the rhs of the factored equation
c  everything at i+1/2,j,k+1/2
c
c    points inside the flowfield
c
c
c     do kc=1,n3
c       do jc=1,n2
c         do ic=1,n1 
c           iadd=ic+(jc-1)*n1+(kc-1)*n1*n2
c           rhs(iadd) = 0.
c         end do
c       end do 
c     end do 
c
c
      n3mm=n3m-1
      ugkk=4./3.
        udx1q =  dx1q*iaxsy
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3mm,kmv,kpv,n2m,jpv,jmv,n1m,ipv,imv,n2)
!$OMP$ SHARED(q2,udx1q,ap2j,ac2j,am2j,ap2k,ac2k,am2k,pr,udx1)
!$OMP$ SHARED(ga,ro,dph,ru2,alre,visc,dt,udx2c,al)
!$OMP$ PRIVATE(km,kp,jp,jm,ip,im,d11q2,d22q2,d33q2,dcq2,dpx22,iadd)
!$OMP$ PRIVATE(aap,aac,aam)
        do 26 kc=2,n3mm
          km=kmv(kc)
          kp=kpv(kc)
          do 26 jc=2,n2m
            jm=jc-1
            jp=jc+1
            aap=ap2j(jc)
            aac=ac2j(jc)
            aam=am2j(jc)
             do 26 ic=1,n1m
               im=imv(ic)
               ip=ipv(ic)
c
c   11 second derivative of q2
c
            d11q2=(q2(ip,jc,kc)-2.*q2(ic,jc,kc)+q2(im,jc,kc))*
     %            udx1q
c
c   22 second derivative of q2
c
            d22q2=q2(ic,jp,kc)*aap+q2(ic,jc,kc)*aac+q2(ic,jm,kc)*aam
c
c   33 second derivative of q2
c
            d33q2=q2(ic,jc,kp)*ap2k(kc)
     %           +q2(ic,jc,kc)*ac2k(kc)
     %           +q2(ic,jc,km)*am2k(kc)
c
            dcq2 = d11q2 + d22q2 + d33q2
c
c   component of grad(pr) along 2 direction
c
            dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*udx2c(jc)*al
c
c
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
c
            rhs(iadd)=+(ga*dph(ic,jc,kc)+ro*ru2(ic,jc,kc)
     %                    +alre*dcq2*visc(kc)-dpx22)*dt
c
            ru2(ic,jc,kc)=dph(ic,jc,kc)
   26   continue
!$OMP  END PARALLEL DO
c

c       LOWER WALL
c
          kc=1
          kp=kc+1
          codea = dx3q/(g3c(kp)*g3m(kc))
          codeb = dx3q/(g3c(kc)*g3m(kc))
          do 101 jc=2,n2m
             jm=jc-1
             jp=jc+1
            aap=ap2j(jc)
            aac=ac2j(jc)
            aam=am2j(jc)
             do 101 ic=1,n1m
                im=imv(ic)
                ip=ipv(ic)
c
c   11 second derivative of q2
c
            d11q2=(q2(ip,jc,kc)-2.*q2(ic,jc,kc)+q2(im,jc,kc))*
     %            udx1q
c
c   22 second derivative of q2
c
            d22q2=q2(ic,jp,kc)*aap+q2(ic,jc,kc)*aac+q2(ic,jm,kc)*aam
c
c   33 second derivative of q2
c
             d33q2=((q2(ic,jc,kp)-q2(ic,jc,kc))*codea - 
     %              (q2(ic,jc,kc)-qb2s(ic,jc) )*codeb*2. )
c
            dcq2 = d11q2 + d22q2 + d33q2
c
c   grad(pr) along 2
c
              dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*udx2c(jc)*al
c
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
c
            rhs(iadd)=+(ga*dph(ic,jc,kc)+ro*ru2(ic,jc,kc)
     %                    +alre*dcq2*visc(kc)-dpx22)*dt
c
              ru2(ic,jc,kc)=dph(ic,jc,kc)
  101   continue


c
c       UPPER WALL
c
          kc=n3m
          km=kc-1
          kp=kc+1
          codea = dx3q/(g3c(kp)*g3m(kc))
          codeb = dx3q/(g3c(kc)*g3m(kc))
            do 102 jc=2,n2m
              jm=jc-1
              jp=jc+1
            aap=ap2j(jc)
            aac=ac2j(jc)
            aam=am2j(jc)
              do 102 ic=1,n1m
                 im=imv(ic)
                 ip=ipv(ic)
c
c   11 second derivative of q2
c
            d11q2=(q2(ip,jc,kc)-2.*q2(ic,jc,kc)+q2(im,jc,kc))*
     %            udx1q
c
c   22 second derivative of q2
c
            d22q2=q2(ic,jp,kc)*aap+q2(ic,jc,kc)*aac+q2(ic,jm,kc)*aam
c
c   33 second derivative of q2
c
             d33q2=((qb2n(ic,jc)-q2(ic,jc,kc))*codea*2. -
     %              (q2(ic,jc,kc)-q2(ic,jc,km))*codeb ) 
c
            dcq2 = d11q2 + d22q2 + d33q2
c
c   grad(pr) along 2
c
              dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*udx2c(jc)*al
c
c
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
c
            rhs(iadd)=+(ga*dph(ic,jc,kc)+ro*ru2(ic,jc,kc)
     %                    +alre*dcq2*visc(kc)-dpx22)*dt
c
              ru2(ic,jc,kc)=dph(ic,jc,kc)
  102   continue
c
       if(infig.ne.-1) then
         usaldto = 1./aldto
!$OMP PARALLEL DO
!$OMP$ SHARED(npunr,indgeo,indgeoe,rhs,distb)
!$OMP$ SHARED(q2,q2bo,al,dt,usaldto,aldto)
!$OMP$ PRIVATE(i,j,k,ie,je,ke,q2e,iadd)

         do n=1,npunr
             i=indgeo(2,n,1)
             j=indgeo(2,n,2)
             k=indgeo(2,n,3)
             ie=indgeoe(2,n,1)
             je=indgeoe(2,n,2)
             ke=indgeoe(2,n,3)

C          Modifica nuova versione NS
             dist = distb(2,n)  !new
             x3b=(x3m(k)-x3m(ke)*dist)/(1.-dist) !new
             x2b=(x2m(j)-x2m(je)*dist)/(1.-dist) !new             
             velvert = u0y-wx*(x3b-x30) !new
             q2e=((al*dt+aldto)*q2(ie,je,ke)-al*dt*q2bo(n))*usaldto
c
c      punti sul contorno esterno
c
             iadd=i+(j-1)*n1m+(k-1)*n1m*n2
             rhs(iadd) = -q2(i,j,k) + q2e*distb(2,n)
     %            +velvert*(1.-distb(2,n)) !new
c             if (n.eq.2000)write(98,*)' x3b - velvert ',x3b,velvert
c             write(99,*)x3b,x2b,velvert
C          Modifica nuova versione NS

             q2bo(n)= q2(ie,je,ke)

         end do
!$OMP END PARALLEL DO
       end if
c
c
      timsol=etime(t)
      call solq2i(q1)
      timsol2=etime(t)
      call solq2j
      timsol3=etime(t)
      call solq2k(q2)
      timsol4=etime(t)
C     print *,'        in invtr2',-timinv+timsol4,timsol2-timsol,
C    c                 timsol3-timsol2,timsol4-timsol3

      return
      end
************************************************************************
c                       SUBROUTINE INVTR3
c   This subroutine performs the computation of Q~~ for the q3 momentum 
c   equation (axial direction) by a factored implicit scheme.
c   Viscous terms are treated implicitly, nonlinear terms explicitly.
c   For details see the introduction of INVTR1 
c
      subroutine invtr3(q3,pr,ru3,q1,aldto)
      include'param.f'
      REAL   pr(m1,m2,m3)
      REAL ru3(m1,m2,m3)
      REAL q3(m1,m2,m3)
      REAL q1(m1,m2,m3)
      real etime,t(2)
      timinv=etime(t)
      alre=al/ren
c
c  compute the rhs of the factored equation
c  everything at i+1/2,j+1/2,k
c     
c
C     write(6,*) ' I1'
c     do kc=1,n3
c       do jc=1,n2
c         do ic=1,n1 
c           iadd=ic+(jc-1)*n1+(kc-1)*n1*n2
c           rhs(iadd) = 0.
c         end do
c       end do 
c     end do 
c
c
      udx1q = dx1q
      n2mm=n2m-1
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3mm,kmv,kpv,n2m,jpv,jmv,n1m,ipv,imv,n2)
!$OMP$ SHARED(q3,udx1q,ap3j,ac3j,am3j,ap3k,ac3k,am3k,pr)
!$OMP$ SHARED(ga,ro,qcap,ru3,alre,visc,dt,udx3c,al)
!$OMP$ PRIVATE(km,kp,jp,jm,ip,im,dq31,dq32,dq33,dcq3,dpx33,iadd)
!$OMP$ PRIVATE(aap,aac,aam,bap,bam,bac)

      do 18 kc=2,n3m
        km=kc-1
        kp=kc+1
          bap=ap3k(kc)
          bam=am3k(kc)
          bac=ac3k(kc)
c       do 18 jc=1,n2m
        do 18 jc=2,n2mm
c         jm=jmv(jc)
c         jp=jpv(jc)
          jm=jc-1
          jp=jc+1
          aap=ap3j(jc)
          aam=am3j(jc)
          aac=ac3j(jc)
          do 18  ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
c
c   11 second derivatives of q3
c
            dq31=(q3(ip,jc,kc)-2.*q3(ic,jc,kc)+q3(im,jc,kc))*udx1q
c      if (jc.eq.40.and.kc.eq.60.and.ic.eq.n1m)then
c         write(*,*)'invtr3 - 1',q3(ip,jc,kc),q3(ic,jc,kc),q3(im,jc,kc)
c     $        ,q3(im-1,jc,kc),q3(im-2,jc,kc)
c      endif
c
c   22 second derivatives of q3
c
            dq32=q3(ic,jm,kc)*aam
     %          +q3(ic,jc,kc)*aac
     %          +q3(ic,jp,kc)*aap
c
c
c   33 second derivatives of q3
c
            dq33=q3(ic,jc,kp)*bap
     %          +q3(ic,jc,kc)*bac
     %          +q3(ic,jc,km)*bam
c
c   viscous terms
c
            dcq3=dq31+dq32+dq33
c
c  component of grad(pr) along x3 direction
c
            dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*udx3c(kc)*al
c
c
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c
c  right hand side of the momentum equation
c
            rhs(iadd)=+(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc)
     %                    +alre*dcq3*visc(kc)-dpx33)*dt
c
c  updating of the non-linear terms
c
            ru3(ic,jc,kc)=qcap(ic,jc,kc)
   18 continue
!$OMP  END PARALLEL DO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      jc=1
      jp=jc+1
      codea = dx2q/(g2c(jp)*g2m(jc))
      codeb = dx2q/(g2c(jc)*g2m(jc))
      do 19 kc=2,n3m
        km=kc-1
        kp=kc+1
          bap=ap3k(kc)
          bam=am3k(kc)
          bac=ac3k(kc)
c       do 18 jc=1,n2m
c        do 18 jc=2,n2mm
c         jm=jmv(jc)
c         jp=jpv(jc)
c         jm=jc-1
c         jp=jc+1
c         aap=ap3j(jc)
c         aam=am3j(jc)
c         aac=ac3j(jc)
          do 19  ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
c
c   11 second derivatives of q3
c
            dq31=(q3(ip,jc,kc)-2.*q3(ic,jc,kc)+q3(im,jc,kc))*udx1q
c
c   22 second derivatives of q3
c
c           dq32=q3(ic,jm,kc)*aam
c    %          +q3(ic,jc,kc)*aac
c    %          +q3(ic,jp,kc)*aap
            dq32=((q3(ic,jp,kc)-q3(ic,jc,kc))*codea - 
     %           (q3(ic,jc,kc)-qb3dn(ic,kc))*codeb*2. )

c
c
c   33 second derivatives of q3
c
            dq33=q3(ic,jc,kp)*bap
     %          +q3(ic,jc,kc)*bac
     %          +q3(ic,jc,km)*bam
c
c   viscous terms
c
            dcq3=dq31+dq32+dq33
c
c  component of grad(pr) along x3 direction
c
            dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*udx3c(kc)*al
c
c
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c
c  right hand side of the momentum equation
c
            rhs(iadd)=+(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc)
     %                    +alre*dcq3*visc(kc)-dpx33)*dt
c
c  updating of the non-linear terms
c
            ru3(ic,jc,kc)=qcap(ic,jc,kc)
   19 continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      jc=n2m
      jm=jc-1
      jp=jc+1
      codea = dx2q/(g2c(jp)*g2m(jc))
      codeb = dx2q/(g2c(jc)*g2m(jc))
      do 20 kc=2,n3m
        km=kc-1
        kp=kc+1
          bap=ap3k(kc)
          bam=am3k(kc)
          bac=ac3k(kc)
c       do 18 jc=1,n2m
c        do 18 jc=2,n2mm
c         jm=jmv(jc)
c         jp=jpv(jc)
c         jm=jc-1
c         jp=jc+1
c         aap=ap3j(jc)
c         aam=am3j(jc)
c         aac=ac3j(jc)
          do 20  ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
c
c   11 second derivatives of q3
c
            dq31=(q3(ip,jc,kc)-2.*q3(ic,jc,kc)+q3(im,jc,kc))*udx1q
c
c   22 second derivatives of q3
c
c           dq32=q3(ic,jm,kc)*aam
c    %          +q3(ic,jc,kc)*aac
c    %          +q3(ic,jp,kc)*aap
            dq32=((qb3up(ic,kc)-q3(ic,jc,kc))*codea*2. -
     %              (q3(ic,jc,kc)-q3(ic,jm,kc))*codeb ) 

c
c
c   33 second derivatives of q3
c
            dq33=q3(ic,jc,kp)*bap
     %          +q3(ic,jc,kc)*bac
     %          +q3(ic,jc,km)*bam
c
c   viscous terms
c
            dcq3=dq31+dq32+dq33
c
c  component of grad(pr) along x3 direction
c
            dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*udx3c(kc)*al
c
c
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c
c  right hand side of the momentum equation
c
            rhs(iadd)=+(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc)
     %                    +alre*dcq3*visc(kc)-dpx33)*dt
c
c  updating of the non-linear terms
c
            ru3(ic,jc,kc)=qcap(ic,jc,kc)
   20 continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       if(infig.ne.-1) then
         usaldto = 1./aldto
!$OMP PARALLEL DO
!$OMP$ SHARED(npunz,indgeo,indgeoe,rhs,distb)
!$OMP$ SHARED(q3,q3bo,al,dt,usaldto,aldto)
!$OMP$ PRIVATE(i,j,k,ie,je,ke,q3e,iadd)
         do n=1,npunz
           i=indgeo(3,n,1)
           j=indgeo(3,n,2)
           k=indgeo(3,n,3)
           ie=indgeoe(3,n,1)
           je=indgeoe(3,n,2)
           ke=indgeoe(3,n,3)

c          Modifica nuova versione NS
           dist = distb(3,n)   !new
           x3b=(x3m(k)-x3m(ke)*dist)/(1.-dist) !new
           x2b=(x2m(j)-x2m(je)*dist)/(1.-dist)   !new
           velhor = u0z + wx*(x2b-x20) !new
           q3e=((al*dt+aldto)*q3(ie,je,ke)-al*dt*q3bo(n))*usaldto
c
c      punti sul contorno esterno
c
           iadd=i+(j-1)*n1m+(k-1)*n1m*n2m
           rhs(iadd) = -q3(i,j,k) + q3e*distb(3,n)
     $          +velhor*(1.-distb(3,n))   !new
c           if (n.eq.2000)write(98,*)' x2b - velhor ',x2b,velhor
c           write(99,*)x3b,x2b,velhor
c          Modifica nuova versione NS

           q3bo(n)= q3(ie,je,ke)
C          rhsm = max(rhsm,abs(rhs(iaddi)))
         end do
!$OMP END PARALLEL DO
       end if
c
      timsol=etime(t)
      call solq3i(q1)

      timsol2=etime(t)
      call solq3j

      timsol3=etime(t)
      call solq3k(q3)
      timsol4=etime(t)
C     print *,'        in invtr3',-timinv+timsol4,timsol2-timsol,
C    c                 timsol3-timsol2,timsol4-timsol3

      return
      end
c
c***********************************************************************
c                                                                      * 
c    subroutine generating random numbers                              *
c                                                                      *
c***********************************************************************
      subroutine gerand(idum,npa,r)
      include 'param.f'
      dimension r(10000),x(10000)
      data iff /0/
      mk1=259200
      ia1=7141
      ic1=54773
      rm1=3.8580247e-6
      mk2=134456
      ia2=8121
      ic2=28411
      rm2=7.4373773e-6
      mk3=243000
      ia3=4561
      ic3=51349
      if (idum.lt.0.or.iff.eq.0) then
        iff=1
        ix1=mod(ic1-idum,mk1)
        ix1=mod(ia1*ix1+ic1,mk1)
        ix2=mod(ix1,mk2)
        ix1=mod(ia1*ix1+ic1,mk1)
        ix3=mod(ix1,mk3)
        do 11 j=1,npa
          ix1=mod(ia1*ix1+ic1,mk1)
          ix2=mod(ia2*ix2+ic2,mk2)
          r(j)=(float(ix1)+float(ix2)*rm2)*rm1
11      continue
        idum=1
      endif
      ix1=mod(ia1*ix1+ic1,mk1)
      ix2=mod(ia2*ix2+ic2,mk2)
      ix3=mod(ia3*ix3+ic3,mk3)
      j=1+(npa*ix3)/mk3
      if(j.gt.npa.or.j.lt.1) pause
      ran1=r(j)
      ll=0
      do 12 l=j,npa
      ll=ll+1
      x(ll)=r(l)
   12 continue
      do 13 l=1,j-1
      ll=ll+1
      x(ll)=r(l)
   13 continue
      r(j)=(float(ix1)+float(ix2)*rm2)*rm1
      return
      end
c************************************************************************
c                                                                       *
c ****************************** subrout coetar  ********************** *
c                                                                       *
c    this subroutine calculates the coefficients for the                *
c    integration in the radial direction with non-uniform coor. trasf.  *
c                                                                       *
c************************************************************************
      subroutine coetar
      include 'param.f'
c
c  ***********  coefficients for q1   inner points
c
c
      do 151 jc=2,n2m-1
      jp=jc+1
      jm=jc-1
      a22=dx2q/g2m(jc)
      a22p= +a22/g2c(jp)
      a22m= +a22/g2c(jc)
      ap1j(jc)=a22p
      am1j(jc)=a22m
      ac1j(jc)=-(a22p+a22m)
  151 continue
c
c    external bound. conditions
c
      jc=n2m
      jp=jc+1
      jm=jc-1
      a22=dx2q/g2m(jc)
      apjjc=a22/g2c(jp)*2.
      acjjc=a22/g2c(jc)
      amjjc=a22/g2c(jc)
       if(inslww.eq.1) then
        ap1j(jc)=-apjjc
       else
        ap1j(jc)=0.
       end if
      ac1j(jc)=-acjjc
      am1j(jc)=amjjc
c
c    internal boundary conditions
c
      jc=1
      jp=jc+1
      a22=dx2q/g2m(jc)
      a22p=a22/g2c(jp)
      a22m=a22/g2c(jc)*2.
      ap1j(jc)=a22p
      ac1j(jc)=-a22p
       if(inslwe.eq.1) then
      am1j(jc)=-a22m
       else
      am1j(jc)=0.
       end if
c
c  ***********  coefficients for q2   
c
      am2j(1)=0.
      ap2j(1)=0.
      ac2j(1)=1.
      am2j(n2)=0.
      ap2j(n2)=0.
      ac2j(n2)=1.
      do 2 jc=2,n2m
      jm=jc-1
      jp=jc+1
      a22=dx2q/g2c(jc)
      a22p=1./g2m(jc)
      a22m=1./g2m(jm)
      ap2j(jc)=a22*a22p
      am2j(jc)=a22*a22m
      ac2j(jc)=-(a22*a22p+a22*a22m)
    2 continue
c
c
c
c  ***********  coefficients for q3   inner points
c  *********** are equal to those for dens
c
      do 51 jc=2,n2m-1
      jp=jc+1
      jm=jc-1
      a22=dx2q/g2m(jc)
      a22p= +a22/g2c(jp)
      a22m= +a22/g2c(jc)
      ap3j(jc)=a22p
      am3j(jc)=a22m
      ac3j(jc)=-(a22p+a22m)
      apscj(jc)=a22p
      amscj(jc)=a22m
      acscj(jc)=-(a22p+a22m)
   51 continue
c     
c    external bound. conditions  q3
c     
      jc=n2m
      jp=jc+1
      jm=jc-1
      a22=dx2q/g2m(jc)
      apjjc=a22/g2c(jp)*2.
      acjjc=a22/g2c(jc)
      amjjc=a22/g2c(jc)
       if(inslww.eq.1) then
        ap3j(jc)=-apjjc
       else
        ap3j(jc)=0.
       end if
      ac3j(jc)=-acjjc
      am3j(jc)=amjjc
c
c    internal boundary conditions
c
      jc=1
      jp=jc+1
      a22=dx2q/g2m(jc)
      a22p=a22/g2c(jp)
      a22m=a22/g2c(jc)*2.
      ap3j(jc)=a22p
      ac3j(jc)=-a22p
       if(inslwe.eq.1) then
      am3j(jc)=-a22m
       else
      am3j(jc)=0.
       end if
c
c
c  ***********  coefficients for q3   for x3 differentation
c  c means centered that is at k location
c
      am3k(1)=0.
      ap3k(1)=0.
      ac3k(1)=1.
      am3k(n3)=0.
      ap3k(n3)=0.
      ac3k(n3)=1.
      do kc=2,n3m
      km=kc-1
      kp=kc+1
      a33=dx3q/g3c(kc)
      a33p=1./g3m(kc)
      a33m=1./g3m(km)
      ap3k(kc)=a33*a33p
      am3k(kc)=a33*a33m
      ac3k(kc)=-(ap3k(kc)+am3k(kc))
      enddo
c
c
c  **coefficients for q1, q2 ap1/2k,am1/2k,ac1/2k
c   ap3ssk,am3ssk,ac3ssk, psc   for x3 differentation
c  s means staggered that is at k+1/2 location
c
      do kc=2,n3m-1
      kp=kc+1
      km=kc-1
      a33=dx3q/g3m(kc)
      a33p= +a33/g3c(kp)
      a33m= +a33/g3c(kc)
      ap2k(kc)=a33p
      am2k(kc)=a33m
      ac2k(kc)=-(ap2k(kc)+am2k(kc))
      ap1k(kc)=a33p
      am1k(kc)=a33m
      ac1k(kc)=-(ap1k(kc)+am1k(kc))
      ap3ssk(kc)=ap2k(kc)
      am3ssk(kc)=am2k(kc)
      ac3ssk(kc)=ac2k(kc)
      enddo
c    
c    lower wall  bound. conditions  indicated by inslws
c    differemtiation of sec. derivative at 1+1/2
c    
      kc=1
      kp=kc+1
      a33=dx3q/g3m(kc)
      a33p= +a33/g3c(kp)
      a33m= +a33/g3c(kc)
      ap2k(kc)=a33p
      am2k(kc)=0.
      ac2k(kc)=-(a33p+inslws*a33m*2.)
      ap1k(kc)=a33p
      am1k(kc)=0.
      ac1k(kc)=-(a33p+inslws*a33m*2.)
      ap3ssk(kc)=ap2k(kc)
      am3ssk(kc)=am2k(kc)
      ac3ssk(kc)=-(a33p)
c    
c    upper wall  bound. conditions  indicated by inslws
c    differentiation of sec. derivative at n3-1/2
c    
      kc=n3m
      kp=kc+1
      a33=dx3q/g3m(kc)
      a33p= +a33/g3c(kp)
      a33m= +a33/g3c(kc)
      am2k(kc)=a33m
      ap2k(kc)=0.
      ac2k(kc)=-(a33m+inslwn*a33p*2.)
      am1k(kc)=a33m
      ap1k(kc)=0.
      ac1k(kc)=-(a33m+inslwn*a33p*2.)
      ap3ssk(kc)=ap2k(kc)
      am3ssk(kc)=am2k(kc)
      ac3ssk(kc)=-(a33m)
C
C     write(6,*) ' centred coeff. K '
C     do k=1,n3
C     write(6,*) k, am3k(k), ac3k(k), ap3k(k)
C     end do
C     pause
C     write(6,*) ' staggered coeff. K '
C     do k=1,n3m
C     write(6,*) k, am2k(k), ac2k(k), ap2k(k)
C     end do
C     pause
C     write(6,*) ' staggered coeff. K (2)'
C     do k=1,n3m
C     write(6,*) k, am3ssk(k), ac3ssk(k), ap3ssk(k)
C     end do
C     pause
c
c     additional coefficients
c
      do j=1,n2m
        udx2m(j) = dx2/g2m(j)
        udx2c(j) = dx2/g2c(j)
      end do
      udx2c(n2) = dx2/g2c(n2)
c
      do k=1,n3m
        udx3m(k) = dx3/g3m(k)
        udx3c(k) = dx3/g3c(k)
      end do
      udx3c(n3) = dx3/g3c(n3)
      return
      end

c
c***********************************************************************
c                                                                      *
c                       CONDIZIONI INIZIALI                            *
c                                                                      *
c***********************************************************************
      subroutine inqpr(q1,q2,q3,pr)                         
      include'param.f'
      REAL q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)
      REAL   pr(m1,m2,m3)
      pi=2.*asin(1.)                 
      do 810 k=1,n3                 
         do 810 j=1,n2             
            do 810 i=1,n1         
              q3(i,j,k)=0.
              q2(i,j,k)=0.
c             q3(i,j,k)=0.45473*cos(15./180.*pi)+0.5*sin(15./180.*pi)
c             q2(i,j,k)=0.45473*sin(15./180.*pi)-0.5*cos(15./180.*pi)
              q1(i,j,k)=0.
              pr(i,j,k)=0.
  810 continue                
c            write(*,*) 'q3ini= ',q3(4,40,10),'q2ini= ',q2(4,40,10)
      call ini
c       write(6,*) ' dopo ini '
      return                                                            
      end                                                               
c                                                                       
c  **************  function ftir                                        
c                                                                       
      function  ftir(tina)
      if(tina.ge.1.) then
      ftir=1.
                     else
      ftir=3.*tina**2-2.*tina**3
                     endif
      return
      end
      
      function dftir(tina)
      if(tina.ge.1.) then
      dftir=0.
                     else
      dftir=6.*tina-6.*tina**2
                     endif
      return
      end
      
      function ddftir(tina)
      if(tina.ge.1.) then
      ddftir=0.
                     else
      ddftir=6.-12.*tina**2
                     endif
      return
      end
c                                                                       
c
c  time dependent inflow axial velocity profile for the
c     RING FORMATION
c
      subroutine ringfo(tin,tino)
      include 'param.f'
      pi=2.*asin(1.)
      if(tin.le.tau2) then
         tina=tin/tau1
         tinao=tino/tau1
         ft=ftir(tina)
         fto=ftir(tinao)
         dft=ft-fto
         dftdt=dftir(tina)
         ddftdt=ddftir(tina)
      else
         if(tin.lt.tau3) then
            tina=(tin-tau2)/tau1
            tinao=(tino-tau2)/tau1
            ft=1.-ftir(tina)
            fto=1.-ftir(tinao)
            dft=(ft-fto)
            dftdt=dftir(tina)
            ddftdt=ddftir(tina)
         else
            ft=0. 
            dft=0. 
            dftdt=0.
            ddftdt=0.
         endif
      endif


c     write(50,*) tin,ft,dft
c     write(6,*) ' RINGFO ',tin,tino,ft,dft,dftdt,ddftdt
      return
      end
c***********************************************************************
c                                                                      *
c  **********************  inizialize inflow profiles      *********   *
c                                                                      *
c***********************************************************************
      subroutine ini
      include 'param.f'
      dimension uve(m2)
      dimension va(m2)
c                                                                       
      pi=2.*asin(1.)                                                    
c                                                                        
      rad = x2c(n2)
c     rad = radinf
      qint = 0.
      dr = 1./dx2
      do j=1,n2m
         xi=x2m(j)
         if(xi.gt.rad) go to 501
         y=rad-xi
         qint = qint +
     %   tanh(all*y)/tanh(all*(rad))*dr*g2m(j)
      end do
 501  continue
      qintv = 0.
      do j=1,n2m
         xi=x2m(j)
         if(xi.gt.rad) go to 502
         y=rad-xi
         qintv = qintv +
     %   tanh(all*y)/tanh(all*(rad))*dr*g2m(j)
      end do
 502  continue
      fasca = qint/qintv
c     write(6,*) ' SCALE FACTOR FOR VEL. PROF. ',fasca
c     fasca = 0.5*rpis**2*velvalp/qintv
c     fasca = 2.* velval
c     fasca= 1.
c
c     uinflw=1./(sigma*amp)
      rad = alx2
      jv3 = -10
      do i=1,n1m                                                    
        do  j=1,n2m
          xi=x2m(j)
          if(xi.le.rad) then
            y=rad-xi
            jv3 = j
c           uinfth(i,j)=0.45473 
            uinfth(i,j)= uzero
c           uinfth(i,j)= 0.
c    %                  (1.-x2m(j)**2/rext**2)*fasca
c    %      tanh(all*y)/tanh(all*(rad))*fasca
c    %      tanh(all*y)/tanh(all*(rad))
          else
            uinfth(i,j)= uzero
c            uinfth(i,j)=0.45473
          endif
          uve(j)=uinfth(i,j)
c        if(i.eq.1) write(21,*) j,xi,y,uve(j)
        enddo
        vamax=0.
        amo=0.
        if(i.eq.1) then
        do  j=1,n2m
         amo=amo+(1.-uve(j))*uve(j)/dx2*g2m(j)
         write(21,*) j,uve(j),g2m(j)
        enddo
c       write(6,*) fasca,velvalp,uve(1)
        amo=1./amo
        do  j=2,n2m
          xi=x2c(j)  
          va(j)=-(uve(j)-uve(j-1))*dx2/g2c(j)
          vamax=max(abs(va(j)),vamax)
c         if(i.eq.1) write(13,1000) j,xi,uve(j),va(j)
        enddo
c       write(6,*) ' JV3 = ',jv3
        va(1)=va(2)
        va(n2)=va(n2m)
      end if
      enddo
 1000 format(1x,i3,4(2x,e12.6))
       close(21)
      return                                                            
      end                                                               
c
c************************************************************************
c
      subroutine brtrk(aa,bb,cc,r,m1,m2,m3,forclo,mm1,mm2,mm3)
c     a() sottodiagonale
c     b() diagonale
c     c() sopradiagonale
c     r() vettore termine noti e soluzione
c     m1,m2,m3   dimensioni della matrice
      dimension aa(1),bb(1),cc(1)
      REAL r(1)
      parameter (mg=1500)
      dimension ggm(mg),a(mg),b(mg),c(mg)
      dimension qqq(mg)
      REAL forclo(mm1,mm2,mm3)
c***
      if(m2.gt.mg) go to 8888
c***
      do j=1,m2
        do i=1,m1
          do k=1,m3
            a(k) = aa(k)*(1.-forclo(i,j,k))
            b(k) = bb(k)*(1.-forclo(i,j,k))+1.*forclo(i,j,k)
            c(k) = cc(k)*(1.-forclo(i,j,k))
          end do
          bet=b(1)
          subet=1./bet
          iadd=i+(j-1)*m1
          qqq(1)=r(iadd)*subet
c***
          do k=2,m3
            ggm(k)=c(k-1)/bet
            bet=b(k)-a(k)*ggm(k)
            subet=1./bet
            iaddc=i+(j-1)*m1+(k-1)*m1*m2
            qqq(k)=(r(iaddc)-a(k)*qqq(k-1))*subet
          end do
          do k=m3-1,1,-1
            qqq(k)=qqq(k)-ggm(k+1)*qqq(k+1)
          end do
          do k=1,m3
            iaddc=i+(j-1)*m1+(k-1)*m1*m2
            r(iaddc)=qqq(k)
          end do
        end do
      end do
c***
      return
8888  print*,'mg too small'
      stop
      end
c
c************************************************************************
c
      subroutine brtrj(aa,bb,cc,r,m1,m2,m3,forclo,mm1,mm2,mm3,visc)
c     a() sottodiagonale
c     b() diagonale
c     c() sopradiagonale
c     r() vettore termine noti e soluzione
c     m1,m2,m3   dimensioni della matrice
      dimension aa(1),bb(1),cc(1)
      REAL r(1)
      parameter (mg=1500)
      dimension ggm(mg),a(mg),b(mg),c(mg)
      dimension qqq(mg)
      REAL forclo(mm1,mm2,mm3)
      dimension   visc(mm3)
c***
      if(m2.gt.mg) go to 8888
c***
      do k=1,m3
        do i=1,m1
          do j=1,m2
            a(j) = aa(j)*(1.-forclo(i,j,k))*visc(k)
            b(j) = ((bb(j)-1.)*visc(k)+1)*(1.-forclo(i,j,k))
     %               +1.*forclo(i,j,k)
c           b(j) = bb(j)*(1.-forclo(i,j,k))+1.*forclo(i,j,k)
            c(j) = cc(j)*(1.-forclo(i,j,k))*visc(k)
          end do
          bet=b(1)
          subet=1./bet
          iadd=i+(k-1)*m1*m2
          qqq(1)=r(iadd)*subet
c***
          do j=2,m2
            ggm(j)=c(j-1)/bet
            bet=b(j)-a(j)*ggm(j)
            subet=1./bet
            iaddc=i+(j-1)*m1+(k-1)*m1*m2
            qqq(j)=(r(iaddc)-a(j)*qqq(j-1))*subet
          end do
          do j=m2-1,1,-1
            qqq(j)=qqq(j)-ggm(j+1)*qqq(j+1)
          end do
          do j=1,m2
            iaddc=i+(j-1)*m1+(k-1)*m1*m2
            r(iaddc)=qqq(j)
          end do
        end do
      end do
c***
      return
8888  print*,'mg too small'
      stop
      end

************************************************************************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c       direction x3
c
      subroutine solq3k(q3)
      include'param.f'
      REAL q3(m1,m2,m3)
      dimension amkl(m1,m3),apkl(m1,m3),ackl(m1,m3),fkl(m1,m3)
c  ********* compute the dq3* sweeping in the x3 direction
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3,n2,n1,forclo)

      do k=1,n3
        do j=1,n2
          do i=1,n1
            forclo(i,j,k)=1.
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
       if(infig.ne.-1) then
!$OMP  PARALLEL DO
!$OMP$ SHARED(indgeo,npunz,forclo)
      do n=1,npunz
        i=indgeo(3,n,1)
        j=indgeo(3,n,2)
        k=indgeo(3,n,3)
        forclo(i,j,k)= 0.
      end do
!$OMP  END PARALLEL DO
      end if
c
      betadx=beta*al
c

!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,n3)
!$OMP$ SHARED(q3,forclo,dqb3s,dqb3n,am3k,betadx,visc,ac3k,ap3k,rhs)
!$OMP$ PRIVATE(amkl,apkl,ackl,fkl,fclo,iadd)

      do j=1,n2m
         do i=1,n1m
            amkl(i,1)=0.
            apkl(i,1)=0.
            ackl(i,1)=1.
            amkl(i,n3)=0.
            apkl(i,n3)=0.
            ackl(i,n3)=1.
            fkl(i,1)=dqb3s(i,j)
            fkl(i,n3)=dqb3n(i,j)
         end do
         do k=2,n3m
            do i=1,n1m
               fclo=forclo(i,j,k)*betadx*visc(k)
               amkl(i,k)=-am3k(k)*fclo
               ackl(i,k)=1.-ac3k(k)*fclo
               apkl(i,k)=-ap3k(k)*fclo
               iadd=i+(j-1)*n1m+(k-1)*n1m*n2m
               fkl(i,k)=rhs(iadd)
            end do
         end do
c     
         call trvb(amkl,ackl,apkl,fkl,m1,m3,n3,n1m)
c     
         do k=2,n3m
            do i=1,n1m
               iadd=i+(j-1)*n1m+(k-1)*n1m*n2m
               q3(i,j,k)=q3(i,j,k) + fkl(i,k)
               rhs(iadd) = 0.
            end do
         end do
      end do
!$OMP  END PARALLEL DO

c
      return
      end

c
c  ****************************** subrout trvb  **********************
c
      subroutine trvb(am,ac,ap,f,mm,nn,n,m)
c
      dimension am(mm,nn),ac(mm,nn),ap(mm,nn),f(mm,nn)
c
c  ******** reduction of trid. matrix to an upper rigth matrix
c
      do 1 i=2,n
      do 2 k=1,m
      ac(k,i)=ac(k,i)*ac(k,i-1)-ap(k,i-1)*am(k,i)
      ap(k,i)=ap(k,i)*ac(k,i-1)
      f(k,i)=f(k,i)*ac(k,i-1)-f(k,i-1)*am(k,i)
    2 continue
    1 continue
c  ******** calculation of the unknown by backward elimination
c
      do 3 k=1,m
      f(k,n)=f(k,n)/ac(k,n)
    3 continue
      nm=n-1
      do 10 ii=1,nm
      i=n-ii
      do 11 k=1,m
      f(k,i)=(f(k,i)-ap(k,i)*f(k,i+1))/ac(k,i)
   11 continue
   10 continue
      return
      end

************************************************************************
c   this subroutine performs the inversion of the q2 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q2
c   are treated implicitly
c   in the first part the rhs is calculated
c        direction x3
c
      subroutine solq2k(q2)
      include'param.f'
      REAL q2(m1,m2,m3)
      dimension amkl(m1,m3),apkl(m1,m3),ackl(m1,m3),fkl(m1,m3)

c
c  ************ compute dq2 sweeping along the x3 direction
c               periodic
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3,n2,n1,forclo)
      do k=1,n3
        do j=1,n2
          do i=1,n1
            forclo(i,j,k)=1.
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
       if(infig.ne.-1) then
!$OMP  PARALLEL DO
!$OMP$ SHARED(indgeo,npunr,forclo)

      do n=1,npunr
        i=indgeo(2,n,1)
        j=indgeo(2,n,2)
        k=indgeo(2,n,3)
        forclo(i,j,k)= 0.
      end do
!$OMP  END PARALLEL DO
      end if
c
      betadx=beta*al
      ugkks = 1./(g3m(1)*g3c(1))*dx3q*2.*betadx
      ugkkn = 1./(g3m(n3m)*g3c(n3))*dx3q*2.*betadx

!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,n2,ugkks,ugkkn)
!$OMP$ SHARED(q2,forclo,dqb2s,dqb2n,am1k,betadx,visc,ac1k,ap1k,rhs)
!$OMP$ PRIVATE(amkl,apkl,ackl,fkl,fclo,iadd,fclo1,fclon3m)

c      do j=1,n2m
       do j=1,n2  !con questo funziona ma non ha senso
c      do j=2,n2m

          do k=2,n3m-1
            do i=1,n1m
            fclo=forclo(i,j,k)*betadx*visc(k)
            amkl(i,k)=-am2k(k)*fclo
            ackl(i,k)=1.-ac2k(k)*fclo
            apkl(i,k)=-ap2k(k)*fclo
            iadd=i+(j-1)*n1m+(k-1)*n1m*n2
            fkl(i,k)=rhs(iadd)
          end do
        end do

c k=1 and k=n3m

        do i=1,n1m
         iadd=i+(j-1)*n1m
            fclo1=forclo(i,j,1)*betadx*visc(1)
         fkl(i,1)= rhs(iadd) + dqb2s(i,j)*ugkks
            amkl(i,1)=-am2k(1)*fclo1
            ackl(i,1)=1.-ac2k(1)*fclo1
            apkl(i,1)=-ap2k(1)*fclo1
         iadd=i+(j-1)*n1m+(n3m-1)*n1m*n2
            fclon3m=forclo(i,j,n3m)*betadx*visc(n3m)
         fkl(i,n3m)= rhs(iadd) + dqb2n(i,j)*ugkkn
            amkl(i,n3m)=-am2k(n3m)*fclon3m
            ackl(i,n3m)=1.-ac2k(n3m)*fclon3m
            apkl(i,n3m)=-ap2k(n3m)*fclon3m
        end do
c
      call trvb(amkl,ackl,apkl,fkl,m1,m3,n3m,n1m)
c
      do k=1,n3m
        do i=1,n1m
         iadd=i+(j-1)*n1m+(k-1)*n1m*n2
         q2(i,j,k)=q2(i,j,k) + fkl(i,k)
         rhs(iadd) = 0.
        end do
       end do
      end do
!$OMP  END PARALLEL DO
c
      return
      end

************************************************************************
c
c   this subroutine performs the inversion of the q1 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q1
c   are treated implicitly
c       direction x3
c
      subroutine solq1k(q1)
      include'param.f'
      REAL q1(m1,m2,m3)
      dimension amkl(m1,m3),apkl(m1,m3),ackl(m1,m3),fkl(m1,m3)

c
c   compute  from dq** sweeping along the x3 periodic direction
c   
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3,n2,n1,forclo)
      do k=1,n3
        do j=1,n2
          do i=1,n1
            forclo(i,j,k)=1.
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
       if(infig.ne.-1) then
!$OMP  PARALLEL DO
!$OMP$ SHARED(indgeo,npunt,forclo)
      do n=1,npunt
          i=indgeo(1,n,1)
          j=indgeo(1,n,2)
          k=indgeo(1,n,3)
          forclo(i,j,k)= 0.
      end do
!$OMP  END PARALLEL DO
      end if
c
      betadx=beta*al
       ugkks = 1./(g3m(1)*g3c(1))*dx3q*2.*betadx
       ugkkn = 1./(g3m(n3m)*g3c(n3))*dx3q*2.*betadx
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,ugkks,ugkkn)
!$OMP$ SHARED(q1,forclo,dqb1s,dqb1n,am1k,betadx,visc,ac1k,ap1k,rhs)
!$OMP$ PRIVATE(amkl,apkl,ackl,fkl,fclo,iadd,fclo1,fclon3m)
      do j=1,n2m

          do k=2,n3m-1
            do i=1,n1m
            fclo=forclo(i,j,k)*betadx*visc(k)
            amkl(i,k)=-am1k(k)*fclo
            ackl(i,k)=1.-ac1k(k)*fclo
            apkl(i,k)=-ap1k(k)*fclo
            iadd=i+(j-1)*n1m+(k-1)*n1m*n2m
            fkl(i,k)=rhs(iadd)
          end do
        end do

c k=1 and k=n3m

        do i=1,n1m
         iadd=i+(j-1)*n1m
            fclo1=forclo(i,j,1)*betadx*visc(1)
         fkl(i,1)= rhs(iadd) + dqb1s(i,j)*ugkks
            amkl(i,1)=-am1k(1)*fclo1
            ackl(i,1)=1.-ac1k(1)*fclo1
            apkl(i,1)=-ap1k(1)*fclo1
         iadd=i+(j-1)*n1m+(n3m-1)*n1m*n2m
            fclon3m=forclo(i,j,n3m)*betadx*visc(n3m)
         fkl(i,n3m)= rhs(iadd) + dqb1n(i,j)*ugkkn
            amkl(i,n3m)=-am1k(n3m)*fclon3m
            ackl(i,n3m)=1.-ac1k(n3m)*fclon3m
            apkl(i,n3m)=-ap1k(n3m)*fclon3m
        end do

      call trvb(amkl,ackl,apkl,fkl,m1,m3,n3m,n1m)
c
      do kc=1,n3m
        do ic=1,n1m
         iadd=ic+(j-1)*n1m+(kc-1)*n1m*n2m
         q1(ic,j,kc)=q1(ic,j,kc) + fkl(ic,kc)
         rhs(iadd) = 0.
        end do
       end do
      end do
!$OMP  END PARALLEL DO

c
c
      return
      end
c ************************************************************************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c       direction x2
c
      subroutine solq3j
      include'param.f'
      dimension amjl(m1,m2),apjl(m1,m2),acjl(m1,m2),fjl(m1,m2)

c
c  ************ compute dq3 sweeping along the x2 direction
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3,n2,n1,forclo)
      do k=1,n3
        do j=1,n2
          do i=1,n1
            forclo(i,j,k)=1.
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
       if(infig.ne.-1) then
!$OMP  PARALLEL DO
!$OMP$ SHARED(indgeo,npunz,forclo)
      do n=1,npunz
        i=indgeo(3,n,1)
        j=indgeo(3,n,2)
        k=indgeo(3,n,3)
        forclo(i,j,k)= 0.
      end do
!$OMP  END PARALLEL DO
      end if
c
      betadx=beta*al
      ugkkdn=1./(g2m(1)*g2c(1))*dx2q*2.*betadx
      ugkkup=1./(g2m(n2m)*g2c(n2m))*dx2q*2.*betadx
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,betadx,forclo,rhs)
!$OMP$ PRIVATE(amjl,apjl,acjl,fjl,fclo,iadd)
      do k=2,n3m
        do j=2,n2m-1
          do i=1,n1m
            fclo=forclo(i,j,k)*betadx
            amjl(i,j)=-am3j(j)*fclo
            acjl(i,j)=1.-ac3j(j)*fclo
            apjl(i,j)=-ap3j(j)*fclo
            iadd=i+(j-1)*n1m+(k-1)*n1m*n2m
            fjl(i,j)=rhs(iadd)
          end do
        end do

          do i=1,n1m
            fclo1=forclo(i,1,k)*betadx
            amjl(i,1)=-am3j(1)*fclo1
            acjl(i,1)=1.-ac3j(1)*fclo1
            apjl(i,1)=-ap3j(1)*fclo1
            iadd=i+(k-1)*n1m*n2m
            fjl(i,1)=rhs(iadd) + dqb3dn(i,k)*ugkkdn
            fclon2m=forclo(i,n2m,k)*betadx
            amjl(i,n2m)=-am3j(n2m)*fclon2m
            acjl(i,n2m)=1.-ac3j(n2m)*fclon2m
            apjl(i,n2m)=-ap3j(n2m)*fclon2m
            iadd=i+(n2m-1)*n1m+(k-1)*n1m*n2m
            fjl(i,n2m)=rhs(iadd) + dqb3up(i,k)*ugkkup
          end do

c
      call trvb(amjl,acjl,apjl,fjl,m1,m2,n2m,n1m)
c
      do j=1,n2m
        do i=1,n1m
          iadd=i+(j-1)*n1m+(k-1)*n1m*n2m
          rhs(iadd) = fjl(i,j)
          end do
        end do
      end do
!$OMP  END PARALLEL DO

      return
      end

c
c     it inverts in 1 direction, vectors of n1m length
c     from 1, n2m volte
c
c     call trvpki(amil,acil,apil,fil,1,n1m,1,n2m,m1,m2)
c
c  ****************************** subrout trvpki  **********************
c
      subroutine trvpki (a,b,c,f,j1,j2,mi,mf,mm,nn)
c
      dimension    a(nn,mm),b(nn,mm),c(nn,mm),f(nn,mm)
     %             ,q(nn,mm),qe(nn,mm),s(nn,mm),fn(nn)
      ja = j1 + 1
      jj = j1 + j2
      do 20 k=mi,mf
      q(k,j1) = -c(k,j1)/b(k,j1)
      s(k,j1) = - a(k,j1)/b(k,j1)
      fn(k) = f(k,j2)
      f(k,j1) = f(k,j1)/b(k,j1)
   20 continue
c
c     forward elimination sweep
c
      do 10 j=ja,j2
      do 21 k=mi,mf
      p =1./( b(k,j) + a(k,j)*q(k,j-1))
      q(k,j) = - c(k,j)*p
      s(k,j) = - a(k,j)*s(k,j-1)*p
      f(k,j) = ( f(k,j) - a(k,j)*f(k,j-1))*p
   21 continue
   10 continue
c
c     backward pass
c
      do 22 k=mi,mf
      s(k,j2) = 1.
      qe(k,j2) = 0.
   22 continue
      do 11 i=ja,j2
      j = jj - i
      s(k,j) = s(k,j) + q(k,j)*s(k,j+1)
      do 23 k=mi,mf
      qe(k,j) = f(k,j) + q(k,j)*qe(k,j+1)
   23 continue
   11 continue
      do 24 k=mi,mf
      f(k,j2)=(fn(k) - c(k,j2)*qe(k,j1) - a(k,j2)*qe(k,j2-1))
     &       /(c(k,j2)*s(k,j1) + a(k,j2)*s(k,j2-1)  +b(k,j2))
   24 continue
c
c     backward elimination pass
c
      do 12 i=ja,j2
      j = jj -i
      do 25 k=mi,mf
      f(k,j) = f(k,j2)*s(k,j) + qe(k,j)
   25 continue
   12 continue
      return
      end
