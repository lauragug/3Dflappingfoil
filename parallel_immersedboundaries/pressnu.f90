!************************************************************************
!  this subroutine perform the calculation of trigz for temperton fft
!
      subroutine fftqua
      include'param.f'

      pi=2.*asin(1.)
!
!     wave number definition
!
      n1mh=n1mglob/2+1
      n1mp=n1mh+1
      do 18 i=1,n1mh
        ap(i)=(i-1)*2.*pi
   18 continue
      do 19 i=n1mp,n1mglob
        ap(i)=-(n1mglob-i+1)*2.*pi
   19 continue
      nx1fft=n1mglob
      if(n1mglob.ne.1) call fftfax(nx1fft,ifx1,trigx1)
!
!   modified wave number
!
      do 28 i=1,n1mglob
        ak1(i)=2.*(1.-cos(ap(i)/n1mglob))*dx1q
   28 continue
!
      return
      end subroutine fftqua

!***********************************************************************
!  this subroutine perform the calculation of dph , periodic direction
!  along x3 and x1to use the real fourier transform
!
      subroutine phcalc(comm)
      include'param.f'
      include "mpif.h"

      real yfisl(m3,m2),wl1(m3*m2),wl2(m3*m2)
      real amjl(m2),acjl(m2),apjl(m2)  !inutilizzati
      real amkl(m3),ackl(m3),apkl(m3)
!      dimension xr(m1m+2,m2m),work(m1,m2m)   !xr dovrebber essere sostituito da temp
      real   etime,t(2)  !t(2)inutilizzato
      integer sy, ey , snew, enew, snew1, enew1
      integer comm, ierr

      real temp1(0:n1,n2m)
      real,    dimension(:,:), allocatable :: atrig
      real,    dimension(:,:), allocatable :: temp,work
      real,    dimension(:,:,:), allocatable :: qcap_q,dph_d
      real dphno

      my = m3
      pi = 2.*asin(1.)

      call MPE_DECOMP1D( n2m, numprocs, myid, sy, ey )
      allocate(temp(n1mglob+2,ey-sy+1))
      allocate(work(n1mglob+1,ey-sy+1))
!
!   fft applied to the x1 direction to the
!   complex coeff. from cos fft
!   from physical to wave number space
!
      phtim1 = etime(t)
!      n1mh=n1m/2+1
!      n1md=n1m+2
      call MPE_DECOMP1D( n1mglob/2+1, numprocs, myid, snew, enew )
      n1mh= enew - snew + 1            !numero armoniche per proc
      snew1 = snew*2 -1 
      enew1 = enew*2
      allocate(atrig(enew1-snew1+1,n2m))
      allocate(qcap_q(n1mh*2,n2m,n3m))
      allocate(dph_d(n1mh*2,n2m,n3m))

      do k = 1,n3m
         fftloop1=etime(t)
         do j = 1,n2m
            do i= 1,n1m
               temp1(i,j)=qcap(i,j,k)
            enddo
         enddo
         call exchng2(temp1,temp,n2m,sy,ey,n1mglob,isx1,iex1,1,comm)
!         do j=1,n2m               !parte seriale
!            xr(1,j)=qcap(n1m,j,k)
!            xr(n1md,j)=qcap(1,j,k)
!            do i=1,n1m
!               is=i+1
!               xr(is,j)=qcap(i,j,k)
!            end do
!         end do
        
         fftin1=etime(t)
!         call fft99(xr,work,trigx1,ifx1,1,m1+1,n1m,n2m,-1) !parte seriale
         call fft99(temp,work,trigx1,ifx1,1,m1+1,n1mglob,ey-sy+1,-1)
         fftout1=etime(t)

         call exchng_wave(atrig,temp,n2m,sy,ey,n1mglob,snew1,enew1,-1,comm)
         
!         do j=1,n2m                    !parte seriale
!            do i=1,n1mh
!               ip=2*i
!               id=2*i-1
!               qcap(i,j,k)=xr(id,j)
!               dph(i,j,k)=xr(ip,j)
!            end do
!         end do
         do j = 1,n2m
            do i = 1,n1mh
               ip=2*i
               id=2*i-1
               qcap_q(i,j,k) = atrig(id,j)
               dph_d(i,j,k) = atrig(ip,j)
            enddo
         enddo  
         fftloop2=etime(t)
      end do
!     
!     inversion of the matrix in the x2 and x3 directions (real part)
!     
      phtim2=etime(t)
!      do k=1,n3m
!         amkl(k)=amphk(k)
!         apkl(k)=apphk(k)
!      enddo
      do i=1,n1mh
         blkloop1=etime(t)
         wl1=w
         wl2=w
         ii = i + snew -1
         do k=1,n3m
            ackl(k)=acphkk(ii,k)! modifica nel caso piu' processori
            amkl(k)=amphk(k)
            apkl(k)=apphk(k)
            do j=1,n2m
               yfisl(k,j)=qcap_q(i,j,k)
!    amjl(j)=amphj(j)
!    acjl(j)=acphj(j)
!    apjl(j)=apphj(j)
            end do
         end do

         blkin1=etime(t)
         
         call blktri(1,1,n2m,amphj,acphj,apphj,1,n3m,amkl,ackl,apkl,m3,yfisl,ierror,wl1)
!     call blktri(1,1,n2m,amjl,acjl,apjl
!     %             ,1,n3m,amkl,ackl,apkl,m3,yfisl,ierror,wl1)     !seriale

         blkout1=etime(t)

         do k=1,n3m
            do j=1,n2m
               qcap_q(i,j,k)=yfisl(k,j)
            end do
         end do
!    
!     inversion of the matrix in the x2 and x3 directions (imaginary part)
!     
         do k=1,n3m
            do j=1,n2m
               yfisl(k,j)=dph_d(i,j,k)
            end do
         end do

         blkin2=etime(t)
         call blktri(1,1,n2m,amphj,acphj,apphj,1,n3m,amkl,ackl,apkl,m3,yfisl,ierror,wl2)
         blkout2=etime(t)
!       call blktri(1,1,n2m,amjl,acjl,apjl
!    %             ,1,n3m,amkl,ackl,apkl,m3,yfisl,ierror,wl2)   !seriale

         do k=1,n3m
            do j=1,n2m
               dph_d(i,j,k)=yfisl(k,j)
            end do
         end do
         blkloop2=etime(t)
      end do
      phtim3=etime(t)
!     print *,'maxval qcap',maxval(qcap(:,:,:)),minval(qcap(:,:,:))
!     print *,'maxval dph',maxval(dph(:,:,:)),minval(dph(:,:,:))
!    
!   inverse fft applied to the phi x1 direction
!   from wave number space to physical space
!    
!      n1mu=n1m-1
      do k=1,n3m
         fftloop3=etime(t)
         do j=1,n2m
            do i=1,n1mh
               ip=2*i
               id=2*i-1
!              xr(id,j)=qcap(i,j,k)          !seriale
!              xr(ip,j)=dph(i,j,k)
               atrig(id,j) = qcap_q(i,j,k)
               atrig(ip,j) = dph_d(i,j,k)
            enddo
         enddo

         call exchng_wave(atrig,temp,n2m,sy,ey,n1mglob,snew1,enew1,1,comm)         
         fftin2=etime(t)

         call fft99(temp,work,trigx1,ifx1,1,m1+1,n1mglob,ey-sy+1,+1)
         fftout2=etime(t)

         call exchng2(temp1,temp,n2m,sy,ey,n1mglob,isx1,iex1,-1,comm)        

!         do j=1,n2m               !seriale
!            dph(n1m,j,k)=xr(1,j)
!            do i=1,n1mu
!               is=i+1
!               dph(i,j,k)=xr(is,j)
!            enddo
!         enddo
         do j = 1,n2m
            do i=1,n1m
               dph(i,j,k) = temp1(i,j)
            enddo
         enddo
         fftloop4=etime(t)
      enddo
!
      deallocate(qcap_q)
      deallocate(dph_d)
      deallocate(atrig)
      deallocate(temp)
      deallocate(work)
!      dphno = dph(1,1,1)
      if (myid.eq.0)dphno=dph(1,1,1)
      call MPI_BCAST(dphno,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      do k=1,n3m
         do j=1,n2m
            do i=1,n1m
               dph(i,j,k)=dph(i,j,k)-dphno
            enddo
         enddo
      enddo
      phtim4=etime(t)
!      print *,' in phcalc =',phtim2-phtim1,-phtim2+phtim3,phtim4-phtim3
!      print *,' in fft99 -1,loop,+1,loop=',fftout1-fftin1,fftloop2-fftloop1,fftout2-fftin2,fftloop4-fftloop3 
!      print *,' in blktri r,i,loop =',blkout1-blkin1,blkout2-blkin2,blkloop2-blkloop1

      return
      end subroutine phcalc

!***********************************************************************
!   in this subr the coefficients of the poisson eq. for dph
!   are calculated this subr. is called only at the beginning
!
      subroutine phini
      include'param.f'
!      real,    dimension(:), allocatable :: ww

      call fftqua
      n1mh=n1mglob/2+1
!
!   tridiagonal matrix coefficients at each k and i
!   x1 and x3 cartesian coordinates
!
      do jc=1,n2m
        jp=jpv(jc)
        a22icc=jmc(jc)*dx2q/g2c(jc)
        a22icp=jpc(jc)*dx2q/g2c(jp)
        ac2=-(a22icc+a22icp)
        ugmmm=1./g2m(jc)
        amphj(jc)=a22icc*ugmmm
        apphj(jc)=a22icp*ugmmm
        acphj(jc)=-(amphj(jc)+apphj(jc))
      end do
!
      do kc=1,n3m
        km=kmv(kc)
        kp=kpv(kc)
        a33icc=kmc(kc)*dx3q/g3c(kc)
        a33icp=kpc(kc)*dx3q/g3c(kp)
        ugmmm=1./g3m(kc)
        amphk(kc)=a33icc*ugmmm
        apphk(kc)=a33icp*ugmmm
        acphk(kc)=-(amphk(kc)+apphk(kc))
      enddo
      do ic=1,n1mh
        do kc=1,n3m
          acphkk(ic,kc)=acphk(kc)-ak1(ic)*iaxsy
        enddo
      enddo
      an3=n2m
      ax=alog(an3)/alog(2.)
      k=ax+1
      L=2**(K+1)
      mw=(K-2)*L+K+5+MAX(2*n2m,6*n3m)
      if(mw.gt.(m3*m2)) then
        write(6,*) ' WARNING THE ARRAY W HAS INSUFFICIENT '
        write(6,*) ' DIMENSIONS : mw = ',mw
        pause
      end if
     
      np=1
      mp=1
      do j=1,n2m
        do k=1,n3m
          yfis(k,j)=1.
        enddo
      enddo
      my=m3
      call blktri(0,np,n2m,amphj,acphj,apphj,mp,n3m,amphk,acphk,apphk,my,yfis,ierror,w)
!     write(6,*) ' 1 IERROR = ',ierror 

      return
      end subroutine phini
