************************************************************************
c  this subroutine perform the calculation of trigz for temperton fft
c
      subroutine fftqua
      include'param.f'
      pi=2.*asin(1.)
c
c     wave number definition
c
      n1mh=n1m/2+1
      n1mp=n1mh+1
      do 18 i=1,n1mh
        ap(i)=(i-1)*2.*pi
   18 continue
      do 19 i=n1mp,n1m
        ap(i)=-(n1m-i+1)*2.*pi
   19 continue
      nx1fft=n1m
      if(n1m.ne.1) call fftfax(nx1fft,ifx1,trigx1)
c
c   modified wave number
c
      do 28 i=1,n1m
        ak1(i)=2.*(1.-cos(ap(i)/n1m))*dx1q
   28 continue
c
      return
      end
************************************************************************
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 and x1to use the real fourier transform
c
      subroutine phcalc
      include'param.f'
      dimension yfisl(m3,m2),wl1(m3*m2),wl2(m3*m2)
      dimension amjl(m2),acjl(m2),apjl(m2)
      dimension amkl(m3),ackl(m3),apkl(m3)

      dimension xr(m1m+2,m2m),work(m1,m2m)
      real*4 etime,t(2)
      my=m3
c
      pi=2.*asin(1.)

c
c
c   fft applied to the x1 direction to the
c   complex coeff. from cos fft
c   from physical to wave number space
c
c     write(6,*) ' aa'
      phtim=etime(t)
      n1mh=n1m/2+1
      n1md=n1m+2
!$OMP  PARALLEL DO
!$OMP$ SHARED(qcap,dph,n1m,n3m,n2m,n1mh,trigx1,ifx1)
!$OMP$ PRIVATE(xr,work,is,ip,id)
        do k=1,n3m
          do j=1,n2m
            xr(1,j)=qcap(n1m,j,k)
            xr(n1md,j)=qcap(1,j,k)
            do i=1,n1m
              is=i+1
              xr(is,j)=qcap(i,j,k)
            end do
          end do

          call fft99(xr,work,trigx1,ifx1,1,m1+1,n1m,n2m,-1)

          do j=1,n2m
            do i=1,n1mh
              ip=2*i
              id=2*i-1
              qcap(i,j,k)=xr(id,j)
              dph(i,j,k)=xr(ip,j)
            end do
          end do
        end do
!$OMP  END PARALLEL DO
c     write(6,*) ' ab'
c    
c     inversion of the matrix in the x2 and x3 directions (real part)
c    
      phtim2=etime(t)
!$OMP  PARALLEL DO
!$OMP$ SHARED(acphkk,amphk,apphk,amphj,acphj,apphj)
!$OMP$ SHARED(qcap,dph,n1mh,n3m,n2m,w) 
!$OMP$ PRIVATE(yfisl,ackl,amkl,apkl,amjl,acjl,apjl,ierror)
!$OMP$ PRIVATE(wl1,wl2)

      do i=1,n1mh

      wl1=w
      wl2=w
 
        do k=1,n3m
          ackl(k)=acphkk(i,k)
          amkl(k)=amphk(k)
          apkl(k)=apphk(k)
          do j=1,n2m
            yfisl(k,j)=qcap(i,j,k)
c           amjl(j)=amphj(j)
c           acjl(j)=acphj(j)
c           apjl(j)=apphj(j)
          end do
        end do

c      print *,'calling blktri real from proc',i,OMP_GET_THREAD_NUM()

        call blktri(1,1,n2m,amphj,acphj,apphj
     %             ,1,n3m,amkl,ackl,apkl,m3,yfisl,ierror,wl1)
c       call blktri(1,1,n2m,amjl,acjl,apjl
c    %             ,1,n3m,amkl,ackl,apkl,m3,yfisl,ierror,wl1)

        do k=1,n3m
          do j=1,n2m
            qcap(i,j,k)=yfisl(k,j)
          end do
        end do
c    
c     inversion of the matrix in the x2 and x3 directions (imaginary part)
c    
          do k=1,n3m
            do j=1,n2m
              yfisl(k,j)=dph(i,j,k)
            end do
          end do

c      print *,'calling blktri imag. from proc',i,OMP_GET_THREAD_NUM()

        call blktri(1,1,n2m,amphj,acphj,apphj
     %             ,1,n3m,amkl,ackl,apkl,m3,yfisl,ierror,wl2)

c       call blktri(1,1,n2m,amjl,acjl,apjl
c    %             ,1,n3m,amkl,ackl,apkl,m3,yfisl,ierror,wl2)

          do k=1,n3m
            do j=1,n2m
              dph(i,j,k)=yfisl(k,j)
            end do
          end do
      end do
!$OMP END PARALLEL DO
      phtim3=etime(t)
c     print *,'maxval qcap',maxval(qcap(:,:,:)),minval(qcap(:,:,:))
c     print *,'maxval dph',maxval(dph(:,:,:)),minval(dph(:,:,:))

c     write(6,*) ' ae'
c    
c   inverse fft applied to the phi x1 direction
c   from wave number space to physical space
c    
        n1mu=n1m-1
!$OMP  PARALLEL DO
!$OMP$ SHARED(qcap,dph,n1m,n1mu,n3m,n2m,n1mh,trigx1,ifx1)
!$OMP$ PRIVATE(xr,work,is,ip,id)
        do k=1,n3m
          do j=1,n2m
            do i=1,n1mh
              ip=2*i
              id=2*i-1
              xr(id,j)=qcap(i,j,k)
              xr(ip,j)=dph(i,j,k)
            enddo
          enddo
          call fft99(xr,work,trigx1,ifx1,1,m1+1,n1m,n2m,+1)
          do j=1,n2m
            dph(n1m,j,k)=xr(1,j)
            do i=1,n1mu
              is=i+1
              dph(i,j,k)=xr(is,j)
            enddo
          enddo
        enddo
!$OMP  END PARALLEL DO
c
        dphno = dph(1,1,1)
!$OMP  PARALLEL DO
!$OMP$ SHARED(dph,n1m,n3m,n2m,dphno)
        do k=1,n3m
          do j=1,n2m
            do i=1,n1m
              dph(i,j,k)=dph(i,j,k)-dphno
            enddo
          enddo
        enddo
!$OMP  END PARALLEL DO
c
c     write(6,*) ' af'
c
      phtim4=etime(t)
C     print *,' in phcalc =',-phtim+phtim4,phtim2-phtim,
C    %                       -phtim2+phtim3,phtim4-phtim3
      return
      end
************************************************************************
c   in this subr the coefficients of the poisson eq. for dph
c   are calculated this subr. is called only at the beginning
c
      subroutine phini
      include'param.f'
      call fftqua
      n1mh=n1m/2+1
c
c   tridiagonal matrix coefficients at each k and i
c   x1 and x3 cartesian coordinates
c
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
c
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
      call blktri(0,np,n2m,amphj,acphj,apphj
     %             ,mp,n3m,amphk,acphk,apphk,my,yfis,ierror,w)
c     write(6,*) ' 1 IERROR = ',ierror 
      return
      end
