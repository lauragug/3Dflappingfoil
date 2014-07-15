c************************************************************************
c                                                                       *
c  ****************************** subrout mkmov  ********************** *
c                                                                       *
c************************************************************************
      subroutine mkmov(time,q1,q2,q3,pr)
      include 'param.f'
      REAL*4 q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)
      REAL   pr(m1,m2,m3)
      REAL*4 voc1(m2,m3),voc2(m2,m3),voc3(m2,m3),prc(m2,m3)
      character*70 namevec,namfi3
      character*70 stringa
      character*5 ipfi,ipre
c
c     form the name of the file
c
      tprfi = 1000.
      irep=re
      itime=nint(time*tprfi)
      write(ipre,98)irep
      write(ipfi,82)itime
   98 format(i5.5)
   82 format(i5.5)
c
c     SELECT THE 2D SLICE
c
      ic=n1m/2+1                                                     
      im=imv(ic)
      do kc=1,n3m
        km=kmv(kc)
        do jc=1,n2m
          jm=jmv(jc)
          voc2(jc,kc)=(q2(ic,jc,kc)+q2(im,jc,kc)+
     %                    q2(ic,jc,km)+q2(im,jc,km))*0.25
          voc1(jc,kc)=(q1(ic,jc,kc)+q1(ic,jm,kc)+
     %                    q1(ic,jc,km)+q1(ic,jm,km))*0.25
          voc3(jc,kc)=(q3(ic,jc,kc)+q3(im,jc,kc)+
     %                    q3(ic,jm,kc)+q3(im,jm,kc))*0.25
          prc(jc,kc)=(pr(ic,jc,kc)+pr(im,jc,kc)
     %                  +pr(ic,jm,kc)+pr(im,jm,kc)
     %                  +pr(ic,jc,km)+pr(im,jc,km)
     %                  +pr(ic,jm,km)+pr(im,jm,km))*0.125
        end do
      end do
c
c
      namfi3='fv'//ipfi//'.dat'
      open(59,file=namfi3,form='unformatted')
      rewind 59
      write(59) n2,n3,4
      write(59) ((voc1(j,k),j=1,n2),k=1,n3)            
     %         ,((voc2(j,k),j=1,n2),k=1,n3)            
     %         ,((voc3(j,k),j=1,n2),k=1,n3)            
     %         ,((prc(j,k),j=1,n2),k=1,n3)            
      close(59)
      return                                                            
      end                                                               
