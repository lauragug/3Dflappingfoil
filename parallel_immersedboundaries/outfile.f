
      subroutine outpf1(time,comm)
      include'param.f'
      include "mpif.h"

      integer comm,ierr
      REAL    stime
      character*20 namfi3
      character*5 ipfi1
      character*2 ipfi2

      stime = time
      itime=nint(1000.*time)
      write(ipfi1,99)itime
      write(ipfi2,98)myid
 99   format(i5.5)
 98   format(i2.2)
      namfi3='debug1'//ipfi1//'-'//ipfi2//'.dat'
      open(99,file=namfi3,form='formatted',status='unknown')
      j = n2m/2
      do i = 1,n1m
         do k = 1,n3m
            write(99,100)dq(i,j,k),dph(i,j,k),qcap(i,j,k)
         enddo
      enddo
      close(99)
 100  format(3(e12.6,1x))

      return
      end
      
c ************************************************************************************

      subroutine outpf2(time,comm)
      include'param.f'
      include "mpif.h"

      integer comm,ierr
      REAL    stime
      character*20 namfi3
      character*5 ipfi1
      character*2 ipfi2

      stime = time
      itime=nint(1000.*time)
      write(ipfi1,99)itime
      write(ipfi2,98)myid
 99   format(i5.5)
 98   format(i2.2)
      namfi3='debug2'//ipfi1//'-'//ipfi2//'.dat'
      open(99,file=namfi3,form='formatted',status='unknown')
      j = n2m/2
      do i = 1,n1m
         do k = 1,n3m
            write(99,100)q3(i,j,k),q2(i,j,k),q1(i,j,k)
         enddo
      enddo
      close(99)
 100  format(3(e12.6,1x))

      return
      end

c ************************************************************************************

      subroutine outpf3(time,comm)
      include'param.f'
      include "mpif.h"

      integer comm,ierr
      REAL    stime
      character*20 namfi3
      character*5 ipfi1
      character*2 ipfi2

      stime = time
      itime=nint(1000.*time)
      write(ipfi1,99)itime
      write(ipfi2,98)myid
 99   format(i5.5)
 98   format(i2.2)
      namfi3='debug3'//ipfi1//'-'//ipfi2//'.dat'
      open(99,file=namfi3,form='formatted',status='unknown')
c      write(99,*) 'VARIABLES = "X","Y","Z","Vx","Vy","Vz",',
c     $     '"Pres"'
c      write(99,*) 'ZONE I = ',n3,', J = 1, K = ',n1, ', F=POINT'
c      j = n2m/2
c      do i = 1, n1
c        igeo = isx1 + i - 1
c        do k=1,n3
c          write(99,*)
c     $           x3c(k),x2c(j),x1c(igeo),
c     $       q3(i,j,k),q2(i,j,k),q1(i,j,k),pr(i,j,k)
c        enddo
c      enddo
      j = n2m/2
      do i = 1,n1m
         do k = 1,n3m
            write(99,100)qcap(i,j,k)
         enddo
      enddo
      close(99)
 100  format(3(e12.6,1x))

      return
      end

c ************************************************************************************

      subroutine outpf4(time,comm)
      include'param.f'
      include "mpif.h"

      integer comm,ierr
      REAL    stime
      character*20 namfi3
      character*5 ipfi1
      character*2 ipfi2

      stime = time
      itime=nint(1000.*time)
      write(ipfi1,99)itime
      write(ipfi2,98)myid
 99   format(i5.5)
 98   format(i2.2)
      namfi3='debug4'//ipfi1//'-'//ipfi2//'.dat'
      open(99,file=namfi3,form='formatted',status='unknown')
c      write(99,*) 'VARIABLES = "X","Y","Z","Vx","Vy","Vz",',
c     &'"Pres"'
c      write(99,*) 'ZONE I = ',n3,', J = 1, K = ',n1, ', F=POINT'
c      j = n2m/2
c      do i = 1, n1
c        igeo = isx1 + i - 1
c        do k=1,n3
c          write(99,*)
c     $           x3c(k),x2c(j),x1c(igeo),
c     $       q3(i,j,k),q2(i,j,k),q1(i,j,k),pr(i,j,k)
c        enddo
c      enddo
      j = n2m/2
      do i = 1,n1m
         do k = 1,n3m
            write(99,100)dph(i,j,k)
         enddo
      enddo
      close(99)
 100  format(3(e12.6,1x))

      return
      end

c ************************************************************************************

      subroutine outpf5(time,comm)
      include'param.f'
      include "mpif.h"

      integer comm,ierr
      REAL    stime
      character*20 namfi3
      character*5 ipfi1
      character*2 ipfi2

      stime = time
      itime=nint(1000.*time)
      write(ipfi1,99)itime
      write(ipfi2,98)myid
 99   format(i5.5)
 98   format(i2.2)
      namfi3='debug5'//ipfi1//'-'//ipfi2//'.dat'
      open(99,file=namfi3,form='formatted',status='unknown')
c      write(99,*) 'VARIABLES = "X","Y","Z","Vx","Vy","Vz",',
c     $'"Pres"'
c      write(99,*) 'ZONE I = ',n3,', J = 1, K = ',n1, ', F=POINT'
c      j = n2m/2
c      do i = 1, n1
c        igeo = isx1 + i - 1
c        do k=1,n3
c          write(99,*)
c     $           x3c(k),x2c(j),x1c(igeo),
c     $       q3(i,j,k),q2(i,j,k),q1(i,j,k),pr(i,j,k)
c        enddo
c      enddo
      j = n2m/2
      do i = 1,n1m
         do k = 1,n3m
            write(99,100)q3(i,j,k),q2(i,j,k),q1(i,j,k)
         enddo
      enddo
      close(99)
 100  format(3(e12.6,1x))

      return
      end

c ************************************************************************************

      subroutine outpf6(time,comm)
      include'param.f'
      include "mpif.h"

      integer comm,ierr
      REAL    stime
      character*20 namfi3
      character*5 ipfi1
      character*2 ipfi2

      stime = time
      itime=nint(1000.*time)
      write(ipfi1,99)itime
      write(ipfi2,98)myid
 99   format(i5.5)
 98   format(i2.2)
      namfi3='debug6'//ipfi1//'-'//ipfi2//'.dat'
      open(99,file=namfi3,form='formatted',status='unknown')
c      write(99,*) 'VARIABLES = "X","Y","Z","Vx","Vy","Vz",',
c     &'"Pres"'
c      write(99,*) 'ZONE I = ',n3,', J = 1, K = ',n1, ', F=POINT'
c      j = n2m/2
c      do i = 1, n1
c        igeo = isx1 + i - 1
c        do k=1,n3
c          write(99,*)
c     $           x3c(k),x2c(j),x1c(igeo),
c     $       q3(i,j,k),q2(i,j,k),q1(i,j,k),pr(i,j,k)
c        enddo
c      enddo
      j = n2m/2
      do i = 1,n1m
         do k = 1,n3m
            write(99,100)pr(i,j,k)
         enddo
      enddo
      close(99)
 100  format(3(e12.6,1x))

      return
      end

c ************************************************************************************

      subroutine outpf7(time,comm)
      include'param.f'
      include "mpif.h"

      integer comm,ierr
      REAL    stime
      character*20 namfi3
      character*5 ipfi1
      character*2 ipfi2

      stime = time
      itime=nint(1000.*time)
      write(ipfi1,99)itime
      write(ipfi2,98)myid
 99   format(i5.5)
 98   format(i2.2)
      namfi3='debug7'//ipfi1//'-'//ipfi2//'.dat'
      open(99,file=namfi3,form='formatted',status='unknown')
      do i = 1,m1*m2*m3
         write(99,100)rhs(i)
      enddo
      close(99)
 100  format(3(e12.6,1x))

      return
      end

c ************************************************************************************
