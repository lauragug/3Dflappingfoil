      program media
      parameter (m2=145,m3=222,ndd=2,ndv=3)
      real*4 q1(m2,m3),q2(m2,m3),q3(m2,m3)
      real*4 pr(m2,m3)
      real*4 q1ave(m2,m3),q2ave(m2,m3),q3ave(m2,m3)
     %    ,prave(m2,m3)
      real*4 x2c(m2),x3c(m3)
      character*70 namefile(1000)
      character*70 namfi4
      real*4  epsil,re,time
      real*4  stime
      nfield = 0
      open(unit=53,file='fieldlist',form='formatted')
      do n=1,1000
        read(53,1001,end=2000) namefile(n)
        nfield = nfield + 1
      end do
      close(53)
 1001 format(a70)
 2000  print *,'total number of fields = ',nfield

      do l=1,nfield           !{ the do loop on the fields starts
      write(6,201) namefile(l)
  201 format(10x,'reading from ',a70)

      coef=1./dble(nfield)

        open(99,file=namefile(l),form='unformatted')
        rewind 99
        read(99) n2,n3,nnn
        read(99)
     %   ((q1(j,k),j=1,n2),k=1,n3),
     %   ((q2(j,k),j=1,n2),k=1,n3),
     %   ((q3(j,k),j=1,n2),k=1,n3),
     %   ((pr(j,k),j=1,n2),k=1,n3)
       close (99)
c  
        do k=1,n3
         do j=1,n2
          q1ave(j,k)=q1ave(j,k)+q1(j,k)
          q2ave(j,k)=q2ave(j,k)+q2(j,k)
          q3ave(j,k)=q3ave(j,k)+q3(j,k)
          prave(j,k)=prave(j,k)+pr(j,k)
         end do
        end do
        end do
c 
        do k=1,n3
         do j=1,n2
          q1ave(j,k)=q1ave(j,k)*coef
          q2ave(j,k)=q2ave(j,k)*coef
          q3ave(j,k)=q3ave(j,k)*coef
          prave(j,k)=prave(j,k)*coef
         end do
        end do

        open(99,file='avframe.dat',form='unformatted')
        rewind 99
        write(99) n3,n2,4
        write(99)
     %   ((q1ave(j,k),k=1,n3),j=1,n2),
     %   ((q2ave(j,k),k=1,n3),j=1,n2),
     %   ((q3ave(j,k),k=1,n3),j=1,n2),
     %   ((prave(j,k),k=1,n3),j=1,n2) 
       close (99)
       stop
       end 
