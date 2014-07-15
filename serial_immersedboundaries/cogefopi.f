c
c***********************************************************************
c                                                                      *
c   ************** TOPOGR ************************                     *
c                                                                      *
c***********************************************************************
      subroutine topogr
c
c     This routine finds the indices in the computational grid
c     close to the physical location of the body.
c
      include 'param.f'
      common/parsta/ npunr_st, npunz_st
      common/parleb/ p_es,p_en,d_ls,p_pl,p_cs
c      common/inter/ indgeoi(3,mpun,3)
      common/inter/ indgeoi(3,1600000,3)
c     write(6,* ) ' a '
c
      open(60,file='GRID.data',form='formatted',status='old')
      read(60,*) nz,ny,nx
      n3=nz
      n2=ny
      n1=nx
      n1m=n1-1
      n2m=n2-1
      n3m=n3-1
      if(n1.eq.1) n1m = 1
      if(n1m.eq.1) then
        iaxsy=0
      else
        iaxsy=1
      endif
      if(nz.gt.m3) then
        write(6,*) ' WARNING: nz gt m3,  nz = ',nz,' m3 = ',m3
        pause
      end if
      if(ny.gt.m2) then
        write(6,*) ' WARNING: ny gt m2 '
        pause
      end if
      if(nx.gt.m1) then
        write(6,*) ' WARNING: nx gt m1 '
        pause
      end if
c
      read(60,*) (x3c(k),k=1,nz)
      read(60,*) (x2c(j),j=1,ny)
      read(60,*) (x1c(i),i=1,nx)
c     do k=1,nz
c        read(60,*) x3c(k)
c        write(6,*) x3c(k)
c     enddo
c     do j=1,ny
c        read(60,*) x2c(j)
c        write(6,*) x2c(j)
c     enddo
c     do i=1,nx
c        read(60,*) x1c(i)
c        write(6,*) x1c(i)
c     enddo
c
c     q3
c
      read(60,*) npunz_i
      write(6,*) ' Number of internal points for q3 = ',npunz_i
      if(npunz_i.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for q3',mpun
c     read(60,*) (indgeoi(3,n,3),indgeoi(3,n,2),indgeoi(3,n,1)
c    %            ,n=1,npunz_i)
      do n=1,npunz_i
      read(60,*) indgeoi(3,n,3),indgeoi(3,n,2),indgeoi(3,n,1)
      enddo            
      read(60,*) npunz_e
      alx3=x3c(nz)-x3c(1)
      write(6,*) ' Number of external points for q3 = ',npunz_e
      if(npunz_e.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for q3',mpun
c     read(60,*) (indgeo(3,n,3),indgeo(3,n,2),indgeo(3,n,1),
c    %            indgeoe(3,n,3),indgeoe(3,n,2),indgeoe(3,n,1),
c    %            distb(3,n),n=1,npunz_e)
      do n=1,npunz_e 
      read(60,*) indgeo(3,n,3),indgeo(3,n,2),indgeo(3,n,1),
     %           indgeoe(3,n,3),indgeoe(3,n,2),indgeoe(3,n,1),
     %           distb(3,n)
      enddo
c
      npunz_i = 0  ! PRESSURE
      npunz = npunz_e + npunz_i
      if(npunz.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for q3',mpun
      do n=1,npunz_i
        ll = npunz_e+n 
        indgeo(3,ll,1) = indgeoi(3,n,1)
        indgeo(3,ll,2) = indgeoi(3,n,2)
        indgeo(3,ll,3) = indgeoi(3,n,3)
        indgeoe(3,ll,1) = indgeoi(3,n,1)
        indgeoe(3,ll,2) = indgeoi(3,n,2)
        indgeoe(3,ll,3) = indgeoi(3,n,3)
        distb(3,ll) = 0.
      end do
c
c     q2
c
      read(60,*) npunr_i
      write(6,*) ' Number of internal points for q2 = ',npunr_i
      if(npunr_i.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for q2',mpun
c     read(60,*) (indgeoi(2,n,3),indgeoi(2,n,2),indgeoi(2,n,1)
c    %            ,n=1,npunr_i)
      do n=1,npunr_i
      read(60,*) indgeoi(2,n,3),indgeoi(2,n,2),indgeoi(2,n,1)
      enddo
      read(60,*) npunr_e
      write(6,*) ' Number of external points for q2 = ',npunr_e
      if(npunr_e.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for q2',mpun
c     read(60,*) (indgeo(2,n,3),indgeo(2,n,2),indgeo(2,n,1),
c    %            indgeoe(2,n,3),indgeoe(2,n,2),indgeoe(2,n,1),
c    %            distb(2,n),n=1,npunr_e)
      do n=1,npunr_e
      read(60,*) indgeo(2,n,3),indgeo(2,n,2),indgeo(2,n,1),
     %            indgeoe(2,n,3),indgeoe(2,n,2),indgeoe(2,n,1),
     %            distb(2,n)
      enddo
c
      npunr_i = 0 !PRESSURE
      npunr = npunr_e + npunr_i
      if(npunr.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for q2',mpun
      do n=1,npunr_i
        ll = npunr_e+n 
        indgeo(2,ll,1) = indgeoi(2,n,1)
        indgeo(2,ll,2) = indgeoi(2,n,2)
        indgeo(2,ll,3) = indgeoi(2,n,3)
        indgeoe(2,ll,1) = indgeoi(2,n,1)
        indgeoe(2,ll,2) = indgeoi(2,n,2)
        indgeoe(2,ll,3) = indgeoi(2,n,3)
        distb(2,ll) = 0.
      end do
c
c     q1
c
      read(60,*) npunt_i
      write(6,*) ' Number of internal points for q1 = ',npunt_i
      if(npunt_i.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for q1',mpun
c     read(60,*) (indgeoi(1,n,3),indgeoi(1,n,2),indgeoi(1,n,1)
c    %            ,n=1,npunt_i)
      do n=1,npunt_i
      read(60,*) indgeoi(1,n,3),indgeoi(1,n,2),indgeoi(1,n,1)
      enddo
      read(60,*) npunt_e
      write(6,*) ' Number of external points for q1 = ',npunt_e
      if(npunt_e.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for q1',mpun
c     read(60,*) (indgeo(1,n,3),indgeo(1,n,2),indgeo(1,n,1),
c    %            indgeoe(1,n,3),indgeoe(1,n,2),indgeoe(1,n,1),
c    %            distb(1,n),n=1,npunt_e)
      do n=1,npunt_e
      read(60,*) indgeo(1,n,3),indgeo(1,n,2),indgeo(1,n,1),
     %            indgeoe(1,n,3),indgeoe(1,n,2),indgeoe(1,n,1),
     %            distb(1,n)
      enddo
c
      npunt_i = 0 !PRESSURE
      npunt = npunt_e + npunt_i
      if(npunt.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for q1',mpun
      do n=1,npunt_i
        ll = npunt_e+n 
        indgeo(1,ll,1) = indgeoi(1,n,1)
        indgeo(1,ll,2) = indgeoi(1,n,2)
        indgeo(1,ll,3) = indgeoi(1,n,3)
        indgeoe(1,ll,1) = indgeoi(1,n,1)
        indgeoe(1,ll,2) = indgeoi(1,n,2)
        indgeoe(1,ll,3) = indgeoi(1,n,3)
        distb(1,ll) = 0.
      end do
c
      close(60)
c
      return
      end
c
