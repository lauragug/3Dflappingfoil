	subroutine exchng1( a, n1, n2, n3,comm, nbrleft, nbrright )
	include 'mpif.h'

	REAL a(0:n1,n2,n3)
	integer comm, nbrleft, nbrright, myid
	integer status(MPI_STATUS_SIZE), ierr
	integer stridetype
c
c        write(*,*)'exchng1-IN',myid
	call MPI_COMM_RANK( comm, myid, ierr )
	call MPI_TYPE_VECTOR(n2*n3,1,n1+1,MPI_DOUBLE_PRECISION,
     $     stridetype,ierr)
	call MPI_TYPE_COMMIT(stridetype,ierr)

	call MPI_SENDRECV( 
     &            a( n1-1 ,1,1),1,stridetype, nbrright, 0, 
     &            a( 0 ,1,1),1,stridetype, nbrleft, 0, 
     &            comm, status, ierr )

	call MPI_SENDRECV(      
     &            a( 1 ,1,1),1,stridetype, nbrleft, 0, 
     &            a( n1 ,1,1),1,stridetype, nbrright, 0, 
     &            comm, status, ierr )

        call MPI_TYPE_FREE(stridetype,ierr)
c        write(*,*)'exchng1-OUT',myid
	return
	end

c*************************************************************

      subroutine MPE_DECOMP1D( n, numprocs, myid, s, e )
      integer n, numprocs, myid, s, e
      integer nlocal
      integer deficit
c
      nlocal  = n / numprocs
      s	      = myid * nlocal + 1
      deficit = mod(n,numprocs)
      s	      = s + min(myid,deficit)
      if (myid .lt. deficit) then
          nlocal = nlocal + 1
      endif
      e = s + nlocal - 1
      if (e .gt. n .or. myid .eq. numprocs-1) e = n
      return
      end

c*************************************************************

      subroutine exchng2(a,temp,ny,sy,ey,nxglob,s,e,dir,comm)
      include 'mpif.h'

      parameter(maxnum=100)
      integer nxglob, ny, s, e, nx, srec, dir
      integer ii,jj,i,j
      integer sy,ey,sysend,eysend
      integer myid, numprocs,comm
      double precision a(0:e-s+2,ny), temp(nxglob+2,ey-sy+1)
      integer status(MPI_STATUS_SIZE,maxnum),ierr,req(maxnum)
c      integer status(MPI_STATUS_SIZE),ierr
      integer stridetype,stridetype1,irec

c      write(*,*)'exchng2-IN',myid
      call MPI_COMM_RANK( comm, myid, ierr )   
      call MPI_COMM_SIZE( comm, numprocs, ierr )

c      call MPI_TYPE_VECTOR(ey-sy+1,e-s+1,e-s+3,MPI_DOUBLE_PRECISION,
c     $     stridetype,ierr)
      call MPI_TYPE_VECTOR(ey-sy+1,e-s+1,nxglob+2,MPI_DOUBLE_PRECISION,
     $     stridetype1,ierr)

c      call MPI_TYPE_COMMIT(stridetype,ierr)
      call MPI_TYPE_COMMIT(stridetype1,ierr)

      irec=0
      nx=e-s+1
c      write(*,*)'inside exchng2_I',myid
      if(dir.eq.1)then
c      write(35+myid,*)'a',((a(i,j),i=s,e),j=1,ny)
         do jj = 1,ey-sy+1
            do ii = 1,e-s+1
               temp(ii+s,jj) = a(ii,jj-1+sy)      
            enddo         
         enddo
c         write(51+myid,*)'temp',((temp(ii,jj),ii=1,nxglob+2),jj=sy:ey)
         do 10 ii=0,numprocs-1
            if (ii.eq.myid) goto 10
            srec = 1 + ii*int(nxglob/numprocs)
c            sysend = 1 + ii*int(ny/numprocs)
            call MPE_DECOMP1D( ny, numprocs, ii, sysend, eysend )
c            write(*,*)'MODIFICA=',myid,ii,ny,srec,sysend,eysend,ey,sy
            call MPI_TYPE_VECTOR(eysend-sysend+1,e-s+1,e-s+3,
     $           MPI_DOUBLE_PRECISION,stridetype,ierr)
            call MPI_TYPE_COMMIT(stridetype,ierr)
c            call MPI_SENDRECV(a(s,sysend),1,stridetype,ii,0,
c     $           temp(srec+1,1),1, stridetype1,ii,0,
c     $           comm,status,ierr)
            irec=irec+1
            call MPI_ISEND(a(1,sysend),1,stridetype,ii,0,
     $           comm,req(irec),ierr)
            irec=irec+1
            call MPI_IRECV(temp(srec+1,1),1, stridetype1,ii,0,
     $           comm,req(irec),ierr)

            call MPI_TYPE_FREE(stridetype,ierr)
 10      continue

         call MPI_WAITALL(2*numprocs-2,req,status,ierr )

         do jj=1,ey-sy+1
            temp(1,jj)=temp(nxglob+1,jj)
            temp(nxglob+2,jj)=temp(2,jj)
         enddo
c        write(45+myid,*)'temp',myid,((temp(i,j),i=1,nxglob+2),j=sy,ey)
      else
c         write(21+myid,*)'temp',myid,((temp(i,j),i=1,nxglob+2),j=sy,ey)
         do jj = 1,ey-sy+1
            do ii = 1,e-s+1
               a(ii,jj-1+sy) = temp(ii+s,jj)    
            enddo
         enddo
         do 11 ii=0,numprocs-1
            if (ii.eq.myid) goto 11
            srec = 1 + ii*int(nxglob/numprocs)
c            sysend = 1 + ii*int(ny/numprocs)
            call MPE_DECOMP1D( ny, numprocs, ii, sysend, eysend )
c            write(*,*)'MODIFICA=',myid,ii,ny,srec,sysend,eysend,ey,sy
            call MPI_TYPE_VECTOR(eysend-sysend+1,e-s+1,e-s+3,
     $           MPI_DOUBLE_PRECISION,stridetype,ierr)
            call MPI_TYPE_COMMIT(stridetype,ierr)
c            call MPI_SENDRECV(temp(srec+1,1),1,stridetype1,ii,
c     $           0,a(s,sysend),1, stridetype,ii,0,
c     $           comm,status,ierr)
            irec=irec+1
            call MPI_ISEND(temp(srec+1,1),1,stridetype1,ii,0,
     $           comm,req(irec),ierr)
            irec=irec+1
            call MPI_IRECV(a(1,sysend),1, stridetype,ii,0,
     $           comm,req(irec),ierr)
            call MPI_TYPE_FREE(stridetype,ierr)
 11      continue
c         write(25+myid,*)'a',((a(i,j),i=s,e),j=1,ny)
         call MPI_WAITALL(2*numprocs-2,req,status,ierr )
      endif

c      write(61+myid,*)'myid',myid,((temp(i,j),i=1,nxglob),j=sy,ey)
	write(*,*)' VERIFICA - irec',irec

c      call MPI_TYPE_FREE(stridetype,ierr)
      call MPI_TYPE_FREE(stridetype1,ierr)
c      WRITE(*,*)'OUT OF HERE',MYID
c      write(*,*)'end exchng2_I',myid
c      write(*,*)'exchng2-OUT',myid
      return
      end

c*************************************************************

      subroutine exchng_wave(a,temp,ny,sy,ey,nxglob,sn,en,dir,comm)
      include 'mpif.h'

      parameter (maxnum = 100)
      integer nx, ny, sn, en, nxglob, dir
      integer ii,jj,i,j
      integer sy,ey,sysend,eysend,ssend,esend
      integer myid, numprocs,comm
      double precision a(en-sn+1,ny), temp(nxglob+2,ey-sy+1)
c      integer status(MPI_STATUS_SIZE),ierr
      integer status(MPI_STATUS_SIZE,maxnum),ierr,req(maxnum)
      integer stridetype,stridetype1,irec

c      write(*,*)'exchng_WAVE-IN',myid
      call MPI_COMM_RANK( comm, myid, ierr )   
      call MPI_COMM_SIZE( comm, numprocs, ierr )

c      call MPI_TYPE_VECTOR(ey-sy+1,en-sn+1,en-sn+1,MPI_DOUBLE_PRECISION,
c     $     stridetype,ierr)

c      call MPI_TYPE_COMMIT(stridetype,ierr)
      irec = 0
      if(dir.eq.1)then
c      write(75+myid,*)'a',dir,((a(i,j),i=1,en-sn+1),j=1,ny)
         do jj = 1,ey-sy+1
            do ii = 1,en-sn+1
               temp(sn+ii-1,jj) = a(ii,jj-1+sy) 
            enddo         
         enddo
         do 10 ii=0,numprocs-1
            if (ii.eq.myid) goto 10
c            sysend = 1 + ii*int(ny/numprocs)
            call MPE_DECOMP1D( ny, numprocs, ii, sysend, eysend )
            call MPI_TYPE_VECTOR(eysend-sysend+1,en-sn+1,en-sn+1,
     $           MPI_DOUBLE_PRECISION,stridetype,ierr)
            call MPI_TYPE_COMMIT(stridetype,ierr)
            call MPE_DECOMP1D( nxglob/2+1, numprocs, ii, ssend, esend )
            ssend = ssend*2 -1 
            esend = esend*2
            call MPI_TYPE_VECTOR(ey-sy+1,esend-ssend+1,nxglob+2,
     $           MPI_DOUBLE_PRECISION,stridetype1,ierr)
            call MPI_TYPE_COMMIT(stridetype1,ierr)
C            call MPI_SENDRECV(a(1,sysend),1,stridetype,ii,0,
C     $           temp(ssend,1),1, stridetype1,ii,0,
C     $           comm,status,ierr)

            irec = irec +1
            call MPI_ISEND(a(1,sysend),1,stridetype,ii,0,
     $           comm,req(irec),ierr)
            irec = irec +1
            call MPI_IRECV(temp(ssend,1),1, stridetype1,ii,0,
     $           comm,req(irec),ierr)

            call MPI_TYPE_FREE(stridetype,ierr)
            call MPI_TYPE_FREE(stridetype1,ierr)
 10      continue
         call MPI_WAITALL(2*numprocs-2,req,status,ierr )
c        write(65+myid,*)'temp',dir,
c     $        ((temp(i,j),i=1,nxglob+2),j=1,ey-sy+1)
      else
c         write(65+myid,*)'temp',dir,
c     $        ((temp(i,j),i=1,nxglob+2),j=1,ey-sy+1)
         do jj = 1,ey-sy+1
            do ii = 1,en-sn+1
               a(ii,jj-1+sy) = temp(sn+ii-1,jj)    
            enddo
         enddo
         do 11 ii=0,numprocs-1
            if (ii.eq.myid) goto 11
c            sysend = 1 + ii*int(ny/numprocs)
            call MPE_DECOMP1D( ny, numprocs, ii, sysend, eysend )
            call MPI_TYPE_VECTOR(eysend-sysend+1,en-sn+1,en-sn+1,
     $           MPI_DOUBLE_PRECISION,stridetype,ierr)
            call MPI_TYPE_COMMIT(stridetype,ierr)

            call MPE_DECOMP1D( nxglob/2+1, numprocs, ii, ssend, esend )
            ssend = ssend*2 -1 
            esend = esend*2
            call MPI_TYPE_VECTOR(ey-sy+1,esend-ssend+1,nxglob+2,
     $           MPI_DOUBLE_PRECISION,stridetype1,ierr)
            call MPI_TYPE_COMMIT(stridetype1,ierr)
c            write(*,*)'WAVE',myid,ii,ssend,esend,sn,en,
c     $           sysend,eysend,sy,ey
C            call MPI_SENDRECV(temp(ssend,1),1,stridetype1,ii,
C     $           0,a(1,sysend),1, stridetype,ii,0,
C     $           comm,status,ierr)

            irec = irec +1
            call MPI_ISEND(temp(ssend,1),1,stridetype1,ii,0,
     $           comm,req(irec),ierr)
            irec = irec +1
            call MPI_IRECV(a(1,sysend),1, stridetype,ii,0,
     $           comm,req(irec),ierr)

            call MPI_TYPE_FREE(stridetype,ierr)
            call MPI_TYPE_FREE(stridetype1,ierr)
 11      continue
         call MPI_WAITALL(2*numprocs-2,req,status,ierr )
c         write(75+myid,*)'a',dir,((a(i,j),i=1,en-sn+1),j=1,ny)
      endif
c      write(*,*)'exchng_WAVE-out',myid
c      call MPI_TYPE_FREE(stridetype,ierr)

      return
      end

c*************************************************************

      subroutine exchng3(a,temp,ny,sy,ey,nxglob,s,e,dir,comm)
      include 'mpif.h'

      parameter(maxnum=100)
      integer nxglob, ny, s, e, nx, srec, dir
      integer ii,jj,i,j
      integer sy,ey,sysend,eysend
      integer myid, numprocs,comm
      double precision a(ny+1,e-s+1), temp(ey-sy+2,nxglob+1)
      integer status(MPI_STATUS_SIZE,maxnum),ierr,req(maxnum)
      integer stridetype,stridetype1,irec

c      write(*,*)'exchng3-IN',myid 
      call MPI_COMM_RANK( comm, myid, ierr )   
      call MPI_COMM_SIZE( comm, numprocs, ierr )

c      write(*,*)' error ',myid,nxglob,e,s,ny,ey,sy
      call MPI_TYPE_VECTOR(e-s+1,ey-sy+1,ey-sy+2,MPI_DOUBLE_PRECISION,
     $     stridetype1,ierr)

      call MPI_TYPE_COMMIT(stridetype1,ierr)

      irec=0
      nx=e-s+1
      if(dir.eq.1)then

         do jj = 1,ey-sy+1
            do ii = 1,e-s+1
               temp(jj,s+ii-1) = a(jj-1+sy,ii)      
            enddo         
         enddo

         do 10 ii=0,numprocs-1
            if (ii.eq.myid) goto 10
            srec = 1 + ii*int(nxglob/numprocs)

            call MPE_DECOMP1D( ny, numprocs, ii, sysend, eysend )
C            write(*,*)' error1 ',myid,e,s,eysend,sysend,ny
            call MPI_TYPE_VECTOR(e-s+1,eysend-sysend+1,ny+1,
     $           MPI_DOUBLE_PRECISION,stridetype,ierr)
            call MPI_TYPE_COMMIT(stridetype,ierr)
            irec=irec+1
            call MPI_ISEND(a(sysend,1),1,stridetype,ii,0,
     $           comm,req(irec),ierr)
            irec=irec+1
            call MPI_IRECV(temp(1,srec),1, stridetype1,ii,0,
     $           comm,req(irec),ierr)

            call MPI_TYPE_FREE(stridetype,ierr)
 10      continue

         call MPI_WAITALL(2*numprocs-2,req,status,ierr )

      else

         do jj = 1,ey-sy+1
            do ii = 1,e-s+1
               a(jj-1+sy,ii) = temp(jj,s+ii-1)    
            enddo
         enddo
         do 11 ii=0,numprocs-1
            if (ii.eq.myid) goto 11
            srec = 1 + ii*int(nxglob/numprocs)

            call MPE_DECOMP1D( ny, numprocs, ii, sysend, eysend )

            call MPI_TYPE_VECTOR(e-s+1,eysend-sysend+1,ny+1,
     $           MPI_DOUBLE_PRECISION,stridetype,ierr)
            call MPI_TYPE_COMMIT(stridetype,ierr)

            irec=irec+1
            call MPI_ISEND(temp(1,srec),1,stridetype1,ii,0,
     $           comm,req(irec),ierr)
            irec=irec+1
            call MPI_IRECV(a(sysend,1),1, stridetype,ii,0,
     $           comm,req(irec),ierr)
            call MPI_TYPE_FREE(stridetype,ierr)
 11      continue

         call MPI_WAITALL(2*numprocs-2,req,status,ierr )
      endif

      call MPI_TYPE_FREE(stridetype1,ierr)
c      write(*,*)'exchng3-OUT',myid 
      return
      end

c*************************************************************
