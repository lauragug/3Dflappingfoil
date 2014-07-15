!************************************************************************
!
!   this subroutine performs the inversion of the q1 momentum equation
!   by a factored implicit scheme, only the derivatives 11,22,33 of q1
!   are treated implicitly
!         direction x1
!
      subroutine solq1i(comm)
      include'param.f'
      
      real amil(m2,max_m1-1),acil(m2,max_m1-1),apil(m2,max_m1-1),fil(m2,max_m1-1)
      real,    dimension(:,:), allocatable :: amilP,acilP,apilP,filP
      integer comm, ey, sy
!
      call MPE_DECOMP1D( n2m, numprocs, myid, sy, ey )
      allocate(amilP(ey-sy+2,n1glob))
      allocate(acilP(ey-sy+2,n1glob))
      allocate(apilP(ey-sy+2,n1glob))
      allocate(filP(ey-sy+2,n1glob))

      do k=1,n3
        do j=1,n2
          do i=1,n1
            forclo(i,j,k)=1.
          end do
        end do
      end do

      if(infig.ne.-1) then
       do n=1,npunt
!          i=indgeo(1,n,1)                     ! PARALL         
          iglob=indgeo(1,n,1)                     ! PARALL         
!!!          if (iglob.ge.isx1.and.i.le.iex1)then    ! PARALL         
          if (iglob.ge.isx1.and.iglob.le.iex1)then    ! PARALL         
             i = iglob - myid*n1m                    ! PARALL
             j=indgeo(1,n,2)
             k=indgeo(1,n,3)
             forclo(i,j,k)= 0.
          endif
       end do
      end if

      betadx=beta*dx1q*al

      do k=1,n3m
         do j=1,n2m
            do i=1,n1m
               forcl = betadx*forclo(i,j,k)
               amil(j,i)=-forcl
!               acil(j,i)=1.+2.*forcl
!               apil(j,i)=-forcl
               iadd = i+(j-1)*n1m+(k-1)*n1m*n2m
               fil(j,i)=rhs(iadd)
            end do
         end do
!
!         write(*,*)' solq1 ',myid,n1mglob,iex1,isx1,n2m,ey,sy
         call exchng3(amil,amilP,n2m,sy,ey,n1mglob,isx1,iex1,1,comm)   !PARALL
         call exchng3(fil,filP,n2m,sy,ey,n1mglob,isx1,iex1,1,comm)     !PARALL
         do j=1,ey-sy+1
            do i=1,n1mglob
               acilP(j,i) = 1.-2.*amilP(j,i)
               apilP(j,i) = amilP(j,i)
            end do
         end do
!
!         call trvpki(amil,acil,apil,fil,1,n1m,1,n2m,m1,m2)            !PARALL
         call trvpki(amilP,acilP,apilP,filP,1,n1mglob,1,ey-sy+1,m1,ey-sy+2)
!
         call exchng3(fil,filP,n2m,sy,ey,n1mglob,isx1,iex1,-1,comm)    !PARALL
         do j=1,n2m
            do i=1,n1m
               iadd = i+(j-1)*n1m+(k-1)*n1m*n2m
               rhs(iadd)=fil(j,i)
            end do
         end do
      end do
!
      deallocate(amilP)
      deallocate(acilP)
      deallocate(apilP)
      deallocate(filP)
!
      return
      end subroutine solq1i
!
!***********************************************************************
!   this subroutine performs the inversion of the q2 momentum equation
!   by a factored implicit scheme, only the derivatives 11,22,33 of q2
!   are treated implicitly
!   in the first part the rhs is calculated
!        direction x1
!
!     subroutine solq2i(q1)
      subroutine solq2i(comm)
      include'param.f'
      
      real amil(m2,max_m1-1),acil(m2,max_m1-1),apil(m2,max_m1-1),fil(m2,max_m1-1)
      real,    dimension(:,:), allocatable :: amilP,acilP,apilP,filP
      integer comm, ey, sy
!
!  ************ compute dq2** sweeping along the x1 direction
!
      call MPE_DECOMP1D( n2m, numprocs, myid, sy, ey )
      allocate(amilP(ey-sy+2,n1glob))
      allocate(acilP(ey-sy+2,n1glob))
      allocate(apilP(ey-sy+2,n1glob))
      allocate(filP(ey-sy+2,n1glob))
      
      do k=1,n3
        do j=1,n2
          do i=1,n1
            forclo(i,j,k)=1.
          end do
        end do
      end do

      if(infig.ne.-1) then
         do n=1,npunr
!            i=indgeo(2,n,1)                     ! PARALL         
            iglob=indgeo(2,n,1)                     ! PARALL         
!!!            if (iglob.ge.isx1.and.i.le.iex1)then    ! PARALL         
            if (iglob.ge.isx1.and.iglob.le.iex1)then    ! PARALL         
               i = iglob - myid*n1m                    ! PARALL
               j=indgeo(2,n,2)
               k=indgeo(2,n,3)
               forclo(i,j,k)= 0.
            endif
         end do
      end if
      betadx=beta*dx1q*al
      do k=1,n3m
         do j=1,n2m
            do i=1,n1m
               forcl = forclo(i,j,k)*betadx
               amil(j,i)=-forcl
!               acil(j,i)=1.+2.*forcl
!               apil(j,i)=-forcl
               iadd = i+(j-1)*n1m+(k-1)*n1m*n2
               fil(j,i)=rhs(iadd)
            end do
         end do
 !
         call exchng3(amil,amilP,n2m,sy,ey,n1mglob,isx1,iex1,1,comm)   !PARALL
         call exchng3(fil,filP,n2m,sy,ey,n1mglob,isx1,iex1,1,comm)     !PARALL
         do j=1,ey-sy+1
            do i=1,n1mglob
               acilP(j,i) = 1.-2.*amilP(j,i)
               apilP(j,i) = amilP(j,i)
            end do
         end do
!
!         call trvpki(amil,acil,apil,fil,1,n1m,1,n2m,m1,m2)     ! PARALL
         call trvpki(amilP,acilP,apilP,filP,1,n1mglob,1,ey-sy+1,m1,ey-sy+2)
!
         call exchng3(fil,filP,n2m,sy,ey,n1mglob,isx1,iex1,-1,comm)    !PARALL
!
         do j=1,n2m
            do i=1,n1m
               iadd = i+(j-1)*n1m+(k-1)*n1m*n2
               rhs(iadd)=fil(j,i)
            end do
         end do
      end do
!
      deallocate(amilP)
      deallocate(acilP)
      deallocate(apilP)
      deallocate(filP)
!
      return
      end subroutine solq2i
!
!***********************************************************************
!   this subroutine performs the inversion of the q3 momentum equation
!   by a factored implicit scheme, only the derivatives 11,22,33 of q3
!   are treated implicitly
!       direction x1
!
!      subroutine solq3i(q1)        ! PARALL
      subroutine solq3i(comm)
      include'param.f'

      real amil(m2,max_m1-1),acil(m2,max_m1-1),apil(m2,max_m1-1),fil(m2,max_m1-1)
      real,    dimension(:,:), allocatable :: amilP,acilP,apilP,filP
      integer comm, ey, sy
!
!  ************ compute dq3** sweeping along the x1 direction
!               inflow outflow
!
      call MPE_DECOMP1D( n2m, numprocs, myid, sy, ey )
      allocate(amilP(ey-sy+2,n1glob))
      allocate(acilP(ey-sy+2,n1glob))
      allocate(apilP(ey-sy+2,n1glob))
      allocate(filP(ey-sy+2,n1glob))

      do k=1,n3
        do j=1,n2
          do i=1,n1
            forclo(i,j,k)=1.
          end do
        end do
      end do

      if(infig.ne.-1) then
         do n=1,npunz
!            i=indgeo(3,n,1)    ! PARALL
            iglob=indgeo(3,n,1)
!!!            if (iglob.ge.isx1.and.i.le.iex1)then    ! PARALL         
            if (iglob.ge.isx1.and.iglob.le.iex1)then    ! PARALL         
               i = iglob - myid*n1m                    ! PARALL
               j=indgeo(3,n,2)
               k=indgeo(3,n,3)
               forclo(i,j,k)= 0.
            endif
         end do
      end if

      betadx=beta*dx1q*al

      do k=1,n3m
         do j=1,n2m
            do i=1,n1m
               forcl = betadx*forclo(i,j,k)
               amil(j,i)=-forcl
!               acil(j,i)=1.+2.*forcl
!               apil(j,i)=-forcl
               iadd = i+(j-1)*n1m+(k-1)*n1m*n2m
               fil(j,i)=rhs(iadd)
            end do
         end do
!
         call exchng3(amil,amilP,n2m,sy,ey,n1mglob,isx1,iex1,1,comm)   !PARALL
         call exchng3(fil,filP,n2m,sy,ey,n1mglob,isx1,iex1,1,comm)     !PARALL
         do j=1,ey-sy+1
            do i=1,n1mglob
               acilP(j,i) = 1.-2.*amilP(j,i)
               apilP(j,i) = amilP(j,i)
            end do
         end do
!
!         call trvpki(amil,acil,apil,fil,1,n1m,1,n2m,m1,m2)
         call trvpki(amilP,acilP,apilP,filP,1,n1mglob,1,ey-sy+1,m1,ey-sy+2)
!
         call exchng3(fil,filP,n2m,sy,ey,n1mglob,isx1,iex1,-1,comm)    !PARALL
         do j=1,n2m
            do i=1,n1m
               iadd = i+(j-1)*n1m+(k-1)*n1m*n2m
               rhs(iadd)=fil(j,i)
            end do
         end do
      end do
!
      deallocate(amilP)
      deallocate(acilP)
      deallocate(apilP)
      deallocate(filP)
!
      return
      end subroutine solq3i
!
!***********************************************************************
