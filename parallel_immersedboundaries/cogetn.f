c************************************************************************
      subroutine gcurv     
      include "mpif.h"                                         
      include'param.f'
c
c     Code for the computation of three-dimensional incompressible flows    
c     in cylindrical polar coordinates.       
                                      
c                                                                       
c     This code solves flow fields bounded in the x3 (axial) and 
c     x2 (radial) directions. All boundaries can be no-slip or free-slip
c     by setting the appropriate indices in the input file.
c     The geometry  (a cylindrical can) is given in the subroutine cordi.
c     The discretization is uniform in the axial (3) and azimuthal (1)
c     directions, while can be non-uniform in the radial direction.
c
c     The equations for the following variables
c                                                                       
c      q1=v(theta)    q2=v(r)*r     q3=v(zeta)                          
c                                                                       
c     are discretized by finite-difference schemes.                    
c     The introduction of the variable q2 is necessary to avoid the
c     problem of the singularity of the equations at the axis of symmetry
c     (r=0).
c
c     All spatial derivatives are discretized by central second-order
c     accurate finite-difference schemes including the non linear terms.  
c
c     The nonlinear terms are treated explicitly, while the viscous terms
c     are computed implicitly. This would lead to the inversion of a
c     large banded matrix, that however is avoided by introducing
c     a factored scheme bringing to the solution of three tridiagonal
c     matrices for each velocity component (subroutine INVTR*).
c                              
c     In time a fractional-step procedure is used in the version of 
c     Nagi Mansour introducing the pressure in the first step.           
c
c     The non-linear terms and the cross derivatives of the viscous terms
c     are discretized by explicit  Adams-Bashfort or 3rd order Runge-Kutta
c     method (A. Wray, personal communication).                      
c
c     The scalar quantity Phi, which projects the provisional velocity field
c     onto  a divergence free field, is solved by a direct method. 
c     For the axial and azimuthal directions modified wave numbers coupled 
c     with trigonometric expansions (FFTs) are used. The equation is then
c     solved by simply inverting a tridiagonal matrix for the radial direction.
c     No explicit boundary conditions are necessary for this Poisson equation. 
     
c                                                                       
c     Other details of the scheme are given in the introduction of the  
c     subroutine TSCHEM
c    
      integer comm1d, ierr 
      real etime,t(2)
      character(50) namfile
      common/vvtt/ visma,vismi
      pi=2.*asin(1.)                                                    
c                                                                       
c     grid information                                                 
c       
      call MPI_CART_CREATE( MPI_COMM_WORLD, 1, numprocs,.true., 
     $     .true., comm1d, ierr )
      call MPI_COMM_RANK( comm1d, myid, ierr )
      call MPI_Cart_shift( comm1d, 0, 1, nbrleft, nbrright, ierr )

      isx1 = 1 + myid*int(n1mglob/numprocs) ! if nx/numprocs is int
      iex1 = isx1 + int(n1mglob/numprocs) - 1 
c      call initia(q1,ru1,q2,ru2,q3,ru3,pr,dq,qcap,dph,isx1,iex1)   !PARALL
      do ii=0,numprocs-1
         if (myid.eq.ii)then
            call topogr
            call MPI_BARRIER( comm1d, ierr )         
         endif
      enddo

      call indic                                                        
      call meshes
      call cordin
      
      if (myid.eq.0)then
         write(6,754)n1,n2,n3                                              
         write(32,754)n1,n2,n3                                             
 754     format(/,5x,'numero di punti griglia: ',' n1= ',i3,' n2= ',i3,
     %        ' n3= ',i3/)                       
         if(iaxsy.eq.0) write(6,856)
         if(iaxsy.eq.0) write(32,856)
 856     format(/,5x,'calcolo con condizione di assialsimmetria',i3/)
         write(6,755) dx1,dx2,dx3,dt,ntst                               
         write(32,755) dx1,dx2,dx3,dt,ntst                              
 755     format(/,2x,'dx1=',e9.3,' dx2=',e9.3,' dx3=',e9.3,'  dt=',e9.3,
     %        '  ntst=',i5,/)
         write(6,*)
         write(32,*)
      endif
      time=0.
      ntii=0   
   
c     if(nson.ne.0) then                         ! PARALL
c         call indson
c         do i = 1,nson
c            write(6,123) i,(mason(i,j),j=1,3)
c            write(32,123) i,(mason(i,j),j=1,3)
c         end do
c 123     format(2x,'sonda n= ',i2,' : i=',i3,', j=',i3,', k=',i3)
c         do i = 1,nson
c            write(6,923) i,(coson(i,j),j=1,3)
c            write(32,923) i,(coson(i,j),j=1,3)
c         end do
c      endif
c 923  format(2x,'sonda n= ',i2,' : th=',f8.2,', r=',f8.2,', x=',f8.2)
c     
c     read or create initial fields                                       
c
      call coetar
c     call initia(q1,ru1,q2,ru2,q3,ru3,pr,dq)
      call initia               !PARALL
      call phini
c     
      do 22 l=1,ndv                                                   
         vmax(l)=0.
 22   continue                                                        
c     
c     create the initial conditions
c     
      if(nread.eq.0) then                                               
         
         if (myid.eq.0)write(6,*)' nread=0 -> c. i. create nel codice'     
         if (myid.eq.0)write(32,*)' nread=0 -> c. i. create nel codice'
         ntii=0                                                         
         ntime=0                                                         
         time=0.
c     call topogr(time)
         cflm=0.   
c     call inqpr(q1,q2,q3,pr)        !PARALL                                
         call inqpr                                        
         call divgck(dmax,dtot,comm1d)                             
c     if(istat.eq.1) call stacle        !COMMENTATO PARALL 
         if (myid.eq.0)write(6,*)' initial divg dmax,dtot  ',dmax,dtot  
         if (myid.eq.0)write(32,*)' initial divg dmax,dtot  ',dmax,dtot 
c     call outpf(time,q1,q2,q3,pr,dq) 
         
      else                                                              
         
         if (myid.eq.0)write(6,*)' nread=1 -> c. i. lette da file'      
         if (myid.eq.0)write(32,*)' nread=1 -> c. i. lette da file'      
         call ini
         call inirea(ntii,time,comm1d)
         call divgck(dmax,dtot,comm1d)                             
         if (myid.eq.0)write(6,*)' initial divg dmax,dtot  ',dmax,dtot                 
         if (myid.eq.0)write(32,*)' initial divg dmax,dtot ',dmax,dtot                
         call vort(voamax,vormax,vozmax,comm1d)
         call cfl(cflm,comm1d)                                       
         call vmaxv(time,comm1d)                                            
         if (myid.eq.0)call outh(time,dt,voamax,vormax,vozmax)
         
c     if(nson.ne.0) call wson(time,q1,q2,q3,pr)       !COMMENTATO PARALL 
c     if(istat.eq.1.and.ireset.eq.1) then
c     call stacle
c     end if
c     if(istat.eq.1.and.ireset.eq.0) then
c     call starea
c     C         call stacle
c     end if

      endif                                                             

      ntstf=ntii+ntst                                                   
      ntii=ntii+1                                                       
      if (myid.eq.0)then
         write(6,711)tprint,ntii,ntstf,tpin                                  
         write(32,711)tprint,ntii,ntstf,tpin   
      endif                              
 711  format(3x,'check in cond : tprint =',f8.4,'  ntii =',i8,            
     $' ntstf =',i8,2x,'tpin =',f8.4//)                          
      if(idtv.eq.1) then
         if (myid.eq.0)then     ! PARALL
      write(6,*)'ntime,time,voamax,vormax,vozmax,vmax(1),imax(i),jmaxv(1
     %),kmaxv(1),vmax(2),imaxv(2),jmaxv(2),kmaxv(2),vmax(3),imaxv(3),jma
     %xv(3),kmaxv(3),dt,qmax,iqm,jqm,kqm,densm,densmax,densmin'
      write(32,*)'ntime,time,voamax,vormax,vozmax,vmax(1),imax(i),jmaxv(
     %1),kmaxv(1),vmax(2),imaxv(2),jmaxv(2),kmaxv(2),vmax(3),imaxv(3),jm
     %axv(3),kmaxv(3),dt,qmax,iqm,jqm,kqm,densm,densmax,densmin'
c         call outh(time,dt,voamax,vormax,vozmax,q3,pr)! PARALL
         call outh(time,dt,voamax,vormax,vozmax) ! PARALL
         endif                     ! PARALL
c         if(nson.ne.0) call wson(time,q1,q2,q3,pr)    ! PARALL
      else
         if (myid.eq.0)then     ! PARALL
      write(6,*)'ntime,time,voamax,vormax,vozmax,vmax(1),imax(i),jmaxv(1
     %),kmaxv(1),vmax(2),imaxv(2),jmaxv(2),kmaxv(2),vmax(3),imaxv(3),jma
     %xv(3),kmaxv(3),cflm,qmax,iqm,jqm,kqm,densm,densmax,densmin'
      write(32,*)'ntime,time,voamax,vormax,vozmax,vmax(1),imax(i),jmaxv(
     %1),kmaxv(1),vmax(2),imaxv(2),jmaxv(2),kmaxv(2),vmax(3),imaxv(3),jm
     %axv(3),kmaxv(3),cflm,qmax,iqm,jqm,kqm,densm,densmax,densmin'
         endif                     ! PARALL
         cflm=cflm*dt
c     call outh(time,cflm,voamax,vormax,vozmax,q3,pr)
c     if(nson.ne.0) call wson(time,q1,q2,q3,pr)
      endif
c
      call rotation(time)
      call ringfo(time,time)
c      call boucqt(q1,q2,q3,time)           ! PARALL
      call boucqt(time)

      call exchng1(q1,n1,n2,n3,comm1d,nbrleft,nbrright)
      call exchng1(q2,n1,n2,n3,comm1d,nbrleft,nbrright)
      call exchng1(q3,n1,n2,n3,comm1d,nbrleft,nbrright)
      call exchng1(pr,n1,n2,n3,comm1d,nbrleft,nbrright)

      namfile = 'surf3n.stl'
      pitch1 = 0.
      do ii=0,numprocs-1
         if (myid.eq.ii)then
            call readgeo3ddim(nb,namfile,pitch1) ! PARALL
            write(6,*)' STL with ',nb,' surface triangles'
            call readgeo3d(time,namfile,pitch1)
!           call readgeo3d(time,xyzb,norb,nb,dim3_min,dim3_max,dim2_min, ! PARALL
!     $     dim2_max,dim1_min,dim1_max,namfile,pitch1)       ! PARALL
            call MPI_BARRIER( comm1d, ierr )         
         endif
      enddo
c  
c  ***************************************************
c  ****** starts the time dependent calculation ******                  
c  ***************************************************
c                                                                       
      ntiia = ntii
      do 350 ntime=ntii,ntstf                                           
c
c     the calculation stops if the velocities are diverging for numerica
c     stability conditions (courant number restrictions)                
c
            call MPI_BARRIER( MPI_COMM_WORLD, ierr )
            t1loop = MPI_WTIME 

c            call cfl(cflm,q1,q2,q3)           ! PARALL
            call cfl(cflm,comm1d)                                     
            if(idtv.eq.1.and.cflm.gt.1.e-05) then
               if(ntime.gt.1) then
                  dt=cflmax/cflm
                  if (myid.eq.0)write(6,*)'dt, cflm',dt, cflm !laura
                  if(dt.gt.dtmax) dt=dtmax
               endif
               if(dt.lt..000001) go to 166
            else
               if (myid.eq.0)write(6,*)'dt, cflm',dt, cflm !laura
               cflm=cflm*dt
               if(cflm.gt.cfllim) go to 165
            endif
c
c     In this code the nonlinear terms are computed explicitly while the
c     viscous terms are computed implicitly. This implies that the stability
c     limit for the time integration is given by the CFL condition.
c     When an LES simulation is carried out, however, the additional LES
c     viscous terms are computed explicitly and this introduces an 
c     additional stability constrain. It turns out that generally this new
c     constrain is less restrictive than the CFL condition (B. Cabot) although
c     intermittently for very limited time periods the latter limitation
c     might yield a time step smaller than the CFL condition. In the
c     following this check is performed and the time step is computed
c     from the most restrictive of the two conditions.
c
c           R.V.  08/10/98
c

c   PARALL COMMENTATA PARTE LARGE EDDIES
c            if(iles.gt.0.and.ntime.gt.1) then
c              tim1=etime(t)
c              fvis = 10.
c              usren = 1./ren
c              del1 = 1./(dx1*dx2*dx3)
c
c              do k=1,n3m
c                del3 = g3m(k)
c                do j=1,n2m
c                  del2 = g2m(j)
c                  do i=1,n1m
c                    deltsq = (del1*del2*del3)**0.66666666
c                    usrece = max(usren ,visct(i,j,k)*usren)
c                    fvis = min(fvis,0.167*deltsq/usrece)
c                  end do
c                end do
c              end do
c              dt = min(dt,0.8*fvis)
c              tim2=etime(t)
c            end if
            beta=dt/ren*0.5
c            call tschem(q1,q2,q3,ru1,ru2,ru3,pr,dq,time)

            call MPI_BARRIER( MPI_COMM_WORLD, ierr )
            t1 = MPI_WTIME
            call tschem(time,comm1d)       
            call MPI_BARRIER( MPI_COMM_WORLD, ierr )
            t2 = MPI_WTIME
            if(myid.eq.0)print *,time,'*****Secs for time step ',t2-t1

            time=time+dt                                                      
            if(ntime.eq.1) go to 306       

c   PARALL COMMENTATA FORZE E SONDE
            if(amod(time,(tpin/5)).lt.dt) 
c     %      call dralical(time,q1,q2,q3,pr) 
c     $           call forces(time,q1,q2,q3,pr) 
     $           call forces(time,comm1d) 
c            if(nson.ne.0) then
c              if(amod(time,(tpin/5)).lt.dt) 
c     %        call wson(time,q1,q2,q3,pr)
c            end if
c   PARALL COMMENTATA STATISTICHE
C            if(istat.eq.1) then
C              if(mod(time,tpin/5.).lt.dt) then
C                if(ntime.gt.5) then 
C                  call stacal(q1,q2,q3,pr,dq)
C                  call stawri(time)
C                end if
C              end if
C            end if
            if(amod(time,tpin).lt.dt) go to 306                                
            go to 305                                                         
  306 continue                                                          
c     call vmaxv(time,q1,q2,q3) ! commentato PARALL  
c     if(vmax(1).gt.1000.and.vmax(2).gt.1000) go to 266  ! commentato PARALL                 
      call vmaxv(time,comm1d) 
      if(vmaxall(1).gt.1000.and.vmaxall(2).gt.1000) go to 266                 
c     call cfl(cflm,q1,q2,q3) !                 commentato PARALL                           
      call cfl(cflm,comm1d)                                     
c     call vort(voamax,vormax,vozmax,q1,q2,q3) !           commentato PARALL 
      call vort(voamax,vormax,vozmax,comm1d)
c     call divgck(dmax,dtot,q1,q2,q3)                                 
      call divgck(dmax,dtot,comm1d) !           commentato PARALL                            
      if(idtv.eq.1) then
         if(myid.eq.0)call outh(time,dt,voamax,vormax,vozmax) ! PARALL
c     call outh(time,dt,voamax,vormax,vozmax,q3,pr)
c             if(nson.ne.0) call wson(time,q1,q2,q3,pr)  NON PARALL
      else
         cflm=cflm*dt
         if(myid.eq.0)call outh(time,dt,voamax,vormax,vozmax) ! PARALL
c     call outh(time,cflm,voamax,vormax,vozmax,q3,pr)
CRV           if(nson.ne.0) call wson(time,q1,q2,q3,pr)  NON PARALL
      endif
      if(dmax.gt.resid) go to 169                                       
c                                                                       
  305 continue                                                          
c                                                                       
c     write the flow field                                              
c     
      if(amod(time,tprint).lt.dt) then                 
         if(nwrit.eq.1) then 
c     call continua(time,q1,q2,q3,pr) ! PARALL
            call continua(time,comm1d)
         endif                                                  
         call outpf(time,comm1d) 
      endif                                                             
      if(time.gt.tmax) go to 167
C     if(imovie.eq.1) then               COMMENTATO PARALL
C         if(amod(time,tframe).lt.dt) then
C            call mkmov(time,q1,q2,q3,pr)
C         end if
C      end if

      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      t2loop = MPI_WTIME 
      if(myid.eq.0)print *,time,'*****Secs for loop ',t2loop-t1loop
      
  350 continue

      go to 167                                                         
  165 continue                                                          
      if (myid.eq.0)then 
         write(6,164)                                                      
         write(32,164)
      endif                                                     
  164 format(10x,'cfl troppo grande   ')                                
      go to 167                                                         
  166 continue                                                          
      if (myid.eq.0)then 
         write(6,168)                                                      
         write(32,168) dt   
      endif                                               
  168 format(10x,'dt troppo piccolo  DT = ',e14.7)                                
      go to 167                                                         
  266 continue   
      if (myid.eq.0)then 
         write(6,268)                                                      
         write(32,268)         
      endif                                            
  268 format(10x,'velocities diverged')                                 
      go to 167                                                         
  169 continue
      if (myid.eq.0)then                                                           
         write(6,178) dmax,imxq,jmxq,kmxq                                  
         write(32,178) dmax,imxq,jmxq,kmxq
      endif                                 
  178 format(10x,'too large local residue for mass conservation : '     
     1       ,e12.5,' at ',3(1x,i3) )
  167 continue      
      call continua(time,comm1d)                                                    
      return                                                            
      end                                                               
c
c
************************************************************************
c
c           SUBROUTINE  TSCHEM
c
c   This subroutine manages the whole integration scheme.
c   The following equations are solved:          
c   
c    ~~     n
c   Q  -  Q                n         n       n-1   alp       2  ~~   n 
c  --------- = -alp*grad (P ) + gam*H + rho*H   + ----- nabla ( Q + Q )
c    d t                                          2 Re
c
c          i                           i               i
c   where H  are the nonlinear terms, P  the pressure Q  the velocities
c       ~~
c   and Q  the provisional non solenoidal velocity field.
c   The superscripts (~~, n, n-1) indicate the time step level. 
c                        n
c   The nonlinear terms H  are computed in the routines HDNL*, while
c   in the routines INVTR* are computed the remaining terms, updated
c   the non linear terms and inverted the equation to find the provisional
c   field at the new time step.
c       ~~
c   The Q  velocity field is projected onto a solenoidal field by a 
c   scalar Phi computed through the equation
c
c                         2            1          ~~
c                    nabla (Phi ) =  ------ div ( Q  )
c                                    alp dt
c
c   The right hand side of this equation is computed in the routine
c   DIVG, while the equation is solved in PHCALC.
c
c   In the routine UPDVP the solenoidal velocity field at the new time
c   step is then computed through
c
c                n+1  ~~
c               Q   = Q  - alt*dt grad (Phi)
c
c   Finally in the routine PRCALC is updated the pressure field
c
c                n+1   n        alp dt      2
c               P   = P + Phi - ------ nabla (Phi)
c                                2 Re
c
c   When the scalar field is computed (density, concentration,
c   temperature) the routines HDNLRO and INVTRRO are used. The same
c   strategy at the velocity field is used, except that the scalar
c   field does not need any correction.
c
c   All variables are located on a staggered grid with the velocities
c   on the faces of the computational cell and all the scalars at the
c   centre. This is important when terms belonging to different equations
c   are avaluated.
c
c   Further details of the scheme can be found in the paper
c   "A finite-difference scheme for three-dimensional incompressible
c    flows in cylindrical coordinates" by R. Verzicco and P. Orlandi
c    J. of Comp. Phys. 1996.
c
c
c      subroutine tschem(q1,q2,q3,ru1,ru2,ru3,pr,dq,time)    ! PARALL
      subroutine tschem(time,comm)
c 
      include'param.f'
      include "mpif.h"

      integer comm,ierr
      real etime,t(2)
c
c   TIME INTEGRATION : implicit viscous, 3rd order RK (Adams Bashfort)  
c

      tin=time
      do 2000 ns=1,nsst
        if(ntime.eq.1) then
          aldto = alm(1)*dt
        else
          aldto = al*dt
        end if
        al=alm(ns)
        ga=gam(ns)                                                        
        ro=rom(ns)                                                        
        tino=tin
        tin=tin+al*dt
c
c   axial velocity time dependency
c
        tinfl=etime(t) 
        call ringfo(tin,tino)
        call rotation (tin)                                       
        call boucdq(tino,comm) 

c
c       LES computation; the turbulent viscosity visct is calculated only
c       at the beginning of the three substeps of the 3rd order R-K scheme.
c
c       if(ns.eq.1.and.iles.gt.0) then                       ! COMMENTATO PARALL
c          tmod=etime(t) 
c          if(iles.eq.1) call dynamrz(time,q1,q2,q3,dq)
c          if(iles.eq.2) call smagorz(time,q1,q2,q3,dq)
c          tmod2=etime(t) 
c       end if
c
c       computation of the LES terms 2*nu*S_ij
c
C        if(iles.gt.0) then
C           tles=etime(t) 
C           call strain(q1,q2,q3)
C           tles1=etime(t) 
C           usren = 1./ren
C           do kc=1,n3m
C              do jc=1,n2m
C                 do ic=1,n1m
C               visk=visct(ic,jc,kc)*usren
C               st(ic,jc,kc,1)=2.*st(ic,jc,kc,1)*visk
C               st(ic,jc,kc,2)=2.*st(ic,jc,kc,2)*visk
C               st(ic,jc,kc,3)=2.*st(ic,jc,kc,3)*visk
C               st(ic,jc,kc,4)=2.*st(ic,jc,kc,4)*visk
C               st(ic,jc,kc,5)=2.*st(ic,jc,kc,5)*visk
C               st(ic,jc,kc,6)=2.*st(ic,jc,kc,6)*visk
C                 end do
C              end do
C          end do
C          tles2=etime(t) 
C        end if
c
c       if iles > 0 the viscous LES terms  div(2*nu*S_ij) are computed
c       explicitly inside the nonlinear terms routines
c

        call hdnl1                 !! NON LINEAR TERMS
        call hdnl2                 !! NON LINEAR TERMS
        call hdnl3                 !! NON LINEAR TERMS
c        call outpf1(tin,commd)

        call invtr1  (aldto,tin,comm)     !! DQ1HAT=Q1HAT-Q1(N)
        call invtr2  (aldto,comm)   !! DQ2HAT=Q2HAT-Q2(N)
        call invtr3  (aldto,comm)   !! DQ3HAT=Q3HAT-Q3(N)
c        call outpf2(tin,commd)
        
        call exchng1(q1,n1,n2,n3,comm,nbrleft,nbrright)
        call exchng1(q2,n1,n2,n3,comm,nbrleft,nbrright)
        call exchng1(q3,n1,n2,n3,comm,nbrleft,nbrright)
       
        call divg             !! DIVG(QHAT)
c        call outpf3(tin,commd)

        call phcalc(comm)        !! PRESSURE (FFT & TRID)
        call exchng1(dph,n1,n2,n3,comm,nbrleft,nbrright)
c        call outpf4(tin,commd)

        call updvp               !! SOLENOIDAL VEL FIELD
        call exchng1(q1,n1,n2,n3,comm,nbrleft,nbrright)
        call exchng1(q2,n1,n2,n3,comm,nbrleft,nbrright)
        call exchng1(q3,n1,n2,n3,comm,nbrleft,nbrright)
c        if (myid.eq.0)write(*,*)q1(9,9,9),q2(9,9,9),q3(9,9,9)
c        call outpf5(tin,commd)

        call prcalc                          !! PRESSURE FIELD
        call exchng1(pr,n1,n2,n3,comm,nbrleft,nbrright)
c        if (myid.eq.0)write(*,*)pr(9,9,9)
c        call outpf6(tin,commd)

        call boucqt(tin)

 

 2000 end do
      return                                                            
      end                                                               


C********************************************************

      subroutine out_qcap
      include'param.f'
      character*20 namfile

      namfile='qcap.dat'
      open(99,file=namfile,form='formatted',status='unknown')
      do 4 k=1,n3   
        do 4 j=1,n2
          do 4 i=1,n1 
             write(99,*)qcap(i,j,k)
    4 continue
      close(99)

      namfile='i-qcap.dat'
      open(99,file=namfile,form='formatted',status='unknown')
      write(99,*) 'VARIABLES = "X","Y","Z","qcap"'
      write(99,*) 'ZONE I = ',n3,', J = ',n2,', K = ',
     $     1,', F=POINT'
      i = 3
      do j = 1,n2
         do k = 1,n3
            write(99,*)x3c(k),x2c(j),x1c(i),qcap(i,j,k)
         enddo
      enddo

      close(99)
      end

C********************************************************

      subroutine out_dph
      include'param.f'
      character*20 namfile

      namfile='dph.dat'
      open(99,file=namfile,form='formatted',status='unknown')
      do 4 k=1,n3   
        do 4 j=1,n2
          do 4 i=1,n1 
             write(99,*)dph(i,j,k)
    4 continue
      close(99)

      namfile='i-dph.dat'
      open(99,file=namfile,form='formatted',status='unknown')
      write(99,*) 'VARIABLES = "X","Y","Z","qcap"'
      write(99,*) 'ZONE I = ',n3,', J = ',n2,', K = ',
     $     1,', F=POINT'
      i = 3
      do j = 1,n2
         do k = 1,n3
            write(99,*)x3c(k),x2c(j),x1c(i),dph(i,j,k)
         enddo
      enddo

      close(99)
      end

C********************************************************
