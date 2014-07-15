c************************************************************************
      subroutine gcurv                                                  
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
      REAL q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)                  
      REAL ru1(m1,m2,m3),ru2(m1,m2,m3),ru3(m1,m2,m3)               
      REAL dq(m1,m2,m3)
      REAL   pr(m1,m2,m3)
      real etime,t(2)
      character(50) namfile
      common/vvtt/ visma,vismi
      pi=2.*asin(1.)                                                    
c                                                                       
c     grid information                                                 
c                                                                       
      call initia(q1,ru1,q2,ru2,q3,ru3,pr,dq)

      call topogr
      call indic                                                        
      call meshes
      call cordin                                                       
      write(6,754)n1,n2,n3                                              
      write(32,754)n1,n2,n3                                             
 754  format(/,5x,'numero di punti griglia: ',' n1= ',i3,' n2= ',i3,
     %     ' n3= ',i3/)                       
      if(iaxsy.eq.0) write(6,856)
      if(iaxsy.eq.0) write(32,856)
 856  format(/,5x,'calcolo con condizione di assialsimmetria   '
     %     ,i3/)                       
      write(6,755) dx1,dx2,dx3,dt,ntst                               
      write(32,755) dx1,dx2,dx3,dt,ntst                              
 755  format(/,2x,'dx1=',e9.3,' dx2=',e9.3,' dx3=',e9.3,'  dt=',e9.3,
     % '  ntst=',i5,/)
      time=0.
      ntii=0                                                            
      if(nson.ne.0) then
         call indson
         do i = 1,nson
            write(6,123) i,(mason(i,j),j=1,3)
            write(32,123) i,(mason(i,j),j=1,3)
         end do
  123 format(2x,'sonda n= ',i2,' : i=',i3,', j=',i3,', k=',i3)
      do i = 1,nson
         write(6,923) i,(coson(i,j),j=1,3)
         write(32,923) i,(coson(i,j),j=1,3)
      end do
      endif
 923  format(2x,'sonda n= ',i2,' : th=',f8.2,', r=',f8.2,', x=',f8.2)
      write(6,*)
      write(32,*)
c     
c     read or create initial fields                                       
c     
      call coetar
      call initia(q1,ru1,q2,ru2,q3,ru3,                        
     %     pr,dq)  
      call phini                                                        
c     
      do 22 l=1,ndv                                                   
         vmax(l)=0.
 22   continue                                                        
c     
c     create the initial conditions
c     
      if(nread.eq.0) then                                               
         write(6,*)' nread=0 ---> cond. iniz. create nel codice'     
         write(32,*)' nread=0 ---> cond. iniz. create nel codice'     
         ntii=0                                                          
         ntime=0                                                         
         time=0.
c     call topogr(time)
         cflm=0.   
         call inqpr(q1,q2,q3,pr)                                        
         write(6,*) ' dopo inqpr '
         call divgck(dmax,dtot,q1,q2,q3)                                 
         if(istat.eq.1) call stacle
         write(6,*)' initial divg dmax,dtot  ',dmax,dtot                 
         write(32,*)' initial divg dmax,dtot  ',dmax,dtot                
c     call outpf(time,q1,q2,q3,pr,dq) 
      else                                                              
         write(6,*)' nread=1 ---> cond. iniz. lette da file'      
         write(32,*)' nread=1 ---> cond. iniz. lette da file'      
         call ini
         call inirea(ntii,time,q1,q2,q3,pr)
         call divgck(dmax,dtot,q1,q2,q3)                                 
         write(6,*)' initial divg dmax,dtot  ',dmax,dtot                 
         write(32,*)' initial divg dmax,dtot  ',dmax,dtot                
         call vort(voamax,vormax,vozmax,q1,q2,q3)                        
         call cfl(cflm,q1,q2,q3)                                         
         call vmaxv(time,q1,q2,q3)                                            
         call outh(time,dt,voamax,vormax,vozmax,q3,pr)
         if(nson.ne.0) call wson(time,q1,q2,q3,pr)
         if(istat.eq.1.and.ireset.eq.1) then
            call stacle
         end if
         if(istat.eq.1.and.ireset.eq.0) then
            call starea
C     C         call stacle
         end if
      endif                                                             
c     
      ntstf=ntii+ntst                                                   
      ntii=ntii+1                                                       
      write(6,711)tprint,ntii,ntstf,tpin                                  
      write(32,711)tprint,ntii,ntstf,tpin                                 
 711  format(3x,'check in cond : tprint =',f8.4,'  ntii =',i8,            
     1     '  ntstf =',i8,2x,'tpin =',f8.4//)                          
      if(idtv.eq.1) then
         write(6,*)'ntime,time,voamax,vormax,vozmax,vmax(1),imax(i),jmaxv(1
     %        ),kmaxv(1),vmax(2),imaxv(2),jmaxv(2),kmaxv(2),vmax(3),imaxv(3),jma
     %        xv(3),kmaxv(3),dt,qmax,iqm,jqm,kqm,densm,densmax,densmin'
         write(32,*)'ntime,time,voamax,vormax,vozmax,vmax(1),imax(i),jmaxv(
     %        1),kmaxv(1),vmax(2),imaxv(2),jmaxv(2),kmaxv(2),vmax(3),imaxv(3),jm
     %        axv(3),kmaxv(3),dt,qmax,iqm,jqm,kqm,densm,densmax,densmin'
         call outh(time,dt,voamax,vormax,vozmax,q3,pr)
         if(nson.ne.0) call wson(time,q1,q2,q3,pr)
      else
      write(6,*)'ntime,time,voamax,vormax,vozmax,vmax(1),imax(i),jmaxv(1
     %    ),kmaxv(1),vmax(2),imaxv(2),jmaxv(2),kmaxv(2),vmax(3),imaxv(3),jma
     %    xv(3),kmaxv(3),cflm,qmax,iqm,jqm,kqm,densm,densmax,densmin'
      write(32,*)'ntime,time,voamax,vormax,vozmax,vmax(1),imax(i),jmaxv(
     %    1),kmaxv(1),vmax(2),imaxv(2),jmaxv(2),kmaxv(2),vmax(3),imaxv(3),jm
     %    axv(3),kmaxv(3),cflm,qmax,iqm,jqm,kqm,densm,densmax,densmin'
         cflm=cflm*dt
c     call outh(time,cflm,voamax,vormax,vozmax,q3,pr)
c     if(nson.ne.0) call wson(time,q1,q2,q3,pr)
      endif
c
      call rotation(time)
      call ringfo(time,time)
      call boucqt(q1,q2,q3,time)
      namfile = 'naca0030-832-0p1.stl'
      pitch1 = 0.
      call readgeo3ddim(nb,namfile,pitch1)
c      call readgeo3d(time,xyzb,norb,nb,dim3_min,dim3_max,dim2_min,
c     $     dim2_max,dim1_min,dim1_max,namfile,pitch1)
      call readgeo3d(time,namfile,pitch1)
c                                                                       
c  ********* starts the time dependent calculation ***                  
c                                                                       
      ntiia = ntii
      do 350 ntime=ntii,ntstf                                           
c
c     the calculation stops if the velocities are diverging for numerica
c     stability conditions (courant number restrictions)                
c
            call cfl(cflm,q1,q2,q3)
            if(idtv.eq.1.and.cflm.gt.1.e-05) then
               if(ntime.gt.1) then
                  if (time.gt.0.7)cflmax = 0.7
                  if (time.gt.1.)cflmax = 1.
                  dt=cflmax/cflm
                  write(6,*)'dt, cflm',dt, cflm !laura
                  if(dt.gt.dtmax) dt=dtmax
              endif
              if(dt.lt..000001) go to 166
            else
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
            if(iles.gt.0.and.ntime.gt.1) then
            tim1=etime(t)
              fvis = 10.
              usren = 1./ren
              del1 = 1./(dx1*dx2*dx3)

              do k=1,n3m
                del3 = g3m(k)
                do j=1,n2m
                  del2 = g2m(j)
                  do i=1,n1m
                    deltsq = (del1*del2*del3)**0.66666666
                    usrece = max(usren ,visct(i,j,k)*usren)
                    fvis = min(fvis,0.167*deltsq/usrece)
                  end do
                end do
              end do
c             write(6,*) ' DT = ',dt,' FVIS = ',fvis
              dt = min(dt,0.8*fvis)
            tim2=etime(t)
C           print *,' in delta=',tim2-tim1,dt,fvis
            end if
            beta=dt/ren*0.5
            timetsch=etime(t)
            call tschem(q1,q2,q3,ru1,ru2,ru3,            
     %                   pr,dq,time)       
            print *,time,'*****Secs for time step ',-timetsch+etime(t)
c           write(6,*) ' a '
            time=time+dt                                                      
            if(ntime.eq.1) go to 306       
CRV         if(nson.ne.0) then
CRV         if(mod(ntime,5).eq.0) 
CRV  %      call wson(time,q1,q2,q3,pr)
CRV         end if
            timeforces=etime(t)
c            if(amod(time,(tpin/5)).lt.dt)then
c     %      call dralical(time,q1,q2,q3,pr) 
               call forces(time,q1,q2,q3,pr) 
c            print *,'*****Secs for forces calcul ',-timeforces+etime(t)
c            end if
            if(nson.ne.0) then
              if(amod(time,(tpin/5)).lt.dt) 
     %        call wson(time,q1,q2,q3,pr)
            end if
            if(istat.eq.1) then
              if(mod(time,tpin/5.).lt.dt) then
                if(ntime.gt.5) then 
                  call stacal(q1,q2,q3,pr,dq)
                  call stawri(time)
                end if
              end if
            end if
            if(amod(time,tpin).lt.dt) go to 306                                
            go to 305                                                         
  306 continue                                                          
            call vmaxv(time,q1,q2,q3) 
            if(vmax(1).gt.1000.and.vmax(2).gt.1000) go to 266                 
            call cfl(cflm,q1,q2,q3)                                           
            call vort(voamax,vormax,vozmax,q1,q2,q3)                        
            call divgck(dmax,dtot,q1,q2,q3)                                 
            if(idtv.eq.1) then
              call outh(time,dt,voamax,vormax,vozmax,q3,pr)
c             if(nson.ne.0) call wson(time,q1,q2,q3,pr)
            else
              cflm=cflm*dt
              call outh(time,cflm,voamax,vormax,vozmax,q3,pr)
CRV           if(nson.ne.0) call wson(time,q1,q2,q3,pr)
            endif
            if(dmax.gt.resid) go to 169                                       
c                                                                       
  305 continue                                                          
c                                                                       
c     write the flow field                                              
c     
      if(amod(time,tprint).lt.dt) then                 
         if(nwrit.eq.1) then    
            call continua(time,q1,q2,q3,pr)
         endif                                                  
        call outpf(time,q1,q2,q3,pr,dq) 
      endif                                                             
      if(time.gt.tmax) go to 167
      if(imovie.eq.1) then
         if(amod(time,tframe).lt.dt) then
            call mkmov(time,q1,q2,q3,pr)
         end if
      end if
            
C           write(6,*) ' G '

  350 continue
      go to 167                                                         
  165 continue                                                          
      write(6,164)                                                      
      write(32,164)                                                     
  164 format(10x,'cfl troppo grande   ')                                
      go to 167                                                         
  166 continue                                                          
      write(6,168)                                                      
      write(32,168) dt                                                  
  168 format(10x,'dt troppo piccolo  DT = ',e14.7)                                
      go to 167                                                         
  266 continue                                                          
      write(6,268)                                                      
      write(32,268)                                                     
  268 format(10x,'velocities diverged')                                 
      go to 167                                                         
  169 continue                                                          
      write(6,178) dmax,imxq,jmxq,kmxq                                  
      write(32,178) dmax,imxq,jmxq,kmxq                                 
  178 format(10x,'too large local residue for mass conservation : '     
     1       ,e12.5,' at ',3(1x,i3) )
  167 continue                                                          
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
      subroutine tschem(q1,q2,q3,ru1,ru2,ru3,
     %                  pr,dq,time)
c 
      include'param.f'
      REAL q1(m1,m2,m3),  q2(m1,m2,m3),  q3(m1,m2,m3)                  
      REAL ru1(m1,m2,m3), ru2(m1,m2,m3), ru3(m1,m2,m3)               
      REAL dq(m1,m2,m3)
      REAL   pr(m1,m2,m3)
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
c       write(6,*)'u0x= ',u0x
c  boundary conditions outflow                                          
        call boucdq(q1,q2,q3,tino) 
        tinfl2=etime(t) 
C       print *,' inflow setup ',tinfl2-tinfl
c
c       LES computation; the turbulent viscosity visct is calculated only
c       at the beginning of the three substeps of the 3rd order R-K scheme.
c
        if(ns.eq.1.and.iles.gt.0) then
        tmod=etime(t) 
          if(iles.eq.1) call dynamrz(time,q1,q2,q3,dq)
          if(iles.eq.2) call smagorz(time,q1,q2,q3,dq)
        tmod2=etime(t) 
C       print *,' computation of les model ',tmod2-tmod
        end if
c
c       computation of the LES terms 2*nu*S_ij
c
        if(iles.gt.0) then
        tles=etime(t) 
          call strain(q1,q2,q3)
        tles1=etime(t) 
          usren = 1./ren
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,usren,visct,st)
!$OMP$ PRIVATE(visk)
            do kc=1,n3m
              do jc=1,n2m
                do ic=1,n1m
               visk=visct(ic,jc,kc)*usren
               st(ic,jc,kc,1)=2.*st(ic,jc,kc,1)*visk
               st(ic,jc,kc,2)=2.*st(ic,jc,kc,2)*visk
               st(ic,jc,kc,3)=2.*st(ic,jc,kc,3)*visk
               st(ic,jc,kc,4)=2.*st(ic,jc,kc,4)*visk
               st(ic,jc,kc,5)=2.*st(ic,jc,kc,5)*visk
               st(ic,jc,kc,6)=2.*st(ic,jc,kc,6)*visk
                end do
              end do
            end do
!$OMP  END PARALLEL DO
        tles2=etime(t) 
C       print *,' in strain and 2nuSij',tles2-tles1,tles1-tles
        end if
c
c       if iles > 0 the viscous LES terms  div(2*nu*S_ij) are computed
c       explicitly inside the nonlinear terms routines
c
        
        timehd=etime(t) 
c        call hdnltot   (q1,q2,q3,dq)    !! NON LINEAR TERMS
         call hdnl1   (q1,q2,q3,dq)    !! NON LINEAR TERMS
        timehd2=etime(t) 
         call hdnl2   (q1,q2,q3)    !!  "    "      "
        timehd3=etime(t) 
         call hdnl3   (q1,q2,q3)        !!  "    "      "
        timehd4=etime(t) 
c       print *,' in hdnltot',timehd4-timehd
C       print *,' in hdnltot',timehd2-timehd,timehd3-timehd2,
C    %                     timehd4-timehd3
c       print *,'maxval dq',maxval(dq(:,:,:)),minval(dq(:,:,:))
c       print *,'maxval qcap',maxval(qcap(:,:,:)),minval(qcap(:,:,:))
c       print *,'maxval dph',maxval(dph(:,:,:)),minval(dph(:,:,:))

c    
        timehd=etime(t) 
        call invtr1  (q1,pr,ru1,dq,aldto)     !! DQ1HAT=Q1HAT-Q1(N)
        timehd2=etime(t) 
        call invtr2  (q2,pr,ru2,q1,aldto)   !! DQ2HAT=Q2HAT-Q2(N)
        timehd3=etime(t) 
        call invtr3  (q3,pr,ru3,q1,aldto)   !! DQ3HAT=Q3HAT-Q3(N)
        timehd4=etime(t) 

C       print *,' in invtr',timehd2-timehd,timehd3-timehd2,
C    %                     timehd4-timehd3
c
        call divg    (q1,q2,q3)          !! DIVG(QHAT)
C       print *,' in divg',etime(t)-timehd4

        call phcalc                      !! PRESSURE (FFT & TRID)
        timeuvp=etime(t) 

        call updvp   (q1,q2,q3) !! SOLENOIDAL VEL FIELD
        timepr=etime(t) 
        call prcalc  (pr)                   !! PRESSURE FIELD
        timepr2=etime(t) 
        call boucqt(q1,q2,q3,tin)
        tinfl2=etime(t) 

C       print *,' in updvp',timepr-timeuvp
C       print *,' in prcalc',timepr2-timepr
C       print *,' in boucqt= ',tinfl2-timepr2
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
      i = n1m/2
      do j = 1,n2
         do k = 1,n3
            write(99,*)x3c(k),x2c(j),x1c(i),qcap(i,j,k)
         enddo
      enddo

      namfile='all-qcap.dat'
      open(99,file=namfile,form='formatted',status='unknown')
      write(99,*) 'VARIABLES = "X","Y","Z","qcap"'
      write(99,*) 'ZONE I = ',n3,', J = ',n2,', K = ',
     $     n1,', F=POINT'
      do i = 1,n1
         do j = 1,n2
            do k = 1,n3
               write(99,*)x3c(k),x2c(j),x1c(i),qcap(i,j,k)
            enddo
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
      i = n1/2
      do j = 1,n2
         do k = 1,n3
            write(99,*)x3c(k),x2c(j),x1c(i),dph(i,j,k)
         enddo
      enddo

      namfile='all-dph.dat'
      open(99,file=namfile,form='formatted',status='unknown')
      write(99,*) 'VARIABLES = "X","Y","Z","qcap"'
      write(99,*) 'ZONE I = ',n3,', J = ',n2,', K = ',
     $     n1,', F=POINT'
      do i = 1,n1
         do j = 1,n2
            do k = 1,n3
               write(99,*)x3c(k),x2c(j),x1c(i),dph(i,j,k)
            enddo
         enddo
      enddo

      close(99)
      end

C********************************************************
