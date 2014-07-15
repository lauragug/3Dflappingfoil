c     this code is made for simulating three-dimensional flows in polar
c     ( cilindrical ) coordinates.
c     boundary condition are slip-walls ( =0 ) or non-slip wall ( =1)       
c     in the radial and the vertical direction.                  
c                                                                       
c     this code allows , by changing suitable indices , to introduce
c     the buoiancy term ( ibuo=1 ) the coriolis terms ( icorio=1 )
c     and to consider the transport of a scalar ( isca=1 )
c     if the scalar is passive, the term for the centrifug force can
c     be included in the pressure so the index ifugo can be set =0,
c     otherwise it has to be set =1 (if the scalar is not a passiv one)
c
c      navier-stokes equations are solved by a fractional step method   
c     ( Kim and Moin ) with the pressure in the first step.              
c            
c     The time advancement of the solution is obtained by a
c     Runge-Kutta 3rd order low storage scheme (Wray) or a 2nd   
c     order Adams-Bashfort scheme.                                      
c
c     The Poisson  equation for the pressure is solved directly                
c     introducing FFT in the azimutal and vertical direction.
c     Because of the presence of the walls in the vertical direction
c     the cosFFT in that direction is required
c                                                                       
c                 Roberto Verzicco and Paolo Orlandi                       
c                 dipartimento di meccanica ed aeronautica              
c                 universita' la sapienza di roma                 
c                                                                       
c                                                                       
c                                                                       
c     All variables are calculated in a staggered grid:
c   
c        Instead of velocities, "flux" variables are introduced to avoid
c        the singularity for r=0
c        q1=vtheta, q2= r*vr, q3= vz       
c
c        dq1,dq2,dq3 : velocity correction                              
c
c        qcap :divergence of the  non free divergent velocity field    
c
c        non linear terms:
c
c        ru1, ru2, ru3, ruro : old step      
c        h1,h2,h3, hro : new step           
c        dens : density
c        pr : pressure                                                  
c        dph : pressure correction                                      
c       pressure solver is in the package press.f                      
c       non linear terms are calculated in hdnl routines                
c       the invertions of momentum equations is performed in invtr      
c       routines                                             
c                                                                       
      program papero
      include'param.f'
      character*4 dummy
      open(unit=15,file='coge.in',status='old')
        read(15,301) dummy
        read(15,*) n1,n2,n3,nsst,nwrit,nread
        read(15,301) dummy
        read(15,*) n1p,n2p,n3p
        read(15,301) dummy
        read(15,*) alx3old,istr3,str3,rmed31,etdp3,strb3
        read(15,301) dummy
        read(15,*) alx2,istr,strr,rmed,etdp,strb
        read(15,301) dummy
        read(15,*) alx1
        read(15,301) dummy
        read(15,*) inslws,inslwn,inslwe,inslww,ifugo,icorio,ibuo,isca
        read(15,301) dummy
        read(15,*) ren,ri,ros,pec
        read(15,301) dummy
        read(15,*) ntst,tprint,tpin,tmax,ireset,ikick
        read(15,301) dummy
        read(15,*) idtv,dtmax,dt,cflmax,cfllim,resid
        read(15,301) dummy
        read(15,*) nini,nfin,nstri
        read(15,301) dummy
        read(15,*)radinf,tau1,tau2,all,alt,cou
        read(15,301) dummy
        read(15,*)infig,r1,r2,z1,z2,hn
        read(15,301) dummy
        read(15,*)imovie,tframe
        read(15,301) dummy
        read(15,*)iles,nson
        read(15,301) dummy
        read(15,*)istat,imed,jmed
        read(15,301) dummy
        read(15,*)igext
        read(15,301) dummy
        read(15,*) uzero,amplit,thetamax
301     format(a4)                
      close(15)

      uzero_max = uzero
      uzero = 0.
      tau3 = tau1+tau2
      if(nson.ne.0) then
      open(unit=15,file='sonda.in',status='old')
        do i = 1,nson
          read(15,*) (coson(i,j),j=1,3)
        end do
      close(15)
      end if
c                                                                       
c
      call openfi
c
      pi=2.*asin(1.)                                                    
      tfini=dt*ntst                                                     
      tfin=tfini                                                        
      n1m=n1-1                                                          
      if(n1.eq.1) n1m = 1
      n2m=n2-1                                                          
      n3m=n3-1                                                          
      n3mh=n3m/2+1                                                      
      if(n1m.eq.1) then
        iaxsy=0
      else
        iaxsy=1
      endif
c                                                                       
c                                                                       
c     assign coefficients for time marching schemes                     
c
      if(nsst.gt.1) then   
        gam(1)=8./15.                                                     
        gam(2)=5./12.                                                      
        gam(3)=3./4.                                                       
        rom(1)=0.                                                         
        rom(2)=-17./60.                                                    
        rom(3)=-5./12.                                                     
        write(6,100) (gam(n),n=1,nsst),(rom(n),n=1,nsst)                    
  100   format(/,5x,'The time scheme is a III order Runge-Kutta'
     %  ,4x,'gam= ',3f8.3,4x,'ro= ',3f8.3)                                
      else                                                              
        gam(1)=1.5                                                         
        gam(2)=0.                                                         
        gam(3)=0.                                                         
        rom(1)=-0.5                                                       
        rom(2)=0.                                                         
        rom(3)=0.                                                         
        write(6,110) gam(1),rom(1)                                          
  110   format(/,5x,'The time scheme is the Adams-Bashfort',4x,
     %   'gam= ',f8.3,4x,'ro= ',f8.3)                                  
      endif                                                             
      do 10 ns=1,nsst
        alm(ns)=(gam(ns)+rom(ns))
   10 continue
c      write(6,112)alx3,rext,h1                                              
c      write(32,112)alx3,rext,h1                                             
  112 format(//,20x,'P A P E R O ',//,10x,
     % 'Perturbazioni Azimutali in PEntola ROtante',//,5x,
     % 'Cilindro:  H0=',f4.2,' R0=',f4.2,' H1=',e10.3,' R1=1.')
      write(6,201) epsil,lamb,rper
      write(32,201) epsil,lamb,rper
  201 format(/,5x,'perturb (seno): eps= ',e10.4,' lamb= ',i3,/,
     % ' perturb (random) rper= ',e10.4)
      write(6,202) ren,ri,ros,pec,delro
      write(32,202) ren,ri,ros,pec,delro
  202 format(/,5x,'Parametri: ',' Re=',e9.3,' Ri= ',e9.3,' Ro= ',
     % e9.3,' Pe=',e9.3,' delta_rho= ',e9.3)
      write(6,203) nlev,dism,nn
      write(32,203) nlev,dism,nn
  203 format(/,5x,'parametri di smoothing: nlev= ',i3,' dism= ',
     % e10.4,' nn= ',i3)
      if(idtv.eq.1) then
         write(6,204) cflmax
         write(32,204) cflmax
  204 format(/,5x,'calcolo effettuato con dt variabile e cfl= ',
     % e10.4,/ )            
      else 
         write(6,205) dtmax,cfllim
         write(32,205) dtmax,cfllim
  205 format(/,5x,'calcolo effettuato con dt= ',e10.4,' e cflmax=
     %  ',e10.4,/ )            
      endif
c
c     starts the solution of the problem
c
      call gcurv                                                        
c
      stop                                                              
      end                                                               
c                                                                       
c***********************************************************************
      subroutine openfi
      include'param.f'
      character*14 filth,filte,filba,filen,filso
      character*14 nome
      character*2 nums
      character*5 fili
      open(46,file='nfcoge')
      read(46,'(a)')filth
      read(46,'(a)')filte
      read(46,'(a)')filba
      read(46,'(a)')filen
      read(46,'(a)')filso
      open(32,file=filth)
      open(34,file=filba)
      open(35,file='filbaper.out')
      open(39,file=filen)
      open(40,file='errcal.out')
      open(65,file='forces.out',form='formatted')
      open(66,file='coefor.out',form='formatted')
      open(67,file='move.out',form='formatted')
      open(88,file='drag.dat',form='formatted')
      open(95,file='hsep.dat',form='formatted')
      open(97,file=filte)
c     if(nson.ne.0) then
c      nul=50
c      write(fili,200) filso
c200   format(a5)   
c      do i=1,nson
c      write(nums,100) i
c100   format(i2.2)
c      nome=fili//nums//'.out'
c      open(nul,file=nome,status='unknown')
c      rewind(nul)
c      nul=nul+i
c      end do
c     endif
      rewind(32)
      rewind(34)
      rewind(35)
      rewind(39)
      rewind(97)
      close(46)
      return
      end   
c
c***********************************************************************
      subroutine cfl(cflm,q1,q2,q3)                                     
      include'param.f'
      REAL q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)                  
c                                                                       
      cflm=0.                                                          
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,ipv)
!$OMP$ SHARED(q1,q2,q3,dx1,udx2m,udx3m,usdtal)
!$OMP$ PRIVATE(qcf,kp,jp,ip)
!$OMP$ REDUCTION(max:cflm)

      do 7 k=1,n3m 
        kp=k+1   
          do 7 j=1,n2m 
            jp=j+1 
              do 7 i=1,n1m   
                ip=ipv(i)  
                qcf=((
     %                abs((q1(i,j,k)+q1(ip,j,k))*0.5*dx1) 
     %               +abs((q2(i,j,k)+q2(i,jp,k))*0.5*udx2m(j)))
     %               +abs((q3(i,j,k)+q3(i,j,kp))*0.5*udx3m(k)))    
                cflm = max(cflm,qcf)
  7   continue
!$OMP  END PARALLEL DO
      return  
      end                                                               
c                                                                       
c***********************************************************************
      subroutine divgck(qmax,qtot,q1,q2,q3) ! residuo !
      include'param.f'
      REAL q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)                  
c                                                                       
      dvol= 1./(dx1*dx2*dx3)
      qtot=0.                                                          
      qmax=0.                                                          
      do 11 kc=1,n3m                                                  
        usrnu3=dx3/g3m(kc)
        kp=kc+1                                                        
        do 11 jc=1,n2m                                                    
          jp=jc+1                                                           
          udx2 = dx2/g2m(jc)
          do 11 ic=1,n1m                                                   
             ip=ipv(ic)                                                        
             dqcap=(((q1(ip,jc,kc)-q1(ic,jc,kc))*dx1*iaxsy
     %             +(q2(ic,jp,kc)-q2(ic,jc,kc))*udx2)
     %             +(q3(ic,jc,kp)-q3(ic,jc,kc))*usrnu3)
             if(abs(dqcap).gt.qmax) then
                imxq=ic
                jmxq=jc
                kmxq=kc
             endif                             
             qmax = max(abs(dqcap),qmax)
             qtot=qtot+dqcap*g2m(jc)
   11 continue      
      qqmax=qmax
      qtot=qtot*dvol
      return     
      end         
c                  
c***********************************************************************
      subroutine meshes
      include'param.f'
      dx1=(x1c(n1)-x1c(1))/float(n1m)                                              
      dx2=1./float(n2m)                                                 
      dx3=1./float(n3m)                                               
      write(6,*)'dx1, dx2, dx3  ',dx1,dx2,dx3                           
      dx1=1./dx1                                                        
      dx2=1./dx2                                                        
      dx3=1./dx3                                                        
      dx1q=dx1*dx1                                                      
      dx2q=dx2*dx2                                                      
      dx3q=dx3*dx3                                                      
      return                                                            
      end                                                               
c                                                                       
c***********************************************************************
      subroutine vmaxv(time,q1,q2,q3) 
      include'param.f'
      REAL q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)
c
      write(*,*)' VMAX-IN'
        vmax(1)=0.  
        vmax(2)=0.  
        vmax(3)=0.  
        do 312 k=1,n3m   
          do 312 j=2,n2m 
            do 312 i=1,n1m 
               vca1=abs(q1(i,j,k))
               vca2=abs(q2(i,j,k))
               vca3=abs(q3(i,j,k)) 
              if(vca1.ge.vmax(1)) then 
                vmax(1)=vca1
                imaxv(1)=i 
                jmaxv(1)=j 
                kmaxv(1)=k 
              end if
              if(vca2.ge.vmax(2)) then 
                vmax(2)=vca2 
                imaxv(2)=i 
                jmaxv(2)=j 
                kmaxv(2)=k 
              end if
              if(vca3.ge.vmax(3)) then 
                vmax(3)=vca3
                imaxv(3)=i 
                jmaxv(3)=j 
                kmaxv(3)=k 
              end if
  312   continue 
        
        write(*,*)'VMAX',vmax(1),imaxv(1),jmaxv(1),kmaxv(1),
     $       vmax(2),imaxv(2),jmaxv(2),kmaxv(2),vmax(3),
     $       imaxv(3),jmaxv(3),kmaxv(3)

      return   
      end     
c                                                                       
c                                                                       
c***********************************************************************
      subroutine continua(time,q1,q2,q3,pr)             
      include'param.f'
      character*19 filcnw
      REAL q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)                  
      REAL   pr(m1,m2,m3)

      filcnw = 'continua.dat'
      open(13,file=filcnw,form='unformatted')
        rewind(13)                                                      
        write(13)n1,n2,n3
        write(13) epsil,nusf,re,time 
c        write(13) q1,q2,q3,pr,csrz_t
        do i=1,n1
           do j=1,n2
              do k=1,n3
                 write(13)q1(i,j,k),q2(i,j,k),q3(i,j,k),
     $                pr(i,j,k),csrz_t(i,j,k)
              enddo
           enddo
        enddo
c        write(13)(((q1(i,j,k),i=1,n1),j=1,n2),k=1,n3)           
c        write(13)(((q2(i,j,k),i=1,n1),j=1,n2),k=1,n3)        
c        write(13)(((q3(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c        write(13)(((pr(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c        write(13)(((csrz_t(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c        write(13) dq1x2o,dq2x2o,dq3x2o
        do i=1,n1
           do j=1,n2
              write(13)dq1x2o(i,j),dq2x2o(i,j),dq3x2o(i,j)
           enddo
        enddo
c        write(13)((dq1x2o(i,j),i=1,n1),j=1,n2)
c        write(13)((dq2x2o(i,j),i=1,n1),j=1,n2)
c        write(13)((dq3x2o(i,j),i=1,n1),j=1,n2)
c        write(13) qb1s,qb2s,qb3s,qb1n,qb2n,qb3n,qb2dn,qb2up,
c     %            qb3dn,qb3up
        do i=1,n1
           do j=1,n2
              write(13)qb1s(i,j),qb2s(i,j),qb3s(i,j),
     $             qb1n(i,j),qb2n(i,j),qb3n(i,j)
           enddo
        enddo
        do i=1,n1
           do k=1,n3
              write(13)qb2dn(i,k),qb2up(i,k),qb3dn(i,k),qb3up(i,k)
           enddo
        enddo
c        write(13)((qb1s(i,j),i=1,n1),j=1,n2)
c        write(13)((qb2s(i,j),i=1,n1),j=1,n2)
c        write(13)((qb3s(i,j),i=1,n1),j=1,n2)
c        write(13)((qb1n(i,j),i=1,n1),j=1,n2)
c        write(13)((qb2n(i,j),i=1,n1),j=1,n2)
c        write(13)((qb3n(i,j),i=1,n1),j=1,n2)
c        write(13)((qb2dn(i,k),i=1,n1),k=1,n3)
c        write(13)((qb2up(i,k),i=1,n1),k=1,n3)
c        write(13)((qb3dn(i,k),i=1,n1),k=1,n3)
c        write(13)((qb3up(i,k),i=1,n1),k=1,n3)
      close(13)
      return                                                            
      end                                                               
c                                                                       
c***********************************************************************
      subroutine inirea(ntii,time,q1,q2,q3,pr)
      include'param.f'
      character*19 filcnw
      REAL q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)                  
      REAL   pr(m1,m2,m3)
      filcnw = 'continua.dat'
      open(13,file=filcnw,form='unformatted')
      read(13)n1,n2,n3
      write(6,*) n1,n2,n3
      read(13) epsil,nusf,re,time
      write(6,*) epsil,nusf,re,time 
c        read(13) q1,q2,q3,pr,csrz_t
c      read(13)(((q1(i,j,k),i=1,n1),j=1,n2),k=1,n3)           
c      read(13)(((q2(i,j,k),i=1,n1),j=1,n2),k=1,n3)        
c      read(13)(((q3(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c      read(13)(((pr(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c      read(13)(((csrz_t(i,j,k),i=1,n1),j=1,n2),k=1,n3)
      do i=1,n1
         do j=1,n2
            do k=1,n3
               read(13)q1(i,j,k),q2(i,j,k),q3(i,j,k),
     $              pr(i,j,k),csrz_t(i,j,k)
            enddo
         enddo
      enddo
      write(6,*) ' AAA '
c     read(13) dq1x2o,dq2x2o,dq3x2o
c      read(13)((dq1x2o(i,j),i=1,n1),j=1,n2)
c      read(13)((dq2x2o(i,j),i=1,n1),j=1,n2)
c      read(13)((dq3x2o(i,j),i=1,n1),j=1,n2)
c     read(13) qb1s,qb2s,qb3s,qb1n,qb2n,qb3n,qb2dn,qb2up,
c     %           qb3dn,qb3up
      do i=1,n1
         do j=1,n2
            read(13)dq1x2o(i,j),dq2x2o(i,j),dq3x2o(i,j)
         enddo
      enddo
      write(6,*) ' BBB '
c      read(13)((qb1s(i,j),i=1,n1),j=1,n2)
c      read(13)((qb2s(i,j),i=1,n1),j=1,n2)
c      read(13)((qb3s(i,j),i=1,n1),j=1,n2)
c      read(13)((qb1n(i,j),i=1,n1),j=1,n2)
c      read(13)((qb2n(i,j),i=1,n1),j=1,n2)
c      read(13)((qb3n(i,j),i=1,n1),j=1,n2)
c      read(13)((qb2dn(i,k),i=1,n1),k=1,n3)
c      read(13)((qb2up(i,k),i=1,n1),k=1,n3)
c      read(13)((qb3dn(i,k),i=1,n1),k=1,n3)
c      read(13)((qb3up(i,k),i=1,n1),k=1,n3)
      do i=1,n1
         do j=1,n2
            read(13)qb1s(i,j),qb2s(i,j),qb3s(i,j),
     $           qb1n(i,j),qb2n(i,j),qb3n(i,j)
         enddo
      enddo
      do i=1,n1
         do k=1,n3
            read(13)qb2dn(i,k),qb2up(i,k),qb3dn(i,k),qb3up(i,k)
         enddo
      enddo
      write(6,*) ' CCC '
      close(13)
      if (ireset.eq.1) then                                             
      ihist=0                                                           
      ntii=0                                                            
      time=0.                                                          
      endif                                                             
      return                                                            
      end                                                               
c                                                                       
c                                                                       
c***********************************************************************
      subroutine initia(q1,ru1,q2,ru2,q3,ru3,pr,dq)
      include'param.f'
      REAL   pr(m1,m2,m3)
      REAL q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)
      REAL dq(m1,m2,m3)
      REAL ru1(m1,m2,m3),ru2(m1,m2,m3),ru3(m1,m2,m3)
!$OMP  PARALLEL
!$OMP$ SHARED(n3,n2,n1)
!$OMP$ SHARED(q1,q2,q3,pr,dq,ru1,ru2,ru3,qcap,dph,st)
!$OMP$ SHARED(forclo,visct,csrz_t)
!$OMP$ PRIVATE(k,j,i)
!$OMP DO

      do 4 k=1,n3   
        do 4 j=1,n2
          do 4 i=1,n1  
            pr(i,j,k)=0.  
            q1(i,j,k)=0.   
            q2(i,j,k)=0.  
            q3(i,j,k)=0. 
c            q2(i,j,k)=uzero*sin(thetamax*pi/180)
c            q3(i,j,k)=uzero*cos(thetamax*pi/180) 
            dq(i,j,k)=0. 
            ru1(i,j,k)=0.
            ru2(i,j,k)=0.
            ru3(i,j,k)=0.
            qcap(i,j,k)=0.
            dph(i,j,k)=0.
            st(i,j,k,1)=0.
            st(i,j,k,2)=0.
            st(i,j,k,3)=0.
            st(i,j,k,4)=0.
            st(i,j,k,5)=0.
            st(i,j,k,6)=0.
            forclo(i,j,k)=0.
            visct(i,j,k)=0.
            csrz_t(i,j,k)=0.
    4 continue
!$OMP END DO
!$OMP  END PARALLEL

c
      do i=1,n1                                                     
        do k=1,n3
        qb2dn(i,k)=0.
        qb2up(i,k)=0.
        qb3dn(i,k)=0.
        qb3up(i,k)=0.
        dqb2dn(i,k)=0.
        dqb2up(i,k)=0.
        dqb3dn(i,k)=0.
        dqb3up(i,k)=0.
        enddo
        do j=1,n2                                                
          uinfth(i,j)=0.
          qb1s(i,j)=0.
c          qb2s(i,j) = uzero*sin(thetamax*pi/180)
c          qb3s(i,j) = uzero*cos(thetamax*pi/180) 
          qb2s(i,j) = 0.
          qb3s(i,j) = 0.
          qb1n(i,j) = 0.
c          qb2n(i,j) = uzero*sin(thetamax*pi/180)
c          qb3n(i,j) = uzero*cos(thetamax*pi/180)
          qb2n(i,j) = 0.
          qb3n(i,j) = 0.
          dqb1s(i,j)=0.
          dqb2s(i,j)=0.
          dqb3s(i,j)=0.
          dqb1n(i,j)=0.
          dqb2n(i,j)=0.
          dqb3n(i,j)=0.
          dq1x2o(i,j)=0.
          dq2x2o(i,j)=0.
          dq3x2o(i,j)=0.
        end do
      end do
      do n=1,npunr
         q2bo(n)=0.
      end do
      do n=1,npunz
         q3bo(n)=0.
      end do
      return 
      end   
c          
c 
c***********************************************************************
      subroutine outh(time,cflm,voamax,vormax,vozmax,q3,pr)
      include'param.f'
      common/profi/  yupm(m3),yupc(m3),ylom(m3),yloc(m3),
     %               ymec(m3),ymem(m3),zz_min,zz_max,
     %               zm_min,zm_max,pitch
      common/vvtt/ visma,vismi
      REAL q3(m1,m2,m3)
      REAL   pr(m1,m2,m3)
      character*20 namprl,nampru
      character*6 ipfi
c                                                                       
      dpi=2.*2.*asin(1.)                                                    
      write(6,159)ntime,time,voamax,vormax,vozmax, 
     % vmax(1),imaxv(1),jmaxv(1),kmaxv(1),
     % vmax(2),imaxv(2),jmaxv(2),kmaxv(2),
     % vmax(3),imaxv(3),jmaxv(3),kmaxv(3),
     % cflm,qqmax,imxq,jmxq,kmxq,densm,denmax,denmin
      write(32,159)ntime,time,voamax,vormax,vozmax, 
     % vmax(1),imaxv(1),jmaxv(1),kmaxv(1),
     % vmax(2),imaxv(2),jmaxv(2),kmaxv(2),
     % vmax(3),imaxv(3),jmaxv(3),kmaxv(3),
     % cflm,qqmax,imxq,jmxq,kmxq,densm,denmax,denmin
  159 format(1x,i5,2x,e10.4,4x,3(e10.3,1x),/,4x,
     % 3(e9.3,1x,i3,1x,i3,1x,i3,1x),/,4x,e9.3,
     % 2x,e9.3,1x,3(i3,1x),3(e12.6,1x))
       write(6,1234) visma/ren,vismi/ren
 1234  format(' Nu_M = ',e10.4,' Nu_m = ',e10.4)
c
       write(96,*) time, vmax(1), vmax(2), vmax(3)
      return                                                            
      end                                                               
c
c***********************************************************************
c                                                                      *
c   ************** WSON ************************                       *
c                                                                      *
c***********************************************************************
      subroutine wson(time,q1,q2,q3,pr)
      include 'param.f'
      character*3 ipfis
      character*70 namson
      REAL q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)
      REAL   pr(m1,m2,m3)
c
c     This is to refer the pressure to a fixed value
c
        ice = n1m/2+1
        jce = n2m/2+1
        prsub = pr(ice,jce,2)
c
        nul=50
      do i = 1,nson
        ir = mason(i,1)
        jr = mason(i,2)
        kr = mason(i,3)
        irm = imv(ir)
        jrm = jmv(jr)
        krm = kmv(kr)
        udx3=dx3/g3c(kr)
        udx2=dx2/g2c(jr)
        q2sur = q2(ir,jr,kr)
c
c       azimuthal vorticity
c
              dq2x3=(q2(ir,jr,kr)-q2(ir,jr,krm))*udx3
              dq3x2=(q3(ir,jr,kr)-q3(ir,jrm,kr))*udx2
              voraz=dq2x3-dq3x2
c
c       radial vorticity
c
              dq3x1=(q3(ir,jr,kr)-q3(irm,jr,kr))*dx1*iaxsy
              dq1x3=(q1(ir,jr,kr)-q1(ir,jr,krm))*dx3/g3c(kr)
              vorr=dq3x1-dq1x3
c
c       axial vorticity
c
              dq1x2=(q1(ir,jr,kr)-q1(ir,jrm,kr))
     %           *udx2
              dq2x1=(q2(ir,jr,kr)-q2(irm,jr,kr))*dx1*iaxsy
              vorz=(dq1x2-dq2x1)
        write(ipfis,99) nul
   99 format(i3.3)
      namson='son'//ipfis//'.out'
      open(99,file=namson,form='formatted',status='unknown')
      if(ntime.eq.1) go to 30
 20   read(99,100,end=10) aa,aa,aa,aa,aa,aa,aa,aa
      go to 20
 10   continue
      backspace(99)
 30     write(99,100)
     %   time,q1(ir,jr,kr),q2sur,q3(ir,jr,kr),
     %   pr(ir,jr,kr)-prsub,voraz,vorr,vorz
      close(99)
        nul=nul+1
      end do
 100  format(2x,8(e12.5,2x))
c
      return
      end     
c***********************************************************************
c                                                                      *
c   ************** INDSON ************************                     *
c                                                                      *
c***********************************************************************
      subroutine indson
c
c     This routnine finds the indices in the computational grid
c     close to the physical location of the probe.
c     This is necessary since the location of the probe is given
c     in terms of physical coordinates (sonda.in)
c
      include 'param.f'
c
      pi = 2. * asin(1.)
      do ll = 1,nson
         write(*,*)'sonda',ll
c         ierr=0
c         if((coson(ll,1).lt.0.).or.(coson(ll,1).gt.alx1)) ierr=1
c         if((coson(ll,2).lt.0.).or.(coson(ll,2).gt.alx2)) ierr=1
c         if((coson(ll,3).lt.0.).or.(coson(ll,3).gt.alx3)) ierr=1
c      if(ierr.eq.1) then
c         write(6,*) 'WARNING: la sonda ',ll,' e'' fuori dal dominio '
c         write(32,*) 'WARNING: la sonda ',ll,' e'' fuori dal dominio '
c         mason(ll,1)=1
c         mason(ll,2)=1
c         mason(ll,3)=1
c      else
c
c     SEARCHING FOR THE INDICES IN THE GRID CLOSE TO THE PROBE
c
         do j=1,n2m
            jp = j+1
            if((x2c(j).le.coson(ll,2)).and.
     %         (x2c(jp).gt.coson(ll,2))) then
               mason(ll,2)=j
               write(*,*)'trovato',mason(ll,2),j,jp
            end if
         end do
         write(*,*)'trovato',mason(ll,2)
         do k=1,n3m
            kp = k+1
            if((x3c(k).le.coson(ll,3)).and.
     %         (x3c(kp).gt.coson(ll,3))) then
               mason(ll,3)=k
               write(*,*)'trovato',mason(ll,3),k,kp
            end if
         end do
         write(*,*)'trovato',mason(ll,3)
         do i=1,n1m
            ip = ipv(i)
c            if(coson(ll,1).gt.1.e-02) then
               if((x1c(i).le.coson(ll,1)).and.
     %            (x1c(ip).gt.coson(ll,1))) then
                  mason(ll,1)=i
               write(*,*)'trovato',mason(ll,1),i,ip
               end if
c            else
c               mason(ll,1)=1
c            end if
         end do
c      end if
         write(*,*)'trovato',mason(ll,1)
      end do
c
      return
      end     
c                                                                       
c                                                                       
c***********************************************************************
      subroutine vort(voamax,vormax,vozmax,q1,q2,q3)                    
      include'param.f'
      REAL q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)                  
c                                                                       
      voamax=0.   ! vorticit… azimutale !   
      vmmx=0.                 
      do 2 kc=1,n3m           
        km=kmv(kc)           
        do 2 jc=2,n2m       
          jm=jmv(jc)       
          udx3=dx3/g3c(kc)
          udx2=dx2/g2c(jc)
            do 2 ic=1,n1m 
              dq2x3=(q2(ic,jc,kc)-q2(ic,jc,km))*udx3
              dq3x2=(q3(ic,jc,kc)-q3(ic,jm,kc))*udx2    
              voraz=dq2x3-dq3x2                 
              if(abs(voraz).gt.vmmx) then      
              jzm=jc                
              kzm=kc               
              vmmx=abs(voraz)  
              if(vmmx.gt.voamax) then  
               voamax=abs(voraz)      
               iazm=ic               
               jazm=jc              
               kazm=kc             
              endif               
              endif              
    2 continue                  
c                              
      vormax=0.      ! vorticit… radiale !  
      do 4 kc=1,n3m                         
        km=kmv(kc)                         
        do 4 jc=1,n2m                     
        do 4 ic=1,n1m                    
          im=imv(ic)                    
          dq3x1=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1*iaxsy
          dq1x3=(q1(ic,jc,kc)-q1(ic,jc,km))*dx3/g3c(kc)             
          vorr=dq3x1-dq1x3                                                  
          if(abs(vorr).gt.vormax) then 
              irrm=ic                 
              jrrm=jc                
              krrm=kc               
              vormax=abs(vorr)     
          endif                   
    4 continue                   
c                               
      vozmax=0.     ! vorticit… verticale !
      do 3 kc=2,n3m                         
        do 3 jc=2,n2m                         
          jm=jmv(jc)                           
          udx2 = dx2/g2c(jc)                               
            do 3 ic=1,n1m                       
            im=imv(ic)                           
            dq1x2=(q1(ic,jc,kc)-q1(ic,jm,kc)) 
     %           *udx2
            dq2x1=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1*iaxsy
            vorz=(dq1x2-dq2x1)
            if(abs(vorz).gt.vozmax) then 
              irzm=ic 
              jrzm=jc 
              krzm=kc 
              vozmax=abs(vorz)
            endif            
    3 continue              
      return               
      end                 
c 
cccccccccccccccccccccccccccccccccccccccccccccc
c      Max of velocity on boundary           c
cccccccccccccccccccccccccccccccccccccccccccccc

       subroutine errcal(q1,q2,q3,time)
       include 'param.f'
       REAL q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)

       umax = 0.
       uhat = 0.
       uhat1 = 0.
       suptot = 0.
        do i=1,n1m
         do l=1,npuns
c         j=indgeos(l,1)
c         k=indgeos(l,2)
          usquare=(q1(i,j,k))**2+q2(i,j,k)**2+(q3(i,j,k))**2
          umod=sqrt(usquare)
          umax=max(umax,umod)
          uhat1=uhat1+usquare*g2m(j)*g3m(k)/(dx2*dx3)
          suptot = suptot + g2m(j)*g3m(k)/(dx2*dx3)
         end do
        end do
       uhat=sqrt(uhat1/float(n1m)/suptot)
       umaxt = 0.
       do k=1,n3m
         do j=2,n2m
           do i=1,n1m
             usquare=(q1(i,j,k))**2+q2(i,j,k)**2+(q3(i,j,k))**2
             umod=sqrt(usquare)
             umaxt=max(umaxt,umod)
           end do
         end do
       end do
       write(40,*) time, umax, uhat,umaxt
        umax = 0.
        do i=1,n1m
         do l=1,npuns
c         j=indgeos(l,1)
c         k=indgeos(l,2)
          usquare=(q1(i,j,k))**2+q2(i,j,k)**2+(q3(i,j,k))**2
          umod=sqrt(usquare)
          if(umod.gt.umax) then
            umax=umod
            jvm = j
            kvm = k
          end if
         end do
        end do
        write(6,*) ' UMAX = ',umax,' at J = ',jvm,' K = ',kvm
        write(6,*) ' rms U_b = ',uhat,' UMAX_t = ',umaxt
       return
       end 
c
c***********************************************************************
      subroutine tripvmy(n1i,n1f,n2i,n2f,n3i,n3f,nn1,nn2,iequa,q1)                                    
      include'param.f'
c                                                                       
c     vectorized for right hand side and coefficients                   
c                                                                       
      dimension fn(m2),p(m2)                                            
      REAL q1(m1,m2,m3)
      ia = n1i + 1                                                       
      ii = n1i + n1f                                                      
      do k=n3i,n3f
      kc = k
c
c   COEFFICIENTS FOR TRIDIAGONAL INVERSION
c

c
c   Q_1
c
      if(iequa.eq.1) then
      betadx=beta*dx1q*al
      do jc=1,n2m
        ugmmm=betadx
        do ic=1,n1m
            im = imv(ic)
            ip = ipv(ic)
            api(ic,jc)=-ugmmm *(1.-forclo(ic,jc,kc))*visc(kc)
            ami(ic,jc)=-ugmmm *(1.-forclo(ic,jc,kc))*visc(kc)
            aci(ic,jc)=1.+ugmmm*2.*(1.-forclo(ic,jc,kc))*visc(kc)
         end do
       end do
       end if
c
c  Q_2
c
      if(iequa.eq.2) then
      betadx=beta*dx1q*al
      do jc=2,n2m
        jm = jc - 1
        jp = jc + 1
        ugmmm=betadx
         do ic=1,n1m
            im = imv(ic)
            ip = ipv(ic)
            api(ic,jc)=-ugmmm*(1.-forclo(ic,jc,kc))*visc(kc)
            ami(ic,jc)=-ugmmm*(1.-forclo(ic,jc,kc))*visc(kc)
            aci(ic,jc)=1.+ugmmm*2.*(1.-forclo(ic,jc,kc))*visc(kc)
         end do
       end do
      do ic=1,n1m
            api(ic,1)=0.
            ami(ic,1)=0.
            aci(ic,1)=1.
            api(ic,n2)=0.
            ami(ic,n2)=0.
            aci(ic,n2)=1.
       end do
       end if
c
c  Q_3
c
      if(iequa.eq.3) then
      km = kmv(kc)
      betadx=beta*dx1q*al
      do jc=1,n2m
        ugmmm=betadx
        do ic=1,n1m
            im = imv(ic)
            ip = ipv(ic)
            api(ic,jc)=-ugmmm*(1.-forclo(ic,jc,kc))*visc(kc)
            ami(ic,jc)=-ugmmm*(1.-forclo(ic,jc,kc))*visc(kc)
            aci(ic,jc)=1.+ugmmm*2.*(1.-forclo(ic,jc,kc))*visc(kc)
         end do
       end do
       end if
c
c  RHO
c
      if(iequa.eq.4) then
      betadx=0.5*dx1q*al*dt/pec
      betanl=dt*dx1*al*0.25
      do jc=1,n2m
        ugmmm=betadx
        do ic=1,n1m
            im = imv(ic)
            ip = ipv(ic)
            api(ic,jc)=-ugmmm*visc(kc)
            ami(ic,jc)=-ugmmm*visc(kc)
            aci(ic,jc)=1.+ugmmm*2.*visc(kc)
         end do
       end do
       end if
c
c  THE INVERSION STARTS
c
      do 20 j=n2i,n2f                                                     
        q(j,n1i) = -api(n1i,j)/aci(n1i,j)                                        
        s(j,n1i) = -ami(n1i,j)/aci(n1i,j)                                        
        iaddf = n1f+(j-1)*nn1+(k-1)*nn1*nn2
        iaddi = n1i+(j-1)*nn1+(k-1)*nn1*nn2
        fn(j) = rhs(iaddf)                                                   
        rhs(iaddi) = rhs(iaddi)/aci(n1i,j)                                         
   20 continue                                                          
c                                                                       
c     forward elimination sweep                                         
c                                                                       
      do 10 i=ia,n1f                              
        do 21 j=n2i,n2f                             
          iadd = i+(j-1)*nn1+(k-1)*nn1*nn2
          iaddm = i-1+(j-1)*nn1+(k-1)*nn1*nn2
          p(j) =1./( aci(i,j) + ami(i,j)*q(j,i-1))
          q(j,i) = - api(i,j)*p(j)                     
          s(j,i) = - ami(i,j)*s(j,i-1)*p(j)                                   
          rhs(iadd) = ( rhs(iadd) - ami(i,j)*rhs(iaddm))*p(j)
   21   continue          
   10 continue             
c                                                                       
c     backward pass                                                     
c                                                                       
      do 22 j=n2i,n2f      
        s(j,n1f) = 1.   
        fei(j,n1f) = 0.
   22 continue        
      do 11 l=ia,n1f       
        i = ii - l       
        do 23 j=n2i,n2f   
          iadd = i+(j-1)*nn1+(k-1)*nn1*nn2
          s(j,i) = s(j,i) + q(j,i)*s(j,i+1)        
          fei(j,i) = rhs(iadd) + q(j,i)*fei(j,i+1)  
   23   continue                    
   11 continue                     
      do 24 j=n2i,n2f               
        iaddf = n1f+(j-1)*nn1+(k-1)*nn1*nn2
        rhs(iaddf)=(fn(j)-api(i,j)*fei(j,n1i) - 
     %      ami(i,j)*fei(j,n1f-1))/(api(i,j)*s(j,n1i) +
     %      ami(i,j)*s(j,n1f-1)+aci(i,j))          
   24 continue          
c                                                                       
c     backward elimination pass                                         
c                                                                       
      do 12 l=ia,n1f         
        i = ii -l          
        do 25 j=n2i,n2f     
          iadd = i+(j-1)*nn1+(k-1)*nn1*nn2
          iaddf = n1f+(j-1)*nn1+(k-1)*nn1*nn2
          rhs(iadd) = rhs(iaddf)*s(j,i) + fei(j,i)  
   25   continue                               
   12 continue                                
      end do
      return                                 
      end                           
c
c**************************************************************
c
      subroutine cordin                                                 
      include'param.f'
      dimension eta(m2),etaz(m3)
      dimension xt31(m3),xt32(m3)
      pi=2.*asin(1.)
c
c     AZIMUTHAL COORDINATE DEFINITION
c
      if(igext.ne.1) then
      do i=1,n1                                                    
        x1c(i)= float(i-1)/dx1                                     
      end do
      end if
      do i=1,n1m                                                    
      ip=i+1
      x1m(i)= (x1c(i)+x1c(ip))*.5
c       x1m(i)= (float(i-1)+0.5)/dx1 
      end do
c
c     RADIAL COORDINATE DEFINITION
c
      rint = 0.
      r0 = alx2 
      if(igext.ne.1) then
      tstr= tanh(strr)
      tstr2=tanh(strr*0.5)
      tstr3=tanh(strr*(etdp-0.5))
      tstr4=tanh(strb*(etdp-1.))
      xn11= rmed * tanh(strr*0.5)/tstr3
c
c     UNIFORM GRID
c
      if (istr.eq.-1) then
        do  j=1,n2
         x2=float(j-1)/float(n2m)
         x2c(j)= r0*x2
        end do
      end if
c
c     CLUSTERING AT THE EXTERNAL RADIAL WALL
c
      if (istr.eq.0) then
        tstr2=tanh(strr)
        do  j=1,n2
          x2=float(j-1)/float(n2m)
          eta(j)=r0*(1.+tanh(strr*(x2-1.))/tstr2)
          x2c(j)=eta(j)
        enddo
      end if
c
c     CLUSTERING AT THE AXIS
c
      if (istr.eq.1) then
        do j=1,n2
          x2dp=float(j-1)/float(n2m)
          rcdp=tanh(strr*x2dp)/tstr
          x2c(j)=rcdp*alx2
        end do
      end if
c
c     VARIABLE CLUSTERING
c
      if (istr.eq.2) then
        do j=1,n2
          x2=(j-1)/dx2
          x22 = 0.5 + x2*0.5
          xx1= rmed * tanh(strr*(x22- 0.5))/tstr3
          xx2= 1./xn11 + (1. - 1./xn11)* tanh(strb*(x22-1.))/tstr4
          eta(j)= xx1*xx2
          x2c(j)=rint+eta(j)*(alx2-rint)
        end do

c       open(66,file='pippo.dat',form='formatted')
c       do j=1,n2
c         x2=(j-1)/dx2
c         x2=(j-1)/85.
c         x22 = 0.5 + x2*0.5
c         xx1= rmed * tanh(strr*(x22- 0.5))/tstr3
c         xx2= 1./xn11 + (1. - 1./xn11)* tanh(strb*(x22-1.))/tstr4
c         eta(j)= xx1*xx2
c         x2c(j)=rint+eta(j)*(alx2-rint)
c         write(66,*) j,x2c(j)
c       end do
c       close(66)

CCC   Aggiunta

      end if
C
C    MESH FROM AN EXTERNAL INPUT FILE
C
      if (istr.eq.3) then
        open(66,file='gridy.data',form='formatted')
        do j=1,n2
          read(66,*) x2,x2c(j)
        end do
        close(66)
      end if
      x2c(1)=rint
c
c     PETROS PROBLEM
c
      if (istr.eq.4) then
        dx2t=n2m/2
        alx2h = alx2*0.5
        n2pe = n2m/2+1
c       tstr2=tanh(strr)
c       do  j=1,n2pe
c         x2=float(j-1)/dx2t
c         eta(j)=alx2h*(1.+tanh(strr*(x2-1.))/tstr2)
c         x2c(j)=eta(j)
c       enddo
        tstr= tanh(strr)
        do j=1,n2pe
          x2dp=float(j-1)/dx2t
          rcdp=tanh(strr*x2dp)/tstr
          x2c(j)=rcdp*alx2h
        end do
        do j=n2pe+1,n2
          ji = j - n2pe
          x2c(j)= x2c(n2pe) + (x2c(n2pe)-x2c(n2pe-ji))
        end do
      end if
c
c     CHANNEL
c
      if (istr.eq.5) then
      tstr2=tanh(strr*0.5)
      do j=1,n2
      x2=float(j-1)/float(n2m)
      eta(j)=0.5*(1.+tanh(strr*(x2-0.5))/tstr2)
      end do
      do j=1,n2
      x2c(j)=+eta(j)*alx2
      end do
      end if
      end if
c
c     STAGGERED COORDINATES AND
c     METRIC QUANTITIES
c
      do j=1,n2m
        x2m(j)=(x2c(j)+x2c(j+1))*0.5
        g2m(j)=(x2c(j+1)-x2c(j))*dx2
      end do
      do j=2,n2m
        g2c(j)=(x2c(j+1)-x2c(j-1))*dx2*0.5
      end do
      g2c(1)=(x2c(2)-x2c(1))*dx2
      g2c(n2)=(x2c(n2)-x2c(n2m))*dx2
c
c     AXIAL COORDINATE DEFINITION
c
      if(igext.ne.1) then
      tstr3=tanh(str3)
c
c     UNIFORM GRID
c
      if (istr3.lt.0) then
        do k=1,n3
          x3=float(k-1)/float(n3m)
          etaz(k)=alx3*x3
          x3c(k)=etaz(k)
        enddo
      endif
c
c     CLUSTERING AT THE SOUTH WALL
c
      if (istr3.eq.0) then
        do k=1,n3
          x3=float(k-1)/float(n3m)
          etaz(k)=alx3*(1.+tanh(str3*(x3-1.0))/tstr3)
          x3c(k)=etaz(k)
        enddo
      endif
c
c     CLUSTERING AT THE NORTH WALL
c
      if (istr3.eq.1) then
        do k=1,n3
          x3=float(k-1)/float(n3m)
          etaz(k)=alx3*tanh(str3*x3)/tstr3
          x3c(k)=etaz(k)
        enddo
      endif
c
c     VARIABLE CLUSTERING
c
      if (istr3.eq.2) then
        dx3t=n3m
        tstr3=tanh(str3*etdp3)
        do k=1,n3
          x3=float(k-1)/dx3t
          xt31(k)=rmed31/alx3*tanh(str3*x3)/tstr3
        end do
        do k=1,n3
          x3=(k-1)/dx3t
          xt32(k)=1./xt31(n3)+(1.-1./xt31(n3))
     %    *tanh(strb3*(x3-1.))/tanh(strb3*(etdp3-1.))
          etaz(k)= xt31(k)*xt32(k)
          x3c(k)=etaz(k)*alx3
        end do
      end if
c
c     PETROS PROBLEM
c
      if (istr3.eq.4) then
        dx3t=n3m/2
        alx3h = alx3*0.5
        n3pe = n3m/2+1
        tstr3=tanh(str3*etdp3)
        do k=1,n3pe
          x3=float(k-1)/dx3t
          xt31(k)=rmed31/alx3h*tanh(str3*x3)/tstr3
        end do
        do k=1,n3pe
          x3=(k-1)/dx3t
          xt32(k)=1./xt31(n3pe)+(1.-1./xt31(n3pe))
     %    *tanh(strb3*(x3-1.))/tanh(strb3*(etdp3-1.))
          etaz(k)= xt31(k)*xt32(k)
          x3c(k)=etaz(k)*alx3h
        end do
        do k=n3pe+1,n3
          ki = k - n3pe
          x3c(k)= x3c(n3pe) + (x3c(n3pe)-x3c(n3pe-ki))
        end do
      end if
C
C    MESH FROM AN EXTERNAL INPUT FILE
C
      if (istr3.eq.3) then
        open(66,file='gridz.data',form='formatted')
        do k=1,n3
          read(66,*) x3,x3c(k)
        end do
        close(66)
      end if
      end if
c
c     STAGGERED COORDINATES AND
c     METRIC QUANTITIES
c
      do k=1,n3m
        x3m(k)=(x3c(k)+x3c(k+1))*0.5
        g3m(k)=(x3c(k+1)-x3c(k))*dx3
      enddo
      do k=2,n3m
        g3c(k)=(x3c(k+1)-x3c(k-1))*dx3*0.5
      enddo
      g3c(1)=(x3c(2)-x3c(1))*dx3
      g3c(n3)= (x3c(n3)-x3c(n3m))*dx3
c
c     WRITE GRID INFORMATION
c
      open(unit=98,file='radcor.out',status='unknown')
      do j=1,n2
        write(98,348) j,x2c(j),x2m(j),g2c(j),g2m(j)
      end do
      close(98)
      open(unit=78,file='axicor.out',status='unknown')
      do k=1,n3
        write(78,348) k,x3c(k),x3m(k),g3c(k),g3m(k)
      end do
      close(78)
 348  format(2x,i3,4(2x,e14.7))
c
c     LES length scales
c
      if(iles.gt.0) then
        do kc=1,n3m
          deltaz1=g3m(kc)/dx3
          deltaz2=2.*deltaz1
          do jc=1,n2m
            deltar1=g2m(jc)/dx2
            deltar2=2.*deltar1
            deltat1=1./dx1
            deltat2=2.*deltat1
            ell1(jc,kc)=(deltat1*deltaz1*deltar1)**.66667
            ell2(jc,kc)=(deltat2*deltaz2*deltar2)**.66667
          end do
        end do
      end if
c
      return                                                            
      end                                                               
c                                                                       
c***********************************************************************

      subroutine outpf(time,q1,q2,q3,pr,dq)
      include'param.f'
      REAL q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)                  
      REAL      pr(m1,m2,m3)
      REAL    dq(m1,m2,m3)
      REAL    stime
      character*20 namfi3
      character*20 namfi4
      character*20 namfi5
      character*20 namfi6
      character*6 ipfi
      do k=1,n3                                                  
        do j=1,n2                                                    
          do i=1,n1                                                   
            dq(i,j,k) = 0.
          end do
        end do
      end do
      ice = n1m/2+1
      jce = n2m/2+1
      prsub = pr(ice,jce,2)
      do k=1,n3m
        do j=1,n2m                                                   
          do i=1,n1m                                                  
            dq(i,j,k) = pr(i,j,k) - prsub
          end do
        end do
      end do
      stime = time
      itime=nint(1000.*time)
      write(ipfi,99) itime
   99 format(i6.6)
      namfi3='coge'//ipfi//'.dat'
      write(6,201)time,namfi3
  201 format(10x,'At t=',e10.3,' write data on ',a20)
      n1pp=(n1-1)/n1p+1                                                 
      n2pp=(n2-1)/n2p+1                                                 
      n3pp=(n3-1)/n3p+1                                                 
      open(99,file=namfi3,form='unformatted')
        rewind 99
        write(99) n1pp,n2pp,n3pp 
        write(99) stime,lamb,stime,stime 
        write(99) 
     %   ((((q1(i,j,k)),i=1,n1,n1p),j=1,n2,n2p),k=1,n3,n3p),
     %   ((((q2(i,j,k)),i=1,n1,n1p),j=1,n2,n2p),k=1,n3,n3p),
     %   ((((q3(i,j,k)),i=1,n1,n1p),j=1,n2,n2p),k=1,n3,n3p),
     %   ((((dq(i,j,k)),i=1,n1,n1p),j=1,n2,n2p),k=1,n3,n3p)
      close(99)
      
c      namfi3='pressione'//ipfi//'.dat'
c      open(99,file=namfi3,form='formatted')
c        rewind 99
c      write(99,*) n3,n2,n1,1
c      write(99,*)(((dq(i,j,k),k=1,n3),j=1,n2),i=1,n1)
c      close(99)

      namfi3='slOUT'//ipfi//'.dat'
      open(99,file=namfi3,form='formatted',status='unknown')
      write(99,*) 'VARIABLES = "X","Y","Z","Vx","Vy","Vz",
     &  "dVx","dVy","dVz"'
      write(99,*) 'ZONE I = 1, J = ',n2, ', K = ',n1, ' F=POINT'
      do i=1,n1
        do j=1,n2
          write(99,*)
     $           x3c(n3),x2c(j),x1c(i),
     $           qb3n(i,j),qb2n(i,j),qb1n(i,j),
     $           dqb3n(i,j),dqb2n(i,j),dqb1n(i,j)
        enddo
      enddo
      close(99)
c      namfi3='slIN'//ipfi//'.dat'
c      open(99,file=namfi3,form='formatted',status='unknown')
c      write(99,*) 'VARIABLES = "X","Y","Z","Vx","Vy","Vz",
c     &"dVx","dVy","dVz"'
c      write(99,*) 'ZONE I = 1, J = ',n2, ', K = ',n1, ' F=POINT'
c      do i=1,n1
c        do j=1,n2
c          write(99,*)
c     $           x3c(n3),x2c(j),x1c(i),
c     $           qb3s(i,j),qb2s(i,j),qb1s(i,j),
c     $           dqb3s(i,j),dqb2s(i,j),dqb1s(i,j)
c        enddo
c      enddo
c      close(99)
c      
c      namfi3='slDN-UP'//ipfi//'.dat'
c      open(99,file=namfi3,form='formatted',status='unknown')
c      write(99,*) 'VARIABLES = "X","Y","Z","VxD","VyD","dVxD","dVyD",
c     &"VxU","VyU","dVxU","dVyU"'
c      write(99,*) 'ZONE I = ',n3,', J = 1, K = ',n1, ', F=POINT'
c      do i=1,n1
c        do k=1,n3
c          write(99,*)
c     $           x3c(k),x2c(n2),x1c(i),
c     $           qb3dn(i,j),qb2dn(i,j),dqb3dn(i,j),dqb2dn(i,j),
c     $           qb3up(i,j),qb2up(i,j),dqb3up(i,j),dqb2up(i,j)
c        enddo
c      enddo
c      close(99)

      namfi3='sl-i9'//ipfi//'.dat'
      open(99,file=namfi3,form='formatted',status='unknown')
      write(99,*) 'VARIABLES = "X","Y","Z","Vx","Vy","Vz",
     &"VxAss","VyAss","Pres"'
      write(99,*) 'ZONE I = ',n3, ', J = ',n2, ', K = 1, F=POINT'
      do j=1,n2
        do k=1,n3
          write(99,*)
     $       x3c(k),x2c(j),x1c(9),
     $       q3(9,j,k),q2(9,j,k),q1(9,j,k),
     $       ( q3(9,j,k)*cos(alphax) + q2(9,j,k)*sin(alphax) ),
     $       ( -q3(9,j,k)*sin(alphax) + q2(9,j,k)*cos(alphax) ),
     $       pr(9,j,k)
         enddo
      enddo
      close(99)

c      namfi3='sl-j94'//ipfi//'.dat'
c      open(99,file=namfi3,form='formatted',status='unknown')
c      write(99,*) 'VARIABLES = "X","Y","Z","Vx","Vy","Vz",
c     $     "VxAss","VyAss","Pres"'
c      write(99,*) 'ZONE I = ',n3,', J = 1, K = ',n1, ', F=POINT'
c      j=94
c      do i=1,n1
c        do k=1,n3
c          write(99,*)
c     $           x3c(k),x2c(j),x1c(i),
c     $       q3(i,j,k),q2(i,j,k),q1(i,j,k),
c     $       ( q3(i,j,k)*cos(alphax) + q2(i,j,k)*sin(alphax) ),
c     $       ( -q3(i,j,k)*sin(alphax) + q2(i,j,k)*cos(alphax) ),
c     $       pr(i,j,k)
c        enddo
c      enddo
c      close(99)

c      namfi3='sl-j109'//ipfi//'.dat'
c      open(99,file=namfi3,form='formatted',status='unknown')
c      write(99,*) 'VARIABLES = "X","Y","Z","Vx","Vy","Vz",
c     $     "VxAss","VyAss","Pres"'
c      write(99,*) 'ZONE I = ',n3,', J = 1, K = ',n1, ', F=POINT'
c      j=109
c      do i=1,n1
c        do k=1,n3
c          write(99,*)
c     $           x3c(k),x2c(j),x1c(i),
c     $       q3(i,j,k),q2(i,j,k),q1(i,j,k),
c     $       ( q3(i,j,k)*cos(alphax) + q2(i,j,k)*sin(alphax) ),
c     $       ( -q3(i,j,k)*sin(alphax) + q2(i,j,k)*cos(alphax) ),
c     $       pr(i,j,k)
c        enddo
c      enddo
c      close(99)
c      namfi3='sl-j80'//ipfi//'.dat'
c      open(99,file=namfi3,form='formatted',status='unknown')
c      write(99,*) 'VARIABLES = "X","Y","Z","Vx","Vy","Vz",
c     $     "VxAss","VyAss","Pres"'
c      write(99,*) 'ZONE I = ',n3,', J = 1, K = ',n1, ', F=POINT'
c      j=80
c      do i=1,n1
c        do k=1,n3
c          write(99,*)
c     $           x3c(k),x2c(j),x1c(i),
c     $       q3(i,j,k),q2(i,j,k),q1(i,j,k),
c     $       ( q3(i,j,k)*cos(alphax) + q2(i,j,k)*sin(alphax) ),
c     $       ( -q3(i,j,k)*sin(alphax) + q2(i,j,k)*cos(alphax) ),
c     $       pr(i,j,k)
c        enddo
c      enddo
c      close(99)
c      namfi3='sl-j169'//ipfi//'.dat'
c      open(99,file=namfi3,form='formatted',status='unknown')
c      write(99,*) 'VARIABLES = "X","Y","Z","Vx","Vy","Vz",
c     $     "VxAss","VyAss","Pres"'
c      write(99,*) 'ZONE I = ',n3,', J = 1, K = ',n1, ', F=POINT'
c      j=169
c      do i=1,n1
c        do k=1,n3
c          write(99,*)
c     $           x3c(k),x2c(j),x1c(i),
c     $       q3(i,j,k),q2(i,j,k),q1(i,j,k),
c     $       ( q3(i,j,k)*cos(alphax) + q2(i,j,k)*sin(alphax) ),
c     $       ( -q3(i,j,k)*sin(alphax) + q2(i,j,k)*cos(alphax) ),
c     $       pr(i,j,k)
c        enddo
c      enddo
c      close(99)
c      namfi3='u2-9j1.dat'
c      open(99,file=namfi3,form='formatted',status='unknown')
c      write(99,*) 'VARIABLES = "X2c","X2m","Q2","Q3"'
c      do j=1,n2
c         write(99,*)
c     $       x2c(j),x2m(j),q2(9,j,1),q2(9,j,2),q3(9,j,1),q3(9,j,2)
c      enddo
c      close(99)
c      namfi3='u2-9j1W.dat'
c      open(99,file=namfi3,form='formatted',status='unknown')
c      write(99,*) 'VARIABLES = "X2c","X2m","Q2","Q3"'
c      do j=1,n2
c         write(99,*)
c     $       x2c(j),x2m(j),q2(9,j,n3m),q2(9,j,n3m-1),q3(9,j,n3),
c     $        q3(9,j,n3m)
c      enddo
c      close(99)     
      if(iles.ne.0) then 
        namfi3='visct'//ipfi//'.dat'
        n1pp=(n1-1)/n1p+1                                                 
        n2pp=(n2-1)/n2p+1                                                 
        n3pp=(n3-1)/n3p+1                                                 
        open(99,file=namfi3,form='unformatted')
        rewind 99
        write(99) n1pp,n2pp,n3pp ,1
        write(99) 
     %  ((((visct(i,j,k)),i=1,n1,n1p),j=1,n2,n2p),k=1,n3,n3p)
        close(99)
      end if
      return                                                            
      end                                                               
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c     Drag Lift and Side force coefficient calculation
c
c 
c***********************************************************************
      subroutine dralical(time,q1,q2,q3,pr)
      include'param.f'
      common/parlea/ p_ew,p_ee,d_lw,p_w,p_d,p_s
      common/forcec/alunt,aliun,dunt,drun,
     % dunto1,dunto2,dunto3,alunto1,alunto2,alunto3,
     % siunto1,siunto2,siunto3,sifo1,sifo2,sifo3
     % drag1,drag2,drag3,alift1,alift2,alift3
      REAL q3(m1,m2,m3)
      REAL q2(m1,m2,m3)
      REAL q1(m1,m2,m3)
      REAL   pr(m1,m2,m3)
      do k=1,n3   
        do j=1,n2
          do i=1,n1
            forclo(i,j,k) = 1.
          end do
        end do
      end do
      do n=1,npunr
        i=indgeo(2,n,1)
        j=indgeo(2,n,2)
        k=indgeo(2,n,3)
        ie=indgeoe(2,n,1)
        je=indgeoe(2,n,2)
        ke=indgeoe(2,n,3)
        forclo(i,j,k) = 0.
        forclo(ie,je,ke) = 0.
      end do
c                                                                       
      if(ntime.eq.ntiia) timeo = 0.
c
      do l=1,3
      if(l.eq.1) then
      kin = n3m/10
      kfi = n3m-n3m/10
      end if
      if(l.eq.2) then
      kin = n3m/12
      kfi = n3m-n3m/12
      end if
      if(l.eq.3) then
      kin = n3m/14
      kfi = n3m-n3m/14
      end if
      jin = 1
      jfi = n2m
      drag = 0.
      alift = 0.
      sifo = 0.
c
c
c     integration along vertical box walls
c
      do j=jin,jfi
      jp = jpv(j)
      jm = jmv(j)
      do i=1,n1m
      im = imv(i)
      ip = ipv(i)
      drag = drag + 2.*(
     %   q3(i,j,kfi)**2-q3(i,j,kin)**2
     % -(q3(i,j,kfi+1)-q3(i,j,kfi-1))*2.*dx3/(g3c(kfi)*ren)
     % +(q3(i,j,kin+1)-q3(i,j,kin-1))*2.*dx3/(g3c(kin)*ren)
     % + (pr(i,j,kfi) + pr(i,j,kfi-1))*0.5
     % - (pr(i,j,kin) + pr(i,j,kin-1))*0.5
     %  )*g2m(j)/(dx2*dx1)
      alift = alift + 2.*(
     % (q3(i,j,kfi)+q3(i,jm,kfi))*
     % (q2(i,j,kfi)+q2(i,j,kfi-1))*0.25
     %-(q3(i,j,kin)+q3(i,jm,kin))*
     % (q2(i,j,kin)+q2(i,j,kin-1))*0.25
     % -((q3(i,j,kfi)-q3(i,jm,kfi))*dx2/(g2c(j)*ren)
     %  +(q2(i,j,kfi)-q2(i,j,kfi-1))*dx3/(g3c(kfi)*ren))*0.5
     % +((q3(i,j,kin)-q3(i,jm,kin))*dx2/(g2c(j)*ren)
     %  +(q2(i,j,kin)-q2(i,j,kin-1))*dx3/(g3c(kin)*ren))*0.5
     %  )*g2c(j)/(dx2*dx1) 
      sifo = sifo + 2.*(
     % (q3(i,j,kfi)+q3(im,j,kfi))*
     % (q1(i,j,kfi)+q1(i,j,kfi-1))*0.25
     %-(q3(i,j,kin)+q3(im,j,kin))*
     % (q1(i,j,kin)+q1(i,j,kin-1))*0.25
     % -((q1(i,j,kfi)-q1(i,j,kfi-1))*dx3/(g3c(kfi)*ren)
     %  +(q3(i,j,kfi)-q3(im,j,kfi))*dx1/ren)*0.5
     % +((q1(i,j,kin)-q1(i,j,kin-1))*dx3/(g3c(kin)*ren)
     %  +(q3(i,j,kfi)-q3(im,j,kfi))*dx1/ren)*0.5
     %  )*g2c(j)/(dx2*dx1) 
      end do
      end do
c
c     integration along horizontal box walls
c
      ice = n1m/2+2
      jce = n2m/2+1
      prsub = pr(ice,jce,2)
      do k=kin,kfi
      do i=1,n1m
      alift = alift + 2.*(
     % +(pr(i,jfi,k)-prsub)
     % -(pr(i,jin,k)-prsub)*forclo(i,jin,k)
     %  )*g3m(k)/(dx3*dx1)
      drag = drag + 2.*(
     % +q3(i,1,k)*2.*dx2/(g2c(1)*ren)*forclo(i,jin,k)
     %  )*g3m(k)/(dx3*dx1)
      sifo = sifo + 2.*(
     % +q1(i,1,k)*2.*dx2/(g2c(1)*ren)*forclo(i,jin,k)
     %  )*g3m(k)/(dx3*dx1)
      end do
      end do
c
c     unsteady terms
c
      dunt = 0.
      alunt = 0.
      siunt = 0.
      do k=kin,kfi
        do j=jin,jfi
          do i=1,n1m
            dunt = dunt + 2.*q3(i,j,k)
     %      *g2m(j)*g3c(k)/(dx1*dx2*dx3)
            alunt = alunt + 2.*q2(i,j,k)
     %      *g2c(j)*g3m(k)/(dx1*dx2*dx3)
            siunt = siunt + 2.*q1(i,j,k)
     %      *g2m(j)*g3m(k)/(dx1*dx2*dx3)
          end do
        end do
      end do
      if(l.eq.1) then
      if(ntime.ne.ntiia) then
        drun = (dunt-dunto1)/(time-timeo)
        aliun = (alunt-alunto1)/(time-timeo)
        sifun = (siunt-siunto1)/(time-timeo)
      else
        drun = 0.
        aliun = 0.
        sifun = 0.
      end if
      dunto1 = dunt
      alunto1 = alunt
      siunto1 = siunt
      end if
      if(l.eq.2) then
      if(ntime.ne.ntiia) then
        drun = (dunt-dunto2)/(time-timeo)
        aliun = (alunt-alunto2)/(time-timeo)
        sifun = (siunt-siunto2)/(time-timeo)
      else
        drun = 0.
        aliun = 0.
        sifun = 0.
      end if
      dunto2 = dunt
      alunto2 = alunt
      siunto2 = siunt
      end if
      if(l.eq.3) then
      if(ntime.ne.ntiia) then
        drun = (dunt-dunto3)/(time-timeo)
        aliun = (alunt-alunto3)/(time-timeo)
        sifun = (siunt-siunto3)/(time-timeo)
      else
        drun = 0.
        aliun = 0.
        sifun = 0.
      end if
      dunto3 = dunt
      alunto3 = alunt
      siunto3 = siunt
      end if
      drag = drag + drun
      alift = alift + aliun
      sifo = sifo + sifun
      if(l.eq.1) then
      drag1 = drag
      alift1 = alift
      sifo1 = sifo
      end if
      if(l.eq.2) then
      drag2 = drag
      alift2 = alift
      sifo2 = sifo
      end if
      if(l.eq.3) then
      drag3 = drag
      alift3 = alift
      sifo3 = sifo
      end if
      end do
c
      timeo = time
c
      write(65,1000) time,drag1,alift1,sifo1,
     % drag2,alift2,sifo2,drag3,alift3,sifo3
      write(6,1100) time,drag1,alift1,sifo1
      write(6,1200) time,drag2,alift2,sifo2
      write(6,1300) time,drag3,alift3,sifo3
1000  format(10(2x,e12.5))
1100  format(' (1) ',4(2x,e12.5))
1200  format(' (2) ',4(2x,e12.5))
1300  format(' (3) ',4(2x,e12.5))
      return                                                            
      end                                                               
