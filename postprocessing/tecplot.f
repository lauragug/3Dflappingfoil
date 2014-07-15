cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c          DESTAggerizzazione della griglia per plot3d                 c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c
      program desta              ! careful: il file di dati e' diverso !
      include 'parade_new.f'

      character*4 dummy

      open(unit=15,file='coge.in',status='old')
        read(15,301) dummy                                                
        read(15,*) n1,n2,n3,nsst,nwrit,nread
        read(15,301) dummy                                                
        read(15,*) n1p,n2p,n3p                                            
        read(15,301) dummy                                                
        read(15,*) alx3,istr3,str3,rmed31,etdp3,strb3
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
        read(15,*)iii,iii
        read(15,301) dummy
        read(15,*)iii,iii,iii   
        read(15,301) dummy
        read(15,*)igext
        read(15,301) dummy
        read(15,*)uzero,amplit,thetamax
301     format(a4)                
      close(15)
c
      rext = alx2
      n1m=n1-1    
      if(n1.eq.1) n1m=1
      n2m=n2-1   
      n3m=n3-1  

c     assegnazione grandezze moto libero corpo MOTION      
      pi=2.0*asin(1.0)
      osci = -0.*pi/180.
      phi = 90.*pi/180.
      sigma = pi*0.35/0.75
      x3b_0_in = 0.333
      x2b_0_in = 0.
      ampp = 0.
c     assegnazione grandezze moto libero corpo MOTION      
      call gcurv 

      stop      
      end   
c                                                                       
c  **************  function ftir                                        
c                                                                       
      function  ftir(tina)
      if(tina.ge.1.) then
      ftir=1.
                     else
      ftir=3.*tina**2-2.*tina**3
                     endif
      return
      end
      
      function dftir(tina)
      if(tina.ge.1.) then
      dftir=0.
                     else
      dftir=6.*tina-6.*tina**2
                     endif
      return
      end
      
      function ddftir(tina)
      if(tina.ge.1.) then
         ddftir=0.
      else
         ddftir=6.-12.*tina**2
      endif
      return
      end
c
c***********************************************************************
c
      subroutine gcurv    
      include 'parade_new.f'

      dimension q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)   
      dimension vo1(m1,m2,m3),vo2(m1,m2,m3),vo3(m1,m2,m3)
      dimension voc1(m1,m2,m3),voc2(m1,m2,m3),voc3(m1,m2,m3)
      dimension eig1(m1,m2,m3),eig2(m1,m2,m3)
      dimension pr(m1,m2,m3),prc(m1,m2,m3)
      dimension dens(m1,m2,m3),densc(m1,m2,m3)
      dimension vora(m1,m2,m3) 

      character*5 nfil



      nvort = 0
      neigen = 0
      write(6,*)' 1 per calcolare la vorticita'' '
      read(5,*) nvort
      write(6,*)' 1 per calcolare lambda_2 '
      write(6,*)' 2 per calcolare lambda_ci e lambda_cr/lambda_ci '
      read(5,*) neigen
     
      call topogr
      call meshes
      call cordin
      call indic
      write(6,754)n1,n2,n3            
  754 format(/,4x,'Griglia :',2x,i3,'x',i3,'x',i3,//) 
      call input(time,q1,q2,q3,pr,dens,nfil)    

C     MOTION
      angle = osci*sin( sigma*time + phi )
      omegak = osci*sigma*cos( sigma*time + phi )
      x3b_0 = x3b_0_in
      x2b_0 = x2b_0_in + ampp*sin( sigma*time )  
C     MOTION

      call rotation (time)
 
      if(nvort.eq.1)call vort(vora,q1,q2,q3,vo1,vo2,vo3,time)

      if(neigen.eq.1.or.neigen.eq.2)
     $     call eigen(neigen,q1,q2,q3,eig1,eig2,time)

      write(6,*) ' 1 se vuoi calcolare le forze sul corpo '
      read(5,*)nforze
      if(nforze.eq.1) 
     $     call svelt(q1,q2,q3,pr,dens,voc1,voc2,voc3,prc,
     $        densc,time) 
c      goto 6789   
      call outpf1(time,vo1,vo2,vo3,voc1,voc2,voc3)  
c      call outpf(time,vo1,vo2,vo3,voc1,voc2,voc3,pr,densc,
c     %                               vora,eig1)
c 6789 continue       
      write(6,*)' 1 se vuoi stampare la griglia '
      read(5,*)iprintG
      if(iprintG.eq.1)call pricor(dens,time)
      
      if(neigen.eq.1.or.neigen.eq.2)
     $     call outeigen(neigen,eig1,eig2,time)


      return
      end  
c         
c***********************************************************************
c
      subroutine indic  
      include 'parade_new.f'

      n1mm=n1m-1     
      do 1 ic=1,n1m 
        ip=ic+1    
        imv(ic)=ic-1 
        if(ic.eq.1) imv(ic)=n1m   
        ipv(ic)=ic+1             
        if(ic.eq.n1m) ipv(ic)=1 
    1 continue                 
      do 4 kc=1,n3m           
        kmv(kc)=kc-1         
        kpv(kc)=kc+1        
        if(kc.eq.1) kmv(kc)=kc   
        if(kc.eq.n3m) kpv(kc)=kc
    4 continue 
      do 3 jc=1,n2m    
        jp=jc+1       
        jmv(jc)=jc-1 
        jpv(jc)=jc+1
        if(jc.eq.1) jmv(jc)=jc    
        if(jc.eq.n2m) jpv(jc)=jc 
    3 continue                  
      do 15 jc=1,n2m           
        jpc(jc)=jpv(jc)-jc    
        jmc(jc)=jc-jmv(jc)   
        jup(jc)=1-jpc(jc)   
        jum(jc)=1-jmc(jc)  
   15 continue            
      return             
      end               
c                      
c***********************************************************************
c
      subroutine meshes
      include 'parade_new.f'
                     
c      dx1=alx1/float(n1m)                                              
      dx1 = (x1c(n1)-x1c(1))/float(n1m)     
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
c
      subroutine input(time,q1,q2,q3,pr,dens,nfilp)
      include 'parade_new.f'

      dimension q1(m1,m2,m3), q2(m1,m2,m3), q3(m1,m2,m3),   
     %                      pr(m1,m2,m3), dens(m1,m2,m3)   
      real*8 qq1(m1_loc,m2,m3),qq2(m1_loc,m2,m3), 
     %     qq3(m1_loc,m2,m3),ppr(m1_loc,m2,m3)
      real*8 rtime

      character*20 namfi3
      character*70 namfi6
      character*5 nfilp,ipfi
      character*2 ipfi1

      write(6,*) ' CHOOSE THE REFERENCE FRAME '
      write(6,*) ' 1) Absolute'
      write(6,*) ' 2) Relative'
      read(5,*)  refsys

      write(6,*) ' ENTER THE NUMBER OF INPUT FILES '
      read(5,*) numprocs
      n1m_loc = int(n1m/numprocs)
      n1_loc = n1m_loc + 1

      write(6,*) ' FILE LIST: '
      write(6,*) ' ---------------------------------'
      call system("ls *.dat")
      write(6,*) ' ---------------------------------'
      write(6,*) ' ENTER the root of the INPUT FILE NAME '
      read(5,1000) namfi3
      write(6,1000) namfi3
 1000 format(a20)

      do imyid = 0,numprocs-1
         time1 = etime(t)

         write(ipfi1,98)imyid
 98      format(i2.2)
         namfi6 = namfi3//'-'//ipfi1//'.dat'
         write(6,201) namfi6
 201     format(10x,'leggo da ',a20)
         read(*,*)yes
         n3pp=(n3-1)/n3p+1           
         n2pp=(n2-1)/n2p+1          
         n1pp=(n1_loc-1)/n1p+1         
         open(59,file=namfi6,form='unformatted',status='unknown')
         read(59) n1pp,n2pp,n3pp               
         write(6,*) ' GRID: ',n1pp,'x',n2pp,'x',n3pp
         read(59) a,a,a,rtime
         time = rtime
         write(6,*) ' TIME: ',time
         read(59) 
     %        (((qq1(i,j,k),i=1,n1_loc,n1p),j=1,n2,n2p),k=1,n3,n3p), 
     %        (((qq2(i,j,k),i=1,n1_loc,n1p),j=1,n2,n2p),k=1,n3,n3p),
     %        (((qq3(i,j,k),i=1,n1_loc,n1p),j=1,n2,n2p),k=1,n3,n3p),
     %        (((ppr(i,j,k),i=1,n1_loc,n1p),j=1,n2,n2p),k=1,n3,n3p)
         time2 = etime(t)
         WRITE(*,*)' LETTURA CONCLUSA ',time2-time1
         do k = 1,n3,n3p    
            do j = 1,n2,n2p
               do i = 1,n1_loc,n1p
                  iglob = n1m_loc*imyid + i
                  q1(iglob,j,k)=qq1(i,j,k)
                  q2(iglob,j,k)=qq2(i,j,k)
                  q3(iglob,j,k)=qq3(i,j,k)
                  pr(iglob,j,k)=ppr(i,j,k)
               enddo
            enddo
         enddo

         time3 = etime(t)
         WRITE(*,*)' RISCRITTO IN SINGOLA PRECISIONE ',
     $     time3-time2
c     read(59)qq1,qq2,qq3,ppr
         close(59)
      enddo
      time2 = etime(t)
      WRITE(*,*)' LETTURA CONCLUSA ',time2-time1
c      do k = 1,n3,n3p    
c         do j = 1,n2,n2p
c            do i = 1,n1,n1p
c               q1(i,j,k)=qq1(i,j,k)
c               q2(i,j,k)=qq2(i,j,k)
c               q3(i,j,k)=qq3(i,j,k)
c               pr(i,j,k)=ppr(i,j,k)
c            enddo
c         enddo
c      enddo
c      q1 = qq1
c      q2 = qq2
c      q3 = qq3
c      pr = ppr

c      do j=1,n2
c         do k=1,n3
cc            q1(0,j,k)=q1(1,j,k)
cc            q2(0,j,k)=q2(1,j,k)
cc            q3(0,j,k)=q3(1,j,k)
cc            pr(0,j,k)=pr(1,j,k)
c            q1(n1+1,j,k)=q1(n1,j,k)   ! a che serve?
c            q2(n1+1,j,k)=q2(n1,j,k)
c            q3(n1+1,j,k)=q3(n1,j,k)
c            pr(n1+1,j,k)=pr(n1,j,k)
c         enddo
c      enddo

      itime=nint(1000.*time)
      write(ipfi,99) itime
   99 format(i5.5)
      namfi6='pressioneP'//ipfi//'.dat'
      open(99,file=namfi6,form='unformatted')
        rewind 99
      write(99) n3,n2,n1,1
      write(99)(((pr(i,j,k),k=1,n3),j=1,n2),i=1,n1)
      close(99)

c      n1mh=(n1-1)/2
c      open(69,file='middleplane.tec',form='formatted',status='unknown')
c      rewind(69)
c      write(69,*) 'VARIABLES = "X","Y","Z","Vx","Vy","Vz","Pres"'
c      write(69,*) 'ZONE I = ',n3, ', J = ',n2, ', K = 1, F=POINT'
c      do 5 j=1,n2
c        do 5 k=1,n3
c          write(69,*)
c     $       x3c(k),x2c(j),x1c(n1mh),
c     $       q3(n1mh,j,k),q2(n1mh,j,k),q1(n1mh,j,k),
c     $       pr(n1mh,j,k)
c 5    continue
c      close(69)

c      n2mh=(n2+1)/2
c      open(70,file='slice1-j84.tec',form='formatted',status='unknown')
c      rewind(70)
c      write(70,*) 'VARIABLES = "X","Y","Z","Vx","Vy","Vz","Pres"'
c      write(70,*) 'ZONE I = ',n3,', J = 1, K = ',n1,', F=POINT'
c      do 6 i=1,n1
c        do 6 k=1,n3
c          write(70,*)
c     $           x3c(k),x2c(n2mh),x1c(i),
c     $           q3(i,n2mh,k),q2(i,n2mh,k),q1(i,n2mh,k),
c     $           pr(i,n2mh,k)
c 6    continue
c      close(70)
      
c     do i=1,n1
c      do j=1,n2
c       do k=1,n3
c        q1(i,j,k)=q1(i,j,k)+u0x+wy*(x3m(k)-x30)-wz*(x2m(j)-x20)
c        q2(i,j,k)=q2(i,j,k)+u0y+wz*(x1m(i)-x10)-wx*(x3m(k)-x30)
c        q3(i,j,k)=q3(i,j,k)+u0z+wx*(x2m(j)-x20)-wy*(x1m(i)-x10)
c       enddo
c      enddo
c     enddo 
      go to 800
      do n=1,npunt
        i=indgeo(1,n,1)
        j=indgeo(1,n,2)
        k=indgeo(1,n,3)
        q1(i,j,k)=0.
c        q1(i,j,k)=u0x+wy*(x3m(k)-x30)-wz*(x2m(j)-x20)
        pr(i,j,k)=0.
      end do
      do n=1,npunr
        i=indgeo(2,n,1)
        j=indgeo(2,n,2)
        k=indgeo(2,n,3)
        q2(i,j,k)=0.
c        q2(i,j,k)=u0y+wz*(x1m(i)-x10)-wx*(x3m(k)-x30)
        pr(i,j,k)=0.
      end do
      do n=1,npunz
        i=indgeo(3,n,1)
        j=indgeo(3,n,2)
        k=indgeo(3,n,3)
        q3(i,j,k)=0.
c        q3(i,j,k)=u0z+wx*(x2m(j)-x20)-wy*(x1m(i)-x10)
        pr(i,j,k)=0.
      end do
c 800  continue
      q1max = -10
      q2max = -10
      q3max = -10
      prmax = -10
      demax = -10
      q1max_abs = -10
      q2max_abs = -10
      q3max_abs = -10
      prmax_abs = -10
      demax_abs = -10
      q1min = -10
      q2min = -10
      q3min = -10
      prmin = -10
      demin = -10
      do k=1,n3
         do j=1,n2
            do i=1,n1
               write(99,*)i,j,k
               if (i.eq.15.and.j.eq.60.and.k.eq.1)write(99,*)q1(i,j,k),
     $              q2(i,j,k),q3(i,j,k),pr(i,j,k),dens(i,j,k)
               if (i.eq.15.and.j.eq.60.and.k.eq.1)write(99,*)q1max,
     $              q2max,q3max,prmax,demax,q1max_abs,q2max_abs,
     $              q3max_abs,prmax_abs,demax_abs
               if (i.eq.15.and.j.eq.60.and.k.eq.1)write(99,*)q1min
     $              ,q2min,q3min
c     $              prmin,demin
               q1max_abs = max(q1max_abs,abs(q1(i,j,k)))
               q2max_abs = max(q2max_abs,abs(q2(i,j,k)))
               q3max_abs = max(q3max_abs,abs(q3(i,j,k)))
               prmax_abs = max(prmax_abs,abs(pr(i,j,k)))
               demax_abs = max(demax_abs,abs(dens(i,j,k)))
               q1max = max(q1max,(q1(i,j,k)))
               q2max = max(q2max,(q2(i,j,k)))
               q3max = max(q3max,(q3(i,j,k)))
               prmax = max(prmax,(pr(i,j,k)))
               demax = max(demax,(dens(i,j,k)))
               q1min = min(q1min,(q1(i,j,k)))
               q2min = min(q2min,(q2(i,j,k)))
               q3min = min(q3min,(q3(i,j,k)))
c               prmin = min(prmin,(pr(i,j,k)))
               demin = min(demin,(dens(i,j,k)))
            end do
         end do
      end do

      write(6,*) ' Valori assoluti massimi di input:'
      write(6,*) ' Q_1 = ',q1max_abs,' Q_2 = ',q2max_abs,
     $     ' Q_3 = ',q3max_abs
      write(6,*) ' Pr  = ',prmax_abs,' Rho = ',demax_abs
      write(6,*) ' Valori massimi di input:'
      write(6,*) ' Q_1 = ',q1max,' Q_2 = ',q2max,' Q_3 = ',q3max
      write(6,*) ' Pr  = ',prmax,' Rho = ',demax
      write(6,*) ' Valori minimi di input:'
      write(6,*) ' Q_1 = ',q1min,' Q_2 = ',q2min,' Q_3 = ',q3min
      write(6,*) ' Pr  = ',prmin,' Rho = ',demin

 234  format(i6,4(2x,e12.5))
 800  continue
      return  
      end
c                                                                       
c***********************************************************************
c
      subroutine outpf(time,vo1,vo2,vo3,voc1,voc2,voc3,prc,densc,
     $     vora,eig1)
      include 'parade_new.f'

      dimension vo1(m1,m2,m3),vo2(m1,m2,m3),vo3(m1,m2,m3)
      dimension eig1(m1,m2,m3),eig2(m1,m2,m3)  !laura arrows
      dimension vo321(m1,m2,m3)   !Laura
      dimension voc1(m1,m2,m3),voc2(m1,m2,m3),voc3(m1,m2,m3)
     %          ,voc31(m1,m2,m3),voc321(m1,m2,m3)
      dimension voc1ass(m1,m2,m3),voc2ass(m1,m2,m3),voc3ass(m1,m2,m3)
      dimension prc(m1,m2,m3),densc(m1,m2,m3),vora(m1,m2,m3)
      dimension vocmax(5),iocm(5),jocm(5),kocm(5) 
      dimension vocmin(5),iocl(5),jocl(5),kocl(5) 
      dimension vo1l(m1*6),vo2l(m1*6),vo3l(m1*6),vo4l(m1*6)
      dimension vo5l(m1*6)
      character*70 namfi3
      character*5 nfilp,ipfi

      itime=nint(1000.*time)
      write(ipfi,99)itime
      write(*,*) 'ITIME ',itime
 99   format(i5.5)

c      q1max = -10
c      q2max = -10
c      q3max = -10
c      prmax = -10
c      demax = -10
c      do k=1,n3
c         do j=1,n2
c            do i=1,n1
c               q1max = max(q1max,abs(voc1(i,j,k)))
c               q2max = max(q2max,abs(voc2(i,j,k)))
c               q3max = max(q3max,abs(voc3(i,j,k)))
c               prmax = max(prmax,abs(prc(i,j,k)))
c            end do
c         end do
c      end do
c      write(6,*) ' Valori massimi di input:'
c      write(6,*) ' Q_1 = ',q1max,' Q_2 = ',q2max,' Q_3 = ',q3max
c      write(6,*) ' Pr  = ',prmax,' Rho = ',demax
      
      write(6,*)' 1 se vuoi stampare il campo di velocita'' '
      read(5,*)iprintV
      if(iprintV.eq.1)then
         namfi3='velocita'//ipfi//'.dat'
         write(6,201) namfi3
 201     format(10x,'scrivo su ',a70)
         n3pp=(n3-1)/n3p+1  
         n2pp=(n2-1)/n2p+1 
         open(59,file=namfi3,form='unformatted',status='unknown')
         rewind(59)
         if(refsys.eq.1) then
            write(6,*)' x10,x20,x30 ',x10,x20,x30
            write(6,*)'wx= ',wx,'u0x= ',u0x,'du0x= ',du0x
            write(6,*)'u0z= ',u0z,'du0z= ',du0z
            do i=1,n1
               do j=1,n2
                  do k=1,n3
c     voc1(i,j,k)=voc1(i,j,k)+u0x+ros*
c     $                    (wy*(x3m(k)-x30)-wz*(x2m(j)-x20))
c     voc2(i,j,k)=voc2(i,j,k)+u0y+ros*
c     $                    (wz*(x1m(i)-x10)-wx*(x3m(k)-x30))
c     voc3(i,j,k)=voc3(i,j,k)+u0z+ros*
c     $                    (wx*(x2m(j)-x20)-wy*(x1m(i)-x10))
                     voc1ass(i,j,k) = voc1(i,j,k)
                     voc2ass(i,j,k) = 
c     $                    uhy +                                 ! nuove incognite
     $                    voc2(i,j,k)*cos(alphax)-
     $                    voc3(i,j,k)*sin(alphax) 
c     $                    - wx*(x3m(k)-x30)*cos(alphax) -       ! nuove incognite
c     $                    wx*(x2m(j)-x20)*sin(alphax)           ! nuove incognite
                     voc3ass(i,j,k) = voc2(i,j,k)*sin(alphax) +
     $                    voc3(i,j,k)*cos(alphax)
c     $                    - wx*(x3m(k)-x30)*sin(alphax) +       ! nuove incognite
c     $                    wx*(x2m(j)-x20)*cos(alphax)           ! nuove incognite
                     voc31(i,j,k)=sqrt(voc2ass(i,j,k)**2.
     $                    +voc3ass(i,j,k)**2.)
                     voc321(i,j,k)=sqrt(voc1ass(i,j,k)**2.
     $                    +voc3ass(i,j,k)**2.+voc2ass(i,j,k)**2.)
                  enddo
               enddo
            enddo 
         else
            do i=1,n1
               do j=1,n2
                  do k=1,n3
                     voc1ass(i,j,k) = voc1(i,j,k)
                     voc2ass(i,j,k) = voc2(i,j,k)
                     voc3ass(i,j,k) = voc3(i,j,k)
                     voc31(i,j,k)=sqrt(voc2(i,j,k)**2.+voc3(i,j,k)**2.)
                     voc321(i,j,k)=sqrt(voc1(i,j,k)**2.+voc3(i,j,k)**2.
     %                    +voc2(i,j,k)**2.)
                  enddo
               enddo
            enddo 
         endif
         write(59) n3,n2,n1,4
c     write(59) epsil, lamb, re, time
         write(59) 
     %        (((voc3ass(i,j,k),k=1,n3),j=1,n2),i=1,n1) !Laura
     %        ,(((voc2ass(i,j,k),k=1,n3),j=1,n2),i=1,n1) !Laura
     %        ,(((voc1ass(i,j,k),k=1,n3),j=1,n2),i=1,n1) !Laura
         close(59)
         
         namfi3='streamline'//ipfi//'.dat'
         write(6,201) namfi3
         open(59,file=namfi3,form='formatted',status='unknown')
         n1mh=180
         rewind(59)
         write(59,*) 'VARIABLES = "X","Y","Z","Vx","Vy","Vz","W_z"'
         write(59,*) 'ZONE I = ',n3, ', J = ',n2, ', K = 1, F=POINT'
         do 8 j=1,n2
            do 8 k=1,n3
               write(59,*)( x30 + (x3c(k)-x30)*cos(alphax)
     $              + (x2c(j)-x20)*sin(alphax) ),
     $              ( hy + x20 - (x3c(k)-x30)*sin(alphax)
     $              + (x2c(j)-x20)*cos(alphax) ),
     $              x1c(n1mh)
     %              ,voc3ass(n1mh,j,k),voc2ass(n1mh,j,k),
     $              voc1ass(n1mh,j,k),vo1(n1mh,j,k)-wx             
 8       continue 
         close(59)
         write(6,202) namfi3
 202     format(10x,' scritto su ',a70)

      endif
      if (nvort.eq.1)then 
         do i=1,n1
c     if (x1c(i).eq.0.)then
            do k=1,n3
               do j=1,n2
                  vo1(i,j,k)=vo1(i,j,k)-wx
c     vo1(i,j,k)=sqrt(vo1(i,j,k)**2.+vo2(i,j,k)**2.+vo3(i,j,k)**2.)   !Laura
                  vo321(i,j,k)=sqrt(vo1(i,j,k)**2.+vo2(i,j,k)**2.
     $                 +vo3(i,j,k)**2.) !Laura
c     raggio=(- 0.5 + 
c     $                    ( x2c(j)*x2c(j)+x3c(k)*x3c(k)+x1c(i)*x1c(i) )**0.5)
c     if (raggio.gt.2.5d-2.and.raggio.lt.5.d-2)then
c     if (x2c(j).gt.0.)then
c     write(45,43)x2c(j),x3c(k),vo321(i,j,k),
c     $                          vo1(i,j,k),vo2(i,j,k),vo3(i,j,k),prc(i,j,k)
c     else
c     write(55,43)x2c(j),x3c(k),vo321(i,j,k),
c     $                          vo1(i,j,k),vo2(i,j,k),vo3(i,j,k),prc(i,j,k)
c     endif
c     endif
               enddo
            enddo
c     endif
         enddo 
c     do j=1,n2
c     if (x2c(j).gt.0.and.x2c(j).lt.1.d-2)then
c     do k=1,n3
c     do i=1,n1
c     vo1(i,j,k)=vo1(i,j,k)-wx
c     c     vo1(i,j,k)=sqrt(vo1(i,j,k)**2.+vo2(i,j,k)**2.+vo3(i,j,k)**2.)   !Laura
c     vo321(i,j,k)=sqrt(vo1(i,j,k)**2.+vo2(i,j,k)**2.
c     $                    +vo3(i,j,k)**2.) !Laura
c                     raggio=(- 0.5 + 
c     $          ( x2c(j)*x2c(j)+x3c(k)*x3c(k)+x1c(i)*x1c(i) )**0.5)
c                     if (raggio.gt.2.5d-2.and.raggio.lt.5.d-2)then
c                        if (x1c(i).ge.0.)then
c                           write(42,43)x1c(i),x3c(k),vo321(i,j,k),
c     $              vo1(i,j,k),vo2(i,j,k),vo3(i,j,k),prc(i,j,k)
c                        else
c                           write(52,43)x1c(i),x3c(k),vo321(i,j,k),
c     $              vo1(i,j,k),vo2(i,j,k),vo3(i,j,k),prc(i,j,k)
c                        endif
c                     endif
c                  enddo
c               enddo
c            endif
c         enddo 
c 43      format(7(2x,e12.5))

c         namfi3='vorticita'//ipfi//'.dat'
         namfi3='vorticita.dat'
         write(6,201) namfi3
         n3pp=(n3-1)/n3p+1  
         n2pp=(n2-1)/n2p+1 
         open(59,file=namfi3,form='unformatted',status='unknown')
         rewind(59)
         write(59) n3,n2,n1,4
         write(59) 
     %        (((vo3(i,j,k),k=1,n3),j=1,n2),i=1,n1)
     %        ,(((vo2(i,j,k),k=1,n3),j=1,n2),i=1,n1)  
     %        ,(((vo1(i,j,k),k=1,n3),j=1,n2),i=1,n1) 
     %        ,(((eig1(i,j,k),k=1,n3),j=1,n2),i=1,n1) 
         close(59)
      
         namfi3='intvortic'//ipfi//'.dat'
         write(6,201) namfi3
         open(59,file=namfi3,form='unformatted',status='unknown')
         rewind(59)
         write(59) n3,n2,n1,1
         write(59) 
     %        (((vo321(i,j,k),k=1,n3),j=1,n2),i=1,n1)
         close(59)
         
      endif
      goto 1979

      n1mh=6
      open(69,file='slice1.tec',form='formatted',status='unknown')
      rewind(69)
      write(69,*) 'VARIABLES = "X","Y","Z","Vx","Vy","Vz","Pres"'
      write(69,*) 'ZONE I = ',n3, ', J = ',n2, ', K = 1, F=POINT'
      do 5 j=1,n2
         do 5 k=1,n3
            write(69,*) x3c(k),x2c(j),x1c(n1mh)
     %           ,voc3(n1mh,j,k),voc2(n1mh,j,k),voc1(n1mh,j,k),
     $           prc(n1mh,j,k)       
 5    continue 
      close(69)
      
      n2mh=33
      open(70,file='slice2.tec',form='formatted',status='unknown')
      rewind(70)
      write(70,*) 'VARIABLES = "X","Y","Z","Vx","Vy","Vz","Pres"'
      write(70,*) 'ZONE I = ',n3, ', J = 1, K = ',n1, ', F=POINT'
      do 6 i=1,n1
         do 6 k=1,n3
            write(70,*) x3c(k),x2c(n2mh),x1c(i)
     %           ,voc3(i,n2mh,k),voc2(i,n2mh,k),voc1(i,n2mh,k),
     $           prc(i,n2mh,k)     
 6    continue 
      close(70)
      
      n3mh=100
      open(71,file='slice3.tec',form='formatted',status='unknown')
      rewind(71)
      write(71,*) 'VARIABLES = "X","Y","Z","Vx","Vy","Vz","Pres"'
      write(71,*) 'ZONE I = 1, J = ',n2, ', K = ',n1, ', F=POINT'
      do 7 i=1,n1
         do 7 j=1,n2
            write(71,*) x3c(n3mh),x2c(j),x1c(i)
     %           ,voc3(i,j,n3mh),voc2(i,j,n3mh),voc1(i,j,n3mh),
     $           prc(i,j,n3mh)     
 7    continue 
      close(71)
      
c      open(72,file='prof_vel.tec',form='formatted',status='unknown')
c      write(72,*) 'Variables = "X" , "Y" '
c      write(72,*) 'ZONE I = ',99,',   F=POINT'
c      do n=2,100
c         write(72,*) voc3(5,n,5),x2c(n)
c      end do
c      close(72)
 1979 continue

      close(42)
      close(45)
      close(52)
      close(55)

      return
      end
c     
c***********************************************************************
c
      subroutine pricor(dens,time)
      include 'parade_new.f'

      dimension yp1(m1,m2),yp2(m1,m2)
      dimension dens(m1,m2,m3)

      character*16 namfile
      character*5 ipfi
      dimension x2cfix(m2,m3),x3cfix(m2,m3)

      if(refsys.eq.1)then
         itime=nint(time*1000.)
         write(ipfi,82)itime
 82      format(i5.5)   
         namfile='bou'//ipfi//'.dat'
         write(*,*)' pricor - hy',hy
         do j=1,n2
            do k=1,n3                !Laura
c     x2c(j)=x2c(j)+hy               !Laura
               x3cfix(j,k) = x30 + (x3c(k)-x30)*cos(alphax)
     $              + (x2c(j)-x20)*sin(alphax) !Laura
               x2cfix(j,k) = hy + x20 - (x3c(k)-x30)*sin(alphax)
     $              + (x2c(j)-x20)*cos(alphax) !Laura               
c     write(6,*) j,x2c(j)
            end do !Laura
         end do
      else
         namfile='bou'//'co'//'.dat'
      endif

      n3pp=(n3-1)/n3p+1   
      n2pp=(n2-1)/n2p+1   
      n1pp=(n1-1)/n1p+1   
      open(18,file=namfile,form='unformatted',status='unknown')
      rewind(18)
      write(18) n3pp,n2pp,n1pp 
      if(refsys.eq.1)then     !Laura
         write(18)                !Laura
     %        (((x3cfix(j,k),k=1,n3,n3p),j=1,n2,n2p),i=1,n1,n1p),    !Laura
     %        (((x2cfix(j,k),k=1,n3,n3p),j=1,n2,n2p),i=1,n1,n1p),    !Laura  
     %        (((x1c(i),k=1,n3,n3p),j=1,n2,n2p),i=1,n1,n1p)        !Laura
      else                        !Laura
         write(18) 
     %        (((x3c(k),k=1,n3,n3p),j=1,n2,n2p),i=1,n1,n1p),
     %        (((x2c(j),k=1,n3,n3p),j=1,n2,n2p),i=1,n1,n1p),
     %        (((x1c(i),k=1,n3,n3p),j=1,n2,n2p),i=1,n1,n1p) 
c     % ,(((nint(1.-dens(i,j,k)),i=1,n1,n1p),j=1,n2,n2p),k=1,n3,n3p) 
      endif                        !Laura
      close(18)
      return
      end
c 
c***********************************************************************
c
      subroutine svelt(q1,q2,q3,pr,dens,voc1,voc2,
     $     voc3,prc,densc,time)  
      include 'parade_new.f'

      dimension q1(m1,m2,m3), q2(m1,m2,m3), q3(m1,m2,m3),   
     %     pr(m1,m2,m3), dens(m1,m2,m3)
      dimension voc1(m1,m2,m3),voc2(m1,m2,m3),voc3(m1,m2,m3),
     %     prc(m1,m2,m3),densc(m1,m2,m3),vmax(4)
      dimension isym(m1)
c     
      dimension indint(5,mpun,3) ! x punti interni
      dimension q3c3(m1,m2,m3),q3c2(m1,m2,m3),q3c1(m1,m2,m3)
      dimension q2c3(m1,m2,m3),q2c2(m1,m2,m3),q2c1(m1,m2,m3)
      dimension q1c3(m1,m2,m3),q1c2(m1,m2,m3),q1c1(m1,m2,m3)
      dimension prc3(m1,m2,m3),prc2(m1,m2,m3),prc1(m1,m2,m3)
c     dimension blankv(3000000),blankp(3000000)
      dimension vert(45000,3),var(45000),ielem(15000,3)
      dimension tauw(55000),tauw1(55000),tauw2(55000),tauw3(55000)
     &     ,prsw(55000)
      dimension vett32(3),vett21(3),vettnorm(3),vettrif(3),versnorm(3)
      dimension defvel(3,3)
c
      character*50 namefile     !LAURA
      character*16 namfile
      character*5 ipfi
      
      real*8 norb(17000,3)      !LAURA
      real*8 xyzb(17000,9)      !LAURA
      integer nb                !LAURA
      real sup_tot,pitch        !LAURA
      
      itime=nint(time*1000.)
      write(ipfi,82)itime
 82   format(i5.5)
      goto 2004
      do k=1,n3
         psirv(1,k) = 0.
         do j=2,n2
            psirv(j,k) = psirv(j-1,k) + q3(1,j-1,k)/dx2*g2m(j-1)
         enddo
      enddo
c     do j=1,n2
c     psirv(j,n3) = psirv(j,1)
c     enddo
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       do k=1,n3m
        do j=1,n2m
         do i=1,n1m
          q3c3(i,j,k)=x3c(k)
          q3c2(i,j,k)=(x2c(j)+x2c(j+1))/2.0
          q3c1(i,j,k)=(x1c(i)+x1c(i+1))/2.0 
c
          q2c3(i,j,k)=(x3c(k)+x3c(k+1))/2.0
          q2c2(i,j,k)=x2c(j)
          q2c1(i,j,k)=(x1c(i)+x1c(i+1))/2.0
c  
          q1c3(i,j,k)=(x3c(k)+x3c(k+1))/2.0
          q1c2(i,j,k)=(x2c(j)+x2c(j+1))/2.0
          q1c1(i,j,k)=x1c(i)

          prc3(i,j,k)=(x3c(k)+x3c(k+1))/2.0
          prc2(i,j,k)=(x2c(j)+x2c(j+1))/2.0
          prc1(i,j,k)=(x1c(i)+x1c(i+1))/2.0
         end do
        end do
       end do
c
ccccccccccccccccccccccccc
c
       do n=1,npunz_i
        kz=indgeoi(3,n,3)
        jy=indgeoi(3,n,2)
        ix=indgeoi(3,n,1)
        q3(ix,jy,kz)=0.
       end do
c
       do n=1,npunz_e
         indint(3,n,3)=indgeo(3,n,3)*2-indgeoe(3,n,3)
         indint(3,n,2)=indgeo(3,n,2)*2-indgeoe(3,n,2)
         indint(3,n,1)=indgeo(3,n,1)*2-indgeoe(3,n,1)
         ki=indint(3,n,3)
         ji=indint(3,n,2)
         ii=indint(3,n,1)
         ke=indgeo(3,n,3)
         je=indgeo(3,n,2)
         ie=indgeo(3,n,1)
         kee=indgeoe(3,n,3)
         jee=indgeoe(3,n,2)
         iee=indgeoe(3,n,1)
         if(distb(3,n).eq.0.) write(*,*) 'p3',ie,je,ke,q3(ie,je,ke) 
         if(distb(3,n).eq.0.) goto 100
c
         diste_ee=sqrt((q3c3(iee,jee,kee)-q3c3(ie,je,ke))**2+
     %                 (q3c2(iee,jee,kee)-q3c2(ie,je,ke))**2+
     %                 (q3c1(iee,jee,kee)-q3c1(ie,je,ke))**2)
         disti_e=sqrt((q3c3(ii,ji,ki)-q3c3(ie,je,ke))**2+
     %                (q3c2(ii,ji,ki)-q3c2(ie,je,ke))**2+
     %                (q3c1(ii,ji,ki)-q3c1(ie,je,ke))**2)  
         distb_e=distb(3,n)/(1-distb(3,n))*diste_ee
         q3(ii,ji,ki)=-q3(ie,je,ke)*distb_e/(disti_e-distb_e)
       if (abs(q3(ii,ji,ki)).gt.2.) q3(ii,ji,ki)=0.  
 100   end do
c
ccccccccccccccccccccccccc
c
       do n=1,npunr_i
        kz=indgeoi(2,n,3)
        jy=indgeoi(2,n,2)
        ix=indgeoi(2,n,1)
        q2(ix,jy,kz)=0.
       end do
c   
       do n=1,npunr_e

         indint(2,n,3)=indgeo(2,n,3)*2-indgeoe(2,n,3)
         indint(2,n,2)=indgeo(2,n,2)*2-indgeoe(2,n,2)
         indint(2,n,1)=indgeo(2,n,1)*2-indgeoe(2,n,1)
         ki=indint(2,n,3)
         ji=indint(2,n,2)
         ii=indint(2,n,1)
         ke=indgeo(2,n,3)
         je=indgeo(2,n,2)
         ie=indgeo(2,n,1)
         kee=indgeoe(2,n,3)
         jee=indgeoe(2,n,2)
         iee=indgeoe(2,n,1)
         if(distb(2,n).eq.0.) write(*,*) 'p2',ie,je,ke,q2(ie,je,ke) 
         if(distb(2,n).eq.0.) goto 200
c
         diste_ee=sqrt((q2c3(iee,jee,kee)-q2c3(ie,je,ke))**2+
     %                 (q2c2(iee,jee,kee)-q2c2(ie,je,ke))**2+
     %                 (q2c1(iee,jee,kee)-q2c1(ie,je,ke))**2)
         disti_e=sqrt((q2c3(ii,ji,ki)-q2c3(ie,je,ke))**2+
     %                (q2c2(ii,ji,ki)-q2c2(ie,je,ke))**2+
     %                (q2c1(ii,ji,ki)-q2c1(ie,je,ke))**2)
         distb_e=distb(2,n)/(1-distb(2,n))*diste_ee
         q2(ii,ji,ki)=-q2(ie,je,ke)*distb_e/(disti_e-distb_e)
       if (abs(q2(ii,ji,ki)).gt.2.) q2(ii,ji,ki)=0.  
 200   end do
c   
ccccccccccccccccccccccccc
c
       do n=1,npunt_i
        kz=indgeoi(1,n,3)
        jy=indgeoi(1,n,2)
        ix=indgeoi(1,n,1)
        q1(ix,jy,kz)=0.
       end do
c
       do n=1,npunt_e
         indint(1,n,3)=indgeo(1,n,3)*2-indgeoe(1,n,3)
         indint(1,n,2)=indgeo(1,n,2)*2-indgeoe(1,n,2)
         indint(1,n,1)=indgeo(1,n,1)*2-indgeoe(1,n,1)
         ki=indint(1,n,3)
         ji=indint(1,n,2)
         ii=indint(1,n,1)
         ke=indgeo(1,n,3)
         je=indgeo(1,n,2)
         ie=indgeo(1,n,1)
         kee=indgeoe(1,n,3)
         jee=indgeoe(1,n,2)
         iee=indgeoe(1,n,1)
         if(distb(1,n).eq.0.) write(*,*) 'p1',ie,je,ke,q1(ie,je,ke) 
         if(distb(1,n).eq.0.) goto 300
c
         diste_ee=sqrt((q1c3(iee,jee,kee)-q1c3(ie,je,ke))**2+
     %                 (q1c2(iee,jee,kee)-q1c2(ie,je,ke))**2+
     %                 (q1c1(iee,jee,kee)-q1c1(ie,je,ke))**2)
         disti_e=sqrt((q1c3(ii,ji,ki)-q1c3(ie,je,ke))**2+
     %                (q1c2(ii,ji,ki)-q1c2(ie,je,ke))**2+
     %                (q1c1(ii,ji,ki)-q1c1(ie,je,ke))**2)
         distb_e=distb(1,n)/(1-distb(1,n))*diste_ee
         q1(ii,ji,ki)=-q1(ie,je,ke)*distb_e/(disti_e-distb_e)
       if (abs(q1(ii,ji,ki)).gt.2.) q1(ii,ji,ki)=0.  
 300   end do
c  
ccccccccccccccccccccc
c
       do n=1,npunp_i
        kz=indgeoi(4,n,3)
        jy=indgeoi(4,n,2)       
        ix=indgeoi(4,n,1)
c        ip=ix + jy*(n1-1) + kz*(n1-1)*(n2-1)
c        blankp(ip)=1.
        pr(ix,jy,kz)=0.
       end do
c
       kkki=0
       kkke=0
       kkkee=0
       kkkeee=0
       kkkint=0
       kkkest=0
       do n=1,npunp_e
         ki=indgeo(4,n,3)*2-indgeoe(4,n,3)
         ji=indgeo(4,n,2)*2-indgeoe(4,n,2)
         ii=indgeo(4,n,1)*2-indgeoe(4,n,1)
         ke=indgeo(4,n,3)
         je=indgeo(4,n,2)
         ie=indgeo(4,n,1)
         kee=indgeoe(4,n,3)
         jee=indgeoe(4,n,2)
         iee=indgeoe(4,n,1)
         keee=indgeoe(4,n,3)*2-indgeo(4,n,3)
         jeee=indgeoe(4,n,2)*2-indgeo(4,n,2)
         ieee=indgeoe(4,n,1)*2-indgeo(4,n,1)
         if(abs(pr(ii,ji,ki)).gt.1.) kkki=kkki+1
c         if(abs(pr(ii,ji,ki)).gt.1.) write(*,*) n,ii,ji,ki,ie,je,ke
         if(abs(pr(ie,je,ke)).gt.1.) kkke=kkke+1
c         if(abs(pr(ie,je,ke)).gt.1.) write(*,*) n,ie,je,ke
         if(abs(pr(iee,jee,kee)).gt.1.) kkkee=kkkee+1
c         if(abs(pr(iee,jee,kee)).gt.1.) write(*,*)'ee',n,iee,jee,kee
c         if(abs(pr(iee,jee,kee)).gt.1.) pr(iee,jee,kee)=0.
         if(abs(pr(ieee,jeee,keee)).gt.1.) kkkeee=kkkeee+1
c      if(abs(pr(ieee,jeee,keee)).gt.1.) write(*,*)'eee',n,ieee,jeee,keee
c         if(abs(pr(ieee,jeee,keee)).gt.1.) pr(ieee,jeee,keee)=0.
c
         distee_eee=sqrt((prc3(iee,jee,kee)-prc3(ieee,jeee,keee))**2+
     %                   (prc2(iee,jee,kee)-prc2(ieee,jeee,keee))**2+
     %                   (prc1(iee,jee,kee)-prc1(ieee,jeee,keee))**2)
         diste_ee=sqrt((prc3(iee,jee,kee)-prc3(ie,je,ke))**2+
     %                 (prc2(iee,jee,kee)-prc2(ie,je,ke))**2+
     %                 (prc1(iee,jee,kee)-prc1(ie,je,ke))**2)
         disti_e=sqrt((prc3(ii,ji,ki)-prc3(ie,je,ke))**2+
     %                (prc2(ii,ji,ki)-prc2(ie,je,ke))**2+
     %                (prc1(ii,ji,ki)-prc1(ie,je,ke))**2)
c      pr(ie,je,ke)=pr(iee,jee,kee)+(pr(iee,jee,kee)-pr(ieee,jeee,keee))*
c     %                diste_ee/distee_eee 
c      pr(ii,ji,ki)=pr(ie,je,ke)+(pr(ie,je,ke)-pr(iee,jee,kee))*
c     %                disti_e/diste_ee
c         if(abs(pr(ie,je,ke)).gt.1.) kkkest=kkkest+1  
c         if(abs(pr(ie,je,ke)).gt.1.) write(*,*) n,ie,je,ke
c         if(abs(pr(ii,ji,ki)).gt.1.) kkkint=kkkint+1
c         if(abs(pr(ii,ji,ki)).gt.1.) write(*,*) n,ii,ji,ki
       end do 
 
       write(*,*) 'pres_i>1........',kkki
       write(*,*) 'pres_e>1........',kkke
       write(*,*) 'pres_ee>1.......',kkkee
       write(*,*) 'pres_eee>1......',kkkeee
c       write(*,*) 'pres_est>1......',kkkest
c       write(*,*) 'pres_int>1......',kkkint
c

 2004 continue
C     parte di lettura del file surface.dat      
c      open(30,file='surface.dat',form='formatted',status='old')
c      read(30,*) nv,ne
c      read(30,*) vert(1,3),vert(1,2),vert(1,1)
c      vert3max=vert(1,3)
c      vert2max=vert(1,2)
c      vert1max=vert(1,1)
c      vert3min=vert(1,3)
c      vert2min=vert(1,2)
c      vert1min=vert(1,1)
c      do l=2,nv
c        read(30,*) vert(l,3),vert(l,2),vert(l,1)
c        vert3max=max(vert(l,3),vert3max)
c        vert2max=max(vert(l,2),vert2max)
c        vert1max=max(vert(l,1),vert1max)
c        vert3min=min(vert(l,3),vert3min)
c        vert2min=min(vert(l,2),vert2min)
c        vert1min=min(vert(l,1),vert1min)
c      end do

c      do l=1,ne
c        read(30,*) ielem(l,3),ielem(l,2),ielem(l,1)
c      enddo 
c      close(30)
c      write(*,*)' surface.dat file with',ne,'elements'

c     call maxmin(vert3,nv,vert3max,vert3min)
c     call maxmin(vert2,nv,vert2max,vert2min)
c     call maxmin(vert1,nv,vert1max,vert1min)
c      refp3=vert3max-vert3min
c      refp2=vert2max-vert2min
c      refp1=vert1max-vert1min
c      write(6,*)'refp3= ',refp3
c      write(6,*)'refp2= ',refp2
c      write(6,*)'refp1= ',refp1
c     do l=1,ne
c      do n=1,3
c      ntmp=ntmp+1
c      ielem(l,n) = ntmp
c       read(30,*) ielem(l,3),ielem(l,2),ielem(l,1)
c      enddo 
c     end do

c      do l=1,nv  
c        do k=1,n3
c          if(vert(l,3).gt.x3c(k).and.vert(l,3).le.x3c(k+1)) then
c            do j=1,n2
c              if(vert(l,2).gt.x2c(j).and.vert(l,2).le.x2c(j+1)) then
c                do i=1,n1
c                  if(vert(l,1).gt.x1c(i).and.vert(l,1).le.x1c(i+1))
c     %            var(l)=pr(i,j,k)
c                end do
c              end if
c            end do
c          end if       
c        end do
c      end do 

c     if(refsys.eq.1) then
c     namfile='surface'//ipfi//'.tec'
c     else
c     namfile='surface.tec'
c     endif
c     open(31,file=namfile,form='formatted',status='unknown')
c     write(31,*) 'Variables = X, Y, Z, Pres'
c     write(31,*) 'ZONE N = ',nv,',E = ',ne,',F=FEPOINT, ET=TRIANGLE'
c     write(31,*) 'DT=(SINGLE SINGLE SINGLE SINGLE)'
c     if(refsys.eq.1) then
c     do l=1,nv
c       write(31,*) vert(l,3),vert(l,2)+hy,vert(l,1)+hx,var(l)
c     end do
c     else
c     do l=1,nv
c       write(31,*) vert(l,3),vert(l,2),vert(l,1),var(l)
c     end do
c     endif
c     do l=1,ne
c       write(31,*) ielem(l,3),ielem(l,2),ielem(l,1)
c     end do
c     close(31)


c
c     goto 1222	
      
      risul3_v=0.
      risul2_v=0.
      risul1_v=0.
      risul3_p=0.
      risul2_p=0.
      risul1_p=0.
      risulM_x=0.
      sup_tot=0.0
      
      namefile = 'surf3n.stl'
      write(*,*)'read the geometry from file ',namefile
      pitch = 0.
      call readgeo3ddim(nb,namefile,pitch)
      call readgeo3d(time,xyzb,norb,nb,
     $     x1,x2,y1,y2,z1,z2,namefile,pitch)
      write(*,*)' triangles which describe the surface = ',nb
      do itr=1,nb
c     if (norb(itr,3).gt.0.5)then 
         write(32,*)' tr',itr,
     $        'normale x,y,z = ',norb(itr,1),norb(itr,2),norb(itr,3)
c               goto 987
c     endif
         v13 = xyzb(itr,1)
         v12 = xyzb(itr,2)
         v11 = xyzb(itr,3)
         v23 = xyzb(itr,4)
         v22 = xyzb(itr,5)
         v21 = xyzb(itr,6)
         v33 = xyzb(itr,7)
         v32 = xyzb(itr,8)
         v31 = xyzb(itr,9)
c     v33=vert(ielem(itr,3),3)      
c     v32=vert(ielem(itr,3),2)
c     v31=vert(ielem(itr,3),1)
c     v23=vert(ielem(itr,2),3)
c     v22=vert(ielem(itr,2),2)
c     v21=vert(ielem(itr,2),1)
c     v13=vert(ielem(itr,1),3)
c     v12=vert(ielem(itr,1),2)
c     v11=vert(ielem(itr,1),1)
         d32=sqrt((v33-v23)**2+(v32-v22)**2+(v31-v21)**2)
         d21=sqrt((v23-v13)**2+(v22-v12)**2+(v21-v11)**2)
         d13=sqrt((v13-v33)**2+(v12-v32)**2+(v11-v31)**2)
         sp = (d32+d21+d13)/2.
         area = sqrt(sp*(sp-d32)*(sp-d21)*(sp-d13))
         sup_tot = sup_tot + area
         bar3=(v33+v23+v13)/3.
         bar2=(v32+v22+v12)/3.
         bar1=(v31+v21+v11)/3.
         
c     calcolo del vettore normale dalle coordinate vertici
C     vett32(3)=v33-v23
C     vett32(2)=v32-v22
C     vett32(1)=v31-v21
C     vett21(3)=v23-v13
C     vett21(2)=v22-v12
C     vett21(1)=v21-v11
C     c     calcolo analitico del vettore normale al triangolo       
C     vettnorm(3)=  vett32(1)*vett21(2)-vett21(1)*vett32(2)
C     vettnorm(2)=-(vett32(1)*vett21(3)-vett21(1)*vett32(3))
C     vettnorm(1)=  vett32(2)*vett21(3)-vett21(2)*vett32(3)
C     c      if (itr.le.5) write(*,*) '3',vettnorm(3)
C     C     if (itr.le.5) write(*,*) '2',vettnorm(2)
C     C     if (itr.le.5) write(*,*) '1',vettnorm(1)
C     c     vettrif(3)=refp3-bar3
C     c     vettrif(2)=refp2-bar2
C     c     vettrif(1)=refp1-bar1  MODIFICA LAURA SFERA
C     vettrif(3)=0.-bar3
C     vettrif(2)=0.-bar2
C     vettrif(1)=0.-bar1
C     prsc=vettnorm(3)*vettrif(3)+vettnorm(2)*vettrif(2)+
C     %           vettnorm(1)*vettrif(1)
C     if(prsc.gt.0.)then
C     vettnorm(3)=-vettnorm(3)
C     vettnorm(2)=-vettnorm(2)
C     vettnorm(1)=-vettnorm(1)
C     endif
         
         vettnorm(3) = norb(itr,1) !!!asse x stl e' asse x3
         vettnorm(2) = norb(itr,2)
         vettnorm(1) = norb(itr,3) !!!asse z stl e' asse x1
         vmod=sqrt(vettnorm(3)**2+vettnorm(2)**2+vettnorm(1)**2)
         if (vmod.le.0.001) then               
            write(*,*)' tr',itr,
     $           'normale x,y,z = ',norb(itr,1),norb(itr,2),norb(itr,3)
            goto 987
         endif
         vettrif(3)=0.2-bar3
         vettrif(2)=0.-bar2
c     vettrif(1)=0.-bar1
         prsc=vettnorm(3)*vettrif(3)+vettnorm(2)*vettrif(2)
c     %           +vettnorm(1)*vettrif(1)
         if(prsc.gt.0.)then
            vettnorm(3)=-vettnorm(3)
            vettnorm(2)=-vettnorm(2)
c     vettnorm(1)=-vettnorm(1)
         endif
         versnorm(3) = vettnorm(3)/vmod
         versnorm(2) = vettnorm(2)/vmod
         versnorm(1) = vettnorm(1)/vmod
C     cfr calcolo normali con le due tecniche 
C     if (abs(versnorm(1)-norb(itr,3)).gt.3.d-3)
C     $           write(*,*)'versnorm(1),norb(itr,3)',itr,versnorm(1),
C     $           norb(itr,3)
C     if (abs(versnorm(2)-norb(itr,2)).gt.3.d-3)
C     $           write(*,*)'versnorm(2),norb(itr,2)',itr,versnorm(2),
C     $           norb(itr,2)
C     if (abs(versnorm(3)-norb(itr,1)).gt.3.d-3)
C     $           write(*,*)'versnorm(3),norb(itr,1)',itr,versnorm(3),
C     $           norb(itr,1)
         vmodvers=sqrt(versnorm(3)**2+versnorm(2)**2+versnorm(1)**2)
         delta=0.
         do k=1,n3m
            
            if(bar3.gt.x3c(k).and.bar3.le.x3c(k+1)) then
               write(98,*)'itr',itr,'if3'
               do j=1,n2m
                  
                  if(bar2.gt.x2c(j).and.bar2.le.x2c(j+1)) then
                     write(98,*)'itr',itr,'if2'
                     do i=1,n1m
                        
                        if(bar1.gt.x1c(i).and.bar1.le.x1c(i+1)) then
                           delta=1.5*sqrt((x1c(i)-x1c(i+1))**2.+
     &                          (x2c(j)-x2c(j+1))**2.+
     &                          (x3c(k)-x3c(k+1))**2.)
                           write(98,*)'itr',itr,'if3',bar3,bar2,bar1
                        endif
                     end do
                  end if
               end do
            end if       
         end do
         if (delta.eq.0.)then
            write(97,*)'itr',itr,'out',bar3,bar2,bar1
            goto 987
         endif

         delta3=delta*versnorm(3)+bar3
         delta2=delta*versnorm(2)+bar2
         delta1=delta*versnorm(1)+bar1
         distdb3=delta3-bar3
         distdb2=delta2-bar2
         distdb1=delta1-bar1
         kint=1
         jint=1
         iint=1
c     se componenti del vettore normale negative
         if(distdb3.lt.0.) kint=-1
         if(distdb2.lt.0.) jint=-1
         if(distdb1.lt.0.) iint=-1
         taumod=0.
         tau1=0.
         tau2=0.
         tau3=0.
         prsint=0.
         forza1_v=0.
         forza2_v=0.
         forza3_v=0.
         forza1_p=0.
         forza2_p=0.
         forza3_p=0.
         do k=4,n3m-1
            if(delta3.gt.x3c(k-1).and.delta3.le.x3c(k)) then
               do j=2,n2m
                  if(delta2.gt.x2c(j-1).and.delta2.le.x2c(j)) then
c     do i=4,n1m-1
                     do i=2,n1
                        if(delta1.gt.x1c(i-1).and.delta1.le.x1c(i))then
                           km=kmv(k)
                           jm=jmv(j)
                           im=imv(i)
                           ks=km+kint
                           js=jm+jint
                           is=im+iint
                           if ((x1c(i) - x1c(im)).eq.0.)
     $                          write(*,*)k,j,i,x1c(i),x1c(im)
                           if ((x1m(is)-x1m(im)).eq.0.)
     $                          write(*,*)k,j,i,x1m(is),x1m(im)
C     INTERPOLAZIONE PRIMA COMPONENTE 
                           q1int1=(q1(i ,jm,km)  * (delta1-x1c(im))
     &                          -q1(im,jm,km) * (delta1-x1c(i )))
     $                          /(x1c(i) - x1c(im))
c     $                             *dx1
                           
                           q1int2=(q1(i ,js,km)  * (delta1-x1c(im))
     &                          -q1(im,js,km) * (delta1-x1c(i )))
     $                          /(x1c(i) - x1c(im))
c     $                             *dx1
                           
                           q1int3=(q1int2        * (delta2-x2m(jm))
     &                          -q1int1        * (delta2-x2m(js)))/
     &                          (x2m(js)-x2m(jm))
                           
                           q1int1=(q1(i ,jm,ks)  * (delta1-x1c(im))
     &                          -q1(im,jm,ks) * (delta1-x1c(i )))
     $                          /(x1c(i) - x1c(im))
c     $                             *dx1
                           
                           q1int2=(q1(i ,js,ks)  * (delta1-x1c(im))
     &                          -q1(im,js,ks) * (delta1-x1c(i )))
     $                          /(x1c(i) - x1c(im))
C     *dx1
                           q1int4=(q1int1        * (delta2-x2m(js))
     &                          -q1int2        * (delta2-x2m(jm)))/
     &                          (x2m(jm)-x2m(js))
c     &                             (x2m(js)-x2m(jm))
                           
                           q1int =(q1int4        * (delta3-x3m(km))
     &                          -q1int3        * (delta3-x3m(ks)))/
     &                          (x3m(ks)-x3m(km))
c     &                             (x3m(js)-x3m(jm))
                           
C     FINE  INTERPOLAZIONE PRIMA COMPONENTE 
                           
C     INTERPOLAZIONE SECONDA COMPONENTE 
                           q2int1=(q2(im,j ,km)  * (delta2-x2c(jm))
     &                          -q2(im,jm,km)  * (delta2-x2c(j )))/
     &                          (x2c(j)-x2c(jm))
                           
                           q2int2=(q2(im,j ,ks)  * (delta2-x2c(jm))
     &                          -q2(im,jm,ks)  * (delta2-x2c(j )))/
     &                          (x2c(j)-x2c(jm))
                           
                           q2int3=(q2int2        * (delta3-x3m(km))
     &                          -q2int1        * (delta3-x3m(ks)))/
     &                          (x3m(ks)-x3m(km))
                           
                           q2int1=(q2(is,j ,km)  * (delta2-x2c(jm))
     &                          -q2(is,jm,km)  * (delta2-x2c(j )))/
     &                          (x2c(j)-x2c(jm))
                           
c     q2int2=(q2(is,j ,km)  * (delta2-x2c(jm))
c     &                             -q2(is,jm,km)  * (delta2-x2c(j )))/
                           q2int2=(q2(is,j ,ks)  * (delta2-x2c(jm))
     &                          -q2(is,jm,ks)  * (delta2-x2c(j )))/
     &                          (x2c(j)-x2c(jm))
                           
                           q2int4=(q2int1        * (delta3-x3m(ks))
     &                          -q2int2        * (delta3-x3m(km)))/
     &                          (x3m(km)-x3m(ks))
c     &                             (x3m(ks)-x3m(km))
                           
                           q2int =(q2int4        * (delta1-x1m(im))
     &                          -q2int3        * (delta1-x1m(is)))/
     &                          (x1m(is)-x1m(im))
C     FINE  INTERPOLAZIONE SECONDA COMPONENTE 
                           
C     INTERPOLAZIONE TERZA COMPONENTE 
                           q3int1=(q3(im,jm,k )  * (delta3-x3c(km))
     &                          -q3(im,jm,km)  * (delta3-x3c(k )))/
     &                          (x3c(k)-x3c(km))
                           
                           q3int2=(q3(is,jm,k )  * (delta3-x3c(km))
     &                          -q3(is,jm,km)  * (delta3-x3c(k )))/
     &                          (x3c(k)-x3c(km))
                           
                           q3int3=(q3int2        * (delta1-x1m(im))
     &                          -q3int1        * (delta1-x1m(is)))/
     &                          (x1m(is)-x1m(im))
                           
                           q3int1=(q3(im,js,k )  * (delta3-x3c(km))
     &                          -q3(im,js,km)  * (delta3-x3c(k )))/
     &                          (x3c(k)-x3c(km))
                           
c     q3int2=(q3(im,js,k )  * (delta3-x3c(km))
c     &                             -q3(im,js,km)  * (delta3-x3c(k )))/
                           q3int2=(q3(is,js,k )  * (delta3-x3c(km))
     &                          -q3(is,js,km)  * (delta3-x3c(k )))/
     &                          (x3c(k)-x3c(km))
                           
                           q3int4=(q3int1        * (delta1-x1m(is))
     &                          -q3int2        * (delta1-x1m(im)))/
     &                          (x1m(im)-x1m(is))
c     &                             (x1m(is)-x1m(im))
                           
                           q3int =(q3int4        * (delta2-x2m(jm))
     &                          -q3int3        * (delta2-x2m(js)))/
     &                          (x2m(js)-x2m(jm))
C     FINE  INTERPOLAZIONE TERZA COMPONENTE 
                           
C     INTERPOLAZIONE PRESSIONE
                              prsint1=(pr(im,jm,ks)  * (delta3-x3m(km))
     &                             - pr(im,jm,km)  * (delta3-x3m(ks)))/
     &                             (x3m(ks)-x3m(km))
                         
                              prsint2=(pr(is,jm,ks)  * (delta3-x3m(km))
     &                             - pr(is,jm,km)  * (delta3-x3m(ks)))/
     &                             (x3m(ks)-x3m(km))
                              
                              prsint3=(prsint2        * (delta1-x1m(im))
     &                             - prsint1        * (delta1-x1m(is)))/
     &                             (x1m(is)-x1m(im))
                              
                              prsint1=(pr(im,js,ks)  * (delta3-x3m(km))
     &                             - pr(im,js,km)  * (delta3-x3m(ks)))/
     &                             (x3m(ks)-x3m(km))
                              
c                              prsint2=(pr(im,js,ks)  * (delta3-x3m(km))
c     &                             - pr(im,js,km)  * (delta3-x3m(ks)))/
                              prsint2=(pr(is,js,ks)  * (delta3-x3m(km))
     &                             - pr(is,js,km)  * (delta3-x3m(ks)))/
     &                             (x3m(ks)-x3m(km))
                              
                              prsint4=(prsint1        * (delta1-x1m(is))
     &                             - prsint2        * (delta1-x1m(im)))/
     &                             (x1m(im)-x1m(is))
c     &                             (x1m(is)-x1m(im))
                              
                              prsint =(prsint4        * (delta2-x2m(jm))
     &                             -prsint3        * (delta2-x2m(js)))/
     &                             (x2m(js)-x2m(jm))
C     FINE  INTERPOLAZIONE PRESSIONE 
C     ESTRAPOLAZIONE PRESSIONE A BAR1,BAR2,BAR3
c                           prsint1=(pr(im,jm,ks)  * (bar3-x3m(km))
c     &                          - pr(im,jm,km)  * (bar3-x3m(ks)))/
c     &                          (x3m(ks)-x3m(km))
c                           
c                           prsint2=(pr(is,jm,ks)  * (bar3-x3m(km))
c     &                          - pr(is,jm,km)  * (bar3-x3m(ks)))/
c     &                          (x3m(ks)-x3m(km))
c                           
c                           prsint3=(prsint2        * (bar1-x1m(im))
c     &                          - prsint1        * (bar1-x1m(is)))/
c     &                          (x1m(is)-x1m(im))
c                           
c                           prsint1=(pr(im,js,ks)  * (bar3-x3m(km))
c     &                          - pr(im,js,km)  * (bar3-x3m(ks)))/
c     &                          (x3m(ks)-x3m(km))
c                          
cc     prsint2=(pr(im,js,ks)  * (delta3-x3m(km))
cc     &                             - pr(im,js,km)  * (delta3-x3m(ks)))/
c                           prsint2=(pr(is,js,ks)  * (bar3-x3m(km))
c     &                          - pr(is,js,km)  * (bar3-x3m(ks)))/
c     &                          (x3m(ks)-x3m(km))
c                           
c                           prsint4=(prsint1        * (bar1-x1m(is))
c     &                          - prsint2        * (bar1-x1m(im)))/
c     &                          (x1m(im)-x1m(is))
cc     &                             (x1m(is)-x1m(im))
c                           
c                           prsint =(prsint4        * (bar2-x2m(jm))
c     &                          -prsint3        * (bar2-x2m(js)))/
c     &                          (x2m(js)-x2m(jm))
C     FINE  ESTRAPOLAZIONE PRESSIONE 
                           gradpr = duhy*cos(alphax)*versnorm(2)-
     $                             duhy*sin(alphax)*versnorm(3)
     $           +dwx*(-versnorm(2)*(bar3-x30)+versnorm(3)*(bar2-x20))
     $           -wx*wx*(versnorm(2)*(bar2-x20)+versnorm(3)*(bar3-x30))
                           prsurf = prsint + gradpr*delta
C     PROIEZIONE DI QINT LUNGO N     
                           qintn = q1int*versnorm(1)+q2int*versnorm(2)+
     &                          q3int*versnorm(3)
C     VETTORE TANGENTE ALLA SUPERFICIE
                           qtang1 = q1int-qintn*versnorm(1)
                           qtang2 = q2int-qintn*versnorm(2)
                           qtang3 = q3int-qintn*versnorm(3)                           
C     velocita' del corpo
                           vsurf1 = 0.
                           vsurf2 = uhy*cos(alphax)-wx*(delta3-x30)
                           vsurf3 = - uhy*sin(alphax)+wx*(delta2-x20)
                           vsurfn = vsurf1*versnorm(1)
     $                          +vsurf2*versnorm(2)+vsurf3*versnorm(3)
                           vsurf1_tang = vsurf1 - vsurfn*versnorm(1)
                           vsurf2_tang = vsurf2 - vsurfn*versnorm(2)
                           vsurf3_tang = vsurf3 - vsurfn*versnorm(3)
C     SFORZO TANGENZIALE
                           tau1=(qtang1 - vsurf1_tang)/delta/ren
                           tau2=(qtang2 - vsurf2_tang)/delta/ren
                           tau3=(qtang3 - vsurf3_tang)/delta/ren
c     forza1=(tau1-prsint*versnorm(1))*area
c     forza2=(tau2-prsint*versnorm(2))*area
c     forza3=(tau3-prsint*versnorm(3))*area
                           forza1_v=tau1*area
                           forza2_v=tau2*area
                           forza3_v=tau3*area
                           forza1_p=(-prsurf*versnorm(1))*area
                           forza2_p=(-prsurf*versnorm(2))*area
                           forza3_p=(-prsurf*versnorm(3))*area
                           momento_x = 
     $                          (forza3_p+forza3_v)*(delta2-x20)
     $                          - (forza2_p+forza2_v)*(delta3-x30)
c     if (bar3.lt.0.055.and.bar3.gt.0.045)then 
                           delta3=delta*versnorm(3)+bar3
                           delta2=delta*versnorm(2)+bar2
                           delta1=delta*versnorm(1)+bar1
                           write(40,43)itr,bar1,bar2,bar3,
     $                          qtang1,qtang2,qtang3
c     $                          ,qintn
     $                          ,vsurf1_tang,vsurf2_tang,vsurf3_tang
                           write(41,43)itr,bar1,bar2,bar3,
     $                          q1int,q2int,q3int
c     $                          ,prsint
     $                          ,vsurf1,vsurf2,vsurf3
                           write(45,43)itr,bar1,bar2,bar3,
     $                          delta1,delta2,delta3,delta
                           write(46,43)itr,bar1,bar2,bar3,
     $                          forza1_v,forza2_v,forza3_v,delta
                           write(47,43)itr,bar1,bar2,bar3,versnorm(1),
     $                          versnorm(2),versnorm(3)
                           write(48,43)itr,bar1,bar2,bar3,
     $             qtang2*sqrt(qtang3**2.+qtang2**2.)/abs(qtang2)/delta, 
     $             qtang3*sqrt(qtang3**2.+qtang2**2.)/abs(qtang3)/delta

c               endif
c               if (bar1.gt.-0.015.and.bar1.lt.0.015.and.bar2.gt.0.)then 
c
c                     write(50,43)itr,bar1,bar2,bar3,
c     $                                qtang1,qtang2,qtang3,qintn
c                     write(51,43)itr,bar1,bar2,bar3,
c     $                       q1int,q2int,q3int,prsint
c                     write(56,43)itr,bar1,bar2,bar3,
c     $                    delta1,delta2,delta3,delta
c                     write(57,43)itr,bar1,bar2,bar3,versnorm(1),
c     $                    versnorm(2),versnorm(3)
c                     write(58,43)itr,bar1,bar2,bar3,
c     $            qtang3*sqrt(qtang2**2.+qtang3**2.)/abs(qtang3)/delta, 
c     $            qtang2*sqrt(qtang2**2.+qtang3**2.)/abs(qtang2)/delta
c
c               endif
 43                        format(i4,10(2x,e12.5))
                           
                           taumod= sqrt(tau1**2.+tau2**2.+tau3**2.)
                           
                        endif
                     end do
                  end if
               end do
            end if       
         end do
         do m=1,3
            ntau=int(itr*3-3+m)
            tauw(ntau)=taumod
            tauw1(ntau)=tau1
            tauw2(ntau)=tau2
            tauw3(ntau)=tau3
            prsw(ntau)=prsurf
         enddo
c     pause
c     stop
c     do k=1,n3
c     if(bar3.gt.x3c(k).and.bar3.le.x3c(k+1)) then
c     do j=1,n2
c     if(bar2.gt.x2c(j).and.bar2.le.x2c(j+1)) then
c     do i=1,n1
c     if(bar1.gt.x1c(i).and.bar1.le.x1c(i+1)) then
c                    forza3=pr(i,j,k)*area*versnorm(3)
c     forza2=pr(i,j,k)*area*versnorm(2)
c     forza1=pr(i,j,k)*area*versnorm(1)
c     endif
c     end do
c     end if
c     end do
c     end if       
c     end do
c     if((abs(forza3).gt.1).or.(abs(forza2).gt.1)
c     &       .or.(abs(forza1).gt.1)) then
c     write(*,*) 'SBALLATO'
c     write(*,*) itr,forza3,forza2,forza1
c     write(*,*) itr,risul3,risul2,risul1
c     goto 987
c     endif
c     write(*,*) itr,forza3,forza2,forza1
         risul3_v = risul3_v + forza3_v
         risul2_v = risul2_v + forza2_v
         risul1_v = risul1_v + forza1_v
         risul3_p = risul3_p + forza3_p
         risul2_p = risul2_p + forza2_p
         risul1_p = risul1_p + forza1_p
         risulM_x = risulM_x + momento_x
c     write(*,*) itr,risul3,risul2,risul1
 987     continue
      end do
      
      write(*,*)'forza risultante'
      write(*,*)'3 comp visc, press',risul3_v,risul3_p
      write(*,*)'2 comp visc, press',risul2_v,risul2_p    
      write(*,*)'1 comp visc, press',risul1_v,risul1_p
      write(*,*)'3 risultante_x',(risul3_v+risul3_p)
      write(*,*)'2 risultante_y',(risul2_v+risul2_p)
      write(*,*)'1 risultante_z',(risul1_v+risul1_p)
      risult3_hat = (risul2_v+risul2_p)*sin(alphax) +
     $     (risul3_v+risul3_p)*cos(alphax)
      risult2_hat = (risul2_v+risul2_p)*cos(alphax) - 
     $     (risul3_v+risul3_p)*sin(alphax)
      write(*,*)'3 risultante_xhat',risult3_hat
      write(*,*)'2 risultante_yhat',risult2_hat
      write(*,*)' sup_tot, ne ',sup_tot,nb
      sup_rif=3.00              ! RONDINE
c     cx=2.0*risul3/sup_rif
c     cy=2.0*risul2/sup_rif
c     cz=2.0*risul1/sup_rif
      cx_v = 8.0*risul3_v/sup_tot
      cy_v = 8.0*risul2_v/sup_tot
      cz_v = 8.0*risul1_v/sup_tot
      cx_p = 8.0*risul3_p/sup_tot
      cy_p = 8.0*risul2_p/sup_tot
      cz_p = 8.0*risul1_p/sup_tot
      
      write(*,*)'coefficienti'
      write(*,*)'cx v - p ',cx_v,cx_p
      write(*,*)'cy v - p ',cy_v,cy_p 
      write(*,*)'cz v - p ',cz_v,cz_p
      
      itime=nint(time*1000.)
      write(ipfi,83)itime
 83   format(i5.5)         
      namfile='geo'//ipfi//'.dat'
      open(12,file=namfile)
      write(12,*)'Variables = X, Y, Z'
      write(12,*)'Zone N = ',3*nb,', E= ',nb,
     $     ', F=FEPOINT, ET=TRIANGLE'
      do k = 1,nb
         if (refsys.eq.1)then
            v13 = x30 +( xyzb(k,2) - x20 )*sin(alphax) + 
     $           ( xyzb(k,1) - x30 )*cos(alphax)
            v12 = x20 +( xyzb(k,2) - x20 )*cos(alphax) - 
     $           ( xyzb(k,1) - x30 )*sin(alphax) + hy
            v11 = xyzb(k,3)
            v23 = x30 +( xyzb(k,5) - x20 )*sin(alphax) + 
     $           ( xyzb(k,4) - x30 )*cos(alphax)
            v22 = x20 +( xyzb(k,5) - x20 )*cos(alphax) - 
     $           ( xyzb(k,4) - x30 )*sin(alphax) + hy
            v21 = xyzb(k,6)
            v33 = x30 +( xyzb(k,8) - x20 )*sin(alphax) + 
     $           ( xyzb(k,7) - x30 )*cos(alphax)
            v32 = x20 +( xyzb(k,8) - x20 )*cos(alphax) - 
     $           ( xyzb(k,7) - x30 )*sin(alphax) + hy
            v31 = xyzb(k,9)
            write(12,*)v13,v12,v11
            write(12,*)v23,v22,v21
            write(12,*)v33,v32,v31
         else
            write(12,*)xyzb(k,1),xyzb(k,2),xyzb(k,3)
            write(12,*)xyzb(k,4),xyzb(k,5),xyzb(k,6)
            write(12,*)xyzb(k,7),xyzb(k,8),xyzb(k,9)
         endif
      enddo
      j = 1
      do k = 1,nb
         write(12,*)j,j+1,j+2
         j=j+3
      enddo
      close(12)
      
      namfile='norm.tec'
      open(31,file=namfile,form='formatted',status='unknown')
      write(31,*) 'Variables = X, Y, Z, Norm1, Norm2, Norm3'
      write(31,*)'ZONE N = ',3*nb,',E = ',nb,',F=FEPOINT,ET=TRIANGLE'
      write(31,*)'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)'
      do k = 1,nb                 
         write(31,*)xyzb(k,1),xyzb(k,2),xyzb(k,3),
     $        norb(k,1),norb(k,2),norb(k,3)
         write(31,*)xyzb(k,4),xyzb(k,5),xyzb(k,6),
     $        norb(k,1),norb(k,2),norb(k,3)
         write(31,*)xyzb(k,7),xyzb(k,8),xyzb(k,9),
     $        norb(k,1),norb(k,2),norb(k,3)
      enddo
      j = 1
      do k = 1,nb
         write(31,*)j,j+1,j+2
         j=j+3
      enddo
      close(31)
      namfile='tau.tec'
c     c     endif
      open(31,file=namfile,form='formatted',status='unknown')
      write(31,*) 'Variables = X, Y, Z, Tau(module), TauX, TauY,
     &     TauZ,Pressure'
      write(31,*)'ZONE N = ',3*nb,',E = ',nb,',F=FEPOINT,ET=TRIANGLE'
c     write(31,*)'ZONE N = ',3*nb,',E = ',nb,',F=FEPOINT, ET=TRIANGLE'
      write(31,*)'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE
     &     SINGLE SINGLE)'
      do k = 1,nb                 
         write(31,*)xyzb(k,1),xyzb(k,2),xyzb(k,3),
     $        tauw(k*3-3+1),tauw3(k*3-3+1),tauw2(k*3-3+1),
     $        tauw1(k*3-3+1),prsw(k*3-3+1)
         write(31,*)xyzb(k,4),xyzb(k,5),xyzb(k,6),
     $        tauw(k*3-3+2),tauw3(k*3-3+2),tauw2(k*3-3+2),
     $        tauw1(k*3-3+2),prsw(k*3-3+2)
         write(31,*)xyzb(k,7),xyzb(k,8),xyzb(k,9),
     $        tauw(k*3-3+3),tauw3(k*3-3+3),tauw2(k*3-3+3),
     $        tauw1(k*3-3+3),prsw(k*3-3+3)
      enddo
      j = 1
      do k = 1,nb
         write(31,*)j,j+1,j+2
         j=j+3
      enddo
      close(31)
         
cc     if(refsys.eq.1) then
cc     do l=1,nv
cc     write(31,*) vert(l,3),vert(l,2)+hy,vert(l,1)+hx,var(l)
cc     end do
cc     else
c         do l=1,nv
c            write(31,*) vert(l,3),vert(l,2),vert(l,1),tauw(l),tauw3(l),
c     &           tauw2(l),tauw1(l),prsw(l)
c         end do
cc     endif
c         do l=1,ne
c            write(31,*) ielem(l,3),ielem(l,2),ielem(l,1)
c         end do
c         close(31)

c 1222 continue
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       do n=1,npunv_i
c         kz=indgeoi(5,n,3)
c         jy=indgeoi(5,n,2)       
c         ix=indgeoi(5,n,1)
c         iv=ix + jy*n1 + kz*n1*n2
c         blankv(iv)=1.
c       enddo
c
      do  kc=1,n3m 
        km=kmv(kc)                                                    
        do  jc=1,n2m                                                     
          jm=jmv(jc)                                                        
          do  ic=1,n1m                                                     
            im=imv(ic)
            if(abs(pr(ic,jc,kc)).gt.3) pr(ic,jc,kc)=0.
          enddo
        enddo
      enddo  
            
      do 4 kc=1,n3m 
        km=kmv(kc)                                                    
        do 4 jc=1,n2m                                                     
          jm=jmv(jc)                                                        
          do 4 ic=1,n1m                                                     
            im=imv(ic)

       iq3=4
       if(q3(ic,jc,kc).eq.0.) iq3=iq3-1
       if(q3(im,jc,kc).eq.0.) iq3=iq3-1
       if(q3(ic,jm,kc).eq.0.) iq3=iq3-1
       if(q3(im,jm,kc).eq.0.) iq3=iq3-1
       if(iq3.eq.0) then
          voc3(ic,jc,kc)=0.
       else  
          voc3(ic,jc,kc)=(q3(ic,jc,kc)+q3(im,jc,kc)+
     %                    q3(ic,jm,kc)+q3(im,jm,kc))/iq3
       endif 
c      if(abs(voc3(ic,jc,kc)).gt.2.)
c    %    write(*,*)'v3',ic,jc,kc,voc3(ic,jc,kc)    
c
c
       iq2=4
       if(q2(ic,jc,kc).eq.0.) iq2=iq2-1
       if(q2(im,jc,kc).eq.0.) iq2=iq2-1
       if(q2(ic,jc,km).eq.0.) iq2=iq2-1
       if(q2(im,jc,km).eq.0.) iq2=iq2-1   
       if(iq2.eq.0) then 
          voc2(ic,jc,kc)=0. 
       else                                          
          voc2(ic,jc,kc)=(q2(ic,jc,kc)+q2(im,jc,kc)+
     %                    q2(ic,jc,km)+q2(im,jc,km))/iq2
       endif
c      if(abs(voc2(ic,jc,kc)).gt.1.)
c    %    write(*,*)'v2',ic,jc,kc,voc2(ic,jc,kc) 
c
c
       iq1=4
       if(q1(ic,jc,kc).eq.0.) iq1=iq1-1
       if(q1(ic,jm,kc).eq.0.) iq1=iq1-1
       if(q1(ic,jc,km).eq.0.) iq1=iq1-1
       if(q1(ic,jm,km).eq.0.) iq1=iq1-1   
       if(iq1.eq.0) then 
          voc1(ic,jc,kc)=0. 
       else                                          
          voc1(ic,jc,kc)=(q1(ic,jc,kc)+q1(ic,jm,kc)+
     %                    q1(ic,jc,km)+q1(ic,jm,km))/iq1
       endif
c      if(abs(voc1(ic,jc,kc)).gt.1.)
c    %    write(*,*)'v1',ic,jc,kc,voc1(ic,jc,kc)
c
c 
       ipr=8
       if(pr(ic,jc,kc).eq.0.) ipr=ipr-1
       if(pr(im,jc,kc).eq.0.) ipr=ipr-1
       if(pr(ic,jm,kc).eq.0.) ipr=ipr-1
       if(pr(im,jm,kc).eq.0.) ipr=ipr-1   
       if(pr(ic,jc,km).eq.0.) ipr=ipr-1
       if(pr(im,jc,km).eq.0.) ipr=ipr-1
       if(pr(ic,jm,km).eq.0.) ipr=ipr-1
       if(pr(im,jm,km).eq.0.) ipr=ipr-1 
       if(ipr.eq.0) then 
          prc(ic,jc,kc)=0. 
       else                                          
          prc(ic,jc,kc)=(pr(ic,jc,kc)+pr(im,jc,kc)
     %                  +pr(ic,jm,kc)+pr(im,jm,kc)
     %                  +pr(ic,jc,km)+pr(im,jc,km)
     %                  +pr(ic,jm,km)+pr(im,jm,km))/ipr
       endif
c      if (abs(prc(ic,jc,kc)).gt.1.) write(*,*)'p',ic,jc,kc,prc(ic,jc,kc)  
c
c
            densc(ic,jc,kc)=(dens(ic,jc,kc)+dens(im,jc,kc)
     %                      +dens(ic,jm,kc)+dens(im,jm,kc)
     %                      +dens(ic,jc,km)+dens(im,jc,km)
     %                      +dens(ic,jm,km)+dens(im,jm,km))*0.125

    4 continue 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
      goto 9876
      ind_up=0
      ind_lw=0
      do n=1,npunv_e   
         ke=indgeo(5,n,3)
         je=indgeo(5,n,2)
         ie=indgeo(5,n,1) 
         if (ie.eq.5) then
           if (je.gt.34) ind_up=ind_up+1
         end if 
         if (ie.eq.5) then
           if (je.lt.34) ind_lw=ind_lw+1
         end if 
      end do
c
      open(32,file='cp_up.tec',form='formatted',status='unknown')
      write(32,*) 'Variables = "X" , "Y" '
      write(32,*) 'ZONE I = ',ind_up,',   F=POINT'
      open(33,file='cp_dw.tec',form='formatted',status='unknown')
      write(33,*) 'Variables = "X" , "Y" '
      write(33,*) 'ZONE I = ',ind_lw,',   F=POINT'
       do n=1,npunv_e   
         ke=indgeo(5,n,3)
         je=indgeo(5,n,2)
         ie=indgeo(5,n,1) 
         if (ie.eq.5) then
           cp=2.0*prc(ie,je,ke)  
           if (je.gt.34) write(32,*) x3c(ke),cp
         end if 
         if (ie.eq.5) then
           cp=2.0*prc(ie,je,ke) 
           if (je.lt.34) write(33,*) x3c(ke),cp  
         end if 
      end do
      close(32)
      close(33)

      open(34,file='cp_wall.tec',form='formatted',status='unknown')
      write(34,*) 'Variables = "X" , "Y" '
      write(34,*) 'ZONE I = ',n3m,',   F=POINT'
        do n=1,n3m
           cp=2.0*prc(5,1,n)
           write(34,*) x3c(n),cp
        end do
      close(34)
 9876 continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do 6 kc=1,n3m
        km=kmv(kc)                                                    
        do 61 ic=1,n1m                                                     
          if(inslww.eq.1) then
          voc1(ic,n2,kc)=0.d0
          voc2(ic,n2,kc)=0.d0
          voc3(ic,n2,kc)=0.d0
          else
          voc1(ic,n2,kc)=voc1(ic,n2m,kc)
          voc2(ic,n2,kc)=voc2(ic,n2m,kc)
          voc3(ic,n2,kc)=voc3(ic,n2m,kc)
          endif
          prc(ic,n2,kc)=prc(ic,n2m,kc)
          densc(ic,n2,kc)= densc(ic,n2m,kc)
          if(inslwe.eq.1) then
          voc1(ic,1,kc)=0.d0
          voc2(ic,1,kc)=0.d0
          voc3(ic,1,kc)=0.d0
          else
          voc1(ic,1,kc)=voc1(ic,2,kc)
          voc2(ic,1,kc)=voc2(ic,2,kc)
          voc3(ic,1,kc)=voc3(ic,2,kc)
          endif
          prc(ic,1,kc)=prc(ic,2,kc)
          densc(ic,1,kc)=densc(ic,2,kc)
   61   continue
    6 continue  
c
      do 15 jc=1,n2                                                    
        do 15 ic=1,n1
          voc3(ic,jc,n3)=q3(ic,jc,n3)
          voc3(ic,jc,1)=q3(ic,jc,1)
          prc(ic,jc,1)=prc(ic,jc,2)
          prc(ic,jc,n3)=prc(ic,jc,n3m)
          densc(ic,jc,1)=densc(ic,jc,2)
          densc(ic,jc,n3)=densc(ic,jc,n3m)
c         if(inslws.eq.1) then
c           voc2(ic,jc,1)=0.d0            
c           voc1(ic,jc,1)=0.d0            
c         else
            voc2(ic,jc,1)=voc2(ic,jc,2)
            voc1(ic,jc,1)=voc1(ic,jc,2)
c         endif
c         if(inslwn.eq.1) then
c           voc2(ic,jc,n3)=0.d0            
c           voc1(ic,jc,n3)=0.d0            
c         else
            voc2(ic,jc,n3)=voc2(ic,jc,n3m)
            voc1(ic,jc,n3)=voc1(ic,jc,n3m)
c         endif
   15 continue  
      do kc=1,n3                                                    
        do jc=1,n2
          voc1(n1,jc,kc)=voc1(1,jc,kc)
          voc2(n1,jc,kc)=voc2(1,jc,kc)
          voc3(n1,jc,kc)=voc3(1,jc,kc)
          prc(n1,jc,kc)=prc(1,jc,kc)
          densc(n1,jc,kc)=densc(1,jc,kc)
        end do
      end do
      return
      end
c                                                                       
c                                                                       
c  ****************************** subrout vort  **********************  
c                                                                       
c     this subroutine calculates the azimuthal vorticity               
c                                                                       
      subroutine vort(vora,q1,q2,q3,voc1,voc2,voc3,time)                    
      include 'parade_new.f'

      dimension q1(m1,m2,m3), q2(m1,m2,m3), q3(m1,m2,m3)
      dimension vora(m1,m2,m3) 
      dimension voc1(m1,m2,m3),voc2(m1,m2,m3),voc3(m1,m2,m3)  
      dimension voc1m(m2,m3)
      dimension vmax(3),vocmax(4)  

c     go to 1234
c
c  

c
c                                                                       
c  ***********  compute the azimuthal vorticity component               
c                                                                       
      vo1max=0.
      vc1max=-1000.
      vc1min=1000.
      do 1 jc=2,n2m                                               
         jm=jmv(jc)                                                  
         udx2=dx2/g2c(jc)
         do 1 kc=2,n3m                                               
            udx3=dx3/g3c(kc)
            km=kmv(kc)                                                  
            do 1 ic=1,n1m                                               
               im=imv(ic)
               dq2x3=(q2(ic,jc,kc)-q2(ic,jc,km))*udx3
               dq2x3m=(q2(im,jc,kc)-q2(im,jc,km))*udx3
               dq3x2=(q3(ic,jc,kc)-q3(ic,jm,kc))*udx2                        
               dq3x2m=(q3(im,jc,kc)-q3(im,jm,kc))*udx2                    
               voraz=dq2x3-dq3x2                                            
               vorazm=dq2x3m-dq3x2m                                        
               voc1(ic,jc,kc)=(voraz+vorazm)*0.5 
c               vo1max=max(vo1max,abs(voc1(ic,jc,kc))) 
c               vc1max=max(vc1max,voc1(ic,jc,kc)) 
c               vc1min=min(vc1min,voc1(ic,jc,kc)) 
               if(vo1max.lt.abs(voc1(ic,jc,kc)))then
                  vo1max = abs(voc1(ic,jc,kc))
                  vo1max_x1 = x1c(ic)
                  vo1max_x2 = x2c(jc)
                  vo1max_x3 = x3c(kc)
               endif
               if(vc1max.lt.voc1(ic,jc,kc))then
                  vc1max = voc1(ic,jc,kc)
                  vc1max_x1 = x1c(ic)
                  vc1max_x2 = x2c(jc)
                  vc1max_x3 = x3c(kc)
               endif
               if(vc1min.gt.voc1(ic,jc,kc))then
                  vc1min = voc1(ic,jc,kc)
                  vc1min_x1 = x1c(ic)
                  vc1min_x2 = x2c(jc)
                  vc1min_x3 = x3c(kc)
               endif
 1    continue      
      write(*,*)' vo1max ',vo1max,' vo1max_x1 ',vo1max_x1,
     $     ' vo1max_x2 ',vo1max_x2,' vo1max_x3 ',vo1max_x3
      write(*,*)' vc1max ',vc1max,' vc1max_x1 ',vc1max_x1,
     $     ' vc1max_x2 ',vc1max_x2,' vc1max_x3 ',vc1max_x3
      write(*,*)' vc1min ',vc1min,' vc1min_x1 ',vc1min_x1,
     $     ' vc1min_x2 ',vc1min_x2,' vc1min_x3 ',vc1min_x3
c                                                                  
c  ***********  compute the radial  vorticity component            
c                                                                  
      vo2max=0.
      vc2max=-1000.
      vc2min=1000.
      do 2 ic=1,n1m                                                
         im=imv(ic)                                                  
         do 2 jc=1,n2m                                               
            jm=jmv(jc)
            do 2 kc=2,n3m                                               
               km=kmv(kc)                                                   
               dq3x1=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1
               dq3x1m=(q3(ic,jm,kc)-q3(im,jm,kc))*dx1
               dq1x3=(q1(ic,jc,kc)-q1(ic,jc,km))*dx3/g3c(kc)
               dq1x3m=(q1(ic,jm,kc)-q1(ic,jm,km))*dx3/g3c(kc)
               vorr=dq3x1-dq1x3                                            
               vorrm=dq3x1m-dq1x3m                                          
               voc2(ic,jc,kc)=(vorr+vorrm)*0.5 
c               vo2max=max(vo2max,abs(voc2(ic,jc,kc))) 
c               vc2max=max(vc2max,voc2(ic,jc,kc)) 
c               vc2min=min(vc2min,voc2(ic,jc,kc)) 
               if(vo2max.lt.abs(voc2(ic,jc,kc)))then
                  vo2max = abs(voc2(ic,jc,kc))
                  vo2max_x1 = x1c(ic)
                  vo2max_x2 = x2c(jc)
                  vo2max_x3 = x3c(kc)
               endif
               if(vc2max.lt.voc2(ic,jc,kc))then
                  vc2max = voc2(ic,jc,kc)
                  vc2max_x1 = x1c(ic)
                  vc2max_x2 = x2c(jc)
                  vc2max_x3 = x3c(kc)
               endif
               if(vc2min.gt.voc2(ic,jc,kc))then
                  vc2min = voc2(ic,jc,kc)
                  vc2min_x1 = x1c(ic)
                  vc2min_x2 = x2c(jc)
                  vc2min_x3 = x3c(kc)
               endif
 2    continue     
      write(*,*)' vo2max ',vo2max,' vo2max_x1 ',vo2max_x1,
     $     ' vo2max_x2 ',vo2max_x2,' vo2max_x3 ',vo2max_x3
      write(*,*)' vc2max ',vc2max,' vc2max_x1 ',vc2max_x1,
     $     ' vc2max_x2 ',vc2max_x2,' vc2max_x3 ',vc2max_x3
      write(*,*)' vc2min ',vc2min,' vc2min_x1 ',vc2min_x1,
     $     ' vc2min_x2 ',vc2min_x2,' vc2min_x3 ',vc2min_x3
c                                                                  
c  ***********  compute the vertical vorticity component        
c                                
      vo3max=0.                                       
      vc3max=-1000.
      vc3min=1000.
      do 3 ic=1,n1m                                               
         im=imv(ic)                                                  
         do 3 jc=1,n2m                                               
            jm=jmv(jc)                                                   
            udx1=dx1
            udx2=dx2/g2c(jc)
            do 3 kc=2,n3m                                                
               km=kmv(kc)
               dq1x2=(q1(ic,jc,kc)-q1(ic,jm,kc))              
     $              *udx2                                                   
               dq1x2m=(q1(ic,jc,km)-q1(ic,jm,km))          
     $              *udx2                                                   
               dq2x1=(q2(ic,jc,kc)-q2(im,jc,kc))*udx1
               dq2x1m=(q2(ic,jc,km)-q2(im,jc,km))*udx1
               vorz=(dq1x2-dq2x1)
               vorzm=(dq1x2m-dq2x1m)
               voc3(ic,jc,kc)=(vorz+vorzm)*0.5 
c               vo3max=max(vo3max,abs(voc3(ic,jc,kc))) 
c               vc3max=max(vc3max,voc3(ic,jc,kc)) 
c               vc3min=min(vc3min,voc3(ic,jc,kc)) 
               if(vo3max.lt.abs(voc3(ic,jc,kc)))then
                  vo3max = abs(voc3(ic,jc,kc))
                  vo3max_x1 = x1c(ic)
                  vo3max_x2 = x2c(jc)
                  vo3max_x3 = x3c(kc)
               endif
               if(vc3max.lt.voc3(ic,jc,kc))then
                  vc3max = voc3(ic,jc,kc)
                  vc3max_x1 = x1c(ic)
                  vc3max_x2 = x2c(jc)
                  vc3max_x3 = x3c(kc)
               endif
               if(vc3min.gt.voc3(ic,jc,kc))then
                  vc3min = voc3(ic,jc,kc)
                  vc3min_x1 = x1c(ic)
                  vc3min_x2 = x2c(jc)
                  vc3min_x3 = x3c(kc)
               endif
 3    continue     
      write(*,*)' vo3max ',vo3max,' vo3max_x1 ',vo3max_x1,
     $     ' vo3max_x2 ',vo3max_x2,' vo3max_x3 ',vo3max_x3
      write(*,*)' vc3max ',vc3max,' vc3max_x1 ',vc3max_x1,
     $     ' vc3max_x2 ',vc3max_x2,' vc3max_x3 ',vc3max_x3
      write(*,*)' vc3min ',vc3min,' vc3min_x1 ',vc3min_x1,
     $     ' vc3min_x2 ',vc3min_x2,' vc3min_x3 ',vc3min_x3
c
c     ****************** vorticities at the upper and lower walls
c 
      do 234 jc=2,n2m                                               
         jm=jmv(jc)                                                  
         do 234 kc=1,n3,n3m                                            
            km=kmv(kc)                                                  
            do 234 ic=1,n1m                                               
      im=imv(ic)
      if(kc.eq.1) then     
         if(inslws.eq.0) then
            voc3(ic,jc,kc)=voc3(ic,jc,kc+1)
            dq2x3=0.
            dq2x3m=0.
            dq1x3=0.                                             
            dq1x3m=0.                                             
         else 
            voc3(ic,jc,kc)=0.              
            dq1x3=(q1(ic,jc,kc)*2.          )*dx3/g3c(kc)       
            dq1x3m=(q1(ic,jm,kc)*2.          )*dx3/g3c(kc)       
            dq2x3=(q2(ic,jc,kc)*2.          )*dx3
            dq2x3m=(q2(im,jc,kc)*2.          )*dx3
         endif
         dq3x1=0.                                                    
         dq3x1m=0.                                                    
         dq3x2=0.                                                     
         dq3x2m=0.                                                 
         voraz=dq2x3-dq3x2                                            
         vorazm=dq2x3m-dq3x2m                                        
         voc1(ic,jc,kc)=(voraz+vorazm)*0.5 
         vorr=dq3x1-dq1x3                                            
         vorrm=dq3x1m-dq1x3m                                          
         voc2(ic,jc,kc)=(vorr+vorrm)*0.5 
      else
         if(inslwn.eq.0) then
            voc3(ic,jc,kc)=voc3(ic,jc,kc-1)
            dq2x3=0.
            dq2x3m=0.
            dq1x3=0.                                             
            dq1x3m=0.                                             
         else 
            voc3(ic,jc,kc)=0.              
            dq1x3=-(q1(ic,jc,kc)*2.          )*dx3/g3c(kc)       
            dq1x3m=-(q1(ic,jm,kc)*2.          )*dx3/g3c(kc)       
            dq2x3=-(q2(ic,jc,kc)*2.          )*dx3/g3c(kc)        
            dq2x3m=-(q2(im,jc,kc)*2.          )*dx3/g3c(kc)   
         endif
         dq3x1=0.                                                    
         dq3x1m=0.                                                    
         dq3x2=0.                                                     
         dq3x2m=0.                                                 
         voraz=dq2x3-dq3x2                                            
         vorazm=dq2x3m-dq3x2m                                        
         voc1(ic,jc,kc)=(voraz+vorazm)*0.5 
         vorr=dq3x1-dq1x3                                            
         vorrm=dq3x1m-dq1x3m                                          
         voc2(ic,jc,kc)=(vorr+vorrm)*0.5 
      endif 
 234  continue    
c                                                                  
      jc=n2
      do 5 ic=1,n1m                                               
         im=imv(ic)                                                  
         do 5 kc=1,n3 
            km=kmv(kc)                                                 
            voc2(ic,jc,kc)=0.
            voc1(ic,jc,kc)=0.
            voc3(ic,jc,kc)=0.
 5    continue                             
c     
      jc=1
      do 55 ic=1,n1m                                               
         im=imv(ic)                                                  
         do 55 kc=1,n3 
            km=kmv(kc)                                                 
            voc2(ic,jc,kc)=0.
            voc1(ic,jc,kc)=0.
            voc3(ic,jc,kc)=0.
 55   continue                             
      do 11 l=1,3
        vmax(l) = 0.
         do 12 ic=1,n1m                                              
            do 12 jc=1,n2                                               
               do 12 kc=1,n3
                  if(l.eq.1) vca=q1(ic,jc,kc)  
                  if(l.eq.2) vca=q2(ic,jc,kc)  
                  if(l.eq.3) vca=q3(ic,jc,kc)
                  if(abs(vca).ge.vmax(l)) vmax(l)=abs(vca)  
 12      continue                                                     
 11   continue                                                     
c     
      do 21 l=1,3
         vocmax(l)=0.
         do 7 ic=1,n1m                                                
            do 7 jc=1,n2m                                                     
               do 7 kc=1,n3m
                  if(l.eq.1) voca=voc1(ic,jc,kc)  
                  if(l.eq.2) voca=voc2(ic,jc,kc)  
                  if(l.eq.3) voca=voc3(ic,jc,kc)
                  if(abs(voca).ge.vocmax(l)) vocmax(l)=abs(voca)  
 7       continue                                                          
 21   continue                                                          
c     
      vcamax=0.                 !FRANC non def
      vcamin=0.                 !FRANC non def
      voamax=0.                 !FRANC non def
      write(6,277)vc1max,vc2max,vc3max,vcamax
      write(32,277)vc1max,vc2max,vc3max,vcamax
 277  format(2x,'max vc1,vc2,vc3,vca',2x,4e12.3)
      write(6,278)vc1min,vc2min,vc3min,vcamin
      write(32,278)vc1min,vc2min,vc3min,vcamin
 278  format(2x,'min vc1,vc2,vc3,vca',2x,4e12.3)
      write(6,177)vo1max,vo2max,vo3max,voamax
      write(32,177)vo1max,vo2max,vo3max,voamax
  177 format(2x,'vo1m,vo2m,vo3m,voa',2x,4e12.3)
      write(6,176) (vocmax(l),l=1,4) 
      write(32,176) (vocmax(l),l=1,4) 
  176 format(2x,' voc1m,voc2m,voc3m,voram',2x,4e12.3)
      write(32,178) (vmax(l),l=1,3) 
      write(6,178) (vmax(l),l=1,3) 
  178 format(2x,' v1m,v2m,v3m',2x,3e12.3)
c     do 14 jc=1,n2                                                    
c     do 14 kc=1,n3m 
c     voc3(n1,jc,kc)=voc3(1,jc,kc)
c     voc2(n1,jc,kc)=voc2(1,jc,kc)
c     voc1(n1,jc,kc)=voc1(1,jc,kc)
c  14 continue                                                          
c
c                                                                       
c1234 continue
      do kc=2,n3m                                               
        udx3=dx3/g3c(kc)
        udx3m=dx3/g3m(kc)
        do jc=2,n2m                                               
          jm=jmv(jc)                                                  
          udx2=dx2/g2c(jc)
          udx2m=dx2/g2m(jc)
          km=kmv(kc)                                                  
          do ic=1,n1m                                               
            im=imv(ic)
            dq1x3=(q1(ic,jc,kc)-q1(ic,jc,km))*udx3
            dq2x3=(q2(ic,jc,kc)-q2(ic,jc,km))*udx3
c           dq3x3=(q3(ic,jc,kc)-q3(ic,jc,km))*udx3m
            dq1x2=(q1(ic,jc,kc)-q1(ic,jm,kc))*udx2                        
c           dq2x2=(q2(ic,jc,kc)-q2(ic,jm,kc))*udx2m
            dq3x2=(q3(ic,jc,kc)-q3(ic,jm,kc))*udx2                        
            dq3x1=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1                        
            dq2x1=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1                        
c           dq1x1=(q1(ic,jc,kc)-q1(im,jc,kc))*dx1                        
            vora(ic,jc,kc)= -(dq2x1*dq1x2+dq2x3*dq3x2+dq3x1*dq1x3)
          end do
        end do
      end do
      return                                                            
      end                                                               

c                                                                       
c
c   ********************************************************************
c 
      subroutine cordin
      include 'parade_new.f'

      dimension eta(m2),etaz(m3)
      dimension xt31(m3),xt32(m3)

c
      pi=2.*asin(1.)
c
c     AZIMUTHAL COORDINATE DEFINITION
c
      if(igext.ne.1) then
      do i=1,n1
        x1c(i)= float(i-1)/dx1
      end do
      x1c(0)=-1/dx1    !laura
c      x1c(n1+1)=x1c(n1)+1/dx1   !laura!prova laura 24gennaio06
      end if
c      do i=1,n1m                    !laura
      write(*,*)(x1c(i),i=1,n1)
      do i=0,n1                 !laura                                                   
         ip=i+1
c     x1m(i)= (x1c(i)+x1c(ip))*.5  !prova laura 24gennaio06
         x1m(i)= (float(i-1)+0.5)/dx1 !prova laura 24gennaio06
      end do

c
c     RADIAL COORDINATE DEFINITION
c
      rint = 0.
      r0 = rext 
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
          x2c(j)=rcdp*rext
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
          x2c(j)=rint+eta(j)*(rext-rint)
        end do
      end if
      x2c(1)=rint
C
C    MESH FROM AN EXTERNAL INPUT FILE
C
      if (istr.eq.3) then
        open(66,file='gridy.data',form='formatted',status='unknown')
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
        alx2h = rext*0.5
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
      x2c(j)=+eta(j)*rext
      end do
      end if
      end if

CCC   Aggiunta
c
c     STAGGERED COORDINATES AND
c     METRIC QUANTITIES
c
      do j=1,n2m
        x2m(j)=(x2c(j)+x2c(j+1))*0.5
        g2m(j)=(x2c(j+1)-x2c(j))*dx2
      end do
      x2m(n2)=x2c(n2)
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
C
C    MESH FROM AN EXTERNAL INPUT FILE
C
      if (istr3.eq.3) then
        open(66,file='gridz.data',form='formatted',status='unknown')
        do k=1,n3
          read(66,*) x3,x3c(k)
        end do
        close(66)
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
      end if
c
c     STAGGERED COORDINATES AND
c     METRIC QUANTITIES
c
      do k=1,n3m
        x3m(k)=(x3c(k)+x3c(k+1))*0.5
        g3m(k)=(x3c(k+1)-x3c(k))*dx3
      enddo
      x3m(n3)=x3c(n3)
      do k=2,n3m
        g3c(k)=(x3c(k+1)-x3c(k-1))*dx3*0.5
      enddo
      g3c(1)=(x3c(2)-x3c(1))*dx3
      g3c(n3)= (x3c(n3)-x3c(n3m))*dx3
c
c     WRITE GRID INFORMATION
c
c      open(unit=98,file='radcor.outp',status='unknown')
c      do j=1,n2
c        write(98,*) j,float(j)/float(n2m),x2c(j),x2m(j),g2c(j),g2m(j)
c      end do
c      close(98)
c      open(unit=78,file='axicor.outp',status='unknown')
c      do k=1,n3
c        write(78,*) k,float(k)/float(n3m),x3c(k),g3c(k),g3m(k)
c      end do
c      close(78)
      return                                                            
      end                                                               
c                                                                       
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
      include 'parade_new.f'

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
      write(6,*)nx,ny,nz
      if(nz.gt.m3) then
         write(6,*)nx,ny,nz
         write(6,*) ' WARNING: nz gt m3 '
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
c
      read(60,*) (x3c(k),k=1,nz)
      read(60,*) (x2c(j),j=1,ny)
      read(60,*) (x1c(i),i=1,nx)
c     do k=1,nz
c        read(60,*) x3c(k)
c     enddo
c     do j=1,ny
c        read(60,*) x2c(j)
c     enddo
c     do i=1,nx
c        read(60,*) x1c(i)
c     enddo
c     q3
c
      read(60,*) npunz_i
      write(6,*) ' Number of internal points for q3 = ',npunz_i
      if(npunz_i.gt.1100000) 
     %  write(*,*) ' WARNING: dimension mpun too small for q3'
c     read(60,*) (indgeoi(3,n,3),indgeoi(3,n,2),indgeoi(3,n,1)
c    %            ,n=1,npunz_i)
      do n=1,npunz_i
        read(60,*) indgeoi(3,n,3),indgeoi(3,n,2),indgeoi(3,n,1)
      enddo            
      read(60,*) npunz_e
      write(6,*) ' Number of external points for q3 = ',npunz_e
      if(npunz_e.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for q3'
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
     %  write(*,*) ' WARNING: dimension mpun too small for q3'
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
     %  write(*,*) ' WARNING: dimension mpun too small for q2'
c     read(60,*) (indgeoi(2,n,3),indgeoi(2,n,2),indgeoi(2,n,1)
c    %            ,n=1,npunr_i)
      do n=1,npunr_i
      read(60,*) indgeoi(2,n,3),indgeoi(2,n,2),indgeoi(2,n,1)
      enddo
      read(60,*) npunr_e
      write(6,*) ' Number of external points for q2 = ',npunr_e
      if(npunr_e.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for q2'
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
     %  write(*,*) ' WARNING: dimension mpun too small for q2'
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
     %  write(*,*) ' WARNING: dimension mpun too small for q1'
c     read(60,*) (indgeoi(1,n,3),indgeoi(1,n,2),indgeoi(1,n,1)
c    %            ,n=1,npunt_i)
      do n=1,npunt_i
      read(60,*) indgeoi(1,n,3),indgeoi(1,n,2),indgeoi(1,n,1)
      enddo
      read(60,*) npunt_e
      write(6,*) ' Number of external points for q1 = ',npunt_e
      if(npunt_e.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for q1'
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
     %  write(*,*) ' WARNING: dimension mpun too small for q1'
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
c     pres
c
      read(60,*) npunp_i
      write(6,*) ' Number of internal points for pres = ',npunp_i
      if(npunp_i.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for pres'
c     read(60,*) (indgeoi(4,n,3),indgeoi(4,n,2),indgeoi(4,n,1)
c    %            ,n=1,npunp_i)
      do n=1,npunp_i
      read(60,*) indgeoi(4,n,3),indgeoi(4,n,2),indgeoi(4,n,1)
      enddo
      read(60,*) npunp_e
      write(6,*) ' Number of external points for pres = ',npunp_e
      if(npunp_e.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for pres'
c     read(60,*) (indgeo(4,n,3),indgeo(4,n,2),indgeo(4,n,1),
c    %            indgeoe(4,n,3),indgeoe(4,n,2),indgeoe(4,n,1),
c    %            distb(4,n),n=1,npunp_e)
      do n=1,npunp_e
      read(60,*) indgeo(4,n,3),indgeo(4,n,2),indgeo(4,n,1),
     %            indgeoe(4,n,3),indgeoe(4,n,2),indgeoe(4,n,1),
     %            distb(4,n)
      enddo
c
      npunp_i = 0 !PRESSURE
      npunp = npunp_e + npunp_i
      if(npunp.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for pres'
      do n=1,npunp_i
        ll = npunp_e+n 
        indgeo(4,ll,1) = indgeoi(4,n,1)
        indgeo(4,ll,2) = indgeoi(4,n,2)
        indgeo(4,ll,3) = indgeoi(4,n,3)
        indgeoe(4,ll,1) = indgeoi(4,n,1)
        indgeoe(4,ll,2) = indgeoi(4,n,2)
        indgeoe(4,ll,3) = indgeoi(4,n,3)
        distb(4,ll) = 0.
      end do
c
c     vert
c
      read(60,*) npunv_i
      write(6,*) ' Number of internal points for vert = ',npunv_i
      if(npunv_i.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for vert'
c     read(60,*) (indgeoi(5,n,3),indgeoi(5,n,2),indgeoi(5,n,1)
c    %            ,n=1,npunv_i)
      do n=1,npunv_i
      read(60,*) indgeoi(5,n,3),indgeoi(5,n,2),indgeoi(5,n,1)
      enddo
      read(60,*) npunv_e
      write(6,*) ' Number of external points for vert = ',npunv_e
      if(npunv_e.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for vert'
c     read(60,*) (indgeo(5,n,3),indgeo(5,n,2),indgeo(5,n,1),
c    %            indgeoe(5,n,3),indgeoe(5,n,2),indgeoe(5,n,1),
c    %            distb(5,n),n=1,npunv_e)
      do n=1,npunv_e
      read(60,*) indgeo(5,n,3),indgeo(5,n,2),indgeo(5,n,1),
     %            indgeoe(5,n,3),indgeoe(5,n,2),indgeoe(5,n,1),
     %            distb(5,n)
      enddo
c
      npunv_i  = 0 !PRESSURE
      npunv = npunv_e + npunv_i
      if(npunv.gt.mpun) 
     %  write(*,*) ' WARNING: dimension mpun too small for vert'
      do n=1,npunv_i
        ll = npunv_e+n 
        indgeo(5,ll,1) = indgeoi(5,n,1)
        indgeo(5,ll,2) = indgeoi(5,n,2)
        indgeo(5,ll,3) = indgeoi(5,n,3)
        indgeoe(5,ll,1) = indgeoi(5,n,1)
        indgeoe(5,ll,2) = indgeoi(5,n,2)
        indgeoe(5,ll,3) = indgeoi(5,n,3)
        distb(5,ll) = 0.
      end do
c
      close(60)
c
c     DEBUGGING
c
c      n1mh = n1m/2+1
c      open(71,file='debugz.out',form='formatted',status='unknown')
c      do n=1,npunz
c        if(indgeo(2,n,1).eq.n1mh) write(71,71)
c     %  indgeo(3,n,3),indgeo(3,n,2),distb(3,n)
c      end do
c      close(71)
c 71   format(2(2x,i3),2x,f14.7)
      return
      end

C************************************************************************************

      subroutine rotation (tin)
      include 'parade_new.f'

c    
c    uhx e' la velocita relativa alla oscillazione di
c    traslazione nella direzione x spanwise nel sistema di
c    riferimento assoluto
c
      write(*,*)' ROTATION '
      pi=2.*asin(1.)
      sigma=1.
      amp=amplit
      alphaxmax=thetamax*pi/180 
      alphaymax=0.

      !velocita' di traslazione
      hy=1./amp*sin(amp*tin)
      uhy=cos(amp*tin)
      duhy=-amp*sin(amp*tin)
      write(*,*)' hy = ',hy
      write(*,*)' uhy = ',uhy
      write(*,*)' duhy = ',duhy
      hx=0.
      uhx=0.
      duhx=0.
c    velocita' di inflow 
c   ATTENZIONE E' NEGATIVA POICHE' E' LA VELOCITA' DEL CORPO
      uinflow=0.
      
      alphax = alphaxmax*cos(amp*tin)
      wx = - alphaxmax*amp*sin(amp*tin)
      dwx = - alphaxmax*amp**2*cos(amp*tin)
      write(*,*)' alphax(deg)= ',alphax/pi*180.
      write(*,*)' wx= ',wx
      write(*,*)' dwx= ',dwx

      alphay=-alphaymax*cos(sigma*tin)
      wy=alphaymax*sigma*sin(sigma*tin)
      dwy=alphaymax*sigma**2.*cos(sigma*tin)
      write(*,*)' alphay(deg)= ',alphay/pi*180.
      write(*,*)' wy= ',wy
      write(*,*)' dwy= ',dwy

      wz=0.
      dwz=0.
      u0x=uhx*cos(alphay)-uinflow*sin(alphay)
c     write(*,*)' u0x= ',u0x
      du0x=duhx*cos(alphay)-wy*uhx*sin(alphay)-wy*uinflow*cos(alphay)
c     write(*,*)' du0x= ',du0   
      u0y=+uhy*cos(alphax)
      du0y=+duhy*cos(alphax)-wx*uhy*sin(alphax)
      u0z=uhx*sin(alphay)+uinflow*cos(alphay)-uhy*sin(alphax)
      du0z=duhx*sin(alphay)+wy*uhx*cos(alphay)-uinflow*wy*sin(alphay)
     %     -duhy*sin(alphax)-uhy*wx*cos(alphax)  

c     Centro di rotazione
      x10=0
      x20=0 
      x30=-0.16667
      write(*,*)' Posizione centro di rotazione: ',x10,x20,x30
c
       write(6,*)'wx= ',wx,'u0x= ',u0x,'du0x= ',du0x
       write(6,*)'u0z= ',u0z,'du0z= ',du0z
      return
      end

c                                                                       
c                                                                       
c  ****************************** subrout eigen  **********************  
c                                                                       
c     this subroutine calculates lambda_2, lambda_ci, lambda_cr
c                                                                       
      subroutine eigen(neigen,q1,q2,q3,eig1,eig2,time)
      include 'parade_new.f'

      dimension q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)  
      dimension eig1(m1,m2,m3),eig2(m1,m2,m3)
      dimension tenso(3,3),defvel(3,3), rotvel(3,3)

      real*4 pten,qten,rten

      do kc=2,n3m    
         km=kmv(kc)                                    
         udx3=dx3/g3c(kc)
         udx3m=dx3/g3m(kc)
         do jc=2,n2m                                               
            jm=jmv(jc)                                                  
            udx2=dx2/g2c(jc)
            udx2m=dx2/g2m(jc)            
            do ic=1,n1m                                               
               im=imv(ic)
c     calcola derivate della velocita'
               dq1x3=(q1(ic,jc,kc)-q1(ic,jc,km))*udx3
               dq2x3=(q2(ic,jc,kc)-q2(ic,jc,km))*udx3
               dq3x3=(q3(ic,jc,kc)-q3(ic,jc,km))*udx3m
               dq1x2=(q1(ic,jc,kc)-q1(ic,jm,kc))*udx2                        
               dq2x2=(q2(ic,jc,kc)-q2(ic,jm,kc))*udx2m
               dq3x2=(q3(ic,jc,kc)-q3(ic,jm,kc))*udx2                        
               dq3x1=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1                        
               dq2x1=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1                        
               dq1x1=(q1(ic,jc,kc)-q1(im,jc,kc))*dx1
c     riempie matrice D*D+OMEGA*OMEGA
               if(neigen.eq.1)then 
                  defvel(1,1)= dq1x1
                  defvel(1,2)= (dq1x2+dq2x1)*0.5
                  defvel(2,1)= (dq1x2+dq2x1)*0.5
                  defvel(1,3)= (dq1x3+dq3x1)*0.5
                  defvel(3,1)= (dq1x3+dq3x1)*0.5
                  defvel(2,2)= dq2x2
                  defvel(2,3)= (dq2x3+dq3x2)*0.5
                  defvel(3,2)= (dq2x3+dq3x2)*0.5
                  defvel(3,3)= dq3x3
                  rotvel(1,1)= 0.
                  rotvel(1,2)= (dq1x2-dq2x1)*0.5
                  rotvel(2,1)= -(dq1x2-dq2x1)*0.5
                  rotvel(1,3)= (dq1x3-dq3x1)*0.5
                  rotvel(3,1)= -(dq1x3-dq3x1)*0.5
                  rotvel(2,2)= 0.
                  rotvel(2,3)= (dq2x3-dq3x2)*0.5
                  rotvel(3,2)= -(dq2x3-dq3x2)*0.5
                  defvel(3,3)= 0.
                  do nn=1,3
                     do mm=1,3
                        tenso(nn,mm)=0.
                        do ll=1,3
                           tenso(nn,mm)=tenso(nn,mm)+
     $           defvel(nn,ll)*defvel(ll,mm)+rotvel(nn,ll)*rotvel(ll,mm)
                        enddo
                     enddo
                  enddo
c     riempie matrice grad ( V )                 
               elseif(neigen.eq.2)then
                  tenso(1,1)=dq1x1
                  tenso(2,1)=dq1x2
                  tenso(3,1)=dq1x3 
                  tenso(1,2)=dq2x1
                  tenso(2,2)=dq2x2
                  tenso(3,2)=dq2x3 
                  tenso(1,3)=dq3x1
                  tenso(2,3)=dq3x2
                  tenso(3,3)=dq3x3
               endif
c               write(6,*)'tenso',((tenso(i,j),i=1,3),j=1,3)
               call calcPQR(tenso,pten,qten,rten)
               call calcEIG(eig1(ic,jc,kc),eig2(ic,jc,kc),
     $              pten,qten,rten,neigen)
c               write(6,*)'ic,jc,kc,eig1(ic,jc,kc)',
c     $              ic,jc,kc,eig1(ic,jc,kc)
            end do
         end do
      end do
      
      return
      end

      subroutine calcPQR(tenso,pten,qten,rten)
      implicit real*4(a-h,o-z)
      
      dimension tenso(3,3),doten(3,3),trten(3,3) 
      real*4 pten,qten,rten      
      
      do nn=1,3
         do mm=1,3
            doten(nn,mm)=0.
            do ll=1,3
               doten(nn,mm)=doten(nn,mm)+
     $              tenso(nn,ll)*tenso(ll,mm)
            enddo
         enddo
      enddo
      do nn=1,3
         do mm=1,3
            trten(nn,mm)=0.
            do ll=1,3
               trten(nn,mm)=trten(nn,mm)+
     $              doten(nn,ll)*tenso(ll,mm)
            enddo
         enddo
      enddo
      pten=-(tenso(1,1)+tenso(2,2)+tenso(3,3))
      qten=( pten*pten-(doten(1,1)+doten(2,2)+doten(3,3)))*0.5
      rten=( -pten*pten*pten+3.*pten*qten-          
     $     ( trten(1,1)+trten(2,2)+trten(3,3) )   )/3.
      
      return
      end

      subroutine calcEIG(evalue1,evalue2,pten,qten,rten,neigen)
      implicit real*4(a-h,o-z)

      PARAMETER (TRIPOW= 0.333333333)      
      DIMENSION EIVAL(3,2),TRIP(2,2),CONS(2,2)
      real*4 pten,qten,rten 
      REAL*4 REDEL,IMDEL
      real*4 modu1, modu2

      pi=2.*asin(1.)
      radtre=sqrt(3.)
      
      tildar=rten + 2.*pten*pten*pten/27. - pten*qten/3.
      tildaq=qten - pten*pten/3.
      
      delta = tildar*tildar*0.5*0.5+tildaq*tildaq*tildaq/27.
      IF(DELTA.GE.0.0)THEN
         TRIP(1,1) = -tildar/2.+SQRT(DELTA)
         IF(TRIP(1,1).GE.0.0)THEN
            TRIP(1,1) = (TRIP(1,1))**TRIPOW   
         ELSE 
            TRIP(1,1) =-(-1.0*TRIP(1,1))**TRIPOW
         ENDIF
         TRIP(2,1) = -tildar/2.-SQRT(DELTA)
         IF(TRIP(2,1).GE.0.0)THEN
            TRIP(2,1) = (TRIP(2,1))**TRIPOW   
         ELSE
            TRIP(2,1) =-(-1.0*TRIP(2,1))**TRIPOW
         ENDIF
         TRIP(1,2) = 0.0
         TRIP(2,2) = 0.0
      ELSEIF(DELTA.LT.0.0)THEN
         REDEL = 0.0
         IMDEL = SQRT(-1.*DELTA)
         
         TRIP(1,1) = -tildar/2.+REDEL
         TRIP(1,2) =  IMDEL
         TRIP(2,1) = -tildar/2.-REDEL
         TRIP(2,2) = -IMDEL
         
         MODU1 = SQRT(TRIP(1,1)**2.+TRIP(1,2)**2.)
         MODU2 = SQRT(TRIP(2,1)**2.+TRIP(2,2)**2.)
         
         IF(TRIP(1,1).EQ.0.0)THEN                  
            IF(TRIP(1,2).LT.0.0)ATAN1=-3./2.*PI
            IF(TRIP(1,2).GT.0.0)ATAN1= 1./2.*PI
         ELSE 
            ATAN1 = ATAN(TRIP(1,2)/TRIP(1,1))
         ENDIF
         
         IF(TRIP(2,1).EQ.0.0)THEN                  
            IF(TRIP(2,2).LT.0.0)ATAN2=-3./2.*PI
            IF(TRIP(2,2).GT.0.0)ATAN2= 1./2.*PI
         ELSE
            ATAN2 = ATAN(TRIP(2,2)/TRIP(2,1))
         ENDIF 
         
         TRIP(1,1) = MODU1**TRIPOW*COS(ATAN1/3.0)	
         TRIP(1,2) = MODU1**TRIPOW*SIN(ATAN1/3.0)
         TRIP(2,1) = MODU2**TRIPOW*COS(ATAN2/3.0)
         TRIP(2,2) = MODU2**TRIPOW*SIN(ATAN2/3.0)
         
      ENDIF
      
      EIVAL(1,1)= TRIP(1,1)+TRIP(2,1)-Pten/3.0 
      EIVAL(1,2)= TRIP(1,2)+TRIP(2,2)
      
      EIVAL(2,1)= -0.5*( trip(1,1)+trip(2,1) )
     $     -0.5*radtre*( trip(1,2)-trip(2,2) )
     $     -Pten/3.0
      EIVAL(2,2)= -0.5*( trip(1,2)+trip(2,2) )
     $       +0.5*radtre*( trip(1,1)-trip(2,1) )
      EIVAL(3,1)= -0.5*( trip(1,1)+trip(2,1) )
     $     +0.5*radtre*( trip(1,2)-trip(2,2) )
     $     -Pten/3.0
      EIVAL(3,2)= -0.5*( trip(1,2)+trip(2,2) )
     $     -0.5*radtre*( trip(1,1)-trip(2,1) )
      
      if(neigen.eq.1)then
         evalue1 =  EIVAL(1,1)+EIVAL(2,1)+EIVAL(3,1)
     $        - MAX(EIVAL(1,1),EIVAL(2,1),EIVAL(3,1))
     $        - MIN(EIVAL(1,1),EIVAL(2,1),EIVAL(3,1))
         evalue2 = 0.      
      elseif(neigen.eq.2)then
         evalue1 = EIVAL(2,2)*EIVAL(2,2) !lamba_ci
         if(abs(EIVAL(2,2)).gt.0.)then
            evalue2 = EIVAL(2,1)/EIVAL(2,2)
c     $           *EIVAL(2,2))
         else
            evalue2 = 0.
         endif
c         if(evalue2.gt.1.d5)write(6,*)'spiralling',evalue2,EIVAL(2,1),
c     $        EIVAL(2,2)
      endif

      END
      
      subroutine  outeigen(neigen,eig1,eig2,time)
      include 'parade_new.f'

      dimension eig1(m1,m2,m3),eig2(m1,m2,m3) 
      character*16 namfile
      character*5 ipfi
      
      itime=nint(time*1000.)
      write(ipfi,82)itime
 82   format(i5.5)
      if(neigen.eq.1)then
         namfile='lam2-'//ipfi//'.dat'
      elseif(neigen.eq.2)then
         namfile='swirl'//ipfi//'.dat'
      endif
      open(19,file=namfile,form='unformatted',status='unknown')
      rewind(19)
      if(neigen.eq.1)then
      write(19) n3,n2,n1,1
         write(19)
     %        (((eig1(i,j,k),k=1,n3),j=1,n2),i=1,n1)
      elseif(neigen.eq.2)then
      write(19) n3,n2,n1,1
         write(19)
     %        (((eig1(i,j,k),k=1,n3),j=1,n2),i=1,n1)
c     %        ,(((eig2(i,j,k),k=1,n3),j=1,n2),i=1,n1)
      endif
      close(19)

      end


      SUBROUTINE MAXMIN(IARRAY, NSIZE, MAX, MIN)
      DIMENSION IARRAY(NSIZE)
      MAX=IARRAY(1)
      MIN=IARRAY(1)
      DO 10 J=2,10
        IF (MAX.LT.IARRAY(J)) MAX=IARRAY(J) 
        IF (MIN.GT.IARRAY(J)) MIN=IARRAY(J)
10    CONTINUE
      RETURN
      END


!--------------------------------------------------------------------
       subroutine readgeo3ddim(nb,namefile,pitch)
!--------------------------------------------------------------------
!  scope: Read Input Geometry (STL) dimension
!--------------------------------------------------------------------
c       implicit none
       character(50)  namefile,string
       integer nb
       real pitch
       open(11,file=namefile)
       nb = 0
 10    read(11,*,err=98)string
       if (string.eq.'solid') goto 10
       if (string.eq.'endsolid') goto 99
       if (string.eq.'facet') goto 10
       if (string.eq.'endfacet') goto 10
       if (string.eq.'outer') then
          nb = nb+1
          read(11,*)string
          read(11,*)string
          read(11,*)string
          goto 10
       endif
       if (string.eq.'endloop') goto 10
 98    write(*,*)' Error in STL file ',string,nb+1
 99    close(11)
       if (pitch.ne.0) nb = 3*nb
       write(*,*)' STL with ',nb,' surface triangles '
!
       return
       end

!--------------------------------------------------------------------
      subroutine readgeo3d(time,xyzb,norb,nb,
     $     x1,x2,y1,y2,z1,z2,namefile,pitch)
c     subroutine readgeo3d(nb,x1,x2,y1,y2,z1,z2,namefile,pitch)
!--------------------------------------------------------------------
!  scope: Read Input Geometry (STL)
!  input: -
!  output: xyzc Input traingle Vertices (nb)
!--------------------------------------------------------------------
      include 'parade_new.f'
      
c      real x1,x2,y1,y2,z1,z2,small,pi,pitch
c      real v11,v12,v13,v21,v22,v23,v31,v32,v33,sp,area,d32,d21,d13,
c     $     sup_tot
c       integer j,nb,i,nbtmp,ntmp,isg
       real*8 norb(17000,3)
       real*8 xyzb(17000,9)
       character*50  namefile,string,string1
c       real*4  time,x3b_0,x2b_0,x1b_0,x3b_0_in,x2b_0_in,
c     %                angle,omegak,osci,phi,ampp,sigma
c       real*4 xcen,ycen,angtime,radpun,angpun,xpun,ypun
!
       x1 = 10000.
       x2 =-10000.
       y1 = 10000.
       y2 =-10000.
       z1 = 10000.
       z2 =-10000.

       pi=2.0*asin(1.0)
    
       xcen = x3b_0
       ycen = x2b_0
       angtime = angle
       write(*,*)'readgeo',time,xcen,ycen,angtime*180/pi,sigma
!
!      read geometry and compute boundaing box
!
       open(11,file=namefile)
       j = 0
       sup_tot = 0.
 10    read(11,*,err=98)string
c       if (string.eq.'solid') goto 10
c       if (j.eq.nb+1)goto 99
       if (string.eq.'solid')then
          j = j+1
          read(11,*)string,string1,(norb(j,i),i=1,3) 
          goto 10
       endif
       if (string.eq.'endsolid') goto 99
       if (string.eq.'facet') goto 10
c     if (string.eq.'endfacet') goto 10
       if (string.eq.'endfacet')then
          j = j+1
          if (j.eq.nb+1)goto 99
          read(11,*)string,string1,(norb(j,i),i=1,3) 
          goto 10
       endif
       if (string.eq.'outer') then
c          j = j+1
          read(11,*)string,(xyzb(j,i),i=1,3)
          xpun=xyzb(j,1)-x3b_0_in
          ypun=xyzb(j,2)-x2b_0_in
          radpun=sqrt(xpun**2.+ypun**2.)
          angpun=atan(ypun/xpun)
          if (xpun.lt.0.0) angpun=angpun+pi
          xyzb(j,1)=radpun*cos(angpun+angtime)+xcen
          xyzb(j,2)=radpun*sin(angpun+angtime)+ycen
          
          read(11,*)string,(xyzb(j,i),i=4,6)
          xpun=xyzb(j,4)-x3b_0_in
          ypun=xyzb(j,5)-x2b_0_in
          radpun=sqrt(xpun**2.+ypun**2.)
          angpun=atan(ypun/xpun)
          if (xpun.lt.0.0) angpun=angpun+pi
          xyzb(j,4)=radpun*cos(angpun+angtime)+xcen
          xyzb(j,5)=radpun*sin(angpun+angtime)+ycen

          read(11,*)string,(xyzb(j,i),i=7,9)
          xpun=xyzb(j,7)-x3b_0_in
          ypun=xyzb(j,8)-x2b_0_in
          radpun=sqrt(xpun**2.+ypun**2.)
          angpun=atan(ypun/xpun)
          if (xpun.lt.0.0) angpun=angpun+pi
          xyzb(j,7)=radpun*cos(angpun+angtime)+xcen
          xyzb(j,8)=radpun*sin(angpun+angtime)+ycen
c          v11 = xyzb(j,1)
c          v12 = xyzb(j,2)
c          v13 = xyzb(j,3)
c          v21 = xyzb(j,4)
c          v22 = xyzb(j,5)
c          v23 = xyzb(j,6)
c          v31 = xyzb(j,7)
c          v32 = xyzb(j,8)
c          v33 = xyzb(j,9)
c
c          d32=sqrt((v33-v23)**2+(v32-v22)**2+(v31-v21)**2)
c          d21=sqrt((v23-v13)**2+(v22-v12)**2+(v21-v11)**2)
c          d13=sqrt((v13-v33)**2+(v12-v32)**2+(v11-v31)**2)
cc     if (itr.le.5) write(*,*) '32',d32
cc     if (itr.le.5) write(*,*) '21',d21
cc     if (itr.le.5) write(*,*) '13',d13
c          sp=(d32+d21+d13)/2.
c          area=sqrt(sp*(sp-d32)*(sp-d21)*(sp-d13))
          x1 = min(x1,xyzb(j,1),xyzb(j,4),xyzb(j,7))
          x2 = max(x2,xyzb(j,1),xyzb(j,4),xyzb(j,7))
          y1 = min(y1,xyzb(j,2),xyzb(j,5),xyzb(j,8))
          y2 = max(y2,xyzb(j,2),xyzb(j,5),xyzb(j,8))
          z1 = min(z1,xyzb(j,3),xyzb(j,6),xyzb(j,9))
          z2 = max(z2,xyzb(j,3),xyzb(j,6),xyzb(j,9))
          
c          sup_tot = sup_tot + area

c          if (j.le.5.or.j.gt.3960) write(*,*)'sp =',sp,'area=',
c     $         area,'sup_tot',sup_tot
c            if (j.le.5.or.j.gt.3960) write(*,*)'vv',j,(xyzb(j,i),i=1,9)
          goto 10
       endif
       if (string.eq.'endloop') goto 10
 98    write(*,*)' Error in STL file '
 99    close(11)

       if (pitch.eq.0..and.j.ne.nb  ) write(*,*)' Error in Reading 
     $      STL file '
       if (pitch.ne.0..and.j.ne.nb/3) write(*,*)' Error in Reading 
     $      STL file '

c       if (pitch.ne.0) then
c          do j=1,nb/3
c             xyzb(   nb/3 +j,1:9) = xyzb(j,1:9)
c             xyzb(2*(nb/3)+j,1:9) = xyzb(j,1:9)
c          enddo
c       endif
!

c       write(*,*)' readgeom3d - sup_tot = ',sup_tot

       return
       end
c                                                                       
c***********************************************************************
c
      subroutine outpf1(time,vo1,vo2,vo3,voc1,voc2,voc3)
      include 'parade_new.f'

      dimension vo1(m1,m2,m3),vo2(m1,m2,m3),vo3(m1,m2,m3)
      dimension vo321(m1,m2,m3)
      dimension voc1(m1,m2,m3),voc2(m1,m2,m3),voc3(m1,m2,m3)
      character*70 namfi3
      character*5 ipfi

      itime=nint(1000.*time)
      write(ipfi,99)itime
      write(*,*) 'ITIME ',itime
 99   format(i5.5)
      
      write(6,*)' 1 se vuoi stampare il campo di velocita'' '
      read(5,*)iprintV
      if(iprintV.eq.1)then
         namfi3='velocita'//ipfi//'.dat'
         write(6,201) namfi3
 201     format(10x,'scrivo su ',a70)
         n3pp=(n3-1)/n3p+1  
         n2pp=(n2-1)/n2p+1 
         open(59,file=namfi3,form='unformatted',status='unknown')
         rewind(59)
         if(refsys.eq.1) then
            write(6,*)' x10,x20,x30 ',x10,x20,x30
            write(6,*)'wx= ',wx,'u0x= ',u0x,'du0x= ',du0x
            write(6,*)'u0z= ',u0z,'du0z= ',du0z
            write(59) n3,n2,n1,3
            sinalph = sin(alphax)
            cosalph = cos(alphax)
            write(59) 
     %           (((voc2(i,j,k)*sinalph+voc3(i,j,k)*cosalph,
     $           k=1,n3),j=1,n2),i=1,n1),
     %           (((voc2(i,j,k)*cosalph-voc3(i,j,k)*sinalph,
     $           k=1,n3),j=1,n2),i=1,n1),
     %           (((voc1(i,j,k),k=1,n3),j=1,n2),i=1,n1)
         else
            write(59) 
     %           (((voc3(i,j,k),k=1,n3),j=1,n2),i=1,n1),
     %           (((voc2(i,j,k),k=1,n3),j=1,n2),i=1,n1),
     %           (((voc1(i,j,k),k=1,n3),j=1,n2),i=1,n1)
         endif
         close(59)
      endif
      
      if (nvort.eq.1)then 
         do i=1,n1
            do k=1,n3
               do j=1,n2
                  vo1(i,j,k)=vo1(i,j,k)-wx
                  vo321(i,j,k)=sqrt(vo1(i,j,k)**2.+vo2(i,j,k)**2.
     $                 +vo3(i,j,k)**2.)
               enddo
            enddo
         enddo 

         namfi3='vorticita'//ipfi//'.dat'
         write(6,201) namfi3
         n3pp=(n3-1)/n3p+1  
         n2pp=(n2-1)/n2p+1 
         open(59,file=namfi3,form='unformatted',status='unknown')
         rewind(59)
         write(59) n3,n2,n1,3
         write(59) 
     %        (((vo3(i,j,k),k=1,n3),j=1,n2),i=1,n1)
     %        ,(((vo2(i,j,k),k=1,n3),j=1,n2),i=1,n1)  
     %        ,(((vo1(i,j,k),k=1,n3),j=1,n2),i=1,n1) 
         close(59)
      
         namfi3='intvortic'//ipfi//'.dat'
         write(6,201) namfi3
         open(59,file=namfi3,form='unformatted',status='unknown')
         rewind(59)
         write(59) n3,n2,n1,1
         write(59) 
     %        (((vo321(i,j,k),k=1,n3),j=1,n2),i=1,n1)
         close(59)
         
      endif

      return
      end

!--------------------------------------------------------------------
