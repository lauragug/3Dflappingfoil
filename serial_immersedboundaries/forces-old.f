c 
c***********************************************************************
c
      subroutine forces(time,q1,q2,q3,pr)  
      include 'param.f'

      REAL q3(m1,m2,m3)
      REAL q2(m1,m2,m3)
      REAL q1(m1,m2,m3)
      REAL pr(m1,m2,m3)
c      real norb(17000,3)          
c      real xyzb(17000,9) 
c      dimension tauw(55000),tauw1(55000),tauw2(55000),tauw3(55000)
c     &          ,prsw(55000)
      real vett32(3),vett21(3),vettnorm(3),vettrif(3),versnorm(3)

      real risul3_v,risul2_v,risul1_v,risul3_p,risul2_p,risul1_p,
     $     risulM_x,sup_tot   
c      real dim3_min,dim3_max,dim2_min,dim2_max,dim1_min,dim1_max
      real v13,v12,v11,v23,v22,v21,v33,v32,v31,dist32,dist21,dist13,sp,
     $     area,bar3,bar2,bar1,vmod,prsc
      real delta,delta3,delta2,delta1,distdb3,distdb2,distdb1
      real tauvisc1,tauvisc2,tauvisc3,tauviscmod,forza1_v,forza2_v,
     $     forza3_v,forza1_p,forza2_p,forza3_p,momento_x
      real q1int,q1int1,q1int2,q1int3,q1int4,q2int,q2int1,q2int2,q2int3,
     $     q2int4,q3int,q3int1,q3int2,q3int3,q3int4,prsint,prsint1,
     $     prsint2,prsint3,prsint4
      real inv_delta,vrel1,vrel2,vrel3,vreln
      real cx_v,cy_v,cz_v,cx_p,cy_p,cz_p,cM_x
      real time
c     $     pitch

      integer itr,ntau,kint,jint,iint,k,j,i,km,jm,im,ks,js,is
c     ,nb
      character(50) namfile
      character*20 namfi3
      character*5 ipfi

      if(amod(time,tprint).lt.dt) then                 
         itime=nint(1000.*time)
         write(ipfi,99) itime
 99      format(i5.5)
         namfi3='prof'//ipfi//'.dat'
         write(6,201)time,namfi3
 201     format(10x,'At t=',e10.3,' write data on ',a20)
         open(99,file=namfi3,form='formatted')
      endif          

      risul3_v=0.
      risul2_v=0.
      risul1_v=0.
      risul3_p=0.
      risul2_p=0.
      risul1_p=0.
      risulM_x=0.
      sup_tot=0.0
c      namfile = 'n30-370dT.stl'
c      write(*,*)'read the geometry from file ',namfile
c      pitch1 = 0.
c      call readgeo3ddim(nb,namfile,pitch1)
c      call readgeo3d(time,xyzb,norb,nb,dim3_min,dim3_max,dim2_min,
c     $     dim2_max,dim1_min,dim1_max,namfile,pitch1)
      write(*,*)' triangles which describe the surface = ',nb
      do itr=1,nb
c     if (norb(itr,3).gt.0.5)then 
c         write(33,*)' tr',itr,
c     $        'normale x,y,z = ',norb(itr,1),norb(itr,2),norb(itr,3)
c     goto 987
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
         dist32=sqrt((v33-v23)**2+(v32-v22)**2+(v31-v21)**2)
         dist21=sqrt((v23-v13)**2+(v22-v12)**2+(v21-v11)**2)
         dist13=sqrt((v13-v33)**2+(v12-v32)**2+(v11-v31)**2)
         sp = (dist32+dist21+dist13)/2.
         area = sqrt(sp*(sp-dist32)*(sp-dist21)*(sp-dist13))
         bar3=(v33+v23+v13)/3.
         bar2=(v32+v22+v12)/3.
         bar1=(v31+v21+v11)/3.
         
c     calcolo del vettore normale dalle coordinate vertici
c         vett32(3)=v33-v23
c         vett32(2)=v32-v22
c         vett32(1)=v31-v21
c         vett21(3)=v23-v13
c         vett21(2)=v22-v12
c         vett21(1)=v21-v11
cc     calcolo analitico del vettore normale al triangolo       
c         vettnorm(3)=  vett32(1)*vett21(2)-vett21(1)*vett32(2)
c         vettnorm(2)=-(vett32(1)*vett21(3)-vett21(1)*vett32(3))
c         vettnorm(1)=  vett32(2)*vett21(3)-vett21(2)*vett32(3)
         
         vettnorm(3) = norb(itr,1) !!!asse x stl e' asse x3
         vettnorm(2) = norb(itr,2)
         vettnorm(1) = norb(itr,3) !!!asse z stl e' asse x1
         vmod=sqrt(vettnorm(3)**2+vettnorm(2)**2+vettnorm(1)**2)
c         if (vmod.le.0.001.or.vettnorm(1).gt.0.05) then               
c            write(*,*)' tr',itr,
c     $           'normale x,y,z = ',norb(itr,1),norb(itr,2),norb(itr,3)
c            goto 987
c         endif
         vettrif(3)=0.25-bar3       ! PARTE DELICATA
         vettrif(2)=0.-bar2     ! ADATTARE AD OGNI CORPO
         vettrif(1)=0.-bar1     !
         prsc=vettnorm(3)*vettrif(3)+vettnorm(2)*vettrif(2)
     %        +vettnorm(1)*vettrif(1)
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
         delta=0.
         do k=1,n3m
            if(bar3.gt.x3c(k).and.bar3.le.x3c(k+1)) then
               do j=1,n2m
                  if(bar2.gt.x2c(j).and.bar2.le.x2c(j+1)) then
                     do i=1,n1m-1  !PROVA ELIMIN PUNTI SPARSI
                        if(bar1.gt.x1c(i).and.bar1.le.x1c(i+1)) then
                           delta=1.1*sqrt((x1c(i)-x1c(i+1))**2.+
     &                          (x2c(j)-x2c(j+1))**2.+
     &                          (x3c(k)-x3c(k+1))**2.)
                        endif
                     end do
                  end if
               end do
            end if       
         end do
         
         inv_delta = 1./(delta*ren)
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
         tauviscmod=0.          !inizializzate nel caso non vengano calcolate 
         tauvisc1=0.            !dentro ciclo do - variabili che escono
         tauvisc2=0.
         tauvisc3=0.
         forza1_v=0.
         forza2_v=0.
         forza3_v=0.
         forza1_p=0.
         forza2_p=0.
         forza3_p=0.
         momento_x=0.
         do k=4,n3m-1 
            if(delta3.gt.x3c(k-1).and.delta3.le.x3c(k)) then
               do j=2,n2m
                  if(delta2.gt.x2c(j-1).and.delta2.le.x2c(j)) then
                     do i=2,n1m !PROVA ELIMIN PUNTI SPARSI
                        if(delta1.gt.x1c(i-1).and.delta1.le.x1c(i))then
                           km=kmv(k)
                           jm=jmv(j)
                           im=imv(i)
                           ks=km+kint
                           js=jm+jint
                           is=im+iint
C     INTERPOLAZIONE PRIMA COMPONENTE 
                           q1int1=(q1(i ,jm,km)  * (delta1-x1c(im))
     &                          -q1(im,jm,km) * (delta1-x1c(i )))
     $                          /(x1c(i) - x1c(im))
                           
                           q1int2=(q1(i ,js,km)  * (delta1-x1c(im))
     &                          -q1(im,js,km) * (delta1-x1c(i )))
     $                          /(x1c(i) - x1c(im))
                           
                           q1int3=(q1int2        * (delta2-x2m(jm))
     &                          -q1int1        * (delta2-x2m(js)))/
     &                          (x2m(js)-x2m(jm))
                           
                           q1int1=(q1(i ,jm,ks)  * (delta1-x1c(im))
     &                          -q1(im,jm,ks) * (delta1-x1c(i )))
     $                          /(x1c(i) - x1c(im))
                           
                           q1int2=(q1(i ,js,ks)  * (delta1-x1c(im))
     &                          -q1(im,js,ks) * (delta1-x1c(i )))
     $                          /(x1c(i) - x1c(im))

                           q1int4=(q1int1        * (delta2-x2m(js))
     &                          -q1int2        * (delta2-x2m(jm)))/
     &                          (x2m(jm)-x2m(js))
                           
                           q1int =(q1int4        * (delta3-x3m(km))
     &                          -q1int3        * (delta3-x3m(ks)))/
     &                          (x3m(ks)-x3m(km))
                           
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
                           
                           q2int2=(q2(is,j ,ks)  * (delta2-x2c(jm))
     &                          -q2(is,jm,ks)  * (delta2-x2c(j )))/
     &                          (x2c(j)-x2c(jm))
                           
                           q2int4=(q2int1        * (delta3-x3m(ks))
     &                          -q2int2        * (delta3-x3m(km)))/
     &                          (x3m(km)-x3m(ks))
                           
                           q2int =(q2int4        * (delta1-x1m(im))
     &                          -q2int3        * (delta1-x1m(is)))/
     &                          (x1m(is)-x1m(im))
                           
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
                           
                           q3int2=(q3(is,js,k )  * (delta3-x3c(km))
     &                          -q3(is,js,km)  * (delta3-x3c(k )))/
     &                          (x3c(k)-x3c(km))
                           
                           q3int4=(q3int1        * (delta1-x1m(is))
     &                          -q3int2        * (delta1-x1m(im)))/
     &                          (x1m(im)-x1m(is))
                           
                           q3int =(q3int4        * (delta2-x2m(jm))
     &                          -q3int3        * (delta2-x2m(js)))/
     &                          (x2m(js)-x2m(jm))
                           
C     INTERPOLAZIONE PRESSIONE
                           prsint1=(pr(im,jm,ks)  * (delta3-x3m(km))
     &                          - pr(im,jm,km)  * (delta3-x3m(ks)))/
     &                          (x3m(ks)-x3m(km))
                           
                           prsint2=(pr(is,jm,ks)  * (delta3-x3m(km))
     &                          - pr(is,jm,km)  * (delta3-x3m(ks)))/
     &                          (x3m(ks)-x3m(km))
                           
                           prsint3=(prsint2        * (delta1-x1m(im))
     &                          - prsint1        * (delta1-x1m(is)))/
     &                          (x1m(is)-x1m(im))
                           
                           prsint1=(pr(im,js,ks)  * (delta3-x3m(km))
     &                          - pr(im,js,km)  * (delta3-x3m(ks)))/
     &                          (x3m(ks)-x3m(km))
                           
                           prsint2=(pr(is,js,ks)  * (delta3-x3m(km))
     &                          - pr(is,js,km)  * (delta3-x3m(ks)))/
     &                          (x3m(ks)-x3m(km))
                           
                           prsint4=(prsint1        * (delta1-x1m(is))
     &                          - prsint2        * (delta1-x1m(im)))/
     &                          (x1m(im)-x1m(is))

                           prsint =(prsint4        * (delta2-x2m(jm))
     &                          -prsint3        * (delta2-x2m(js)))/
     &                          (x2m(js)-x2m(jm))
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
c                           prsint2=(pr(is,js,ks)  * (bar3-x3m(km))
c     &                          - pr(is,js,km)  * (bar3-x3m(ks)))/
c     &                          (x3m(ks)-x3m(km))
c                           
c                           prsint4=(prsint1        * (bar1-x1m(is))
c     &                          - prsint2        * (bar1-x1m(im)))/
c     &                          (x1m(im)-x1m(is))
c                           
c                           prsint =(prsint4        * (bar2-x2m(jm))
c     &                          -prsint3        * (bar2-x2m(js)))/
c     &                          (x2m(js)-x2m(jm))
                           gradpr = duhy*cos(alphax)*versnorm(2)-
     $                          duhy*sin(alphax)*versnorm(3)
     $           +dwx*(-versnorm(2)*(bar3-x30)+versnorm(3)*(bar2-x20))
     $           -wx*wx*(versnorm(2)*(bar2-x20)+versnorm(3)*(bar3-x30))
                           prsurf = prsint + gradpr*delta
CC     PROIEZIONE DI QINT LUNGO N     
C                           qintn = q1int*versnorm(1)+q2int*versnorm(2)+
C     &                          q3int*versnorm(3)
CC     VETTORE TANGENTE ALLA SUPERFICIE
C                           qtang1 = q1int-qintn*versnorm(1)
C                           qtang2 = q2int-qintn*versnorm(2)
C                           qtang3 = q3int-qintn*versnorm(3) 
C     velocita' del corpo
                           vsurf1 = 0.
                           vsurf2 = uhy*cos(alphax)-wx*(delta3-x30)
                           vsurf3 = - uhy*sin(alphax)+wx*(delta2-x20)
C                           vsurfn = vsurf1*versnorm(1)
C     $                          +vsurf2*versnorm(2)+vsurf3*versnorm(3)
C                           vsurf1_tang = vsurf1 - vsurfn*versnorm(1)
C                           vsurf2_tang = vsurf2 - vsurfn*versnorm(2)
C                           vsurf3_tang = vsurf3 - vsurfn*versnorm(3)
CC     SFORZO TANGENZIALE
C                           tauvisc1=(qtang1 - vsurf1_tang)/delta/ren
C                           tauvisc2=(qtang2 - vsurf2_tang)/delta/ren
C                           tauvisc3=(qtang3 - vsurf3_tang)/delta/ren
                           vrel1 = q1int - vsurf1
                           vrel2 = q2int - vsurf2
                           vrel3 = q3int - vsurf3
                           vreln = vrel1*versnorm(1)
     $                          +vrel2*versnorm(2)+vrel3*versnorm(3)
                           tauvisc1=(vrel1 - vreln*versnorm(1))
     $                          *inv_delta
                           tauvisc2=(vrel2 - vreln*versnorm(3))
     $                          *inv_delta
                           tauvisc3=(vrel3 - vreln*versnorm(3))
     $                          *inv_delta
                           
                           forza1_v=tauvisc1*area
                           forza2_v=tauvisc2*area
                           forza3_v=tauvisc3*area
                           forza1_p=(-prsurf*versnorm(1))*area
                           forza2_p=(-prsurf*versnorm(2))*area
                           forza3_p=(-prsurf*versnorm(3))*area
                           momento_x = 
     $                          (forza3_p+forza3_v)*(delta2-x20)
     $                          - (forza2_p+forza2_v)*(delta3-x30)
c                           write(44,42)itr,bar1,bar2,bar3,
c     $                          vrel1,vrel2,vrel3,
c     $                          forza1_p,forza2_p,forza3_p
c                           write(41,42)itr,bar1,bar2,bar3,
c     $                          q1int,q2int,q3int
c     $                          ,vsurf1,vsurf2,vsurf3
c                           write(45,43)itr,bar1,bar2,bar3,
c     $                          delta1,delta2,delta3,delta
c                           write(46,43)itr,bar1,bar2,bar3,
c     $                          forza1_p,forza2_p,forza3_p
c                           write(47,43)itr,bar1,bar2,bar3,versnorm(1),
c     $                          versnorm(2),versnorm(3)
                           if(amod(time,tprint).lt.dt) then
                              write(99,42)itr,bar1,bar2,bar3,
     $             vrel2*sqrt(vrel3**2.+vrel2**2.)/abs(vrel2)/delta, 
     $             vrel3*sqrt(vrel3**2.+vrel2**2.)/abs(vrel3)/delta,
     $                          prsurf,prsint,gradpr
                           endif
                           tauviscmod= sqrt(tauvisc1**2.+tauvisc2**2.
     $                          +tauvisc3**2.)
                        endif
                     end do
                  end if
               end do
            end if       
         end do
 42      format(i4,10(2x,e12.5))
 43      format(i4,8(2x,e12.5))
c         do m=1,3
c            ntau=int(itr*3-3+m)
c            tauw(ntau)=tauviscmod
c            tauw1(ntau)=tauvisc1
c            tauw2(ntau)=tauvisc2
c            tauw3(ntau)=tauvisc3
c            prsw(ntau)=prsurf
c         enddo
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
         if (tauviscmod.gt.0.)then
            sup_tot = sup_tot + area
c            write(99,1002)itr,area,sup_tot,bar3,bar2,bar1
c         else
c            write(66,1002)itr,area,sup_tot,bar3,bar2,bar1
c            if (abs(bar1).gt.0.2)write(*,1002)itr,area,sup_tot,
c     $           bar3,bar2,bar1
         endif
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
      
c      write(*,*)'forza risultante'
c      write(*,*)'3 comp visc, press',risul3_v,risul3_p
c      write(*,*)'2 comp visc, press',risul2_v,risul2_p    
c      write(*,*)'1 comp visc, press',risul1_v,risul1_p
c      write(*,*)'3 risultante_x',(risul3_v+risul3_p)
c      write(*,*)'2 risultante_y',(risul2_v+risul2_p)
c      write(*,*)'1 risultante_z',(risul1_v+risul1_p)
      risult3_hat = (risul2_v+risul2_p)*sin(alphax) +
     $     (risul3_v+risul3_p)*cos(alphax)
      risult2_hat = (risul2_v+risul2_p)*cos(alphax) - 
     $     (risul3_v+risul3_p)*sin(alphax)
c      write(*,*)'3 risultante_xhat',risult3_hat
c      write(*,*)'2 risultante_yhat',risult2_hat
c      write(*,*)' sup_tot, ne ',sup_tot,nb
      write(65,1000) time,risul1_v,risul1_p,risul2_v,risul2_p,
     $     risul3_v,risul3_p,(risul1_v+risul1_p),risult2_hat,
     $     risult3_hat,risulM_x
      cx_v = 20.0*risul3_v
      cy_v = 20.0*risul2_v
      cz_v = 20.0*risul1_v
      cx_p = 20.0*risul3_p
      cy_p = 20.0*risul2_p
      cz_p = 20.0*risul1_p
      cM_x = 20.0*risulM_x
      write(66,1001) time,cx_v,cy_v,cz_v,cx_p,cy_p,cz_p,cM_x
c      write(*,*)'coefficienti'
c      write(*,*)'cx v - p ',cx_v,cx_p
c      write(*,*)'cy v - p ',cy_v,cy_p 
c      write(*,*)'cz v - p ',cz_v,cz_p
 1000 format(11(2x,e12.5))
 1001 format(8(2x,e12.5))
 1002 format(i5,5(2x,e12.5))
c      namfile='norm.tec'
c      open(31,file=namfile,form='formatted',status='unknown')
c      write(31,*) 'Variables = X, Y, Z, Norm1, Norm2, Norm3'
c      write(31,*)'ZONE N = ',3*nb,',E = ',nb,',F=FEPOINT,ET=TRIANGLE'
c      write(31,*)'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)'
c      do k = 1,nb                 
c         write(31,*)xyzb(k,1),xyzb(k,2),xyzb(k,3), 
c     $        norb(k,1),norb(k,2),norb(k,3)
c         write(31,*)xyzb(k,4),xyzb(k,5),xyzb(k,6),
c     $        norb(k,1),norb(k,2),norb(k,3)
c         write(31,*)xyzb(k,7),xyzb(k,8),xyzb(k,9),
c     $        norb(k,1),norb(k,2),norb(k,3)
c      enddo
c      j = 1
c      do k = 1,nb
c         write(31,*)j,j+1,j+2
c         j=j+3
c      enddo
c      close(31)
c      namfile='tau.tec'
cc     c     endif
c      open(31,file=namfile,form='formatted',status='unknown')
c      write(31,*) 'Variables = X, Y, Z, Tau(module), TauX, TauY,
c     &     TauZ,Pressure'
c      write(31,*)'ZONE N = ',3*nb,',E = ',nb,',F=FEPOINT,ET=TRIANGLE'
cc     write(31,*)'ZONE N = ',3*nb,',E = ',nb,',F=FEPOINT, ET=TRIANGLE'
c      write(31,*)'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE
c     &     SINGLE SINGLE)'
c      do k = 1,nb                 
c         write(31,*)xyzb(k,1),xyzb(k,2),xyzb(k,3),
c     $        tauw(k*3-3+1),tauw3(k*3-3+1),tauw2(k*3-3+1),
c     $        tauw1(k*3-3+1),prsw(k*3-3+1)
c         write(31,*)xyzb(k,4),xyzb(k,5),xyzb(k,6),
c     $        tauw(k*3-3+2),tauw3(k*3-3+2),tauw2(k*3-3+2),
c     $        tauw1(k*3-3+2),prsw(k*3-3+2)
c         write(31,*)xyzb(k,7),xyzb(k,8),xyzb(k,9),
c     $        tauw(k*3-3+3),tauw3(k*3-3+3),tauw2(k*3-3+3),
c     $        tauw1(k*3-3+3),prsw(k*3-3+3)
c      enddo
c      j = 1
c      do k = 1,nb
c         write(31,*)j,j+1,j+2
c         j=j+3
c      enddo
c      close(31)
         
      return
      end

!--------------------------------------------------------------------
       subroutine readgeo3ddim(nb,namefile,pitch)
!--------------------------------------------------------------------
!  scope: Read Input Geometry (STL) dimension
!--------------------------------------------------------------------
       implicit none
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
      subroutine readgeo3d(time,xyzb,norb,nb,dim3_min,dim3_max,dim2_min,
     $     dim2_max,dim1_min,dim1_max,namefile,pitch)
c     subroutine readgeo3d(nb,x1,x2,y1,y2,z1,z2,namefile,pitch)
!--------------------------------------------------------------------
!  scope: Read Input Geometry (STL)
!  input: -
!  output: xyzc Input traingle Vertices (nb)
!--------------------------------------------------------------------
c       include 'param.f'
      real time,xcen,ycen,angtime,radpun,angpun,xpun,ypun
      real dim3_min,dim3_max,dim2_min,dim2_max,dim1_min,dim1_max
      real small,pi
c     $      ,pitch
      real v11,v12,v13,v21,v22,v23,v31,v32,v33,sp,area,d32,d21,d13,
     $     sup_tot
      real norb(10000,3)
      real xyzb(10000,9)
      integer j,i,nbtmp,ntmp,isg,nb
      character(50) namefile,string,string1
!     
      dim3_min = 10000.
      dim3_max = -10000.
      dim2_min = 10000.
      dim2_max = -10000.
      dim1_min = 10000.
      dim1_max = -10000.
      
      pi=2.0*asin(1.0)
      
c     xcen = x30
c     ycen = hy
c     angtime = alphax
      xcen = 0.
      ycen = 0.
      angtime = 0.
c     write(*,*)'readgeo',time,xcen,ycen,angtime*180/pi
!     
!     read geometry and compute bounding box
!     modificata rispetto quella del code_grid per 
!     poter leggere anche le normali
      open(11,file=namefile)
      j = 0
      sup_tot = 0.
 10   read(11,*,err=98)string
c     if (string.eq.'solid') goto 10
c     if (j.eq.nb)goto 99
      if (string.eq.'solid')then
         j = j+1       
         read(11,*)string,string1,(norb(j,i),i=1,3) 
         goto 10
      endif
      if (string.eq.'endsolid') then
         goto 99
      endif
      if (string.eq.'facet') goto 10
c     if (string.eq.'endfacet') goto 10
      if (string.eq.'endfacet')then
         j = j+1
         if (j.gt.nb)then
c            write(*,*)'yes'
            goto 99
         endif
         read(11,*)string,string1,(norb(j,i),i=1,3) 
         goto 10
      endif
      if (string.eq.'outer') then
c     j = j+1
         read(11,*)string,(xyzb(j,i),i=1,3)
c     xpun=xyzb(j,1)-x30
c     ypun=xyzb(j,2)-x20
c     radpun=sqrt(xpun**2.+ypun**2.)
c          angpun=atan(ypun/xpun)
c          if (xpun.lt.0.0) angpun=angpun+pi
c          xyzb(j,1)=radpun*cos(angpun+angtime)+xcen
c          xyzb(j,2)=radpun*sin(angpun+angtime)+ycen
          
         read(11,*)string,(xyzb(j,i),i=4,6)
c          xpun=xyzb(j,4)-x30
c          ypun=xyzb(j,5)-x20
c          radpun=sqrt(xpun**2.+ypun**2.)
c          angpun=atan(ypun/xpun)
c          if (xpun.lt.0.0) angpun=angpun+pi
c          xyzb(j,4)=radpun*cos(angpun+angtime)+xcen
c          xyzb(j,5)=radpun*sin(angpun+angtime)+ycen
         
         read(11,*)string,(xyzb(j,i),i=7,9)
c          xpun=xyzb(j,7)-x30
c          ypun=xyzb(j,8)-x20
c          radpun=sqrt(xpun**2.+ypun**2.)
c          angpun=atan(ypun/xpun)
c          if (xpun.lt.0.0) angpun=angpun+pi
c          xyzb(j,7)=radpun*cos(angpun+angtime)+xcen
c          xyzb(j,8)=radpun*sin(angpun+angtime)+ycen
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
c     area=sqrt(sp*(sp-d32)*(sp-d21)*(sp-d13))
         dim3_min = min(dim3_min,xyzb(j,1),xyzb(j,4),xyzb(j,7))
         dim3_max = max(dim3_max,xyzb(j,1),xyzb(j,4),xyzb(j,7))
         dim2_min = min(dim3_min,xyzb(j,2),xyzb(j,5),xyzb(j,8))
         dim2_max = max(dim3_max,xyzb(j,2),xyzb(j,5),xyzb(j,8))
         dim1_min = min(dim3_min,xyzb(j,3),xyzb(j,6),xyzb(j,9))
         dim1_max = max(dim3_max,xyzb(j,3),xyzb(j,6),xyzb(j,9))
         
c          sup_tot = sup_tot + area

         goto 10
      endif
      if (string.eq.'endloop') goto 10
 98   write(*,*)' Error in STL file ',string,j
 99   close(11)

c       if (pitch.eq.0..and.j.ne.nb  ) write(*,*)' Error in Reading 
c     $      STL file ',pitch,j
c       if (pitch.ne.0..and.j.ne.nb/3) write(*,*)' Error in Reading 
c     $      STL file ',pitch,j

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
!--------------------------------------------------------------------
