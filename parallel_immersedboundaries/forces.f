c***********************************************************************
c
c      subroutine forces(time,q1,q2,q3,pr)  !PARALL
      subroutine forces(time,comm)  
      include 'param.f'
      include "mpif.h"       !PARALL   

c      dimension tauw(55000),tauw1(55000),tauw2(55000),tauw3(55000)
c&          ,prsw(55000)
      real vett32(3),vett21(3),vettnorm(3),vettrif(3),versnorm(3)

      real risul3_v,risul2_v,risul1_v,risul3_p,risul2_p,risul1_p,
     $     risulM1,risulM2,risulM3,sup_tot   
      real v13,v12,v11,v23,v22,v21,v33,v32,v31,dist32,dist21,dist13,sp,
     $     area,bar3,bar2,bar1,vmod,prsc
      real delta,delta3,delta2,delta1,distdb3,distdb2,distdb1
      real tauvisc1,tauvisc2,tauvisc3,tauviscmod,forza1_v,forza2_v,
     $    forza3_v,forza1_p,forza2_p,forza3_p,momento1,momento2,momento3
      real q1int,q1int1,q1int2,q1int3,q1int4,q2int,q2int1,q2int2,q2int3,
     $     q2int4,q3int,q3int1,q3int2,q3int3,q3int4,prsint,prsint1,
     $     prsint2,prsint3,prsint4
      real inv_delta,vrel1,vrel2,vrel3,vreln,inv_uzero
      real cx_v,cy_v,cz_v,cx_p,cy_p,cz_p,cM_x
      real time

      integer itr,ntau,kint,jint,iint,k,j,i,km,jm,im,ks,js,is
      integer comm,ierr    !PARALL      

      character(50) namfile
      character(20) namfi2,namfi3
      character(6)  ipfi

c      if(amod(time,tprint).lt.dt) then                 
c         itime=nint(1000.*time)
c         write(ipfi,99) itime
c 99      format(i6.6)
c         namfi2='tau'//ipfi//'.dat'
c         namfi3='prof'//ipfi//'.dat'
c         write(6,201)time,namfi3
c 201     format(10x,'At t=',e10.3,' write data on ',a20)
c         open(99,file=namfi3,form='formatted')
c      endif          

      risul3_v=0.
      risul2_v=0.
      risul1_v=0.
      risul3_p=0.
      risul2_p=0.
      risul1_p=0.
      risulM1=0.
      risulM2=0.
      risulM3=0.
      sup_tot=0.0
c      inv_uzero=1./(alx1*uzero*uzero) !CIL
c      if (uzero.eq.0.)inv_uzero=1./alx1   !CIL
c      namfile = 'n30-370dT.stl'
c      write(*,*)'read the geometry from file ',namfile
c      pitch1 = 0.
c      call readgeo3ddim(nb,namfile,pitch1)
c      call readgeo3d(time,xyzb,norb,nb,dim3_min,dim3_max,dim2_min,
c     $     dim2_max,dim1_min,dim1_max,namfile,pitch1)
c      write(*,*)' triangles which describe the surface = ',nb
      itr = 0
      itrcount = 0
      do itr=1,nb
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
         vettrif(3)=0.  - bar3    ! PARTE DELICATA
         vettrif(2)=0.  - bar2     ! ADATTARE AD OGNI CORPO
         vettrif(1)=0.  - bar1     !
         prsc=vettnorm(3)*vettrif(3)+vettnorm(2)*vettrif(2)
     %        +vettnorm(1)*vettrif(1)
         if(prsc.gt.0.)then
            vettnorm(3)=-vettnorm(3)
            vettnorm(2)=-vettnorm(2)
            vettnorm(1)=-vettnorm(1)    
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
c         do k=1,n3m
         do k = kb_min, kb_max
            if(bar3.gt.x3c(k).and.bar3.le.x3c(k+1)) then
c               do j=1,n2m
               do j = jb_min, jb_max
                  if(bar2.gt.x2c(j).and.bar2.le.x2c(j+1)) then
c                     do i=1,n1m-1  !PROVA ELIMIN PUNTI SPARSI
                     do i=1,n1m  ! PARALL
                      iglob = myid*n1m + i
                   if(bar1.gt.x1c(iglob).and.bar1.le.x1c(iglob+1))delta=
     $                     1.5*sqrt((x1c(iglob)-x1c(iglob+1))**2.+
     &                          (x2c(j)-x2c(j+1))**2.+
     &                          (x3c(k)-x3c(k+1))**2.)
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
c         do k=4,n3m-1 
         do k = kb_min, kb_max
            if(delta3.gt.x3c(k-1).and.delta3.le.x3c(k)) then
c               do j=2,n2m
               do j = jb_min,jb_max
                  if(delta2.gt.x2c(j-1).and.delta2.le.x2c(j)) then
C                     do i=2,n1m !PROVA ELIMIN PUNTI SPARSI
                     do i=2,n1m+1 ! PARALL

                 iglob = myid*n1m + i
                 if(delta1.gt.x1c(iglob-1).and.delta1.le.x1c(iglob))then
                           itrcount = itrcount + 1
                           km=kmv(k)
                           jm=jmv(j)
                           im=imv(i)
                           imglob=imv(iglob)! PARALL
                           ks=km+kint
                           js=jm+jint
                           is=im+iint
                           isglob=imglob+iint! PARALL
C     INTERPOLAZIONE PRIMA COMPONENTE 
                           q1int1=
     $                          (q1(i ,jm,km)  * (delta1-x1c(imglob))! PARALL
     &                          -q1(im,jm,km) * (delta1-x1c(iglob )))! PARALL
     $                          /(x1c(iglob) - x1c(imglob))! PARALL
                           
                           q1int2=(q1(i ,js,km)  * (delta1-x1c(imglob))! PARALL
     &                          -q1(im,js,km) * (delta1-x1c(iglob )))! PARALL
     $                          /(x1c(iglob) - x1c(imglob))! PARALL
                           
                           q1int3=(q1int2        * (delta2-x2m(jm))
     &                          -q1int1        * (delta2-x2m(js)))/
     &                          (x2m(js)-x2m(jm))
                           
                           q1int1=(q1(i ,jm,ks)  * (delta1-x1c(imglob))! PARALL
     &                          -q1(im,jm,ks) * (delta1-x1c(iglob )))! PARALL
     $                          /(x1c(iglob) - x1c(imglob))! PARALL
                           
                           q1int2=(q1(i ,js,ks)  * (delta1-x1c(imglob))! PARALL
     &                          -q1(im,js,ks) * (delta1-x1c(iglob )))! PARALL
     $                          /(x1c(iglob) - x1c(imglob))! PARALL

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
                           
                           q2int =(q2int4        * (delta1-x1m(imglob))! PARALL
     &                          -q2int3        * (delta1-x1m(isglob)))/! PARALL
     &                          (x1m(isglob)-x1m(imglob))! PARALL
                           
C     INTERPOLAZIONE TERZA COMPONENTE 
                           q3int1=(q3(im,jm,k )  * (delta3-x3c(km))
     &                          -q3(im,jm,km)  * (delta3-x3c(k )))/
     &                          (x3c(k)-x3c(km))
                           
                           q3int2=(q3(is,jm,k )  * (delta3-x3c(km))
     &                          -q3(is,jm,km)  * (delta3-x3c(k )))/
     &                          (x3c(k)-x3c(km))
                           
                           q3int3=(q3int2        * (delta1-x1m(imglob))! PARALL
     &                          -q3int1        * (delta1-x1m(isglob)))/! PARALL
     &                          (x1m(isglob)-x1m(imglob))! PARALL
                           
                           q3int1=(q3(im,js,k )  * (delta3-x3c(km))
     &                          -q3(im,js,km)  * (delta3-x3c(k )))/
     &                          (x3c(k)-x3c(km))
                           
                           q3int2=(q3(is,js,k )  * (delta3-x3c(km))
     &                          -q3(is,js,km)  * (delta3-x3c(k )))/
     &                          (x3c(k)-x3c(km))
                           
                           q3int4=(q3int1        * (delta1-x1m(isglob))! PARALL
     &                          -q3int2        * (delta1-x1m(imglob)))/! PARALL
     &                          (x1m(imglob)-x1m(isglob))! PARALL
                           
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
                           
                           prsint3=(prsint2     * (delta1-x1m(imglob))! PARALL
     &                          - prsint1       * (delta1-x1m(isglob)))/! PARALL
     &                          (x1m(isglob)-x1m(imglob))! PARALL
                           
                           prsint1=(pr(im,js,ks)  * (delta3-x3m(km))
     &                          - pr(im,js,km)  * (delta3-x3m(ks)))/
     &                          (x3m(ks)-x3m(km))
                           
                           prsint2=(pr(is,js,ks)  * (delta3-x3m(km))
     &                          - pr(is,js,km)  * (delta3-x3m(ks)))/
     &                          (x3m(ks)-x3m(km))
                           
                           prsint4=(prsint1     * (delta1-x1m(isglob))! PARALL
     &                          - prsint2       * (delta1-x1m(imglob)))/! PARALL
     &                          (x1m(imglob)-x1m(isglob))! PARALL

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
                           momento1 = 
     $                          (forza3_p+forza3_v)*(delta2-x20)
     $                          - (forza2_p+forza2_v)*(delta3-x30)
                           momento2 = -
     $                          (forza3_p+forza3_v)*(delta1-x10)
     $                          + (forza1_p+forza1_v)*(delta3-x30)
                           momento3 = 
     $                          (forza2_p+forza2_v)*(delta1-x10)
     $                          - (forza1_p+forza1_v)*(delta2-x20)

c      write(53+myid,44)itr,delta1,delta2,delta3,q1int,q2int,q3int,
c     $     tauvisc1,tauvisc2,tauvisc3,prsint,prsurf


c         write(*,*)'VEL',itr,vrel1,vrel2,vrel3,tauvisc1,tauvisc2,
c     $    tauvisc3,prsurf,versnorm(1),versnorm(2),versnorm(3)
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
                           if (bar2.gt.0.)then
                              vortic = vrel3*sqrt(vrel3**2.+vrel2**2.)/
     $                             abs(vrel3)/delta
                           else
                              vortic = -vrel3*sqrt(vrel3**2.+vrel2**2.)/
     $                             abs(vrel3)/delta
                           endif
c                           if(amod(time,tprint).lt.dt) then
c                              write(99,42)itr,bar1,bar2,bar3,
c     $                             prsurf,q2int,q3int
c                           endif
                           tauviscmod= sqrt(tauvisc1**2.+tauvisc2**2.
     $                          +tauvisc3**2.)
                        endif
                     end do
                  end if
               end do
            end if       
         end do
 42      format(i4,7(2x,e12.5))
 43      format(i4,8(2x,e12.5))
 44      format(i4,12(2x,e12.5))
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
         risulM1 = risulM1 + momento1
         risulM2 = risulM2 + momento2
         risulM3 = risulM3 + momento3
         write(43+myid,44)itr,(forza1_p+forza1_v),(forza2_p+forza2_v),
     $       (forza3_p+forza3_v),(delta1-x10),(delta2-x20),(delta3-x30),
     $        momento1,momento2,momento3,risulM1,risulM2,risulM3
 987     continue            
      end do 
      
      risul3_hat = (risul2_v+risul2_p)*sin(alphax) +
     $     (risul3_v+risul3_p)*cos(alphax)
      risul2_hat = (risul2_v+risul2_p)*cos(alphax) - 
     $     (risul3_v+risul3_p)*sin(alphax)

      call MPI_ALLREDUCE(risul3_v,risul3_v_all,1,MPI_DOUBLE_PRECISION,
     $     MPI_SUM,comm,ierr)
      call MPI_ALLREDUCE(risul2_v,risul2_v_all,1,MPI_DOUBLE_PRECISION,
     $     MPI_SUM,comm,ierr)
      call MPI_ALLREDUCE(risul1_v,risul1_v_all,1,MPI_DOUBLE_PRECISION,
     $     MPI_SUM,comm,ierr)
      call MPI_ALLREDUCE(risul3_p,risul3_p_all,1,MPI_DOUBLE_PRECISION,
     $     MPI_SUM,comm,ierr)
      call MPI_ALLREDUCE(risul2_p,risul2_p_all,1,MPI_DOUBLE_PRECISION,
     $     MPI_SUM,comm,ierr)
      call MPI_ALLREDUCE(risul1_p,risul1_p_all,1,MPI_DOUBLE_PRECISION,
     $     MPI_SUM,comm,ierr)
      call MPI_ALLREDUCE(risulM1,risulM1_all,1,MPI_DOUBLE_PRECISION,
     $     MPI_SUM,comm,ierr)
      call MPI_ALLREDUCE(risulM2,risulM2_all,1,MPI_DOUBLE_PRECISION,
     $     MPI_SUM,comm,ierr)
      call MPI_ALLREDUCE(risulM3,risulM3_all,1,MPI_DOUBLE_PRECISION,
     $     MPI_SUM,comm,ierr)
      call MPI_ALLREDUCE(risul3_hat,risul3_hat_all,1,
     $     MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      call MPI_ALLREDUCE(risul2_hat,risul2_hat_all,1,
     $     MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      risul3_v = risul3_v_all
      risul2_v = risul2_v_all
      risul1_v = risul1_v_all
      risul3_p = risul3_p_all
      risul2_p = risul2_p_all
      risul1_p = risul1_p_all
      risulM1 = risulM1_all
      risulM2 = risulM2_all
      risulM3 = risulM3_all
      risul3_hat = risul3_hat_all
      risul2_hat = risul2_hat_all
      if (myid.eq.0)then
         write(65,1000)time,risul1_v,risul1_p,risul2_v,risul2_p,
     $        risul3_v,risul3_p,(risul1_v+risul1_p),risul2_hat,
     $        risul3_hat,risulM1,risulM2,risulM3
      endif
c      cx_v = 2.0*risul3_v*inv_uzero
c      cy_v = 2.0*risul2_v*inv_uzero
c      cz_v = 2.0*risul1_v*inv_uzero
c      cx_p = 2.0*risul3_p*inv_uzero
c      cy_p = 2.0*risul2_p*inv_uzero
c      cz_p = 2.0*risul1_p*inv_uzero
c      cM_x = 2.0*risulM_x*inv_uzero
c      write(66,1001) time,cx_v,cy_v,cz_v,cx_p,cy_p,cz_p,cM_x
 1000 format(13(2x,e12.5))
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
c     c     endif
c      if(amod(time,tprint).lt.dt) then
c      open(31,file=namfi2,form='formatted',status='unknown')
c      write(31,*) 'Variables = X, Y, Z, Tau(module), TauX, TauY,
c     &     TauZ,Pressure'
c      write(31,*)'ZONE N = ',3*nb,',E = ',nb,',F=FEPOINT,ET=TRIANGLE'
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
c      endif
         
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
c      subroutine readgeo3d(time,xyzb,norb,nb,dim3_min,dim3_max,dim2_min,
c     $     dim2_max,dim1_min,dim1_max,namefile,pitch)     
      subroutine readgeo3d(time,namefile,pitch)

c     subroutine readgeo3d(nb,x1,x2,y1,y2,z1,z2,namefile,pitch)
!--------------------------------------------------------------------
!  scope: Read Input Geometry (STL)
!  input: -
!  output: xyzc Input traingle Vertices (nb)
!--------------------------------------------------------------------
      include 'param.f'

      real time
c      real xcen,ycen,angtime,radpun,angpun,xpun,ypun
      real dim3_min,dim3_max,dim2_min,dim2_max,dim1_min,dim1_max
      real small,pi
c     $      ,pitch
c      real v11,v12,v13,v21,v22,v23,v31,v32,v33,sp,area,d32,d21,d13,
c     $     sup_tot
c      real norb(10000,3)
c      real xyzb(10000,9)
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
            goto 99
         endif
         read(11,*)string,string1,(norb(j,i),i=1,3) 
         goto 10
      endif
      if (string.eq.'outer') then
c          j = j+1
         read(11,*)string,(xyzb(j,i),i=1,3)
c          xpun=xyzb(j,1)-x30
c          ypun=xyzb(j,2)-x20
c          radpun=sqrt(xpun**2.+ypun**2.)
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

         dim3_min = min(dim3_min,xyzb(j,1),xyzb(j,4),xyzb(j,7))
         dim3_max = max(dim3_max,xyzb(j,1),xyzb(j,4),xyzb(j,7))
         dim2_min = min(dim2_min,xyzb(j,2),xyzb(j,5),xyzb(j,8))
         dim2_max = max(dim2_max,xyzb(j,2),xyzb(j,5),xyzb(j,8))
         dim1_min = min(dim1_min,xyzb(j,3),xyzb(j,6),xyzb(j,9))
         dim1_max = max(dim1_max,xyzb(j,3),xyzb(j,6),xyzb(j,9))
         
         goto 10
      endif
      if (string.eq.'endloop') goto 10
 98   write(*,*)' Error in STL file ',string,j
 99   close(11)

      if (myid.eq.0)write(*,*)' Dimensioni corpo - i,j,k '
     $     ,dim1_min,dim1_max,
     $     dim2_min,dim2_max,dim3_min,dim3_max

c      do ic = 1,n1m      !PARALL
      do ic = 1,n1mglob
         if ( x1c(ic).lt.dim1_min )ib_min = ic
         if ( x1c(ic).lt.dim1_max )ib_max = ic
      enddo
      ib_min = ib_min - 5
      ib_max = ib_max + 6
      do jc = 1,n2m
         if ( x2c(jc).lt.dim2_min )jb_min = jc
         if ( x2c(jc).lt.dim2_max )jb_max = jc
      enddo
      jb_min = jb_min - 5
      jb_max = jb_max + 5
      do kc = 1,n3m
         if ( x3c(kc).lt.dim3_min )kb_min = kc
         if ( x3c(kc).lt.dim3_max )kb_max = kc
      enddo
      kb_min = kb_min - 5
      kb_max = kb_max + 5
      if (myid.eq.0)write(*,*)' Posizione corpo - i,j,k '
     $     ,ib_min,ib_max,
     &  jb_min,jb_max,kb_min,kb_max

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

      return
      end
!--------------------------------------------------------------------