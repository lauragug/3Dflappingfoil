c
c***********************************************************************
c
c***********************************************************************
c                                                                      *
c   ************** TOPOGR ************************                     *
c                                                                      *
c***********************************************************************
      subroutine topogr
c
c     This routnine finds the indices in the computational grid
c     close to the physical location of the sloping bottom.
c
c     implicit real (a-h,o-z)
      include 'param.f'
      real xx(m1,m2),yy(m1,m2),zpl(m3),aaa
       parameter (one=1.)
c 
c ciclo di lettura delle metriche
c     
       
     
      if(infig.ne.-1) then
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pi = 2.*asin(1.)
c       all1 = 1.
c       all2 = 0.909090909
c       r_in = 1.2
C       r_ex = rext   !Questo va' bene per la marmitta

c Il raggio di uscita dell' uscita dell' ugello e' posto
c uguale a 1 in accordo con la definizione del numero di
c Reynolds in questa simulazione.
        r_ex = 1.  !Raggio uscita parete interna
        r_in_1 = 5.0* r_ex    !Raggio entrata parete esterna
        ltot = all1 + all2 + all3 !Lunghezza totale ugello
       do l = 1,3 !{ start do over the 3 velocity components
      n=0
c     l = 1   Q_1 vel. component
c     l = 2   Q_2 vel. component
c     l = 3   Q_3 vel. component
c
      do k=1,n3m
        km=kmv(k)
        kp=kpv(k)
        ze=zm(k)
        zem=zm(km)
        zep=zm(kp)
        if(l.eq.3) then
          ze=zz(k)
          zem=zz(km)
          zep=zz(kp)
        end if
        do j=1,n2m
          jm=jmv(j)
          jp=jpv(j)
          ye=rm(j)
          yem=rm(jm)
          yep=rm(jp)
          if(l.eq.2) then
            ye=rc(j)
            yem=rc(jm)
            yep=rc(jp)
          end if
          do i=1,n1m
            im = imv(i)
            ip = ipv(i)
            xe=thetam(i)
            xem=thetam(im)
            xep=thetam(ip)
            if(l.eq.1) then
              xe=thetac(i)
              xem=thetac(im)
              xep=thetac(ip)
            end if
c
c         FUNZIONI CUBICHE PER L'UGELLO
c
c FF all1 e all2 corrispondono a L1 e L2 della marmitta
          if((ze.ge.all1).and.(ze.le.(all1+all2))) then
           rnoz_in = r_in + 3.*(r_ex-r_in)*((ze-all1)/(all2))**2.
     %              +2.*(r_in-r_ex)*((ze-all1)/(all2))**3.
            rnoz_ex = r_in_1 + 3.*(1.1*r_ex-r_in_1)*
     %       ((ze-all1)/(all2))**2.+ 2.*(r_in_1-1.1*r_ex)
     %        *((ze-all1)/(all2))**3.
          end if
c
c         TRATTO UNIFORME
c
          if(ze.lt.all1) then
           rnoz_in = r_in 
           rnoz_ex = r_in_1
          end if
          
          if((ze.gt.(all1+all2)).and.(ze.le.ltot)) then
           rnoz_in = r_ex
           rnoz_ex = 1.1*r_ex
          end if
c
c         USCITA UGELLO
c
          if((ze.gt.ltot).and.(zem.le.ltot).and.
     %       (ye.ge.rnoz_in).and.(ye.le.rnoz_ex)) then
                n=n+1
                indgeo(l,n,1)=i
                indgeo(l,n,2)=j
                indgeo(l,n,3)=k
                indgeoe(l,n,1)=i 
                indgeoe(l,n,2)=j
                indgeoe(l,n,3)=kp
                delta1x=(zep-ze)
                delta2x=(ze-ltot)
                distb(l,n)= delta2x/(delta1x+delta2x)
              end if
c         
c        PARETE LATERALE UGELLO
c
           if(ze.le.(ltot)) then
              if((ye.lt.rnoz_in).and.(yep.ge.rnoz_in)) then
                n=n+1
                indgeo(l,n,1)=i
                indgeo(l,n,2)=j
                indgeo(l,n,3)=k
                indgeoe(l,n,1)=i 
                indgeoe(l,n,2)=jm
                indgeoe(l,n,3)=k
                delta1x=(ye-yem)
                delta2x=(rnoz_in-ye)
                distb(l,n)= delta2x/(delta1x+delta2x)
              end if
              if((yem.le.rnoz_ex).and.(ye.gt.rnoz_ex)) then
                n=n+1
                indgeo(l,n,1)=i
                indgeo(l,n,2)=j
                indgeo(l,n,3)=k
                indgeoe(l,n,1)=i 
                indgeoe(l,n,2)=jp 
                indgeoe(l,n,3)=k
                delta1x=(yep-ye)
                delta2x=-(rnoz_ex-ye)
                distb(l,n)= delta2x/(delta1x+delta2x)
              end if
c
c            PUNTI INTERNI
c
              if((ye.ge.rnoz_in).and.(ye.le.rnoz_ex)) then
                n=n+1
                indgeo(l,n,1)=i
                indgeo(l,n,2)=j
                indgeo(l,n,3)=k
                indgeoe(l,n,1)=i 
                indgeoe(l,n,2)=j
                indgeoe(l,n,3)=k
                distb(l,n)= 0.
              end if
            end if
          end do
        end do
      end do
      if(l.eq.1) then
        if(n.gt.mpun)
     %  write(*,*) 'Dim max di indgeot e'' stata superata n=',n
        npunt= n
        write(6,*) ' for Q_1 N = ',npunt
      end if
      if(l.eq.2) then
        if(n.gt.mpun)
     %  write(*,*) 'Dim max di indgeor e'' stata superata n=',n
        npunr= n
        write(6,*) ' for Q_2 N = ',npunr
      end if
      if(l.eq.3) then
        if(n.gt.mpun)
     %  write(*,*) 'Dim max di indgeoz e'' stata superata n=',n
        npunz= n
        write(6,*) ' for Q_3 N = ',npunz
      end if
      end do   !} end do over the 3 velocity components return
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         do   k=1,n3
         do j=1,n2
          do i=1,n1
            forclo(i,j,k)=0.
           end do
         end do
        end do
 
        do n=1,npunz
          i=indgeo(3,n,1)
          j=indgeo(3,n,2)
          k=indgeo(3,n,3)
           forclo(i,j,k)=1.-distb(3,n)
        end do
 
        do n=1,npunr
          i=indgeo(2,n,1)
          j=indgeo(2,n,2)
          k=indgeo(2,n,3)
           qcap(i,j,k)=1.-distb(2,n)
        end do
 
        do n=1,npunt
          i=indgeo(1,n,1)
          j=indgeo(1,n,2)
          k=indgeo(1,n,3)
           qcap(i,j,k)=1.-distb(1,n)
        end do

        do k = 1,n3m
          do j = 1,n2m
            do i = 1,n1m
              forclo(i,j,k) = 
     %      ( forclo(i,j,k)+dph(i,j,k)+qcap(i,j,k) )*0.333333
            end do
          end do
        end do
 
      n1pp=(n1-1)/n1p+1
      n2pp=(n2-1)/n2p+1
      n3pp=(n3-1)/n3p+1
c
      do 1 i=1,n1
        do 1 j=1,n2
          xx(i,j)=rc(j)*cos(thetac(i))
          yy(i,j)=rc(j)*sin(thetac(i))
    1 continue
      open(93,file='coplo.dat',form='unformatted')
      rewind 93
      write(93) n1pp,n2pp,n3pp
      write(93)
     %          (((xx(i,j),i=1,n1,n1p),j=1,n2,n2p),k=1,n3,n3p),
     %          (((yy(i,j),i=1,n1,n1p),j=1,n2,n2p),k=1,n3,n3p),
     %          (((zz(k),i=1,n1,n1p),j=1,n2,n2p),k=1,n3,n3p)
       close(93)
      open(93,file='forclo.dat',form='unformatted')
      rewind 93
      aaa = 1.
      write(93) n1pp,n2pp,n3pp,1
c      write(93) aaa,lamb,aaa,aaa
       write(93) 
     %     (((forclo(i,j,k),i=1,n1,n1p),j=1,n2,n2p),k=1,n3,n3p)
       close(93)

C prova formatted
c     open(93,file='forclo.tst',form='formatted')
c      rewind 93
c        do j=1,n2
c         do k=1,n3
c          if (forclo(1,j,k).ne.0) then
c       write(93,*) forclo(1,j,k),j,k
c          end if
c         end do
c          end do 
c       close(93) 

      open(93,file='indgeoz.dat',form='formatted')
        rewind 93
        do n=1,npunz
      write(93,1000) indgeo(3,n,2),indgeo(3,n,3),distb(3,n)
     % ,rm(indgeo(3,n,2)),zz(indgeo(3,n,3))
         end do
         close(93) 
       
      open(94,file='indgeor.dat',form='formatted')
        rewind 94
         do n=1,npunr
      write(94,1000) indgeo(2,n,2),indgeo(2,n,3),distb(2,n) 
     % ,rc(indgeo(2,n,2)),zm(indgeo(2,n,3))
         end do
        close(94)
C   Fine prova formatted
 1000 format(2(1x,i4),3(1x,f14.7))
       
      end if
      return
      end
c
