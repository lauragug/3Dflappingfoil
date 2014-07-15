c***********************************************
c     Smagorinsky subgrid-scale stress model 
c     constant function of r and z
c***********************************************
      subroutine smagorz(time,q1,q2,q3,dq)
      include 'param.f'
      REAL*4 dq(m1,m2,m3)
      REAL*4 q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)
c
c      the rate of strain tensor S_ij at the grid filter
c      level is computed
c
      call strain(q1,q2,q3)
c                    _         _   _
c       dq contains |S|=sqrt(2*S_ij*S_ij) 
c                  _
c       dq stores |S| for the computation of viscosity at the end
c     of the subroutine
c
      do k=1,n3m
        do j=1,n2m
          do i=1,n1m
            dq(i,j,k)= sqrt(
     %         2.*(st(i,j,k,1)**2+st(i,j,k,2)**2+st(i,j,k,3)**2)+ 
     %         4.*(st(i,j,k,4)**2+st(i,j,k,5)**2+st(i,j,k,6)**2) 
     %                             )
          end do
        end do
      end do
       visma=-10000.
       vismi=+10000.
c
c      Smagorinsky coefficient: this might change from flow to flow
c
       ccsrz = 0.01
       do k=1,n3m
         do j=1,n2m
           do i=1,n1m
             visct(i,j,k)=ccsrz*ell1(j,k)*dq(i,j,k)*ren
             visma=max(visct(i,j,k),visma)
             vismi=min(visct(i,j,k),vismi)
           end do
         end do
       end do
       write(29,1000) time,visma/ren,vismi/ren
       write(6,*) 
     %         ' Nu_M = ',visma/ren,' Nu_m = ',vismi/ren
 1000  format(4(2x,e10.4))
       return
       end 
c***********************************************
c     dynamic subgrid-scale stress model 
c     constant function of r and z
c***********************************************
      subroutine dynamrz(time,q1,q2,q3,dq)
      include 'param.f'
      REAL*4 dq(m1,m2,m3)
      REAL*4 q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)
      REAL*4 amij(m1,m2,m3,6)
      dimension alijmrz(m2,m3),amijmrz(m2,m3)
      common/vvtt/ visma,vismi
      common/tuvis/ accrz(m2,m3),csrzo(m2,m3)
c
c      the rate of strain tensor S_ij at the grid filter
c      level is computed
c
      call strain(q1,q2,q3)
c                    _         _   _
c     dq   contains |S|=sqrt(2*S_ij*S_ij) 
c                  _
c       dq stores |S| for the computation of viscosity at the end
c     of the subroutine
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,dq,st)
!$SGI$ NEST(k,j,i)

      do k=1,n3m
        do j=1,n2m
          do i=1,n1m
            dq(i,j,k)= sqrt(
     %         2.*(st(i,j,k,1)**2+st(i,j,k,2)**2+st(i,j,k,3)**2)+ 
     %         4.*(st(i,j,k,4)**2+st(i,j,k,5)**2+st(i,j,k,6)**2) )
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c                    _ _
c  calculates   dph=|S|S_ij and filter this quantity
c                                  ___
c                                 /_ _\
c     m_ij contains the filtered (|S|S_ij)
c
      do n=1,6
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,dq,st,dph,forclo,amij)
        do k=1,n3m
          do j=1,n2m
            do i=1,n1m
              dph(i,j,k)=st(i,j,k,n)*dq(i,j,k)
c                         /_\
c   compute   the filtered S_ij
              forclo(i,j,k)=st(i,j,k,n)
            end do
          end do
        end do
!$OMP  END PARALLEL DO
        call filter(dph,amij(1,1,1,n))
        call filters(forclo,st(1,1,1,n))
      end do
c                                /_\
c     st_ij contains the filtered S_ij
c
c                           /_\       /_\  /_\
c     dph contains filtered |S|=sqrt(2*S_ij*S_ij) 
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,dph,st)
!$SGI$ NEST(k,j,i)
      do k=1,n3m
        do j=1,n2m
          do i=1,n1m
            dph(i,j,k)= sqrt(
     %           2.*(st(i,j,k,1)**2+st(i,j,k,2)**2+st(i,j,k,3)**2)+
     %           4.*(st(i,j,k,4)**2+st(i,j,k,5)**2+st(i,j,k,6)**2) 
     %                          )
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c                  __               ____
c                 /__\/_\/_\   __  /_  _\
c    compute M_ij= l2*|S| S_ij-l1*(|S| S_ij)
c
              rell= 6.
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,rell,dph,st,amij)
!$SGI$ NEST(k,j,i)
        do k=1,n3m
          do j=1,n2m
            do i=1,n1m
c     
c             according to Lund's derivation when the 
c             test filter is twice the grid filter
c             the square of the ratio of the grid filter size
c             is 6 and not 4 (according to the coefficients of
c             the trapezoidal rule)
c
              amij(i,j,k,1)=rell*dph(i,j,k)*st(i,j,k,1)-amij(i,j,k,1)
              amij(i,j,k,2)=rell*dph(i,j,k)*st(i,j,k,2)-amij(i,j,k,2)
              amij(i,j,k,3)=rell*dph(i,j,k)*st(i,j,k,3)-amij(i,j,k,3)
              amij(i,j,k,4)=rell*dph(i,j,k)*st(i,j,k,4)-amij(i,j,k,4)
              amij(i,j,k,5)=rell*dph(i,j,k)*st(i,j,k,5)-amij(i,j,k,5)
              amij(i,j,k,6)=rell*dph(i,j,k)*st(i,j,k,6)-amij(i,j,k,6)
            end do
          end do
        end do
!$OMP  END PARALLEL DO
c
c   contraction of  M_ij M_ij           
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,dph,amij)
      do k=1,n3m
        do j=1,n2m
          do i=1,n1m
       dph(i,j,k)=(amij(i,j,k,1)**2+amij(i,j,k,2)**2+amij(i,j,k,3)**2)+
     %         2.*(amij(i,j,k,4)**2+amij(i,j,k,5)**2+amij(i,j,k,6)**2)
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
c     averaging M_ij M_ij in the azimuthal direction and among the
c     closest neighbours
c     M_ij M_ij is temporarily stored into rhs
c
      call filter(dph,forclo)
c                                                  _____
c                                       /_\ /_\   /_   _\
c     computation of the tensor L_ij = - q_i q_j + q_i q_j
c                                               _
c    applying the test filter to the velocities q_i
c      _
c    ( q_i must be defined at the cell centre)
c                   _
c     st(.,4)<------q_t
c                   _
c     st(.,5)<------q_r
c                   _
c     st(.,6)<------q_x
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,kpv,jpv,ipv,q1,q2,q3,st)
!$OMP$ PRIVATE(kp,jp,ip)
      do k=1,n3m
        kp=kpv(k)
        do j=1,n2m
          jp = jpv(j)
          do i=1,n1m
            ip=ipv(i)
            st(i,j,k,4)=.5*(q1(ip,j,k)+q1(i,j,k))
            st(i,j,k,5)=.5*(q2(i,jp,k)+q2(i,j,k))
            st(i,j,k,6)=.5*(q3(i,j,kp)+q3(i,j,k))
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c                             _    /_\
c     st(.,1)<------ filtered q_t = q_t
c                             _    /_\
c     st(.,2)<------ filtered q_r = q_r
c                             _    /_\
c     st(.,3)<------ filtered q_x = q_x
c
      call filters(st(1,1,1,4),st(1,1,1,1))        
      call filters(st(1,1,1,5),st(1,1,1,2))        
      call filters(st(1,1,1,6),st(1,1,1,3))        
c
c  computation of L_ij  and contraction with M_ij
c
c  qcap=L_ij M_ij
c
c  _tt component
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,dph,st)
      do k=1,n3m
        do j=1,n2m
          do i=1,n1m
            dph(i,j,k)=st(i,j,k,4)*st(i,j,k,4)
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
c     filtered quantities are temporarily stored into visct(
c
      call filter(dph,visct)        
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,qcap,visct,amij,dph,st)
      do k=1,n3m
        do j=1,n2m
          do i=1,n1m
            qcap(i,j,k)=(visct(i,j,k)-st(i,j,k,1)*st(i,j,k,1))*
     %                   amij(i,j,k,1)
            dph(i,j,k)=st(i,j,k,5)*st(i,j,k,5)
          end do
        end do
      end do
!$OMP END  PARALLEL DO
c
      call filter(dph,visct)        
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,qcap,visct,amij,dph,st)
      do k=1,n3m
        do j=1,n2m
          do i=1,n1m
            qcap(i,j,k)=(visct(i,j,k)-st(i,j,k,2)*st(i,j,k,2))*
     %                   amij(i,j,k,2) + qcap(i,j,k)        
            dph(i,j,k)=st(i,j,k,6)*st(i,j,k,6)
          end do
        end do
      end do
!$OMP END  PARALLEL DO
c
      call filter(dph,visct)        
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,qcap,visct,amij,dph,st)
      do k=1,n3m
        do j=1,n2m
          do i=1,n1m
            qcap(i,j,k)=(visct(i,j,k)-st(i,j,k,3)*st(i,j,k,3))*
     %                   amij(i,j,k,3) + qcap(i,j,k)        
            dph(i,j,k)=st(i,j,k,4)*st(i,j,k,5)
          end do
        end do
      end do
!$OMP END  PARALLEL DO
c
      call filter(dph,visct)        
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,qcap,visct,amij,dph,st)
      do k=1,n3m
        do j=1,n2m
          do i=1,n1m
            qcap(i,j,k)=(visct(i,j,k)-st(i,j,k,1)*st(i,j,k,2))*
     %                   amij(i,j,k,6)*2. + qcap(i,j,k)        
            dph(i,j,k)=st(i,j,k,4)*st(i,j,k,6)
          end do
        end do
      end do
!$OMP END  PARALLEL DO
c
      call filter(dph,visct)        
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,qcap,visct,amij,dph,st)
      do k=1,n3m
        do j=1,n2m
          do i=1,n1m
            qcap(i,j,k)=(visct(i,j,k)-st(i,j,k,1)*st(i,j,k,3))*
     %                   amij(i,j,k,5)*2. + qcap(i,j,k)        
            dph(i,j,k)=st(i,j,k,5)*st(i,j,k,6)
          end do
        end do
      end do
!$OMP END  PARALLEL DO
      call filter(dph,visct)        
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,qcap,visct,amij,dph,st)
      do k=1,n3m
        do j=1,n2m
          do i=1,n1m
            qcap(i,j,k)=(visct(i,j,k)-st(i,j,k,2)*st(i,j,k,3))*
     %                   amij(i,j,k,4)*2. + qcap(i,j,k)        
          end do
        end do
      end do
!$OMP END  PARALLEL DO
c
c     qcap  contains L_ij M_ij
c     averaging among the closest neighbours
c     L_ij M_ij is temporarily stored into dph
c
c     call filters(qcap,st(1,1,1,1))        
      av=1./64.
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,kpv,kmv,jpv,jmv,ipv,imv,av,uf,uv)
!$OMP$ PRIVATE(kp,km,jp,jm,ip,im)
      do kc=1,n3m
        kp=kpv(kc)
        km=kmv(kc)
        do jc=1,n2m
          jp=jpv(jc)
          jm=jmv(jc)
          do ic=1,n1m
            ip=ipv(ic)
            im=imv(ic)
      dph(ic,jc,kc)= qcap(ic,jc,kc)
c     dph(ic,jc,kc)=av*(
c    %         (  (qcap(im,jm,kp)+2.*qcap(ic,jm,kp)+qcap(ip,jm,kp)+
c    %         2.*(qcap(im,jm,kc)+2.*qcap(ic,jm,kc)+qcap(ip,jm,kc))+
c    %            qcap(im,jm,km)+2.*qcap(ic,jm,km)+qcap(ip,jm,km)) )
c    %     +2.*(  (qcap(im,jc,kp)+2.*qcap(ic,jc,kp)+qcap(ip,jc,kp)+
c    %         2.*(qcap(im,jc,kc)+2.*qcap(ic,jc,kc)+qcap(ip,jc,kc))+
c    %            qcap(im,jc,km)+2.*qcap(ic,jc,km)+qcap(ip,jc,km)) )
c    %       + (  (qcap(im,jp,kp)+2.*qcap(ic,jp,kp)+qcap(ip,jp,kp)+
c    %       2.*(qcap(im,jp,kc)+2.*qcap(ic,jp,kc)+qcap(ip,jp,kc))+
c    %          qcap(im,jp,km)+2.*qcap(ic,jp,km)+qcap(ip,jp,km)) ) )
          end do
        end do
      end do
!$OMP  END PARALLEL DO

c
c       computation of Smagorinsky-like constant
c
c      initialization of the constant at the beginning
c      (this must be modified if the computation starts
c       from a restart file)
c
       if((ntime.eq.ntiia).and.(ireset.eq.1)) then
         do k=1,n3m
           do j=1,n2m
             do i=1,n1m
               csrz_t(i,j,k)= 0.
             end do
           end do
         end do
         nusf = 0
       end if 
c
       nusf = nusf + 1
       usnusf = 1./float(nusf)
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,forclo,csrz_t,dph,usnusf,visct)
       do k=1,n3m
         do j=1,n2m
           do i=1,n1m
             if(forclo(i,j,k).gt.1.e-8) then
               csrz_t(i,j,k)= csrz_t(i,j,k) 
     %                      -.5*dph(i,j,k)/forclo(i,j,k)
             else 
               csrz_t(i,j,k)= csrz_t(i,j,k)
             end if
c      The constant is averaged in time and temporarily stored in visct
             visct(i,j,k)= csrz_t(i,j,k)*usnusf
           end do
         end do
       end do
!$OMP  END PARALLEL DO
c
c
c      print *,"maxval(qcap)=",maxval(qcap(:,:,:))
c      print *,"minval(qcap)=",minval(qcap(:,:,:))
c      print *,"maxval(dph)=",maxval(st(:,:,:,1))
c      print *,"minval(dph)=",minval(st(:,:,:,1))
c      print *,"maxval(forclo)=",maxval(forclo(:,:,:))
c      print *,"minval(forclo)=",minval(forclo(:,:,:))
c      print *,"maxval(visct)=",maxval(visct(:,:,:))
c      print *,"minval(visct)=",minval(visct(:,:,:))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       visma=-10000.
       vismi=+10000.
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,visct,dq,ren)
!$OMP$ REDUCTION(max:visma)
!$OMP$ REDUCTION(min:vismi)
       do k=1,n3m
         do j=1,n2m
           do i=1,n1m
             visct(i,j,k)=visct(i,j,k)*dq(i,j,k)*ren
c
c     clipping all the points with negative total viscosity
c
             if((visct(i,j,k)+1.).lt.0.) then
               visct(i,j,k) = 0.
             end if
             visma=max(visct(i,j,k),visma)
             vismi=min(visct(i,j,k),vismi)
           end do
         end do
       end do
!$OMP  END PARALLEL DO
c      print *,
c    % ' T = ',time,' Nu_M = ',visma/ren,' Nu_m = ',vismi/ren
c1000  format(4(2x,e10.4))
       return
       end 
c  ****************************** subrout filter  **********************  
      subroutine filter(uv,uf)
      include 'param.f'
c     calculates filtered function using a box filter in physical space. 
c     The filtering is performed in all directions and the grid can be
c     non uniform the the radial and axial directions (2 and 3).
c     Each contribution is weighted according to volume of the corresponding
c     cell.
c     Modified by R.V. 07/13/98
c     The volume weighted filter turns out to be quite heavy to be computed.
c     J. Jimenez, suggests to forget about the cell volume and to compute the
c     filter as if the grid were uniform. Using an LES model an approximation
c     is made anyway.
c     R.V. 07/14/98
      dimension uv(m1,m2,m3)
      REAL*4 uf(m1,m2,m3)
      av=1./64.
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,kpv,kmv,jpv,jmv,ipv,imv,av,uf,uv)
!$OMP$ PRIVATE(kp,km,jp,jm,ip,im)

      do kc=1,n3m
        kp=kpv(kc)
        km=kmv(kc)
        do jc=1,n2m
          jp=jpv(jc)
          jm=jmv(jc)
          do ic=1,n1m
            ip=ipv(ic)
            im=imv(ic)
      uf(ic,jc,kc)=av*(
     %             (  (uv(im,jm,kp)+2.*uv(ic,jm,kp)+uv(ip,jm,kp)+
     %             2.*(uv(im,jm,kc)+2.*uv(ic,jm,kc)+uv(ip,jm,kc))+
     %                 uv(im,jm,km)+2.*uv(ic,jm,km)+uv(ip,jm,km)) )
     %         +2.*(  (uv(im,jc,kp)+2.*uv(ic,jc,kp)+uv(ip,jc,kp)+
     %             2.*(uv(im,jc,kc)+2.*uv(ic,jc,kc)+uv(ip,jc,kc))+
     %                 uv(im,jc,km)+2.*uv(ic,jc,km)+uv(ip,jc,km)) )
     %           + (  (uv(im,jp,kp)+2.*uv(ic,jp,kp)+uv(ip,jp,kp)+
     %             2.*(uv(im,jp,kc)+2.*uv(ic,jp,kc)+uv(ip,jp,kc))+
     %                 uv(im,jp,km)+2.*uv(ic,jp,km)+uv(ip,jp,km)) ) )
          end do
        end do
      end do
!$OMP  END PARALLEL DO
      return
      end
c  ****************************** subrout filters  **********************  
      subroutine filters(uv,uf)
      include 'param.f'
c
c     this routine is identical to "filter", the only difference is that
c     here two REAL arrays are processed. THis is just to save some
c     memory
c     R.V. 10/20/99
      REAL*4 uv(m1,m2,m3)
      REAL*4 uf(m1,m2,m3)
C
      av=1./64.
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,kpv,kmv,jpv,jmv,ipv,imv,av,uf,uv)
!$OMP$ PRIVATE(kp,km,jp,jm,ip,im)
      do kc=1,n3m
        kp=kpv(kc)
        km=kmv(kc)
        do jc=1,n2m
          jp=jpv(jc)
          jm=jmv(jc)
          do ic=1,n1m
            ip=ipv(ic)
            im=imv(ic)
      uf(ic,jc,kc)=av*(
     %             (  (uv(im,jm,kp)+2.*uv(ic,jm,kp)+uv(ip,jm,kp)+
     %             2.*(uv(im,jm,kc)+2.*uv(ic,jm,kc)+uv(ip,jm,kc))+
     %                 uv(im,jm,km)+2.*uv(ic,jm,km)+uv(ip,jm,km)) )
     %         +2.*(  (uv(im,jc,kp)+2.*uv(ic,jc,kp)+uv(ip,jc,kp)+
     %             2.*(uv(im,jc,kc)+2.*uv(ic,jc,kc)+uv(ip,jc,kc))+
     %                 uv(im,jc,km)+2.*uv(ic,jc,km)+uv(ip,jc,km)) )
     %           + (  (uv(im,jp,kp)+2.*uv(ic,jp,kp)+uv(ip,jp,kp)+
     %             2.*(uv(im,jp,kc)+2.*uv(ic,jp,kc)+uv(ip,jp,kc))+
     %                 uv(im,jp,km)+2.*uv(ic,jp,km)+uv(ip,jp,km)) ) )
          end do
        end do
      end do
!$OMP  END PARALLEL DO
      return
      end
c
c
c  ****************************** subrout strain  **********************
c
c     this subroutine computes the rate of strain tensor in Cartesian 
c     coordinates. After computation each element of the tensor is
c     averaged at the cell centre where the turbulent viscosity is
c     eventually computed. The routine is essentially the one used by
c     M. Fatica during his PhD with some minor modifications to allow
c     for non uniform meshes in the radial and axial directions.
c     
c     07/13/98     R. Verzicco
c
c     simplified from cylindrical to Cartesian coordinates
c
c     08/10/99     R. Verzicco
c
      subroutine strain(q1,q2,q3)
      include 'param.f'
      REAL*4 q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)
      real*4 st1(m1,m2),st2(m1,m3),st3(m2,m3)
      dimension dvr(m2),vror(m2)
      dimension strma(6)
c
c     S_tt --->  st(1
c     S_rr --->  st(2
c     S_xx --->  st(3
c     S_rx --->  st(4
c     S_tx --->  st(5
c     S_tr --->  st(6
c
c      DIAGONAL TERMS
c
c
c      st(l,i,j,k)  =S_ll(i,j,k)   at i+1/2,j+1/2,k+1/2
c
c
C     deftma=0.
!$OMP  PARALLEL
!$OMP$ SHARED(n3m,n2m,n1m,st,q1,q2,q3)
!$OMP$ SHARED(kpv,kmv,jmv,ipv,imv,dx1,udx2m,udx3m,udx3c,udx2c)
!$OMP$ PRIVATE(kp,km,jm,jp,ip,im,dq2x3,dq3x2,dq3x1,dq1x3,dq1x2,dq2x1)
!$OMP  DO
      do kc=2,n3m
        kp=kpv(kc)
        km=kmv(kc)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jc+1
          do ic=1,n1m
            ip=ipv(ic)
            im=imv(ic)
c
            st(ic,jc,kc,1)=(q1(ip,jc,kc)-q1(ic,jc,kc))*dx1
            st(ic,jc,kc,2)=(q2(ic,jp,kc)-q2(ic,jc,kc))*udx2m(jc)
            st(ic,jc,kc,3)=(q3(ic,jc,kp)-q3(ic,jc,kc))*udx3m(kc)
c
            dq2x3=(q2(ic,jc,kc)-q2(ic,jc,km))*udx3c(kc)
            dq3x2=(q3(ic,jc,kc)-q3(ic,jm,kc))*udx2c(jc)
            st(ic,jc,kc,4)=(dq2x3+dq3x2)*0.5
c
            dq3x1=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1
            dq1x3=(q1(ic,jc,kc)-q1(ic,jc,km))*udx3c(kc)
            st(ic,jc,kc,5)=(dq3x1+dq1x3)*0.5
c
            dq1x2=(q1(ic,jc,kc)-q1(ic,jm,kc))*udx2c(jc)
            dq2x1=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1
            st(ic,jc,kc,6)=(dq1x2+dq2x1)*0.5
          end do
        end do
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c   Average the off-diagonal terms at the cell centre
c
c   S_rx
c
!$OMP  PARALLEL
!$OMP$ SHARED(n3m,n2m,n1m,st,qcap,rhs,forclo,kpv,ipv)
!$OMP$ PRIVATE(kp,jp,ip,iadd)
!$OMP DO
      do kc=2,n3m
        kp=kpv(kc)
        do jc=1,n2m
          jp=jc+1
          do ic=1,n1m
            iadd=ic+(jc-1)*n1m+(kc-1)*(n1m*n2m)
            ip=ipv(ic)
            qcap(ic,jc,kc)=(st(ic,jc,kc,4)+st(ic,jp,kc,4)+
     %                     st(ic,jc,kp,4)+st(ic,jp,kp,4))*0.25
            forclo(ic,jc,kc)=(st(ic,jc,kc,5)+st(ic,jc,kp,5)+
     %                     st(ip,jc,kc,5)+st(ip,jc,kp,5))*0.25
            rhs(iadd)=(st(ic,jc,kc,6)+st(ic,jp,kc,6)+
     %                     st(ip,jc,kc,6)+st(ip,jp,kc,6))*0.25
          end do
        end do
      end do
!$OMP END DO
!$OMP DO
      do kc=2,n3m
        do jc=1,n2m
          do ic=1,n1m
            iadd=ic+(jc-1)*n1m+(kc-1)*(n1m*n2m)
            st(ic,jc,kc,4)=qcap(ic,jc,kc)
            st(ic,jc,kc,5)=forclo(ic,jc,kc)
            st(ic,jc,kc,6)=rhs(iadd)
          end do
        end do
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c    rate of strain tensor at the boundaries
c
        do jc=1,n2m
          do ic=1,n1m
            st(ic,jc,n3,1) = st(ic,jc,n3m,1)
            st(ic,jc,n3,2) = st(ic,jc,n3m,2)
            st(ic,jc,n3,3) = st(ic,jc,n3m,3)
            st(ic,jc,n3,4) = st(ic,jc,n3m,4)
            st(ic,jc,n3,5) = st(ic,jc,n3m,5)
            st(ic,jc,n3,6) = st(ic,jc,n3m,6)
            st(ic,jc,1,1)  = st(ic,jc,2,1)
            st(ic,jc,1,2)  = st(ic,jc,2,2)
            st(ic,jc,1,3)  = st(ic,jc,2,3)
            st(ic,jc,1,4)  = st(ic,jc,2,4)
            st(ic,jc,1,5)  = st(ic,jc,2,5)
            st(ic,jc,1,6)  = st(ic,jc,2,6)
          end do
        end do
c     
      return
      end
      
c  ****************************** subrout filters  **********************
      subroutine filterd(uv,uf)
      include 'param.f'
c
c     this routine is identical to "filter", the only difference is that
c     here two REAL arrays are processed. THis is just to save some
c     memory
c     R.V. 10/20/99
      REAL*4 uv(m1,m2,m3)
      REAL*8 uf(m1,m2,m3)
C
c     av=1./64.
      av=.015625
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,kpv,kmv,jpv,jmv,ipv,imv,av,uf,uv)
!$OMP$ PRIVATE(kp,km,jp,jm,ip,im)
      do kc=1,n3m
        kp=kpv(kc)
        km=kmv(kc)
        do jc=1,n2m
          jp=jpv(jc)
          jm=jmv(jc)
          do ic=1,n1m
            ip=ipv(ic)
            im=imv(ic)
      uf(ic,jc,kc)=av*(
     %             (  (uv(im,jm,kp)+2.*uv(ic,jm,kp)+uv(ip,jm,kp)+
     %             2.*(uv(im,jm,kc)+2.*uv(ic,jm,kc)+uv(ip,jm,kc))+
     %                 uv(im,jm,km)+2.*uv(ic,jm,km)+uv(ip,jm,km)) )
     %         +2.*(  (uv(im,jc,kp)+2.*uv(ic,jc,kp)+uv(ip,jc,kp)+
     %             2.*(uv(im,jc,kc)+2.*uv(ic,jc,kc)+uv(ip,jc,kc))+
     %                 uv(im,jc,km)+2.*uv(ic,jc,km)+uv(ip,jc,km)) )
     %           + (  (uv(im,jp,kp)+2.*uv(ic,jp,kp)+uv(ip,jp,kp)+
     %             2.*(uv(im,jp,kc)+2.*uv(ic,jp,kc)+uv(ip,jp,kc))+
     %                 uv(im,jp,km)+2.*uv(ic,jp,km)+uv(ip,jp,km)) ) )
          end do
        end do
      end do
!$OMP  END PARALLEL DO
      return
      end
