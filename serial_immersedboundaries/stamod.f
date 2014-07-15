c***********************************************************************
c    
      subroutine stacle
      include 'param.f'
      common/rmsvalv/ q1rmsv(m2,m3), q2rmsv(m2,m3), q3rmsv(m2,m3),
     %                      prrmsv(m2,m3)
      common/avevalv/ q1avev(m2,m3), q2avev(m2,m3), q3avev(m2,m3),
     %                      pravev(m2,m3),nfield
      common/rmsvalh/ q1rmsh(m1,m3), q2rmsh(m1,m3), q3rmsh(m1,m3),
     %                      prrmsh(m1,m3)
      common/avevalh/ q1aveh(m1,m3), q2aveh(m1,m3), q3aveh(m1,m3),
     %                      praveh(m1,m3)
      nfield = 0
      do k=1,n3
        do j=1,n2
          q1rmsv(j,k) = 0.
          q2rmsv(j,k) = 0.
          q3rmsv(j,k) = 0.
          prrmsv(j,k) = 0.
          q1avev(j,k) = 0.
          q2avev(j,k) = 0.
          q3avev(j,k) = 0.
          pravev(j,k) = 0.
        end do
      end do
      do k=1,n3
        do i=1,n1
          q1rmsh(i,k) = 0.
          q2rmsh(i,k) = 0.
          q3rmsh(i,k) = 0.
          prrmsh(i,k) = 0.
          q1aveh(i,k) = 0.
          q2aveh(i,k) = 0.
          q3aveh(i,k) = 0.
          praveh(i,k) = 0.
        end do
      end do
      return
      end
c***********************************************************************
      subroutine stacal(q1,q2,q3,pr,dq)
      include 'param.f'
      real*4   dq(m1,m2,m3)
      real*4   q1(m1,m2,m3), q2(m1,m2,m3), q3(m1,m2,m3) 
      real                  pr(m1,m2,m3)
      dimension q1mev(m2,m3), q2mev(m2,m3), q3mev(m2,m3),
     %                      prmev(m2,m3)
      dimension q1meh(m1,m3), q2meh(m1,m3), q3meh(m1,m3),
     %                      prmeh(m1,m3)
      common/rmsvalv/ q1rmsv(m2,m3), q2rmsv(m2,m3), q3rmsv(m2,m3),
     %                      prrmsv(m2,m3)
      common/avevalv/ q1avev(m2,m3), q2avev(m2,m3), q3avev(m2,m3),
     %                      pravev(m2,m3),nfield
      common/rmsvalh/ q1rmsh(m1,m3), q2rmsh(m1,m3), q3rmsh(m1,m3),
     %                      prrmsh(m1,m3)
      common/avevalh/ q1aveh(m1,m3), q2aveh(m1,m3), q3aveh(m1,m3),
     %                      praveh(m1,m3)
      nfield = nfield + 1
c
c     averaged AVE and RMS vertical slices
c
      do k=1,n3m
        do j=1,n2m
c
          q1avev(j,k) = q1avev(j,k) + q1(imed,j,k)*amasv(j,k)
          q2avev(j,k) = q2avev(j,k) + q2(imed,j,k)*amasv(j,k)
          q3avev(j,k) = q3avev(j,k) + q3(imed,j,k)*amasv(j,k)
c         pravev(j,k) = pravev(j,k) + pr(imed,j,k)*amasv(j,k)
c
          q1rmsv(j,k) = q1rmsv(j,k) + 
     %   (q1(imed,j,k)-q1avev(j,k)/float(nfield))**2*amasv(j,k)
          q2rmsv(j,k) = q2rmsv(j,k) + 
     %   (q2(imed,j,k)-q2avev(j,k)/float(nfield))**2*amasv(j,k)
          q3rmsv(j,k) = q3rmsv(j,k) + 
     %   (q3(imed,j,k)-q3avev(j,k)/float(nfield))**2*amasv(j,k)
c         prrmsv(j,k) = prrmsv(j,k) + 
c    %   (pr(imed,j,k)-pravev(j,k)/float(nfield))**2*amasv(j,k)
c
        end do
      end do
c
c     averaged AVE and RMS horizontal slices
c
      do k=1,n3m
        do i=1,n1m
c
          q1aveh(i,k) = q1aveh(i,k) + q1(i,jmed,k)*amash(i,k)
          q2aveh(i,k) = q2aveh(i,k) + q2(i,jmed,k)*amash(i,k)
          q3aveh(i,k) = q3aveh(i,k) + q3(i,jmed,k)*amash(i,k)
c         praveh(i,k) = praveh(i,k) + pr(i,jmed,k)*amash(i,k)
c
          q1rmsh(i,k) = q1rmsh(i,k) + 
     %   (q1(i,jmed,k)-q1aveh(i,k)/float(nfield))**2*amash(i,k)
          q2rmsh(i,k) = q2rmsh(i,k) + 
     %   (q2(i,jmed,k)-q2aveh(i,k)/float(nfield))**2*amash(i,k)
          q3rmsh(i,k) = q3rmsh(i,k) + 
     %   (q3(i,jmed,k)-q3aveh(i,k)/float(nfield))**2*amash(i,k)
c         prrmsh(i,k) = prrmsh(i,k) + 
c    %   (pr(i,jmed,k)-praveh(i,k)/float(nfield))**2*amash(i,k)
c
        end do
      end do
c
      return  
      end
c    
c***********************************************************************
      subroutine stawri(time)
      include 'param.f'
      common/rmsvalv/ q1rmsv(m2,m3), q2rmsv(m2,m3), q3rmsv(m2,m3),
     %                      prrmsv(m2,m3)
      common/avevalv/ q1avev(m2,m3), q2avev(m2,m3), q3avev(m2,m3),
     %                      pravev(m2,m3),nfield
      common/rmsvalh/ q1rmsh(m1,m3), q2rmsh(m1,m3), q3rmsh(m1,m3),
     %                      prrmsh(m1,m3)
      common/avevalh/ q1aveh(m1,m3), q2aveh(m1,m3), q3aveh(m1,m3),
     %                      praveh(m1,m3)
      open(59,file='stafield.dat',form='unformatted')
      rewind(59)
      write(59) nfield,n2,n3               
      write(59) a,a,a,time
      write(59) 
     %           q1avev,
     %           q2avev,
     %           q3avev,
     %           pravev,
     %           q1rmsv,
     %           q2rmsv,
     %           q3rmsv,
     %           prrmsv,
     %           q1aveh,
     %           q2aveh,
     %           q3aveh,
     %           praveh, 
     %           q1rmsh,
     %           q2rmsh,
     %           q3rmsh,
     %           prrmsh 
      close(59)
c 
      return  
      end
c 
c***********************************************************************
      subroutine starea
      include 'param.f'
      common/rmsvalv/ q1rmsv(m2,m3), q2rmsv(m2,m3), q3rmsv(m2,m3),
     %                      prrmsv(m2,m3)
      common/avevalv/ q1avev(m2,m3), q2avev(m2,m3), q3avev(m2,m3),
     %                      pravev(m2,m3),nfield
      common/rmsvalh/ q1rmsh(m1,m3), q2rmsh(m1,m3), q3rmsh(m1,m3),
     %                      prrmsh(m1,m3)
      common/avevalh/ q1aveh(m1,m3), q2aveh(m1,m3), q3aveh(m1,m3),
     %                      praveh(m1,m3)
      real*4 aaa
      open(59,file='stafield.dat',form='unformatted')
      read(59) nfield,n2pp,n3pp               
      read(59) a,a,a,aaa
      read(59) 
     %           q1avev,
     %           q2avev,
     %           q3avev,
     %           pravev,
     %           q1rmsv,
     %           q2rmsv,
     %           q3rmsv,
     %           prrmsv,
     %           q1aveh,
     %           q2aveh,
     %           q3aveh,
     %           praveh, 
     %           q1rmsh,
     %           q2rmsh,
     %           q3rmsh,
     %           prrmsh 
      close(59)
      return  
      end
