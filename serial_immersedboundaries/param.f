c     parameter (m1=9,m2=259,m3=359,ndd=2,ndv=3)
c     parameter (m1=121,m2=109,m3=192,ndd=2,ndv=3)
c     parameter (m1=151,m2=189,m3=192,ndd=2,ndv=3)
c     parameter (m1=181,m2=181,m3=301,ndd=2,ndv=3)
c     parameter (m1=121,m2=91,m3=151,ndd=2,ndv=3)
c     parameter (m1=241,m2=121,m3=201,ndd=2,ndv=3)
c     parameter (m1=321,m2=136,m3=201,ndd=2,ndv=3)
      parameter (m1=25,m2=251,m3=351,ndd=2,ndv=3)
c     parameter (m1=401,m2=161,m3=161,ndd=2,ndv=3)
c     parameter (m1=91,m2=91,m3=150,ndd=2,ndv=3)
c     parameter (m1=9,m2=270,m3=210,ndd=2,ndv=3)
c     parameter (m1=41,m2=109,m3=192,ndd=2,ndv=3)
C     parameter (m1=193,m2=88,m3=163,ndd=2,ndv=3)
c     parameter (m1=1,m2=257,m3=257,ndd=2,ndv=3)
c     parameter (m1=97,m2=100,m3=193,ndd=2,ndv=3)
c     parameter (m1=49,m2=65,m3=191,ndd=2,ndv=3)
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
c     parameter (m1m=1,m2m=m2-1,m3m=m3-1)
      parameter (m3mh=m3m/2+1)
      parameter (m12=2*m1m+1,m3c=2*m3m,m32=2*m3c)
      parameter (m3p=m3c+2)
      parameter (m1p=m1m+2)
      parameter (m1c=m1m)
      parameter (mpun=100000)
      REAL*4   visct(m1,m2,m3),st(m1,m2,m3,6),forclo(m1,m2,m3)
     %      ,rhs(m1*m2*m3),csrz_t(m1,m2,m3)
      real norb
      common /codin/  csrz_t,nusf
      common /averou/ iav
      common /axsym/  iaxsy
      common /cft/    nx3fft
      common /cordvo/ x1c(m1),x1m(m1),x2c(m2),x2m(m2),x3c(m3),x3m(m3)
      common /ctrdph/ amph(m2),acph(m2),apph(m2)
      common /d1/     re,tfin,eps
      common /d13/    alx1,alx3
      common /d2/     nstop,ntst,npstf,ireset,ntime,ikick
      common /dd2/    tprint,tpin,tmax
      common /d4/     ntiia,ntstf
      common /dim/    n1,n1m,n2,n2m,n3,n3m
      common /dimens/ alx3d
      common /div/    qcap(m1,m2,m3)
      common /dphc/   dph(m1,m2,m3)
      common /enerw/  nini,nfin,nstri   
      common /fftcm1/ ifx1(13)
      common /fftcm3/ trigx1(m12)
      common /force/  ifugo,icorio,ibuo,isca
      common /idtopt/ cflmax,dtmax,cfllim,idtv
      common /indat/  h0,h1,delro
      common /indbo/  imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
     %               ,isym(m1)
      common /inqca/  imxq,jmxq,kmxq
      common /inslip/ inslws,inslwn,inslwe,inslww
      common /invoma/ iazm,jazm,kazm,irrm,jrrm,krrm
     %               ,irzm,jrzm,krzm
      common /inwal/  jpc(m2),jup(m2),jmc(m2),jum(m2)
      common /inwalk/ kpc(m3),kup(m3),kmc(m3),kum(m3)
      common /maxdi/  qqmax,qqtot
      common /mesh/   dx1,dx1q,dx2,dx2q,dx3,dx3q
      common /mmdens/ denmax,denmin,densm
      common /njump/  n1p,n2p,n3p
      common /parcoo/ r0
      common /peckle/ pec
      common /posmax/ imaxv(ndv),jmaxv(ndv),kmaxv(ndv)
      common /rando/  rnd(10000)
      common /residue/resid
      common /rhsc/   rhs
      common /richar/ ri
      common /rossby/ ros
      common /sonda/  coson(65,3),mason(65,3),nson
      common /trici/  ami(m1,m2),aci(m1,m2),api(m1,m2),
     %                fei(m2,m1),q(m2,m1),s(m2,m1)
      common /tricj/  amj(m2),acj(m2),apj(m2)
      common /trick/  amk(m3),ack(m3),apk(m3)
      common /tscoe/  ga,ro,al,nsst
      common /tscoev/ gam(3),rom(3),alm(3)
      common /tstep/  dt,beta,ren
      common /velmax/ vmax(ndv),vmaxo(ndv)
      common /vopa/   vmx,sig,yc2mo,yc3mo
      common /vperin/ vper,epsil,rper,nsect,lamb
      common /wrre/   nwrit,nread
      common /waves/  an(m3),ap(m1),ak3(m3),ak1(m1)
C      NUOVE AGGIUNTE
      common /metrr/  g2m(m2),g2c(m2)
      common /nonunif/strr,alx2,rint,rmed,etdp,strb,istr
      common /cor1j/  ap1j(m2),ac1j(m2),am1j(m2)
      common /cor2j/  ap2j(m2),ac2j(m2),am2j(m2)
      common /cor2je/ ap2je(m2),ac2je(m2),am2je(m2)
      common /cor3j/  ap3j(m2),ac3j(m2),am3j(m2)
      common /corscj/ apscj(m2),acscj(m2),amscj(m2)
      common /anusse/ anusslow,anussupp
c
c
      common/forcoll/ forclo
      common/diffrot/ vdiff(m2),omdiff
      common/geom   / r1,r2,z1,z2,hn
c************************************************************************
      common/dqx2ol/dq2x2o(m1,m2),dq3x2o(m1,m2),dq1x2o(m1,m2) 
C    %  ,ddex2o(m1,m2)
      common/bcodq1/dqb1s(m1,m2),dqb1n(m1,m2)
      common/bcodq2/dqb2s(m1,m2),dqb2n(m1,m2)
     %             ,dqb2dn(m1,m3),dqb2up(m1,m3)
      common/bcodq3/dqb3s(m1,m2),dqb3n(m1,m2)
     %             ,dqb3dn(m1,m3),dqb3up(m1,m3)
      common/bconq1/qb1s(m1,m2),qb1n(m1,m2)
      common/bconq2/qb2s(m1,m2),qb2n(m1,m2)
     %             ,qb2dn(m1,m3),qb2up(m1,m3)
      common/bconq3/qb3s(m1,m2),qb3n(m1,m2)
     %             ,qb3dn(m1,m3),qb3up(m1,m3)
      common/velcon/cou,cor
      common/taus/tau1,tau2,tau3,tau4
      common/tinfl/dft,uinf(m2),deinf(m2),ft,fto,dftdt,ddftdt
      common/inflth/uinfth(m1,m2),uinflw
      common/jbou/jv3,infig,npun,
     % npuns,indgeo(3,mpun,3),npunr,
     % npunz,npunt,npunp,indgeoe(3,mpun,3),
     % npunpl,npunpr
      common/abou/radinf,all,alt,
     % distb(3,mpun)
      common/bbou/ q2bo(mpun),q3bo(mpun),q1bo(mpun)
c************************************************************************
      common/appo/fej(m3,m2),fek(m2,m3),fj(m3,m2),fk(m2,m3)
      common/anim/tframe,imovie
      common/phcok/amphk(m3),acphk(m3),apphk(m3),acphkk(m1,m3)
      common/quafi/yfis(m3,m2),w(m2*m3)
      common/cofisn/mw,np,mp
      common/phcoj/amphj(m2),
     %   acphj(m2),apphj(m2),qsbph(m1,m2),fphj(m1,m2)
      common /cor3ck/  ap3k(m3),ac3k(m3),am3k(m3)
      common /cor2ck/  ap2k(m3),ac2k(m3),am2k(m3)
      common /cor1ck/  ap1k(m3),ac1k(m3),am1k(m3)
      common /cor3ss/  ap3ssk(m3),ac3ssk(m3),am3ssk(m3)
      common /ametrz/  g3c(m3),g3m(m3),str3
      common/nonunz/   rmed31,etdp3,strb3,istr3
      common/piston/   velval1,velval2,velvalp
      common/geocon/   dims_0,dimss,pisth,dimr,dimc,dimh,dimd,
     %                 dimg,dimra,dimrb,dimrf,dimrv,dimrw,dimrs
      common/visc/   visc(m3)
      common/les/      ell1(m2,m3),ell2(m2,m3),iles
      common/less/     visct,st
      common/addc/     udx2m(m2),udx2c(m2),udx3m(m3),udx3c(m3)
      common /aveonl/  istat,imed,jmed,amasv(m2,m3),amash(m1,m3),igext
c*************************************************************************
      common /rotat/ wx,dwx,alphax,alphay,wy,dwy,wz,dwz,
     %    hy,uhy,duhy,u0x,du0x,u0y,du0y,u0z,du0z,x10,x20,x30
      common /oscilinput/ uzero,uzero_max,amplit,thetamax
      common /readstl/ norb(10000,3),xyzb(10000,9),nb
      common /dimbody/ ib_min,ib_max,jb_min,jb_max,kb_min,kb_max
