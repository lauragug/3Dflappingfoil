      parameter (m1=385,m2=224,m3=274,mpun=77000,mpuni=600000)
      parameter (m1_loc=33)
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1,m3md=m3m+2,m3mp=m3m+4) 
      implicit real*4(a-h,o-z)
      
      common/dim/n1,n1m,n2,n2m,n3,n3m 
      common/njump/n1p,n2p,n3p  
      common/parcoo/r0,alx1
      common/d13/alx3
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q 
      common/nonunif/strr,rext,rint,rmed,etdp,strb,istr
      common/cordvo/x2c(m2),x2m(m2),x3c(m3),x3m(m3)
      common/corrt/x1c(0:m1),x1m(0:m1)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/inwal/jpc(m2),jup(m2),jmc(m2),jum(m2) 
      common/metrr/g2m(m2),g2c(m2)
      common/ametrz/g3c(m3),g3m(m3),str3
      common/nonunz/rmed31,istr3,etdp3,strb3
      common/aveonl/igext  

      common/vperin/vper,epsil,lamb   
      common/tstep/dt,beta,ren   
      common/d1/re,tfin,eps  
      common/inslnsw/inslws,inslwn,inslww,inslwe
      common/section/nsect 
                    
      common/stream/psirv(m2,m3)
      common/circol/gamz,gamr,gamt
      common/impuls/aipz,aipr,aipt
      common/mimpul/ammz,ammr,ammt
      common/parsta/ npunr_st, npunz_st
      common/parleb/ p_es,p_en,d_ls,p_pl,p_cs

      common/refsys/ refsys
      common/replic/irepli
      common/nvort/ nvort

      common/inter/ indgeoi(5,mpuni,3),
c     &  indgeo(5,mpuni,3),   
     %              indgeo(5,mpun,3),    !prova SIGSEGV
     %              indgeoe(5,mpun,3),distb(5,mpun),
     %              npunz,npunz_e,npunz_i,
     %              npunr,npunr_e,npunr_i,
     %              npunt,npunt_e,npunt_i,
     %              npunp,npunp_e,npunp_i,
     %              npunv,npunv_e,npunv_i

      common /oscilinput/ uzero,amplit,thetamax      
      common/tau/tau1,tau2
      common/emotion/ x3b_0,x2b_0,x1b_0,x3b_0_in,x2b_0_in,
     %                angle,omegak,osci,phi,ampp,sigma
      common/rotat/ wx,dwx,alphax,alphay,wy,dwy,wz,dwz
     %,u0x,du0x,u0y,du0y,u0z,du0z,x10,x20,x30,hx,hy,uhy,duhy,ros




