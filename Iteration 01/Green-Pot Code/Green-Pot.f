ccADJHELASTIC, GREEN-POT.  ELASTIC - ELASTIC SCATTERING OF ELECTRONS FROM ADJH0000
cc1   IONS AND ATOMS.  D.R. SCHULTZ, C.O. REINHOLD.                       ADJH0000
cc///////////////////////////////////////////////////////////////////////
cc    Calculation of elastic scattering cross sections of electrons 
cc                      at ions and atoms
cc
cc The interaction potential is contained in the routine ''pot'' and
cc can be set to any desired form which rises less rapidly than the 
cc centrifugal potential at the origin. Presently implemented is the
cc parameterization of the potential following Garvey, Jackman, and
cc Green (Phys. Rev. 12, 1144 (1975)).  Use the auxiliary code Green_pot
cc to determine these parameters, or use the tabulation and formulae given
cc by Garvey et al. to compute them.
cc
cc The asymptotic behavior of the wavefunction at large distances is
cc matched to Coulomb functions to include treated scattering from
cc ions.  If the Coulomb functions diverge, increase the input
cc variable ''rend0'' (matching radius).  For zero asymptotic charges
cc (atoms), the Coulomb functions authomatically behave as the 
cc appropriate Bessel functions.
cc                                    
cc The integration of the radial Schrodinger equation is performed based 
cc on Johnson's log derivative method (J. Comput. Phys. 13, 445 (1973)). 
cc                                                   
cc The evaluation of the phase shifts is vectorized.
cc
cc See routine readit at the bottom for a brief explanation of the 
cc input data.
cc
cc D.R. Schultz and C.O. Reinhold, version 5/15/98
cc Physics Division, Oak Ridge National Laboratory
cc///////////////////////////////////////////////////////////////////////
c
c      program elastic
c
c      implicit real*8(a-h,o-z)
c      parameter (nl = 1000)
c      complex*16 cz, cgamma, cg, ci, cratio, cdexpon, cfc, cfnc, cf
c      real*8 logder
c      dimension fl(nl), z11(nl), phase(nl), sigma(nl)
c      dimension f(0:nl), g(0:nl), fp(0:nl), gp(0:nl)
c      common veloc, en
c      common /modpot/ zmod, xi, eta, rmod, zasy, nmod
c      common /data  / rmu,rstart,rend0,spac,lmin,lmax,lspc
c      common /unit  / iout
c      common /born  / iborn
c
cc/// open output  
c      iout=6
c      open(iout,file='elastic.out',status='unknown')
c
cc/// read and echo input
c      call readit
c 
cc/// definitions 
c      pi = 3.14159265358979323846d0
c      ci = (0.d0,1.d0)
c 
c      veloc = sqrt(2.d0*en/rmu)
c      rmu2 = 2.d0*rmu
c      p1 =  sqrt(en*rmu2)
c 
c      etahyp = zasy / p1 + 1.e-30
c 
c      con = 0.125d0
c      lmin1 = lmin + 1
c      lmmx = lmax + 1
c 
c      nos = 2*int ((rend0-rstart)/2.d0/spac)+2
c      space = (rend0-rstart)/ dfloat(nos)
c      del1 = space*space*rmu2/3.d0
c      del2 = 2.d0*del1
c      del4 = 4.d0*del1
c 
c      do  ll = lmin1,lmmx,lspc
c         l = ll-1
c         fl(ll) =  dfloat (l*(l+1))/rmu2
c         z11(ll) = 1.e20
c      end do
c 
cc/// Integrate radial Schrodinger equation for log derivative
c      do ll = lmin1,lmmx,lspc
c         r=rstart
c         do nr = 1,nos,2
c            r=r+space
c            call pot(v,r)
c            v11=(en-v)*del4
c            cen=1.d0/(r**2)*del4
c            r=r+space
c            call pot(v,r)
c            vv11=(en-v)*del2
c            ccen=1.d0/(r**2)*del2
c            y11=z11(ll)
c            w11=v11-cen*fl(ll)
c            den=1.d0+con*w11
c            u11=w11/den
c            den=1.d0+y11
c            y11=y11/den-u11
c            u11=vv11-ccen*fl(ll)
c            den=1.d0+y11
c            z11(ll)=y11/den-u11
c         end do
c      end do 
c 
c      do ll=lmin1,lmmx,lspc
c         u11=vv11-ccen*fl(ll)
c         z11(ll)=(z11(ll)+u11/2.d0)/space
c      end do
c  
cc/// calculation of the phase shifts 
c      rho = rend0 * p1
c      tol = 1.0e-04
c      max = 50
c 
c      call cwf(lmax, rho, etahyp, tol, max, f, g, fp, gp)
c 
c      hcon = 4.d0*pi/p1/p1
c      tcs=0.d0
c      do ll=lmin1,lmmx,lspc
c         l = ll - 1
c         logder = z11(ll) / p1
c         fhyp = f(l)
c         ghyp = g(l)
c         dfhyp = fp(l)
c         dghyp = gp(l)
c         phase(ll) = atan( (logder*fhyp-dfhyp)/(dghyp-ghyp*logder) )
c         cg=(1.d0,0.d0) 		
c         cz = l + 1.d0 + ci*etahyp
c         cg = cgamma(cz)  
c         sigma(ll) = atan2 ( dimag(cg), dreal(cg) )
c         write(iout,6999) l, phase(ll), sigma(ll)
c         tcs = tcs + hcon * (2.d0*l+1) * sin(phase(ll))**2
c      end do
c 
c6999  format(' ','l = ',i3,2x,'phase shift =',1pe15.7,2x,
c     @        'Coulomb phase shift=',1pe15.7)
c
cc/// total cross section for the case of neutral targets (i.e. atoms)
c      if (abs(zasy).lt.1.e-5) then
c      write(iout,*)
c      write(iout,*)'*************************************'
c      write(iout,*)'*** TCS = ', tcs, ' a.u.'
c      write(iout,*)'*************************************'
c      end if
c 
cc/// differential cross sections 
c      degrad = pi/180.
c      angdeg = -0.999
c 
c      write(iout,*)
c      write(iout,7012)
c      if(iborn.eq.1)then
c         write(iout,7011)
c      else
c         write(iout,7010)
c      endif
c7012  format(//,10X,'*** ELASTIC CROSS SECTIONS ***',/)
c7011  format(1x,'theta(deg)    EXACT',11x,'BORN',' a.u.',/,50('='))
c7010  format(1x,'theta(deg)   dcs(a.u.)     Rutherford-q  ',
c     1          '  Rutherford-Z      ',/, 60('='))
c 
c      do i=0,180
c         angdeg = angdeg + 1.d0
c         angrad = angdeg * degrad
c         sint22 = sin(angrad/2.d0)**2
c         ruther = -zasy / ( 2.d0 * p1*p1 * sint22 )
c         ruthzm =  zmod / ( 2.d0 * p1*p1 * sint22 )
c         cratio = cgamma( 1.d0+ci*etahyp ) / cgamma( 1.d0-ci*etahyp )
c         cdexpon = cdexp( -ci*etahyp*dlog(sint22) )
c         cfc = ruther * cratio * cdexpon
c         cost = cos(angrad)
c         cfnc = (0.d0,0.d0)
c         do ll=lmin1,lmmx,lspc
c            l = ll - 1
c            cfnc = cfnc + (2*l+1.d0) * cdexp(2.d0*ci*sigma(ll)) *
c     &                    (cdexp(2.d0*ci*phase(ll))-1.d0) * pl(l,cost)
c         end do
c         cfnc = cfnc / (2.d0*ci*p1)
c         cf = cfc + cfnc
c         dcs = cdabs(cf)**2
c         if(iborn.eq.1)then
c            pmod = 2.d0*p1*sin(angrad/2.d0)
c            dcsb1 = 4.d0*((-zasy+zeff(pmod))/pmod**2)**2
c            write(iout,8001) angdeg, dcs, dcsb1
c         else
c            write(iout,8001) angdeg, dcs, ruther**2, ruthzm**2
c         endif
c      end do
c8001  format(1x,f6.1,6x,1pe10.4,6x,e10.4,6x,e10.4)
c
c      close(iout)
c      stop
c      end
c 
cc///////////////////////////////////////////////////////////////////////
cc Central potential for elastic scattering.  Implemented for the
cc Garvey et al. parameterized Hartree-Fock model potential. 
c
c      subroutine pot(v,r)
c
c      implicit real*8(a-h,o-z)
c      common /modpot/ zmod, xi, eta, rmod, zasy, nmod
c      dimension fact(8)
c      data alp1,alp2,d2/1.38414,2.32,0.11578/
c      data r0,hla/2.373,10.6172/
c      data fact/1.,2.,6.,24.,120.,720.,5040.,40320./
c 
c      o(r) = 1.d0/((eta/xi)*exp(xi*r)-(eta/xi)+1.d0)
c      z(r) = ( dfloat(nmod-1) * (1.d0-o(r)) ) - zmod
c      vmod(r) = z(r) / r
c      if (r.lt.rmod) then
c         vtemp = vmod(r)
c      else
c         vtemp = zasy / r
c      end if
c
c      v = vtemp 
c 
c      return
c      end
c 
cc/////////////////////////////////////////////////////////////////////// 
cc This function subprogram evaluates the complex gamma function,
cc using Stirling's formula (after algorithm in Numerical Recipes, W.H. 
cc Press et al. (Cambridge University Press, Cambridge, 1989)
c
c      function cgamma(cz)
c
c      implicit real*8(a-b,d-h,o-z)
c      implicit complex*16(c)
c      real*8 cof(6)
c      data cof / 76.18009173d0, -86.50532033d0, 24.01409822d0,
c     &           -1.231739516d0, 1.20858003d-3, -5.36382d-6 /
c      data stp / 2.50662827465d0 /
c 
cc/// if 0 < Re(z) < 1 use reflection formula
c 
c      if ( dreal(cz).lt.1. ) then
c         cx = cz
c         ctmp = (cx+5.5d0)**(cx+0.5d0) * cdexp(-(cx+5.5d0))
c         cser = 1.d0
c         do 100 i = 1,6
c            cx = cx + 1.d0
c            cser = cser + cof(i)/cx
c100      continue
c         cgamma = ctmp * stp * cser / cz
c      else
c         cx = cz - 1.d0
c         ctmp = (cx+5.5d0)**(cx+0.5d0) * cdexp(-(cx+5.5d0))
c         cser = 1.d0
c         do 200 i = 1,6
c            cx = cx + 1.d0
c            cser = cser + cof(i)/cx
c200      continue
c         cgamma = ctmp * stp * cser
c      end if
c 
c      return
c      end
c 
cc/////////////////////////////////////////////////////////////////////// 
cc Legendre polynomial of order lin and argument x
c
c      function pl(lin, x)
c
c      implicit real*8(a-h,o-z)
c 
c      p0 = 1.d0
c      p1 = x
c 
c      if(lin.eq.0) then
c         pl = p0
c         return
c      end if
c      if(lin.eq.1) then
c         pl = p1
c         return
c      end if
c 
c      plm2 = p0
c      plm1 = p1
c 
c      do 100 l=2,lin
c 
c         pl = ( (2*l-1)*x*plm1 - (l-1)*plm2 ) / l
c         plm2 = plm1
c         plm1 = pl
c 
c100   continue
c 
c      return
c      end
c 
cc/////////////////////////////////////////////////////////////////////// 
cc this subroutine evaluates the regular and irregular hypergeometric
cc functions in the asymptotic regime, after the method of Abramowitz,
cc as well as their derivatives (see Abramowitz and Stegun, Handbook of
cc Mathematical Functions, (Dover, New York, 1972) chapter 14.5)
c
c      subroutine hyper(l, rho, eta, sigmal, tol, max,
c     &                 fhyp, ghyp, dfhyp, dghyp)
c      implicit real*8(a-h,o-z)
c 
c      pi = 3.14159265358979323846d0
c 
c      thetal = rho - eta*dlog(2*rho) - l*pi*0.5d0 + sigmal
c 
c      f1k = 1.d0
c      g1k = 0.d0
c      f2k = 0.d0
c      g2k = 1.d0 - eta/rho
c 
c      f1 = f1k
c      g1 = g1k
c      f2 = f2k
c      g2 = g2k
c 
c      kount = 0
c 
c      ak = eta * 0.5d0 / rho
c      bk = (l*(l+1) + eta*eta) * 0.5d0 / rho
c 
c      kmin = 1
c      kmax = 10
c 
c1     f1old = f1
c      g1old = g1
c      f2old = f2
c      g2old = g2
c 
c      kount = kount + 1
c 
c      do 100 k=kmin,kmax
c 
c         f1kp1 = ak*f1k - bk*g1k
c         g1kp1 = ak*g1k + bk*f1k
c         f2kp1 = ak*f2k - bk*g2k - f1kp1/rho
c         g2kp1 = ak*g2k + bk*f2k - g1kp1/rho
c 
c         f1 = f1 + f1kp1
c         g1 = g1 + g1kp1
c         f2 = f2 + f2kp1
c         g2 = g2 + g2kp1
c 
c         ak = (2.d0*k+1) * eta / (2.d0*k+2.d0) / rho
c         bk = (l*(l+1) - k*(k+1) + eta*eta) / (2.d0*k+2.d0) / rho
c 
c         f1k = f1kp1
c         g1k = g1kp1
c         f2k = f2kp1
c         g2k = g2kp1
c 
c100   continue
c
cc///diagnostics 
cc     errf1 = abs( (f1 - f1old) / f1 )
cc     errg1 = abs( (g1 - g1old) / g1 )
cc     errf2 = abs( (f2 - f2old) / f2 )
cc     errg2 = abs( (g2 - g2old) / g2 )
cc
cc     errmax = amax1(errf1, errg1, errf2, errg2)
c 
c      wronk = g1*f2-f1*g2
c 
c      errmax = abs ( wronk + 1.d0 )
c 
c      if ( errmax.gt.tol ) then
c         if ( kount.gt.max ) then
c            write(6,*) 'Warning: Hypergeometric convergence failure',
c     &            ' kount=',kount
c            write(6,*) 'Wronskian (-1) = ',wronk
c            write(6,*)
c            stop
c         end if
c         kmin=kmin+10
c         kmax=kmax+10
c         goto 1
c      end if
c 
c      s = sin(thetal)
c      c = cos(thetal)
c 
c      fhyp = g1 * c + f1 * s
c      ghyp = f1 * c - g1 * s
c      dfhyp = g2 * c + f2 * s
c      dghyp = f2 * c - g2 * s
c 
c      return
c      end
c
cc/////////////////////////////////////////////////////////////////////// 
cc This subroutine computes the regular and irregular Coulomb wave  
cc functions and their derivatives for l's from zero to lmax. The first   
cc two are calculated and then the rest found by recursion.  (The method
cc is from Abramowitz and Stegun,Handbook of Mathematical Functions, 
cc (Dover, New York, 1972)).   
c
c      subroutine cwf(lmax, rho, eta, tol, max,f, g, fp, gp)
c
c      implicit real*8(a-h,o-z)
c      parameter ( nl = 1000 )
c      complex*16 ci, cgamma, cg, cz
c      dimension f(0:nl), g(0:nl), fp(0:nl), gp(0:nl)
c 
c      if (lmax.gt.nl) stop 'lmax is greater than the dimensions'
c 
c      ci = (0.d0,1.d0)
c      eta2 = eta*eta
c 
cc/// first get regular and irregular hypergeometric functions f and g
cc     and their derivatives fp and gp for l = 0 and 1
c 
c      do 100 l = 0,1
c 
c         cz = l + 1.d0 + ci*eta
c         cg = cgamma(cz)
c         sigmal = atan2 ( dimag(cg), dreal(cg) )
c 
c         call hyper(l, rho, eta, sigmal, tol, max,
c     &              fhyp, ghyp, dfhyp, dghyp)
c 
c         f(l) = fhyp
c         fp(l) = dfhyp
c         g(l) = ghyp
c         gp(l) = dghyp
c 
c100   continue
c 
cc/// recur for higher l's
c 
c      do 200 l=1, lmax-1
c 
c         lm1 = l - 1
c         lp1 = l + 1
c 
c         term1 = ( l + lp1 ) * ( eta + ( l * lp1 / rho ) )
c         term2 = lp1 * sqrt( l * l + eta2 )
c         term3 = 1.d0 / ( l * sqrt( lp1 * lp1 + eta2 ) )
c 
c         f(lp1) = ( term1 * f(l) - term2 * f(lm1) ) * term3
c         g(lp1) = ( term1 * g(l) - term2 * g(lm1) ) * term3
c 
c         term1 = sqrt( lp1*lp1 + eta2 ) / lp1
c         term2 = ( lp1*lp1 / rho + eta) / lp1
c 
c         fp(lp1) = term1 * f(l) - term2 * f(lp1)
c         gp(lp1) = term1 * g(l) - term2 * g(lp1)
c 
c         wronk = fp(lp1)*g(lp1) - f(lp1)*gp(lp1)
c         if ( abs(wronk-1.d0) .gt. tol) then
c            write(6,*) 'Warning: Wronskian is not correct, l= ',l
c            write(6,*) 'Wronskian (+1) = ',wronk
c            write(6,*)
c            stop
c         end if
c 
c200   continue
c 
c      return
c      end
c
cc/////////////////////////////////////////////////////////////////////// 
cc effective charge computed as the Fourier Transform of the potential 
cc  for use in the Born approximation
c
c      function zeff(p)
c      implicit real*8(a-h,o-z)
c      external fint
c      common /zef1/ pm
c      common /modpot/ zmod, xi, eta, rmod, zasy, nmod
c 
c      pm=p
c      rmin=1.d-10
c      rmax=rmod
c      call romba(rmin,rmod,fint,zef,1.d-5,ndim,ier)
c      if(ier.ne.0)then
c          write(6,*)'ndim=',ndim,'  ier=',ier
c          stop
c      endif
c      zeff = -pm * zef
c 
c      return
c      end
c 
cc///////////////////////////////////////////////////////////////////////
cc integrand for Fourier Transform integral used in the Born approximation
c
c      function fint(r)
c      implicit real*8(a-h,o-z)
c      common /zef1/ pm
c      common /modpot/ zmod, xi, eta, rmod, zasy, nmod
c        call pot(v,r)
c        fint = (r*v-zasy)*sin(pm*r)
c      return
c      end
c 
cc///////////////////////////////////////////////////////////////////////
cc Romberg integration routine
c
c      subroutine romba(xl,xu,cfct,cresul,er,ndim,ier)
c      implicit real*8(a-h,o-z)
c      external cfct
c      dimension c(20,20)
c      data c0/0.d0/
c 
c      h=xu-xl
c      i=1
c      h0=h
c      caux1= 0.5d0 *(cfct(xl)+cfct(xu))
c      c(1,1)=caux1
c      cr1=caux1
c 
c      i=2
c      hn=h0/2.d0
c      x=xl+hn
c      p=2.d0
c      caux2=(caux1+cfct(x))/p
c      c(1,2)=caux2
c      ert1= abs(caux2-caux1)/ abs(caux2)
c 
c      c(2,1)=c(1,2)+(c(1,2)-c(1,1))/3.d0
c      cr2=c(2,1)
c      err1= abs(cr2-cr1)/ abs(cr2)
c 
c  100 continue
c 
c      i=i+1
c      p=p*2.d0
c      ndim=i
c      h0=hn
c      hn= 0.5d0 *h0
c 
c      x=xl+hn
c      cs=c0
c      jj=2**(i-2)
c      do 3 l=1,jj
c      cs=cs+cfct(x)
c    3 x=x+h0
c 
c      caux3=0.5d0*caux2+cs/p
c      c(1,i)=caux3
c      ert2= abs(caux3-caux2)/ abs(caux3)
c      if(i.gt.19) go to 99
c 
c      jcross=i-1
c      q=1.d0
c      do 4 j=1,jcross
c      ij=i-j
c      ii=j+1
c      q=q+q
c      q=q+q
c      c(ii,ij)=c(ii-1,ij+1) + (c(ii-1,ij+1)-c(ii-1,ij))/(q-1.d0)
c    4 continue
c      cr3=c(i,1)
c      err2= abs(cr3-cr2)/ abs(cr3)
c      if(err1.lt.2.d0*er.and.err2.lt.er) go to 90
c 
c      err1=err2
c      cr2=cr3
c      ert1=ert2
c      caux2=caux3
c      go to 100
c 
c   90 continue
c      ier=0
c      cresul=cr3*h
c      return
c 
c   99 continue
c      ier=1
c      cresul=h*caux3
c      return
c 
c      end
c 
cc///////////////////////////////////////////////////////////////////////
cc read the input and set other related variables
c
c      subroutine readit
c
c      implicit real*8(a-h,o-z)
c      common veloc, en
c      common /modpot/ zmod, xi, eta, rmod, zasy, nmod
c      common /data  / rmu,rstart,rend0,spac,lmin,lmax,lspc
c      common /unit  / iout
c      common /born  / iborn
c
cc/// open read unit
c      open(unit=1,file='elastic.in',status='old')
c
cc/// electron energy in eV, convert to a.u.
c      read(1,*) en      
c      en=en/27.2
c
cc/// reduced mass
c      read(1,*) rmu
c
cc/// starting point for the numerical integration of the Schrodinger equation
c      read(1,*) rstart
c    
cc/// Ending point for the numerical integration of the Schrodinger equation
cc    rend0 should have two properties: (1) the potential at this point
cc    must reach its aymptotic form (e.g. COulomb or negligible) so that
cc    the wavefunction can be expressed as a linear combination of the
cc    regular and irregular Coulomb functions  (2) the value of rend0
cc    must be large enough so that the routine used to evaluate the
cc    hypergeometric functions is correct (note that only an asymptotic
cc    series is used).  Rule of thumb: make (k*rend0) > 10a.u.  (where
cc    k is the electron momentum).
c      read(1,*) rend0    
c 
cc/// step in a.u. for the numerical integration of the Scroedinger Eq
c      read(1,*) spac
c
cc/// minimum, maximum angular momenta and step
c      read(1,*) lmin, lmax, lspc
c
cc/// data for atomic potential of Garvey et al.
c      read(1,*) qmod, zmod, xi, eta, rmod
c      nmod = int(zmod-qmod+1.0001d0)
c      zasy =-(zmod-dfloat(nmod)+1.d0) + 1.d-10
c
cc/// use first Born approximation (iborn=1) or numerical integration of 
cc    Schrodinger equation (iborn=0)?
c      read(1,*) iborn
c
cc/// echo input
c      write(iout,*) 
c     &'Elastic cross section by the Johnson method for potentials ',
c     &'             of the Garvey et al. form'
c      write(iout,*) 
c     &'==========================================================='
c      write(iout,*)
c      write(iout,1000) en
c1000  format(' ','impact energy: ',1pe12.5,' a.u.')
c      write(iout,1001) rmu
c1001  format(' ','reduced mass: ',1pe12.5)
c      write(iout,1002) rstart,rend0,spac
c1002  format(' ','rstart: ',1pe12.5,' rend0: ',1pe12.5,
c     1        ' spac: ',1pe12.5,' a.u.')
c      write(iout,1003) lmin,lmax,lspc
c1003  format(' ','lmin: ',i3,' lmax: ',i3,' lspc: ',i3)
c      write(iout,*)
c      write(iout,1004) nmod,zmod,xi,eta
c1004  format(' ','nmod= ',i3,' zmod= ',1pe12.5,' xi= ',1pe12.5,
c     &           ' eta= ',1pe12.5)
c      write(iout,1005) rmod,qmod
c1005  format(' ','rmod= ',1pe12.5,' qmod= ',1pe12.5)
c      write(iout,*)
c 
c      close(1)
c      return
c      end
cc/////////////////////////////////////////////////////////////////////// 
c--CUT------------------------------------------------------------------
c1000.				! en, electron energy in eV	
c1.				! rmu, electron-projectile reduced mass
c0.				! rstart, beginning of radial integ. (a.u.)
c14.1 	  			! rend0, radius ending numerical integ. (a.u.)
c1.e-3				! space, numerical integration step (a.u.)
c0,100,1				! lmin, lmax, lspc
c0.,54.,1.044,5.101,14.1		! qmod, zmod, xi, eta, rmod
c0				! iborn, =0 numerical integ., =1 Born approx.
c--CUT------------------------------------------------------------------
cc///////////////////////////////////////////////////////////////////////
cc
cc This interactive program computes the independent particle model
cc  potential of Garvey, Jackman, and Green, Phys. Rev. A 12, 1144 (1975),
cc  which is based on the Hartree-Fock approach.  The resulting potential 
cc  for any (Hartree-Fock) ground state ion with nuclear charge Z.ge.2 and 
cc  Z.le.54 was parameterized by Garvey et al. and is tabulated with this 
cc  code by inputing two numbers; q, the ionic charge and Z, the nuclear
cc  charge of the ion or atom. Some examples of possible input are given 
cc  below.  Since the Garvey et al. fitting parameters vary smoothly, we 
cc  have extrapolated them to allow for ions and atoms with Z.le.100. 
cc  The user should be aware of the limitations of these approaches as 
cc  discussed and implied in the work by Garvey et al. and the accompanying 
cc  CPC article describing this program.
cc
cc The intended purpose of the program is to generate a plottable output
cc  file for displaying the Garvey et al. potential in comparison with
cc  the Coulomb potential for the inputted q.  Thus, comparison is made
cc  for a given ion or atom of the potential experienced by another charge
cc  outside a bare charge Z, or screened by the inputted number of electrons.
cc  In addition, comparison of these two potential outputs allows the user
cc  to determine the "range" of the screened potential.  This range is
cc  input to other programs, such as e_elastic.f, which computes the 
cc  elastic scattering cross section for scattering from a screened ion
cc  or from an atom.  Since the Garvey et al. potential formula contains
cc  an exponential, it is often numerically convenient to use this range
cc  to know when to switch from the Garvey et al. form to a Coulomb form.
cc
cc Inputs:
cc  q -- The ionic charge of the ion or atom (=0) under consideration
cc  Z -- The nuclear charge of the ion or atom under consideration
cc Outputs: (screen, unit 6)
cc  q, Z -- echo
cc  xi, eta -- Garvey et al. potential parameters
cc  range -- estimated range of the Garvey et al potential to which it 
cc           differs from the Coulomb potential
cc Outputs: (output file unit 7)
cc  r  -- radius (since it is a central potential, in atomic units
cc        the code has a default set of radial points at which the
cc        potentials will be tabulated
cc  vm -- the Garvey et al. potential at a given r, in atomic units
cc  vc -- the Coulomb potential at a given r, i.e. vc(r) = -q/r, in atomic
cc        units
cc  d  -- the difference vm(r) - vc(r), to help judge the "range" as
cc        discussed above (Note: program computes but does not print
cc        the derivative of the Garvey et al. program, so the user could
cc        add this to the output along with the derivative of the Coulomb
cc        potential, to have an additional measure of the convergence of
cc        the Garvey et al. and Coulomb potentials)
cc
cc  Examples: 
cc   1) for the potential experienced by an electron colliding with C^2+
cc      q=2., Z=6. (vm will asymptotically go to vc for a charge of 2)
cc   2) for the potential experienced by an electron colliding with Xe^22+
cc      q=22., Z=54.
cc   3) for the potential experienced by an electron colliding with He 
cc      q=0., Z=2. (vm will asymptotically be zero, vc will be zero; this
cc      illustrates an extrapolation of the Garvey et al. potential for 
cc      neutrals)
cc   4) Garvey et al.'s work does not include a potential for an electron
cc      colliding with atomic hydrogen.  We provide
cc      parameters which mimic the exact, electrostatic 
cc      potential of an electron colliding with the ground state of hydrogen
cc      (neglecting exchange). This user may wish to attempt to improve this.  
cc   5) For Z.gt.54 we provide parameters which are an extrapolation of
cc      the values provided by Garvey et al.
cc
cc D.R. Schultz and C.O. Reinhold, version 5/15/98
cc Physics Division, Oak Ridge National Laboratory
cc///////////////////////////////////////////////////////////////////////

      program green_pot

      implicit real*8 (a-h,o-z)
      dimension delr(5)

c/// these statement functions define the Garvey et al. potential, vm,
c     and its derivative dvm

      om(x) = 1. / ( (eta/xi) * (exp(xi*x)-1.) + 1. )
      zm(x) = ( float(n-1) * (1.-om(x)) - z )
      vm(x) = zm(x) / x

      dom(x) = -1. * om(x)**2 * eta * exp(xi*x)
      dzm(x) = -1. * (float(n-1)) * dom(x)
      dvm(x) = zm(x) * (-1./(x**2)) + dzm(x)/x

c/// open the output file

      open(unit=7,file='green-pot.out',status='unknown')

c/// take input from the console 

      write(6,*) 'Enter q, Z: '
      read(5,*) q,z
      n = int(z-q+1.0001d0)

c/// look up, or interpolate, to obtain the parameters of the Garvey
c     et al. potential, xi and eta

      call param(z,n,xi,eta)

c/// write parameters to screen

      write(6,*)
      write(6,*) 'q, Z, xi, eta:'
      write(6,*) q,z,xi,eta

c/// define the asymptotic Coulomb charge

      zc = z - n + 1

c/// define radial mesh for r in atomic units

      delr(1) = 0.01
      delr(2) = 0.1
      delr(3) = 0.5
      delr(4) = 1.
      delr(5) = 5.

c/// set tolerance for difference between vm and vc and their derivative 
c     at which range can be set

      tol = 1.e-6
      told = 1.e-6
      range = 10.
      iskip = 0

c/// loop over radial mesh points, print potentials, etc.

      write(7,*)
     &'   r        vm           vc           vm-vc        dvm-dvc'
      r = 0.
      do i=1,5
      do j=1,10
         r = r + delr(i)
	 diff = vm(r) + zc/r
	 diffd = dvm(r) - zc/r/r
c        write(6,11) r, vm(r), -zc/r, diff, diffd
         write(7,11) r, vm(r), -zc/r, diff, diffd
         if ((iskip.eq.0)
     &       .and.(abs(diff).le.tol).and.(abs(diffd).le.told)) then
	    range = r
	    iskip = 999
         end if
c         if (-vm(r).gt.1.e-60) write(7,*) r, -vm(r)
      end do
      end do

c/// write range estimate

      write(6,*)
      write(6,*) 'Estimated range:'
      write(6,*) range

c/// finished

      close(7)
 11   format(1x,f7.4,2x,1p,4(e11.4,2x),0p)
      stop
      end

c///////////////////////////////////////////////////////////////////////
      subroutine param(z,n,xi,eta)

      implicit real*8 (a-h,o-z)
      dimension xi0(2:54), xi1(2:54), eta0(2:54), eta1(2:54)

c/// fitting parameters from Garvey et al., note values for N= 38, 40,
c     43, 45, 47, 49, 51, and 53 have been linearly interpolated from
c     their Table 1 as they suggest 

      data xi0 /2.625, 2.164, 1.300, 1.031, 1.065, 1.179, 1.360, 
     &          1.508, 1.792, 1.712, 1.492, 1.170, 1.012, 0.954,
     &          0.926, 0.933, 0.957, 0.964, 0.941, 0.950, 0.998,
     &          1.061, 1.138, 1.207, 1.308, 1.397, 1.455, 1.520,
     &          1.538, 1.541, 1.512, 1.492, 1.460, 1.407, 1.351,
     &          1.286, 1.208, 1.129, 1.134, 1.139, 1.136, 1.167,
     &          1.197, 1.222, 1.246, 1.225, 1.205, 1.168, 1.130,
     &          1.090, 1.050, 1.047, 1.044/
      data xi1 /1.2996, 0.9764, 0.6465, 0.4924, 0.4800, 0.4677,
     &          0.4613, 0.4602, 0.4515, 0.3923, 0.3452, 0.3191,
     &          0.2933, 0.2659, 0.2478, 0.2368, 0.2165, 0.2151,
     &          0.2248, 0.2324, 0.2345, 0.2243, 0.2291, 0.2408,
     &          0.2391, 0.2462, 0.2397, 0.2246, 0.2106, 0.1988,
     &          0.1914, 0.1990, 0.1857, 0.1897, 0.1872, 0.1686,
     &          0.1745, 0.1784, 0.1743, 0.1702, 0.1694, 0.1648, 
     &          0.1601, 0.1594, 0.1587, 0.1493, 0.1358, 0.1377, 
     &          0.1395, 0.1374, 0.1354, 0.1231, 0.1107/
      data eta0 /1.770, 1.750, 1.880, 2.000, 2.130, 2.270, 2.410,
     &           2.590, 2.710, 2.850, 3.010, 3.170, 3.260, 3.330,
     &           3.392, 3.447, 3.500, 3.516, 3.570, 3.627, 3.667,
     &           3.709, 3.745, 3.803, 3.840, 3.891, 3.973, 4.000,
     &           4.050, 4.110, 4.182, 4.230, 4.290, 4.369, 4.418,
     &           4.494, 4.556, 4.618, 4.649, 4.680, 4.749, 4.759,
     &           4.769, 4.799, 4.829, 4.867, 4.904, 4.947, 4.990,
     &           5.020, 5.050, 5.076, 5.101/ 
      data eta1 /1.1402, 0.6821, 0.5547, 0.4939, 0.4434, 0.4143,
     &           0.3925, 0.3755, 0.3671, 0.3469, 0.3269, 0.3087,
     &           0.2958, 0.2857, 0.2739, 0.2633, 0.2560, 0.2509,
     &           0.2404, 0.2328, 0.2238, 0.2171, 0.2187, 0.2090,
     &           0.2088, 0.2048, 0.1925, 0.1985, 0.1878, 0.2001,
     &           0.1897, 0.1782, 0.1772, 0.1686, 0.1611, 0.1619,
     &           0.1564, 0.1509, 0.1497, 0.1485, 0.1412, 0.1424, 
     &           0.1435, 0.1416, 0.1397, 0.1406, 0.1414, 0.1369, 
     &           0.1324, 0.1319, 0.1314, 0.1315, 0.1316/

c/// data for extrapolated values

      data xi0_60/1.13/, xi0_70/1.10/, xi0_80/1.08/, xi0_90/1.05/, 
     &     xi0_100/1.02/
      data xi1_60/0.11/, xi1_70/0.082/, xi1_80/0.072/, xi1_90/0.060/,
     &     xi1_100/0.047/
      data eta0_60/5.35/, eta0_70/5.65/, eta0_80/5.83/, eta0_90/6.00/,
     &     eta0_100/6.15/
      data eta1_60/0.120/, eta1_70/0.112/, eta1_80/0.108/,
     &     eta1_90/0.102/, eta1_100/0.100/

c/// if N=2, Z=1 (electron colliding with atomic hydrogen), use parameters
c     chosen to mimic the exact electrostatic potential (allow some variance
c     in z so that the value entered numerically may not equal 1.

      if ((n.eq.2).and.(z.ge.0.95.and.z.le.1.05)) then
         xi = 1.78
         eta = 1.0
         goto 100
      end if

c/// if N.le.54 use tabulated values of Garvey et al.
     
      if (n.le.54) then
 	 zn = z-float(n)
         xi = xi0(n) + zn*xi1(n)
 	 eta = eta0(n) + zn*eta1(n)
         goto 100
      end if

c/// if N.gt.54 use visually extrapolated values of Garvey et al.
c     parameters (i.e. plot tabulated values and extrapolate). The
c     user should note the uncertainty in using this approach.

      if ((n.gt.54).and.(n.le.60)) then
         xn=float(n)
         xi0e=(xn-60.)/(-6.)*xi0(54) + (xn-54.)/6.*xi0_60
         xi1e=(xn-60.)/(-6.)*xi1(54) + (xn-54.)/6.*xi1_60
         eta0e=(xn-60.)/(-6.)*eta0(54) + (xn-54.)/6.*eta0_60
         eta1e=(xn-60.)/(-6.)*eta1(54) + (xn-54.)/6.*eta1_60
 	 zn = z-float(n)
         xi = xi0e + zn*xi1e
 	 eta = eta0e + zn*eta1e
         goto 100
      end if

      if ((n.gt.60).and.(n.le.70)) then
         xn=float(n)
         xi0e=(xn-70.)/(-10.)*xi0_60 + (xn-60.)/10.*xi0_70
         xi1e=(xn-70.)/(-10.)*xi1_60 + (xn-60.)/10.*xi1_70
         eta0e=(xn-70.)/(-10.)*eta0_60 + (xn-60.)/10.*eta0_70
         eta1e=(xn-70.)/(-10.)*eta1_60 + (xn-60.)/10.*eta1_70
 	 zn = z-float(n)
         xi = xi0e + zn*xi1e
 	 eta = eta0e + zn*eta1e
         goto 100
      end if

      if ((n.gt.70).and.(n.le.80)) then
         xn=float(n)
         xi0e=(xn-80.)/(-10.)*xi0_70 + (xn-70.)/10.*xi0_80
         xi1e=(xn-80.)/(-10.)*xi1_70 + (xn-70.)/10.*xi1_80
         eta0e=(xn-80.)/(-10.)*eta0_70 + (xn-70.)/10.*eta0_80
         eta1e=(xn-80.)/(-10.)*eta1_70 + (xn-70.)/10.*eta1_80
 	 zn = z-float(n)
         xi = xi0e + zn*xi1e
 	 eta = eta0e + zn*eta1e
         goto 100
      end if

      if ((n.gt.80).and.(n.le.90)) then
         xn=float(n)
         xi0e=(xn-90.)/(-10.)*xi0_80 + (xn-80.)/10.*xi0_90
         xi1e=(xn-90.)/(-10.)*xi1_80 + (xn-80.)/10.*xi1_90
         eta0e=(xn-90.)/(-10.)*eta0_80 + (xn-80.)/10.*eta0_90
         eta1e=(xn-90.)/(-10.)*eta1_80 + (xn-80.)/10.*eta1_90
 	 zn = z-float(n)
         xi = xi0e + zn*xi1e
 	 eta = eta0e + zn*eta1e
         goto 100
      end if

      if ((n.gt.90).and.(n.le.100)) then
         xn=float(n)
         xi0e=(xn-100.)/(-10.)*xi0_90 + (xn-90.)/10.*xi0_100
         xi1e=(xn-100.)/(-10.)*xi1_90 + (xn-90.)/10.*xi1_100
         eta0e=(xn-100.)/(-10.)*eta0_90 + (xn-90.)/10.*eta0_100
         eta1e=(xn-100.)/(-10.)*eta1_90 + (xn-90.)/10.*eta1_100
 	 zn = z-float(n)
         xi = xi0e + zn*xi1e
 	 eta = eta0e + zn*eta1e
         goto 100
      end if

c/// if N.gt.100 write message

      write(6,*) 'N chosen to be greater than 100, no parameters 
     & available'

c/// return here

 100  return
      end
                                                                            ****
