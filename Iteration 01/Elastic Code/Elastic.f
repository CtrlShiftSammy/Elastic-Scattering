cADJHELASTIC, GREEN-POT.  ELASTIC - ELASTIC SCATTERING OF ELECTRONS FROM ADJH0000
c1   IONS AND ATOMS.  D.R. SCHULTZ, C.O. REINHOLD.                       ADJH0000
c///////////////////////////////////////////////////////////////////////
c    Calculation of elastic scattering cross sections of electrons 
c                      at ions and atoms
c
c The interaction potential is contained in the routine ''pot'' and
c can be set to any desired form which rises less rapidly than the 
c centrifugal potential at the origin. Presently implemented is the
c parameterization of the potential following Garvey, Jackman, and
c Green (Phys. Rev. 12, 1144 (1975)).  Use the auxiliary code Green_pot
c to determine these parameters, or use the tabulation and formulae given
c by Garvey et al. to compute them.
c
c The asymptotic behavior of the wavefunction at large distances is
c matched to Coulomb functions to include treated scattering from
c ions.  If the Coulomb functions diverge, increase the input
c variable ''rend0'' (matching radius).  For zero asymptotic charges
c (atoms), the Coulomb functions authomatically behave as the 
c appropriate Bessel functions.
c                                    
c The integration of the radial Schrodinger equation is performed based 
c on Johnson's log derivative method (J. Comput. Phys. 13, 445 (1973)). 
c                                                   
c The evaluation of the phase shifts is vectorized.
c
c See routine readit at the bottom for a brief explanation of the 
c input data.
c
c D.R. Schultz and C.O. Reinhold, version 5/15/98
c Physics Division, Oak Ridge National Laboratory
c///////////////////////////////////////////////////////////////////////

      program elastic

      implicit real*8(a-h,o-z)
      parameter (nl = 1000)
      complex*16 cz, cgamma, cg, ci, cratio, cdexpon, cfc, cfnc, cf
      real*8 logder
      dimension fl(nl), z11(nl), phase(nl), sigma(nl)
      dimension f(0:nl), g(0:nl), fp(0:nl), gp(0:nl)
      common veloc, en
      common /modpot/ zmod, xi, eta, rmod, zasy, nmod
      common /data  / rmu,rstart,rend0,spac,lmin,lmax,lspc
      common /unit  / iout
      common /born  / iborn

c/// open output  
      iout=6
      open(iout,file='elastic.out',status='unknown')

c/// read and echo input
      call readit
 
c/// definitions 
      pi = 3.14159265358979323846d0
      ci = (0.d0,1.d0)
 
      veloc = sqrt(2.d0*en/rmu)
      rmu2 = 2.d0*rmu
      p1 =  sqrt(en*rmu2)
 
      etahyp = zasy / p1 + 1.e-30
 
      con = 0.125d0
      lmin1 = lmin + 1
      lmmx = lmax + 1
 
      nos = 2*int ((rend0-rstart)/2.d0/spac)+2
      space = (rend0-rstart)/ dfloat(nos)
      del1 = space*space*rmu2/3.d0
      del2 = 2.d0*del1
      del4 = 4.d0*del1
 
      do  ll = lmin1,lmmx,lspc
         l = ll-1
         fl(ll) =  dfloat (l*(l+1))/rmu2
         z11(ll) = 1.e20
      end do
 
c/// Integrate radial Schrodinger equation for log derivative
      do ll = lmin1,lmmx,lspc
         r=rstart
         do nr = 1,nos,2
            r=r+space
            call pot(v,r)
            v11=(en-v)*del4
            cen=1.d0/(r**2)*del4
            r=r+space
            call pot(v,r)
            vv11=(en-v)*del2
            ccen=1.d0/(r**2)*del2
            y11=z11(ll)
            w11=v11-cen*fl(ll)
            den=1.d0+con*w11
            u11=w11/den
            den=1.d0+y11
            y11=y11/den-u11
            u11=vv11-ccen*fl(ll)
            den=1.d0+y11
            z11(ll)=y11/den-u11
         end do
      end do 
 
      do ll=lmin1,lmmx,lspc
         u11=vv11-ccen*fl(ll)
         z11(ll)=(z11(ll)+u11/2.d0)/space
      end do
  
c/// calculation of the phase shifts 
      rho = rend0 * p1
      tol = 1.0e-04
      max = 50
 
      call cwf(lmax, rho, etahyp, tol, max, f, g, fp, gp)
 
      hcon = 4.d0*pi/p1/p1
      tcs=0.d0
      do ll=lmin1,lmmx,lspc
         l = ll - 1
         logder = z11(ll) / p1
         fhyp = f(l)
         ghyp = g(l)
         dfhyp = fp(l)
         dghyp = gp(l)
         phase(ll) = atan( (logder*fhyp-dfhyp)/(dghyp-ghyp*logder) )
         cg=(1.d0,0.d0) 		
         cz = l + 1.d0 + ci*etahyp
         cg = cgamma(cz)  
         sigma(ll) = atan2 ( dimag(cg), dreal(cg) )
         write(iout,6999) l, phase(ll), sigma(ll)
         tcs = tcs + hcon * (2.d0*l+1) * sin(phase(ll))**2
      end do
 
6999  format(' ','l = ',i3,2x,'phase shift =',1pe15.7,2x,
     @        'Coulomb phase shift=',1pe15.7)

c/// total cross section for the case of neutral targets (i.e. atoms)
      if (abs(zasy).lt.1.e-5) then
      write(iout,*)
      write(iout,*)'*************************************'
      write(iout,*)'*** TCS = ', tcs, ' a.u.'
      write(iout,*)'*************************************'
      end if
 
c/// differential cross sections 
      degrad = pi/180.
      angdeg = -0.999
 
      write(iout,*)
      write(iout,7012)
      if(iborn.eq.1)then
         write(iout,7011)
      else
         write(iout,7010)
      endif
7012  format(//,10X,'*** ELASTIC CROSS SECTIONS ***',/)
7011  format(1x,'theta(deg)    EXACT',11x,'BORN',' a.u.',/,50('='))
7010  format(1x,'theta(deg)   dcs(a.u.)     Rutherford-q  ',
     1          '  Rutherford-Z      ',/, 60('='))
 
      do i=0,180
         angdeg = angdeg + 1.d0
         angrad = angdeg * degrad
         sint22 = sin(angrad/2.d0)**2
         ruther = -zasy / ( 2.d0 * p1*p1 * sint22 )
         ruthzm =  zmod / ( 2.d0 * p1*p1 * sint22 )
         cratio = cgamma( 1.d0+ci*etahyp ) / cgamma( 1.d0-ci*etahyp )
         cdexpon = cdexp( -ci*etahyp*dlog(sint22) )
         cfc = ruther * cratio * cdexpon
         cost = cos(angrad)
         cfnc = (0.d0,0.d0)
         do ll=lmin1,lmmx,lspc
            l = ll - 1
            cfnc = cfnc + (2*l+1.d0) * cdexp(2.d0*ci*sigma(ll)) *
     &                    (cdexp(2.d0*ci*phase(ll))-1.d0) * pl(l,cost)
         end do
         cfnc = cfnc / (2.d0*ci*p1)
         cf = cfc + cfnc
         dcs = cdabs(cf)**2
         if(iborn.eq.1)then
            pmod = 2.d0*p1*sin(angrad/2.d0)
            dcsb1 = 4.d0*((-zasy+zeff(pmod))/pmod**2)**2
            write(iout,8001) angdeg, dcs, dcsb1
         else
            write(iout,8001) angdeg, dcs, ruther**2, ruthzm**2
         endif
      end do
8001  format(1x,f6.1,6x,1pe10.4,6x,e10.4,6x,e10.4)

      close(iout)
      stop
      end
 
c///////////////////////////////////////////////////////////////////////
c Central potential for elastic scattering.  Implemented for the
c Garvey et al. parameterized Hartree-Fock model potential. 

      subroutine pot(v,r)

      implicit real*8(a-h,o-z)
      common /modpot/ zmod, xi, eta, rmod, zasy, nmod
      dimension fact(8)
      data alp1,alp2,d2/1.38414,2.32,0.11578/
      data r0,hla/2.373,10.6172/
      data fact/1.,2.,6.,24.,120.,720.,5040.,40320./
 
      o(r) = 1.d0/((eta/xi)*exp(xi*r)-(eta/xi)+1.d0)
      z(r) = ( dfloat(nmod-1) * (1.d0-o(r)) ) - zmod
      vmod(r) = z(r) / r
      if (r.lt.rmod) then
         vtemp = vmod(r)
      else
         vtemp = zasy / r
      end if

      v = vtemp 
 
      return
      end
 
c/////////////////////////////////////////////////////////////////////// 
c This function subprogram evaluates the complex gamma function,
c using Stirling's formula (after algorithm in Numerical Recipes, W.H. 
c Press et al. (Cambridge University Press, Cambridge, 1989)

      function cgamma(cz)

      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      real*8 cof(6)
      data cof / 76.18009173d0, -86.50532033d0, 24.01409822d0,
     &           -1.231739516d0, 1.20858003d-3, -5.36382d-6 /
      data stp / 2.50662827465d0 /
 
c/// if 0 < Re(z) < 1 use reflection formula
 
      if ( dreal(cz).lt.1. ) then
         cx = cz
         ctmp = (cx+5.5d0)**(cx+0.5d0) * cdexp(-(cx+5.5d0))
         cser = 1.d0
         do 100 i = 1,6
            cx = cx + 1.d0
            cser = cser + cof(i)/cx
100      continue
         cgamma = ctmp * stp * cser / cz
      else
         cx = cz - 1.d0
         ctmp = (cx+5.5d0)**(cx+0.5d0) * cdexp(-(cx+5.5d0))
         cser = 1.d0
         do 200 i = 1,6
            cx = cx + 1.d0
            cser = cser + cof(i)/cx
200      continue
         cgamma = ctmp * stp * cser
      end if
 
      return
      end
 
c/////////////////////////////////////////////////////////////////////// 
c Legendre polynomial of order lin and argument x

      function pl(lin, x)

      implicit real*8(a-h,o-z)
 
      p0 = 1.d0
      p1 = x
 
      if(lin.eq.0) then
         pl = p0
         return
      end if
      if(lin.eq.1) then
         pl = p1
         return
      end if
 
      plm2 = p0
      plm1 = p1
 
      do 100 l=2,lin
 
         pl = ( (2*l-1)*x*plm1 - (l-1)*plm2 ) / l
         plm2 = plm1
         plm1 = pl
 
100   continue
 
      return
      end
 
c/////////////////////////////////////////////////////////////////////// 
c this subroutine evaluates the regular and irregular hypergeometric
c functions in the asymptotic regime, after the method of Abramowitz,
c as well as their derivatives (see Abramowitz and Stegun, Handbook of
c Mathematical Functions, (Dover, New York, 1972) chapter 14.5)

      subroutine hyper(l, rho, eta, sigmal, tol, max,
     &                 fhyp, ghyp, dfhyp, dghyp)
      implicit real*8(a-h,o-z)
 
      pi = 3.14159265358979323846d0
 
      thetal = rho - eta*dlog(2*rho) - l*pi*0.5d0 + sigmal
 
      f1k = 1.d0
      g1k = 0.d0
      f2k = 0.d0
      g2k = 1.d0 - eta/rho
 
      f1 = f1k
      g1 = g1k
      f2 = f2k
      g2 = g2k
 
      kount = 0
 
      ak = eta * 0.5d0 / rho
      bk = (l*(l+1) + eta*eta) * 0.5d0 / rho
 
      kmin = 1
      kmax = 10
 
1     f1old = f1
      g1old = g1
      f2old = f2
      g2old = g2
 
      kount = kount + 1
 
      do 100 k=kmin,kmax
 
         f1kp1 = ak*f1k - bk*g1k
         g1kp1 = ak*g1k + bk*f1k
         f2kp1 = ak*f2k - bk*g2k - f1kp1/rho
         g2kp1 = ak*g2k + bk*f2k - g1kp1/rho
 
         f1 = f1 + f1kp1
         g1 = g1 + g1kp1
         f2 = f2 + f2kp1
         g2 = g2 + g2kp1
 
         ak = (2.d0*k+1) * eta / (2.d0*k+2.d0) / rho
         bk = (l*(l+1) - k*(k+1) + eta*eta) / (2.d0*k+2.d0) / rho
 
         f1k = f1kp1
         g1k = g1kp1
         f2k = f2kp1
         g2k = g2kp1
 
100   continue

c///diagnostics 
c     errf1 = abs( (f1 - f1old) / f1 )
c     errg1 = abs( (g1 - g1old) / g1 )
c     errf2 = abs( (f2 - f2old) / f2 )
c     errg2 = abs( (g2 - g2old) / g2 )
c
c     errmax = amax1(errf1, errg1, errf2, errg2)
 
      wronk = g1*f2-f1*g2
 
      errmax = abs ( wronk + 1.d0 )
 
      if ( errmax.gt.tol ) then
         if ( kount.gt.max ) then
            write(6,*) 'Warning: Hypergeometric convergence failure',
     &            ' kount=',kount
            write(6,*) 'Wronskian (-1) = ',wronk
            write(6,*)
            stop
         end if
         kmin=kmin+10
         kmax=kmax+10
         goto 1
      end if
 
      s = sin(thetal)
      c = cos(thetal)
 
      fhyp = g1 * c + f1 * s
      ghyp = f1 * c - g1 * s
      dfhyp = g2 * c + f2 * s
      dghyp = f2 * c - g2 * s
 
      return
      end

c/////////////////////////////////////////////////////////////////////// 
c This subroutine computes the regular and irregular Coulomb wave  
c functions and their derivatives for l's from zero to lmax. The first   
c two are calculated and then the rest found by recursion.  (The method
c is from Abramowitz and Stegun,Handbook of Mathematical Functions, 
c (Dover, New York, 1972)).   

      subroutine cwf(lmax, rho, eta, tol, max,f, g, fp, gp)

      implicit real*8(a-h,o-z)
      parameter ( nl = 1000 )
      complex*16 ci, cgamma, cg, cz
      dimension f(0:nl), g(0:nl), fp(0:nl), gp(0:nl)
 
      if (lmax.gt.nl) stop 'lmax is greater than the dimensions'
 
      ci = (0.d0,1.d0)
      eta2 = eta*eta
 
c/// first get regular and irregular hypergeometric functions f and g
c     and their derivatives fp and gp for l = 0 and 1
 
      do 100 l = 0,1
 
         cz = l + 1.d0 + ci*eta
         cg = cgamma(cz)
         sigmal = atan2 ( dimag(cg), dreal(cg) )
 
         call hyper(l, rho, eta, sigmal, tol, max,
     &              fhyp, ghyp, dfhyp, dghyp)
 
         f(l) = fhyp
         fp(l) = dfhyp
         g(l) = ghyp
         gp(l) = dghyp
 
100   continue
 
c/// recur for higher l's
 
      do 200 l=1, lmax-1
 
         lm1 = l - 1
         lp1 = l + 1
 
         term1 = ( l + lp1 ) * ( eta + ( l * lp1 / rho ) )
         term2 = lp1 * sqrt( l * l + eta2 )
         term3 = 1.d0 / ( l * sqrt( lp1 * lp1 + eta2 ) )
 
         f(lp1) = ( term1 * f(l) - term2 * f(lm1) ) * term3
         g(lp1) = ( term1 * g(l) - term2 * g(lm1) ) * term3
 
         term1 = sqrt( lp1*lp1 + eta2 ) / lp1
         term2 = ( lp1*lp1 / rho + eta) / lp1
 
         fp(lp1) = term1 * f(l) - term2 * f(lp1)
         gp(lp1) = term1 * g(l) - term2 * g(lp1)
 
         wronk = fp(lp1)*g(lp1) - f(lp1)*gp(lp1)
         if ( abs(wronk-1.d0) .gt. tol) then
            write(6,*) 'Warning: Wronskian is not correct, l= ',l
            write(6,*) 'Wronskian (+1) = ',wronk
            write(6,*)
            stop
         end if
 
200   continue
 
      return
      end

c/////////////////////////////////////////////////////////////////////// 
c effective charge computed as the Fourier Transform of the potential 
c  for use in the Born approximation

      function zeff(p)
      implicit real*8(a-h,o-z)
      external fint
      common /zef1/ pm
      common /modpot/ zmod, xi, eta, rmod, zasy, nmod
 
      pm=p
      rmin=1.d-10
      rmax=rmod
      call romba(rmin,rmod,fint,zef,1.d-5,ndim,ier)
      if(ier.ne.0)then
          write(6,*)'ndim=',ndim,'  ier=',ier
          stop
      endif
      zeff = -pm * zef
 
      return
      end
 
c///////////////////////////////////////////////////////////////////////
c integrand for Fourier Transform integral used in the Born approximation

      function fint(r)
      implicit real*8(a-h,o-z)
      common /zef1/ pm
      common /modpot/ zmod, xi, eta, rmod, zasy, nmod
        call pot(v,r)
        fint = (r*v-zasy)*sin(pm*r)
      return
      end
 
c///////////////////////////////////////////////////////////////////////
c Romberg integration routine

      subroutine romba(xl,xu,cfct,cresul,er,ndim,ier)
      implicit real*8(a-h,o-z)
      external cfct
      dimension c(20,20)
      data c0/0.d0/
 
      h=xu-xl
      i=1
      h0=h
      caux1= 0.5d0 *(cfct(xl)+cfct(xu))
      c(1,1)=caux1
      cr1=caux1
 
      i=2
      hn=h0/2.d0
      x=xl+hn
      p=2.d0
      caux2=(caux1+cfct(x))/p
      c(1,2)=caux2
      ert1= abs(caux2-caux1)/ abs(caux2)
 
      c(2,1)=c(1,2)+(c(1,2)-c(1,1))/3.d0
      cr2=c(2,1)
      err1= abs(cr2-cr1)/ abs(cr2)
 
  100 continue
 
      i=i+1
      p=p*2.d0
      ndim=i
      h0=hn
      hn= 0.5d0 *h0
 
      x=xl+hn
      cs=c0
      jj=2**(i-2)
      do l=1,jj
      cs=cs+cfct(x)
      x=x+h0
      end do

      caux3=0.5d0*caux2+cs/p
      c(1,i)=caux3
      ert2= abs(caux3-caux2)/ abs(caux3)
      if(i.gt.19) go to 99
 
      jcross=i-1
      q=1.d0
      do 4 j=1,jcross
      ij=i-j
      ii=j+1
      q=q+q
      q=q+q
      c(ii,ij)=c(ii-1,ij+1) + (c(ii-1,ij+1)-c(ii-1,ij))/(q-1.d0)
    4 continue
      cr3=c(i,1)
      err2= abs(cr3-cr2)/ abs(cr3)
      if(err1.lt.2.d0*er.and.err2.lt.er) go to 90
 
      err1=err2
      cr2=cr3
      ert1=ert2
      caux2=caux3
      go to 100
 
   90 continue
      ier=0
      cresul=cr3*h
      return
 
   99 continue
      ier=1
      cresul=h*caux3
      return
 
      end
 
c///////////////////////////////////////////////////////////////////////
c read the input and set other related variables

      subroutine readit

      implicit real*8(a-h,o-z)
      common veloc, en
      common /modpot/ zmod, xi, eta, rmod, zasy, nmod
      common /data  / rmu,rstart,rend0,spac,lmin,lmax,lspc
      common /unit  / iout
      common /born  / iborn

c/// open read unit
      open(unit=1,file='elastic.in',status='old')

c/// electron energy in eV, convert to a.u.
      read(1,*) en      
      en=en/27.2

c/// reduced mass
      read(1,*) rmu

c/// starting point for the numerical integration of the Schrodinger equation
      read(1,*) rstart
    
c/// Ending point for the numerical integration of the Schrodinger equation
c    rend0 should have two properties: (1) the potential at this point
c    must reach its aymptotic form (e.g. Coulomb or negligible) so that
c    the wavefunction can be expressed as a linear combination of the
c    regular and irregular Coulomb functions  (2) the value of rend0
c    must be large enough so that the routine used to evaluate the
c    hypergeometric functions is correct (note that only an asymptotic
c    series is used).  Rule of thumb: make (k*rend0) > 10a.u.  (where
c    k is the electron momentum).
      read(1,*) rend0    
 
c/// step in a.u. for the numerical integration of the Scroedinger Eq
      read(1,*) spac

c/// minimum, maximum angular momenta and step
      read(1,*) lmin, lmax, lspc

c/// data for atomic potential of Garvey et al.
      read(1,*) qmod, zmod, xi, eta, rmod
      nmod = int(zmod-qmod+1.0001d0)
      zasy =-(zmod-dfloat(nmod)+1.d0) + 1.d-10

c/// use first Born approximation (iborn=1) or numerical integration of 
c    Schrodinger equation (iborn=0)?
      read(1,*) iborn

c/// echo input
      write(iout,*) 
     &'Elastic cross section by the Johnson method for potentials ',
     &'             of the Garvey et al. form'
      write(iout,*) 
     &'==========================================================='
      write(iout,*)
      write(iout,1000) en
1000  format(' ','impact energy: ',1pe12.5,' a.u.')
      write(iout,1001) rmu
1001  format(' ','reduced mass: ',1pe12.5)
      write(iout,1002) rstart,rend0,spac
1002  format(' ','rstart: ',1pe12.5,' rend0: ',1pe12.5,
     1        ' spac: ',1pe12.5,' a.u.')
      write(iout,1003) lmin,lmax,lspc
1003  format(' ','lmin: ',i3,' lmax: ',i3,' lspc: ',i3)
      write(iout,*)
      write(iout,1004) nmod,zmod,xi,eta
1004  format(' ','nmod= ',i3,' zmod= ',1pe12.5,' xi= ',1pe12.5,
     &           ' eta= ',1pe12.5)
      write(iout,1005) rmod,qmod
1005  format(' ','rmod= ',1pe12.5,' qmod= ',1pe12.5)
      write(iout,*)
 
      close(1)
      return
      end
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
c///////////////////////////////////////////////////////////////////////
c
c This interactive program computes the independent particle model
c  potential of Garvey, Jackman, and Green, Phys. Rev. A 12, 1144 (1975),
c  which is based on the Hartree-Fock approach.  The resulting potential 
c  for any (Hartree-Fock) ground state ion with nuclear charge Z.ge.2 and 
c  Z.le.54 was parameterized by Garvey et al. and is tabulated with this 
c  code by inputing two numbers; q, the ionic charge and Z, the nuclear
c  charge of the ion or atom. Some examples of possible input are given 
c  below.  Since the Garvey et al. fitting parameters vary smoothly, we 
c  have extrapolated them to allow for ions and atoms with Z.le.100. 
c  The user should be aware of the limitations of these approaches as 
c  discussed and implied in the work by Garvey et al. and the accompanying 
c  CPC article describing this program.
c
c The intended purpose of the program is to generate a plottable output
c  file for displaying the Garvey et al. potential in comparison with
c  the Coulomb potential for the inputted q.  Thus, comparison is made
c  for a given ion or atom of the potential experienced by another charge
c  outside a bare charge Z, or screened by the inputted number of electrons.
c  In addition, comparison of these two potential outputs allows the user
c  to determine the "range" of the screened potential.  This range is
c  input to other programs, such as e_elastic.f, which computes the 
c  elastic scattering cross section for scattering from a screened ion
c  or from an atom.  Since the Garvey et al. potential formula contains
c  an exponential, it is often numerically convenient to use this range
c  to know when to switch from the Garvey et al. form to a Coulomb form.
c
c Inputs:
c  q -- The ionic charge of the ion or atom (=0) under consideration
c  Z -- The nuclear charge of the ion or atom under consideration
c Outputs: (screen, unit 6)
c  q, Z -- echo
c  xi, eta -- Garvey et al. potential parameters
c  range -- estimated range of the Garvey et al potential to which it 
c           differs from the Coulomb potential
c Outputs: (output file unit 7)
c  r  -- radius (since it is a central potential, in atomic units
c        the code has a default set of radial points at which the
c        potentials will be tabulated
c  vm -- the Garvey et al. potential at a given r, in atomic units
c  vc -- the Coulomb potential at a given r, i.e. vc(r) = -q/r, in atomic
c        units
c  d  -- the difference vm(r) - vc(r), to help judge the "range" as
c        discussed above (Note: program computes but does not print
c        the derivative of the Garvey et al. program, so the user could
c        add this to the output along with the derivative of the Coulomb
c        potential, to have an additional measure of the convergence of
c        the Garvey et al. and Coulomb potentials)
c
c  Examples: 
c   1) for the potential experienced by an electron colliding with C^2+
c      q=2., Z=6. (vm will asymptotically go to vc for a charge of 2)
c   2) for the potential experienced by an electron colliding with Xe^22+
c      q=22., Z=54.
c   3) for the potential experienced by an electron colliding with He 
c      q=0., Z=2. (vm will asymptotically be zero, vc will be zero; this
c      illustrates an extrapolation of the Garvey et al. potential for 
c      neutrals)
c   4) Garvey et al.'s work does not include a potential for an electron
c      colliding with atomic hydrogen.  We provide
c      parameters which mimic the exact, electrostatic 
c      potential of an electron colliding with the ground state of hydrogen
c      (neglecting exchange). This user may wish to attempt to improve this.  
c   5) For Z.gt.54 we provide parameters which are an extrapolation of
c      the values provided by Garvey et al.
c
c D.R. Schultz and C.O. Reinhold, version 5/15/98
c Physics Division, Oak Ridge National Laboratory
c///////////////////////////////////////////////////////////////////////
c
c      program green_pot
c
c      implicit real*8 (a-h,o-z)
c      dimension delr(5)
c
cc/// these statement functions define the Garvey et al. potential, vm,
cc     and its derivative dvm
c
c      om(x) = 1. / ( (eta/xi) * (exp(xi*x)-1.) + 1. )
c      zm(x) = ( float(n-1) * (1.-om(x)) - z )
c      vm(x) = zm(x) / x
c
c      dom(x) = -1. * om(x)**2 * eta * exp(xi*x)
c      dzm(x) = -1. * (float(n-1)) * dom(x)
c      dvm(x) = zm(x) * (-1./(x**2)) + dzm(x)/x
c
cc/// open the output file
c
c      open(unit=7,file='green-pot.out',status='unknown')
c
cc/// take input from the console 
c
c      write(6,*) 'Enter q, Z: '
c      read(5,*) q,z
c      n = int(z-q+1.0001d0)
c
cc/// look up, or interpolate, to obtain the parameters of the Garvey
cc     et al. potential, xi and eta
c
c      call param(z,n,xi,eta)
c
cc/// write parameters to screen
c
c      write(6,*)
c      write(6,*) 'q, Z, xi, eta:'
c      write(6,*) q,z,xi,eta
c
cc/// define the asymptotic Coulomb charge
c
c      zc = z - n + 1
c
cc/// define radial mesh for r in atomic units
c
c      delr(1) = 0.01
c      delr(2) = 0.1
c      delr(3) = 0.5
c      delr(4) = 1.
c      delr(5) = 5.
c
cc/// set tolerance for difference between vm and vc and their derivative 
cc     at which range can be set
c
c      tol = 1.e-6
c      told = 1.e-6
c      range = 10.
c      iskip = 0
c
cc/// loop over radial mesh points, print potentials, etc.
c
c      write(7,*)
c     &'   r        vm           vc           vm-vc        dvm-dvc'
c      r = 0.
c      do i=1,5
c      do j=1,10
c         r = r + delr(i)
c	 diff = vm(r) + zc/r
c	 diffd = dvm(r) - zc/r/r
cc        write(6,11) r, vm(r), -zc/r, diff, diffd
c         write(7,11) r, vm(r), -zc/r, diff, diffd
c         if ((iskip.eq.0)
c     &       .and.(abs(diff).le.tol).and.(abs(diffd).le.told)) then
c	    range = r
c	    iskip = 999
c         end if
cc         if (-vm(r).gt.1.e-60) write(7,*) r, -vm(r)
c      end do
c      end do
c
cc/// write range estimate
c
c      write(6,*)
c      write(6,*) 'Estimated range:'
c      write(6,*) range
c
cc/// finished
c
c      close(7)
c 11   format(1x,f7.4,2x,1p,4(e11.4,2x),0p)
c      stop
c      end
c
cc///////////////////////////////////////////////////////////////////////
c      subroutine param(z,n,xi,eta)
c
c      implicit real*8 (a-h,o-z)
c      dimension xi0(2:54), xi1(2:54), eta0(2:54), eta1(2:54)
c
cc/// fitting parameters from Garvey et al., note values for N= 38, 40,
cc     43, 45, 47, 49, 51, and 53 have been linearly interpolated from
cc     their Table 1 as they suggest 
c
c      data xi0 /2.625, 2.164, 1.300, 1.031, 1.065, 1.179, 1.360, 
c     &          1.508, 1.792, 1.712, 1.492, 1.170, 1.012, 0.954,
c     &          0.926, 0.933, 0.957, 0.964, 0.941, 0.950, 0.998,
c     &          1.061, 1.138, 1.207, 1.308, 1.397, 1.455, 1.520,
c     &          1.538, 1.541, 1.512, 1.492, 1.460, 1.407, 1.351,
c     &          1.286, 1.208, 1.129, 1.134, 1.139, 1.136, 1.167,
c     &          1.197, 1.222, 1.246, 1.225, 1.205, 1.168, 1.130,
c     &          1.090, 1.050, 1.047, 1.044/
c      data xi1 /1.2996, 0.9764, 0.6465, 0.4924, 0.4800, 0.4677,
c     &          0.4613, 0.4602, 0.4515, 0.3923, 0.3452, 0.3191,
c     &          0.2933, 0.2659, 0.2478, 0.2368, 0.2165, 0.2151,
c     &          0.2248, 0.2324, 0.2345, 0.2243, 0.2291, 0.2408,
c     &          0.2391, 0.2462, 0.2397, 0.2246, 0.2106, 0.1988,
c     &          0.1914, 0.1990, 0.1857, 0.1897, 0.1872, 0.1686,
c     &          0.1745, 0.1784, 0.1743, 0.1702, 0.1694, 0.1648, 
c     &          0.1601, 0.1594, 0.1587, 0.1493, 0.1358, 0.1377, 
c     &          0.1395, 0.1374, 0.1354, 0.1231, 0.1107/
c      data eta0 /1.770, 1.750, 1.880, 2.000, 2.130, 2.270, 2.410,
c     &           2.590, 2.710, 2.850, 3.010, 3.170, 3.260, 3.330,
c     &           3.392, 3.447, 3.500, 3.516, 3.570, 3.627, 3.667,
c     &           3.709, 3.745, 3.803, 3.840, 3.891, 3.973, 4.000,
c     &           4.050, 4.110, 4.182, 4.230, 4.290, 4.369, 4.418,
c     &           4.494, 4.556, 4.618, 4.649, 4.680, 4.749, 4.759,
c     &           4.769, 4.799, 4.829, 4.867, 4.904, 4.947, 4.990,
c     &           5.020, 5.050, 5.076, 5.101/ 
c      data eta1 /1.1402, 0.6821, 0.5547, 0.4939, 0.4434, 0.4143,
c     &           0.3925, 0.3755, 0.3671, 0.3469, 0.3269, 0.3087,
c     &           0.2958, 0.2857, 0.2739, 0.2633, 0.2560, 0.2509,
c     &           0.2404, 0.2328, 0.2238, 0.2171, 0.2187, 0.2090,
c     &           0.2088, 0.2048, 0.1925, 0.1985, 0.1878, 0.2001,
c     &           0.1897, 0.1782, 0.1772, 0.1686, 0.1611, 0.1619,
c     &           0.1564, 0.1509, 0.1497, 0.1485, 0.1412, 0.1424, 
c     &           0.1435, 0.1416, 0.1397, 0.1406, 0.1414, 0.1369, 
c     &           0.1324, 0.1319, 0.1314, 0.1315, 0.1316/
c
cc/// data for extrapolated values
c
c      data xi0_60/1.13/, xi0_70/1.10/, xi0_80/1.08/, xi0_90/1.05/, 
c     &     xi0_100/1.02/
c      data xi1_60/0.11/, xi1_70/0.082/, xi1_80/0.072/, xi1_90/0.060/,
c     &     xi1_100/0.047/
c      data eta0_60/5.35/, eta0_70/5.65/, eta0_80/5.83/, eta0_90/6.00/,
c     &     eta0_100/6.15/
c      data eta1_60/0.120/, eta1_70/0.112/, eta1_80/0.108/,
c     &     eta1_90/0.102/, eta1_100/0.100/
c
cc/// if N=2, Z=1 (electron colliding with atomic hydrogen), use parameters
cc     chosen to mimic the exact electrostatic potential (allow some variance
cc     in z so that the value entered numerically may not equal 1.
c
c      if ((n.eq.2).and.(z.ge.0.95.and.z.le.1.05)) then
c         xi = 1.78
c         eta = 1.0
c         goto 100
c      end if
c
cc/// if N.le.54 use tabulated values of Garvey et al.
c     
c      if (n.le.54) then
c 	 zn = z-float(n)
c         xi = xi0(n) + zn*xi1(n)
c 	 eta = eta0(n) + zn*eta1(n)
c         goto 100
c      end if
c
cc/// if N.gt.54 use visually extrapolated values of Garvey et al.
cc     parameters (i.e. plot tabulated values and extrapolate). The
cc     user should note the uncertainty in using this approach.
c
c      if ((n.gt.54).and.(n.le.60)) then
c         xn=float(n)
c         xi0e=(xn-60.)/(-6.)*xi0(54) + (xn-54.)/6.*xi0_60
c         xi1e=(xn-60.)/(-6.)*xi1(54) + (xn-54.)/6.*xi1_60
c         eta0e=(xn-60.)/(-6.)*eta0(54) + (xn-54.)/6.*eta0_60
c         eta1e=(xn-60.)/(-6.)*eta1(54) + (xn-54.)/6.*eta1_60
c 	 zn = z-float(n)
c         xi = xi0e + zn*xi1e
c 	 eta = eta0e + zn*eta1e
c         goto 100
c      end if
c
c      if ((n.gt.60).and.(n.le.70)) then
c         xn=float(n)
c         xi0e=(xn-70.)/(-10.)*xi0_60 + (xn-60.)/10.*xi0_70
c         xi1e=(xn-70.)/(-10.)*xi1_60 + (xn-60.)/10.*xi1_70
c         eta0e=(xn-70.)/(-10.)*eta0_60 + (xn-60.)/10.*eta0_70
c         eta1e=(xn-70.)/(-10.)*eta1_60 + (xn-60.)/10.*eta1_70
c 	 zn = z-float(n)
c         xi = xi0e + zn*xi1e
c 	 eta = eta0e + zn*eta1e
c         goto 100
c      end if
c
c      if ((n.gt.70).and.(n.le.80)) then
c         xn=float(n)
c         xi0e=(xn-80.)/(-10.)*xi0_70 + (xn-70.)/10.*xi0_80
c         xi1e=(xn-80.)/(-10.)*xi1_70 + (xn-70.)/10.*xi1_80
c         eta0e=(xn-80.)/(-10.)*eta0_70 + (xn-70.)/10.*eta0_80
c         eta1e=(xn-80.)/(-10.)*eta1_70 + (xn-70.)/10.*eta1_80
c 	 zn = z-float(n)
c         xi = xi0e + zn*xi1e
c 	 eta = eta0e + zn*eta1e
c         goto 100
c      end if
c
c      if ((n.gt.80).and.(n.le.90)) then
c         xn=float(n)
c         xi0e=(xn-90.)/(-10.)*xi0_80 + (xn-80.)/10.*xi0_90
c         xi1e=(xn-90.)/(-10.)*xi1_80 + (xn-80.)/10.*xi1_90
c         eta0e=(xn-90.)/(-10.)*eta0_80 + (xn-80.)/10.*eta0_90
c         eta1e=(xn-90.)/(-10.)*eta1_80 + (xn-80.)/10.*eta1_90
c 	 zn = z-float(n)
c         xi = xi0e + zn*xi1e
c 	 eta = eta0e + zn*eta1e
c         goto 100
c      end if
c
c      if ((n.gt.90).and.(n.le.100)) then
c         xn=float(n)
c         xi0e=(xn-100.)/(-10.)*xi0_90 + (xn-90.)/10.*xi0_100
c         xi1e=(xn-100.)/(-10.)*xi1_90 + (xn-90.)/10.*xi1_100
c         eta0e=(xn-100.)/(-10.)*eta0_90 + (xn-90.)/10.*eta0_100
c         eta1e=(xn-100.)/(-10.)*eta1_90 + (xn-90.)/10.*eta1_100
c 	 zn = z-float(n)
c         xi = xi0e + zn*xi1e
c 	 eta = eta0e + zn*eta1e
c         goto 100
c      end if
c
cc/// if N.gt.100 write message
c
c      write(6,*) 'N chosen to be greater than 100, no parameters 
c     & available'
c
cc/// return here
c
c 100  return
c      end
c                                                                            ****
