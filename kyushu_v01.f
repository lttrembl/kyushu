c--------------------------------------------------------------------c
c                                                                    c
c kyushu1.f                                                          c
c                                                                    c
c written by Steve Desch                                             c
c April 9, 2017                                                      c
c code updates August 24, 2017                                       c
c comments updated April 3, 2019                                     c
c                                                                    c
c This code implements the algorithm of Hachisu (1986a),             c
c Astrophysical Journal Supplement 62, 461-499.                      c
c                                                                    c
c It calculates the density structure of a fluid in hydrostatic      c
c  equilibrium in rapid uniform rotation. It is applicable to        c
c  uniform-density Jacobi ellipsoids, or situations where the        c
c  density is a known function of pressure.                          c
c In this version, density is set to rho_ice for pressures below a   c
c  threshold value, rho_core for higher pressures. This is           c
c  adequate to model the dwarf planet Haumea.                        c
c                                                                    c
c The code accepts the mass, rotation period, and long and short     c
c  semi-axes normal to the rotation axis. It then finds the          c
c  hydrostatic equilibrium solution that matches the total mass      c
c  and rotation period.                                              c
c                                                                    c
c Inputs:                                                            c
c  a    = the longer semi-axis perpendicular to the rotation axis    c
c  b    = the shorter semi-axis perpendicular to the rotation axis   c
c  rho_ice = the density of the icy mantle                           c
c  P0   = the observed rotation period                               c
c  M0   = the observed total mass                                    c
c  Nr   = number of radial points in spherical grid                  c
c  Nphi = number of azimuthal angles in spherical grid               c
c  Nmu  = number of polar angles in spherical grid                   c
c                                                                    c
c Outputs:                                                           c
c  c    = semi-axis along rotation axis                              c
c  rhoc = the density of the core material                           c
c  Pcmb = pressure of the rocky core-icy mantle boundary             c
c  rho  = density at all grid points                                 c
c  P    = calculated rotation period (should match P0)               c
c  M    = calculated mass (should match M0)                          c
c                                                                    c
c--------------------------------------------------------------------c

       implicit none
       integer iphi,Nphi,nphimax,imu,Nmu,nmumax,ir,Nr,nrmax
       integer iphiA,imuA,irA,iphiB,imuB,irB
       parameter(nphimax=129,nmumax=129,nrmax=999)
       double precision phi(nphimax),mu(nmumax),r(nrmax),rmax
       double precision rho(nphimax,nmumax,nrmax),rhorock,rhoice
       double precision a,b,c,ac,bc,cc,xx
       double precision rho1,P1,rho2,P2,dPdrho
       integer ip,ic
       double precision M,P,M0,P0
       double precision rhoc,Pcmb
       double precision xlo,xhi,xmd,flo,fhi,fmd
       double precision pi
       parameter(pi=3.141592653589793d0)
       double precision km,hr,MPa
       parameter(km=1.0d+03,hr=3600.,MPa=1.0d+06)

       logical verbose
       verbose = .false.

c--------------------------------------------------------------------c
c Fixed parameters                                                   c
c Note: we impose a and b, but if the grid we create doesn't have    c
c  points lying exactly at r=a or r=b, it will choose the closest    c
c  radial points, assuming they are within 2 km of a and b. So the   c
c  actual axes a and b will be approximations to these imposed       c
c  values.                                                           c
c--------------------------------------------------------------------c

       a = 1000.*km
       b =  852.*km
       P0 = 3.9155*hr
       M0 = 4.006d+21

c--------------------------------------------------------------------c
c Create grid                                                        c
c--------------------------------------------------------------------c

       Nphi = 33
       Nmu  = 33
       Nr   = 391.
       rmax = 1300.*km

       call grid(nphi,phi,nmu,mu,nr,rmax,r)

        irA = -99
        irB = -99
       do ir=1,Nr
        if (dabs(r(ir)-a).lt.(2.*km)) then
         irA = ir
        end if
        if (dabs(r(ir)-b).lt.(2.*km)) then
         irB = ir
        end if
       end do
       if ((irA.lt.0).or.(irB.lt.0)) then
         write(*,*) ' Could not assign points A '
         write(*,*) ' irA = ',irA
         write(*,*) ' irB = ',irB
        stop
       end if

       iphiA = 1
       imuA  = 1
       iphiB = Nphi
       imuB  = 1

       write(*,*) ' Input a = ',a/km,' km'
       write(*,*) ' Input b = ',b/km,' km'
       write(*,*) ' b/a = ',b/a
       write(*,*) ' Using a = ',r(irA)/km,' km'
       write(*,*) ' Using b = ',r(irB)/km,' km'
       write(*,*) ' b/a = ',r(irB)/r(irA)
       if (verbose) then
        write(*,*) ' Point A:'
        write(*,*) ' x = ',r(irA)
     &               *dsqrt(1.0d0-mu(imuA)*mu(imuA))*cos(phi(iphiA))/km
        write(*,*) ' y = ',r(irA)
     &               *dsqrt(1.0d0-mu(imuA)*mu(imuA))*sin(phi(iphiA))/km
        write(*,*) ' z = ',r(irA)*mu(imuA)/km
        write(*,*) ' Point B:'
        write(*,*) ' x = ',r(irB)
     &               *dsqrt(1.0d0-mu(imuB)*mu(imuB))*cos(phi(iphiB))/km
        write(*,*) ' y = ',r(irB)
     &               *dsqrt(1.0d0-mu(imuB)*mu(imuB))*sin(phi(iphiB))/km
        write(*,*) ' z = ',r(irB)*mu(imuB)/km
        write(*,*) ' '
       end if

c Initialize the density, just this first time
       rhorock = 2700.
       rhoice  =  921.
       c    = b
       xx   = 1.0d0
       ac   = a *xx
       bc   = b *xx
       cc   = c *xx
       call gridfill(Nphi,phi,Nmu,mu,Nr,r,
     &                a,b,c,ac,bc,cc,rhorock,rhoice,rho)

c--------------------------------------------------------------------c
c Loop over rhoc with fixed Pcmb=0 until Mass is right               c
c--------------------------------------------------------------------c

c      write(*,*) ' Doing homogeneous (Pcmb=20 MPa) case '
c      Pcmb = 20.*MPa
       write(*,*) ' Doing homogeneous (Pcmb=0 MPa) case '
       Pcmb = 0.*MPa
       rho1 = rhorock

       do ic=1,3

       write(*,*) ' Calling hachisu with '
       write(*,*) ' a (km) = ',a/km
       write(*,*) ' b (km) = ',b/km
       write(*,*) ' Pcmb (MPa) = ',Pcmb/MPa
       write(*,*) ' rho1       = ',rho1

       call hachisu(Nphi,Nmu,Nr,nphimax,nmumax,nrmax,phi,mu,r,
     &               iphiA,imuA,irA,iphiB,imuB,irB,
     &               a,b,c,rho1,Pcmb,M,P,rho)

       write(*,*) ' hachisu returns '
       write(*,*) ' M (kg) = ',M
       write(*,*) ' P (hr) = ',P/hr

       rho1 = rho1 *(M0/M)

       end do

       if (P.lt.P0) then
        write(*,*) ' P < P0 in homogeneous case '
        stop
       end if

       write(*,*) ' Doing Pcmb =  50 MPa case '
       Pcmb =  50.*MPa
       rho1 = rhorock

       do ic=1,10

       write(*,*) ' Calling hachisu with '
       write(*,*) ' a (km) = ',a/km
       write(*,*) ' b (km) = ',b/km
       write(*,*) ' Pcmb (MPa) = ',Pcmb/MPa
       write(*,*) ' rho1       = ',rho1

       call hachisu(Nphi,Nmu,Nr,nphimax,nmumax,nrmax,phi,mu,r,
     &               iphiA,imuA,irA,iphiB,imuB,irB,
     &               a,b,c,rho1,Pcmb,M,P,rho)

       write(*,*) ' hachisu returns '
       write(*,*) ' M (kg) = ',M
       write(*,*) ' P (hr) = ',P/hr

       rho1 = rho1 *(M0/M)**0.5

       end do

       if (P.gt.P0) then
        write(*,*) ' P > P0 in very differentiated case '
        stop
       end if

       write(*,*) ' OK, we have bracketed the solution'
       xlo = 0.
       xhi = Pcmb

       do ip=1,10

       xmd = (xlo+xhi)/2.
       Pcmb = xmd
       rho1 = rhorock

       write(*,*) ' Doing Pcmb = ',xmd/MPa,' MPa case '

       do ic=1,10

       write(*,*) ' Calling hachisu with '
       write(*,*) ' a (km) = ',a/km
       write(*,*) ' b (km) = ',b/km
       write(*,*) ' Pcmb (MPa) = ',Pcmb/MPa
       write(*,*) ' rho1       = ',rho1

       call hachisu(Nphi,Nmu,Nr,nphimax,nmumax,nrmax,phi,mu,r,
     &               iphiA,imuA,irA,iphiB,imuB,irB,
     &               a,b,c,rho1,Pcmb,M,P,rho)

       write(*,*) ' hachisu returns '
       write(*,*) ' M (kg) = ',M
       write(*,*) ' P (hr) = ',P/hr

       rho1 = rho1 *(M0/M)**0.5

       end do

       if (P.lt.P0) then
        xhi = xmd
       else
        xlo = xmd
       end if

       end do


c--------------------------------------------------------------------c
c WRITE OUTPUT                                                       c
c--------------------------------------------------------------------c

         write(*,*) ' '
         write(*,*) 'Solution found:'
         write(*,*) ' Mass = ',M,' kg'
         write(*,*) ' a = ',r(irA)/km
         write(*,*) ' b = ',r(irB)/km
         write(*,*) ' c = ',c/km
         write(*,*) ' rhoavg = ',M*3./(4.*pi*a*b*c),' kg m-3'
         write(*,*) ' Period = ',P/hr,' hr'

         write(*,*) 'Density from center out, along each axis:'
         write(*,*) 'r        a axis        b axis        c axis '
         do ir=1,Nr
          write(*,*) r(ir)/km,
     &      rho(1,1,ir),rho(Nphi,1,ir),rho(1,Nmu,ir)
         end do

       end


c====================================================================c
c                                                                    c
c SUBROUTINE HACHISU                                                 c
c                                                                    c
c INPUTS:                                                            c
c  a    = longest semi-major axis                                    c
c  b    = second longest semi-major axis                             c
c  rhoc = density of core material                                   c
c  Pcmb = pressure at core-mantle boundary                           c
c                                                                    c
c OUTPUTS:                                                           c
c  M    = mass of equilibrium configuration                          c
c  P    = rotation period of equilibrium configuration               c
c                                                                    c
c====================================================================c

       subroutine hachisu(Nphi,Nmu,Nr,nphimax,nmumax,nrmax,phi,mu,r,
     &               iphiA,imuA,irA,iphiB,imuB,irB,
     &               a,b,c,rhoc,Pcmb,M,P,rho)

       implicit none
       integer iphi,Nphi,nphimax,imu,Nmu,nmumax,ir,Nr,nrmax
       integer is,it,iu
       double precision phi(nphimax),mu(nmumax),r(nrmax)
       double precision a,b,c,ac,bc,cc,rhoc,rhoi,rmax
       double precision rho(nphimax,nmumax,nrmax)
       double precision Qone(nmumax,nrmax),Qtwo(nrmax)
       double precision dphi,dmu,dr
       double precision M,xx
       integer il,nl,nlmax,im,nm,nmmax
       parameter(nlmax=32,nmmax=32)
       double precision Phig(nphimax,nmumax,nrmax)
       double precision fl(nrmax),emfac,y,Plm(nmumax)
       double precision f1,f2
       double precision cosmphi(nphimax)
       double precision Done(nmumax,nrmax),Dtwo(nrmax),Dthree
       double precision Psi(nphimax,nmumax,nrmax),Omega2,Const
       integer iphiA,imuA,irA,iphiB,imuB,irB
       double precision H(nphimax,nmumax,nrmax)
       double precision Hold(nphimax,nmumax,nrmax)
       double precision H0,rho0,xlo,xhi,xmd,Hlo,Hhi,Hmd,Htarget
       double precision rhonew(nphimax,nmumax,nrmax)
       double precision Constold,Omega2old,eps1,eps2,eps3
       integer iter,itermax,ic
       double precision P,Rc,Pcmb,Hcmb,rhoavg ,rhorock,rhoice
       double precision pi
       parameter(pi=3.141592653589793d0)
       double precision km,MPa
       parameter(km=1.0d+03,MPa=1.0d+06)
       double precision Gconst
       parameter(Gconst=6.672d-11)

       logical verbose
       verbose = .false.

       rhoice  = 921.
       rhorock = rhoc
       Hcmb = Pcmb / rhoice

       Constold = 0.0d0
       Omega2old = 0.0d0

c Set up Psi
       do iphi=1,Nphi
        do imu=1,Nmu
         do ir=1,Nr
          Psi(iphi,imu,ir) = -0.5*r(ir)*r(ir)*(1.0d0-mu(imu)*mu(imu))
          Hold(iphi,imu,ir) = 0.
         end do
        end do
       end do

c--------------------------------------------------------------------c
c BEGIN ITERATIONS                                                   c
c--------------------------------------------------------------------c

       itermax = 25
       iter = 0
610    continue
       iter = iter + 1
        if (verbose) then
        write(*,*) ' ITERATION #',iter
        end if

c STEP #0
c  INTEGRATION TO GET TOTAL MASS
c  This follows Equation 39 of H86b

       dphi = phi(2)-phi(1)
       dmu  = mu (2)-mu (1)
       dr   = r  (2)-r  (1)
       do imu=1,Nmu
        do ir=1,Nr
         Qone(imu,ir) = 0.
         do iphi=1,Nphi-2,2
          Qone(imu,ir) = Qone(imu,ir)
     &      +(2.0)/(3.0)*(2.*dphi)*
     &    (rho(iphi,imu,ir)+4.*rho(iphi+1,imu,ir)+rho(iphi+2,imu,ir))
         end do
        end do
       end do

       do ir=1,nr
        Qtwo(ir) = 0.
        do imu=1,Nmu-2,2
         Qtwo(ir) = Qtwo(ir)
     &   +(1.0)/(3.0)*(2.*dmu)*
     &   (Qone(imu,ir)+4.*Qone(imu+1,ir)+Qone(imu+2,ir))
        end do
       end do

       M = 0.0d0
       do ir=1,Nr-2,2
        M = M
     &    + (2.*dr)/6.0*
     &    (r(ir)*r(ir)*Qtwo(ir)+4.*r(ir+1)*r(ir+1)*Qtwo(ir+1)
     &    +r(ir+2)*r(ir+2)*Qtwo(ir+2))
       end do


c STEP #1
c  INTEGRATION TO GET GRAVITATIONAL POTENTIAL PHIG
c  Using integral in equation 2, but really equations 33-36

c LMAX
      nl = 16

c Zero out potential
      do iphi=1,Nphi
       do imu=1,Nmu
        do ir=1,Nr
         Phig(iphi,imu,ir) = 0.0
        end do
       end do
      end do

c Loop through all even l, even m

      do il=0,nl,2
       do im=0,il,2

c create Plm,cos(mphi) arrays
        do iphi=1,Nphi
         cosmphi(iphi) = cos(dble(im)*phi(iphi))
        end do
        do imu=1,Nmu
         call plgndr(il,im,mu(imu),y)
         Plm(imu) = y
        end do
c find em,factorials
        call fctrl(il-im,f1)
        call fctrl(il+im,f2)
        emfac = (2.0d0)*(f1)/(f2)
        if (im.eq.0) then
         emfac=1.0d0
        end if

c Begin loop over ir
        do ir=1,nr
c create fl array
          do is=1,Nr
           if (r(is).lt.r(ir)) then
            fl(is) = r(is)*(r(is)/r(ir))**dble(il+1)
           end if
           if (r(is).eq.r(ir)) then
            fl(is) = r(is)
           end if
           if (r(is).gt.r(ir)) then
            fl(is) = r(is)*(r(ir)/r(is))**dble(il)
           end if
          end do

c create Done array
          do it=1,Nmu
           do is=1,Nr
            Done(it,is) = 0.0
            do iu=1,Nphi-2,2
            Done(it,is) = Done(it,is)
     &                    +(2.0)/(3.0)*(2.*dphi)*
     &        (    cosmphi(iu  )*rho(iu,it,is  )
     &         +4.*cosmphi(iu+1)*rho(iu+1,it,is)
     &         +   cosmphi(iu+2)*rho(iu+2,it,is) )
            end do
           end do
          end do

c create Dtwo array
          do is=1,Nr
           Dtwo(is) = 0.0
           do it=1,Nmu-2,2
            Dtwo(is) = Dtwo(is)
     &         +2.*dmu/3.*
     &         (    Plm(it  )*Done(it  ,is)
     &          +4.*Plm(it+1)*Done(it+1,is)
     &          +   Plm(it+2)*Done(it+2,is) )
           end do
          end do

c create Dthree
          Dthree = 0.0
          do is=1,Nr-2,2
           Dthree = Dthree
     &        + 2.*dr/6.*
     &         (     fl(is  )*Dtwo(is  )
     &           +4.*fl(is+1)*Dtwo(is+1)
     &           +   fl(is+2)*Dtwo(is+2) )
          end do

c Augment all Phi(iphi,imu,ir) for this il,im combination

        do iphi=1,Nphi
         do imu=1,Nmu
          Phig(iphi,imu,ir) = Phig(iphi,imu,ir)
     &     -Gconst*emfac*Dthree*Plm(imu)*cosmphi(iphi)
         end do
        end do

        end do

       end do
      end do

c STEP #2
c  DETERMINATION OF ENTHALPY
c  Using equation 8
c  CALCULATION OF omega2n and Cn from Equations 16 and 17
c  have to use the values of PHIg at two points, "A" and "B"
c  H=0 on surface, at two points A and B

       Omega2 = -(Phig(iphiA,imuA,irA) - Phig(iphiB,imuB,irB)) /
     &           ( Psi(iphiA,imuA,irA) - Psi (iphiB,imuB,irB))
       Const  = Phig(iphiA,imuA,irA) + Omega2*Psi(iphiA,imuA,irA)
       if (verbose) then
        write(*,*) '  Omega2 = ',Omega2
        write(*,*) '  Const = ',Const
        write(*,*) ' Period = ',2.*pi/sqrt(Omega2)/3600.,' hr'
       end if

       do iphi=1,Nphi
        do imu=1,Nmu
         do ir=1,Nr
          H(iphi,imu,ir) = Const - Phig(iphi,imu,ir)
     &                      -Omega2*Psi(iphi,imu,ir)
         end do
        end do
       end do
c      write(*,*) ' H(A) = ',H(iphiA,imuA,irA)
c      write(*,*) ' H(B) = ',H(iphiB,imuB,irB)

c STEP #3
c  DETERMINE NEW DENSITY DISTRIBUTION
c This is where the very simplified equation of state is applied.

        do iphi=1,Nphi
         do imu=1,Nmu
          do ir=1,Nr
           Htarget = H(iphi,imu,ir)
           if (Htarget.lt.0.) then
            rho(iphi,imu,ir) = 0.
           else
            if (Htarget.lt.Hcmb) then
             rho(iphi,imu,ir) = rhoice
            else
             rho(iphi,imu,ir) = rhorock
            end if
           end if

          end do
         end do
        end do

c CHECK FOR CONVERGENCE

        eps1 = 0.
        do iphi=1,Nphi
         do imu=1,Nmu
          do ir=1,Nr
           if (H(iphi,imu,ir).gt.(1.)) then
            eps2 = dabs(H(iphi,imu,ir)-Hold(iphi,imu,ir))
     &                / H(iphi,imu,ir)
            eps1 = max(eps1,eps2)
            end if
            Hold(iphi,imu,ir) = H(iphi,imu,ir)
          end do
         end do
        end do

        eps2 = dabs(Omega2-Omega2old)/dabs(Omega2)
        Constold = Const
        Omega2old = Omega2
        if (verbose) then
         write(*,*) ' Change in H: ',eps1
         write(*,*) ' Change in Omega2: ',eps2
        end if

c END ITERATIONS

        if (eps2.le.(1.0d-04)) then
         write(*,*) ' Converged after ',iter,' iterations'
        else
         if (iter.lt.itermax) then
          goto 610
         else
           write(*,*) 'Exceeded maximum # iterations!'
         end if
        end if

c GENERATE OUTPUT

         P = 2.*pi/dsqrt(Omega2)
         iphi=1
         imu = Nmu
         do ir=1,nr-1
          if ( (H(iphi,imu,ir  ).gt.0.).and.
     &         (H(iphi,imu,ir+1).le.0.)      ) then
           ic = ir
          end if
         end do
         c = (H(iphi,imu,ic+1)*r(ic)-H(iphi,imu,ic)*r(ic+1)) /
     &       (H(iphi,imu,ic+1)      -H(iphi,imu,ic)        )

      return
      end

c====================================================================c
c                                                                    c
c SUBROUTINE GRID                                                    c
c                                                                    c
c phi from Equation 28. n=2 because even though not 2-body system,   c
c  has twofold symmetry.  mu from Equation 29.  r from Equation 30.  c
c note rmax should be slightly greater than 1; they use 16/15 times  c
c  long axis.  Note they use dimensionless values; I don't.          c
c                                                                    c
c====================================================================c

       subroutine grid(nphi,phi,nmu,mu,nr,rmax,r)

       implicit none
       integer iphi,nphi,nphimax,imu,nmu,nmumax,ir,nr,nrmax
       parameter(nphimax=129,nmumax=129,nrmax=999)
       double precision phi(nphimax),mu(nmumax),r(nrmax),rmax
       double precision pi
       parameter(pi=3.141592653589793d0)

       do iphi=1,nphi
        phi(iphi) = pi*dble(iphi-1)/dble(nphi-1)/(2.0d0)
       end do
       do imu=1,nmu
        mu(imu) = dble(imu-1)/dble(nmu-1)
       end do
       do ir=1,nr
        r(ir)=rmax*dble(ir-1)/dble(nr-1)
       end do

       return
       end

c====================================================================c
c                                                                    c
c SUBROUTINE GRIDFILL                                                c
c                                                                    c
c We compute points at the intersections of grid lines (see comment  c
c  before Equation 31). Those that are inside triaxial ellipsoid     c
c  defined by ac,bc,cc core radii have rho = rhoc; those outside     c
c  have rho=rhoi.  Jacobi test requires setting rhoc=rhoi.           c
c                                                                    c
c====================================================================c

       subroutine gridfill(nphi,phi,nmu,mu,nr,r,
     &                a,b,c,ac,bc,cc,rhoc,rhoi,rho)

       implicit none
       integer iphi,nphi,nphimax,imu,nmu,nmumax,ir,nr,nrmax
       parameter(nphimax=129,nmumax=129,nrmax=129)
       double precision phi(nphimax),mu(nmumax),r(nrmax)
       double precision a,b,c,ac,bc,cc,rhoc,rhoi,f
       double precision rho(nphimax,nmumax,nrmax)
       double precision x,y,z

       do iphi=1,nphi
        do imu=1,nmu
         do ir=1,nr
          rho(iphi,imu,ir) = 0.0d0
          x = r(ir)*dsqrt(1.0d0-mu(imu)*mu(imu))*cos(phi(iphi))
          y = r(ir)*dsqrt(1.0d0-mu(imu)*mu(imu))*sin(phi(iphi))
          z = r(ir)*mu(imu)
          f = (x/a)*(x/a) + (y/b)*(y/b) + (z/c)*(z/c)
          if (f.le.(1.0d0)) then
           rho(iphi,imu,ir) = rhoi
          end if
          f = (x/ac)*(x/ac) + (y/bc)*(y/bc) + (z/cc)*(z/cc)
          if (f.le.(1.0d0)) then
           rho(iphi,imu,ir) = rhoc
          end if
         end do
        end do
       end do

       return
       end

c====================================================================c
c                                                                    c
c SUBROUTINE PLGNDR                                                  c
c                                                                    c
c====================================================================c

      subroutine plgndr(l,m,x,y)

      implicit none
      integer l,m
      double precision x,y
      integer i,ll
      double precision fact,pll,pmm,pmmp1,somx2

      if ((m.lt.0).or.(m.gt.l).or.(dabs(x).gt.1.0d0)) then
       write(*,*) 'bad arguments in plgndr'
       stop
      end if

      pmm = 1.0d0
      if (m.gt.0) then
       somx2 = dsqrt((1.0d0-x)*(1.0d0+x))
       fact = 1.0d0
       do i=1,m
        pmm=-pmm*fact*somx2
        fact=fact+2.0d0
       end do
      end if
      if (l.eq.m) then
       y = pmm
      else
       pmmp1=x*(2.0d0*dble(m)+1.0d0)*pmm
       if (l.eq.m+1) then
        y = pmmp1
       else
        do ll=m+2,l
         pll= ( x*(2.0d0*dble(ll)-1.0d0)*pmmp1 - dble(ll+m-1)*pmm )
     &          / dble(ll-m)
         pmm= pmmp1
         pmmp1=pll
        end do
        y = pll
       end if
      end if

      return
      end

c====================================================================c
c                                                                    c
c SUBROUTINE FCTRL                                                   c
c                                                                    c
c====================================================================c

      subroutine fctrl(n,f)

      implicit none
      integer n,i
      double precision f

      if (n.lt.0) then
       write(*,*) 'negative input in fctrl'
       stop
      end if

      f = 1.0d0

      if (n.gt.0) then
       do i=1,n
        f = f*dble(i)
       end do
      end if

      return
      end
