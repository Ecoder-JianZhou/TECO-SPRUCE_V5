module vegetation
   !=================================================================================
   !  main subroutines  :  canopy,  respiration, plantgrowth 
   !             canopy => yrday,   xlayers,     Tsoil_simu
   !            xlayers => Radiso,  goudriaan,   agsean_day,  agsean_ngt
   !         agsean_day => photosyn
   !           photosyn => ciandA
   ! functions:
   !           sinbet, esat, VJtemp, fJQres, EnzK
   !=================================================================================
   use datatypes
   implicit none

   contains
   subroutine vegn_canopy(vegn, iforcing)
      implicit none
      type(vegn_tile_type), intent(inout) :: vegn
      type(forcing_data_type), intent(in) :: iforcing
      ! local variables
      integer :: doy, hour, ipft
      real    :: lat, radsol, fbeam, coszen, Radabv(2)

      doy    = iforcing%doy
      hour   = iforcing%hour + 1
      radsol = iforcing%PAR
      radsol = AMAX1(radsol,0.01)
      call yrday(doy, hour, lat, radsol, fbeam)  ! calculate beam fraction in incoming solar radiation
      coszen = sinbet(doy, hour, lat)

      do ipft = 1, vegn%npft
         call xlayers(vegn%allSp(ipft), iforcing, fbeam, coszen, Radabv)
      enddo 

   end subroutine vegn_canopy

   subroutine yrday(doy, hour, lat, radsol, fbeam)  
      ! Jian: This subrontine is used to calculate the fbeam.
      integer, intent(in) :: doy, hour
      real, intent(in)    :: lat, radsol
      real, intent(out)   :: fbeam
      ! parameters in this subrontine 
      real pidiv, slatx, sindec, cosdec
      real a, b, sinbet0, solext, tmprat, tmpR, tmpK, fdiff

      pidiv   = pi/180.0                ! divide the pi value
      slatx   = lat*pidiv               
      sindec  = -sin(23.4*pidiv)*cos(2.0*pi*(doy + 10.0)/365.0)
      cosdec  = sqrt(1.-sindec*sindec)
      a       = sin(slatx)*sindec
      b       = cos(slatx)*cosdec
      sinbet0 = a + b*cos(2*pi*(hour - 12.)/24.)
      solext  = 1370.0*(1.0 + 0.033*cos(2.0*pi*(doy - 10.)/365.0))*sinbet0
      tmprat  = radsol/solext                                              ! radsol: par in forcing
      tmpR    = 0.847 - 1.61*sinbet0 + 1.04*sinbet0*sinbet0
      tmpK    = (1.47 - tmpR)/1.66
      if (tmprat .le. 0.22) fdiff = 1.0
      if (tmprat .gt. 0.22 .and. tmprat .le. 0.35) then
         fdiff = 1.0 - 6.4*(tmprat - 0.22)*(tmprat - 0.22)
      end if
      if (tmprat .gt. 0.35 .and. tmprat .le. tmpK) then
         fdiff = 1.47 - 1.66*tmprat
      end if
      if (tmprat .ge. tmpK) then
         fdiff = tmpR
      end if
      fbeam = 1.0 - fdiff
      if (fbeam .lt. 0.0) fbeam = 0.0
      return
   end subroutine yrday

   subroutine xlayers(spec, iforcing, fbeam, coszen, radabv) 
      ! the multi-layered canopy model developed by
      ! Ray Leuning with the new radiative transfer scheme
      ! implemented by Y.P. Wang (from Sellers 1986)
      ! 12/Sept/96 (YPW) correction for mean surface temperature of sunlit
      ! and shaded leaves
      ! Tleaf,i=sum{Tleaf,i(n)*fslt*Gaussw(n)}/sum{fslt*Gaussw(n)}
      ! ----------------------------------------------------------------------
      type(spec_data_type), intent(inout) :: spec
      type(forcing_data_type), intent(in) :: iforcing
      real, intent(in) :: fbeam, coszen, Radabv(2)
      real :: Acan1, Acan2, Ecan1, Ecan2    ! outputs
      
      ! local vars
      real    :: Rnst1 , Rnst2 , Qcan1 , Qcan2  !net rad, sunlit, vis rad
      real    :: Rcan1 , Rcan2 , Hcan1 , Hcan2  !NIR rad !Sens heat
      real    :: Gbwc1 , Gbwc2 , Gswc1 , Gswc2  !Boundary layer conductance; Canopy conductance
      real    :: Tleaf1 , Tleaf2                !Leaf Temp
      real    :: xphi1, xphi2, funG 
      real    :: pi180, cozen15, cozen45, cozen75, xK15, xK45, xK75
      real    :: transd, extkn
      real    :: scatt(2)
      real    :: rhoc(3, 2)       !Goudriaan
      real    :: rhoch, rhoc15, rhoc45, rhoc75
      real    :: scalex, fslt, fshd
      real    :: RnStL(5), QcanL(5), RcanL(5), AcanL(5), EcanL(5), HcanL(5)
      real    :: layer1(5), layer2(5)
      real    :: FLAIT1, flait
      real    :: Gaussx(5), Gaussw(5), Gaussw_cum(5) 
      real    :: wind, raero, WILTPT, FILDCP
      real    :: TairK
      integer :: nw, ng
      ! transfer to other subroutine
      real    :: extkd, extkb, flai, kpr(3,2), reff(3,2)
      real    :: windUx, Vcmxx, eJmxx
      ! calculate in other subrontines
      real    :: emair, Qabs(3,2), Rnstar(2), grdn
      real    :: Tleaf(2)
      
      data Gaussx/0.0469101, 0.2307534, 0.5, 0.7692465, 0.9530899/        ! 5-point
      data Gaussw/0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635/
      data Gaussw_cum/0.11846, 0.35777, 0.64222, 0.88153, 1.0/

      flait  = spec%LAI                ! Jian: LAI relavant variable from vegetable to energy
      wind   = iforcing%WS
      if (wind .lt. 0.01) wind = 0.01   ! set windspeed to the minimum speed to avoid zero Gb
      ! soil water conditions
      WILTPT = st%wsmin/100.
      FILDCP = st%wsmax/100.
      ! reset the vairables
      Rnst1  = 0.0        !net rad, sunlit
      Rnst2  = 0.0        !net rad, shaded
      Qcan1  = 0.0        !vis rad
      Qcan2  = 0.0
      Rcan1  = 0.0        !NIR rad
      Rcan2  = 0.0
      Acan1  = 0.0        !CO2
      Acan2  = 0.0
      Ecan1  = 0.0        !Evap
      Ecan2  = 0.0
      Hcan1  = 0.0        !Sens heat
      Hcan2  = 0.0
      Gbwc1  = 0.0        !Boundary layer conductance
      Gbwc2  = 0.0
      Gswc1  = 0.0        !Canopy conductance
      Gswc2  = 0.0
      Tleaf1 = 0.0       !Leaf Temp
      Tleaf2 = 0.0  
      ! aerodynamic resistance
      raero = 50./wind
      ! Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13)
      xphi1 = 0.5 - 0.633*spec%xfang - 0.33*spec%xfang*spec%xfang
      xphi2 = 0.877*(1.0 - 2.0*xphi1)
      funG = xphi1 + xphi2*coszen      ! G-function: Projection of unit leaf area in direction of beam
      if (coszen .gt. 0) then          ! check if day or night
         extKb = funG/coszen           ! beam extinction coeff - black leaves
      else
         extKb = 100.
      end if
      ! Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
      ! Effective extinction coefficient for diffuse radiation Goudriaan & van Laar Eq 6.6)
      pi180   = 3.1416/180.
      cozen15 = cos(pi180*15)
      cozen45 = cos(pi180*45)
      cozen75 = cos(pi180*75)
      xK15    = xphi1/cozen15 + xphi2
      xK45    = xphi1/cozen45 + xphi2
      xK75    = xphi1/cozen75 + xphi2
      !--------------------------------------
      transd  = 0.308*exp(-xK15*flait) + 0.514*exp(-xK45*flait) +  &
                &  0.178*exp(-xK75*flait)
      extkd   = (-1./flait)*alog(transd)
      extkn   = extkd                        !N distribution coeff

      ! canopy reflection coefficients (Array indices: first -> 1=VIS,  2=NIR; second -> 1=beam, 2=diffuse)
      do nw = 1, 2  ! nw:1=VIS, 2=NIR
         scatt(nw)   = tauL(nw) + rhoL(nw)                                               ! scattering coeff
         if ((1.-scatt(nw)) < 0.0) scatt(nw) = 0.9999                                    ! Weng 10/31/2008
         kpr(nw, 1)  = extKb*sqrt(1.-scatt(nw))                                          ! modified k beam scattered (6.20)
         kpr(nw, 2)  = extkd*sqrt(1.-scatt(nw))                                          ! modified k diffuse (6.20)
         rhoch       = (1.-sqrt(1.-scatt(nw)))/(1.+sqrt(1.-scatt(nw)))                   ! canopy reflection black horizontal leaves (6.19)
         rhoc15      = 2.*xK15*rhoch/(xK15 + extkd)                                      ! canopy reflection (6.21) diffuse
         rhoc45      = 2.*xK45*rhoch/(xK45 + extkd)
         rhoc75      = 2.*xK75*rhoch/(xK75 + extkd)
         rhoc(nw, 2) = 0.308*rhoc15 + 0.514*rhoc45 + 0.178*rhoc75
         rhoc(nw, 1) = 2.*extKb/(extKb + extkd)*rhoch                                    ! canopy reflection (6.21) beam
         reff(nw, 1) = rhoc(nw, 1) + (rhoS(nw) - rhoc(nw, 1))*exp(-2.*kpr(nw, 1)*flait)  ! effective canopy-soil reflection coeff - beam (6.27)
         reff(nw, 2) = rhoc(nw, 2) + (rhoS(nw) - rhoc(nw, 2))*exp(-2.*kpr(nw, 2)*flait)  ! effective canopy-soil reflection coeff - diffuse (6.27)
      end do
      ! isothermal net radiation & radiation conductance at canopy top - needed to calc emair
      call Radiso(iforcing, fbeam, extkd, flait, flai, Qabs, emair, Rnstar, grdn)           ! Jian: some parameters not initialization.
      TairK = iforcing%Tair + 273.2
      Tleaf1 = 0.
      do ng = 1, 5
         flai = gaussx(ng)*flait
         ! radiation absorption for visible and near infra-red
         call goudriaan(fbeam, coszen, Radabv, kpr, reff, flai, scatt, Qabs)
         ! isothermal net radiation & radiation conductance at canopy top
         call Radiso(iforcing, fbeam, extkd, flait, flai, Qabs, emair, Rnstar, grdn)
         windUx = wind*exp(- st%extkU*flai)     ! windspeed at depth xi
         scalex = exp(-extkn*flai)              ! scale Vcmx0 & Jmax0
         Vcmxx  = spec%Vcmx0*scalex             ! Vcmx0 ---> Vcmax0
         eJmxx  = spec%eJmx0*scalex
         if (radabv(1) .ge. 10.0) then                          !check solar Radiation > 10 W/m2
            ! leaf stomata-photosynthesis-transpiration model - daytime
            call agsean_day()
         else
            call agsean_ngt()
         end if
         fslt      = exp(-extKb*flai)                        !fraction of sunlit leaves
         fshd      = 1.0 - fslt                                !fraction of shaded leaves
         Rnst1     = Rnst1 + fslt*Rnstar(1)*Gaussw(ng)*flait  !Isothermal net rad`
         Rnst2     = Rnst2 + fshd*Rnstar(2)*Gaussw(ng)*flait
         RnstL(ng) = Rnst1 + Rnst2

         Qcan1     = Qcan1 + fslt*Qabs(1, 1)*Gaussw(ng)*flait  !visible
         Qcan2     = Qcan2 + fshd*Qabs(1, 2)*Gaussw(ng)*flait
         QcanL(ng) = Qcan1 + Qcan2

         Rcan1     = Rcan1 + fslt*Qabs(2, 1)*Gaussw(ng)*flait  !NIR
         Rcan2     = Rcan2 + fshd*Qabs(2, 2)*Gaussw(ng)*flait
         RcanL(ng) = Rcan1 + Rcan2

         if (Aleaf(1) .lt. 0.0) Aleaf(1) = 0.0      !Weng 2/16/2006
         if (Aleaf(2) .lt. 0.0) Aleaf(2) = 0.0      !Weng 2/16/2006

         Acan1      = Acan1 + fslt*Aleaf(1)*Gaussw(ng)*flait*stom_n    !amphi/hypostomatous
         Acan2      = Acan2 + fshd*Aleaf(2)*Gaussw(ng)*flait*stom_n
         AcanL(ng)  = Acan1 + Acan2

         layer1(ng) = Aleaf(1)
         layer2(ng) = Aleaf(2)

         Ecan1      = Ecan1 + fslt*Eleaf(1)*Gaussw(ng)*flait
         Ecan2      = Ecan2 + fshd*Eleaf(2)*Gaussw(ng)*flait
         EcanL(ng)  = Ecan1 + Ecan2

         Hcan1     = Hcan1 + fslt*Hleaf(1)*Gaussw(ng)*flait
         Hcan2     = Hcan2 + fshd*Hleaf(2)*Gaussw(ng)*flait
         HcanL(ng) = Hcan1 + Hcan2

         Gbwc1     = Gbwc1 + fslt*gbleaf(1)*Gaussw(ng)*flait*stom_n
         Gbwc2     = Gbwc2 + fshd*gbleaf(2)*Gaussw(ng)*flait*stom_n

         Gswc1     = Gswc1 + fslt*gsleaf(1)*Gaussw(ng)*flait*stom_n
         Gswc2     = Gswc2 + fshd*gsleaf(2)*Gaussw(ng)*flait*stom_n

         Tleaf1    = Tleaf1 + fslt*Tleaf(1)*Gaussw(ng)*flait
         Tleaf2    = Tleaf2 + fshd*Tleaf(2)*Gaussw(ng)*flait
      end do  ! 5 layers

      FLAIT1 = (1.0 - exp(-extKb*flait))/extkb
      Tleaf1 = Tleaf1/FLAIT1
      Tleaf2 = Tleaf2/(flait - FLAIT1)
      ! Soil surface energy and water fluxes
      ! Radiation absorbed by soil
      Rsoilab1 = fbeam*(1.-reff(1, 1))*exp(-kpr(1, 1)*flait)        &
          &         + (1.-fbeam)*(1.-reff(1, 2))*exp(-kpr(1, 2)*flait)          !visible
      Rsoilab2 = fbeam*(1.-reff(2, 1))*exp(-kpr(2, 1)*flait)        &
          &         + (1.-fbeam)*(1.-reff(2, 2))*exp(-kpr(2, 2)*flait)          !NIR
      Rsoilab1 = Rsoilab1*Radabv(1)
      Rsoilab2 = Rsoilab2*Radabv(2)
      Tlk1     = Tleaf1 + 273.2
      Tlk2     = Tleaf2 + 273.2
      QLair    = emair*sigma*(TairK**4)*exp(-extkd*flait)
      QLleaf   = emleaf*sigma*(Tlk1**4)*exp(-extkb*flait)           &
                 &      + emleaf*sigma*(Tlk2**4)*(1.0 - exp(-extkb*flait))
      if(ISNAN(QLleaf)) then
         write(*,*)"QLleaf is NAN1111: ",QLleaf, emleaf, sigma, Tlk1, extKb, flait, Tlk2, extkd, Tleaf1, FLAIT1
         stop
      endif
      QLleaf   = QLleaf*(1.0 - exp(-extkd*flait))
      QLsoil   = emsoil*sigma*(TairK**4)
      Rsoilab3 = (QLair + QLleaf)*(1.0 - rhoS(3)) - QLsoil
      if(ISNAN(QLleaf)) then
         write(*,*)"QLleaf is NAN: ", emleaf, sigma, Tlk1, extKb, flait, Tlk2, extkd
         stop
      endif
      ! Net radiation absorbed by soil
      ! the old version of net long-wave radiation absorbed by soils
      ! (with isothermal assumption)
      ! Rsoil3=(sigma*TairK**4)*(emair-emleaf)*exp(-extkd*flait)         !Longwave
      ! Rsoilab3=(1-rhoS(3))*Rsoil3

      ! Total radiation absorbed by soil
      Rsoilabs = Rsoilab1 + Rsoilab2 + Rsoilab3

      ! thermodynamic parameters for air
      TairK  = Tair + 273.2
      rhocp  = cpair*Patm*AirMa/(Rconst*TairK)
      H2OLv  = H2oLv0 - 2.365e3*Tair
      slope  = (esat(Tair + 0.1) - esat(Tair))/0.1
      psyc   = Patm*cpair*AirMa/(H2OLv*H2OMw)
      Cmolar = Patm/(Rconst*TairK)
      fw1    = AMIN1(AMAX1((FILDCP - wcl(1))/(FILDCP - WILTPT), 0.05), 1.0)
      Rsoil  = 30.*exp(0.2/fw1)
      rLAI   = exp(flait)
      ! latent heat flux into air from soil
      ! Eleaf(ileaf)=1.0*
      ! &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    !2* Weng 0215
      ! &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
      Esoil  = (slope*(Rsoilabs - G) + rhocp*Dair/(raero + rLAI))/       &
               &      (slope + psyc*(rsoil/(raero + rLAI) + 1.))!*omega!*AMIN1(0.5, AMAX1(0.333, 0.333 + omega))
      ! sensible heat flux into air from soil
      Hsoil = Rsoilabs - Esoil - G
      
      return
   end subroutine xlayers

   ! functions 
   real function sinbet(doy, hour, lat)
      integer, intent(in) :: doy, hour
      real, intent(in) :: lat
      real rad, sinlat, coslat, sindec, cosdec, A, B
      ! sin(bet), bet = elevation angle of sun
      ! calculations according to Goudriaan & van Laar 1994 P30
      rad = pi/180.
      ! sine and cosine of latitude
      sinlat = sin(rad*lat)
      coslat = cos(rad*lat)
      ! sine of maximum declination
      sindec = -sin(23.45*rad)*cos(2.0*pi*(doy + 10.0)/365.0)
      cosdec = sqrt(1.-sindec*sindec)
      ! terms A & B in Eq 3.3
      A = sinlat*sindec
      B = coslat*cosdec
      sinbet = A + B*cos(pi*(hour - 12.)/12.)
      return
   end function sinbet

end module vegetation