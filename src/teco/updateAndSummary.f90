! This module is used to summary the values of hourly, daily, monthly and yearly
module update_and_summary
    use datatypes
    implicit NONE
    real convert_g2kg, convert_h2s

    contains

    subroutine updateOutVars(itime, ntime, vegn, outvars)
        implicit none
        integer, intent(in) :: itime, ntime
        type(outvars_data_type), intent(inout) :: outVars
        type(vegn_tile_type), intent(in) :: vegn
        integer :: ipft, npft
        ! integer iTotHourly
        convert_g2kg = 0.001
        convert_h2s  = 1/3600.
        if (allocated(outVars_h%allSpec)) then
            npft = size(outVars_h%allSpec)
            do ipft = 1, npft
                ! carbon fluxes (Kg C m-2 s-1)
                outVars_h%allSpec(ipft)%gpp(itime)      = outVars_h%allSpec(ipft)%gpp(itime)     + &
                                                          vegn%allSp(ipft)%gpp*convert_g2kg*convert_h2s/ntime
                ! outVars_h%allSpec(ipft)%nee           = vegn%allSp(ipft)%NEE
                outVars_h%allSpec(ipft)%npp(itime)      = outVars_h%allSpec(ipft)%npp(itime)     + &
                                                          vegn%allSp(ipft)%npp*convert_g2kg*convert_h2s/ntime
                outVars_h%allSpec(ipft)%nppLeaf(itime)  = outVars_h%allSpec(ipft)%nppLeaf(itime) + &
                                                          vegn%allSp(ipft)%NPP_L*convert_g2kg*convert_h2s/ntime
                outVars_h%allSpec(ipft)%nppWood(itime)  = outVars_h%allSpec(ipft)%nppWood(itime) + &
                                                          vegn%allSp(ipft)%NPP_W*convert_g2kg*convert_h2s/ntime
                outVars_h%allSpec(ipft)%nppStem(itime)  = outVars_h%allSpec(ipft)%nppStem(itime) + &
                                                          vegn%allSp(ipft)%NPP_W*convert_g2kg*convert_h2s/ntime
                outVars_h%allSpec(ipft)%nppRoot(itime)  = outVars_h%allSpec(ipft)%nppRoot(itime) + &
                                                          vegn%allSp(ipft)%NPP_R*convert_g2kg*convert_h2s/ntime
                ! outVars_h%allSpec(ipft)%nppOther      = vegn%allSp(ipft)%nppOther    ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
                outVars_h%allSpec(ipft)%ra(itime)       = outVars_h%allSpec(ipft)%ra(itime)      + &
                                                          vegn%allSp(ipft)%Rauto*convert_g2kg*convert_h2s/ntime
                outVars_h%allSpec(ipft)%raLeaf(itime)   = outVars_h%allSpec(ipft)%raLeaf(itime)  + &
                                                          vegn%allSp(ipft)%RmLeaf*convert_g2kg*convert_h2s/ntime
                outVars_h%allSpec(ipft)%raStem(itime)   = outVars_h%allSpec(ipft)%raStem(itime)  + &
                                                          vegn%allSp(ipft)%RmStem*convert_g2kg*convert_h2s/ntime
                outVars_h%allSpec(ipft)%raRoot(itime)   = outVars_h%allSpec(ipft)%raRoot(itime)  + &
                                                          vegn%allSp(ipft)%RmRoot*convert_g2kg*convert_h2s/ntime
                ! outVars_h%allSpec(ipft)%raOther  = vegn%allSp(ipft)%
                outVars_h%allSpec(ipft)%rMaint(itime)   = outVars_h%allSpec(ipft)%rMaint(itime)  + &
                                                          vegn%allSp(ipft)%Rmain*convert_g2kg*convert_h2s/ntime
                outVars_h%allSpec(ipft)%rGrowth(itime)  = outVars_h%allSpec(ipft)%rGrowth(itime) + &
                                                          vegn%allSp(ipft)%Rgrowth*convert_g2kg*convert_h2s/ntime
                ! outVars_h%allSpec(ipft)%nbp      = vegn%allSp(ipft)%
                ! Carbon Pools  (KgC m-2)
                outVars_h%allSpec(ipft)%cLeaf(itime)    = outVars_h%allSpec(ipft)%cLeaf(itime)  + &
                                                          vegn%allSp(ipft)%QC(1)*convert_g2kg/ntime
                outVars_h%allSpec(ipft)%cStem(itime)    = outVars_h%allSpec(ipft)%cStem(itime)  + &
                                                          vegn%allSp(ipft)%QC(2)*convert_g2kg/ntime
                outVars_h%allSpec(ipft)%cRoot(itime)    = outVars_h%allSpec(ipft)%cRoot(itime)  + &
                                                          vegn%allSp(ipft)%QC(3)*convert_g2kg/ntime
                ! Nitrogen pools (kgN m-2)
                outVars_h%allSpec(ipft)%nLeaf(itime)    = outVars_h%allSpec(ipft)%nLeaf(itime)  + &
                                                          vegn%allSp(ipft)%QN(1)*convert_g2kg/ntime
                outVars_h%allSpec(ipft)%nStem(itime)    = outVars_h%allSpec(ipft)%nStem(itime)  + &
                                                          vegn%allSp(ipft)%QN(2)*convert_g2kg/ntime
                outVars_h%allSpec(ipft)%nRoot(itime)    = outVars_h%allSpec(ipft)%nRoot(itime)  + &
                                                          vegn%allSp(ipft)%QN(3)*convert_g2kg/ntime
                ! water fluxes (kg m-2 s-1)
                outVars_h%allSpec(ipft)%tran(itime)     = outVars_h%allSpec(ipft)%tran(itime)   + &
                                                          vegn%allSp(ipft)%transp*convert_g2kg*convert_h2s/ntime
                ! other
                outVars_h%allSpec(ipft)%lai(itime)      = outVars_h%allSpec(ipft)%lai(itime)    + &
                                                          vegn%allSp(ipft)%LAI/ntime
            enddo
        endif
        ! carbon fluxes (KgC m-2 s-1) Jian: TECO unit is gC m-2 h-1
        outVars_h%gpp(itime)             = outVars_h%gpp(itime)      + vegn%gpp*convert_g2kg*convert_h2s/ntime
        outVars_h%npp(itime)             = outVars_h%npp(itime)      + vegn%npp*convert_g2kg*convert_h2s/ntime
        outVars_h%nppLeaf(itime)         = outVars_h%nppLeaf(itime)  + vegn%NPP_L*convert_g2kg*convert_h2s/ntime
        outVars_h%nppWood(itime)         = outVars_h%nppWood(itime)  + vegn%NPP_W*convert_g2kg*convert_h2s/ntime  
        outVars_h%nppStem(itime)         = outVars_h%nppStem(itime)  + vegn%NPP_W*convert_g2kg*convert_h2s/ntime 
        outVars_h%nppRoot(itime)         = outVars_h%nppRoot(itime)  + vegn%NPP_R*convert_g2kg*convert_h2s/ntime
        outVars_h%nppOther(itime)        = outVars_h%nppOther(itime) + 0*convert_g2kg*convert_h2s/ntime 
        outVars_h%ra(itime)              = outVars_h%ra(itime)       + vegn%Rauto*convert_g2kg*convert_h2s/ntime
        outVars_h%raLeaf(itime)          = outVars_h%raLeaf(itime)   + vegn%Rmleaf*convert_g2kg*convert_h2s/ntime
        outVars_h%raStem(itime)          = outVars_h%raStem(itime)   + vegn%Rmstem*convert_g2kg*convert_h2s/ntime
        outVars_h%raRoot(itime)          = outVars_h%raRoot(itime)   + vegn%Rmroot*convert_g2kg*convert_h2s/ntime
        outVars_h%raOther(itime)         = outVars_h%raOther(itime)  + st%Rnitrogen *convert_g2kg*convert_h2s/ntime
        outVars_h%rMaint(itime)          = outVars_h%rMaint(itime)   + vegn%Rmain *convert_g2kg*convert_h2s/ntime
        outVars_h%rGrowth(itime)         = outVars_h%rGrowth(itime)  + vegn%Rgrowth *convert_g2kg*convert_h2s/ntime 
        outVars_h%rh(itime)              = outVars_h%rh(itime)       + st%Rhetero *convert_g2kg*convert_h2s/ntime 
        outVars_h%nbp(itime)             = outVars_h%nbp(itime)      + &
                                           (vegn%gpp - st%Rhetero - vegn%Rauto) *convert_g2kg*convert_h2s/ntime   
        outVars_h%wetlandCH4(itime)      = outVars_h%wetlandCH4(itime) + st%simuCH4 *convert_g2kg*convert_h2s/ntime   
        outVars_h%wetlandCH4prod(itime)  = outVars_h%wetlandCH4prod(itime) + st%Pro_sum *convert_g2kg*convert_h2s/ntime 
        outVars_h%wetlandCH4cons(itime)  = outVars_h%wetlandCH4cons(itime) + st%Oxi_sum *convert_g2kg*convert_h2s/ntime 
        ! Carbon Pools  (KgC m-2)
        outVars_h%cLeaf(itime)           = outVars_h%cLeaf(itime) + st%QC(1)*convert_g2kg/ntime
        outVars_h%cStem(itime)           = outVars_h%cStem(itime) + st%QC(2)*convert_g2kg/ntime
        outVars_h%cRoot(itime)           = outVars_h%cRoot(itime) + st%QC(3)*convert_g2kg/ntime
        outVars_h%cOther(itime)          = outVars_h%cOther(itime)+ vegn%NSC*convert_g2kg/ntime
        outVars_h%cLitter(itime)         = outVars_h%cLitter(itime) + st%QC(4)*convert_g2kg/ntime
        outVars_h%cLitterCwd(itime)      = outVars_h%cLitterCwd(itime) + st%QC(5)*convert_g2kg/ntime
        outVars_h%cSoil(itime)           = outVars_h%cSoil(itime) + (st%QC(6) + st%QC(7) + st%QC(8))*convert_g2kg/ntime
        outVars_h%cSoilLevels(itime, :)  = outVars_h%cSoilLevels(itime, :) + (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)/ntime
        outVars_h%cSoilFast(itime)       = outVars_h%cSoilFast(itime) + st%QC(6)*convert_g2kg/ntime 
        outVars_h%cSoilSlow(itime)       = outVars_h%cSoilSlow(itime) + st%QC(7)*convert_g2kg/ntime 
        outVars_h%cSoilPassive(itime)    = outVars_h%cSoilPassive(itime) + st%QC(8)*convert_g2kg/ntime 
        outVars_h%CH4(itime, :)          = outVars_h%CH4(itime, :) + st%CH4*convert_g2kg/ntime 
        ! Nitrogen fluxes (kgN m-2 s-1)
        outVars_h%fBNF(itime)            = outVars_h%fBNF(itime) + st%N_fixation*convert_g2kg*convert_h2s/ntime 
        outVars_h%fN2O(itime)            = outVars_h%fN2O(itime) + &
                                           (st%N_transfer+st%N_uptake+st%N_fixation)*convert_g2kg*convert_h2s/ntime
        outVars_h%fNloss(itime)          = outVars_h%fNloss(itime) + &
                                            (vegn%N_leaf+vegn%N_wood+vegn%N_root)*convert_g2kg*convert_h2s/ntime
        outVars_h%fNnetmin(itime)        = outVars_h%fNnetmin(itime) + st%fNnetmin*convert_g2kg*convert_h2s/ntime 
        outVars_h%fNdep(itime)           = outVars_h%fNdep(itime) + st%N_deposit*convert_g2kg*convert_h2s/ntime 
        ! Nitrogen pools (kgN m-2)
        outVars_h%nLeaf(itime)           = outVars_h%nLeaf(itime) + st%QN(1)*convert_g2kg/ntime
        outVars_h%nStem(itime)           = outVars_h%nStem(itime) + st%QN(2)*convert_g2kg/ntime
        outVars_h%nRoot(itime)           = outVars_h%nRoot(itime) + st%QN(3)*convert_g2kg/ntime
        outVars_h%nOther(itime)          = outVars_h%nOther(itime) + vegn%NSN*convert_g2kg/ntime
        outVars_h%nLitter(itime)         = outVars_h%nLitter(itime) + st%QN(4)*convert_g2kg/ntime
        outVars_h%nLitterCwd(itime)      = outVars_h%nLitterCwd(itime) + st%QN(5)*convert_g2kg/ntime
        outVars_h%nSoil(itime)           = outVars_h%nSoil(itime) + (st%QN(6)+st%QN(7)+st%QN(8))*convert_g2kg/ntime
        outVars_h%nMineral(itime)        = outVars_h%nMineral(itime) + st%QNminer*convert_g2kg/ntime 
        ! energy fluxes (W m-2)
        outVars_h%hfls(itime)            = outVars_h%hfls(itime) + st%Hsoil/ntime ! Sensible heat flux;
        outVars_h%hfss(itime)            = outVars_h%hfss(itime) + st%Esoil/ntime ! Latent heat flux;
        outVars_h%SWnet(itime)           = outVars_h%SWnet(itime) + 0/ntime       ! Net shortwave radiation;
        outVars_h%LWnet(itime)           = outVars_h%LWnet(itime) + 0/ntime       ! Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        outVars_h%ec(itime)              = outVars_h%ec(itime)    + 0/ntime       !evap*convert_g2kg*convert_h2s/ntime        ! Canopy evaporation;
        outVars_h%tran(itime)            = outVars_h%tran(itime)  + vegn%transp*convert_g2kg*convert_h2s/ntime      ! Canopy transpiration;
        outVars_h%es(itime)              = outVars_h%es(itime)    + st%evap*convert_g2kg*convert_h2s/ntime ! Soil evaporation
        outVars_h%hfsbl(itime)           = outVars_h%hfsbl(itime) + st%sublim*convert_g2kg*convert_h2s/ntime ! Snow sublimation
        outVars_h%mrro(itime)            = outVars_h%mrro(itime)  + st%runoff*convert_g2kg*convert_h2s/ntime
        ! outVars_h%mrros(itime)         = forcing(iforcing)%Rain    
        outVars_h%mrrob(itime)           = outVars_h%mrrob(itime) + 0/ntime ! Total runoff; Surface runoff; Subsurface runoff
        ! other
        outVars_h%mrso(itime, :)         = outVars_h%mrso(itime, :) + st%liq_water*1000/ntime  ! Kg m-2, soil moisture in each soil layer
        outVars_h%tsl(itime,:)           = outVars_h%tsl(itime,:) + (st%tsoil_layer(1:10)+273.15)/ntime                            ! K, soil temperature in each soil layer Jian: not sure the tsoil_layer is correct or not
        ! outVars_h%tsland(itime)        = forcing(iforcing)%Tair+273.15                                   ! K, surface temperature
        outVars_h%wtd(itime)             = outVars_h%wtd(itime) +  (st%zwt/1000)/ntime                                       ! m, Water table depth
        outVars_h%snd(itime)             = outVars_h%snd(itime) +  (st%snow_depth/100)/ntime                               ! m, Total snow depth, Jian: change from m to cm in code, and now change from cm to m
        outVars_h%lai(itime)             = outVars_h%lai(itime) +  (vegn%LAI)/ntime                                           ! m2 m-2, Leaf area index
    end subroutine updateOutVars

    ! subroutine updateHourly(vegn, itime)
    !     implicit none
    !     type(vegn_tile_type), intent(in) :: vegn
    !     integer, intent(in) :: itime
    !     integer :: ipft, npft
    !     ! integer iTotHourly
    !     convert_g2kg = 0.001
    !     convert_h2s  = 1/3600.
    !     if (allocated(outVars_h%allSpec)) then
    !         npft = size(outVars_h%allSpec)
    !         do ipft = 1, npft
    !             ! carbon fluxes (Kg C m-2 s-1)
    !             outVars_h%allSpec(ipft)%gpp(itime)      = vegn%allSp(ipft)%gpp*convert_g2kg*convert_h2s
    !             ! outVars_h%allSpec(ipft)%nee      = vegn%allSp(ipft)%NEE
    !             outVars_h%allSpec(ipft)%npp(itime)      = vegn%allSp(ipft)%npp*convert_g2kg*convert_h2s
    !             outVars_h%allSpec(ipft)%nppLeaf(itime)  = vegn%allSp(ipft)%NPP_L*convert_g2kg*convert_h2s
    !             outVars_h%allSpec(ipft)%nppWood(itime)  = vegn%allSp(ipft)%NPP_W*convert_g2kg*convert_h2s
    !             outVars_h%allSpec(ipft)%nppStem(itime)  = vegn%allSp(ipft)%NPP_W*convert_g2kg*convert_h2s
    !             outVars_h%allSpec(ipft)%nppRoot(itime)  = vegn%allSp(ipft)%NPP_R*convert_g2kg*convert_h2s
    !             ! outVars_h%allSpec(ipft)%nppOther = vegn%allSp(ipft)%nppOther    ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
    !             outVars_h%allSpec(ipft)%ra(itime)       = vegn%allSp(ipft)%Rauto*convert_g2kg*convert_h2s
    !             outVars_h%allSpec(ipft)%raLeaf(itime)   = vegn%allSp(ipft)%RmLeaf*convert_g2kg*convert_h2s
    !             outVars_h%allSpec(ipft)%raStem(itime)   = vegn%allSp(ipft)%RmStem*convert_g2kg*convert_h2s
    !             outVars_h%allSpec(ipft)%raRoot(itime)   = vegn%allSp(ipft)%RmRoot*convert_g2kg*convert_h2s
    !             ! outVars_h%allSpec(ipft)%raOther  = vegn%allSp(ipft)%
    !             outVars_h%allSpec(ipft)%rMaint(itime)   = vegn%allSp(ipft)%Rmain*convert_g2kg*convert_h2s
    !             outVars_h%allSpec(ipft)%rGrowth(itime)  = vegn%allSp(ipft)%Rgrowth*convert_g2kg*convert_h2s
    !             ! outVars_h%allSpec(ipft)%nbp      = vegn%allSp(ipft)%
    !             ! Carbon Pools  (KgC m-2)
    !             outVars_h%allSpec(ipft)%cLeaf(itime)    = vegn%allSp(ipft)%QC(1)*convert_g2kg
    !             outVars_h%allSpec(ipft)%cStem(itime)    = vegn%allSp(ipft)%QC(2)*convert_g2kg
    !             outVars_h%allSpec(ipft)%cRoot(itime)    = vegn%allSp(ipft)%QC(3)*convert_g2kg
    !             ! Nitrogen pools (kgN m-2)
    !             outVars_h%allSpec(ipft)%nLeaf(itime)    = vegn%allSp(ipft)%QN(1)*convert_g2kg
    !             outVars_h%allSpec(ipft)%nStem(itime)    = vegn%allSp(ipft)%QN(2)*convert_g2kg
    !             outVars_h%allSpec(ipft)%nRoot(itime)    = vegn%allSp(ipft)%QN(3)*convert_g2kg
    !             ! water fluxes (kg m-2 s-1)
    !             outVars_h%allSpec(ipft)%tran(itime)     = vegn%allSp(ipft)%transp*convert_g2kg*convert_h2s
    !             ! other
    !             outVars_h%allSpec(ipft)%lai(itime)      = vegn%allSp(ipft)%LAI
    !         enddo
    !     endif
    !     ! carbon fluxes (KgC m-2 s-1) Jian: TECO unit is gC m-2 h-1
    !     outVars_h%gpp(itime)             = vegn%gpp*convert_g2kg*convert_h2s
    !     outVars_h%npp(itime)             = vegn%npp*convert_g2kg*convert_h2s
    !     outVars_h%nppLeaf(itime)         = vegn%NPP_L*convert_g2kg*convert_h2s
    !     outVars_h%nppWood(itime)         = vegn%NPP_W*convert_g2kg*convert_h2s  
    !     outVars_h%nppStem(itime)         = vegn%NPP_W*convert_g2kg*convert_h2s                   ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues. Jian: TECO has no above ground woody tissues, set to equit wood
    !     outVars_h%nppRoot(itime)         = vegn%NPP_R*convert_g2kg*convert_h2s
    !     outVars_h%nppOther(itime)        = 0*convert_g2kg*convert_h2s                       ! Jian: no other storage, NSC seems different from other NPP.
    !     outVars_h%ra(itime)              = vegn%Rauto*convert_g2kg*convert_h2s
    !     outVars_h%raLeaf(itime)          = vegn%Rmleaf*convert_g2kg*convert_h2s
    !     outVars_h%raStem(itime)          = vegn%Rmstem*convert_g2kg*convert_h2s
    !     outVars_h%raRoot(itime)          = vegn%Rmroot*convert_g2kg*convert_h2s
    !     outVars_h%raOther(itime)         = st%Rnitrogen *convert_g2kg*convert_h2s               ! Total C cost for nitrogen
    !     outVars_h%rMaint(itime)          = vegn%Rmain *convert_g2kg*convert_h2s                   ! maintenance respiration
    !     outVars_h%rGrowth(itime)         = vegn%Rgrowth *convert_g2kg*convert_h2s                 ! growth respiration
    !     outVars_h%rh(itime)              = st%Rhetero *convert_g2kg*convert_h2s                 ! heterotrophic respiration
    !     outVars_h%nbp(itime)             = (vegn%gpp - st%Rhetero - vegn%Rauto) *convert_g2kg*convert_h2s   ! NBP(net biome productivity) = GPP - Rh - Ra - other losses  
    !     outVars_h%wetlandCH4(itime)      = st%simuCH4 *convert_g2kg*convert_h2s                 ! wetland net fluxes of CH4
    !     outVars_h%wetlandCH4prod(itime)  = st%Pro_sum *convert_g2kg*convert_h2s                 ! wetland net fluxes of CH4 production
    !     outVars_h%wetlandCH4cons(itime)  = st%Oxi_sum *convert_g2kg*convert_h2s                 ! wetland net fluxes of CH4 consumption
    !     ! Carbon Pools  (KgC m-2)
    !     outVars_h%cLeaf(itime)           = st%QC(1)*convert_g2kg
    !     outVars_h%cStem(itime)           = st%QC(2)*convert_g2kg
    !     outVars_h%cRoot(itime)           = st%QC(3)*convert_g2kg
    !     outVars_h%cOther(itime)          = vegn%NSC*convert_g2kg                        ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
    !     outVars_h%cLitter(itime)         = st%QC(4)*convert_g2kg                      ! litter (excluding coarse woody debris), Jian: fine litter in TECO?
    !     outVars_h%cLitterCwd(itime)      = st%QC(5)*convert_g2kg                      ! cLitterCwd: carbon in coarse woody debris
    !     outVars_h%cSoil(itime)           = (st%QC(6) + st%QC(7) + st%QC(8))*convert_g2kg    ! cSoil: soil organic carbon (Jian: total soil carbon);
    !     outVars_h%cSoilLevels(itime, :)     = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)                                       ! cSoilLevels(depth-specific soil organic carbon, Jian: depth?);
    !     outVars_h%cSoilFast(itime)       = st%QC(6)*convert_g2kg                      ! cSoilPools (different pools without depth)
    !     outVars_h%cSoilSlow(itime)       = st%QC(7)*convert_g2kg 
    !     outVars_h%cSoilPassive(itime)    = st%QC(8)*convert_g2kg 
    !     outVars_h%CH4(itime, :)          = st%CH4*convert_g2kg                        ! methane concentration
    !     ! Nitrogen fluxes (kgN m-2 s-1)
    !     outVars_h%fBNF(itime)            = st%N_fixation*convert_g2kg*convert_h2s                 ! fBNF: biological nitrogen fixation;
    !     outVars_h%fN2O(itime)            = (st%N_transfer+st%N_uptake+st%N_fixation)*convert_g2kg*convert_h2s                    ! fN2O: loss of nitrogen through emission of N2O;
    !     outVars_h%fNloss(itime)          = (vegn%N_leaf+vegn%N_wood+vegn%N_root)*convert_g2kg*convert_h2s                     ! fNloss:Total loss of nitrogen to the atmosphere and from leaching;
    !     outVars_h%fNnetmin(itime)        = st%fNnetmin*convert_g2kg*convert_h2s !((N_transfer+N_uptake+N_fixation)-(N_leaf+N_wood+N_root))*convert_g2kg*convert_h2s                   ! net mineralizaiton
    !     outVars_h%fNdep(itime)           = st%N_deposit*convert_g2kg*convert_h2s                  ! deposition of N
    !     ! Nitrogen pools (kgN m-2)
    !     outVars_h%nLeaf(itime)           = st%QN(1)*convert_g2kg
    !     outVars_h%nStem(itime)           = st%QN(2)*convert_g2kg
    !     outVars_h%nRoot(itime)           = st%QN(3)*convert_g2kg
    !     outVars_h%nOther(itime)          = vegn%NSN*convert_g2kg         ! other N pool
    !     outVars_h%nLitter(itime)         = st%QN(4)*convert_g2kg
    !     outVars_h%nLitterCwd(itime)      = st%QN(5)*convert_g2kg
    !     outVars_h%nSoil(itime)           = (st%QN(6)+st%QN(7)+st%QN(8))*convert_g2kg
    !     outVars_h%nMineral(itime)        = st%QNminer*convert_g2kg                                 ! nMineral: Mineral nitrogen pool
    !     ! energy fluxes (W m-2)
    !     outVars_h%hfls(itime)            = st%Hsoil ! Sensible heat flux;
    !     outVars_h%hfss(itime)            = st%Esoil ! Latent heat flux;
    !     outVars_h%SWnet(itime)           = 0 ! Net shortwave radiation;
    !     outVars_h%LWnet(itime)           = 0 ! Net longwave radiation
    !     ! water fluxes (kg m-2 s-1)
    !     outVars_h%ec(itime)              = 0!evap*convert_g2kg*convert_h2s        ! Canopy evaporation;
    !     outVars_h%tran(itime)            = vegn%transp*convert_g2kg*convert_h2s      ! Canopy transpiration;
    !     outVars_h%es(itime)              = st%evap*convert_g2kg*convert_h2s ! Soil evaporation
    !     outVars_h%hfsbl(itime)           = st%sublim*convert_g2kg*convert_h2s ! Snow sublimation
    !     outVars_h%mrro(itime)            = st%runoff*convert_g2kg*convert_h2s
    !     ! outVars_h%mrros(itime)           = forcing(iforcing)%Rain    
    !     outVars_h%mrrob(itime)           = 0 ! Total runoff; Surface runoff; Subsurface runoff
    !     ! other
    !     outVars_h%mrso(itime, :)         = st%liq_water*1000                                   ! Kg m-2, soil moisture in each soil layer
    !     outVars_h%tsl(itime,:)             = st%tsoil_layer(1:10)+273.15                            ! K, soil temperature in each soil layer Jian: not sure the tsoil_layer is correct or not
    !     ! outVars_h%tsland(itime)          = forcing(iforcing)%Tair+273.15                                   ! K, surface temperature
    !     outVars_h%wtd(itime)             = st%zwt/1000                                       ! m, Water table depth
    !     outVars_h%snd(itime)             = st%snow_depth/100                                ! m, Total snow depth, Jian: change from m to cm in code, and now change from cm to m
    !     outVars_h%lai(itime)             = vegn%LAI                                           ! m2 m-2, Leaf area index
        
    ! end subroutine updateHourly
end module update_and_summary