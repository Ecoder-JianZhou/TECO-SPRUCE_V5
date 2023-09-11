! This module is used to summary the values of hourly, daily, monthly and yearly
module update_and_summary
    use datatypes
    implicit NONE
    real convert_g2kg, convert_h2s

    contains
    subroutine updateHourly(vegn)
        implicit none
        type(vegn_tile_type), intent(in) :: vegn
        ! integer iTotHourly
        convert_g2kg = 0.001
        convert_h2s  = 1/3600.
        ! carbon fluxes (KgC m-2 s-1) Jian: TECO unit is gC m-2 h-1
        outVars_h%gpp             = vegn%gpp*convert_g2kg*convert_h2s
        outVars_h%npp             = npp*convert_g2kg*convert_h2s
        outVars_h%nppLeaf         = NPP_L*convert_g2kg*convert_h2s
        outVars_h%nppWood         = NPP_W*convert_g2kg*convert_h2s  
        outVars_h%nppStem         = NPP_W*convert_g2kg*convert_h2s                   ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues. Jian: TECO has no above ground woody tissues, set to equit wood
        outVars_h%nppRoot         = NPP_R*convert_g2kg*convert_h2s
        outVars_h%nppOther        = 0*convert_g2kg*convert_h2s                       ! Jian: no other storage, NSC seems different from other NPP.
        outVars_h%ra              = Rauto*convert_g2kg*convert_h2s
        outVars_h%raLeaf          = Rmleaf*convert_g2kg*convert_h2s
        outVars_h%raStem          = Rmstem*convert_g2kg*convert_h2s
        outVars_h%raRoot          = Rmroot*convert_g2kg*convert_h2s
        outVars_h%raOther         = Rnitrogen *convert_g2kg*convert_h2s               ! Total C cost for nitrogen
        outVars_h%rMaint          = Rmain *convert_g2kg*convert_h2s                   ! maintenance respiration
        outVars_h%rGrowth         = Rgrowth *convert_g2kg*convert_h2s                 ! growth respiration
        outVars_h%rh              = Rhetero *convert_g2kg*convert_h2s                 ! heterotrophic respiration
        outVars_h%nbp             = (gpp - Rhetero - Rauto) *convert_g2kg*convert_h2s   ! NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        outVars_h%wetlandCH4      = simuCH4 *convert_g2kg*convert_h2s                 ! wetland net fluxes of CH4
        outVars_h%wetlandCH4prod  = Pro_sum *convert_g2kg*convert_h2s                 ! wetland net fluxes of CH4 production
        outVars_h%wetlandCH4cons  = Oxi_sum *convert_g2kg*convert_h2s                 ! wetland net fluxes of CH4 consumption
        ! Carbon Pools  (KgC m-2)
        outVars_h%cLeaf           = QC(1)*convert_g2kg
        outVars_h%cStem           = QC(2)*convert_g2kg
        outVars_h%cRoot           = QC(3)*convert_g2kg
        outVars_h%cOther          = NSC*convert_g2kg                        ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        outVars_h%cLitter         = QC(4)*convert_g2kg                      ! litter (excluding coarse woody debris), Jian: fine litter in TECO?
        outVars_h%cLitterCwd      = QC(5)*convert_g2kg                      ! cLitterCwd: carbon in coarse woody debris
        outVars_h%cSoil           = (QC(6) + QC(7) + QC(8))*convert_g2kg    ! cSoil: soil organic carbon (Jian: total soil carbon);
        outVars_h%cSoilLevels     = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)                                       ! cSoilLevels(depth-specific soil organic carbon, Jian: depth?);
        outVars_h%cSoilFast       = QC(6)*convert_g2kg                      ! cSoilPools (different pools without depth)
        outVars_h%cSoilSlow       = QC(7)*convert_g2kg 
        outVars_h%cSoilPassive    = QC(8)*convert_g2kg 
        outVars_h%CH4            = CH4*convert_g2kg                        ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        outVars_h%fBNF            = N_fixation*convert_g2kg*convert_h2s                 ! fBNF: biological nitrogen fixation;
        outVars_h%fN2O            = (N_transfer+N_uptake+N_fixation)*convert_g2kg*convert_h2s                    ! fN2O: loss of nitrogen through emission of N2O;
        outVars_h%fNloss          = (N_leaf+N_wood+N_root)*convert_g2kg*convert_h2s                     ! fNloss:Total loss of nitrogen to the atmosphere and from leaching;
        outVars_h%fNnetmin        = ((N_transfer+N_uptake+N_fixation)-(N_leaf+N_wood+N_root))*convert_g2kg*convert_h2s                   ! net mineralizaiton
        outVars_h%fNdep           = N_wood*convert_g2kg*convert_h2s                  ! deposition of N
        ! Nitrogen pools (kgN m-2)
        outVars_h%nLeaf           = QN(1)*convert_g2kg
        outVars_h%nStem           = QN(2)*convert_g2kg
        outVars_h%nRoot           = QN(3)*convert_g2kg
        outVars_h%nOther          = NSN*convert_g2kg         ! other N pool
        outVars_h%nLitter         = QN(4)*convert_g2kg
        outVars_h%nLitterCwd      = QN(5)*convert_g2kg
        outVars_h%nSoil           = (QN(6)+QN(7)+QN(8))*convert_g2kg
        outVars_h%nMineral        = QNminer*convert_g2kg                                 ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        outVars_h%hfls            = Hsoil ! Sensible heat flux;
        outVars_h%hfss            = Esoil ! Latent heat flux;
        outVars_h%SWnet           = 0 ! Net shortwave radiation;
        outVars_h%LWnet           = 0 ! Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        outVars_h%ec              = 0!evap*convert_g2kg*convert_h2s        ! Canopy evaporation;
        outVars_h%tran            = transp*convert_g2kg*convert_h2s      ! Canopy transpiration;
        outVars_h%es              = evap*convert_g2kg*convert_h2s ! Soil evaporation
        outVars_h%hfsbl           = sublim*convert_g2kg*convert_h2s ! Snow sublimation
        outVars_h%mrro            = runoff*convert_g2kg*convert_h2s
        outVars_h%mrros           = forcing(iforcing)%Rain    
        outVars_h%mrrob           = 0 ! Total runoff; Surface runoff; Subsurface runoff
        ! other
        outVars_h%mrso            = liq_water*1000                                   ! Kg m-2, soil moisture in each soil layer
        outVars_h%tsl             = tsoil_layer(1:10)+273.15                            ! K, soil temperature in each soil layer Jian: not sure the tsoil_layer is correct or not
        outVars_h%tsland          = forcing(iforcing)%Tair+273.15                                   ! K, surface temperature
        outVars_h%wtd             = zwt/1000                                       ! m, Water table depth
        outVars_h%snd             = snow_depth/100                                ! m, Total snow depth, Jian: change from m to cm in code, and now change from cm to m
        outVars_h%lai             = LAI                                           ! m2 m-2, Leaf area index
        
    end subroutine updateHourly

    subroutine updateDaily()
        implicit none
        ! carbon fluxes
        outVars_d%gpp             = outVars_d%gpp            + outVars_h%gpp/24
        outVars_d%npp             = outVars_d%npp            + outVars_h%npp/24
        outVars_d%nppLeaf         = outVars_d%nppLeaf        + outVars_h%nppLeaf/24
        outVars_d%nppWood         = outVars_d%nppWood        + outVars_h%nppWood/24
        outVars_d%nppStem         = outVars_d%nppStem        + outVars_h%nppStem/24                      ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        outVars_d%nppRoot         = outVars_d%nppRoot        + outVars_h%nppRoot/24
        outVars_d%nppOther        = outVars_d%nppOther       + outVars_h%nppOther/24
        outVars_d%ra              = outVars_d%ra             + outVars_h%ra/24
        outVars_d%raLeaf          = outVars_d%raLeaf         + outVars_h%raLeaf/24
        outVars_d%raStem          = outVars_d%raStem         + outVars_h%raStem/24
        outVars_d%raRoot          = outVars_d%raRoot         + outVars_h%raRoot/24
        outVars_d%raOther         = outVars_d%raOther        + outVars_h%raOther/24
        outVars_d%rMaint          = outVars_d%rMaint         + outVars_h%rMaint/24                           ! maintenance respiration
        outVars_d%rGrowth         = outVars_d%Rgrowth        + outVars_h%rGrowth/24                         ! growth respiration
        outVars_d%rh              = outVars_d%rh             + outVars_h%rh/24                                   ! heterotrophic respiration
        outVars_d%nbp             = outVars_d%nbp            + outVars_h%nbp/24                                 ! NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        outVars_d%wetlandCH4      = outVars_d%wetlandCH4     + outVars_h%wetlandCH4/24                   ! wetland net fluxes of CH4
        outVars_d%wetlandCH4prod  = outVars_d%wetlandCH4prod + outVars_h%wetlandCH4prod/24           ! wetland net fluxes of CH4 production
        outVars_d%wetlandCH4cons  = outVars_d%wetlandCH4cons + outVars_h%wetlandCH4cons/24           ! wetland net fluxes of CH4 consumption
        ! Carbon Pools  (KgC m-2)
        outVars_d%cLeaf           = outVars_d%cLeaf          + outVars_h%cLeaf/24
        outVars_d%cStem           = outVars_d%cStem          + outVars_h%cStem/24
        outVars_d%cRoot           = outVars_d%cRoot          + outVars_h%cRoot/24
        outVars_d%cOther          = outVars_d%cOther         + outVars_h%cOther/24                                    ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        outVars_d%cLitter         = outVars_d%cLitter        + outVars_h%cLitter/24                                   ! litter (excluding coarse woody debris), Jian: fine litter in TECO?
        outVars_d%cLitterCwd      = outVars_d%cLitterCwd     + outVars_h%cLitterCwd/24                                ! cLitterCwd: carbon in coarse woody debris
        outVars_d%cSoil           = outVars_d%cSoil          + outVars_h%cSoil/24                                     ! cSoil: soil organic carbon (Jian: total soil carbon);
        outVars_d%cSoilLevels     = outVars_d%cSoilLevels    + outVars_h%cSoilLevels/24                               ! cSoilLevels(depth-specific soil organic carbon, Jian: depth?);
        outVars_d%cSoilFast       = outVars_d%cSoilFast      + outVars_h%cSoilFast/24                                 ! cSoilPools (different pools without depth)
        outVars_d%cSoilSlow       = outVars_d%cSoilSlow      + outVars_h%cSoilSlow/24 
        outVars_d%cSoilPassive    = outVars_d%cSoilPassive   + outVars_h%cSoilPassive/24 
        outVars_d%CH4            = outVars_d%CH4           + outVars_h%CH4/24                                      ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        outVars_d%fBNF            = outVars_d%fBNF           + outVars_h%fBNF/24                               ! fBNF: biological nitrogen fixation;
        outVars_d%fN2O            = outVars_d%fN2O           + outVars_h%fN2O/24                               ! fN2O: loss of nitrogen through emission of N2O;
        outVars_d%fNloss          = outVars_d%fNloss         + outVars_h%fNloss/24                           ! fNloss:Total loss of nitrogen to the atmosphere and from leaching;
        outVars_d%fNnetmin        = outVars_d%fNnetmin       + outVars_h%fNnetmin/24                       ! net mineralizaiton
        outVars_d%fNdep           = outVars_d%fNdep          + outVars_h%fNdep/24                             ! deposition of N
        ! Nitrogen pools (kgN m-2)
        outVars_d%nLeaf           = outVars_d%nLeaf          + outVars_h%nLeaf/24
        outVars_d%nStem           = outVars_d%nStem          + outVars_h%nStem/24
        outVars_d%nRoot           = outVars_d%nRoot          + outVars_h%nRoot/24
        outVars_d%nOther          = outVars_d%nOther         + outVars_h%nOther/24
        outVars_d%nLitter         = outVars_d%nLitter        + outVars_h%nLitter/24
        outVars_d%nLitterCwd      = outVars_d%nLitterCwd     + outVars_h%nLitter/24
        outVars_d%nSoil           = outVars_d%nSoil          + outVars_h%nSoil/24
        outVars_d%nMineral        = outVars_d%nMineral       + outVars_h%nMineral/24                                  ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        outVars_d%hfls            = outVars_d%hfls           + outVars_h%hfls/24                               ! Sensible heat flux;
        outVars_d%hfss            = outVars_d%hfss           + outVars_h%hfss/24                               ! Latent heat flux;
        outVars_d%SWnet           = outVars_d%SWnet          + outVars_h%SWnet/24                             ! Net shortwave radiation;
        outVars_d%LWnet           = outVars_d%LWnet          + outVars_h%LWnet/24                             ! Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        outVars_d%ec              = outVars_d%ec             + outVars_h%ec/24
        outVars_d%tran            = outVars_d%tran           + outVars_h%tran/24
        outVars_d%es              = outVars_d%es             + outVars_h%es/24                                   ! Canopy evaporation; Canopy transpiration; Soil evaporation
        outVars_d%hfsbl           = outVars_d%hfsbl          + outVars_h%hfsbl/24                             ! Snow sublimation
        outVars_d%mrro            = outVars_d%mrro           + outVars_h%mrro/24
        outVars_d%mrros           = outVars_d%mrros          + outVars_h%mrros/24
        outVars_d%mrrob           = outVars_d%mrrob          + outVars_h%mrrob/24                              ! Total runoff; Surface runoff; Subsurface runoff
        ! other
        outVars_d%mrso            = outVars_d%mrso           + outVars_h%mrro/24                  ! Kg m-2, soil moisture in each soil layer
        outVars_d%tsl             = outVars_d%tsl            + outVars_h%tsl/24                   ! K, soil temperature in each soil layer
        outVars_d%tsland          = outVars_d%tsland         + outVars_h%tsland/24                ! K, surface temperature
        outVars_d%wtd             = outVars_d%wtd            + outVars_h%wtd/24                   ! m, Water table depth
        outVars_d%snd             = outVars_d%snd            + outVars_h%snd/24                   ! m, Total snow depth
        outVars_d%lai             = outVars_d%lai            + outVars_h%lai/24                   ! m2 m-2, Leaf area index   
    end subroutine updateDaily

    subroutine updateMonthly(hoursOfmonth)
        implicit none
        integer hoursOfmonth
        ! carbon fluxes
        outVars_m%gpp             = outVars_m%gpp            + outVars_h%gpp            /hoursOfmonth
        outVars_m%npp             = outVars_m%npp            + outVars_h%npp            /hoursOfmonth
        outVars_m%nppLeaf         = outVars_m%nppLeaf        + outVars_h%nppLeaf        /hoursOfmonth
        outVars_m%nppWood         = outVars_m%nppWood        + outVars_h%nppWood        /hoursOfmonth
        outVars_m%nppStem         = outVars_m%nppStem        + outVars_h%nppStem        /hoursOfmonth ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        outVars_m%nppRoot         = outVars_m%nppRoot        + outVars_h%nppRoot        /hoursOfmonth
        outVars_m%nppOther        = outVars_m%nppOther       + outVars_h%nppOther       /hoursOfmonth
        outVars_m%ra              = outVars_m%ra             + outVars_h%ra             /hoursOfmonth
        outVars_m%raLeaf          = outVars_m%raLeaf         + outVars_h%raLeaf         /hoursOfmonth
        outVars_m%raStem          = outVars_m%raStem         + outVars_h%raStem         /hoursOfmonth
        outVars_m%raRoot          = outVars_m%raRoot         + outVars_h%raRoot         /hoursOfmonth
        outVars_m%raOther         = outVars_m%raOther        + outVars_h%raOther        /hoursOfmonth
        outVars_m%rMaint          = outVars_m%rMaint         + outVars_h%rMaint         /hoursOfmonth ! maintenance respiration
        outVars_m%rGrowth         = outVars_m%rGrowth        + outVars_h%rGrowth        /hoursOfmonth ! growth respiration
        outVars_m%rh              = outVars_m%rh             + outVars_h%rh             /hoursOfmonth ! heterotrophic respiration
        outVars_m%nbp             = outVars_m%nbp            + outVars_h%nbp            /hoursOfmonth ! NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        outVars_m%wetlandCH4      = outVars_m%wetlandCH4     + outVars_h%wetlandCH4     /hoursOfmonth ! wetland net fluxes of CH4
        outVars_m%wetlandCH4prod  = outVars_m%wetlandCH4prod + outVars_h%wetlandCH4prod /hoursOfmonth ! wetland net fluxes of CH4 production
        outVars_m%wetlandCH4cons  = outVars_m%wetlandCH4cons + outVars_h%wetlandCH4cons /hoursOfmonth ! wetland net fluxes of CH4 consumption
        ! Carbon Pools  (KgC m-2)
        outVars_m%cLeaf           = outVars_m%cLeaf          + outVars_h%cLeaf          /hoursOfmonth     
        outVars_m%cStem           = outVars_m%cStem          + outVars_h%cStem          /hoursOfmonth
        outVars_m%cRoot           = outVars_m%cRoot          + outVars_h%cRoot          /hoursOfmonth      
        outVars_m%cOther          = outVars_m%cOther         + outVars_h%cOther         /hoursOfmonth ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        outVars_m%cLitter         = outVars_m%cLitter        + outVars_h%cLitter        /hoursOfmonth ! litter (excluding coarse woody debris), Jian: fine litter in TECO?
        outVars_m%cLitterCwd      = outVars_m%cLitterCwd     + outVars_h%cLitterCwd     /hoursOfmonth ! cLitterCwd: carbon in coarse woody debris
        outVars_m%cSoil           = outVars_m%cSoil          + outVars_h%cSoil          /hoursOfmonth ! cSoil: soil organic carbon (Jian: total soil carbon);
        outVars_m%cSoilLevels     = outVars_m%cSoilLevels    + outVars_h%cSoilLevels    /hoursOfmonth  ! cSoilLevels(depth-specific soil organic carbon, Jian: depth?);
        outVars_m%cSoilFast       = outVars_m%cSoilFast      + outVars_h%cSoilFast      /hoursOfmonth ! cSoilPools (different pools without depth)
        outVars_m%cSoilSlow       = outVars_m%cSoilSlow      + outVars_h%cSoilSlow      /hoursOfmonth
        outVars_m%cSoilPassive    = outVars_m%cSoilPassive   + outVars_h%cSoilPassive   /hoursOfmonth
        outVars_m%CH4            = outVars_m%CH4           + outVars_h%CH4           /hoursOfmonth ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        outVars_m%fBNF            = outVars_m%fBNF           + outVars_h%fBNF           /hoursOfmonth ! fBNF: biological nitrogen fixation;
        outVars_m%fN2O            = outVars_m%fN2O           + outVars_h%fN2O           /hoursOfmonth ! fN2O: loss of nitrogen through emission of N2O;
        outVars_m%fNloss          = outVars_m%fNloss         + outVars_h%fNloss         /hoursOfmonth ! fNloss:Total loss of nitrogen to the atmosphere and from leaching;
        outVars_m%fNnetmin        = outVars_m%fNnetmin       + outVars_h%fNnetmin       /hoursOfmonth ! net mineralizaiton
        outVars_m%fNdep           = outVars_m%fNdep          + outVars_h%fNdep          /hoursOfmonth ! deposition of N
        ! Nitrogen pools (kgN m-2)
        outVars_m%nLeaf           = outVars_m%nLeaf          + outVars_h%nLeaf          /hoursOfmonth
        outVars_m%nStem           = outVars_m%nStem          + outVars_h%nStem          /hoursOfmonth
        outVars_m%nRoot           = outVars_m%nRoot          + outVars_h%nRoot          /hoursOfmonth
        outVars_m%nOther          = outVars_m%nOther         + outVars_h%nOther         /hoursOfmonth
        outVars_m%nLitter         = outVars_m%nLitter        + outVars_h%nLitter        /hoursOfmonth
        outVars_m%nLitterCwd      = outVars_m%nLitterCwd     + outVars_h%nLitterCwd     /hoursOfmonth
        outVars_m%nSoil           = outVars_m%nSoil          + outVars_h%nSoil          /hoursOfmonth
        outVars_m%nMineral        = outVars_m%nMineral       + outVars_h%nMineral       /hoursOfmonth ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        outVars_m%hfls            = outVars_m%hfls           + outVars_h%hfls           /hoursOfmonth ! Sensible heat flux;
        outVars_m%hfss            = outVars_m%hfss           + outVars_h%hfss           /hoursOfmonth ! Latent heat flux;
        outVars_m%SWnet           = outVars_m%SWnet          + outVars_h%SWnet          /hoursOfmonth ! Net shortwave radiation;
        outVars_m%LWnet           = outVars_m%LWnet          + outVars_h%LWnet          /hoursOfmonth !    Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        outVars_m%ec              = outVars_m%ec             + outVars_h%ec             /hoursOfmonth
        outVars_m%tran            = outVars_m%tran           + outVars_h%tran           /hoursOfmonth
        outVars_m%es              = outVars_m%es             + outVars_h%es             /hoursOfmonth ! Canopy evaporation; Canopy transpiration; Soil evaporation
        outVars_m%hfsbl           = outVars_m%hfsbl          + outVars_h%hfsbl          /hoursOfmonth ! Snow sublimation
        outVars_m%mrro            = outVars_m%mrro           + outVars_h%mrro           /hoursOfmonth
        outVars_m%mrros           = outVars_m%mrros          + outVars_h%mrros          /hoursOfmonth
        outVars_m%mrrob           = outVars_m%mrrob          + outVars_h%mrrob          /hoursOfmonth ! Total runoff; Surface runoff; Subsurface runoff
        ! other
        outVars_m%mrso            = outVars_m%mrso           + outVars_h%mrro           /hoursOfmonth                  ! Kg m-2, soil moisture in each soil layer
        outVars_m%tsl             = outVars_m%tsl            + outVars_h%tsl            /hoursOfmonth                   ! K, soil temperature in each soil layer
        outVars_m%tsland          = outVars_m%tsland         + outVars_h%tsland         /hoursOfmonth                ! K, surface temperature
        outVars_m%wtd             = outVars_m%wtd            + outVars_h%wtd            /hoursOfmonth                   ! m, Water table depth
        outVars_m%snd             = outVars_m%snd            + outVars_h%snd            /hoursOfmonth                   ! m, Total snow depth
        outVars_m%lai             = outVars_m%lai            + outVars_h%lai            /hoursOfmonth                   ! m2 m-2, Leaf area index
    end subroutine updateMonthly

    subroutine updateYearly(hoursOfYear)
        implicit none
        integer hoursOfYear
        ! carbon fluxes
        outVars_y%gpp             = outVars_y%gpp            + outVars_h%gpp            /hoursOfYear
        outVars_y%npp             = outVars_y%npp            + outVars_h%npp            /hoursOfYear
        outVars_y%nppLeaf         = outVars_y%nppLeaf        + outVars_h%nppLeaf        /hoursOfYear
        outVars_y%nppWood         = outVars_y%nppWood        + outVars_h%nppWood        /hoursOfYear
        outVars_y%nppStem         = outVars_y%nppStem        + outVars_h%nppStem        /hoursOfYear        ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        outVars_y%nppRoot         = outVars_y%nppRoot        + outVars_h%nppRoot        /hoursOfYear
        outVars_y%nppOther        = outVars_y%nppOther       + outVars_h%nppOther       /hoursOfYear
        outVars_y%ra              = outVars_y%ra             + outVars_h%ra             /hoursOfYear
        outVars_y%raLeaf          = outVars_y%raLeaf         + outVars_h%raLeaf         /hoursOfYear
        outVars_y%raStem          = outVars_y%raStem         + outVars_h%raStem         /hoursOfYear
        outVars_y%raRoot          = outVars_y%raRoot         + outVars_h%raRoot         /hoursOfYear
        outVars_y%raOther         = outVars_y%raOther        + outVars_h%raOther        /hoursOfYear
        outVars_y%rMaint          = outVars_y%rMaint         + outVars_h%rMaint         /hoursOfYear        ! maintenance respiration
        outVars_y%rGrowth         = outVars_y%rGrowth        + outVars_h%rGrowth        /hoursOfYear        ! growth respiration
        outVars_y%rh              = outVars_y%rh             + outVars_h%rh             /hoursOfYear        ! heterotrophic respiration
        outVars_y%nbp             = outVars_y%nbp            + outVars_h%nbp            /hoursOfYear        ! NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        outVars_y%wetlandCH4      = outVars_y%wetlandCH4     + outVars_h%wetlandCH4     /hoursOfYear        ! wetland net fluxes of CH4
        outVars_y%wetlandCH4prod  = outVars_y%wetlandCH4prod + outVars_h%wetlandCH4prod /hoursOfYear        ! wetland net fluxes of CH4 production
        outVars_y%wetlandCH4cons  = outVars_y%wetlandCH4cons + outVars_h%wetlandCH4cons /hoursOfYear        ! wetland net fluxes of CH4 consumption
        ! Carbon Pools  (KgC m-2)
        outVars_y%cLeaf           = outVars_y%cLeaf          + outVars_h%cLeaf          /hoursOfYear
        outVars_y%cStem           = outVars_y%cStem          + outVars_h%cStem          /hoursOfYear
        outVars_y%cRoot           = outVars_y%cRoot          + outVars_h%cRoot          /hoursOfYear
        outVars_y%cOther          = outVars_y%cOther         + outVars_h%cOther         /hoursOfYear        ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        outVars_y%cLitter         = outVars_y%cLitter        + outVars_h%cLitter        /hoursOfYear        ! litter (excluding coarse woody debris), Jian: fine litter in TECO?
        outVars_y%cLitterCwd      = outVars_y%cLitterCwd     + outVars_h%cLitterCwd     /hoursOfYear        ! cLitterCwd: carbon in coarse woody debris
        outVars_y%cSoil           = outVars_y%cSoil          + outVars_h%cSoil          /hoursOfYear        ! cSoil: soil organic carbon (Jian: total soil carbon);
        outVars_y%cSoilLevels     = outVars_y%cSoilLevels    + outVars_h%cSoilLevels    /hoursOfYear        ! cSoilLevels(depth-specific soil organic carbon, Jian: depth?);
        outVars_y%cSoilFast       = outVars_y%cSoilFast      + outVars_h%cSoilFast      /hoursOfYear        ! cSoilPools (different pools without depth)
        outVars_y%cSoilSlow       = outVars_y%cSoilSlow      + outVars_h%cSoilSlow      /hoursOfYear
        outVars_y%cSoilPassive    = outVars_y%cSoilPassive   + outVars_h%cSoilPassive   /hoursOfYear
        outVars_y%CH4            = outVars_y%CH4           + outVars_h%CH4           /hoursOfYear        ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        outVars_y%fBNF            = outVars_y%fBNF           + outVars_h%fBNF           /hoursOfYear        ! fBNF: biological nitrogen fixation;
        outVars_y%fN2O            = outVars_y%fN2O           + outVars_h%fN2O           /hoursOfYear        ! fN2O: loss of nitrogen through emission of N2O;
        outVars_y%fNloss          = outVars_y%fNloss         + outVars_h%fNloss         /hoursOfYear        ! fNloss:Total loss of nitrogen to the atmosphere and from leaching;
        outVars_y%fNnetmin        = outVars_y%fNnetmin       + outVars_h%fNnetmin       /hoursOfYear        ! net mineralizaiton
        outVars_y%fNdep           = outVars_y%fNdep          + outVars_h%fNdep          /hoursOfYear        ! deposition of N
        ! Nitrogen pools (kgN m-2)
        outVars_y%nLeaf           = outVars_y%nLeaf          + outVars_h%nLeaf          /hoursOfYear
        outVars_y%nStem           = outVars_y%nStem          + outVars_h%nStem          /hoursOfYear
        outVars_y%nRoot           = outVars_y%nRoot          + outVars_h%nRoot          /hoursOfYear
        outVars_y%nOther          = outVars_y%nOther         + outVars_h%nOther         /hoursOfYear
        outVars_y%nLitter         = outVars_y%nLitter        + outVars_h%nLitter        /hoursOfYear
        outVars_y%nLitterCwd      = outVars_y%nLitterCwd     + outVars_h%nLitterCwd     /hoursOfYear
        outVars_y%nSoil           = outVars_y%nSoil          + outVars_h%nSoil          /hoursOfYear
        outVars_y%nMineral        = outVars_y%nMineral       + outVars_h%nMineral       /hoursOfYear        ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        outVars_y%hfls            = outVars_y%hfls           + outVars_h%hfls           /hoursOfYear        ! Sensible heat flux;
        outVars_y%hfss            = outVars_y%hfss           + outVars_h%hfss           /hoursOfYear        ! Latent heat flux;
        outVars_y%SWnet           = outVars_y%SWnet          + outVars_h%SWnet          /hoursOfYear        ! Net shortwave radiation;
        outVars_y%LWnet           = outVars_y%LWnet          + outVars_h%LWnet          /hoursOfYear        ! Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        outVars_y%ec              = outVars_y%ec             + outVars_h%ec             /hoursOfYear
        outVars_y%tran            = outVars_y%tran           + outVars_h%tran           /hoursOfYear
        outVars_y%es              = outVars_y%es             + outVars_h%es             /hoursOfYear        ! Canopy evaporation; Canopy transpiration; Soil evaporation
        outVars_y%hfsbl           = outVars_y%hfsbl          + outVars_h%hfsbl          /hoursOfYear        ! Snow sublimation
        outVars_y%mrro            = outVars_y%mrro           + outVars_h%mrro           /hoursOfYear
        outVars_y%mrros           = outVars_y%mrros          + outVars_h%mrros          /hoursOfYear
        outVars_y%mrrob           = outVars_y%mrrob          + outVars_h%mrrob          /hoursOfYear        ! Total runoff; Surface runoff; Subsurface runoff
        ! other
        outVars_y%mrso            = outVars_y%mrso           + outVars_h%mrro           /hoursOfYear                  ! Kg m-2, soil moisture in each soil layer
        outVars_y%tsl             = outVars_y%tsl            + outVars_h%tsl            /hoursOfYear                   ! K, soil temperature in each soil layer
        outVars_y%tsland          = outVars_y%tsland         + outVars_h%tsland         /hoursOfYear                ! K, surface temperature
        outVars_y%wtd             = outVars_y%wtd            + outVars_h%wtd            /hoursOfYear                   ! m, Water table depth
        outVars_y%snd             = outVars_y%snd            + outVars_h%snd            /hoursOfYear                   ! m, Total snow depth
        outVars_y%lai             = outVars_y%lai            + outVars_h%lai            /hoursOfYear
        ! test_gpp_y        = test_gpp_y       + test_gpp         /hoursOfYear

        ! ! not used in SPRUCE-MIP
        ! rain_yr   = rain_yr   + rain
        ! transp_yr = transp_yr + transp
        ! evap_yr   = evap_yr   + evap
        ! runoff_yr = runoff_yr + runoff
    end subroutine updateYearly

    subroutine summaryHourly(iTotHourly)
        implicit NONE
        integer iTotHourly
        ! summary in the total results
        tot_outVars_h%gpp(iTotHourly)            = outVars_h%gpp         
        tot_outVars_h%npp(iTotHourly)            = outVars_h%npp
        tot_outVars_h%nppLeaf(iTotHourly)        = outVars_h%nppLeaf
        tot_outVars_h%nppWood(iTotHourly)        = outVars_h%nppWood
        tot_outVars_h%nppStem(iTotHourly)        = outVars_h%nppStem
        tot_outVars_h%nppRoot(iTotHourly)        = outVars_h%nppRoot
        tot_outVars_h%nppOther(iTotHourly)       = outVars_h%nppOther                  ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        tot_outVars_h%ra(iTotHourly)             = outVars_h%ra
        tot_outVars_h%raLeaf(iTotHourly)         = outVars_h%raLeaf
        tot_outVars_h%raStem(iTotHourly)         = outVars_h%raStem
        tot_outVars_h%raRoot(iTotHourly)         = outVars_h%raRoot
        tot_outVars_h%raOther(iTotHourly)        = outVars_h%raOther
        tot_outVars_h%rMaint(iTotHourly)         = outVars_h%rMaint
        tot_outVars_h%rGrowth(iTotHourly)        = outVars_h%rGrowth                                             ! maintenance respiration and growth respiration
        tot_outVars_h%rh(iTotHourly)             = outVars_h%rh
        tot_outVars_h%nbp(iTotHourly)            = outVars_h%nbp                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        tot_outVars_h%wetlandCH4(iTotHourly)     = outVars_h%wetlandCH4 
        tot_outVars_h%wetlandCH4prod(iTotHourly) = outVars_h%wetlandCH4prod 
        tot_outVars_h%wetlandCH4cons(iTotHourly) = outVars_h%wetlandCH4cons                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        ! Carbon Pools  (KgC m-2)
        tot_outVars_h%cLeaf(iTotHourly)          = outVars_h%cLeaf
        tot_outVars_h%cStem(iTotHourly)          = outVars_h%cStem
        tot_outVars_h%cRoot(iTotHourly)          = outVars_h%cRoot
        tot_outVars_h%cOther(iTotHourly)         = outVars_h%cOther                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        tot_outVars_h%cLitter(iTotHourly)        = outVars_h%cLitter
        tot_outVars_h%cLitterCwd(iTotHourly)     = outVars_h%cLitterCwd                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        tot_outVars_h%cSoil(iTotHourly)          = outVars_h%cSoil
        tot_outVars_h%cSoilLevels(iTotHourly,:)  = outVars_h%cSoilLevels
        tot_outVars_h%cSoilFast(iTotHourly)      = outVars_h%cSoilFast
        tot_outVars_h%cSoilSlow(iTotHourly)      = outVars_h%cSoilSlow
        tot_outVars_h%cSoilPassive(iTotHourly)   = outVars_h%cSoilPassive                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        tot_outVars_h%CH4(iTotHourly,:)         = outVars_h%CH4                                                        ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        tot_outVars_h%fBNF(iTotHourly)           = outVars_h%fBNF
        tot_outVars_h%fN2O(iTotHourly)           = outVars_h%fN2O
        tot_outVars_h%fNloss(iTotHourly)         = outVars_h%fNloss
        tot_outVars_h%fNnetmin(iTotHourly)       = outVars_h%fNnetmin
        tot_outVars_h%fNdep(iTotHourly)          = outVars_h%fNdep                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        tot_outVars_h%nLeaf(iTotHourly)          = outVars_h%nLeaf
        tot_outVars_h%nStem(iTotHourly)          = outVars_h%nStem
        tot_outVars_h%nRoot(iTotHourly)          = outVars_h%nRoot
        tot_outVars_h%nOther(iTotHourly)         = outVars_h%nOther
        tot_outVars_h%nLitter(iTotHourly)        = outVars_h%nLitter
        tot_outVars_h%nLitterCwd(iTotHourly)     = outVars_h%nLitterCwd
        tot_outVars_h%nSoil(iTotHourly)          = outVars_h%nSoil
        tot_outVars_h%nMineral(iTotHourly)       = outVars_h%nMineral                    ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        tot_outVars_h%hfls(iTotHourly)           = outVars_h%hfls
        tot_outVars_h%hfss(iTotHourly)           = outVars_h%hfss
        tot_outVars_h%SWnet(iTotHourly)          = outVars_h%SWnet
        tot_outVars_h%LWnet(iTotHourly)          = outVars_h%LWnet                               ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        tot_outVars_h%ec(iTotHourly)             = outVars_h%ec
        tot_outVars_h%tran(iTotHourly)           = outVars_h%tran
        tot_outVars_h%es(iTotHourly)             = outVars_h%es                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        tot_outVars_h%hfsbl(iTotHourly)          = outVars_h%hfsbl                                                         ! Snow sublimation
        tot_outVars_h%mrro(iTotHourly)           = outVars_h%mrro
        tot_outVars_h%mrros(iTotHourly)          = outVars_h%mrros
        tot_outVars_h%mrrob(iTotHourly)          = outVars_h%mrrob                                        ! Total runoff; Surface runoff; Subsurface runoff
        ! Other
        tot_outVars_h%mrso(iTotHourly,:)         = outVars_h%mrso                                                   ! Kg m-2, soil moisture in each soil layer
        tot_outVars_h%tsl(iTotHourly,:)          = outVars_h%tsl                                             ! K, soil temperature in each soil layer
        tot_outVars_h%tsland(iTotHourly)         = outVars_h%tsland                                                 ! K, surface temperature
        tot_outVars_h%wtd(iTotHourly)            = outVars_h%wtd                                                 ! m, Water table depth
        tot_outVars_h%snd(iTotHourly)            = outVars_h%snd                                                 ! m, Total snow depth
        tot_outVars_h%lai(iTotHourly)            = outVars_h%lai                     ! m2 m-2, Leaf area index
        ! tot_outVars_h%gdd5(iTotHourly)           = GDD5
        ! tot_outVars_h%onset(iTotHourly)           = onset
        ! tot_outVars_h%storage(iTotHourly)        = storage
        ! tot_outVars_h%add(iTotHourly)            = add
        ! tot_outVars_h%accumulation(iTotHourly)   = accumulation
        ! all_test_h(5*(iTotHourly-1)+1:5*(iTotHourly),:) = test_gpp
        ! all_test_h(iTotHourly,:) = test_gpp

        iTotHourly = iTotHourly+1
    end subroutine summaryHourly
       
    subroutine summaryDaily(iTotDaily)
        implicit none
        integer iTotDaily
        ! Summary 
        ! daily: 
        ! ---------------------------------------------------------------------
        ! carbon fluxes (Kg C m-2 s-1)
        tot_outVars_d%gpp(iTotDaily)            = outVars_d%gpp
        tot_outVars_d%npp(iTotDaily)            = outVars_d%npp
        tot_outVars_d%nppLeaf(iTotDaily)        = outVars_d%nppLeaf
        tot_outVars_d%nppWood(iTotDaily)        = outVars_d%nppWood
        tot_outVars_d%nppStem(iTotDaily)        = outVars_d%nppStem
        tot_outVars_d%nppRoot(iTotDaily)        = outVars_d%nppRoot
        tot_outVars_d%nppOther(iTotDaily)       = outVars_d%nppOther   ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        tot_outVars_d%ra(iTotDaily)             = outVars_d%ra
        tot_outVars_d%raLeaf(iTotDaily)         = outVars_d%raLeaf
        tot_outVars_d%raStem(iTotDaily)         = outVars_d%raStem
        tot_outVars_d%raRoot(iTotDaily)         = outVars_d%raRoot
        tot_outVars_d%raOther(iTotDaily)        = outVars_d%raOther
        tot_outVars_d%rMaint(iTotDaily)         = outVars_d%rMaint
        tot_outVars_d%rGrowth(iTotDaily)        = outVars_d%Rgrowth                                            ! maintenance respiration and growth respiration
        tot_outVars_d%rh(iTotDaily)             = outVars_d%rh
        tot_outVars_d%nbp(iTotDaily)            = outVars_d%nbp                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        tot_outVars_d%wetlandCH4(iTotDaily)     = outVars_d%wetlandCH4
        tot_outVars_d%wetlandCH4prod(iTotDaily) = outVars_d%wetlandCH4prod
        tot_outVars_d%wetlandCH4cons(iTotDaily) = outVars_d%wetlandCH4cons                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        ! Carbon Pools  (KgC m-2)
        tot_outVars_d%cLeaf(iTotDaily)          = outVars_d%cLeaf
        tot_outVars_d%cStem(iTotDaily)          = outVars_d%cStem
        tot_outVars_d%cRoot(iTotDaily)          = outVars_d%cRoot
        tot_outVars_d%cOther(iTotDaily)         = outVars_d%cOther                           ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        tot_outVars_d%cLitter(iTotDaily)        = outVars_d%cLitter
        tot_outVars_d%cLitterCwd(iTotDaily)     = outVars_d%cLitterCwd                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        tot_outVars_d%cSoil(iTotDaily)          = outVars_d%cSoil
        tot_outVars_d%cSoilLevels(iTotDaily,:)    = outVars_d%cSoilLevels
        tot_outVars_d%cSoilFast(iTotDaily)      = outVars_d%cSoilFast
        tot_outVars_d%cSoilSlow(iTotDaily)      = outVars_d%cSoilSlow
        tot_outVars_d%cSoilPassive(iTotDaily)   = outVars_d%cSoilPassive                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        tot_outVars_d%CH4(iTotDaily,:)         = outVars_d%CH4                                                         ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        tot_outVars_d%fBNF(iTotDaily)           = outVars_d%fBNF
        tot_outVars_d%fN2O(iTotDaily)           = outVars_d%fN2O
        tot_outVars_d%fNloss(iTotDaily)         = outVars_d%fNloss
        tot_outVars_d%fNnetmin(iTotDaily)       = outVars_d%fNnetmin
        tot_outVars_d%fNdep(iTotDaily)          = outVars_d%fNdep                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        tot_outVars_d%nLeaf(iTotDaily)          = outVars_d%nLeaf
        tot_outVars_d%nStem(iTotDaily)          = outVars_d%nStem
        tot_outVars_d%nRoot(iTotDaily)          = outVars_d%nRoot
        tot_outVars_d%nOther(iTotDaily)         = outVars_d%nOther
        tot_outVars_d%nLitter(iTotDaily)        = outVars_d%nLitter
        tot_outVars_d%nLitterCwd(iTotDaily)     = outVars_d%nLitterCwd
        tot_outVars_d%nSoil(iTotDaily)          = outVars_d%nSoil
        tot_outVars_d%nMineral(iTotDaily)       = outVars_d%nMineral                    ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        tot_outVars_d%hfls(iTotDaily)           = outVars_d%hfls
        tot_outVars_d%hfss(iTotDaily)           = outVars_d%hfss
        tot_outVars_d%SWnet(iTotDaily)          = outVars_d%SWnet
        tot_outVars_d%LWnet(iTotDaily)          = outVars_d%LWnet                              ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        tot_outVars_d%ec(iTotDaily)             = outVars_d%ec
        tot_outVars_d%tran(iTotDaily)           = outVars_d%tran
        tot_outVars_d%es(iTotDaily)             = outVars_d%es                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        tot_outVars_d%hfsbl(iTotDaily)          = outVars_d%hfsbl                                                         ! Snow sublimation
        tot_outVars_d%mrro(iTotDaily)           = outVars_d%mrro
        tot_outVars_d%mrros(iTotDaily)          = outVars_d%mrros
        tot_outVars_d%mrrob(iTotDaily)          = outVars_d%mrrob                                        ! Total runoff; Surface runoff; Subsurface runoff
        ! Other
        tot_outVars_d%mrso(iTotDaily,:)         = outVars_d%mrso                                          ! Kg m-2, soil moisture in each soil layer
        tot_outVars_d%tsl(iTotDaily,:)          = outVars_d%tsl                                          ! K, soil temperature in each soil layer
        tot_outVars_d%tsland(iTotDaily)         = outVars_d%tsland                                                 ! K, surface temperature
        tot_outVars_d%wtd(iTotDaily)            = outVars_d%wtd                                                 ! m, Water table depth
        tot_outVars_d%snd(iTotDaily)            = outVars_d%snd                                                 ! m, Total snow depth
        tot_outVars_d%lai(iTotDaily)            = outVars_d%lai                                                 ! m2 m-2, Leaf area index

        iTotDaily = iTotDaily+1
    end subroutine summaryDaily

    subroutine summaryMonthly(iTotMonthly)
        implicit none
        integer iTotMonthly
        ! monthly
        ! carbon fluxes (Kg C m-2 s-1)
        tot_outVars_m%gpp(iTotMonthly)            = outVars_m%gpp
        tot_outVars_m%npp(iTotMonthly)            = outVars_m%npp
        tot_outVars_m%nppLeaf(iTotMonthly)        = outVars_m%nppLeaf
        tot_outVars_m%nppWood(iTotMonthly)        = outVars_m%nppWood
        tot_outVars_m%nppStem(iTotMonthly)        = outVars_m%nppStem
        tot_outVars_m%nppRoot(iTotMonthly)        = outVars_m%nppRoot
        tot_outVars_m%nppOther(iTotMonthly)       = outVars_m%nppOther   ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        tot_outVars_m%ra(iTotMonthly)             = outVars_m%ra
        tot_outVars_m%raLeaf(iTotMonthly)         = outVars_m%raLeaf
        tot_outVars_m%raStem(iTotMonthly)         = outVars_m%raStem
        tot_outVars_m%raRoot(iTotMonthly)         = outVars_m%raRoot
        tot_outVars_m%raOther(iTotMonthly)        = outVars_m%raOther
        tot_outVars_m%rMaint(iTotMonthly)         = outVars_m%rMaint
        tot_outVars_m%rGrowth(iTotMonthly)        = outVars_m%rGrowth                                             ! maintenance respiration and growth respiration
        tot_outVars_m%rh(iTotMonthly)             = outVars_m%rh
        tot_outVars_m%nbp(iTotMonthly)            = outVars_m%nbp                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        tot_outVars_m%wetlandCH4(iTotMonthly)     = outVars_m%wetlandCH4
        tot_outVars_m%wetlandCH4prod(iTotMonthly) = outVars_m%wetlandCH4prod
        tot_outVars_m%wetlandCH4cons(iTotMonthly) = outVars_m%wetlandCH4cons                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        ! Carbon Pools  (KgC m-2)
        tot_outVars_m%cLeaf(iTotMonthly)          = outVars_m%cLeaf
        tot_outVars_m%cStem(iTotMonthly)          = outVars_m%cStem
        tot_outVars_m%cRoot(iTotMonthly)          = outVars_m%cRoot
        tot_outVars_m%cOther(iTotMonthly)         = outVars_m%cOther                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        tot_outVars_m%cLitter(iTotMonthly)        = outVars_m%cLitter
        tot_outVars_m%cLitterCwd(iTotMonthly)     = outVars_m%cLitterCwd                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        tot_outVars_m%cSoil(iTotMonthly)          = outVars_m%cSoil
        tot_outVars_m%cSoilLevels(iTotMonthly,:)  = outVars_m%cSoilLevels
        tot_outVars_m%cSoilFast(iTotMonthly)      = outVars_m%cSoilFast
        tot_outVars_m%cSoilSlow(iTotMonthly)      = outVars_m%cSoilSlow
        tot_outVars_m%cSoilPassive(iTotMonthly)   = outVars_m%cSoilPassive                           ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        tot_outVars_m%CH4(iTotMonthly,:)         = outVars_m%CH4                                                      ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        tot_outVars_m%fBNF(iTotMonthly)           = outVars_m%fBNF
        tot_outVars_m%fN2O(iTotMonthly)           = outVars_m%fN2O
        tot_outVars_m%fNloss(iTotMonthly)         = outVars_m%fNloss
        tot_outVars_m%fNnetmin(iTotMonthly)       = outVars_m%fNnetmin
        tot_outVars_m%fNdep(iTotMonthly)          = outVars_m%fNdep                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        tot_outVars_m%nLeaf(iTotMonthly)          = outVars_m%nLeaf
        tot_outVars_m%nStem(iTotMonthly)          = outVars_m%nStem
        tot_outVars_m%nRoot(iTotMonthly)          = outVars_m%nRoot
        tot_outVars_m%nOther(iTotMonthly)         = outVars_m%nOther
        tot_outVars_m%nLitter(iTotMonthly)        = outVars_m%nLitter
        tot_outVars_m%nLitterCwd(iTotMonthly)     = outVars_m%nLitterCwd
        tot_outVars_m%nSoil(iTotMonthly)          = outVars_m%nSoil
        tot_outVars_m%nMineral(iTotMonthly)       = outVars_m%nMineral                    ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        tot_outVars_m%hfls(iTotMonthly)           = outVars_m%hfls
        tot_outVars_m%hfss(iTotMonthly)           = outVars_m%hfss
        tot_outVars_m%SWnet(iTotMonthly)          = outVars_m%SWnet
        tot_outVars_m%LWnet(iTotMonthly)          = outVars_m%LWnet                                ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        tot_outVars_m%ec(iTotMonthly)             = outVars_m%ec
        tot_outVars_m%tran(iTotMonthly)           = outVars_m%tran
        tot_outVars_m%es(iTotMonthly)             = outVars_m%es                                             ! Canopy evaporation; Canopy transpiration; Soil evaporation
        tot_outVars_m%hfsbl(iTotMonthly)          = outVars_m%hfsbl                                                         ! Snow sublimation
        tot_outVars_m%mrro(iTotMonthly)           = outVars_m%mrro
        tot_outVars_m%mrros(iTotMonthly)          = outVars_m%mrros
        tot_outVars_m%mrrob(iTotMonthly)          = outVars_m%mrrob 
        ! Other
        tot_outVars_m%mrso(iTotMonthly,:)         = outVars_m%mrso                                         ! Kg m-2, soil moisture in each soil layer
        tot_outVars_m%tsl(iTotMonthly,:)          = outVars_m%tsl                                         ! K, soil temperature in each soil layer
        tot_outVars_m%tsland(iTotMonthly)         = outVars_m%tsland                                                 ! K, surface temperature
        tot_outVars_m%wtd(iTotMonthly)            = outVars_m%wtd                                                 ! m, Water table depth
        tot_outVars_m%snd(iTotMonthly)            = outVars_m%snd                                                 ! m, Total snow depth
        tot_outVars_m%lai(iTotMonthly)            = outVars_m%lai                                                 ! m2 m-2, Leaf area index

        iTotMonthly = iTotMonthly + 1
    end subroutine summaryMonthly

    subroutine summaryYearly(iTotYearly)
        implicit none
        integer iTotYearly
        ! monthly
        ! carbon fluxes (Kg C m-2 s-1)
        tot_outVars_y%gpp(iTotYearly)            = outVars_y%gpp
        tot_outVars_y%npp(iTotYearly)            = outVars_y%npp
        tot_outVars_y%nppLeaf(iTotYearly)        = outVars_y%nppLeaf
        tot_outVars_y%nppWood(iTotYearly)        = outVars_y%nppWood
        tot_outVars_y%nppStem(iTotYearly)        = outVars_y%nppStem
        tot_outVars_y%nppRoot(iTotYearly)        = outVars_y%nppRoot
        tot_outVars_y%nppOther(iTotYearly)       = outVars_y%nppOther   ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        tot_outVars_y%ra(iTotYearly)             = outVars_y%ra
        tot_outVars_y%raLeaf(iTotYearly)         = outVars_y%raLeaf
        tot_outVars_y%raStem(iTotYearly)         = outVars_y%raStem
        tot_outVars_y%raRoot(iTotYearly)         = outVars_y%raRoot
        tot_outVars_y%raOther(iTotYearly)        = outVars_y%raOther
        tot_outVars_y%rMaint(iTotYearly)         = outVars_y%rMaint
        tot_outVars_y%rGrowth(iTotYearly)        = outVars_y%rGrowth                                             ! maintenance respiration and growth respiration
        tot_outVars_y%rh(iTotYearly)             = outVars_y%rh
        tot_outVars_y%nbp(iTotYearly)            = outVars_y%nbp                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        tot_outVars_y%wetlandCH4(iTotYearly)     = outVars_y%wetlandCH4
        tot_outVars_y%wetlandCH4prod(iTotYearly) = outVars_y%wetlandCH4prod
        tot_outVars_y%wetlandCH4cons(iTotYearly) = outVars_y%wetlandCH4cons                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        ! Carbon Pools  (KgC m-2)
        tot_outVars_y%cLeaf(iTotYearly)          = outVars_y%cLeaf
        tot_outVars_y%cStem(iTotYearly)          = outVars_y%cStem
        tot_outVars_y%cRoot(iTotYearly)          = outVars_y%cRoot
        tot_outVars_y%cOther(iTotYearly)         = outVars_y%cOther                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        tot_outVars_y%cLitter(iTotYearly)        = outVars_y%cLitter
        tot_outVars_y%cLitterCwd(iTotYearly)     = outVars_y%cLitterCwd                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        tot_outVars_y%cSoil(iTotYearly)          = outVars_y%cSoil
        tot_outVars_y%cSoilLevels(iTotYearly,:)  = outVars_y%cSoilLevels
        tot_outVars_y%cSoilFast(iTotYearly)      = outVars_y%cSoilFast
        tot_outVars_y%cSoilSlow(iTotYearly)      = outVars_y%cSoilSlow
        tot_outVars_y%cSoilPassive(iTotYearly)   = outVars_y%cSoilPassive                           ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        tot_outVars_y%CH4(iTotYearly,:)         = outVars_y%CH4                                                      ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        tot_outVars_y%fBNF(iTotYearly)           = outVars_y%fBNF
        tot_outVars_y%fN2O(iTotYearly)           = outVars_y%fN2O
        tot_outVars_y%fNloss(iTotYearly)         = outVars_y%fNloss
        tot_outVars_y%fNnetmin(iTotYearly)       = outVars_y%fNnetmin
        tot_outVars_y%fNdep(iTotYearly)          = outVars_y%fNdep                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        tot_outVars_y%nLeaf(iTotYearly)          = outVars_y%nLeaf
        tot_outVars_y%nStem(iTotYearly)          = outVars_y%nStem
        tot_outVars_y%nRoot(iTotYearly)          = outVars_y%nRoot
        tot_outVars_y%nOther(iTotYearly)         = outVars_y%nOther
        tot_outVars_y%nLitter(iTotYearly)        = outVars_y%nLitter
        tot_outVars_y%nLitterCwd(iTotYearly)     = outVars_y%nLitterCwd
        tot_outVars_y%nSoil(iTotYearly)          = outVars_y%nSoil
        tot_outVars_y%nMineral(iTotYearly)       = outVars_y%nMineral                    ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        tot_outVars_y%hfls(iTotYearly)           = outVars_y%hfls
        tot_outVars_y%hfss(iTotYearly)           = outVars_y%hfss
        tot_outVars_y%SWnet(iTotYearly)          = outVars_y%SWnet
        tot_outVars_y%LWnet(iTotYearly)          = outVars_y%LWnet                                ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        tot_outVars_y%ec(iTotYearly)             = outVars_y%ec
        tot_outVars_y%tran(iTotYearly)           = outVars_y%tran
        tot_outVars_y%es(iTotYearly)             = outVars_y%es                                             ! Canopy evaporation; Canopy transpiration; Soil evaporation
        tot_outVars_y%hfsbl(iTotYearly)          = outVars_y%hfsbl                                                         ! Snow sublimation
        tot_outVars_y%mrro(iTotYearly)           = outVars_y%mrro
        tot_outVars_y%mrros(iTotYearly)          = outVars_y%mrros
        tot_outVars_y%mrrob(iTotYearly)          = outVars_y%mrrob 
        ! Other
        tot_outVars_y%mrso(iTotYearly,:)         = outVars_y%mrso                                         ! Kg m-2, soil moisture in each soil layer
        tot_outVars_y%tsl(iTotYearly,:)          = outVars_y%tsl                                         ! K, soil temperature in each soil layer
        tot_outVars_y%tsland(iTotYearly)         = outVars_y%tsland                                                 ! K, surface temperature
        tot_outVars_y%wtd(iTotYearly)            = outVars_y%wtd                                                 ! m, Water table depth
        tot_outVars_y%snd(iTotYearly)            = outVars_y%snd                                                 ! m, Total snow depth
        tot_outVars_y%lai(iTotYearly)            = outVars_y%lai                                                 ! m2 m-2, Leaf area index

        iTotYearly = iTotYearly + 1
    end subroutine summaryYearly
end module update_and_summary