module datatypes
    implicit none
    ! settings
    character(100) :: teco_configfile  ! the TECO cofig file
    character(50)  :: case_name       ! define the case name
    logical :: do_simu                ! simulation mode
    logical :: do_mcmc                ! MCMC for data assimilation mode
    logical :: do_spinup              ! spinup mode
    logical :: do_matrix              ! whether run matrix or not
    logical :: do_restart             ! whether read restart file or not
    ! simulation selections
    logical :: do_snow                ! do soil snow process or not
    logical :: do_soilphy             ! do soil physics or not
    logical :: do_EBG                 ! run EBG or not based on Ma et al., 2022
    logical :: do_ndep                ! N deposit
    logical :: do_leap                ! judge leap year or not
    ! output selections
    logical :: do_out_hr
    logical :: do_out_day
    logical :: do_out_mon
    logical :: do_out_yr

    integer :: dtimes                 ! 24: hourly simulation
    ! set the input and output path and files
    character(200) :: inDir
    character(200) :: outDir
    ! input files
    character(300) :: climfile
    character(300) :: watertablefile
    character(300) :: snowdepthfile
    character(300) :: in_restartfile
    character(300) :: mcmc_configfile
    character(300) :: spinup_configfile
    ! above settings from the nml file
    ! output path
    character(250) :: outDir_nc     = "results_nc_format"
    character(250) :: outDir_csv    = "results_csv_format"
    character(250) :: outDir_mcmc   = "results_mcmc"
    character(250) :: outDir_spinup = "results_spinup"

    ! experiment settings
    real :: Ttreat   = 0.        ! Temperature treatment, warming in air and soil temperature
    real :: CO2treat = 0.        ! CO2 treatmant, up to CO2treat, not add to Ca. CO2
    real :: N_fert   = 0.        ! 5.6 ! (11.2 gN m-2 yr-1, in spring, Duke Forest FACE)
    ! ---------------------------------------------------------------------------------------
    
    integer, parameter :: nlayers = 10                ! how many
    real,    parameter :: pi      = 3.1415926
    ! physical constants
    real,    parameter :: tauL(3) = (/0.1, 0.425, 0.00/)  ! leaf transmittance for vis, for NIR, for thermal
    real,    parameter :: rhoL(3) = (/0.1, 0.425, 0.00/)  ! leaf reflectance for vis, for NIR, for thermal
    real,    parameter :: emleaf  = 0.96
    real,    parameter :: emsoil  = 0.94
    real,    parameter :: Rconst  = 8.314                 ! universal gas constant (J/mol)
    real,    parameter :: sigma   = 5.67e-8               ! Steffan Boltzman constant (W/m2/K4)
    real,    parameter :: cpair   = 1010.                 ! heat capapcity of air (J/kg/K)
    real,    parameter :: Patm    = 101325. !1.e5         ! atmospheric pressure  (Pa)
    real,    parameter :: Trefk   = 293.2                 ! reference temp K for Kc, Ko, Rd
    real,    parameter :: H2OLv0  = 2.501e6               ! latent heat H2O (J/kg)
    real,    parameter :: AirMa   = 29.e-3                ! mol mass air (kg/mol)
    real,    parameter :: H2OMw   = 18.e-3                ! mol mass H2O (kg/mol)
    real,    parameter :: chi     = 0.93                  ! gbH/gbw
    real,    parameter :: Dheat   = 21.5e-6               ! molecular diffusivity for heat
    ! plant parameters
    real,    parameter :: gsw0    = 1.0e-2                ! g0 for H2O in BWB model
    real,    parameter :: theta   = 0.9
    real,    parameter :: wleaf   = 0.01                  ! leaf width (m)
    ! thermodynamic parameters for Kc and Ko (Leuning 1990)
    real,    parameter :: conKc0  = 302.e-6               ! mol mol^-1
    real,    parameter :: conKo0  = 256.e-3               ! mol mol^-1
    real,    parameter :: Ekc     = 59430.                ! J mol^-1
    real,    parameter :: Eko     = 36000.                ! J mol^-1
    ! Erd = 53000.                                        ! J mol^-1
    real,    parameter :: o2ci    = 210.e-3               ! mol mol^-1
    ! thermodynamic parameters for Vcmax & Jmax (Eq 9, Harley et al, 1992; #1392)
    real,    parameter :: Eavm    = 116300.               ! J/mol  (activation energy)
    real,    parameter :: Edvm    = 202900.               ! J/mol  (deactivation energy)
    real,    parameter :: Eajm    = 79500.                ! J/mol  (activation energy) 
    real,    parameter :: Edjm    = 201000.               ! J/mol  (deactivation energy)
    ! parameters for temperature dependence of gamma* (revised from von Caemmerer et al 1993)
    real,    parameter :: gam0    = 28.0e-6               ! mol mol^-1 @ 20C = 36.9 @ 25C
    real,    parameter :: gam1    = .0509
    real,    parameter :: gam2    = .0010
    real,    parameter :: times_storage_use=3*720.        ! 720 hours, 30 days
    real :: rhoS(3) = (/0.1, 0.3,   0.00/)                ! soil reflectance for vis, for NIR, for thermal, update in vegetation?
    ! end of consts parameters -------------------------------------------------------------------

    ! climate data type
    type forcing_data_type
        integer :: year
        integer :: doy
        integer :: hour
        real    :: Tair
        real    :: Tsoil
        real    :: RH                   ! Jian: RH seems confused in forcing and soil respiration
        real    :: VPD
        real    :: Rain
        real    :: WS
        real    :: PAR
        real    :: CO2
        real    :: PBOT                 ! unit patm Pa dynamic atmosphere pressure
        real    :: Ndep
    end type forcing_data_type
    type(forcing_data_type), allocatable, save :: forcing(:)
    real :: co2ca

    type site_data_type   ! data use in this model, but the common parameters in this site
        integer:: lat, lon
        real :: wsmax
        real :: wsmin
        real :: extkU
        real :: G
        real :: GDD5
        real :: rdepth
        real :: Q10pro
        real :: r_me
        real :: Toxi
        real :: Omax
        real :: kCH4
        real :: CH4_thre
        real :: Tveg
        real :: f
        real :: bubprob
        real :: Vmaxfraction  
        ! some environmental variables
        real :: Dair, raero, Tsoill(10), Tsnow, depth_ex, ta, rain_d
        real :: tsoil_layer(11)
        real :: Twater, Tice, snow_dsim, dcount, sublim, fsub
        real :: resht, fa, decay_m, rho_snow
        real :: dpatm
        ! state variables
        real :: depth(nlayers), THKSL(nlayers), FRLEN(10)
        real :: wcl(10), ice(10), wsc(10), liq_water(10), zwt
        real :: fwsoil, topfws, omega, scalW
        real :: Rsoilab1, Rsoilab2, Rsoilab3, Rsoilabs
        real :: melt, infilt, runoff
        real :: evap
        real :: diff_snow, condu_snow, shcap_snow, diff_s, condu_b
        real :: water_tw, ice_tw, snow_depth, albedo_snow
        real :: thd_snow_depth, dcount_soil, sftmp
        ! energy
        real :: Esoil, Hsoil
        ! soil flux
        real :: Rhetero, Rh_pools(5) 
        ! methane
        real :: simuCH4, Pro_sum, CH4(nlayers), Oxi_sum, CH4_V(nlayers)
        real :: Tpro_me, pwater(nlayers), presP(nlayers), Vp(nlayers) 
        real :: methanebP(nlayers), methaneP(nlayers)
        real :: bubble_methane_tot, Nbub
    end type site_data_type
    type(site_data_type) :: st

    ! pft data type for different species
    type spec_data_type
        real :: LAI
        ! special paramaters for different species
        real :: LAImin
        real :: LAImax
        real :: SLA
        real :: xfang
        real :: Vcmx0 
        real :: Vcmax0    ! fixed value from namelist
        real :: eJmx0
        real :: gddonset
        real :: stom_n
        real :: alpha
        real :: Ds0
        real :: Rl0
        real :: Rs0
        real :: Rr0
        real :: Q10
        real :: SapS, SapR
        real :: hmax            ! in plant growth hmax = 24.19   ! m
        real :: hl0             ! in plant growth hl0  = 0.00019  ! m2/kg C
        real :: LAIMAX0         ! in plant growth LAIMAX0 = 8.    ! maybe the LAImax
        real :: la0             ! in plant growht la0     = 0.2
        real :: GLmax, GSmax, GRmax
        ! flux 
        real :: gpp, npp
        real :: transp
        real :: evap
        real :: RmLeaf, RmStem, RmRoot, Rmain
        real :: Rgrowth
        ! states
        real :: QC1,  QC2,  QC3,  QC4,  QC5             ! leaf, stem, root, Litm, Lits
        real :: CN01, CN02, CN03, CN04, CN05
        real :: bmleaf, bmstem, bmroot, bmplant 
        real :: StemSap,RootSap, NSC, NSCmax, add
        real :: storage, stor_use, accumulation, store
        real :: L_fall
        ! special cycle 
        real :: tauCLeaf, tauCStem, tauCRoot
        ! N scalar
        real :: SNvcmax, SNgrowth, SNRauto
        real :: fnsc, NSN
        ! 
        integer :: onset
        real :: alpha_L, alpha_W, alpha_R
        ! energy
        real :: QLleaf
    end type spec_data_type

    type vegn_tile_type
        integer :: npft
        type(spec_data_type), allocatable :: allSp(:)
        real :: LAI
        real :: LAImin
        real :: LAImax
        ! total flux
        real :: gpp
        real :: transp
        real :: evap
        real :: RmLeaf, RmStem, RmRoot, Rmain
        ! states
        real :: bmleaf, bmstem, bmroot, bmplant
        real :: NSC
        ! energy
        real :: QLleaf, Rsoilab1, Rsoilab2
    end type vegn_tile_type

    ! some parameters, may be in species, may be in cycle

    ! outputs
    type spec_outvars_type
        ! carbon fluxes (Kg C m-2 s-1)
        real, allocatable :: gpp(:)
        real, allocatable :: nee(:)
        real, allocatable :: npp(:)
        real, allocatable :: nppLeaf(:)
        real, allocatable :: nppWood(:)
        real, allocatable :: nppStem(:)
        real, allocatable :: nppRoot(:)
        real, allocatable :: nppOther(:)    ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        real, allocatable :: ra(:)
        real, allocatable :: raLeaf(:)
        real, allocatable :: raStem(:)
        real, allocatable :: raRoot(:)
        real, allocatable :: raOther(:)
        real, allocatable :: rMaint(:)
        real, allocatable :: rGrowth(:)
        real, allocatable :: nbp(:)
        ! Carbon Pools  (KgC m-2)
        real, allocatable :: cLeaf(:)
        real, allocatable :: cStem(:)
        real, allocatable :: cRoot(:)
        real, allocatable :: cOther(:)              ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        real, allocatable :: cLitter(:)
        real, allocatable :: cLitterCwd(:)
        ! Nitrogen pools (kgN m-2)
        real, allocatable :: nLeaf(:)
        real, allocatable :: nStem(:)
        real, allocatable :: nRoot(:)
        real, allocatable :: nOther(:)
        real, allocatable :: nLitter(:)
        real, allocatable :: nLitterCwd(:)
        ! water fluxes (kg m-2 s-1)
        real, allocatable :: ec(:)
        real, allocatable :: tran(:)
        real, allocatable :: es(:)   
        ! other
        real, allocatable :: lai(:)                     ! m2 m-2, Leaf area index
    end type spec_outvars_type

    ! total outputs
    type outvars_data_type
        type(spec_data_type), allocatable :: allSpec(:)
        ! carbon fluxes (Kg C m-2 s-1)
        real, allocatable :: gpp(:)
        real, allocatable :: nee(:)
        real, allocatable :: npp(:)
        real, allocatable :: nppLeaf(:)
        real, allocatable :: nppWood(:)
        real, allocatable :: nppStem(:)
        real, allocatable :: nppRoot(:)
        real, allocatable :: nppOther(:)           ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        real, allocatable :: ra(:)
        real, allocatable :: raLeaf(:)
        real, allocatable :: raStem(:)
        real, allocatable :: raRoot(:)
        real, allocatable :: raOther(:)
        real, allocatable :: rMaint(:)
        real, allocatable :: rGrowth(:)            ! maintenance respiration and growth respiration
        real, allocatable :: rh(:)
        real, allocatable :: nbp(:)                ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        real, allocatable :: wetlandCH4(:)
        real, allocatable :: wetlandCH4prod(:)
        real, allocatable :: wetlandCH4cons(:)     ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        ! Carbon Pools  (KgC m-2)
        real, allocatable :: cLeaf(:)
        real, allocatable :: cStem(:)
        real, allocatable :: cRoot(:)
        real, allocatable :: cOther(:)              ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        real, allocatable :: cLitter(:)
        real, allocatable :: cLitterCwd(:)          ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        real, allocatable :: cSoil(:)
        real, allocatable :: cSoilLevels(:, :)
        real, allocatable :: cSoilFast(:)
        real, allocatable :: cSoilSlow(:)
        real, allocatable :: cSoilPassive(:)        ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        real, allocatable :: CH4(:, :)              ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        real, allocatable :: fBNF(:)
        real, allocatable :: fN2O(:)
        real, allocatable :: fNloss(:)
        real, allocatable :: fNnetmin(:)
        real, allocatable :: fNdep(:)               ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        real, allocatable :: nLeaf(:)
        real, allocatable :: nStem(:)
        real, allocatable :: nRoot(:)
        real, allocatable :: nOther(:)
        real, allocatable :: nLitter(:)
        real, allocatable :: nLitterCwd(:)
        real, allocatable :: nSoil(:)
        real, allocatable :: nMineral(:)                ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        real, allocatable :: hfls(:)
        real, allocatable :: hfss(:)
        real, allocatable :: SWnet(:)
        real, allocatable :: LWnet(:)                   ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        real, allocatable :: ec(:)
        real, allocatable :: tran(:)
        real, allocatable :: es(:)                      ! Canopy evaporation; Canopy transpiration; Soil evaporation
        real, allocatable :: hfsbl(:)                   ! Snow sublimation
        real, allocatable :: mrro(:)
        real, allocatable :: mrros(:)
        real, allocatable :: mrrob(:)                   ! Total runoff; Surface runoff; Subsurface runoff
        ! other
        real, allocatable :: mrso(:, :)           ! Kg m-2, soil moisture in each soil layer
        real, allocatable :: tsl(:, :)            ! K, soil temperature in each soil layer
        real, allocatable :: tsland(:)                  ! K, surface temperature
        real, allocatable :: wtd(:)                     ! m, Water table depth
        real, allocatable :: snd(:)                     ! m, Total snow depth
        real, allocatable :: lai(:)                     ! m2 m-2, Leaf area index 
    end type outvars_data_type

contains
    subroutine read_teco_configs()
        implicit none
        integer io
        namelist /nml_teco_settings/ case_name, do_simu, do_mcmc, do_spinup, do_matrix, &
            do_restart, do_snow, do_soilphy, do_EBG, do_ndep, do_leap, do_out_hr,       &
            do_out_day, do_out_mon, do_out_yr, dtimes, inDir, outDir, climfile,         &
            watertablefile, snowdepthfile, in_restartfile, mcmc_configfile, spinup_configfile
        namelist /nml_exps/ Ttreat, CO2treat, N_fert
        print *, "# read TECO config nml file ..."
        open(388, file = teco_configfile)
        read(388, nml  = nml_teco_settings,  iostat=io)
        read(388, nml  = nml_exps,           iostat=io)
        close(388)
    end subroutine read_teco_configs

    subroutine read_parameters_nml(param_nml_file)
        implicit none
        character(*), intent(in) :: param_nml_file
        integer io
        ! ! site special variables that are read from namelist file
        ! type spec_data_type
        real :: lat, lon
        real :: wsmax, wsmin
        real :: LAImin, LAImax
        real :: SLAx, rdepth
        real :: Rootmax, Stemmax
        real :: SapR, SapS
        real :: GLmax, GRmax, Gsmax
        real :: stom_n, a1, Ds0
        real :: Vcmax0                              ! Jian: Vcmax0 and Vcmx0 is same? Vcmax0 is Vcmx0 in consts
        real :: extkU, xfang, alpha               
        real :: Tau_Leaf, Tau_Wood, Tau_Root        ! turnover rate of plant carbon pools : leaf, wood, root  
        real :: Tau_F, Tau_C                        ! turnover rate of litter carbon pools: fine, coarse 
        real :: Tau_Micro, Tau_slowSOM, Tau_Passive ! turnover rate of soil carbon pools  : fast, slow, passive 
        real :: gddonset
        real :: Q10, Q10rh                          ! Q10rh modified from Ma et al.,2023 for aclimate study, change in transfer module of Q10h
        real :: Rl0, Rs0, Rr0
        ! added for parameters in methane module   
        real :: r_me, Q10pro
        real :: kCH4, Omax
        real :: CH4_thre
        real :: Tveg, Toxi 
        real :: Tpro_me
        ! add based on Ma et al., 2022
        real :: f, bubprob, Vmaxfraction  
        ! add based on Ma et al., 2023
        real :: JV, Entrpy              ! J/mol/K (entropy term, for Jmax & Vcmax)
        real :: etaL, etaW, etaR        ! etaL and etaR are not used.
        real :: f_F2M, f_C2M, f_C2S
        real :: f_M2S, f_M2P, f_S2P
        real :: f_S2M, f_P2M
        ! ----------------------------------------------------
        namelist/nml_params/ lat, lon, wsmax, wsmin, LAImax, LAImin, rdepth,    &
            Rootmax, Stemmax, SapR, SapS, SLAx, GLmax, GRmax, Gsmax, stom_n,    &
            a1, Ds0, Vcmax0, extkU, xfang, alpha, Tau_Leaf, Tau_Wood, Tau_Root, &
            Tau_F, Tau_C, Tau_Micro, Tau_SlowSOM, Tau_Passive, gddonset, Q10,   &
            Q10rh, Rl0, Rs0, Rr0, r_me, Q10pro, kCH4, Omax, CH4_thre, Tveg,     &
            Tpro_me, Toxi, f, bubprob, Vmaxfraction, JV, Entrpy, etaL, etaW,    &
            etaR, f_F2M, f_C2M, f_C2S, f_M2S, f_M2P, f_S2P, f_S2M, f_P2M
        ! ------------------------------------------------------------------------
        ! parameters that are needed to be initilized.
        real :: QC(8), CN0(8)               ! leaf,wood,root,fine lit.,coarse lit.,Micr,Slow,Pass
        real :: NSCmin, Storage, nsc
        real :: accumulation, SNvcmax 
        real :: N_deposit, alphaN, NSN
        real :: QNminer, N_deficit
        real :: thksl(10)                   ! thickness of every soil layer
        real :: FRLEN(10)                   ! ratio of roots in every layer, Oak Ridge FACE: Shuang
        real :: liq_water(10)               ! unit m
        real :: fwsoil, topfws, omega       ! update in soilwater module
        real :: zwt, infilt
        real :: sftmp, Tsnow, Twater
        real :: Tice, G, snow_dsim
        real :: dcount, dcount_soil
        real :: ice_tw   
        real :: Tsoill(10), ice(10)
        real :: shcap_snow                  ! tuneice worker better
        real :: condu_snow, condu_b         ! yuanyuan soil thermal version value  ... int: this par is not sensitive to CWE
        real :: depth_ex, diff_s, diff_snow ! int diffusivity of snow not sensitive for ice
        real :: albedo_snow, resht     
        real :: thd_snow_depth, b_bound     ! tuneice  not sensitive for ice
        real :: infilt_rate, fa, fsub
        real :: rho_snow, decay_m           ! aging factor on snow melting
        ! methane module. update: Shuang methane bog species even more shallowly rooted than the tundra. add initials for methane module Shuang version
        real :: CH4_V(10), CH4(10), Vp(10)  ! assume in the very beginning no bubbles exist in the first three layers (30cm)
        real :: bubble_methane_tot, Nbub
        real :: depth_1                     ! calculate soil depth unit cm
        ! -----------------------------------------------------------------
        namelist/nml_initial_values/ QC, CN0, NSCmin, Storage, nsc, accumulation, SNvcmax, &
            N_deposit, alphaN, NSN, QNminer, N_deficit, thksl, FRLEN, liq_water, fwsoil,   &
            topfws, omega, zwt, infilt, sftmp, Tsnow, Twater, Tice, G, snow_dsim, dcount,  &
            dcount_soil, ice_tw, Tsoill, ice, shcap_snow, condu_snow, condu_b, depth_ex,   & 
            diff_s, diff_snow, albedo_snow, resht, thd_snow_depth, b_bound, infilt_rate,   &
            fa, fsub, rho_snow, decay_m, CH4_V, CH4, Vp, bubble_methane_tot, Nbub, depth_1
        ! -----------------------------------------------------------------------------------
        print *, "# read parameters nml file: ", param_nml_file
        open(343, file = param_nml_file)
        read(343, nml  = nml_params,         iostat=io)
        read(343, nml  = nml_initial_values, iostat=io)
        close(343)
    end subroutine read_parameters_nml

end module datatypes