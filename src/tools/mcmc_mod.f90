module mcmc_mod
    ! some functions that both driver module and MCMC module.
    use datatypes
    implicit none

    ! parameters and observation files

    integer npar, nDAsimu, ncov, nRand
    real search_scale
    logical :: do_mc_out_hr, do_mc_out_day, do_mc_out_mon, do_mc_out_yr

    real mc_lat, mc_Longitude, mc_wsmax, mc_wsmin
    real mc_LAIMAX, mc_LAIMIN, mc_rdepth, mc_Rootmax, mc_Stemmax
    real mc_SapR, mc_SapS, mc_SLA, mc_GLmax, mc_GRmax, mc_Gsmax, mc_stom_n
    real mc_a1, mc_Ds0, mc_Vcmx0, mc_extkU, mc_xfang, mc_alpha
    real mc_Tau_Leaf, mc_Tau_Wood, mc_Tau_Root, mc_Tau_F
    real mc_Tau_C,  mc_Tau_Micro, mc_Tau_SlowSOM, mc_Tau_Passive
    real mc_gddonset, mc_Q10, mc_Rl0, mc_Rs0, mc_Rr0
    real mc_r_me, mc_Q10pro, mc_kCH4, mc_Omax, mc_CH4_thre
    real mc_Tveg, mc_Tpro_me, mc_Toxi
    real mc_f, mc_bubprob, mc_Vmaxfraction
    real mc_Q10rh, mc_JV, mc_Entrpy
    real mc_etaL, mc_etaW, mc_etaR
    real mc_f_F2M, mc_f_C2M, mc_f_C2S, mc_f_M2S
    real mc_f_M2P, mc_f_S2P, mc_f_S2M, mc_f_P2M

    type params_mcmc
        real lat, Longitude, wsmax, wsmin
        real LAIMAX, LAIMIN, rdepth, Rootmax, Stemmax
        real SapR, SapS, SLA, GLmax, GRmax, Gsmax, stom_n
        real a1, Ds0, Vcmx0, extkU, xfang, alpha
        real Tau_Leaf, Tau_Wood, Tau_Root, Tau_F
        real Tau_C,  Tau_Micro, Tau_SlowSOM, Tau_Passive
        real gddonset, Q10, Rl0, Rs0, Rr0
        real r_me, Q10pro, kCH4, Omax, CH4_thre
        real Tveg, Tpro_me, Toxi
        real f, bubprob, Vmaxfraction
        real Q10rh, JV, Entrpy
        real etaL, etaW, etaR
        real f_F2M, f_C2M, f_C2S, f_M2S
        real f_M2P, f_S2P, f_S2M, f_P2M
    end type params_mcmc

    type(params_mcmc), allocatable :: parval(:), parmin(:), parmax(:)

    ! type(nml_params_data_type) :: in_params, in_parval, in_parval_min, in_parval_max
    ! real, allocatable :: parval(:), parmin(:), parmax(:)
    character(20), allocatable :: parnames(:)

    ! observational file path
    character(500) :: obsfile_gpp_d, obsfile_nee_d, obsfile_reco_d
    character(500) :: obsfile_gpp_h, obsfile_nee_h, obsfile_reco_h
    ! methane   
    character(500) :: obsfile_ch4_h
    ! c pools
    character(500) :: obsfile_cleaf, obsfile_cwood 

    character(500) :: obsfile_anpp_y, obsfile_bnpp_y 
    character(500) :: obsfile_lai_h, obsfile_npp_y, obsfile_reco_y 

    ! variables for calculating the cost in MCMC processes
    type interCostVariable
        character(300) :: filepath
        logical :: existOrNot
        real, allocatable :: obsData(:,:)
        real, allocatable :: mdData(:,:)
    end type interCostVariable

    type allCostVariables
    ! default variables, you can add the variable names here. (year, doy, hour, value, std.)
        ! carbon flux 
        type(interCostVariable) :: gpp_d
        type(interCostVariable) :: nee_d
        type(interCostVariable) :: reco_d
        type(interCostVariable) :: gpp_h
        type(interCostVariable) :: nee_h
        type(interCostVariable) :: reco_h
        ! methane
        type(interCostVariable) :: ch4_h
        ! c pools
        type(interCostVariable) :: cleaf    ! foliage
        type(interCostVariable) :: cwood
        ! for different species
        type(interCostVariable) :: anpp_y
        type(interCostVariable) :: bnpp_y
        type(interCostVariable) :: lai_h
        type(interCostVariable) :: npp_y
        type(interCostVariable) :: reco_y
    end type allCostVariables

    type(allCostVariables) :: vars4MCMC      ! define a allCostVariables first

    ! variables for marking the cycle number
    integer mc_itime_gpp_d, mc_itime_nee_d, mc_itime_reco_d
    integer mc_itime_gpp_h, mc_itime_nee_h, mc_itime_reco_h
    integer mc_itime_ch4_h, mc_itime_cleaf, mc_itime_cwood
    integer mc_itime_anpp_y, mc_itime_bnpp_y, mc_itime_lai_h
    integer mc_itime_npp_y, mc_itime_reco_y
    integer mc_iyear,  mc_iday, mc_ihour

    contains

    subroutine mcmc_functions_init()
        implicit none
        mc_itime_gpp_d  = 1
        mc_itime_nee_d  = 1 
        mc_itime_reco_d = 1
        mc_itime_gpp_h  = 1
        mc_itime_nee_h  = 1
        mc_itime_reco_h = 1
        mc_itime_ch4_h  = 1
        mc_itime_cleaf  = 1
        mc_itime_cwood  = 1

        mc_itime_anpp_y = 1
        mc_itime_bnpp_y = 1 
        mc_itime_lai_h  = 1
        mc_itime_npp_y  = 1
        mc_itime_reco_y = 1

        mc_iyear = 1
        mc_iday  = 1
        mc_ihour = 1
    end subroutine mcmc_functions_init

    subroutine readConfsNml()
        character(20) :: parnames_1, parnames_2, parnames_3, parnames_4, parnames_5 
        character(20) :: parnames_6, parnames_7, parnames_8, parnames_9, parnames_10

        character(20) :: parnames_11, parnames_12, parnames_13, parnames_14, parnames_15 
        character(20) :: parnames_16, parnames_17, parnames_18, parnames_19, parnames_20

        character(20) :: parnames_21, parnames_22, parnames_23, parnames_24, parnames_25 
        character(20) :: parnames_26, parnames_27, parnames_28, parnames_29, parnames_30

        character(20) :: parnames_31, parnames_32, parnames_33, parnames_34, parnames_35 
        character(20) :: parnames_36, parnames_37, parnames_38, parnames_39, parnames_40

        character(20) :: parnames_41, parnames_42, parnames_43, parnames_44, parnames_45 
        character(20) :: parnames_46, parnames_47, parnames_48, parnames_49, parnames_50

        character(20) :: parnames_51, parnames_52, parnames_53, parnames_54, parnames_55 
        character(20) :: parnames_56, parnames_57, parnames_58, parnames_59, parnames_60 

        namelist /nml_obsfiles/ obsfile_gpp_d, obsfile_nee_d, obsfile_reco_d, &
                obsfile_gpp_h, obsfile_nee_h, obsfile_reco_h, obsfile_ch4_h, &
                obsfile_cleaf, obsfile_cwood, obsfile_anpp_y, obsfile_bnpp_y, & 
                obsfile_lai_h, obsfile_npp_y, obsfile_reco_y 

        namelist /nml_param_names/parnames_1, parnames_2, parnames_3, parnames_4, parnames_5, & 
                parnames_6, parnames_7, parnames_8, parnames_9, parnames_10, &
                parnames_11, parnames_12, parnames_13, parnames_14, parnames_15, & 
                parnames_16, parnames_17, parnames_18, parnames_19, parnames_20, &
                parnames_21, parnames_22, parnames_23, parnames_24, parnames_25, &
                parnames_26, parnames_27, parnames_28, parnames_29, parnames_30, &
                parnames_31, parnames_32, parnames_33, parnames_34, parnames_35, &
                parnames_36, parnames_37, parnames_38, parnames_39, parnames_40, &
                parnames_41, parnames_42, parnames_43, parnames_44, parnames_45, & 
                parnames_46, parnames_47, parnames_48, parnames_49, parnames_50, &
                parnames_51, parnames_52, parnames_53, parnames_54, parnames_55, &
                parnames_56, parnames_57, parnames_58, parnames_59, parnames_60

        namelist /nml_mcmc_settings/ nDAsimu, search_scale, ncov, nRand, &
                do_mc_out_hr, do_mc_out_day, do_mc_out_mon, do_mc_out_yr

        open(145, file=mcmc_conf_file)
        read(145, nml=nml_mcmc_settings)
        read(145, nml=nml_obsfiles)
        read(145, nml=nml_param_names)
        close(145)

        parnames(1)  = parnames_1
        parnames(2)  = parnames_2
        parnames(3)  = parnames_3
        parnames(4)  = parnames_4
        parnames(5)  = parnames_5 
        parnames(6)  = parnames_6
        parnames(7)  = parnames_7
        parnames(8)  = parnames_8
        parnames(9)  = parnames_9
        parnames(10) = parnames_10
        parnames(11)  = parnames_11
        parnames(12)  = parnames_12
        parnames(13)  = parnames_13
        parnames(14)  = parnames_14
        parnames(15)  = parnames_15 
        parnames(16)  = parnames_16
        parnames(17)  = parnames_17
        parnames(18)  = parnames_18
        parnames(19)  = parnames_19
        parnames(20)  = parnames_20
        parnames(21)  = parnames_21
        parnames(22)  = parnames_22
        parnames(23)  = parnames_23
        parnames(24)  = parnames_24
        parnames(25)  = parnames_25 
        parnames(26)  = parnames_26
        parnames(27)  = parnames_27
        parnames(28)  = parnames_28
        parnames(29)  = parnames_29
        parnames(30)  = parnames_30
        parnames(31)  = parnames_31
        parnames(32)  = parnames_32
        parnames(33)  = parnames_33
        parnames(34)  = parnames_34
        parnames(35)  = parnames_35 
        parnames(36)  = parnames_36
        parnames(37)  = parnames_37
        parnames(38)  = parnames_38
        parnames(39)  = parnames_39
        parnames(40)  = parnames_40
        parnames(41)  = parnames_41
        parnames(42)  = parnames_42
        parnames(43)  = parnames_43
        parnames(44)  = parnames_44
        parnames(45)  = parnames_45 
        parnames(46)  = parnames_46
        parnames(47)  = parnames_47
        parnames(48)  = parnames_48
        parnames(49)  = parnames_49
        parnames(50)  = parnames_50        
        parnames(51)  = parnames_51
        parnames(52)  = parnames_52
        parnames(53)  = parnames_53
        parnames(54)  = parnames_54
        parnames(55)  = parnames_55 
        parnames(56)  = parnames_56
        parnames(57)  = parnames_57
        parnames(58)  = parnames_58
        parnames(59)  = parnames_59
        parnames(60)  = parnames_60


        ! give the filepath to each variable
        vars4MCMC%gpp_d%filepath  = adjustl(trim(filepath_in))//"/"//adjustl(trim(obsfile_gpp_d))
        vars4MCMC%nee_d%filepath  = adjustl(trim(filepath_in))//"/"//adjustl(trim(obsfile_nee_d))
        vars4MCMC%reco_d%filepath = adjustl(trim(filepath_in))//"/"//adjustl(trim(obsfile_reco_d))
        vars4MCMC%gpp_h%filepath  = adjustl(trim(filepath_in))//"/"//adjustl(trim(obsfile_gpp_h))
        vars4MCMC%nee_h%filepath  = adjustl(trim(filepath_in))//"/"//adjustl(trim(obsfile_nee_h))
        vars4MCMC%reco_h%filepath = adjustl(trim(filepath_in))//"/"//adjustl(trim(obsfile_reco_h))
        ! methane   
        vars4MCMC%ch4_h%filepath  = adjustl(trim(filepath_in))//"/"//adjustl(trim(obsfile_ch4_h))
        ! c pools
        vars4MCMC%cleaf%filepath  = adjustl(trim(filepath_in))//"/"//adjustl(trim(obsfile_cleaf)) ! foliage
        vars4MCMC%cwood%filepath  = adjustl(trim(filepath_in))//"/"//adjustl(trim(obsfile_cwood))

        vars4MCMC%anpp_y%filepath = adjustl(trim(filepath_in))//"/"//adjustl(trim(obsfile_anpp_y))
        vars4MCMC%bnpp_y%filepath = adjustl(trim(filepath_in))//"/"//adjustl(trim(obsfile_bnpp_y))
        vars4MCMC%lai_h%filepath  = adjustl(trim(filepath_in))//"/"//adjustl(trim(obsfile_lai_h))
        vars4MCMC%npp_y%filepath  = adjustl(trim(filepath_in))//"/"//adjustl(trim(obsfile_npp_y))
        vars4MCMC%reco_y%filepath = adjustl(trim(filepath_in))//"/"//adjustl(trim(obsfile_reco_y))
        return
    end subroutine readConfsNml

    subroutine readParamNml(param_nml_file, in_params, init_params, mc_parval, mc_parmin, mc_parmax)
    ! default nml file name of "TECO_MCMC_configs.nml"
        implicit none
        character(*), intent(in) :: param_nml_file
        type(nml_params_data_type), intent(inout)    :: in_params
        type(nml_initValue_data_type), intent(inout) :: init_params
        type(params_mcmc), intent(inout) :: mc_parval, mc_parmin, mc_parmax
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
        real :: stom_n, a1, Ds0                     ! a1 recalculated in vegn model
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
        real :: etaL, etaW, etaR        ! etaL and etaR are not used. ! the percentage of fine litter of the litters from plant parts
        real :: f_F2M, f_C2M, f_C2S
        real :: f_M2S, f_M2P, f_S2P
        real :: f_S2M, f_P2M
        real :: hmax, hl0, LAIMAX0, la0
        ! ----------------------------------------------------
        namelist/nml_params/ lat, lon, wsmax, wsmin, LAImax, LAImin, rdepth,    &
            Rootmax, Stemmax, SapR, SapS, SLAx, GLmax, GRmax, Gsmax, stom_n,    &
            a1, Ds0, Vcmax0, extkU, xfang, alpha, Tau_Leaf, Tau_Wood, Tau_Root, &
            Tau_F, Tau_C, Tau_Micro, Tau_SlowSOM, Tau_Passive, gddonset, Q10,   &
            Q10rh, Rl0, Rs0, Rr0, r_me, Q10pro, kCH4, Omax, CH4_thre, Tveg,     &
            Tpro_me, Toxi, f, bubprob, Vmaxfraction, JV, Entrpy, etaL, etaW,    &
            etaR, f_F2M, f_C2M, f_C2S, f_M2S, f_M2P, f_S2P, f_S2M, f_P2M, hmax, hl0, LAIMAX0, la0 
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
        real :: thd_snow_depth, b_bound     ! tuneice  not sensitive for ice, b_bound no used
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

        ! -------------------------------------------------------------

        namelist /nml_parval/ mc_lat, mc_Longitude, mc_wsmax, mc_wsmin,            &                                                    
                mc_LAIMAX, mc_LAIMIN, mc_rdepth, mc_Rootmax, mc_Stemmax,           &                        
                mc_SapR, mc_SapS, mc_SLA, mc_GLmax, mc_GRmax, mc_Gsmax, mc_stom_n, &            
                mc_a1, mc_Ds0, mc_Vcmx0, mc_extkU, mc_xfang, mc_alpha,             &                    
                mc_Tau_Leaf, mc_Tau_Wood, mc_Tau_Root, mc_Tau_F, &
                mc_Tau_C,  mc_Tau_Micro, mc_Tau_SlowSOM, mc_Tau_Passive, &            
                mc_gddonset, mc_Q10, mc_Rl0, mc_Rs0, mc_Rr0, &        
                mc_r_me, mc_Q10pro, mc_kCH4, mc_Omax, mc_CH4_thre, &
                mc_Tveg, mc_Tpro_me, mc_Toxi, &
                mc_f, mc_bubprob, mc_Vmaxfraction, &                                    
                mc_Q10rh, mc_JV, mc_Entrpy, &                                
                mc_etaL, mc_etaW, mc_etaR, &
                mc_f_F2M, mc_f_C2M, mc_f_C2S, mc_f_M2S, &
                mc_f_M2P, mc_f_S2P, mc_f_S2M, mc_f_P2M
        namelist /nml_parmin/ mc_lat, mc_Longitude, mc_wsmax, mc_wsmin, &                                                    
                mc_LAIMAX, mc_LAIMIN, mc_rdepth, mc_Rootmax, mc_Stemmax, &                        
                mc_SapR, mc_SapS, mc_SLA, mc_GLmax, mc_GRmax, mc_Gsmax, mc_stom_n, &            
                mc_a1, mc_Ds0, mc_Vcmx0, mc_extkU, mc_xfang, mc_alpha, &                    
                mc_Tau_Leaf, mc_Tau_Wood, mc_Tau_Root, mc_Tau_F, &
                mc_Tau_C,  mc_Tau_Micro, mc_Tau_SlowSOM, mc_Tau_Passive, &            
                mc_gddonset, mc_Q10, mc_Rl0, mc_Rs0, mc_Rr0, &        
                mc_r_me, mc_Q10pro, mc_kCH4, mc_Omax, mc_CH4_thre, &
                mc_Tveg, mc_Tpro_me, mc_Toxi, &
                mc_f, mc_bubprob, mc_Vmaxfraction, &                                    
                mc_Q10rh, mc_JV, mc_Entrpy, &                                
                mc_etaL, mc_etaW, mc_etaR, &
                mc_f_F2M, mc_f_C2M, mc_f_C2S, mc_f_M2S, &
                mc_f_M2P, mc_f_S2P, mc_f_S2M, mc_f_P2M
        namelist /nml_parmax/ mc_lat, mc_Longitude, mc_wsmax, mc_wsmin, &                                                    
                mc_LAIMAX, mc_LAIMIN, mc_rdepth, mc_Rootmax, mc_Stemmax, &                        
                mc_SapR, mc_SapS, mc_SLA, mc_GLmax, mc_GRmax, mc_Gsmax, mc_stom_n, &            
                mc_a1, mc_Ds0, mc_Vcmx0, mc_extkU, mc_xfang, mc_alpha, &                    
                mc_Tau_Leaf, mc_Tau_Wood, mc_Tau_Root, mc_Tau_F, &
                mc_Tau_C,  mc_Tau_Micro, mc_Tau_SlowSOM, mc_Tau_Passive, &            
                mc_gddonset, mc_Q10, mc_Rl0, mc_Rs0, mc_Rr0, &        
                mc_r_me, mc_Q10pro, mc_kCH4, mc_Omax, mc_CH4_thre, &
                mc_Tveg, mc_Tpro_me, mc_Toxi, &
                mc_f, mc_bubprob, mc_Vmaxfraction, &                                    
                mc_Q10rh, mc_JV, mc_Entrpy, &                                
                mc_etaL, mc_etaW, mc_etaR, &
                mc_f_F2M, mc_f_C2M, mc_f_C2S, mc_f_M2S, &
                mc_f_M2P, mc_f_S2P, mc_f_S2M, mc_f_P2M

        
        print *, "# read parameters nml file: ", param_nml_file
        open(343, file = param_nml_file)
        read(343, nml  = nml_params,         iostat=io)
        read(343, nml  = nml_initial_values, iostat=io)
        read(343, nml  = nml_parval,         iostat=io)
        call giveValues2par(mc_parval)
        read(343, nml=nml_parmin,         iostat=io)
        call giveValues2par(mc_parmin)
        read(343, nml=nml_parmax,         iostat=io)
        call giveValues2par(mc_parmax)
        close(343)
        ! update the parameters in in_params and init_params
        in_params%lat         = lat
        in_params%lon         = lon
        in_params%wsmax       = wsmax
        in_params%wsmin       = wsmin
        in_params%LAImin      = LAImin
        in_params%LAImax      = LAImax
        in_params%SLAx        = SLAx
        in_params%rdepth      = rdepth
        in_params%Rootmax     = Rootmax
        in_params%Stemmax     = Stemmax
        in_params%SapR        = SapR
        in_params%SapS        = SapS
        in_params%GLmax       = GLmax
        in_params%GRmax       = GRmax
        in_params%Gsmax       = Gsmax
        in_params%stom_n      = stom_n
        in_params%a1          = a1
        in_params%Ds0         = Ds0
        in_params%Vcmax0      = Vcmax0
        in_params%extkU       = extkU
        in_params%xfang       = xfang
        in_params%alpha       = alpha            
        in_params%Tau_Leaf    = Tau_Leaf
        in_params%Tau_Wood    = Tau_Wood
        in_params%Tau_Root    = Tau_Root 
        in_params%Tau_F       = Tau_F
        in_params%Tau_C       = Tau_C  
        in_params%Tau_Micro   = Tau_Micro
        in_params%Tau_slowSOM = Tau_slowSOM
        in_params%Tau_Passive = Tau_Passive
        in_params%gddonset    = gddonset
        in_params%Q10         = Q10
        in_params%Q10rh       = Q10rh    
        in_params%Rl0         = Rl0
        in_params%Rs0         = Rs0
        in_params%Rr0         = Rr0
        ! added for parameters in methane module   
        in_params%r_me        = r_me
        in_params%Q10pro      = Q10pro
        in_params%kCH4        = kCH4
        in_params%Omax        = Omax
        in_params%CH4_thre    = CH4_thre
        in_params%Tveg        = Tveg
        in_params%Toxi        = Toxi
        in_params%Tpro_me     = Tpro_me
        ! add based on Ma et al., 2022
        in_params%f            = f
        in_params%bubprob      = bubprob
        in_params%Vmaxfraction = Vmaxfraction
        ! add based on Ma et al., 2023
        in_params%JV     = JV
        in_params%Entrpy = Entrpy 
        in_params%etaL   = etaL
        in_params%etaW   = etaW
        in_params%etaR   = etaR
        in_params%f_F2M  = f_F2M
        in_params%f_C2M  = f_C2M
        in_params%f_C2S  = f_C2S
        in_params%f_M2S  = f_M2S
        in_params%f_M2P  = f_M2P
        in_params%f_S2P  = f_S2P
        in_params%f_S2M  = f_S2M
        in_params%f_P2M  = f_P2M
        in_params%hmax   = hmax 
        in_params%hl0    = hl0
        in_params%LAIMAX0 = LAIMAX0
        in_params%la0     = la0
        ! =====================================
        init_params%QC             = QC
        init_params%CN0            = CN0            
        init_params%NSCmin         = NSCmin
        init_params%Storage        = Storage
        init_params%nsc            = nsc
        init_params%accumulation   = accumulation
        init_params%SNvcmax        = SNvcmax
        init_params%N_deposit      = N_deposit
        init_params%alphaN         = alphaN
        init_params%NSN            = NSN
        init_params%QNminer        = QNminer
        init_params%N_deficit      = N_deficit
        init_params%thksl          = thksl 
        init_params%FRLEN          = FRLEN
        init_params%liq_water      = liq_water
        init_params%fwsoil         = fwsoil
        init_params%topfws         = topfws
        init_params%omega          = omega
        init_params%zwt            = zwt
        init_params%infilt         = infilt
        init_params%sftmp          = sftmp
        init_params%Tsnow          = Tsnow
        init_params%Twater         = Twater
        init_params%Tice           = Tice
        init_params%G              = G
        init_params%snow_dsim      = snow_dsim
        init_params%dcount         = dcount
        init_params%dcount_soil    = dcount_soil
        init_params%ice_tw         = ice_tw
        init_params%Tsoill         = Tsoill
        init_params%ice            = ice
        init_params%shcap_snow     = shcap_snow
        init_params%condu_snow     = condu_snow
        init_params%condu_b        = condu_b
        init_params%depth_ex       = depth_ex
        init_params%diff_s         = diff_s
        init_params%diff_snow      = diff_snow
        init_params%albedo_snow    = albedo_snow
        init_params%resht          = resht
        init_params%thd_snow_depth = thd_snow_depth
        init_params%b_bound        = b_bound
        init_params%infilt_rate    = infilt_rate
        init_params%fa             = fa
        init_params%fsub           = fsub
        init_params%rho_snow       = rho_snow
        init_params%decay_m        = decay_m
        ! methane module. update: Shuang methane bog species even more shallowly rooted than the tundra. add initials for methane module Shuang version
        init_params%CH4_V              = CH4_V
        init_params%CH4                = CH4
        init_params%Vp                 = Vp
        init_params%bubble_methane_tot = bubble_methane_tot
        init_params%Nbub               = Nbub
        init_params%depth_1            = depth_1
    end subroutine readParamNml

    subroutine readObsData()
        implicit none
        logical toExistOrNot
        integer toCountLines
        ! existornot, data 

        ! gpp_d
        INQUIRE(FILE=vars4MCMC%gpp_d%filepath, EXIST=toExistOrNot)
        vars4MCMC%gpp_d%existOrNot = toExistOrNot
        if (vars4MCMC%gpp_d%existOrNot) then
            call ReadLineNumFromFile(vars4MCMC%gpp_d%filepath, toCountLines)
            allocate(vars4MCMC%gpp_d%obsData(toCountLines, 5))
            call ReadObsDataFromFile(vars4MCMC%gpp_d%filepath, toCountLines, vars4MCMC%gpp_d%obsData)
            allocate(vars4MCMC%gpp_d%mdData(toCountLines, 4))
        endif
        ! nee_d
        INQUIRE(FILE=vars4MCMC%nee_d%filepath, EXIST=toExistOrNot)
        vars4MCMC%nee_d%existOrNot = toExistOrNot
        if (vars4MCMC%nee_d%existOrNot) then
            call ReadLineNumFromFile(vars4MCMC%nee_d%filepath, toCountLines)
            allocate(vars4MCMC%nee_d%obsData(toCountLines, 5))
            call ReadObsDataFromFile(vars4MCMC%nee_d%filepath, toCountLines, vars4MCMC%nee_d%obsData)
            allocate(vars4MCMC%nee_d%mdData(toCountLines, 4))
        endif
        ! reco_d
        INQUIRE(FILE=vars4MCMC%reco_d%filepath, EXIST=toExistOrNot)
        vars4MCMC%reco_d%existOrNot = toExistOrNot
        if (vars4MCMC%reco_d%existOrNot) then
            call ReadLineNumFromFile(vars4MCMC%reco_d%filepath, toCountLines)
            allocate(vars4MCMC%reco_d%obsData(toCountLines, 5))
            call ReadObsDataFromFile(vars4MCMC%reco_d%filepath, toCountLines, vars4MCMC%reco_d%obsData)
            allocate(vars4MCMC%reco_d%mdData(toCountLines, 4))
        endif
        ! gpp_h
        INQUIRE(FILE=vars4MCMC%gpp_h%filepath, EXIST=toExistOrNot)
        vars4MCMC%gpp_h%existOrNot = toExistOrNot
        if (vars4MCMC%gpp_h%existOrNot) then
            call ReadLineNumFromFile(vars4MCMC%gpp_h%filepath, toCountLines)
            allocate(vars4MCMC%gpp_h%obsData(toCountLines, 5))
            call ReadObsDataFromFile(vars4MCMC%gpp_h%filepath, toCountLines, vars4MCMC%gpp_h%obsData)
            allocate(vars4MCMC%gpp_h%mdData(toCountLines, 4))
        endif
        ! nee_h
        INQUIRE(FILE=vars4MCMC%nee_h%filepath, EXIST=toExistOrNot)
        vars4MCMC%nee_h%existOrNot = toExistOrNot
        if (vars4MCMC%nee_h%existOrNot) then
            call ReadLineNumFromFile(vars4MCMC%nee_h%filepath, toCountLines)
            allocate(vars4MCMC%nee_h%obsData(toCountLines, 5))
            call ReadObsDataFromFile(vars4MCMC%nee_h%filepath, toCountLines, vars4MCMC%nee_h%obsData)
            allocate(vars4MCMC%nee_h%mdData(toCountLines, 4))
        endif
        ! reco_h
        INQUIRE(FILE=vars4MCMC%reco_h%filepath, EXIST=toExistOrNot)
        vars4MCMC%reco_h%existOrNot = toExistOrNot
        if (vars4MCMC%reco_h%existOrNot) then
            call ReadLineNumFromFile(vars4MCMC%reco_h%filepath, toCountLines)
            allocate(vars4MCMC%reco_h%obsData(toCountLines, 5))
            call ReadObsDataFromFile(vars4MCMC%reco_h%filepath, toCountLines, vars4MCMC%reco_h%obsData)
            allocate(vars4MCMC%reco_h%mdData(toCountLines, 4))
        endif
        ! ch4_h
        INQUIRE(FILE=vars4MCMC%ch4_h%filepath, EXIST=toExistOrNot)
        vars4MCMC%ch4_h%existOrNot = toExistOrNot
        if (vars4MCMC%ch4_h%existOrNot) then
            call ReadLineNumFromFile(vars4MCMC%ch4_h%filepath, toCountLines)
            allocate(vars4MCMC%ch4_h%obsData(toCountLines, 5))
            call ReadObsDataFromFile(vars4MCMC%ch4_h%filepath, toCountLines, vars4MCMC%ch4_h%obsData)
            allocate(vars4MCMC%ch4_h%mdData(toCountLines, 4))
        endif
        ! cleaf
        INQUIRE(FILE=vars4MCMC%cleaf%filepath, EXIST=toExistOrNot)
        vars4MCMC%cleaf%existOrNot = toExistOrNot
        if (vars4MCMC%cleaf%existOrNot) then
            call ReadLineNumFromFile(vars4MCMC%cleaf%filepath, toCountLines)
            allocate(vars4MCMC%cleaf%obsData(toCountLines, 5))
            call ReadObsDataFromFile(vars4MCMC%cleaf%filepath, toCountLines, vars4MCMC%cleaf%obsData)
            allocate(vars4MCMC%cleaf%mdData(toCountLines, 4))
        endif
        ! cwood
        INQUIRE(FILE=vars4MCMC%cwood%filepath, EXIST=toExistOrNot)
        vars4MCMC%cwood%existOrNot = toExistOrNot
        if (vars4MCMC%cwood%existOrNot) then
            call ReadLineNumFromFile(vars4MCMC%cwood%filepath, toCountLines)
            allocate(vars4MCMC%cwood%obsData(toCountLines, 5))
            call ReadObsDataFromFile(vars4MCMC%cwood%filepath, toCountLines, vars4MCMC%cwood%obsData)
            allocate(vars4MCMC%cwood%mdData(toCountLines, 4))
        endif

        ! anpp_y
        INQUIRE(FILE=vars4MCMC%anpp_y%filepath, EXIST=toExistOrNot)
        vars4MCMC%anpp_y%existOrNot = toExistOrNot
        if (vars4MCMC%anpp_y%existOrNot) then
            call ReadLineNumFromFile(vars4MCMC%anpp_y%filepath, toCountLines)
            allocate(vars4MCMC%anpp_y%obsData(toCountLines, 5))
            call ReadObsDataFromFile(vars4MCMC%anpp_y%filepath, toCountLines, vars4MCMC%anpp_y%obsData)
            allocate(vars4MCMC%anpp_y%mdData(toCountLines, 4))
        endif

        ! bnpp_y
        INQUIRE(FILE=vars4MCMC%bnpp_y%filepath, EXIST=toExistOrNot)
        vars4MCMC%bnpp_y%existOrNot = toExistOrNot
        if (vars4MCMC%bnpp_y%existOrNot) then
            call ReadLineNumFromFile(vars4MCMC%bnpp_y%filepath, toCountLines)
            allocate(vars4MCMC%bnpp_y%obsData(toCountLines, 5))
            call ReadObsDataFromFile(vars4MCMC%bnpp_y%filepath, toCountLines, vars4MCMC%bnpp_y%obsData)
            allocate(vars4MCMC%bnpp_y%mdData(toCountLines, 4))
        endif

        ! lai_h
        INQUIRE(FILE=vars4MCMC%lai_h%filepath, EXIST=toExistOrNot)
        vars4MCMC%lai_h%existOrNot = toExistOrNot
        if (vars4MCMC%lai_h%existOrNot) then
            call ReadLineNumFromFile(vars4MCMC%lai_h%filepath, toCountLines)
            allocate(vars4MCMC%lai_h%obsData(toCountLines, 5))
            call ReadObsDataFromFile(vars4MCMC%lai_h%filepath, toCountLines, vars4MCMC%lai_h%obsData)
            allocate(vars4MCMC%lai_h%mdData(toCountLines, 4))
        endif

        ! npp_y
        INQUIRE(FILE=vars4MCMC%npp_y%filepath, EXIST=toExistOrNot)
        vars4MCMC%npp_y%existOrNot = toExistOrNot
        if (vars4MCMC%npp_y%existOrNot) then
            call ReadLineNumFromFile(vars4MCMC%npp_y%filepath, toCountLines)
            allocate(vars4MCMC%npp_y%obsData(toCountLines, 5))
            call ReadObsDataFromFile(vars4MCMC%npp_y%filepath, toCountLines, vars4MCMC%npp_y%obsData)
            allocate(vars4MCMC%npp_y%mdData(toCountLines, 4))
        endif

        ! reco_y
        INQUIRE(FILE=vars4MCMC%reco_y%filepath, EXIST=toExistOrNot)
        vars4MCMC%reco_y%existOrNot = toExistOrNot
        if (vars4MCMC%reco_y%existOrNot) then
            call ReadLineNumFromFile(vars4MCMC%reco_y%filepath, toCountLines)
            allocate(vars4MCMC%reco_y%obsData(toCountLines, 5))
            call ReadObsDataFromFile(vars4MCMC%reco_y%filepath, toCountLines, vars4MCMC%reco_y%obsData)
            allocate(vars4MCMC%reco_y%mdData(toCountLines, 4))
        endif

    end subroutine readObsData

    subroutine renewMDpars()
        implicit none


st%extkU        = in_params%extkU
st%rdepth       = in_params%rdepth
st%Q10rh        = in_params%Q10rh
st%Q10pro       = in_params%Q10pro
st%r_me         = in_params%r_me
st%Toxi         = in_params%Toxi
st%Omax         = in_params%Omax
st%kCH4         = in_params%kCH4
st%CH4_thre     = in_params%CH4_thre
st%Tveg         = in_params%Tveg
st%bubprob      = in_params%bubprob
st%Vmaxfraction = in_params%Vmaxfraction
st%etaW         = in_params%etaW
st%f            = in_params%f
st%tauC(1)      = in_params%Tau_Leaf*8760.
st%tauC(2)      = in_params%Tau_Wood*8760.
st%tauC(3)      = in_params%Tau_Root*8760.
st%tauC(4)      = in_params%Tau_F*8760.
st%tauC(5)      = in_params%Tau_C*8760.
st%tauC(6)      = in_params%Tau_Micro*8760.
st%tauC(7)      = in_params%Tau_slowSOM*8760.
st%tauC(8)      = in_params%Tau_Passive*8760.
st%Tpro_me      = in_params%Tpro_me
st%f_F2M        = in_params%f_F2M
st%f_C2M        = in_params%f_C2M
st%f_C2S        = in_params%f_C2S
st%f_M2S        = in_params%f_M2S
st%f_M2P        = in_params%f_M2P
st%f_S2P        = in_params%f_S2P
st%f_S2M        = in_params%f_S2M
st%f_P2M        = in_params%f_P2M
! initilize data
st%G            = init_params%G
st%QC           = init_params%QC
st%CN0          = init_params%CN0
st%CN           = init_params%CN0
st%QN           = init_params%QC/init_params%CN0
st%N_deposit    = init_params%N_deposit/8760.   ! Nitrogen input (gN/h/m2, )
st%alphaN       = init_params%alphaN
st%QNminer      = init_params%QNminer
st%N_deficit    = init_params%N_deficit
st%THKSL        = init_params%thksl
st%FRLEN        = init_params%FRLEN
st%liq_water    = init_params%liq_water
st%fwsoil       = init_params%fwsoil
st%topfws       = init_params%topfws
st%omega        = init_params%omega
st%zwt          = init_params%zwt
st%infilt       = init_params%infilt
st%sftmp        = init_params%sftmp
st%Tsnow        = init_params%Tsnow
st%Twater       = init_params%Twater
st%Tice         = init_params%Tice
st%snow_dsim    = init_params%snow_dsim
st%dcount       = init_params%dcount
st%dcount_soil  = init_params%dcount_soil
st%ice_tw       = init_params%ice_tw
st%Tsoill       = init_params%Tsoill
st%ice          = init_params%ice
st%shcap_snow   = init_params%shcap_snow
st%condu_snow   = init_params%condu_snow
st%condu_b      = init_params%condu_b
st%depth_ex     = init_params%depth_ex
st%diff_s       = init_params%diff_s
st%diff_snow    = init_params%diff_snow
st%albedo_snow  = init_params%albedo_snow
st%resht        = init_params%resht
st%thd_snow_depth = init_params%thd_snow_depth
st%fa           = init_params%fa
st%fsub         = init_params%fsub
st%rho_snow     = init_params%rho_snow
st%decay_m      = init_params%decay_m
st%CH4_V        = init_params%CH4_V
st%CH4          = init_params%CH4
st%Vp           = init_params%Vp
st%bubble_methane_tot = init_params%bubble_methane_tot
st%Nbub         = st%Nbub
st%depth(1)     = init_params%depth_1


spec%LAImin = in_params%LAImin
spec%LAImax = in_params%LAImax
spec%SLA    = in_params%SLAx/10000.          ! Convert unit from cm2/g to m2/g
spec%Rootmax = in_params%Rootmax
spec%Stemmax = in_params%Stemmax
spec%SapR    = in_params%SapR
spec%SapS    = in_params%SapS
spec%GLmax   = in_params%GLmax/8760.          ! growth rates of plant. Jian: per year to per hour ?
spec%GRmax   = in_params%GRmax/8760.
spec%Gsmax   = in_params%Gsmax/8760.
spec%stom_n  = in_params%stom_n
spec%Ds0     = in_params%Ds0
spec%Vcmax0  = in_params%Vcmax0
spec%xfang   = in_params%xfang
spec%alpha   = in_params%alpha
spec%gddonset = in_params%gddonset
spec%Q10     = in_params%Q10
spec%Rl0     = in_params%Rl0
spec%Rs0     = in_params%Rs0
spec%Rr0     = in_params%Rr0
spec%JV      = in_params%JV
spec%Entrpy  = in_params%Entrpy
spec%alphaN  = init_params%alphaN
! 
spec%tauC(1)      = in_params%Tau_Leaf*8760.
spec%tauC(2)      = in_params%Tau_Wood*8760.
spec%tauC(3)      = in_params%Tau_Root*8760.
! st%tauC(4)      = in_params%Tau_F
! st%tauC(5)      = in_params%Tau_C
! st%tauC(6)      = in_params%Tau_Micro
! st%tauC(7)      = in_params%Tau_slowSOM
! st%tauC(8)      = in_params%Tau_Passive
spec%QC      = init_params%QC
spec%CN0     = init_params%CN0
spec%CN      = init_params%CN0
spec%QN      = init_params%QC/init_params%CN0
spec%NSCmin  = init_params%NSCmin
spec%storage = init_params%Storage
spec%stor_use= spec%Storage/times_storage_use
spec%nsc     = init_params%nsc
spec%accumulation = init_params%accumulation
spec%SNvcmax = init_params%SNvcmax
spec%NSN     = init_params%NSN
spec%N_deficit = init_params%N_deficit
spec%LAI     = spec%LAImin
spec%bmleaf  = spec%QC(1)/0.48
spec%bmstem  = spec%QC(2)/0.48
spec%bmroot  = spec%QC(3)/0.48
spec%hmax    = in_params%hmax        ! in plant growth hmax = 24.19   ! m
spec%hl0     = in_params%hl0             ! in plant growth hl0  = 0.00019  ! m2/kg C
spec%LAIMAX0 = in_params%LAIMAX0        ! in plant growth LAIMAX0 = 8.    ! maybe the LAImax
spec%la0     = in_params%la0        ! in plant growht la0     = 0.2


        st%lat         = parval(1)
        st%lon         = parval(2)
        st%wsmax       = parval(3)
        st%wsmin       = parval(4)                                            
        LAIMAX      = parval(5)
        LAIMIN      = parval(6)
        rdepth      = parval(7)
        Rootmax     = parval(8)
        Stemmax     = parval(9)                                    
        SapR        = parval(10)
        SapS        = parval(11)
        SLAx        = parval(12)
        GLmax       = parval(13)
        GRmax       = parval(14)
        Gsmax       = parval(15)
        stom_n      = parval(16)         
        a1          = parval(17)
        Ds0         = parval(18)
        Vcmax0      = parval(19)
        extkU       = parval(20)
        xfang       = parval(21)
        alpha       = parval(22)    
        Tau_Leaf    = parval(23)
        Tau_Wood    = parval(24)
        Tau_Root    = parval(25)
        Tau_F       = parval(26)
        Tau_C       = parval(27)
        Tau_Micro   = parval(28)
        Tau_SlowSOM = parval(29)
        Tau_Passive = parval(30)    
        gddonset    = parval(31)
        Q10         = parval(32)
        Rl0         = parval(33)     
        Rs0         = parval(34)    
        Rr0         = parval(35)                    
        r_me        = parval(36)
        Q10pro      = parval(37)
        kCH4        = parval(38)
        Omax         = parval(39)
        CH4_thre     = parval(40)
        Tveg         = parval(41)
        Tpro_me      = parval(42)
        Toxi         = parval(43)        
        f            = parval(44)
        bubprob      = parval(45)
        Vmaxfraction = parval(46)                                    
        Q10rh        = parval(47)
        JV           = parval(48)
        Entrpy       = parval(49)                    
        etaL         = parval(50)
        etaW         = parval(51)
        etaR         = parval(52)
        f_F2M        = parval(53)
        f_C2M        = parval(54)
        f_C2S        = parval(55)
        f_M2S        = parval(56)
        f_M2P        = parval(57)
        f_S2P        = parval(58)
        f_S2M        = parval(59)
        f_P2M        = parval(60)
    end subroutine renewMDpars


    ! subroutine giveValues2var(filepath, existOrNot, data)
    !     implicit none
    !     character(500) filepath
    !     logical existOrNot
    !     real, allocatable :: data(:, :)
    !     integer count_lines

    !     INQUIRE(FILE=filepath, EXIST=existOrNot)
    !     if(existOrNot)then
    !         call ReadLineNumFromFile(filepath, count_lines)
    !         allocate(data(count_lines, 5))
    !         call ReadObsDataFromFile(filepath, count_lines, data)
    !     end if
    !     return
    ! end subroutine giveValues2var

    subroutine giveValues2par(mc_arr)
        implicit none
        type(params_mcmc), intent(inout) :: mc_arr

        mc_arr%lat  = mc_lat
        mc_arr%Longitude  = mc_Longitude 
        mc_arr%wsmax  = mc_wsmax 
        mc_arr%wsmin  = mc_wsmin                                                      
        mc_arr%LAIMAX  = mc_LAIMAX
        mc_arr%LAIMIN  = mc_LAIMIN    
        mc_arr%rdepth  = mc_rdepth    
        mc_arr%Rootmax  = mc_Rootmax    
        mc_arr%Stemmax  = mc_Stemmax                                            
        mc_arr%SapR = mc_SapR    
        mc_arr%SapS = mc_SapS     
        mc_arr%SLA = mc_SLA        
        mc_arr%GLmax = mc_GLmax    
        mc_arr%GRmax = mc_GRmax    
        mc_arr%Gsmax = mc_Gsmax    
        mc_arr%stom_n = mc_stom_n                                            
        mc_arr%a1 = mc_a1       
        mc_arr%Ds0 = mc_Ds0        
        mc_arr%Vcmx0 = mc_Vcmx0    
        mc_arr%extkU = mc_extkU    
        mc_arr%xfang = mc_xfang    
        mc_arr%alpha = mc_alpha                         
        mc_arr%Tau_Leaf = mc_Tau_Leaf   
        mc_arr%Tau_Wood = mc_Tau_Wood   
        mc_arr%Tau_Root = mc_Tau_Root   
        mc_arr%Tau_F = mc_Tau_F       
        mc_arr%Tau_C = mc_Tau_C       
        mc_arr%Tau_Micro = mc_Tau_Micro   
        mc_arr%Tau_SlowSOM = mc_Tau_SlowSOM 
        mc_arr%Tau_Passive = mc_Tau_Passive                             
        mc_arr%gddonset = mc_gddonset    
        mc_arr%Q10 = mc_Q10         
        mc_arr%Rl0 = mc_Rl0        
        mc_arr%Rs0 = mc_Rs0        
        mc_arr%Rr0 = mc_Rr0                            
        mc_arr%r_me = mc_r_me   
        mc_arr%Q10pro = mc_Q10pro   
        mc_arr%kCH4 = mc_kCH4    
        mc_arr%Omax = mc_Omax   
        mc_arr%CH4_thre = mc_CH4_thre 
        mc_arr%Tveg = mc_Tveg  
        mc_arr%Tpro_me = mc_Tpro_me 
        mc_arr%Toxi = mc_Toxi               
        mc_arr%f = mc_f    
        mc_arr%bubprob = mc_bubprob  
        mc_arr%Vmaxfraction = mc_Vmaxfraction                                        
        mc_arr%Q10rh = mc_Q10rh  
        mc_arr%JV = mc_JV   
        mc_arr%Entrpy = mc_Entrpy                                
        mc_arr%etaL = mc_etaL   
        mc_arr%etaW = mc_etaW  
        mc_arr%etaR = mc_etaR   
        mc_arr%f_F2M = mc_f_F2M   
        mc_arr%f_C2M = mc_f_C2M  
        mc_arr%f_C2S = mc_f_C2S 
        mc_arr%f_M2S = mc_f_M2S  
        mc_arr%f_M2P = mc_f_M2P 
        mc_arr%f_S2P = mc_f_S2P  
        mc_arr%f_S2M = mc_f_S2M  
        mc_arr%f_P2M = mc_f_P2M 
        return
    end subroutine giveValues2par

    subroutine GetSimuData(get_iyear, get_iday, get_ihour, in_vegn)
        implicit none
        type(vegn_tile_type), intent(in) :: in_vegn

        integer get_iyear, get_iday, get_ihour
        integer i
        ! vars4MCMC%
        mc_iyear = get_iyear
        mc_iday  = get_iday
        mc_ihour = get_ihour + 1

        ! do i = 1, 20
        !     write(*,*)vars4MCMC%gpp_d%obsData(i, :)
        !     write(*,*)int(vars4MCMC%gpp_d%obsData(i, 1))
        !     write(*,*)int(vars4MCMC%gpp_d%obsData(i, 2))
        !     write(*,*)int(vars4MCMC%gpp_d%obsData(i, 3))
        ! enddo
        ! stop

        ! gpp_d
        if(vars4MCMC%gpp_d%existOrNot)then
            ! write(*,*) "test here: ", vars4MCMC%gpp_d%obsData
            if(mc_itime_gpp_d<=size(vars4MCMC%gpp_d%obsData, dim=1))then
                do while(vars4MCMC%gpp_d%obsData(mc_itime_gpp_d, 1) .lt. forcing(1)%year)
                    vars4MCMC%gpp_d%mdData(mc_itime_gpp_d, 4) = -9999
                    mc_itime_gpp_d = mc_itime_gpp_d + 1
                enddo

                if(vars4MCMC%gpp_d%obsData(mc_itime_gpp_d, 1) .eq. mc_iyear .and. &
                vars4MCMC%gpp_d%obsData(mc_itime_gpp_d, 2) .eq. mc_iday  .and. &
                vars4MCMC%gpp_d%obsData(mc_itime_gpp_d, 3) .eq. mc_ihour) then
                    vars4MCMC%gpp_d%mdData(mc_itime_gpp_d, 1) = mc_iyear
                    vars4MCMC%gpp_d%mdData(mc_itime_gpp_d, 2) = mc_iday
                    vars4MCMC%gpp_d%mdData(mc_itime_gpp_d, 3) = mc_ihour
                    vars4MCMC%gpp_d%mdData(mc_itime_gpp_d, 4) = outVars_d%gpp*86400000
                    mc_itime_gpp_d = mc_itime_gpp_d + 1
                endif
            endif
        endif
        ! nee_d
        if(vars4MCMC%nee_d%existOrNot)then
            if(mc_itime_nee_d <= size(vars4MCMC%nee_d%obsData,dim=1))then
                do while(vars4MCMC%nee_d%obsData(mc_itime_nee_d, 1) .lt. forcing(1)%year)
                    vars4MCMC%nee_d%mdData(mc_itime_nee_d, 4) = -9999
                    mc_itime_nee_d = mc_itime_nee_d + 1
                enddo
                if(vars4MCMC%nee_d%obsData(mc_itime_nee_d, 1) .eq. mc_iyear .and. &
                vars4MCMC%nee_d%obsData(mc_itime_nee_d, 2) .eq. mc_iday  .and. &
                vars4MCMC%nee_d%obsData(mc_itime_nee_d, 3) .eq. mc_ihour) then
                    vars4MCMC%nee_d%mdData(mc_itime_nee_d, 1) = mc_iyear
                    vars4MCMC%nee_d%mdData(mc_itime_nee_d, 2) = mc_iday
                    vars4MCMC%nee_d%mdData(mc_itime_nee_d, 3) = mc_ihour
                    vars4MCMC%nee_d%mdData(mc_itime_nee_d, 4) = outVars_d%nbp*86400000    ! the same in TECO model
                    mc_itime_nee_d = mc_itime_nee_d + 1
                endif
            endif
        endif
        ! reco_d
        if(vars4MCMC%reco_d%existOrNot)then
            if(mc_itime_reco_d <= size(vars4MCMC%reco_d%obsData, dim=1))then
                do while(vars4MCMC%reco_d%obsData(mc_itime_reco_d, 1) .lt. forcing(1)%year)
                    vars4MCMC%reco_d%mdData(mc_itime_reco_d, 4) = -9999
                    mc_itime_reco_d = mc_itime_reco_d + 1
                enddo

                if(vars4MCMC%reco_d%obsData(mc_itime_reco_d, 1) .eq. mc_iyear .and. &
                vars4MCMC%reco_d%obsData(mc_itime_reco_d, 2) .eq. mc_iday  .and. &
                vars4MCMC%reco_d%obsData(mc_itime_reco_d, 3) .eq. mc_ihour) then
                    vars4MCMC%reco_d%mdData(mc_itime_reco_d, 1) = mc_iyear
                    vars4MCMC%reco_d%mdData(mc_itime_reco_d, 2) = mc_iday
                    vars4MCMC%reco_d%mdData(mc_itime_reco_d, 3) = mc_ihour
                    vars4MCMC%reco_d%mdData(mc_itime_reco_d, 4) = (outVars_d%rh + outVars_d%ra)*86400000
                    mc_itime_reco_d = mc_itime_reco_d + 1
                endif
            endif
        endif
        ! gpp_h
        if(vars4MCMC%gpp_h%existOrNot)then
            if(mc_itime_gpp_h <= size(vars4MCMC%gpp_h%obsData, dim=1) )then
                do while(vars4MCMC%gpp_h%obsData(mc_itime_gpp_h, 1) .lt. forcing(1)%year)
                    vars4MCMC%gpp_h%mdData(mc_itime_gpp_h, 4) = -9999
                    mc_itime_gpp_h = mc_itime_gpp_h + 1
                enddo

                if(vars4MCMC%gpp_h%obsData(mc_itime_gpp_h, 1) .eq. mc_iyear .and. &
                vars4MCMC%gpp_h%obsData(mc_itime_gpp_h, 2) .eq. mc_iday  .and. &
                vars4MCMC%gpp_h%obsData(mc_itime_gpp_h, 3) .eq. mc_ihour) then
                    ! write(*,*) "test gpp hourly", mc_itime_gpp_h
                    ! write(*,*) vars4MCMC%gpp_h%obsData(mc_itime_gpp_h, 1), mc_iyear
                    ! write(*,*) vars4MCMC%gpp_h%obsData(mc_itime_gpp_h, 2), mc_iday
                    ! write(*,*) vars4MCMC%gpp_h%obsData(mc_itime_gpp_h, 3), mc_ihour
                    ! write(*,*) outVars_h%gpp, vars4MCMC%gpp_h%obsData(mc_itime_gpp_h, 4)
                    
                    vars4MCMC%gpp_h%mdData(mc_itime_gpp_h, 1) = mc_iyear
                    vars4MCMC%gpp_h%mdData(mc_itime_gpp_h, 2) = mc_iday
                    vars4MCMC%gpp_h%mdData(mc_itime_gpp_h, 3) = mc_ihour
                    vars4MCMC%gpp_h%mdData(mc_itime_gpp_h, 4) = outVars_h%gpp*3600000
                    mc_itime_gpp_h = mc_itime_gpp_h + 1
                endif
            endif
        endif
        ! nee_h
        if(vars4MCMC%nee_h%existOrNot)then
            if(mc_itime_nee_h <= size(vars4MCMC%nee_h%obsData, dim=1)) then
                do while(vars4MCMC%nee_h%obsData(mc_itime_nee_h, 1) .lt. forcing(1)%year)
                    vars4MCMC%nee_h%mdData(mc_itime_nee_h, 4) = -9999
                    mc_itime_nee_h = mc_itime_nee_h + 1
                enddo
                if(vars4MCMC%nee_h%obsData(mc_itime_nee_h, 1) .eq. mc_iyear .and. &
                vars4MCMC%nee_h%obsData(mc_itime_nee_h, 2) .eq. mc_iday  .and. &
                vars4MCMC%nee_h%obsData(mc_itime_nee_h, 3) .eq. mc_ihour) then
                    vars4MCMC%nee_h%mdData(mc_itime_nee_h, 1) = mc_iyear
                    vars4MCMC%nee_h%mdData(mc_itime_nee_h, 2) = mc_iday
                    vars4MCMC%nee_h%mdData(mc_itime_nee_h, 3) = mc_ihour
                    vars4MCMC%nee_h%mdData(mc_itime_nee_h, 4) = outVars_h%nbp*3600000
                    mc_itime_nee_h = mc_itime_nee_h + 1
                endif
            endif
        endif
        ! reco_h
        if(vars4MCMC%reco_h%existOrNot)then
            if(mc_itime_reco_h <= size(vars4MCMC%reco_h%obsData, dim=1))then
                do while(vars4MCMC%reco_h%obsData(mc_itime_reco_h, 1) .lt. forcing(1)%year)
                    vars4MCMC%reco_h%mdData(mc_itime_reco_h, 4) = -9999
                    mc_itime_reco_h = mc_itime_reco_h + 1
                enddo

                if(vars4MCMC%reco_h%obsData(mc_itime_reco_h, 1) .eq. mc_iyear .and. &
                vars4MCMC%reco_h%obsData(mc_itime_reco_h, 2) .eq. mc_iday  .and. &
                vars4MCMC%reco_h%obsData(mc_itime_reco_h, 3) .eq. mc_ihour) then
                    vars4MCMC%reco_h%mdData(mc_itime_reco_h, 1) = mc_iyear
                    vars4MCMC%reco_h%mdData(mc_itime_reco_h, 2) = mc_iday
                    vars4MCMC%reco_h%mdData(mc_itime_reco_h, 3) = mc_ihour
                    vars4MCMC%reco_h%mdData(mc_itime_reco_h, 4) = (outVars_h%rh + outVars_h%ra)*3600000
                    mc_itime_reco_h = mc_itime_reco_h + 1
                endif
            endif
        endif
        ! ch4_h
        if(vars4MCMC%ch4_h%existOrNot)then
            if(mc_itime_ch4_h <= size(vars4MCMC%ch4_h%obsData, dim=1))then
                do while(vars4MCMC%ch4_h%obsData(mc_itime_ch4_h, 1) .lt. forcing(1)%year)
                    vars4MCMC%ch4_h%mdData(mc_itime_ch4_h, 4) = -9999
                    mc_itime_ch4_h = mc_itime_ch4_h + 1
                enddo

                if(vars4MCMC%ch4_h%obsData(mc_itime_ch4_h, 1) .eq. mc_iyear .and. &
                vars4MCMC%ch4_h%obsData(mc_itime_ch4_h, 2) .eq. mc_iday  .and. &
                vars4MCMC%ch4_h%obsData(mc_itime_ch4_h, 3) .eq. mc_ihour) then
                    vars4MCMC%ch4_h%mdData(mc_itime_ch4_h, 1) = mc_iyear
                    vars4MCMC%ch4_h%mdData(mc_itime_ch4_h, 2) = mc_iday
                    vars4MCMC%ch4_h%mdData(mc_itime_ch4_h, 3) = mc_ihour
                    vars4MCMC%ch4_h%mdData(mc_itime_ch4_h, 4) = sum(in_vegn%CH4)*3600000
                    mc_itime_ch4_h = mc_itime_ch4_h + 1
                endif
            endif
        endif
        ! cleaf
        if(vars4MCMC%cleaf%existOrNot)then
            if(mc_itime_cleaf <= size(vars4MCMC%cleaf%obsData, dim=1))then
                do while(vars4MCMC%cleaf%obsData(mc_itime_cleaf, 1) .lt. forcing(1)%year)
                    vars4MCMC%cleaf%mdData(mc_itime_cleaf, 4) = -9999
                    mc_itime_cleaf = mc_itime_cleaf + 1
                enddo

                if(vars4MCMC%cleaf%obsData(mc_itime_cleaf, 1) .eq. mc_iyear .and. &
                vars4MCMC%cleaf%obsData(mc_itime_cleaf, 2) .eq. mc_iday  .and. &
                vars4MCMC%cleaf%obsData(mc_itime_cleaf, 3) .eq. mc_ihour) then
                    vars4MCMC%cleaf%mdData(mc_itime_cleaf, 1) = mc_iyear
                    vars4MCMC%cleaf%mdData(mc_itime_cleaf, 2) = mc_iday
                    vars4MCMC%cleaf%mdData(mc_itime_cleaf, 3) = mc_ihour
                    vars4MCMC%cleaf%mdData(mc_itime_cleaf, 4) = QC(1) 
                    mc_itime_cleaf = mc_itime_cleaf + 1
                endif
            endif
        endif
        ! cwood
        if(vars4MCMC%cwood%existOrNot)then
            if(mc_itime_cwood <= size(vars4MCMC%cwood%obsData, dim=1)) then
                do while(vars4MCMC%cwood%obsData(mc_itime_cwood, 1) .lt. forcing(1)%year)
                    vars4MCMC%cwood%mdData(mc_itime_cwood, 4) = -9999
                    mc_itime_cwood = mc_itime_cwood + 1
                enddo
                if(vars4MCMC%cwood%obsData(mc_itime_cwood, 1) .eq. mc_iyear .and. &
                vars4MCMC%cwood%obsData(mc_itime_cwood, 2) .eq. mc_iday  .and. &
                vars4MCMC%cwood%obsData(mc_itime_cwood, 3) .eq. mc_ihour) then
                    vars4MCMC%cwood%mdData(mc_itime_cwood, 1) = mc_iyear
                    vars4MCMC%cwood%mdData(mc_itime_cwood, 2) = mc_iday
                    vars4MCMC%cwood%mdData(mc_itime_cwood, 3) = mc_ihour
                    vars4MCMC%cwood%mdData(mc_itime_cwood, 4) = QC(2)
                    mc_itime_cwood = mc_itime_cwood + 1
                endif
            endif
        endif

        ! anpp_y
        if(vars4MCMC%anpp_y%existOrNot)then
            if(mc_itime_anpp_y <= size(vars4MCMC%anpp_y%obsData, dim=1)) then
                do while(vars4MCMC%anpp_y%obsData(mc_itime_anpp_y, 1) .lt. forcing(1)%year)
                    vars4MCMC%anpp_y%mdData(mc_itime_anpp_y, 4) = -9999
                    mc_itime_anpp_y = mc_itime_anpp_y + 1
                enddo
                if (vars4MCMC%anpp_y%obsData(mc_itime_anpp_y, 2) .lt. 0) then 
                    vars4MCMC%anpp_y%obsData(mc_itime_anpp_y, 2) = 365
                endif
                if(vars4MCMC%anpp_y%obsData(mc_itime_anpp_y, 1) .eq. mc_iyear .and. &
                vars4MCMC%anpp_y%obsData(mc_itime_anpp_y, 2) .eq. mc_iday  .and. &
                vars4MCMC%anpp_y%obsData(mc_itime_anpp_y, 3) .eq. mc_ihour) then
                    vars4MCMC%anpp_y%mdData(mc_itime_anpp_y, 1) = mc_iyear
                    vars4MCMC%anpp_y%mdData(mc_itime_anpp_y, 2) = mc_iday
                    vars4MCMC%anpp_y%mdData(mc_itime_anpp_y, 3) = mc_ihour
                    vars4MCMC%anpp_y%mdData(mc_itime_anpp_y, 4) = (outVars_y%nppLeaf + outVars_y%nppStem)*3600000*365*24
                    mc_itime_anpp_y = mc_itime_anpp_y + 1
                endif
            endif
        endif

        ! bnpp_y
        if(vars4MCMC%bnpp_y%existOrNot)then
            if(mc_itime_bnpp_y <= size(vars4MCMC%bnpp_y%obsData, dim=1)) then
                do while(vars4MCMC%bnpp_y%obsData(mc_itime_bnpp_y, 1) .lt. forcing(1)%year)
                    vars4MCMC%bnpp_y%mdData(mc_itime_bnpp_y, 4) = -9999
                    mc_itime_bnpp_y = mc_itime_bnpp_y + 1
                enddo
                if (vars4MCMC%bnpp_y%obsData(mc_itime_bnpp_y, 2) .lt. 0) then 
                    vars4MCMC%bnpp_y%obsData(mc_itime_bnpp_y, 2) = 365
                endif
                if(vars4MCMC%bnpp_y%obsData(mc_itime_bnpp_y, 1) .eq. mc_iyear .and. &
                vars4MCMC%bnpp_y%obsData(mc_itime_bnpp_y, 2) .eq. mc_iday  .and. &
                vars4MCMC%bnpp_y%obsData(mc_itime_bnpp_y, 3) .eq. mc_ihour) then
                    vars4MCMC%bnpp_y%mdData(mc_itime_bnpp_y, 1) = mc_iyear
                    vars4MCMC%bnpp_y%mdData(mc_itime_bnpp_y, 2) = mc_iday
                    vars4MCMC%bnpp_y%mdData(mc_itime_bnpp_y, 3) = mc_ihour
                    vars4MCMC%bnpp_y%mdData(mc_itime_bnpp_y, 4) = outVars_y%nppRoot*3600000*365*24
                    mc_itime_bnpp_y = mc_itime_bnpp_y + 1
                endif
            endif
        endif

        ! lai_h
        if(vars4MCMC%lai_h%existOrNot)then
            if(mc_itime_lai_h <= size(vars4MCMC%lai_h%obsData, dim=1)) then
                do while(vars4MCMC%lai_h%obsData(mc_itime_lai_h, 1) .lt. forcing(1)%year)
                    vars4MCMC%lai_h%mdData(mc_itime_lai_h, 4) = -9999
                    mc_itime_lai_h = mc_itime_lai_h + 1
                enddo
                if(vars4MCMC%lai_h%obsData(mc_itime_lai_h, 1) .eq. mc_iyear .and. &
                vars4MCMC%lai_h%obsData(mc_itime_lai_h, 2) .eq. mc_iday  .and. &
                vars4MCMC%lai_h%obsData(mc_itime_lai_h, 3) .eq. mc_ihour) then
                    vars4MCMC%lai_h%mdData(mc_itime_lai_h, 1) = mc_iyear
                    vars4MCMC%lai_h%mdData(mc_itime_lai_h, 2) = mc_iday
                    vars4MCMC%lai_h%mdData(mc_itime_lai_h, 3) = mc_ihour
                    vars4MCMC%lai_h%mdData(mc_itime_lai_h, 4) = outVars_h%lai
                    mc_itime_lai_h = mc_itime_lai_h + 1
                endif
            endif
        endif

        ! npp_y
        if(vars4MCMC%npp_y%existOrNot)then
            if(mc_itime_npp_y <= size(vars4MCMC%npp_y%obsData, dim=1)) then
                do while(vars4MCMC%npp_y%obsData(mc_itime_npp_y, 1) .lt. forcing(1)%year)
                    vars4MCMC%npp_y%mdData(mc_itime_npp_y, 4) = -9999
                    mc_itime_npp_y = mc_itime_npp_y + 1
                enddo
                if (vars4MCMC%npp_y%obsData(mc_itime_npp_y, 2) .lt. 0) then 
                    vars4MCMC%npp_y%obsData(mc_itime_npp_y, 2) = 365
                endif
                if(vars4MCMC%npp_y%obsData(mc_itime_npp_y, 1) .eq. mc_iyear .and. &
                vars4MCMC%npp_y%obsData(mc_itime_npp_y, 2) .eq. mc_iday  .and. &
                vars4MCMC%npp_y%obsData(mc_itime_npp_y, 3) .eq. mc_ihour) then
                    vars4MCMC%npp_y%mdData(mc_itime_npp_y, 1) = mc_iyear
                    vars4MCMC%npp_y%mdData(mc_itime_npp_y, 2) = mc_iday
                    vars4MCMC%npp_y%mdData(mc_itime_npp_y, 3) = mc_ihour
                    vars4MCMC%npp_y%mdData(mc_itime_npp_y, 4) = outVars_y%npp*3600000*365*24
                    mc_itime_npp_y = mc_itime_npp_y + 1
                endif
            endif
        endif

        ! reco_y
        if(vars4MCMC%reco_y%existOrNot)then
            if(mc_itime_reco_y <= size(vars4MCMC%reco_y%obsData, dim=1)) then
                do while(vars4MCMC%reco_y%obsData(mc_itime_reco_y, 1) .lt. forcing(1)%year)
                    vars4MCMC%reco_y%mdData(mc_itime_reco_y, 4) = -9999
                    mc_itime_reco_y = mc_itime_reco_y + 1
                enddo
                if(vars4MCMC%reco_y%obsData(mc_itime_reco_y, 1) .eq. mc_iyear .and. &
                vars4MCMC%reco_y%obsData(mc_itime_reco_y, 2) .eq. mc_iday  .and. &
                vars4MCMC%reco_y%obsData(mc_itime_reco_y, 3) .eq. mc_ihour) then
                    vars4MCMC%reco_y%mdData(mc_itime_reco_y, 1) = mc_iyear
                    vars4MCMC%reco_y%mdData(mc_itime_reco_y, 2) = mc_iday
                    vars4MCMC%reco_y%mdData(mc_itime_reco_y, 3) = mc_ihour
                    vars4MCMC%reco_y%mdData(mc_itime_reco_y, 4) = QC(2)
                    mc_itime_reco_y = mc_itime_reco_y + 1
                endif
            endif
        endif
           
    end subroutine GetSimuData

    

    ! subroutine ReadLineNumFromFile(filepath, count_lines)
    !     implicit none
    !     character(len=*), intent(in) :: filepath
    !     character(len=100) header, line
    !     integer STAT, count_lines

    !     open(38, file=trim(filepath), status="old", action="read", iostat=STAT) ! open file
    !     read(38, '(a100)') header           ! read the header of the file
    !     count_lines = 0                     ! initilize the count_lines
    !     do while(.TRUE.)
    !         read(38, *, iostat=STAT) line   ! read each line
    !         if(STAT .ne. 0) exit            ! until the end of the file
    !         count_lines = count_lines + 1   ! recording the count of the lines
    !     enddo
    !     return
    ! end subroutine ReadLineNumFromFile

    subroutine ReadObsDataFromFile(filepath, count_lines, resData)
        ! Jian: note that this subroutine is used to read the observational data. 
        ! The observational file must be .txt format, and with 5 columns: year, doy, hour, value, std.
        implicit none
        character(len=*), intent(in) :: filepath
        character(len=100) header
        integer STAT, count_lines, iline, n
        real resData(count_lines, 5), readData(5) ! 5 colunms: year, doy, hour, value, std.

        OPEN(34, FILE=trim(filepath), status='old', ACTION='read', IOSTAT=STAT) ! open file
        read(34, '(a100)') header
        iline = 1
        do
            read(34,*,iostat=STAT, end=567) (readData(n), n = 1, 5)
            if(STAT .ne. 0) exit
            resData(iline, :) = readData
            iline = iline + 1
        end do
567     continue
        close(34)
        return
    end subroutine ReadObsDataFromFile

end module mcmc_mod