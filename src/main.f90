program TECO
    use datatypes
    use driver
    
    implicit none
    integer :: count_mode
    integer num_args, ierr
    type(nml_params_data_type)    :: in_params
    type(nml_initValue_data_type) :: init_params

    print *, ""
    write(*,*) "# -----------------------------------------"
    ! get the count of command
    num_args = COMMAND_ARGUMENT_COUNT()
    ! check if have the file
    if (num_args /= 1) then
        write(*,*) "Usage: ./my_program <config_file>"
        stop
    end if
    ! read the command
    call GET_COMMAND_ARGUMENT(1, config_file, ierr)
    config_file = adjustl(trim("configs/"))//adjustl(trim(config_file))
    write(*,*) "Reading the file of ", config_file

    print *, ""
    write(*,*) "# -----------------------------------------"
    call read_teco_configs()
    ! check the three mode: do_simu; do_mcmc; do_spinup
    count_mode = 0
    if (do_simu)   count_mode = count_mode + 1
    if (do_mcmc)   count_mode = count_mode + 1
    if (do_spinup) count_mode = count_mode + 1

    if (count_mode .eq. 0) then
        print *, "# You must choose a run mode."
        print *, "# Please make sure one of the three modes is True in file of TECO_model_configs.nml"
        print *, "#    *do_simu; *do_mcmc; *do_spinup"
        print *, ""
        stop
    elseif (count_mode .gt. 1) then
        print *, "# You can only select one mode out of the three, please check your file."
        print *, "# Please check the file of TECO_model_configs.nml"
        print *, "#    *do_simu; *do_mcmc; *do_spinup"
        print *, ""
        stop
    else
        continue
    endif

    write(*,*) "# Start to run the case of """, adjustl(trim(simu_name)), """"
    ! update the in-out path and create the revelent ouput paths
    call createNewCase() 

    call get_forcingdata()                      ! read forcing data
    nHours  = nforcing                          
    nDays   = int(nHours/24.)
    nYears  = int(nforcing/(365*24))
    nMonths = nYears*12
    do ipft = 1, count_pft
        call read_parameters_nml(files_pft_params(ipft), in_params, init_params)
    enddo

    if(do_out_hr)  call assign_outVars(nHours,  outVars_h, count_pft)
    if(do_out_day) call assign_outVars(nDays,   outVars_d, count_pft)
    if(do_out_mon) call assign_outVars(nMonths, outVars_m, count_pft)
    if(do_out_yr)  call assign_outVars(nYears,  outVars_y, count_pft)

    if (.not. do_snow) call get_snowdepth()
    
! update the vegn LAImax LAImin
end program TECO


subroutine createNewCase()
    use datatypes
    ! use mcmc_functions
    ! use mod_mcmc
    implicit none
    ! create a new case to run the TECO model
    !   * create the output path

    print *, "# Update and create the output dirs"

    ! update the full path of input file
    climfile        = adjustl(trim(filepath_in))//"/"//adjustl(trim(climfile))       ! climate file name
    snowdepthfile   = adjustl(trim(filepath_in))//"/"//adjustl(trim(snowdepthfile))  ! snow depthfile
    restartfile     = adjustl(trim(filepath_in))//"/"//adjustl(trim(restartfile))    ! restartfile
    watertablefile  = adjustl(trim(filepath_in))//"/"//adjustl(trim(watertablefile)) ! Jian: maybe used when not run soil_physical
    ! check the inputfile
    call check_inputfile(climfile, "climate file")
    if (.not. do_snow)    call check_inputfile(snowdepthfile,  "The file of snowdepth")
    if (do_restart)       call check_inputfile(restartfile,    "The file of restart")
    if (.not. do_soilphy) call check_inputfile(watertablefile, "The file of water table")
    ! create the outdir
    call CreateFolder(adjustl(trim(outdir)))

    ! update and create the output dir of case
    outdir_case = adjustl(trim(outdir))//"/"//adjustl(trim(simu_name))
    call CreateFolder(adjustl(trim(outdir_case)))

    ! update and create the output dir of each format outputs
    outDir_nc  = adjustl(trim(outdir_case))//"/"//adjustl(trim(outDir_nc))
    call CreateFolder(adjustl(trim(outDir_nc)))
    outDir_csv = adjustl(trim(outdir_case))//"/"//adjustl(trim(outDir_csv))
    call CreateFolder(adjustl(trim(outDir_csv)))

    ! update and create the output for each time frequency of nc-format outputs
    outDir_h = adjustl(trim(outdir_nc))//"/"//adjustl(trim(outDir_h))
    outDir_d = adjustl(trim(outdir_nc))//"/"//adjustl(trim(outDir_d))
    outDir_m = adjustl(trim(outdir_nc))//"/"//adjustl(trim(outDir_m))
    outDir_y = adjustl(trim(outdir_nc))//"/"//adjustl(trim(outDir_y))
    
    call CreateFolder(adjustl(trim(outDir_h)))
    call CreateFolder(adjustl(trim(outDir_d)))
    call CreateFolder(adjustl(trim(outDir_m)))
    call CreateFolder(adjustl(trim(outDir_y)))

    if (do_spinup)then
        outDir_spinup = adjustl(trim(outdir_case))//"/"//adjustl(trim(outDir_spinup))
        call CreateFolder(adjustl(trim(outDir_spinup)))
        outfile_restart = adjustl(trim(outDir_spinup))//"/restart.nc"
    endif

    if (do_mcmc)then
        outDir_mcmc = adjustl(trim(outdir_case))//"/"//adjustl(trim(outDir_mcmc))
        call CreateFolder(adjustl(trim(outDir_mcmc)))
        ! if (do_mc_out_hr)then
            outDir_mcmc_h = adjustl(trim(outDir_mcmc))//"/"//adjustl(trim(outDir_mcmc_h))
            call CreateFolder(adjustl(trim(outDir_mcmc_h)))
        ! endif
        ! if (do_mc_out_day) then
            outDir_mcmc_d = adjustl(trim(outDir_mcmc))//"/"//adjustl(trim(outDir_mcmc_d))
            call CreateFolder(adjustl(trim(outDir_mcmc_d)))
        ! endif
        ! if (do_mc_out_mon) then
            outDir_mcmc_m = adjustl(trim(outDir_mcmc))//"/"//adjustl(trim(outDir_mcmc_m))
            call CreateFolder(adjustl(trim(outDir_mcmc_m)))
        ! endif
    endif

    if(do_restart)then
        restartfile = adjustl(trim(filepath_in))//adjustl(trim(restartfile))
    endif
end subroutine createNewCase