module driver
    use datatypes
    contains

    subroutine driver_simu_multi_pft()
        implicit none

        st%Dair = 0
        st%rain_d = 0
    end subroutine driver_simu_multi_pft
end module driver