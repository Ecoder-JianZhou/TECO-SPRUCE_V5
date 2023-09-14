program test
    implicit none
    integer a, io
    namelist /nml_test/ a
    print *, "# read TECO config nml file ...", adjustl(trim("configs/test.nml"))
    open(388, file = adjustl(trim("configs/test.nml")))
    read(388, nml  = nml_test,  iostat=io)
    close(388)
    print*, a
end program test