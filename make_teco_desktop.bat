gfortran -g -fbacktrace -Wall -fcheck=all   ^
    src/teco/datatypes.f90                  ^
    src/teco/updateAndSummary.f90          ^
    src/teco/io_mod.f90            ^
    src/teco/soil.f90                      ^
    src/teco/vegetation.f90                 ^
    src/teco/transfer.f90                  ^
    src/teco/driver.f90                    ^
    src/main.f90                            ^
    -o run_teco                             ^
    -I"C:/Users/jz964/AppData/Local/miniconda3/envs/spruce-da/Library/include/" -L"C:/Users/jz964/AppData/Local/miniconda3/envs/spruce-da/Library/lib/" -lnetcdff -lnetcdf