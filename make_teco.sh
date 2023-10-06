gfortran -g -fbacktrace -Wall -fcheck=all   \
    src/teco/datatypes.f90                  \
    src/teco/updateAndSummary.f90          \
    src/teco/io_mod.f90            \
    src/teco/soil.f90                      \
    src/teco/vegetation.f90                 \
    src/teco/transfer.f90                  \
    src/teco/driver.f90                    \
    src/main.f90                            \
    -o run_teco                             \
    -I/usr/include -L/usr/lib -lnetcdff 

if find "src" -name "*.mod" -print -quit | grep -q '.*'; then
    rm src/*.mod
fi

if find "src/teco" -name "*.mod" -print -quit | grep -q '.*'; then
    rm src/teco/*.mod
fi

if find "src/tools" -name "*.mod" -print -quit | grep -q '.*'; then
    rm src/tools/*.mod
fi

if find "." -name "*.mod" -print -quit | grep -q '.*'; then
    rm *.mod
fi