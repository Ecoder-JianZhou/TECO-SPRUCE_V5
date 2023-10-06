gfortran -g -fbacktrace -Wall -fcheck=all   ^
    src/teco/datatypes.f90                  ^
    src/teco/updateAndSummary.f90          ^
    src/teco/io_mod.f90            ^
    src/teco/soil.f90                      ^
    src/teco/vegetation.f90                 ^
    src/teco/transfer.f90                  ^
    src/teco/driver.f90                    ^
    src/main.f90                            ^
    -o run_teco                             
    @REM -I"C:/Users/jz964/AppData/Local/miniconda3/envs/spruce-da/Library/include/" -L"C:/Users/jz964/AppData/Local/miniconda3/envs/spruce-da/Library/lib/" -lnetcdff -lnetcdf

@echo off
setlocal

@REM set "folder_path=C:\Your\Folder\Path"

@REM :: Navigate to the target folder
@REM cd /d "%folder_path%"

:: Delete all files with the .mod extension
for /r %%i in (*.mod) do (
    del "%%i"
)

echo ".mod files deleted successfully."

run_teco.exe teco_configs_desktop.nml