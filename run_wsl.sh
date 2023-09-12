export LD_LIBRARY_PATH=/home/jz964/miniconda3/lib/:$LD_LIBRARY_PATH 
./make_teco.sh
gdb ./run_teco teco_configs.nml
rm run_teco