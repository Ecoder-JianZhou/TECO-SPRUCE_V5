export LD_LIBRARY_PATH=/home/zhou_j/miniconda3/envs/spruce/lib/:$LD_LIBRARY_PATH 
./make_mac.sh
./run_teco teco_configs_server_P04.nml
rm run_teco
