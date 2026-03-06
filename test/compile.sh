mpicxx=/mnt/gpfs3_amd/share/apps/Intel/mpi/2021.2.0/bin/mpicxx
sz3_inc=/scratch/pji228/gittmp/sz3qp/include/
zstd_inc=/scratch/pji228/gittmp/sz3qp/install/
zstd_lib=/scratch/pji228/gittmp/sz3qp/install/lib64/libzstd.so
$mpicxx -I$sz3_inc \
        -I$zstd_inc\
        -L$zstd_lib -lzstd \
        -Wl,-rpath -Wl,$zstd_lib  \
        parallel.cpp -o parallel --std=c++17 -lstdc++fs -O