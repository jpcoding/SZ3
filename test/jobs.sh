#!/bin/sh
#SBATCH --job-name=sz3qp        # Job name
#SBATCH --output=sz3qp%j.out    # Standard output and error log
#SBATCH --error=sz3qp%j.err       # Separate error log
#SBATCH --ntasks=128              # Number of MPI tasks
#SBATCH --time=1:00:00              # Time limit (hh:mm:ss)
#SBATCH --partition=short          # Change this if needed
#SBATCH --account=coa_xli281_uksr          # Account name 
#SBATCH --exclusive                 # Request nodes exclusively
#SBATCH --mail-type ALL         # Send email when job starts/ends


module load mpi/latest

test=/scratch/pji228/gittmp/sz3qp/build/test/parallel
data_dir=/pscratch/xli281_uksr/shared/rotstrat4096_temp/
out_dir=/pscratch/xli281_uksr/shared/tmp/
on_config=/scratch/pji228/gittmp/sz3qp/test/sz3_on.cfg
off_config=/scratch/pji228/gittmp/sz3qp/test/sz3_off.cfg
mkdir -p $out_dir 

eb_list=(0.00000001 0.0000001  0.000001 0.00001 0.0001 0.001 0.01 0.1)

for eb in ${eb_list[@]}
do
    echo "eb=${eb}" 
    mpirun -n 64 $test $data_dir $out_dir 512 512 512 $on_config $eb
done

for eb in ${eb_list[@]}
do
    echo "eb=${eb}" 
    mpirun -n 64 $test $data_dir $out_dir 512 512 512 $off_config $eb
done


