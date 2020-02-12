#! /bin/bash

#SBATCH -J MPIBWA_32_JOBS
#SBATCH -N 2                            # Ask 2 nodes
#SBATCH -n 4                            # Total number of cores
#SBATCH -c 8                            # use 8 core per mpi job
#SBATCH --tasks-per-node=2              # Ask 2 cores per node
#SBATCH --mem-per-cpu=${MEM}            # See Memory ressources
#SBATCH -t 01:00:00
#SBATCH -o STDOUT_FILE.%j.o
#SBATCH -e STDERR_FILE.%j.e

type mpiBWA

if [[ ! $? -eq 0 ]]
then
    echo "ERROR: mpiBWA not found in your PATH"
    exit 1
fi

output_dir=${HOME}/mpiBWAExample

if [[ ! -d ${output_dir} ]]
then
    mkdir -p ${output_dir}
    echo "INFO: the directory ${output_dir} has been created to store the results generated by mpiBWAExample"
fi

cd ${SLURM_SUBMIT_DIR}

mpirun -n 2 mpiBWA mem -t 2 -o ${HOME}/mpiBWAExample/example.sam data/hg19.small.fa.gz data/HCC1187C_R1_10K.fastq data/HCC1187C_R2_10K.fastq 2> ${output_dir}/mpiBWA.log

echo "INFO: mpiBWA.log is available in ${output_dir}"

