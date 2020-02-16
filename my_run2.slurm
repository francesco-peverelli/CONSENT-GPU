#!/bin/bash -l
#SBATCH -q gpu_preempt
#SBATCH -C gpu
#SBATCH --gres=gpu:1
#SBATCH --exclusive
#SBATCH -L SCRATCH,project
#SBATCH -t 08:00:00
#SBATCH -A m1759 


module load python
source activate fpa
srun -c 80 CONSENT-correct --tmpdir $SCRATCH/reads --in $SCRATCH/reads/d_melanog_20x.fasta  --out $SCRATCH/reads/d_melanog_20x_GPU_out.fasta  --type PB &> errlog3.txt
