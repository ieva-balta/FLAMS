#!/bin/bash -l
#SBATCH --clusters=wice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=batch_long
#SBATCH --time=168:00:00
#SBATCH --job-name=rest_ebi_long
#SBATCH --output=logs/rest_ebi_long%j.out
#SBATCH --error=logs/rest_ebi_long%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ieva.balta@student.kuleuven.be     
#SBATCH --mem=4G
#SBATCH --account=lp_cosubi

cd ${SLURM_SUBMIT_DIR}
mkdir -p logs

#load conda environment
source /data/leuven/368/vsc36826/miniconda3/etc/profile.d/conda.sh
conda activate /data/leuven/368/vsc36826/envs/flams

# check environment works
echo "Which python?"
which python
python --version

OUTPUT_DIR="/data/leuven/368/vsc36826/scratch/IBP/ebi_threads_long"

python ebi_threads.py "${OUTPUT_DIR}" --threads 4