#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=3-20:00:00
#SBATCH --output=%x-%j.out
#SBATCH --mail-user barbara.rentroiapacheco@ucsf.edu
#SBATCH --mail-type BEGIN,END,FAIL


module load CBI miniforge3/24.11.0-0 # This gives access to conda
conda activate intogen
cd /c4/home/barentroia/cga_intogen/EMC_cohort/GROUP/
nextflow run /c4/home/barentroia/softwares/intogen/intogen-plus/intogen.nf -resume -profile local --input /c4/home/barentroia/cga_intogen/EMC_cohort/GROUP/input/ --datasets /c4/home/barentroia/softwares/intogen/intogen-plus/datasets --containers /c4/home/barentroia/softwares/intogen/intogen-plus/containers --output /c4/home/barentroia/cga_intogen/EMC_cohort/GROUP/output/