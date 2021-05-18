#!/bin/bash
#
#PBS -A UQ-QBI
#PBS -l select=1:ncpus=1:mpiprocs=1:mem=10GB:vmem=10GB
#PBS -l walltime=06:00:00
#PBS -o /30days/uqcnolan/super-effects/log/VSL
#PBS -e /30days/uqcnolan/super-effects/log/VSL

# start
PROJPATH=/30days/uqcnolan/super-effects
module load singularity/3.5.0
singularity run -B $PROJPATH/src:/rundir \
                -B $PROJPATH/data:/datadir \
                -B $PROJPATH/out/VSL:/outdir \
                -H $TMPDIR \
                $PROJPATH/images/r-supereffects.sif \
                VSL.R \
                /datadir/total_of_313_subs_VSL_task_trial_level_data.csv \
                /outdir $PBS_ARRAY_INDEX
# OUTFILE=$(printf 'vsl_%02d.zip' $PBS_ARRAYID)
# zip "$PROJPATH/out/VSL/$OUTFILE" "$TMPDIR/*"
