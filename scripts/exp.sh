#!/bin/bash

# set up directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/../cpp/build" 

TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
DATA_DIR="../../data"
EXPERIMENT_DIR="${TIMESTAMP}"
FULL_PATH="${DATA_DIR}/${EXPERIMENT_DIR}"

mkdir -p "${FULL_PATH}"/{sphere,torus,euclidean}

# Define the log file inside the experiment directory
LOG_FILE="${FULL_PATH}/sweep.log"

# default parameters
N=250
V0=1.0
RADIUS=0.02
INT_RAD=0.5
DT=0.005
STEPS=2000
SAVE_EVERY=1

# swept parameters
ds=(0.0 0.2 1.0 5.0)
cs=(0.0 0.2 1.0 5.0)
ss=(10 20 30)

# put everything in the log
{
    echo "------------------------------------------------"
    echo "SWEEP START : $(date)"
    echo "Output Dir  : ${FULL_PATH}"
    echo "Total Sims  : $((${#ds[@]} * ${#cs[@]} * ${#ss[@]} * 3))"
    echo "------------------------------------------------"
    echo "Fixed Params:"
    echo "N           : $N"
    echo "v0          : $V0"
    echo "r (radius)  : $RADIUS"
    echo "dt          : $DT"
    echo "steps       : $STEPS"
    echo "inter_rad   : $INT_RAD"
    echo "------------------------------------------------"
    echo "Swept Params:"
    echo "Metrics     : sphere, torus, euclidean"
    echo "D (diff)    : ${ds[*]}"
    echo "C (clust)   : ${cs[*]}"
    echo "S (speed)   : ${ss[*]}"
    echo "------------------------------------------------"
} | tee -a "$LOG_FILE"

MAX_JOBS=7 

# sweep through parameters
for manifold in sphere torus euclidean; do
    for diffusion in "${ds[@]}"; do
        for coupling_str in "${cs[@]}"; do
            for seed in "${ss[@]}"; do
                
                outfile="d${diffusion}_c${coupling_str}_s${seed}.h5"
                
                echo "Launching: ${manifold} ${outfile}"
                
                ./core \
                    --N "$N" \
                    --v0 "$V0" \
                    --radius "$RADIUS" \
                    --interaction_rad "$INT_RAD" \
                    --dt "$DT" \
                    --steps "$STEPS" \
                    --saveEvery "$SAVE_EVERY" \
                    --manifold "$manifold" \
                    --diffusion "$diffusion" \
                    --coupling_str "$coupling_str" \
                    --seed "$seed" \
                    --output "$outfile" \
                    --expdir "${FULL_PATH}" &

                # job control
                while [ $(jobs -rp | wc -l) -ge $MAX_JOBS ]; do
                    sleep 0.5 
                done
            done
        done
    done
done

wait
echo "All done! Log saved to ${LOG_FILE}"