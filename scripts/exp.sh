#!/bin/bash

# Setup directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/../cpp/build" 

TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
DATA_DIR="../../data"
EXPERIMENT_DIR="${TIMESTAMP}"

mkdir -p "${DATA_DIR}/${EXPERIMENT_DIR}"/{sphere,torus,euclidean}

echo "Created: ${EXPERIMENT_DIR}/{sphere,torus,euclidean}"

# limit the number of concurrent jobs (8 cores, leave one for other tasks)
MAX_JOBS=7 

for manifold in sphere torus euclidean; do
    for diffusion in 0.1 5.0; do   
        for coupling_str in 0.1 5.0; do
            for interaction_rad in 0.2 0.5; do
                for seed in 10 20; do
                    
                    outfile="d${diffusion}_c${coupling_str}_i${interaction_rad}_s${seed}.h5"
                    
                    # launch simulation in the background
                    echo "Launching: ${manifold} ${outfile}"
                    ./core \
                        --manifold $manifold \
                        --diffusion $diffusion \
                        --coupling_str $coupling_str \
                        --interaction_rad $interaction_rad \
                        --seed $seed \
                        --output $outfile \
                        --expdir "${DATA_DIR}/${EXPERIMENT_DIR}" &

                    # check how many jobs are running
                    while [ $(jobs -rp | wc -l) -ge $MAX_JOBS ]; do
                        sleep 1 # Wait a second before checking again
                    done

                done
            done
        done
    done
done

wait
echo "All done!"