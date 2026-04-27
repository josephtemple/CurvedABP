#!/bin/bash

# make cpp/build the working directory so directory references aren't confusing wrt 
# where the bash is and where the cpp is
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/../cpp/build" 

TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
DATA_DIR="../../data"
EXPERIMENT_DIR="${TIMESTAMP}"

mkdir -p "${DATA_DIR}/${EXPERIMENT_DIR}"/{sphere,torus,euclidean}


echo "Created: ${EXPERIMENT_DIR}/{sphere,torus,euclidean}"

for manifold in sphere torus euclidean; do
    for diffusion in 0.05 0.5 5.0; do
        for coupling_str in 0.01 0.1 1.0; do
            for interaction_rad in 0.01 0.1 0.5; do
                for seed in 1 2 3; do
                    outfile="d${diffusion}_c${coupling_str}_i${interaction_rad}_s${seed}.h5"
                    echo "Running: ${manifold} ${outfile}"
                    ./core \
                        --manifold $manifold \
                        --diffusion $diffusion \
                        --coupling_str $coupling_str \
                        --interaction_rad $interaction_rad \
                        --seed $seed \
                        --output $outfile \
                        --expdir "${DATA_DIR}/${EXPERIMENT_DIR}"
                done
            done
        done
    done
done


echo "All done!"