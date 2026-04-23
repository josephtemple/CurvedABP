#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/../cpp/build" 

TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
DATA_DIR="../../data"
EXPERIMENT_DIR="${TIMESTAMP}"

mkdir -p "${DATA_DIR}/${EXPERIMENT_DIR}"/{sphere,torus,euclidean}

echo "Created: ${EXPERIMENT_DIR}/{sphere,torus,euclidean}"


#for manifold in sphere torus euclidean; do
#    for diffusion in 0.05 0.1 0.5 1.0 5.0; do
#        for mobility in 0.1 1.0 5.0 10.0; do
#            for seed in 1 2 3; do
#                outfile="_d${diffusion}_m${mobility}_s${seed}.h5"
#                echo "Running: $outfile"
#                ./core \
#                    --manifold $manifold \
#                    --diffusion $diffusion \
#                    --mobility $mobility \
#                    --seed $seed \
#                    --output $outfile \
#                    --expdir $EXPERIMENT_DIR
#            done
#        done
#    done
#done


for manifold in sphere torus euclidean; do
    for diffusion in 0.05; do
        for mobility in 0.1; do
            for seed in 1 2; do
                outfile="d${diffusion}_m${mobility}_s${seed}.h5"
                echo "Running: ${manifold} ${outfile}"
                ./core \
                    --manifold $manifold \
                    --diffusion $diffusion \
                    --mobility $mobility \
                    --seed $seed \
                    --output $outfile \
                    --expdir "${DATA_DIR}/${EXPERIMENT_DIR}"
            done
        done
    done
done


echo "All done!"