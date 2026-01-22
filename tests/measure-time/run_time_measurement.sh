#!/bin/bash

# Script to run time measurements with all solvers and test files

cd "../../build"

OUTPUT_DIR="testFiles"
RESULTS_DIR="results"
mkdir -p "$RESULTS_DIR"

echo "Running time measurements..."

# Particle counts: 10, 100, 1K, 10K, 100K, 1M
PARTICLE_COUNTS=(10 100 1000 10000 100000 1000000)

# Step counts: 1 10, 100
STEP_COUNTS=(1 10)

for particles in "${PARTICLE_COUNTS[@]}"; do
    for steps in "${STEP_COUNTS[@]}"; do
        testfile="$OUTPUT_DIR/input_p${particles}_s${steps}.in"

        echo "Simulating p=$particles s=$steps"

        # Barnes Hut
        echo "------ BARNES HUT -----"
        # OMP_NUM_THREADS=4 ./simulation "$testfile" BarnesHut 0.2 -200 200 > "$RESULTS_DIR/BH-0.2.p$particles.$steps.log" 2>&1
        ./simulation "$testfile" BarnesHut 0.5 -200 200 > "$RESULTS_DIR/BH-0.5.p$particles.$steps.log" 2>&1
        # OMP_NUM_THREADS=4 ./simulation "$testfile" BarnesHut 0.8 -200 200 > "$RESULTS_DIR/BH-0.8.p$particles.$steps.log" 2>&1

        # # MPI and MPIReduced
        # echo "------ MPI / MPIReduced -----"
        # OMP_NUM_THREADS=2 mpirun -n 2 ./simulation "$testfile" MPI > "$RESULTS_DIR/MPI-n2-omp2.p$particles.$steps.log" 2>&1
        # OMP_NUM_THREADS=2 mpirun -n 2 ./simulation "$testfile" MPIReduced > "$RESULTS_DIR/MPIReduced-n2-omp2.p$particles.$steps.log" 2>&1
        # OMP_NUM_THREADS=4 mpirun -n 2 ./simulation "$testfile" MPI > "$RESULTS_DIR/MPI-n2-omp4.p$particles.$steps.log" 2>&1
        # OMP_NUM_THREADS=4 mpirun -n 2 ./simulation "$testfile" MPIReduced > "$RESULTS_DIR/MPIReduced-n2-omp4.p$particles.$steps.log" 2>&1

        # Serial and SerialReduced
        # echo "------ Serial / SerialReduced -----"
        # OMP_NUM_THREADS=4 ./simulation "$testfile" Serial > "$RESULTS_DIR/Serial.p$particles.$steps.log" 2>&1
        # OMP_NUM_THREADS=4 ./simulation "$testfile" SerialReduced > "$RESULTS_DIR/SerialReduced.p$particles.$steps.log" 2>&1


    done
done

echo "All measurements completed! Results in $RESULTS_DIR"