#!/bin/bash

# Script to run time measurements with all solvers and test files

cd "../../build"

OUTPUT_DIR="testFiles"
RESULTS_DIR="results"
mkdir -p "$RESULTS_DIR"

SOLVERS=("Serial" "SerialReduced" "MPI" "MPIReduced" "BarnesHut")
BARNES_HUT_PARAMS=("0.2 -200 200" "0.5 -200 200" "0.8 -200 200")

echo "Running time measurements..."

# Particle counts: 10, 100, 1K, 10K, 100K, 1M
PARTICLE_COUNTS=(10 100 1000 10000 100000 1000000)

# Step counts: 1 10, 100
STEP_COUNTS=(1 10 100)

for particles in "${PARTICLE_COUNTS[@]}"; do
    for steps in "${STEP_COUNTS[@]}"; do
        testfile="$OUTPUT_DIR/input_p${particles}_s${steps}.in"

        echo "Simulating p=$particles s=$steps"

        # Barnes Hut
        echo "------ BARNES HUT -----"
        OMP_NUM_THREADS=4 ./simulation "$testfile" BarnesHut 0.2 -200 200 > "$RESULTS_DIR/BH-0.2.p$particles.$steps.log" 2>&1
        OMP_NUM_THREADS=4 ./simulation "$testfile" BarnesHut 0.5 -200 200 > "$RESULTS_DIR/BH-0.5.p$particles.$steps.log" 2>&1
        OMP_NUM_THREADS=4 ./simulation "$testfile" BarnesHut 0.8 -200 200 > "$RESULTS_DIR/BH-0.8.p$particles.$steps.log" 2>&1

    done
done

#     for solver in "${SOLVERS[@]}"; do
#         logfile="$RESULTS_DIR/${testname}_${solver}.log"
        
#         echo "Running: $solver with $testname -> $logfile"
        
#         if [ "$solver" == "BarnesHut" ]; then
#             for params in "${BARNES_HUT_PARAMS[@]}"; do
#                 theta=$(echo "$params" | awk '{print $1}')
#                 echo "   BH theta: $theta"
#                 logfile="$RESULTS_DIR/${testname}_${solver}_${theta}.log"
#                 ./simulate "$testfile" "$solver" $params > "$logfile" 2>&1
#             done
#         elif [ "$solver" == "MPI" ] || [ "$solver" == "MPIReduced" ]; then
#                 mpirun -n 4 ./simulate "$testfile" "$solver" > "$logfile" 2>&1
#         else
#             particles=$(echo "$testname" | grep -oP 'p\K[0-9]+')
#             steps=$(echo "$testname" | grep -oP 's\K[0-9]+')

#             # echo "   p $particles s $steps"
            
#             if [ "$particles" -lt 1001 ] && [ "$steps" -lt 1001 ]; then
#                 ./simulate "$testfile" "$solver" > "$logfile" 2>&1
#             else 
#                 echo "Skipped"
#             fi
#         fi
#     done
# done

echo "All measurements completed! Results in $RESULTS_DIR"