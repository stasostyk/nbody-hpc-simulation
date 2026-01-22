#!/bin/bash

# Script to generate test inputs with various particle counts and step counts

cd "../../build"

OUTPUT_DIR="testFiles"
mkdir -p "$OUTPUT_DIR"

# Particle counts: 10, 100, 1K, 10K, 100K, 1M
PARTICLE_COUNTS=(10 100 1000 10000 100000 1000000)

# Step counts: 1 10, 100
STEP_COUNTS=(1 10 100)

# Fixed parameters
DT=0.01
SEED=42

echo "Generating test inputs in $OUTPUT_DIR..."

for particles in "${PARTICLE_COUNTS[@]}"; do
    for steps in "${STEP_COUNTS[@]}"; do
        filename="$OUTPUT_DIR/input_p${particles}_s${steps}.in"
        echo "Generating: $filename (particles=$particles, steps=$steps)"
        ./generate "$filename" "$particles" "$steps" "$DT" "$SEED"
        
        if [ $? -ne 0 ]; then
            echo "Error generating $filename"
            exit 1
        fi
    done
done

echo "All test inputs generated successfully!"