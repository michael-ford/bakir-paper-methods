#!/bin/bash

# Define the data directory
DATA_DIR="../HPRC-assemblies"

# Maximum number of parallel jobs
MAX_JOBS=8

# Current number of jobs
CURRENT_JOBS=0

# Loop over the .fa.gz files in the data directory
find -L "$DATA_DIR" -name "*.fa.gz" | while read -r file; do
    # Extract the subject from the file path
    subject=$(basename "$(dirname "$file")")
    
    # Create an output directory for the subject
    OUTPUT_DIR="${subject}_output"
    mkdir -p "$OUTPUT_DIR"

    # Extract the base name of the file
    base_name=$(basename "$file" .fa.gz)

    # Run kir-annotator in the background
    (
        echo "Processing $file..."
        kir-annotator -o "${OUTPUT_DIR}/${base_name}.yaml" "$file"  > "${OUTPUT_DIR}/${base_name}.stdout" 2> "${OUTPUT_DIR}/${base_name}.stderr"
        echo "$file processed."
    ) &

    # Increment the jobs counter
    ((CURRENT_JOBS++))

    # Wait for some jobs to finish if we've reached the maximum
    if ((CURRENT_JOBS >= MAX_JOBS)); then
        wait -n
        ((CURRENT_JOBS--))
    fi
done

# Wait for all background jobs to finish
wait

echo "KIR-Annotator Workflow completed."
