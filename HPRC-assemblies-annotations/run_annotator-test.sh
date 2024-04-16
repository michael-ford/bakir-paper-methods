#!/bin/bash

# Check if an input file was provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

# Get the input file from the command line
input_file="$1"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: File '$input_file' not found."
    exit 1
fi

# Extract the subject name from the input file path
subject="$(basename "$(dirname "$input_file")")"

# Create a directory based on the subject name
mkdir -p "$subject"

# Extract the base name of the file without the extension
base_name="$(basename "$input_file" .fa.gz)"

# Run kir-annotator and redirect output and error
kir-annotator -o "$subject/$base_name" -m temp "$input_file"

echo "Processing of '$input_file' completed."

