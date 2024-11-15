#!/bin/bash

# Check if at least 4 arguments are provided
if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <file1.fits> <file2.fits> <file3.fits> <file4.fits>"
    exit 1
fi

# Assign inputs to variables for clarity
input_file1="$1"
input_file2="$2"
input_file3="$3"
output_file2="$4"
# Array of input files for convenience
input_files=("$input_file1" "$input_file2" "$input_file3")

# Process each input file
for i in "${!input_files[@]}"; do
    # Compute index starting from 1
    index=$(($i + 1))
    input_file="${input_files[$i]}"
    
    # Check if the temporary output file exists, if so, remove it
    if [ -d "$input_file" ]; then
        rm -r "$input_file" || { echo "Failed to remove existing temporary file $output_file1"; exit 1; }
    fi
    
    # Check if the final output file exists, if so, remove it
    if [ -d "$output_file2" ]; then
        rm -r "$output_file2" || { echo "Failed to remove existing final output file $output_file2"; exit 1; }
    fi
    
    # Perform the xyin operation (assuming a command that works with this syntax)
    fits in="${input_file}.fits" op=xyin out="$input_file" || { echo "Failed to process $input_file with xyin operation"; exit 1; }
done

# Combine the temporary outputs into one
imcomb in="$input_file1","$input_file2","$input_file3" out="$output_file2" options=fqaver || { echo "Failed to combine inputs with imcomb"; exit 1; }

# Convert the combined output back to a FITS file
if [ -f "${input_file2}.fits" ]; then
    rm "${input_file2}.fits"  || { echo "Failed to remove the output file"; exit 1; }
fi
fits in="$output_file2" op=xyout out="${output_file2}.fits" || { echo "Failed to create final FITS file from combined output"; exit 1; }

# Cleanup temporary files
rm -r "$input_file1" "$input_file2" "$input_file3" "$output_file2" || { echo "Failed to clean up temporary files"; exit 1; }

echo "Processing completed successfully."

