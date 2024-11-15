#!/bin/bash

# Define the start and end of the range
start=143
end=233

# File to modify
file="config/figures.yaml"

# Temporary file for holding intermediate results
temp_file=$(mktemp)

# Copy everything before the insertion point
awk '/hcg16_chanmap142_pbc_zoom.pdf/{print; exit} 1' "$file" > "$temp_file"

# Generate new figure names and append them
for i in $(seq $start $end); do
    echo "  - hcg16_chanmap${i}_pbc_zoom.pdf" >> "$temp_file"
done

# Copy everything after the insertion point
awk '/hcg16_chanmap234_pbc_zoom.pdf/,0' "$file" >> "$temp_file"

# Replace the original file with the modified one
mv "$temp_file" "$file"

# Clean up
rm -f "$temp_file"

