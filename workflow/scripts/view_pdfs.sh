#!/bin/bash

# Directory where the PDFs are located
pdf_dir="../../results/publication/figures"

# Temporary directory for converted images
temp_dir="../../results/publication/temp_images"
mkdir -p "$temp_dir"

# Output directory for the montage image
output_dir="../../results/publication/montage_images"
mkdir -p "$output_dir"

# Convert each PDF to an image (first page only) and store in temp_dir
for pdf in "$pdf_dir"/*.pdf; do
	    filename=$(basename "$pdf" .pdf)
	        convert -density 300 "${pdf}[0]" -quality 100 "$temp_dir/${filename}.jpg"
	done

	# Use montage to arrange the images side by side in two rows
	montage "$temp_dir"/*.jpg -tile 2x -geometry +10+10 "$output_dir/montage.jpg"

	# Optionally, remove the temporary image directory
	# rm -r "$temp_dir"

	echo "Montage created in $output_dir/montage.jpg"

