#!/bin/bash

# Directory containing PDF files
PDF_DIR="../../outputs/publication/figures/"

# Temp directory for intermediate files
TEMP_DIR="$PDF_DIR/temp_compression"

# Create temp directory if it doesn't exist
mkdir -p "$TEMP_DIR"

# Loop through all PDF files in the specified directory
find "$PDF_DIR" -iname '*.pdf' -exec bash -c '
  for pdf; do
    echo "Compressing $pdf..."
    # Use Ghostscript to compress PDF, outputting to temp directory
    gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 \
       -dPDFSETTINGS=/ebook -dNOPAUSE -dQUIET -dBATCH \
       -sOutputFile="'$TEMP_DIR'/${pdf##*/}" "$pdf"
    
    # Overwrite original PDF with compressed version
    mv -f "'$TEMP_DIR'/${pdf##*/}" "$pdf"
  done
' bash {} +

