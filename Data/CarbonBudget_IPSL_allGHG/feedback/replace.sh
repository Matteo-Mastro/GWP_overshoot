#!/bin/bash

# Name of the file containing the desired header
header_file="column.csv"

# Loop through all CSV files in the current directory
for file in *.csv; do
    # Skip the header file itself
    if [ "$file" != "$header_file" ]; then
        # Substitute the first row of the current file with the contents of the header file
        sed '1d' "$file" > "$file.tmp"   # Remove the first row from the original file
        cat "$header_file" "$file.tmp" > "$file"   # Concatenate the header and the modified file
        rm "$file.tmp"   # Remove the temporary file
    fi
done

