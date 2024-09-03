#!/bin/bash


##COPY FILES FROM TIM
#Source directory containing the nested folders
src_dir="/kfs2/projects/invpoly/timbernat/polymers_unique"

#Destination directory where you want to copy the files
dest_dir="/kfs2/projects/invpoly/dlazaren/density_study/PRODUCTION/trimer"

##Find and copy .lammps files, excluding directories named "sub_10000_atoms" and "pentamer"
find "$src_dir" \( -type d -name "sub_10000_atoms" -o -type d -name "pentamer" \) -prune -o -type f -name "*.lammps" -print | while read -r file; do
  # Get the base name of the file without the extension
  base_name=$(basename "$file" .lammps)
  
  #Create a directory in the destination with the base name
  new_dir="$dest_dir/$base_name"
  mkdir -p "$new_dir"
  
  #Copy the .lammps file to the new directory
  cp "$file" "$new_dir"
done

##MODIFY DATA FILES
#Find all .lammps files in the destination directory
find "$dest_dir" -type f -name "*.lammps" | while read -r file; do
  # Create a temporary file
  tmp_file=$(mktemp)

  # Process each .lammps file
  awk '
  # Remove the specific line
  !/0.0 0.0 0.0 xy xz yz/ {
    print
  }
  ' "$file" | sed -e 's/harmonic//g' -e 's/fourier//g' > "$tmp_file"

  # Replace the original file with the modified one
  mv "$tmp_file" "$file"
done