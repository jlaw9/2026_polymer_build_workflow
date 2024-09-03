import os
import numpy as np

def extract_density(file_path, start_step):
    densities = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        # Skip the first two lines with headers
        for line in lines[2:]:
            parts = line.strip().split()
            if parts and parts[0].isdigit():
                step = int(parts[0])
                if step >= start_step:
                    try:
                        density = float(parts[1])
                        densities.append(density)
                    except ValueError:
                        continue
    return densities

def calculate_statistics(densities):
    average_density = np.mean(densities)
    std_dev_density = np.std(densities)
    return average_density, std_dev_density

def process_density_directories(main_folder, output_file, start_step):
    with open(output_file, 'w') as out_file:
        for root, dirs, files in os.walk(main_folder):
            for file in files:
                if file == 'density_final-data.txt':  # Update to check for the new file
                    file_path = os.path.join(root, file)
                    densities = extract_density(file_path, start_step)
                    if densities:
                        avg_density, std_dev_density = calculate_statistics(densities)
                        out_file.write(f"Directory: {root}\n")
                        out_file.write(f"Average Density: {avg_density:.6f}\n")
                        out_file.write(f"Standard Deviation: {std_dev_density:.6f}\n")
                        out_file.write("\n")

if __name__ == "__main__":
    main_folder = "/kfs2/projects/invpoly/dlazaren/density_study/PRODUCTION/trimer"  # Your main folder path
    density_output_file = "output.density"
    start_step = 100000  # Your start step for density calculations

    process_density_directories(main_folder, density_output_file, start_step)
