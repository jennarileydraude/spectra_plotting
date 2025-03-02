import numpy as np
import matplotlib.pyplot as plt
import os

def get_dihedral_angle(atom_1_index, atom_2_index, atom_3_index, atom_4_index, atom_array):
    vector_1 = atom_array[atom_2_index] - atom_array[atom_1_index]
    vector_2 = atom_array[atom_3_index] - atom_array[atom_2_index]
    vector_3 = atom_array[atom_4_index] - atom_array[atom_3_index]
    
    cross_12 = np.cross(vector_1, vector_2)
    cross_23 = np.cross(vector_2, vector_3)
    
    term_1 = np.dot(vector_2, np.cross(cross_12, cross_23))
    term_2 = np.linalg.norm(vector_2) * np.dot(cross_12, cross_23)
    
    dihedral_angle = np.arctan2(term_1, term_2)
    return np.rad2deg(dihedral_angle)

def read_xyz_file(file):
    with open(file, 'r') as xyz_file:
        next(xyz_file)
        next(xyz_file)
        atom_data_list = []
        for line in xyz_file:
            if line.strip():  # Skip empty lines
                atom_data = line.split()
                atom_coordinates = [float(coord) for coord in atom_data[1:]]
                atom_data_list.append(atom_coordinates)
    return np.array(atom_data_list)

def process_xyz_files(directory, atom_indices):
    dihedral_angles = []
    file_list = sorted([os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.xyz')])
    
    for file in file_list:
        atom_array = read_xyz_file(file)
        angle = get_dihedral_angle(*atom_indices, atom_array)
        dihedral_angles.append(angle)
    return dihedral_angles

def plot_dihedral_angles(dihedral_angles):
    time = np.arange(1, len(dihedral_angles) + 1)  # 1 Frame = 1 ns ?
    plt.figure(figsize=(10, 6))
    plt.scatter(time, dihedral_angles, marker='o', linestyle='-', label='Dihedral Angles')
    plt.xlabel('Time (ns)')
    plt.ylabel('Dihedral Angle (Degrees)')
    plt.title('Dihedral Angles for Length 6 Oligomer')
    plt.savefig('my_plot_high_dpi.png')

# Call functions
directory_path = 'your_directory'  
atom_indices = (269, 262, 261, 260)  
dihedral_angles = process_xyz_files(directory_path, atom_indices)
plot_dihedral_angles(dihedral_angles)
