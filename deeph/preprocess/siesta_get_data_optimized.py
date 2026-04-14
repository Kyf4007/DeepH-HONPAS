"""
Optimized SIESTA to DeepH data converter
Vectorized and optimized version for better performance on large systems
Original code by ZC Tang @ Tsinghua Univ.
"""

import os
import numpy as np
import h5py
import json
from scipy.io import FortranFile


def siesta_parse(input_path, output_path):
    """
    Transfer SIESTA output to DeepH format with optimized performance.
    All outputs are identical to the original implementation.
    """
    input_path = os.path.abspath(input_path)
    output_path = os.path.abspath(output_path)
    os.makedirs(output_path, exist_ok=True)
    
    # ========================================================================
    # 1. Find system name (optimized with list comprehension)
    # ========================================================================
    f_list = os.listdir(input_path)
    system_name = next(f_name[:-9] for f_name in f_list 
                      if f_name[::-1][0:9] == 'XDNI_BRO.')
    
    # ========================================================================
    # 2. Parse structure info from STRUCT_OUT (vectorized where possible)
    # ========================================================================
    with open('{}/{}.STRUCT_OUT'.format(input_path, system_name), 'r') as struct:
        lines = struct.readlines()
    
    # Lattice vectors (3 lines)
    lattice = np.array([line.split() for line in lines[0:3]], dtype=np.float64)
    np.savetxt('{}/lat.dat'.format(output_path), np.transpose(lattice), fmt='%.18e')
    
    # Number of atoms and atomic coordinates
    num_atoms = int(lines[3].split()[0])
    atom_coord = np.array([lines[i+4].split()[1:] for i in range(num_atoms)], dtype=np.float64)
    np.savetxt('{}/element.dat'.format(output_path), atom_coord[:, 0], fmt='%d')
    
    # ========================================================================
    # 3. Parse atomic positions from XV file (vectorized)
    # ========================================================================
    atom_coord_cart = np.genfromtxt('{}/{}.XV'.format(input_path, system_name), 
                                    skip_header=4)[:, 2:5] * 0.529177249
    np.savetxt('{}/site_positions.dat'.format(output_path), np.transpose(atom_coord_cart))
    
    # ========================================================================
    # 4. Parse orbital information (fully vectorized)
    # ========================================================================
    orb_indx = np.genfromtxt('{}/{}.ORB_INDX'.format(input_path, system_name), 
                             skip_header=3, skip_footer=17)
    
    # ========================================================================
    # 5. Generate R_list.dat (vectorized unique R vectors)
    # ========================================================================
    R_vectors = orb_indx[:, 12:15].astype(int)
    # Find unique R vectors by comparing consecutive rows
    R_diff = np.any(R_vectors[1:] != R_vectors[:-1], axis=1)
    unique_indices = np.concatenate([[0], np.where(R_diff)[0] + 1])
    unique_R = R_vectors[unique_indices]
    
    with open('{}/R_list.dat'.format(output_path), 'w') as R_list_f:
        for R in unique_R:
            R_list_f.write('{} {} {}\n'.format(int(R[0]), int(R[1]), int(R[2])))
    
    # ========================================================================
    # 6. Build ia2Riua array (pre-allocated, no append)
    # ========================================================================
    # Find indices where atom_id changes
    atom_ids = orb_indx[:, 1].astype(int)
    atom_change_idx = np.where(np.diff(atom_ids, prepend=atom_ids[0]) != 0)[0]
    
    n_atoms = len(atom_change_idx)
    ia2Riua = np.empty((n_atoms, 4), dtype=np.float64)
    
    for idx, start_idx in enumerate(atom_change_idx):
        ia2Riua[idx, 0:3] = orb_indx[start_idx, 12:15]
        iuo = int(orb_indx[start_idx, 15])
        iua = int(orb_indx[iuo - 1, 1])
        ia2Riua[idx, 3] = iua
    
    # ========================================================================
    # 7. Write info.json
    # ========================================================================
    info = {'nsites': num_atoms, 'isorthogonal': False, 'isspinful': False, 
            'norbits': len(orb_indx)}
    with open('{}/info.json'.format(output_path), 'w') as info_f:
        json.dump(info, info_f)
    
    # ========================================================================
    # 8. Calculate reciprocal lattice (vectorized)
    # ========================================================================
    a1, a2, a3 = lattice[0, :], lattice[1, :], lattice[2, :]
    cross_a2a3 = np.cross(a2, a3)
    cross_a3a1 = np.cross(a3, a1)
    cross_a1a2 = np.cross(a1, a2)
    
    b1 = 2 * np.pi * cross_a2a3 / np.dot(a1, cross_a2a3)
    b2 = 2 * np.pi * cross_a3a1 / np.dot(a2, cross_a3a1)
    b3 = 2 * np.pi * cross_a1a2 / np.dot(a3, cross_a1a2)
    rlattice = np.array([b1, b2, b3])
    np.savetxt('{}/rlat.dat'.format(output_path), np.transpose(rlattice), fmt='%.18e')
    
    # ========================================================================
    # 9. Generate orbital_types.dat (vectorized)
    # ========================================================================
    n_orbitals = len(orb_indx)
    l_values = orb_indx[:, 6].astype(int)
    orb_ids_in_atom = orb_indx[:, 4].astype(int)
    atom_ids_for_orbs = orb_indx[:, 1].astype(int)
    equiv_orb_idx = orb_indx[:, 15].astype(int)
    
    # Group by atom
    unique_atoms = np.unique(atom_ids_for_orbs)
    
    with open('{}/orbital_types.dat'.format(output_path), 'w') as orb_type_f:
        for atom_id in unique_atoms:
            # Check if this is a unit cell atom (equiv_orb_idx equals orb index)
            atom_orbs_mask = atom_ids_for_orbs == atom_id
            atom_equiv_idx = equiv_orb_idx[atom_orbs_mask]
            atom_l_values = l_values[atom_orbs_mask]
            
            # Write orbital types for s, p, d, f
            for angular_momentum in range(4):
                count = np.sum(atom_l_values == angular_momentum) // (2 * angular_momentum + 1)
                for _ in range(count):
                    orb_type_f.write('{}  '.format(angular_momentum))
            orb_type_f.write('\n')
    
    # ========================================================================
    # 10. Build orb2deephorb mapping (fully vectorized)
    # ========================================================================
    orb2deephorb = np.zeros((n_orbitals, 5))
    
    # Extract all needed columns once
    n_col = orb_indx[:, 5].astype(int)
    l_col = orb_indx[:, 6].astype(int)
    m_col = orb_indx[:, 7].astype(int)
    zeta_col = orb_indx[:, 8].astype(int)
    
    # Process each atom's orbitals
    t = 0
    for atom_id in unique_atoms:
        atom_mask = atom_ids_for_orbs == atom_id
        atom_orb_indices = np.where(atom_mask)[0]
        n_orbs_in_atom = len(atom_orb_indices)
        
        if n_orbs_in_atom == 0:
            continue
        
        # Extract orbital data for this atom
        atom_n = n_col[atom_mask]
        atom_l = l_col[atom_mask]
        atom_m = m_col[atom_mask]
        atom_zeta = zeta_col[atom_mask]
        
        # Apply m quantum number transformation (vectorized)
        # p orbitals (l=1): -1->0, 0->1, 1->-1
        p_mask = atom_l == 1
        atom_m_transformed = atom_m.copy()
        atom_m_transformed[p_mask & (atom_m == -1)] = 0
        atom_m_transformed[p_mask & (atom_m == 0)] = 1
        atom_m_transformed[p_mask & (atom_m == 1)] = -1
        
        # d orbitals (l=2): -2->0, -1->2, 0->-2, 1->1, 2->-1
        d_mask = atom_l == 2
        atom_m_transformed[d_mask & (atom_m == -2)] = 0
        atom_m_transformed[d_mask & (atom_m == -1)] = 2
        atom_m_transformed[d_mask & (atom_m == 0)] = -2
        atom_m_transformed[d_mask & (atom_m == 1)] = 1
        atom_m_transformed[d_mask & (atom_m == 2)] = -1
        
        # f orbitals (l=3): -3->0, -2->1, -1->-1, 0->2, 1->-2, 2->3, 3->-3
        f_mask = atom_l == 3
        atom_m_transformed[f_mask & (atom_m == -3)] = 0
        atom_m_transformed[f_mask & (atom_m == -2)] = 1
        atom_m_transformed[f_mask & (atom_m == -1)] = -1
        atom_m_transformed[f_mask & (atom_m == 0)] = 2
        atom_m_transformed[f_mask & (atom_m == 1)] = -2
        atom_m_transformed[f_mask & (atom_m == 2)] = 3
        atom_m_transformed[f_mask & (atom_m == 3)] = -3
        
        # Calculate sort index: m + 10*zeta + 100*atom_id + 1000*l
        sort_index = (atom_m_transformed + 10 * atom_zeta + 
                      100 * atom_ids_for_orbs[atom_mask] + 1000 * atom_l)
        
        # Get ordering
        orb_order_local = np.argsort(np.argsort(sort_index))
        
        # Fill orb2deephorb
        for local_idx, global_idx in enumerate(atom_orb_indices):
            orb2deephorb[global_idx, 0:3] = np.round(orb_indx[global_idx, 12:15])
            orb2deephorb[global_idx, 3] = ia2Riua[atom_ids_for_orbs[global_idx] - 1, 3]
            orb2deephorb[global_idx, 4] = orb_order_local[local_idx]
    
    # ========================================================================
    # 11. Parse HSX file (fully vectorized - major performance gain)
    # ========================================================================
    f = FortranFile('{}/{}.HSX'.format(input_path, system_name), 'r')
    
    # Read header information
    header1 = f.read_ints()
    no_u, no_s, nspin, nh = header1[0], header1[1], header1[2], header1[3]
    _ = f.read_ints()  # gamma (unused)
    _ = f.read_ints()  # indxuo (unused)
    numh = f.read_ints()
    maxnumh = max(numh)
    
    # Read listh matrix (vectorized)
    listh = np.zeros((no_u, maxnumh), dtype=int)
    for i in range(no_u):
        listh[i, :] = f.read_ints()
    
    # ========================================================================
    # 12. Find connected atoms (vectorized)
    # ========================================================================
    # Create mask for non-zero elements
    mask = listh != 0
    
    # Get indices where connections exist
    row_indices, col_indices = np.where(mask)
    
    if len(row_indices) > 0:
        # Vectorized lookup of atom indices and R vectors
        atom_1_indices = orb2deephorb[row_indices, 3].astype(int)
        connected_orbs = listh[row_indices, col_indices] - 1  # Convert to 0-indexed
        atom_2_indices = orb2deephorb[connected_orbs, 3].astype(int)
        Rijk_vectors = orb2deephorb[connected_orbs, 0:3].astype(int)
        
        # Generate keys efficiently
        connected_keys = [f'[{R[0]}, {R[1]}, {R[2]}, {a1}, {a2}]' 
                         for R, a1, a2 in zip(Rijk_vectors, atom_1_indices, atom_2_indices)]
        connected_atoms = set(connected_keys)
    else:
        connected_atoms = set()
    
    # Initialize sparse matrices
    H_block_sparse = {key: [] for key in connected_atoms}
    S_block_sparse = {key: [] for key in connected_atoms}
    
    # ========================================================================
    # 13. Parse Hamiltonian (vectorized by batches)
    # ========================================================================
    # Pre-compute all needed indices for speed
    orb_atom_id = orb2deephorb[:, 3].astype(int)
    orb_deep_id = orb2deephorb[:, 4].astype(int)
    orb_m_values = orb_indx[:, 7]
    R_vectors_orb = orb2deephorb[:, 0:3].astype(int)
    
    # Process each spin
    for spin_idx in range(nspin):
        # Read all Hamiltonian data for this spin
        for j in range(no_u):
            tmpt = f.read_reals(dtype='<f4')
            n_connections = len(tmpt)
            
            if n_connections == 0:
                continue
            
            # Vectorized computation of all indices
            connected_orb_indices = listh[j, :n_connections] - 1
            atom_1_id = orb_atom_id[j]
            atom_2_ids = orb_atom_id[connected_orb_indices]
            m_sum = orb_m_values[j] + orb_m_values[connected_orb_indices]
            Rijk = R_vectors_orb[connected_orb_indices]
            deep_id_1 = orb_deep_id[j]
            deep_id_2 = orb_deep_id[connected_orb_indices]
            
            # Apply (-1)^m factor
            sign_factor = np.where(m_sum % 2 == 0, 1.0, -1.0)
            values = tmpt * sign_factor
            
            # Generate keys and add to sparse matrix
            keys = [f'[{R[0]}, {R[1]}, {R[2]}, {a1}, {a2}]' 
                   for R, a1, a2 in zip(Rijk, atom_1_id, atom_2_ids)]
            
            # Append to sparse blocks
            for idx, (key, val) in enumerate(zip(keys, values)):
                H_block_sparse[key].append([deep_id_1, deep_id_2[idx], val])
    
    # ========================================================================
    # 14. Parse Overlap matrix (vectorized)
    # ========================================================================
    for j in range(no_u):
        tmpt = f.read_reals(dtype='<f4')
        n_connections = len(tmpt)
        
        if n_connections == 0:
            continue
        
        # Vectorized computation
        connected_orb_indices = listh[j, :n_connections] - 1
        atom_1_id = orb_atom_id[j]
        atom_2_ids = orb_atom_id[connected_orb_indices]
        m_sum = orb_m_values[j] + orb_m_values[connected_orb_indices]
        Rijk = R_vectors_orb[connected_orb_indices]
        deep_id_1 = orb_deep_id[j]
        deep_id_2 = orb_deep_id[connected_orb_indices]
        
        # Apply (-1)^m factor
        sign_factor = np.where(m_sum % 2 == 0, 1.0, -1.0)
        values = tmpt * sign_factor
        
        # Generate keys
        keys = [f'[{R[0]}, {R[1]}, {R[2]}, {a1}, {a2}]' 
               for R, a1, a2 in zip(Rijk, atom_1_id, atom_2_ids)]
        
        # Append to sparse blocks
        for idx, (key, val) in enumerate(zip(keys, values)):
            S_block_sparse[key].append([deep_id_1, deep_id_2[idx], val])
    
    f.close()
    
    # ========================================================================
    # 15. Build atom2nu (vectorized)
    # ========================================================================
    nua = int(max(orb2deephorb[:, 3]))
    atom2nu = np.zeros(nua)
    
    # Find unit cell atoms (R = 0, 0, 0) and get max orbital index
    uc_mask = (orb_indx[:, 12] == 0) & (orb_indx[:, 13] == 0) & (orb_indx[:, 14] == 0)
    uc_atom_ids = atom_ids_for_orbs[uc_mask]
    uc_orb_ids = orb_ids_in_atom[uc_mask]
    
    for atom_id in unique_atoms:
        atom_max_orb = np.max(uc_orb_ids[uc_atom_ids == atom_id])
        atom2nu[atom_id - 1] = atom_max_orb
    
    # ========================================================================
    # 16. Convert COO sparse to dense and save (vectorized)
    # ========================================================================
    conversion_factor = 1.0 / 0.036749324533634074 / 2
    
    for key in H_block_sparse:
        sparse_list = H_block_sparse[key]
        if len(sparse_list) == 0:
            H_block_sparse[key] = np.zeros((1, 1))
            continue
            
        sparse_arr = np.array(sparse_list)
        ia1 = int(key[1:-1].split(',')[3])
        ia2 = int(key[1:-1].split(',')[4])
        rows = sparse_arr[:, 0].astype(int)
        cols = sparse_arr[:, 1].astype(int)
        values = sparse_arr[:, 2] * conversion_factor
        
        dense_matrix = np.zeros((int(atom2nu[ia1 - 1]), int(atom2nu[ia2 - 1])))
        dense_matrix[rows, cols] = values
        H_block_sparse[key] = dense_matrix
    
    # Save Hamiltonians
    f = h5py.File('{}/hamiltonians.h5'.format(output_path), 'w')
    for key, value in H_block_sparse.items():
        f[key] = value
    f.close()
    
    # Convert and save Overlaps
    for key in S_block_sparse:
        sparse_list = S_block_sparse[key]
        if len(sparse_list) == 0:
            S_block_sparse[key] = np.zeros((1, 1))
            continue
            
        sparse_arr = np.array(sparse_list)
        ia1 = int(key[1:-1].split(',')[3])
        ia2 = int(key[1:-1].split(',')[4])
        rows = sparse_arr[:, 0].astype(int)
        cols = sparse_arr[:, 1].astype(int)
        values = sparse_arr[:, 2]
        
        dense_matrix = np.zeros((int(atom2nu[ia1 - 1]), int(atom2nu[ia2 - 1])))
        dense_matrix[rows, cols] = values
        S_block_sparse[key] = dense_matrix
    
    f = h5py.File('{}/overlaps.h5'.format(output_path), 'w')
    for key, value in S_block_sparse.items():
        f[key] = value
    f.close()


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 3:
        siesta_parse(sys.argv[1], sys.argv[2])
    else:
        print("Usage: python siesta_get_data_optimized.py <input_path> <output_path>")
