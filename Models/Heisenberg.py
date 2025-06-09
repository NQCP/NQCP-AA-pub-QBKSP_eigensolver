###################################################################################
############# This file contains the Hubbard model definition  ####################
###################################################################################

# Import Hubbard model from qiskit
from qiskit_nature.second_q.hamiltonians import lattices, HeisenbergModel
from qiskit.quantum_info import SparsePauliOp



def Heisenberg_Hamiltonian(n_sites):
    '''	
    Function to return the Heisenberg Hamiltonian
    n_sites: number of sites

    Returns:
            pauli_op: SparsePauliOp - the Hamiltonian in the Pauli basis
            matrix_op: numpy array - the Hamiltonian in the matrix form
    '''
    # Define a 1D line lattice with n_sites
    line_lattice = lattices.LineLattice(num_nodes=n_sites)
    heisenberg_model = HeisenbergModel(
        line_lattice
    )
    
    qubit_op = heisenberg_model.second_q_op()

    matrix_op = qubit_op.to_matrix()

    pauli_op = SparsePauliOp.from_operator(qubit_op)
    
    return pauli_op, matrix_op



    