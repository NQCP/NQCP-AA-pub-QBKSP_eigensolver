################################################################################
########### This file contains diatomic molecules model definition  ############
################################################################################

# Import libraries and packages
from qiskit_nature.second_q.transformers import FreezeCoreTransformer
from qiskit_nature.second_q.formats.molecule_info import MoleculeInfo
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.circuit.library import  HartreeFock
from qiskit_nature.second_q.drivers import PySCFDriver

def diatomic_molecule_Hamiltonian(symbols,dist, total_num_orbitals=None):
    '''
    Function to return the Hamiltonian of a diatomic molecule
    symbols: list of strings - the symbols of the atoms in the molecule
    dist: float - the distance between the atoms
    total_num_orbitals: int - the total number of spatial orbitals in the molecule to be used

    Returns:
            problem_reduced: the reduced eletronic problem
            qubit_op: SparsePauliOp - the Hamiltonian in the Pauli basis
            matrix_op: numpy array - the Hamiltonian in the matrix form
            coefs: list - the coefficients of the Hamiltonian
            paulis: list - the Pauli operators of the Hamiltonian
            hf_state: HartreeFock - the HartreeFock state of the molecule
    '''

    

    # the freeze-core reduction
    transformer = FreezeCoreTransformer()

    # Define the molecule
    molecule = MoleculeInfo(
        symbols=symbols,
        coords=[(0.0, 0.0, 0.0), (0.0, 0.0, dist)],
        multiplicity=1,  # = 2*spin + 1
        charge=0,

    )
    
    # Define the driver
    driver = PySCFDriver.from_molecule(molecule,basis='sto3g')

    problem = driver.run()

    if total_num_orbitals is None:
        total_num_orbitals = problem.num_spatial_orbitals
        
    transformer.prepare_active_space(molecule,total_num_orbitals)
    problem_reduced = transformer.transform(problem)
    mapper = ParityMapper(problem_reduced.num_particles)
    qubit_op = mapper.map(problem_reduced.second_q_ops()[0])

    matrix_op = qubit_op.to_matrix()
    coefs = qubit_op.coeffs
    paulis = qubit_op.paulis
    a,b =problem_reduced.num_particles
    hf_state = HartreeFock(num_spatial_orbitals=problem_reduced.num_spatial_orbitals,
                        num_particles=(a,b),
                        qubit_mapper=mapper)


    
    return problem_reduced, qubit_op, matrix_op, coefs, paulis, hf_state

