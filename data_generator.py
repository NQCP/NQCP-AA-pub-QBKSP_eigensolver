from qbksp import create_T_S_block_numerical, create_T_S_block_quantum
from utils import generate_random_state_with_overlap
import pickle as pkl
import numpy as np
from scipy import linalg
from Models.hubbard import Hubbard_Hamiltonian
from Models.Heisenberg import Heisenberg_Hamiltonian
from Models.molecules import diatomic_molecule_Hamiltonian
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
from qiskit.circuit.library import StatePreparation
from qiskit_nature.second_q.mappers import ParityMapper

# --- Path to save data ---
Data_path = 'ADD YOUR PATH HERE'


# --- Helper Function ---
def load_pickle(path):
    with open(path, 'rb') as f:
        return pkl.load(f)
    

####################################################################################
######################### Classical expectation ####################################
####################################################################################

# Heisenberg
for sites in [2,3,4,5,6,7,8,9,10,11]:
    _,randm = Heisenberg_Hamiltonian(sites)
    norm = np.linalg.norm(randm,2)
    rand = randm/norm # Normalize it 
    eigvals, eigvecs = linalg.eigh(randm) 
    eigvals = np.unique(np.sort(eigvals))

    dic = {'eigvals': eigvals, 'eigvecs': eigvecs,'norm': norm}
    with open(f'{Data_path}/eigvals_eigvecs_norm_{sites}_heisenberg_line_3evol.pkl', 'wb') as f:
       pkl.dump(dic, f)
print('Heisenberg Hamiltonian eigvals and eigvecs saved')

#####################################################################################
################################ Numerical #########################################
#####################################################################################

################################
#### Heisenberg Hamiltonian ####
################################

for sites in [2,3,4,5,6,7,8,9,10,11]:
  size = 2**sites
  for overlap in [0.1,0.3,0.5,0.7,0.9]:
      for t in [1,2,3]:
        _,randm = Heisenberg_Hamiltonian(sites)
        norm = np.linalg.norm(randm,2)
        #print(norm)
        rand = randm/norm # Normalize it
        exp = linalg.expm(-1j *rand*t)  # Compute the unitary time evolution operator
        dic = load_pickle(f'{Data_path}/eigvals_eigvecs_norm_{sites}_heisenberg_line_3evol.pkl')
        eigvals = dic['eigvals']
        eigvecs = dic['eigvecs']
        np.random.seed(1)

        init_state = generate_random_state_with_overlap(eigvecs[:,0]/linalg.norm(eigvecs[:,0],2),overlap)
        
        state = np.array([init_state])
        T, S = create_T_S_block_numerical(state, exp, 200)
        
        file_t = f'{Data_path}/T_{sites}_heisenberg_line_{t}evol_random_overlap{overlap}.pkl'
        file_s = f'{Data_path}/S_{sites}_heisenberg_line_{t}evol_random_overlap{overlap}.pkl'
        with open(file_t, 'wb') as f:
            pkl.dump(T, f)
        with open(file_s, 'wb') as f:
            pkl.dump(S, f)

        init_state2 = generate_random_state_with_overlap(eigvecs[:,1]/linalg.norm(eigvecs[:,1],2),overlap)
        init_state3 = generate_random_state_with_overlap(eigvecs[:,2]/linalg.norm(eigvecs[:,2],2),overlap)
        T2, S2 = create_T_S_block_numerical(np.array([init_state, init_state2, init_state3]), exp, 200)

        file_t2 = f'{Data_path}/T_{sites}_heisenberg_line_3evol_random_overlap{overlap}_2nd{overlap}_3rd{overlap}.pkl'
        file_s2 = f'{Data_path}/S_{sites}_heisenberg_line_3evol_random_overlap{overlap}_2nd{overlap}_3rd{overlap}.pkl'
        with open(file_t2, 'wb') as f:
            pkl.dump(T2, f)
        with open(file_s2, 'wb') as f:
            pkl.dump(S2, f)

        state = np.array([init_state, init_state2])
        T3, S3 = create_T_S_block_quantum(state, exp, 200)
        
        file_s3 = f'{Data_path}/S_{sites}_heisenberg_line_{t}evol_random_overlap{overlap}_2nd{overlap}.pkl'
        file_t3 = f'{Data_path}/T_{sites}_heisenberg_line_{t}evol_random_overlap{overlap}_2nd{overlap}.pkl'
        with open(file_t3, 'wb') as f:
            pkl.dump(T, f)
        with open(file_s3, 'wb') as f:
            pkl.dump(S, f)
print('Heisenberg Hamiltonian T and S saved')

#####################################################################################
################################ Qiskit Simulations #################################
#####################################################################################

#############################
####### LiH molecule ########
#############################
problem_reduced, qubit_op, matrix_op, coefs, paulis, hf_state = diatomic_molecule_Hamiltonian(['Li', 'H'], 1.6)

# Retrieve the dipole moment property
dipole_property = problem_reduced.properties.electronic_dipole_moment


dipole_operators = dipole_property.second_q_ops()
dip_x = dipole_operators['XDipole']
dip_y = dipole_operators['YDipole']
dip_z = dipole_operators['ZDipole']
    # Initialize a qubit mapper (e.g., Jordan-Wigner)
mapper = ParityMapper(problem_reduced.num_particles)

qubit_op_dipole = mapper.map(dip_x)
qubit_op_dipolem = qubit_op_dipole.to_matrix()
excited_state = qubit_op_dipolem@Statevector(hf_state).data
qubit_op_y = mapper.map(dip_y)
qubit_op_ym= qubit_op_y.to_matrix()

excited_state_y = qubit_op_ym@Statevector(hf_state).data

qubit_op_z = mapper.map(dip_z)
qubit_op_zm = qubit_op_z.to_matrix()
excited_state_z = qubit_op_zm@Statevector(hf_state).data

QC_HF = QuantumCircuit(hf_state.num_qubits)
qcinit = StatePreparation(Statevector(hf_state))
QC_HF.append(qcinit, range(hf_state.num_qubits))
QC_HF.draw()


QC_HFX = QuantumCircuit(hf_state.num_qubits)
qcinit = StatePreparation(excited_state/linalg.norm(excited_state,2))
QC_HFX.append(qcinit, range(hf_state.num_qubits))
QC_HFX.draw()

QC_HFY = QuantumCircuit(hf_state.num_qubits)
qcinit = StatePreparation(excited_state_y/linalg.norm(excited_state_y,2))
QC_HFY.append(qcinit, range(hf_state.num_qubits))

QC_HFZ = QuantumCircuit(hf_state.num_qubits)
qcinit = StatePreparation(excited_state_z/linalg.norm(excited_state_z,2))
QC_HFZ.append(qcinit, range(hf_state.num_qubits))
QC_HFY.draw()

state = [QC_HF,QC_HFX,QC_HFY,QC_HFZ]
ensemble = 1
for n_trotter in [15,10,5]:
        for shots in [100000,1000000]:
    
            T, S = create_T_S_block_quantum(hf_state.num_qubits,qubit_op, state, max_iter=15,ntrotter=n_trotter,t=1,NUM_SHOTS=shots)
            # save T and S
            with open(f'{Data_path}/QISKIT_T_LiH_2ndntrotter{n_trotter}_1evol_hf_{shots}_ensemble{ensemble}vf15_ref4.pkl', 'wb') as f:
                pkl.dump(T, f)
            with open(f'{Data_path}/QISKIT_S_LiH_2ndntrotter{n_trotter}_1evol_hf_{shots}_ensemble{ensemble}vf15_ref4.pkl', 'wb') as f:
                pkl.dump(S, f)
            print('T and S saved', n_trotter, shots)

state = [QC_HF]
for n_trotter in [15,10,5]:
        for shots in [100000,1000000]:
    
            T, S = create_T_S_block_quantum(hf_state.num_qubits,qubit_op, state, max_iter=15,ntrotter=n_trotter,t=1,NUM_SHOTS=shots)
            # save T and S
            with open(f'{Data_path}/QISKIT_T_LiH_2ndntrotter{n_trotter}_1evol_hf_{shots}_ensemble{ensemble}vf15_ref1.pkl', 'wb') as f:
                pkl.dump(T, f)
            with open(f'{Data_path}/QISKIT_S_LiH_2ndntrotter{n_trotter}_1evol_hf_{shots}_ensemble{ensemble}vf15_ref1.pkl', 'wb') as f:
                pkl.dump(S, f)
            print('T and S saved', n_trotter, shots)
