################################################################################
################# This file contains the QBKSP algorithm  ######################
################################################################################

# Imports
from scipy import linalg
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.synthesis.evolution import LieTrotter, SuzukiTrotter
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, transpile
import numpy as np
from qiskit_aer import Aer
from qiskit.quantum_info import SparsePauliOp
import numpy as np
import sys
sys.path.append('Models')

# Defining global variables
NUM_SHOTS = 100000
BACKEND = Aer.get_backend('qasm_simulator')

#####################################################################################
############# Numerical Quantum Block Krylov Subspace Projection (QBKSP) ############
#####################################################################################

def create_T_S_block_numerical(init_states, unitary_operator, max_iter=None, threshold=None):
    """
    Numerical Quantum Block Krylov Subspace Projection (QBKSP)

    Args:
        init_states: ndarray (B, N) - block of reference states
        unitary_operator: ndarray (N, N)
        max_iter: int
        threshold: function to check convergence - if threshold output is 1 convergence is reached

    Returns:
        T: ndarray (B * max_iter, B * max_iter)
        S: ndarray (B * max_iter, B * max_iter)
    """
      
    B, N = init_states.shape
    if max_iter is None:
        max_iter = N // B

    # Normalize initial block
    p0 = []
    for b in range(B):
        norm = np.linalg.norm(init_states[b])
        p0.append(init_states[b] / norm)

    # Precompute the basis states
    p = p0.copy()
    for i in range(max_iter+1):
        for bi in range(B):
            idx_i = i * B + bi
            p.append(unitary_operator @ p[idx_i])
    
    # Precompute all required inner products
    array = []
    for b in range(B):
        for b2 in range(B):
                if b <= b2:
                    array.append([np.vdot(p[b], p[b2])])
                else:
                    array.append([np.conj(array[b2*B+b][0])])
     
   
    T = np.zeros((B * max_iter, B * max_iter), dtype=np.complex128)	
    S = np.zeros((B * max_iter, B * max_iter), dtype=np.complex128)
    for i in range(max_iter):
        for bi in range(B):
            idx_i = i * B + bi
            for bj in range(B):
                if bi <= bj:
                    array[bi*B+bj].append(np.vdot(p[bi], p[(i+1)*B+bj]))
                else:
                    array[bi*B+bj].append(array[bj*B+bi][i+1])
                
                for j in range(i + 1):
                        idx_j = j * B + bj
                        T[idx_j, idx_i] = array[bi*B+bj][i-j+1]
                 
                        S[idx_j, idx_i] = array[bj * B + bi][i - j]
                        if i!=j:
                            T[idx_i, idx_j] = np.conj(array[bi*B+bj][i-j-1])
                            S[idx_i, idx_j] = np.conj(S[idx_j, idx_i]) 

        # Check convergence if threshold is not None:
        if  threshold is not None:
            if threshold(T,S) == 1:
                return T, S                             
    return T, S

#####################################################################################
############## Quantum Block Krylov Subspace Projection (QBKSP) #####################
#####################################################################################


def quantum_expectation_value_aux(size,initial_stateleft, Initial_stateright,Href,time,var='R', backend = BACKEND, shots = NUM_SHOTS, reps = 1):
    '''
    Function to simulate the Hadamard test circuit
    input time: float - time
    input var: str - variable to measure, 'R' for real part, 'I' for imaginary part
    input initial_stateleft: numpy array - initial state of the system (left of the expectation value)
    input Initial_stateright: numpy array - initial state of the system (right of the expectation value)
    input system_size: int - number of qubits in the system
    input Href: qubit op - the model circuit to be simulated
    input shots: int - number of shots
    input backend: the backend to run the simulation
    input reps: int - number of repetitions for the Trotter decomposition
    input size: int - number of qubits in the system

    output result: float - the result of the Hadamard test
    '''
    print(time)
    # Initialize the circuit
    system = QuantumRegister(size)
    ancilla = QuantumRegister(1)

    CLASSICAL = ClassicalRegister(1) # 1 classical bit to measure the ancilla
    QC = QuantumCircuit(ancilla, system, CLASSICAL)

    ###########################################
    ######### Initial State Preparation ########
    ###########################################
    QC.h(ancilla)

 
    control_right = Initial_stateright.control(1)
  
    QC.append(control_right, ancilla[:] + system[:])
    if time!= 0:
        ##############################################
        ######### Hamiltonian Simulation #############
        ##############################################
        trotter = SuzukiTrotter(reps=reps, order=2)
        unitary = PauliEvolutionGate(Href, time = time, synthesis=trotter)
        controlled_unitary = unitary.control(1)
        QC.append(controlled_unitary, ancilla[:]+ system[:])
       
        
   
    ##############################################
    control_left = initial_stateleft.inverse().control(1)
    #control_left = initial_stateleft.control(1)
    QC.append(control_left, ancilla[:] + system[:])
    ##############################################
    ############ Measurement and Output ##########
    ##############################################
    if var == 'I':
        # s dagger gate
        QC.sdg(ancilla)
    
    QC.h(ancilla)
    QC.measure(ancilla,0)

   
    
    job = transpile(QC, backend)
    result = backend.run(job, shots=shots).result()
    counts = result.get_counts()
    
    
    if '0' not in counts.keys():
        return -1
    if '1' not in counts.keys():
        return 1
    return (counts['0'] - counts['1'])/shots

def quantum_expectation_values(size,init_left,init_right,unitary_operator,i,t=1,ntrotter=1,shots = NUM_SHOTS):
    '''
    Function to simulate the Hadamard test circuit <psi_left|U|psi_right> for the given time t
    input t: float - time
    input init_left: numpy array - initial state of the system (left of the expectation value) <psi_left|
    input init_right: numpy array - initial state of the system (right of the expectation value) |psi_right>
    input size: int - number of qubits in the system 
    input unitary_operator: qubit op - the model circuit to be simulated U
    input i: int - index of the krylov iteration
    input ntrotter: int - number of repetitions for the Trotter decomposition
    input shots: int - number of shots
    
    output result: complex - expectation value of the operator - <psi_left|U|psi_right>
    '''
    real = quantum_expectation_value_aux(size,init_left,init_right,unitary_operator,t*i,var='R',backend=BACKEND,reps=ntrotter,shots=shots)
    imag = quantum_expectation_value_aux(size,init_left,init_right,unitary_operator,t*i,var='I',backend=BACKEND,reps=ntrotter,shots=shots)
    
    return real+1j*imag

def create_T_S_block_quantum(size, unitary_operator, qcs_init, max_iter=None,ntrotter=1,t=1,NUM_SHOTS=NUM_SHOTS,threshold=None):
    """
    Block Lanczos T matrix (element-wise, precomputed)

    Args:
        size: int - number of qubits in the system
        unitary_operator: qubit op - the model circuit to be simulated
        qcs_init: list of QuantumCircuits - block of reference states
        max_iter: int - maximum number of iterations
        ntrotter: int - number of repetitions for the Trotter decomposition
        t: float - time
        NUM_SHOTS: int - number of shots for the simulation
        threshold: function to check convergence - if threshold output is 1 convergence is reached
    Returns:
        T: ndarray (B * max_iter, B * max_iter) - T matrix
        S: ndarray (B * max_iter, B * max_iter) - S matrix
    """
   
    B = len(qcs_init)
    
    # Precompute all required inner products
    array = []
    for b in range(B):
        init_left = qcs_init[b]
        for b2 in range(B):
            init_right = qcs_init[b2]
            inner = []
            for i in range(1):
                    val = quantum_expectation_values(size,init_left,init_right,unitary_operator,i,t=t,ntrotter=ntrotter*i,shots=NUM_SHOTS)
                    inner.append(val)
                
            array.append(inner)
      
    T = np.zeros((B * max_iter, B * max_iter), dtype=np.complex128)	
    S = np.zeros((B * max_iter, B * max_iter), dtype=np.complex128)
    for i in range(max_iter):
        for bi in range(B):
            idx_i = i * B + bi
            for bj in range(B):
                array[bi*B+bj].append(quantum_expectation_values(size,qcs_init[bi],qcs_init[bj],unitary_operator,i+1,t=t,ntrotter=ntrotter*(i+1),shots=NUM_SHOTS))
                for j in range(i + 1):           
                    idx_j = j * B + bj
                    T[idx_j, idx_i] = array[bi*B+bj][i-j+1]
                    S[idx_j, idx_i] = array[bj * B + bi][i - j]
                    if i!=j:
                        T[idx_i, idx_j] = np.conj(array[bi*B+bj][i-j-1])
                        S[idx_i, idx_j] = np.conj(S[idx_j, idx_i])

        # Check convergence if threshold is not None:
        if  threshold is not None:
            if threshold(T,S) == 1:
                return T, S
                      
    return T, S
