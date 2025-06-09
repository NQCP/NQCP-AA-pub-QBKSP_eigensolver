################################################################################
#################### This file contains auxiliary functions  ###################
################################################################################

# Imports
import numpy as np

def find_nearest(array, value):
    """
    Find the nearest value in an array to a given value.
    Args:
        array: numpy array
        value: float - the value to find the nearest to
    Returns:
        nearest_value: float - the nearest value in the array
        idx: int - the index of the nearest value in the array
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

def generate_random_state_with_overlap(target_state: np.ndarray, target_overlap: float) -> np.ndarray:
    """
    Generate a random quantum state with specified overlap with target state.
    
    Args:
        target_state (np.ndarray): Target state vector (must be normalized)
        target_overlap (float): Desired overlap between 0 and 1
        
    Returns:
        np.ndarray: New random state with specified overlap
    """
    if not 0 <= target_overlap <= 1:
        raise ValueError("Target overlap must be between 0 and 1")
        
    if abs(np.linalg.norm(target_state) - 1) > 1e-10:
        raise ValueError("Target state must be normalized")
        
    # Generate random perpendicular component
    random_state = np.random.randn(len(target_state)) 
    random_state = random_state - np.vdot(target_state, random_state) * target_state  # Make orthogonal
    random_state /= np.linalg.norm(random_state)  # Normalize
    
    # Combine target and perpendicular components with proper weights
    guess_state = np.sqrt(target_overlap) * target_state + np.sqrt(1 - target_overlap) * random_state

    return guess_state