

from .openqasm3_stdgates_inc import *

OPENQASM3_TO_QSYMPY = {
    'p': P_Gate,
    'x': Pauli_X_Gate,
    'y': Pauli_Y_Gate,
    'z': Pauli_Z_Gate,
    'h': Hadamard_Gate,
    's': None,
    'sdg': None,
    't': None,
    'tdg': None,
    'sx': None,
    'rx': None,
    'ry': None,
    'rz': None,
    'cx': CX_Gate,
    'cy': None,
    'cz': None,
    'cp': None,
    'crx': None,
    'cry': None,
    'crz': None,
    'ch': None,
    'swap': None,
    'ccx': CCX_Gate,
    'cswap': None,
    'cu': CU_Gate
    }