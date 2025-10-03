

from .openqasm3_stdgates_inc import *

OPENQASM3_TO_QSYMPY = {
    'U': U_Gate,
    'p': P_Gate,
    'x': Pauli_X_Gate,
    'y': Pauli_Y_Gate,
    'z': Pauli_Z_Gate,
    'h': Hadamard_Gate,
    's': S_Gate,
    'sdg': Sdg_Gate,
    't': T_Gate,
    'tdg': Tdg_Gate,
    'sx': SX_Gate,
    'rx': RX_Gate,
    'ry': RY_Gate,
    'rz': RZ_Gate,
    'cx': CX_Gate,
    'cy': CY_Gate,
    'cz': CZ_Gate,
    'cp': CP_Gate,
    'crx': CRX_Gate,
    'cry': CRY_Gate,
    'crz': CRZ_Gate,
    'ch': CH_Gate,
    'swap': SWAP_Gate,
    'ccx': CCX_Gate,
    'cswap': None,
    'cu': CU_Gate
    }