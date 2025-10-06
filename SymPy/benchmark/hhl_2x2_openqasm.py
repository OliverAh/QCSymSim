import sys
import pathlib
qsympy_dir = pathlib.Path(__file__).parent.parent
print('qsympy_dir:', qsympy_dir)
sys.path.append(str(qsympy_dir))

import QSymPy
import numpy as np
import sympy as sp
import time

timestamps = {}

timestamps['start'] = time.time()
print(timestamps['start'])

qc = QSymPy.openqasm3_to_qc('../HHL_circuit_Qiskit_2x2.qasm3')

timestamps['finish_build_qc'] = time.time()
print(timestamps['finish_build_qc'])
print('Time to parse QASM3 file and create QuantumCircuit:', timestamps['finish_build_qc'] - timestamps['start'])

qc.assemble_symbolic_unitary(use_alternative_repr=True, replace_symbolic_zeros_and_ones=True)

timestamps['finish_assemble_unitary'] = time.time()
print(timestamps['finish_assemble_unitary'])
print('Time to assemble symbolic unitary:', timestamps['finish_assemble_unitary'] - timestamps['finish_build_qc'])

state = QSymPy.QuantumState(num_qubits=len(qc.qubits))
state.set_state({k: 0 for k in range(len(qc.qubits))})
print('Initial state:')
print(state.state)

timestamps['finished_build_state'] = time.time()
print(timestamps['finished_build_state'])
print('Time to build initial state:', timestamps['finished_build_state'] - timestamps['finish_assemble_unitary'])

state2 = qc.unitary @ state.state

timestamps['finished_apply_unitary'] = time.time()
print(timestamps['finished_apply_unitary'])
print('Time to apply unitary:', timestamps['finished_apply_unitary'] - timestamps['finished_build_state'])

a = qc.gate_collection.collections['CU'][0].matrix22_t_alt[1][1][0,0]
var_for_grad = a.args[0].args[0].args[1]
print('Variable for gradient:', var_for_grad)

grad = sp.diff(state2[32], var_for_grad)

timestamps['finished_compute_grad_state2_32'] = time.time()
print(timestamps['finished_compute_grad_state2_32'])
print('Time to compute gradient of state2[32]:', timestamps['finished_compute_grad_state2_32'] - timestamps['finished_apply_unitary'])
print('Value of gradient of state2[32]:\n', grad)