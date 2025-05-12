import pathlib
import sys

import matplotlib.pyplot as plt
import time
import qiskit.quantum_info
import scipy
import numpy as np
import qiskit

import qiskit_aer

hhl_qiskit_path = pathlib.Path(__file__).parent.parent.parent.parent
hhl_qiskit_path = hhl_qiskit_path.joinpath('QGPE').joinpath('QLinAlg').joinpath('HHL')
print(hhl_qiskit_path)
sys.path.append(str(hhl_qiskit_path))

from HHL_Qiskit import HHL_Qiskit

system_size = 4 # see cases below
device = 'CPU' # 'GPU' or 'CPU'
is_print_figs_circuits = False
is_print_figs_result = True
if system_size == 4:
    c_reg_size = 5 # 5, playaround determined this to work somewhat well
    hhl_instance = HHL_Qiskit(c_reg_size=c_reg_size)
    A,b,alpha = hhl_instance.init_Ab_poisson_first_order_FD(system_size=system_size)
    A_herm, b_herm, b_herm_normalized = hhl_instance.hermitianize_system()



print('Matrix A_herm:\n', hhl_instance.A_herm)
print('Eigvals of A_herm:\n', scipy.linalg.eigvals(hhl_instance.A_herm))

sol_classical_herm = hhl_instance.compute_classical_solutions()[2]
    

# Hamiltonian Parameter --> Defines the total time for the Hamiltonial simulations
#t_ = (1 - 2**(-1*c_reg_size))*n_time*np.pi/4
#t_ = np.pi/max(scipy.linalg.eigvals(A_herm).real)
#t_ = np.pi*3/4
print(hhl_instance.t_hamiltonian)
# Total number of circuit runs
shots = int(1e6)

# Construct the quantum circuit
qc = hhl_instance.construct_quantum_circuit()
qc_wo_meas = hhl_instance.construct_quantum_circuit(include_measurements=False, skip_state_prep=False)
qc_wo_meas_wo_prep = hhl_instance.construct_quantum_circuit(include_measurements=False, skip_state_prep=True)
if is_print_figs_circuits:
    qc.draw('mpl')
    plt.savefig(f'HHL_circuit_Qiskit_{system_size}x{system_size}.png')
print('Available devices to run simulator on:', qiskit_aer.AerSimulator(method='statevector').available_devices())
sim = qiskit_aer.AerSimulator(method='statevector', device=device)
#instructions = [so.name for so in sim.operations if isinstance(so, qiskit.circuit.Instruction)]
#print('Instructions:', instructions)
qc_transpiled = qiskit.transpile(circuits=qc,backend=sim, optimization_level=0)#, basis_gates=instructions)
#print(qc_transpiled.data)

tic = time.time()
job = sim.run(qc_transpiled, shots=shots)
results = job.result()
counts = results.get_counts()
toc_compute_samples = time.time()
print('Time to compute samples:', toc_compute_samples - tic)

y = np.zeros(len(b_herm))
for key in counts:
    counts[key] *= (1/shots)
    counts[key] = np.sqrt(counts[key])
    #print(key, np.sqrt(counts[key]))
for key in counts:
    if key[-1] == '1':
        pos = int(key[0:hhl_instance.b_reg_size],2)
        y[pos] +=counts[key]
        #print('key:', key)
        #print('pos:', pos)


#print('counts:', counts)
print('Quantum solution, ratio elem 0/1:', y, y[0]/y[1])
print('Quantum solution normalized, norm:', y/np.linalg.norm(y), np.linalg.norm(y))
print('Classical solution, ratio elem 0/1:', sol_classical_herm, sol_classical_herm[0]/sol_classical_herm[1])
print('Classical solution normalized, norm:', sol_classical_herm/np.linalg.norm(sol_classical_herm), np.linalg.norm(sol_classical_herm))

    
if is_print_figs_circuits:
    qc_transpiled.draw('mpl')
    plt.savefig(f'HHL_circuit_transpiled_Qiskit_{system_size}x{system_size}.png')

if is_print_figs_result:
    fig, ax = plt.subplots()
    ax.plot(y)
    fig.savefig(f'HHL_result_Qiskit_{system_size}x{system_size}.png')

num_qubits = hhl_instance.c_reg_size + hhl_instance.b_reg_size + 1
bitstring_le = [format(i, '0' + str(num_qubits) + 'b') for i in range(2**num_qubits)]
#print('bitstring_le:', bitstring_le)
print(qc_wo_meas.draw())
op_wo_meas = qiskit.quantum_info.Operator(qc_wo_meas).to_matrix()
op_wo_meas_wo_prep = qiskit.quantum_info.Operator(qc_wo_meas_wo_prep).to_matrix()
#print(op)
z = np.zeros(op_wo_meas.shape[0])
z[0] = 1.
z = op_wo_meas @ z
state_dict = {key:val for key, val in zip(bitstring_le, z)}
#print('state_dict:', state_dict)

y = np.zeros(len(b_herm))
for key in state_dict:
    if key[-1] == '1':
        pos = int(key[0:hhl_instance.b_reg_size],2)
        y[pos] +=np.abs(state_dict[key])
        #print('key:', key)
        #print('pos:', pos)
print('Quantum solution, ratio elem 0/1:', y, y[0]/y[1])





H = 1/np.sqrt(2)*np.array([[1, 1], [1, -1]])
print('H:', H)
HH = np.kron(H, H)
print('HH:', HH)
eye = np.eye(2)
opprep = 1. * H
for i in range(hhl_instance.b_reg_size-1): #-1 because start with H
    opprep = np.kron(opprep, H)
for i in range(hhl_instance.c_reg_size+1): #+1 because ancilla
    opprep = np.kron(opprep, eye)
print('opprep:', opprep)

reconstructed = op_wo_meas_wo_prep @ opprep
print(np.allclose(op_wo_meas, reconstructed))



z = np.zeros(op_wo_meas.shape[0])
z[0] = 1.
z2 = opprep @ z
print('prep state:', z2)
z = op_wo_meas_wo_prep @ opprep @ z
state_dict = {key:val for key, val in zip(bitstring_le, z)}
#print('state_dict:', state_dict)

y = np.zeros(len(b_herm))
for key in state_dict:
    if key[-1] == '1':
        pos = int(key[0:hhl_instance.b_reg_size],2)
        y[pos] +=np.abs(state_dict[key])
        #print('key:', key)
        #print('pos:', pos)
print('Quantum solution, ratio elem 0/1:', y, y[0]/y[1])

np.save(f'HHL_result_Qiskit_{system_size}x{system_size}_c{c_reg_size}.npy', op_wo_meas_wo_prep)


# Qubits
b_reg_size = int(np.log2(system_size))
anc = qiskit.QuantumRegister(1,'anc')
c_reg = qiskit.QuantumRegister(c_reg_size, 'c_reg')
b_reg = qiskit.QuantumRegister(b_reg_size, 'b_reg')
# Classical bits
b_reg_cl = qiskit.ClassicalRegister(b_reg_size, 'b_reg_cl')

qc_plot = qiskit.QuantumCircuit(anc, c_reg, b_reg, b_reg_cl)
qc_plot.barrier(label='state prep')
#qc_plot.h(b_reg[0])
HGate = qiskit.circuit.library.HGate().to_mutable()
HGate.label = '$H_{eI}$'
qc_plot.append(HGate, [b_reg[0]])
qc_plot.h(b_reg[1])
qc_plot.barrier(label='HHL')
op_hhl_wo_meas_wo_prep = qiskit.quantum_info.Operator(qc_wo_meas_wo_prep).to_instruction()
op_hhl_wo_meas_wo_prep.label = '   HHL   '
qargs = [anc]
qargs.extend([c_reg[i] for i in range(c_reg_size)])
qargs.extend([b_reg[i] for i in range(b_reg_size)])
qc_plot.append(op_hhl_wo_meas_wo_prep, qargs=qargs)
qc_plot.barrier(label='measure')
qc_plot.measure(b_reg, b_reg_cl)
#qc_plot.measure(1+c_reg_size, 0)
#qc_plot.measure(1+c_reg_size+1, 0+1)
print(qc_plot)

import matplotlib.pyplot as plt
qc_plot.draw('mpl')
plt.savefig(f'HHL_circuit_Qiskit_{system_size}x{system_size}_plot.png',
            dpi=600, bbox_inches='tight')
#b_reg[:], b_reg_cl[:])