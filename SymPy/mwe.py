import QSymPy as qp
import sympy as sp

qc = qp.QuantumCircuit(num_qubits=3)
qc.add_gate(name='H', qubits_t=[0], step=1)
qc.add_gate(name='CX', qubits_t=[1], qubits_c=[0], step=2)
qc.add_gate(name='CX', qubits_t=[2], qubits_c=[0], step=2)

qc.assemble_symbolic_unitary()

qs = qp.QuantumState(num_qubits=3)

qs2 = qc.unitary @ qs.state

print(qs2)

print(sp.diff(qs2, qc.gate_collection.collections['CX'][0].matrix[1, 1]))