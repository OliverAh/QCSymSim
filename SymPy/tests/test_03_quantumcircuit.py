import pathlib
import sys
import sympy as sp
import numpy as np
sys.path.append(str(pathlib.Path(__file__).parent.parent))

import pytest
import QSymPy as qp

class TestAddGateOrBarrier:
    def test_add_identity_gate(self):
        qc = qp.QuantumCircuit(num_qubits=2)
        qc.add_gate(name='I', qubits_t=[0], step=0)

        assert len(qc.gate_collection.collections['I']) == 1
        assert isinstance(qc.gate_collection.collections['I'][0], qp.Identity_Gate)
        assert qc.steps[0][0].name == 'I'

    def test_add_pauli_x_gate(self):
        qc = qp.QuantumCircuit(num_qubits=2)
        qc.add_gate(name='X', qubits_t=[1], step=1)

        assert len(qc.gate_collection.collections['X']) == 1
        assert isinstance(qc.gate_collection.collections['X'][0], qp.Pauli_X_Gate)
        assert qc.steps[1][0].name == 'X'

    def test_add_pauli_y_gate(self):
        qc = qp.QuantumCircuit(num_qubits=2)
        qc.add_gate(name='Y', qubits_t=[0], step=2)

        assert len(qc.gate_collection.collections['Y']) == 1
        assert isinstance(qc.gate_collection.collections['Y'][0], qp.Pauli_Y_Gate)
        assert qc.steps[2][0].name == 'Y'

    def test_add_pauli_z_gate(self):
        qc = qp.QuantumCircuit(num_qubits=2)
        qc.add_gate(name='Z', qubits_t=[1], step=3)

        assert len(qc.gate_collection.collections['Z']) == 1
        assert isinstance(qc.gate_collection.collections['Z'][0], qp.Pauli_Z_Gate)
        assert qc.steps[3][0].name == 'Z'

    def test_add_hadamard_gate(self):
        qc = qp.QuantumCircuit(num_qubits=2)
        qc.add_gate(name='H', qubits_t=[0], step=4)

        assert len(qc.gate_collection.collections['H']) == 1
        assert isinstance(qc.gate_collection.collections['H'][0], qp.Hadamard_Gate)
        assert qc.steps[4][0].name == 'H'

    def test_add_cx_gate(self):
        qc = qp.QuantumCircuit(num_qubits=2)
        qc.add_gate(name='CX', qubits_t=[1], qubits_c=[0], step=5)

        assert len(qc.gate_collection.collections['CX']) == 1
        assert isinstance(qc.gate_collection.collections['CX'][0], qp.CX_Gate)
        assert qc.steps[5][0].name == 'CX'

    def test_add_barrier(self):
        qc = qp.QuantumCircuit(num_qubits=2)
        qc.add_barrier()

        assert len(qc.barrier_collection) == 1
        assert isinstance(qc.barrier_collection[0], qp.Barrier)
        assert qc.barrier_collection[0].step == 0


        qc = qp.QuantumCircuit(num_qubits=2)
        qc.add_gate(name='I', qubits_t=[0], step=0)
        qc.add_gate(name='I', qubits_t=[0], step=1)
        qc.add_barrier()

        assert len(qc.barrier_collection) == 1
        assert qc.barrier_collection[0].step == 1

        qc = qp.QuantumCircuit(num_qubits=2)
        qc.add_gate(name='I', qubits_t=[0], step=0)
        qc.add_barrier(step=5)

        qc.add_gate(name='I', qubits_t=[0], step=0)
        assert len(qc.barrier_collection) == 1
        assert qc.barrier_collection[0].step == 5

    def test_add_gate_invalid_name(self):
        qc = qp.QuantumCircuit(num_qubits=2)
        with pytest.raises(ValueError, match="Unknown gate name"):
            qc.add_gate(name='INVALID', qubits_t=[0], step=0)

    def test_add_gate_invalid_gate_object(self):
        qc = qp.QuantumCircuit(num_qubits=2)
        with pytest.raises(ValueError, match="Providing a gate object is not yet implemented."):
            qc.add_gate(name=None, gate=qp.Identity_Gate(qubits_t=[0]))

    def test_add_gate_missing_qubits(self):
        qc = qp.QuantumCircuit(num_qubits=2)
        with pytest.raises(ValueError):
            qc.add_gate(name='X', step=0)

    def test_add_gate_control_and_target_same(self):
        qc = qp.QuantumCircuit(num_qubits=2)
        with pytest.raises(AssertionError, match="Control and target qubit must be different."):
            qc.add_gate(name='CX', qubits_t=[0], qubits_c=[0], step=0)

    def test_add_gate_multiple_target_qubits(self):
        qc = qp.QuantumCircuit(num_qubits=2)
        with pytest.raises(AssertionError, match="CX gate can only be applied to a single target qubit."):
            qc.add_gate(name='CX', qubits_t=[0, 1], qubits_c=[1], step=0)

    def test_add_gate_multiple_control_qubits(self):
        qc = qp.QuantumCircuit(num_qubits=2)
        with pytest.raises(AssertionError, match="CX gate can only be applied to a single control qubit."):
            qc.add_gate(name='CX', qubits_t=[0], qubits_c=[1, 0], step=0)

class TestAssemble_symbolic_unitary_2qubits:
    def test_assemble_symbolic_unitary_2qubits(self):
        gates_to_be_tested = ['I', 'X', 'Y', 'Z', 'H']
        for g0 in gates_to_be_tested:
            for g1 in gates_to_be_tested:
                qc = qp.QuantumCircuit(num_qubits=2)
                qc.add_gate(name=g0, qubits_t=[0], step=0)
                qc.add_gate(name=g1, qubits_t=[1], step=0)

                qc.assemble_symbolic_unitary(replace_symbolic_zeros_and_ones=False)
                assert qc.unitary.shape == (4, 4)
                #print(f"Unitary for {g0} and {g1}:")
                if g0 != g1:
                    tmp_mat = sp.physics.quantum.TensorProduct(qc.gate_collection.collections[g1][0].matrix22_t[1][0],
                                                                       qc.gate_collection.collections[g0][0].matrix22_t[0][0])
                else:
                    tmp_mat = sp.physics.quantum.TensorProduct(qc.gate_collection.collections[g1][1].matrix22_t[1][0],
                                                                       qc.gate_collection.collections[g0][0].matrix22_t[0][0])
                
                assert qc.unitary == tmp_mat

class TestSubs_symbolic_zeros_in_symbolic_unitary:
    def test_subs_symbolic_zeros_in_symbolic_unitary_2qubits(self):
        gates_to_be_tested = ['I', 'X', 'Y', 'Z', 'H']
        ids_to_be_zero = {'II': [[0,1],[0,2],[0,3],[1,0],[1,2],[1,3],[2,0],[2,1],[2,3],[3,0],[3,1],[3,2]],#II
                          'IX': [[0,0],[0,2],[0,3],[1,1],[1,2],[1,3],[2,0],[2,1],[2,2],[3,0],[3,1],[3,3]],#IX
                          'IY': [[0,0],[1,1],[0,2],[0,3],[1,1],[1,2],[1,3],[2,0],[2,1],[2,2],[3,0],[3,1],[3,3]],#IY
                          'IZ': [[0,1],[0,2],[0,3],[1,0],[1,2],[1,3],[2,0],[2,1],[2,3],[3,0],[3,1],[3,2]],#IZ
                          'IH': [[0,2],[0,3],[1,2],[1,3],[2,0],[2,1],[3,0],[3,1]],#IH
                          'XI': [[0,0],[0,1],[0,3],[1,0],[1,1],[1,2],[2,1],[2,2],[2,3],[3,0],[3,2],[3,3]],#XI
                          'XX': [[0,0],[0,1],[0,2],[1,0],[1,1],[1,3],[2,0],[2,2],[2,3],[3,1],[3,2],[3,3]],#XX
                          'XY': [[0,0],[0,1],[0,2],[1,0],[1,1],[1,3],[2,0],[2,2],[2,3],[3,1],[3,2],[3,3]],#XY
                          'XZ': [[0,0],[0,1],[0,3],[1,0],[1,1],[1,2],[2,1],[2,2],[2,3],[3,0],[3,2],[3,3]],#XZ
                          'XH': [[0,0],[0,1],[1,0],[1,1],[2,2],[2,3],[3,2],[3,3]],#XH
                          'YI': [[0,0],[0,1],[0,3],[1,0],[1,1],[1,2],[2,1],[2,2],[2,3],[3,0],[3,2],[3,3]],#YI
                          'YX': [[0,0],[0,1],[0,2],[1,0],[1,1],[1,3],[2,0],[2,2],[2,3],[3,1],[3,2],[3,3]],#YX
                          'YY': [[0,0],[0,1],[0,2],[1,0],[1,1],[1,3],[2,0],[2,2],[2,3],[3,1],[3,2],[3,3]],#YY
                          'YZ': [[0,0],[0,1],[0,3],[1,0],[1,1],[1,2],[2,1],[2,2],[2,3],[3,0],[3,2],[3,3]],#YZ
                          'YH': [[0,0],[0,1],[1,0],[1,1],[2,2],[2,3],[3,2],[3,3]],#YH
                          'ZI': [[0,1],[0,2],[0,3],[1,0],[1,2],[1,3],[2,0],[2,1],[2,3],[3,0],[3,1],[3,2]],#ZI
                          'ZX': [[0,0],[1,1],[0,2],[0,3],[1,1],[1,2],[1,3],[2,0],[2,1],[2,2],[3,0],[3,1],[3,3]],#ZX
                          'ZY': [[0,0],[1,1],[0,2],[0,3],[1,1],[1,2],[1,3],[2,0],[2,1],[2,2],[3,0],[3,1],[3,3]],#ZY
                          'ZZ': [[0,1],[0,2],[0,3],[1,0],[1,2],[1,3],[2,0],[2,1],[2,3],[3,0],[3,1],[3,2]],#ZZ
                          'ZH': [[0,2],[0,3],[1,2],[1,3],[2,0],[2,1],[3,0],[3,1]],#ZH
                          'HI': [[0,1],[0,3],[1,0],[1,2],[2,1],[2,3],[3,0],[3,2]],#HI
                          'HX': [[0,0],[0,2],[1,1],[1,3],[2,0],[2,2],[3,1],[3,3]],#HX
                          'HY': [[0,0],[0,2],[1,1],[1,3],[2,0],[2,2],[3,1],[3,3]],#HY
                          'HZ': [[0,1],[0,3],[1,0],[1,2],[2,1],[2,3],[3,0],[3,2]],#HZ
                          'HH': []}#HH
        for i0, g0 in enumerate(gates_to_be_tested):
            for i1, g1 in enumerate(gates_to_be_tested):
                qc = qp.QuantumCircuit(num_qubits=2)
                qc.add_gate(name=g0, qubits_t=[0], step=0)
                qc.add_gate(name=g1, qubits_t=[1], step=0)

                #print(f"Unitary for {g0} and {g1}:")
                qc.assemble_symbolic_unitary()
                qc.subs_symbolic_zerosones_in_symbolic_unitary()

                #print(qc.unitary)
                
                for id in ids_to_be_zero[g1+g0]:
                    if len(id) > 0:
                        #print(id)
                        assert 0 == qc.unitary[id[0],id[1]]
