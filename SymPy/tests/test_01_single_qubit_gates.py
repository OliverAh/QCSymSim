import pathlib
import sys

sys.path.append(str(pathlib.Path(__file__).parent.parent))

import pytest
import QSymPy as qp


def test_identity_gate_initialization():
    gate = qp.Identity_Gate(qubits_t=[0], qubits_c=[1], step=0)
    assert gate.name == 'I'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [0]
    assert gate.qubits_c == [1]
    assert gate.matrix_numeric.tolist() == [[1, 0], [0, 1]]

def test_pauli_x_gate_initialization():
    gate = qp.Pauli_X_Gate(qubits_t=[0], step=1)
    assert gate.name == 'X'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [0]
    assert gate.matrix_numeric.tolist() == [[0, 1], [1, 0]]

def test_pauli_y_gate_initialization():
    gate = qp.Pauli_Y_Gate(qubits_t=[0], step=2)
    assert gate.name == 'Y'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [0]
    assert (gate.matrix_numeric == [[0, -1j], [1j, 0]]).all()

def test_pauli_z_gate_initialization():
    gate = qp.Pauli_Z_Gate(qubits_t=[0], step=3)
    assert gate.name == 'Z'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [0]
    assert gate.matrix_numeric.tolist() == [[1, 0], [0, -1]]

def test_hadamard_gate_initialization():
    gate = qp.Hadamard_Gate(qubits_t=[0], step=4)
    assert gate.name == 'H'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [0]
    assert (gate.matrix_numeric == [[1/2**0.5, 1/2**0.5], [1/2**0.5, -1/2**0.5]]).all()

