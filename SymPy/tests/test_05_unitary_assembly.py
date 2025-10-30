import pathlib
import sys

sys.path.append(str(pathlib.Path(__file__).parent.parent))

import pytest
import numpy as np
import sympy as sp
import QSymPy as qp



def test_cx_gate_initialization():
    qc = qp.QuantumCircuit(num_qubits=2)
    qc.add_gate('CX', qubits_t=[1], qubits_c=[0], step=0)
    qc.assemble_symbolic_unitary(use_alternative_repr=False, replace_symbolic_zeros_and_ones=False)
    qc.create_numeric_unitary_from_symbolic(use_alternative_repr=False)
    a = sp.Matrix([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]])
    assert qc.unitary_numeric == a

def test_cx_gate_again():
    qc = qp.QuantumCircuit(num_qubits=2)
    qc.add_gate('CX', qubits_t=[1], qubits_c=[0], step=0)
    qc.add_gate('CX', qubits_t=[1], qubits_c=[0], step=1)
    qc.assemble_symbolic_unitary(use_alternative_repr=False, replace_symbolic_zeros_and_ones=False)
    qc.create_numeric_unitary_from_symbolic(use_alternative_repr=False)
    a = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    assert qc.unitary_numeric == a
