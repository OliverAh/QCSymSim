import pathlib
import sys

sys.path.append(str(pathlib.Path(__file__).parent.parent))

import pytest
import QSymPy as qp



def test_cx_gate_initialization():
    gate = qp.CX_Gate(qubits_t=[1], qubits_c=[0], step=5)
    assert gate.name == 'CX'
    assert gate.shape == (4, 4)
    assert gate.qubits_t == [1]
    assert gate.qubits_c == [0]
    assert gate.matrix22_t[1][0] is not None
    assert gate.matrix22_c[0][0] is not None
    assert gate.matrix22_t[1][1] is not None
    assert gate.matrix22_c[0][1] is not None
    assert (gate.matrix22_t_numeric[1][0] == [[0, 1], [1, 0]]).all()
    assert (gate.matrix22_c_numeric[0][0] == [[0, 0], [0, 1]]).all()
    assert (gate.matrix22_t_numeric[1][1] == [[1, 0], [0, 1]]).all()
    assert (gate.matrix22_c_numeric[0][1] == [[1, 0], [0, 0]]).all()
    assert gate.matrix_numeric.tolist() == [[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]]

def test_cx_gate_interchanged_qubits():
    gate = qp.CX_Gate(qubits_t=[0], qubits_c=[1], step=5)
    assert gate.name == 'CX'
    assert gate.shape == (4, 4)
    assert gate.qubits_t == [0]
    assert gate.qubits_c == [1]
    assert gate.matrix22_t[0][1] is not None
    assert gate.matrix22_c[1][0] is not None
    assert gate.matrix22_t[0][0] is not None
    assert gate.matrix22_c[1][1] is not None
    assert (gate.matrix22_t_numeric[0][0] == [[0, 1], [1, 0]]).all()
    assert (gate.matrix22_c_numeric[1][0] == [[0, 0], [0, 1]]).all()
    assert (gate.matrix22_t_numeric[0][1] == [[1, 0], [0, 1]]).all()
    assert (gate.matrix22_c_numeric[1][1] == [[1, 0], [0, 0]]).all()
    assert gate.matrix_numeric.tolist() == [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]

def test_cx_gate_invalid_target_and_control():
    with pytest.raises(AssertionError, match="CX gate can only be applied to a single target qubit."):
        qp.CX_Gate(qubits_t=[1, 2], qubits_c=[0], step=5)
    with pytest.raises(AssertionError, match="CX gate can only be applied to a single control qubit."):
        qp.CX_Gate(qubits_t=[1], qubits_c=[0, 2], step=5)
    with pytest.raises(AssertionError, match="Control and target qubit must be different."):
        qp.CX_Gate(qubits_t=[1], qubits_c=[1], step=5)