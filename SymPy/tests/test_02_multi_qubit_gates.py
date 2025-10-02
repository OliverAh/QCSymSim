import pathlib
import sys

sys.path.append(str(pathlib.Path(__file__).parent.parent))

import pytest
import numpy as np
import sympy as sp
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
    assert (gate.matrix22_t_numeric[1][0] == [[1, 0], [0, 1]]).all()
    assert (gate.matrix22_c_numeric[0][0] == [[1, 0], [0, 0]]).all()
    assert (gate.matrix22_t_numeric[1][1] == [[0, 1], [1, 0]]).all()
    assert (gate.matrix22_c_numeric[0][1] == [[0, 0], [0, 1]]).all()
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
    assert (gate.matrix22_t_numeric[0][0] == [[1, 0], [0, 1]]).all()
    assert (gate.matrix22_c_numeric[1][0] == [[1, 0], [0, 0]]).all()
    assert (gate.matrix22_t_numeric[0][1] == [[0, 1], [1, 0]]).all()
    assert (gate.matrix22_c_numeric[1][1] == [[0, 0], [0, 1]]).all()
    assert gate.matrix_numeric.tolist() == [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]

def test_cx_gate_invalid_target_and_control():
    with pytest.raises(AssertionError, match="CX gate can only be applied to a single target qubit."):
        qp.CX_Gate(qubits_t=[1, 2], qubits_c=[0], step=5)
    with pytest.raises(AssertionError, match="CX gate can only be applied to a single control qubit."):
        qp.CX_Gate(qubits_t=[1], qubits_c=[0, 2], step=5)
    with pytest.raises(AssertionError, match="Control and target qubit must be different."):
        qp.CX_Gate(qubits_t=[1], qubits_c=[1], step=5)

def test_cu_gate_initialization():
    parameters = {'theta': np.pi/16, 'phi': np.pi/32, 'lambda': np.pi/64}
    with pytest.raises(AssertionError, match="CU gate needs parameters dict with keys 'theta', 'phi', 'lambda', 'gamma'."):
        qp.CU_Gate(qubits_t=[1], qubits_c=[0], step=15, parameters=parameters)
    parameters = {'theta': np.pi/16, 'phi': np.pi/32, 'lambda': np.pi/64, 'gamma': np.pi/128}
    gate = qp.CU_Gate(qubits_t=[0], qubits_c=[1], step=15, parameters=parameters)
    assert not isinstance(gate, qp.QuantumGateParameterized)
    assert isinstance(gate, qp.QuantumGateMultiQubit)
    assert hasattr(gate, 'parameters')
    assert hasattr(gate, 'matrix_alt')
    assert hasattr(gate, 'matrix22_t_alt')
    assert gate.name == 'CU'
    assert gate.shape == (4, 4)
    assert gate.qubits_t == [0]
    assert gate.qubits_c == [1]
    assert gate.atomics['00'].name == 'CU_qt0_qc1_s15_p00'
    assert gate.atomics_alt['theta'][0].name == 'theta_CU_U_for_CU_Gate_qt0_qc1_s15'
    assert gate.parameters == [parameters]
    assert np.allclose(gate.matrix_numeric, np.array([[1,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,0]]) + 
                       np.exp(1j*parameters['gamma'])*np.array([[0,0,0,0],[0,0,0,0],
                                             [0,0,np.cos(parameters['theta']/2), -np.exp(1j*parameters['lambda'])*np.sin(parameters['theta']/2)],
                                             [0,0,np.exp(1j*parameters['phi'])*np.sin(parameters['theta']/2), np.exp(1j*(parameters['phi']+parameters['lambda']))*np.cos(parameters['theta']/2)]]))
    
def test_cuu_gate_initialization():
    parameters = [{'theta': np.pi/16, 'phi': np.pi/32, 'lambda': np.pi/64}]
    with pytest.raises(AssertionError, match="CUU gate needs list of parameters dicts with keys 'theta', 'phi', 'lambda'."):
        gate = qp.CUU_Gate(qubits_t=[1, 2], qubits_c=[0], step=25, parameters=parameters)
    parameters = [{'theta': np.pi/16, 'phi': np.pi/32, 'lambda': np.pi/64}]*2
    gate = qp.CUU_Gate(qubits_t=[1, 2], qubits_c=[0], step=25, parameters=parameters)
    assert not isinstance(gate, qp.QuantumGateParameterized)
    assert isinstance(gate, qp.QuantumGateMultiQubit)
    assert hasattr(gate, 'parameters')
    assert hasattr(gate, 'matrix_alt')
    assert hasattr(gate, 'matrix22_t_alt')
    assert gate.name == 'CUU'
    assert gate.shape == (8, 8)
    assert gate.qubits_t == [1, 2]
    assert gate.qubits_c == [0]
    assert gate.atomics['00'].name == 'CUU_qt12_qc0_s25_p00'
    assert gate.atomics_alt['theta'][1].name == 'theta_CUU_U_qt2_qc0_s25'
    assert gate.parameters == parameters

def test_cy_gate_initialization():
    gate = qp.CY_Gate(qubits_t=[1], qubits_c=[0], step=35)
    assert gate.name == 'CY'
    assert gate.shape == (4, 4)
    assert gate.qubits_t == [1]
    assert gate.qubits_c == [0]
    assert gate.matrix22_t[1][0] is not None
    assert gate.matrix22_c[0][0] is not None
    assert gate.matrix22_t[1][1] is not None
    assert gate.matrix22_c[0][1] is not None
    assert (gate.matrix22_t_numeric[1][0] == [[1, 0], [0, 1]]).all()
    assert (gate.matrix22_c_numeric[0][0] == [[1, 0], [0, 0]]).all()
    assert (gate.matrix22_t_numeric[1][1] == [[0, -1j], [1j, 0]]).all()
    assert (gate.matrix22_c_numeric[0][1] == [[0, 0], [0, 1]]).all()
    assert gate.matrix_numeric.tolist() == [[1, 0, 0, 0], [0, 0, 0, -1j], [0, 0, 1, 0], [0, 1j, 0, 0]]

def test_swap_gate_initialization():
    gate = qp.SWAP_Gate(qubits_t=[2, 3], step=45)
    assert gate.name == 'SWAP'
    assert gate.shape == (4, 4)
    assert gate.qubits_t == [2, 3]
    assert gate.qubits_c in ([], None)
    assert gate.matrix22_t[2][3] is not None
    assert len(gate.matrix22_t) == 2
    assert np.allclose(gate.matrix22_t_numeric[2][3], 1/np.sqrt(2)*np.array([[1, 0], [0, -1]]))
    assert np.allclose(gate.matrix22_t_numeric[3][2], 1/np.sqrt(2)*np.array([[0, -1j], [1j, 0]]))
    assert gate.matrix_numeric.tolist() == [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]

def test_cswap_gate_initialization():
    with pytest.raises(NotImplementedError, match="CSWAP gate is not implemented yet."):
        gate = qp.CSWAP_Gate(qubits_t=[2, 3], qubits_c=[1], step=55)
