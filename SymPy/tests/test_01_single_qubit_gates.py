import pathlib
import sys

sys.path.append(str(pathlib.Path(__file__).parent.parent))

import pytest
import numpy as np
import sympy as sp
import QSymPy as qp


def test_identity_gate_initialization():
    gate = qp.Identity_Gate(qubits_t=[0], step=0)
    assert isinstance(gate, qp.QuantumGate)
    assert not isinstance(gate, qp.QuantumGateParameterized)
    assert not isinstance(gate, qp.QuantumGateMultiQubit)
    assert gate.parameters is None
    assert not hasattr(gate, 'matrix_alt')
    assert not hasattr(gate, 'matrix22_t_alt')
    assert gate.name == 'I'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [0]
    assert gate.atomics['00'].name == 'I_qt0_qc_s0_p00'
    assert gate.matrix[0,0] == gate.atomics['00']
    assert gate.matrix_numeric.tolist() == [[1, 0], [0, 1]]
    assert gate.num_summands_decomposed == 1
    assert len(gate.matrix22_t) == 1
    assert len(gate.matrix22_t_numeric) == 1
    assert gate.matrix22_c is None
    assert gate.matrix22_c_numeric is None
    assert gate.is_parametric is False

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

def test_hadamard_error_I_gate_initialization():
    gate = qp.Hadamard_error_I_Gate(qubits_t=[0], step=0)
    assert gate.name == 'H_eI'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [0]
    assert gate.matrix_numeric.tolist() == [[1/2**0.5, 1/2**0.5], [1/2**0.5, -1/2**0.5]]

def test_s_gate_initialization():
    gate = qp.S_Gate(qubits_t=[0], step=5)
    assert gate.name == 'S'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [0]
    assert (gate.matrix_numeric == [[1, 0], [0, 1j]]).all()

def test_sdg_gate_initialization():
    gate = qp.Sdg_Gate(qubits_t=[0], step=6)
    assert gate.name == 'Sdg'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [0]
    assert (gate.matrix_numeric == [[1, 0], [0, -1j]]).all()

def test_t_gate_initialization():
    gate = qp.T_Gate(qubits_t=[0], step=7)
    assert gate.name == 'T'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [0]
    assert (gate.matrix_numeric == [[1, 0], [0, np.exp(1j*np.pi/4)]]).all()

def test_tdg_gate_initialization():
    gate = qp.Tdg_Gate(qubits_t=[0], step=8)
    assert gate.name == 'Tdg'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [0]
    assert (gate.matrix_numeric == [[1, 0], [0, np.exp(-1j*np.pi/4)]]).all()

def test_u_gate_initialization():
    parameters = {'theta': np.pi/16, 'phi': np.pi/32, 'lambda': np.pi/64, 'gamma': np.pi/128}
    with pytest.raises(AssertionError, match="U gate requires parameters: {'theta': val, 'phi': val, 'lambda': val}"):
        qp.U_Gate(qubits_t=[0], step=9, parameters={'theta': np.pi/2, 'phi': np.pi/4})
    parameters = {'theta': np.pi/2, 'phi': np.pi/4, 'lambda': np.pi/8}
    gate = qp.U_Gate(qubits_t=[0], step=9, parameters=parameters)
    assert isinstance(gate, qp.QuantumGateParameterized)
    assert not isinstance(gate, qp.QuantumGateMultiQubit)
    assert hasattr(gate, 'parameters')
    assert hasattr(gate, 'matrix_alt')
    assert hasattr(gate, 'matrix22_t_alt')
    assert gate.name == 'U'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [0]
    assert gate.atomics['00'].name == 'U_qt0_qc_s9_p00'
    assert gate.atomics_alt['theta'].name == 'theta_U_qt0_qc_s9'
    assert gate.parameters == parameters
    assert np.allclose(gate.matrix_numeric, np.array([[np.cos(parameters['theta']/2), -np.exp(1j*parameters['lambda'])*np.sin(parameters['theta']/2)],
                                   [np.exp(1j*parameters['phi'])*np.sin(parameters['theta']/2), np.exp(1j*(parameters['phi']+parameters['lambda']))*np.cos(parameters['theta']/2)]]))

def test_gp_gate_initialization():
    parameters = {'phi': np.pi/16}
    with pytest.raises(AssertionError, match="GP gate requires parameters: {'gamma': val}"):
        qp.GP_Gate(qubits_t=[0], step=10, parameters={})
    parameters = {'gamma': np.pi/16 }
    gate = qp.GP_Gate(qubits_t=[0], step=10, parameters=parameters)
    assert isinstance(gate, qp.QuantumGateParameterized)
    assert not isinstance(gate, qp.QuantumGateMultiQubit)
    assert hasattr(gate, 'parameters')
    assert hasattr(gate, 'matrix_alt')
    assert hasattr(gate, 'matrix22_t_alt')
    assert gate.name == 'GP'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [0]
    assert gate.atomics['00'].name == 'GP_qt0_qc_s10_p00'
    assert gate.atomics_alt['gamma'].name == 'gamma_GP_qt0_qc_s10'
    assert gate.parameters == parameters
    assert np.allclose(gate.matrix_numeric, np.exp(1j*parameters['gamma']) * np.array([[1, 0], [0, 1]]))

def test_rx_gate_initialization():
    gate = qp.RX_Gate(qubits_t=[0], step=11, parameters={'theta': np.pi/16})
    assert gate.name == 'RX'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [0]
    assert (gate.matrix_numeric == [[np.cos(np.pi/32), -1j*np.sin(np.pi/32)], [-1j*np.sin(np.pi/32), np.cos(np.pi/32)]]).all()

def test_ry_gate_initialization():
    gate = qp.RY_Gate(qubits_t=[0], step=12, parameters={'theta': np.pi/16})
    assert gate.name == 'RY'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [0]
    assert (gate.matrix_numeric == [[np.cos(np.pi/32), -np.sin(np.pi/32)], [np.sin(np.pi/32), np.cos(np.pi/32)]]).all()

def test_rz_gate_initialization():
    gate = qp.RZ_Gate(qubits_t=[0], step=13, parameters={'lambda': np.pi/16})
    assert gate.name == 'RZ'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [0]
    assert (gate.matrix_numeric == [[np.exp(-1j*np.pi/32), 0], [0, np.exp(1j*np.pi/32)]]).all()

def test_p_gate_initialization():
    gate = qp.P_Gate(qubits_t=[0], step=14, parameters={'lambda': np.pi/16})
    assert gate.name == 'P'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [0]
    assert (gate.matrix_numeric == [[1, 0], [0, np.exp(1j*np.pi/16)]]).all()

def test_uforcu_gate_initialization():
    parameters = {'theta': np.pi/16, 'phi': np.pi/32, 'lambda': np.pi/64}
    with pytest.raises(AssertionError, match="U_for_CU_Gate requires parameters: {'theta': val, 'phi': val, 'lambda': val, 'gamma': val}"):
        qp.U_for_CU_Gate(qubits_t=[1], qubits_c=[0], step=15, parameters=parameters)
    parameters = {'theta': np.pi/16, 'phi': np.pi/32, 'lambda': np.pi/64, 'gamma': np.pi/128}
    gate = qp.U_for_CU_Gate(qubits_t=[1], qubits_c=[0], step=15, parameters=parameters)
    assert isinstance(gate, qp.QuantumGateParameterized)
    assert not isinstance(gate, qp.QuantumGateMultiQubit)
    assert hasattr(gate, 'parameters')
    assert hasattr(gate, 'matrix_alt')
    assert hasattr(gate, 'matrix22_t_alt')
    assert gate.name == 'U'
    assert gate.shape == (2, 2)
    assert gate.qubits_t == [1]
    assert gate.qubits_c == [0]
    assert gate.atomics['00'].name == 'U_qt1_qc0_s15_p00'
    assert gate.atomics_alt['theta'].name == 'theta_U_qt1_qc0_s15'
    assert gate.parameters == parameters
    assert np.allclose(gate.matrix_numeric, np.exp(1j*parameters['gamma']) \
                       * np.array([[np.cos(parameters['theta']/2), -np.exp(1j*parameters['lambda'])*np.sin(parameters['theta']/2)],
                                   [np.exp(1j*parameters['phi'])*np.sin(parameters['theta']/2), np.exp(1j*(parameters['phi']+parameters['lambda']))*np.cos(parameters['theta']/2)]]))