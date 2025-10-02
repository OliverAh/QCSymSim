import numpy as np
import itertools
import sympy as sp
import sympy.physics.quantum as sp_phy_qant

#from .. import gate_bases
from ..gate_bases import QuantumGate, QuantumGateParameterized, QuantumGateMultiQubit, QuantumGateMultiQubit


'''There are lots of gphase(...) used in stdgates.inc. But as far as I checked, they should be neglected.
No idea why they are there.
I guess best would be to verify in the tests (pytest).'''




class Barrier():

    def __init__(self, step: int=0):
        self.step = step
    
    def __str__(self):
        return f'Barrier at step {self.step}'
    def __repr__(self):
        return self.__str__()

class Identity_Gate(QuantumGate):
    '''Class of the identity gate. It is a subclass of QuantumGate.'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'I'
    def __init__(self,name: str='I', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name=name, name_short='I', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 0], [0, 1]])

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class Pauli_X_Gate(QuantumGate):
    '''Class of the Pauli X gate. It is a subclass of QuantumGate.'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'X'
    def __init__(self, name: str='X', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name=name, name_short='X', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[0, 1], [1, 0]])

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class Pauli_Y_Gate(QuantumGate):
    '''Class of the Pauli Y gate. It is a subclass of QuantumGate.'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'Y'
    def __init__(self, name: str='Y', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name=name, name_short='Y', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[0, -1j], [1j, 0]])

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class Pauli_Z_Gate(QuantumGate):
    '''Class of the Pauli Z gate. It is a subclass of QuantumGate.'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'Z'
    def __init__(self, name: str='Z', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name=name, name_short='Z', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 0], [0, -1]])

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class Hadamard_Gate(QuantumGate):
    '''Class of the Hadamard gate. It is a subclass of QuantumGate.'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'H'
    def __init__(self, name: str='H', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name=name, name_short='H', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 1], [1, -1]])/np.sqrt(2)

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class Hadamard_error_I_Gate(QuantumGate):
    '''Class of a Hadamard gate, which fails towards identity with varyaing parameter p. 
    p=0 is 0 error, p=1 means instead of Hadamard Identity is applied. Numeric gates are equal to Hadamard gate as intented way to use is to replace symbolic values.
    It is a subclass of QuantumGate.'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'H_eI'
    def __init__(self, name: str='H_eI', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name=name, name_short='H_eI', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 1], [1, -1]])/np.sqrt(2)

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class S_Gate(QuantumGate):
    '''Class of the S gate. It is a subclass of QuantumGate.'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'S'
    def __init__(self, name: str='S', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name=name, name_short='S', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 0], [0, 1j]])

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class Sdg_Gate(QuantumGate):
    '''Class of the Sdg gate. It is a subclass of QuantumGate.'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'Sdg'
    def __init__(self, name: str='Sdg', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name=name, name_short='Sdg', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 0], [0, -1j]])

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class T_Gate(QuantumGate):
    '''Class of the T gate. It is a subclass of QuantumGate.'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'T'
    def __init__(self, name: str='T', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name=name, name_short='T', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 0], [0, np.exp(1j*np.pi/4)]])

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class Tdg_Gate(QuantumGate):
    '''Class of the Tdg gate. It is a subclass of QuantumGate.'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'Tdg'
    def __init__(self, name: str='Tdg', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name=name, name_short='Tdg', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 0], [0, np.exp(-1j*np.pi/4)]])

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class SX_Gate(QuantumGate):
    '''Class of the SX gate. It is a subclass of QuantumGate.'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'SX'
    def __init__(self, name: str='SX', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name=name, name_short='SX', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = 1/2*np.array([[1+1j, 1-1j], [1-1j, 1+1j]])

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class U_Gate(QuantumGateParameterized):
    '''Class of the OpenQasm3 built-in U gate. It is a subclass of QuantumGate. https://openqasm.com/language/gates.html#built-in-gates 
    parameters must be a dict: {'theta': val, 'phi': val, 'lambda': val}  
    NOTE: There is a difference in the definition between OpenQasm2 and OpenQasm3, OQ3 = e^i(ϕ+λ)/2 OQ2'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'U'
    def __init__(self, name: str='U', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0, parameters: dict=None):
        assert parameters is not None and all(key in parameters for key in ['theta', 'phi', 'lambda']), "U gate requires parameters: {'theta': val, 'phi': val, 'lambda': val}"
        super().__init__(name=name, name_short='U', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.parameters = {'theta': parameters['theta'], 'phi': parameters['phi'], 'lambda': parameters['lambda']} if parameters is not None else {'theta': None, 'phi': None, 'lambda': None}
        self.atomics_alt = {key: sp.symbols(key+'_' + self.atomics['00'].name[:-4]) for key in self.parameters.keys()} # {'theta': 'theta_U_qt0_qc0_s0_p00', ...}
        self.matrix_alt = sp.Matrix([[sp.cos(self.atomics_alt['theta']/2)                                     , -sp.exp(sp.I*self.atomics_alt['lambda'])*sp.sin(self.atomics_alt['theta']/2)],
                                     [sp.exp(sp.I*self.atomics_alt['phi'])*sp.sin(self.atomics_alt['theta']/2), sp.exp(sp.I*(self.atomics_alt['phi']+self.atomics_alt['lambda']))*sp.cos(self.atomics_alt['theta']/2)]])
        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_alt = [[self.matrix_alt]]
        self.matrix_numeric = np.array(self.matrix_alt.subs({self.atomics_alt[key]: val for key, val in self.parameters.items()})).astype(complex)
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric 

class GP_Gate(QuantumGateParameterized):
    '''Class of the OpenQasm3 built-in gphase gate (global phase). It is a subclass of QuantumGate. https://openqasm.com/language/gates.html#built-in-gates 
    parameters must be a dict: {'gamma': val}  
    NOTE: There is a difference in the definition between OpenQasm2 and OpenQasm3, OQ3 = e^iγ/2 OQ2'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'GP'
    def __init__(self, name: str='GP', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0, parameters: dict=None):
        assert parameters is not None and all(key in parameters for key in ['gamma']), "GP gate requires parameters: {'gamma': val}"
        super().__init__(name=name, name_short='GP', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.parameters = {'gamma': parameters['gamma']}
        self.atomics_alt = {key: sp.symbols(key+'_' + self.atomics['00'].name[:-4]) for key in self.parameters.keys()} # {'gamma': 'gamma_GP_qt0_qc0_s0_p00'}
        self.matrix_alt = sp.Matrix([[sp.exp(sp.I*self.atomics_alt['gamma']), 0],
                                     [0, sp.exp(sp.I*self.atomics_alt['gamma'])]])
        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_alt = [[self.matrix_alt]]

        self.matrix_numeric = np.array(self.matrix_alt.subs({self.atomics_alt[key]: val for key, val in self.parameters.items()})).astype(complex)
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class RX_Gate(QuantumGateParameterized):
    '''Class of the RX(theta) gate. It is a subclass of QuantumGateParameterized. 
    Parameter must be a dict: {'theta': val}, for compatibility with gates of more parameters.  
    '''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'RX'
    def __init__(self, name: str='RX', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0, parameters: dict=None):
        super().__init__(name=name, name_short='RX', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.parameters = {'theta': parameters['theta']} if parameters is not None else {'theta': None}
        self.atomics_alt = {key: sp.symbols(key+'_' + self.atomics['00'].name[:-4]) for key in self.parameters.keys()} # {'theta': 'theta_U_qt0_qc0_s0_p00', ...}
        self.matrix_alt = sp.Matrix([[ sp.cos(self.atomics_alt['theta']/2)    , -sp.I*sp.sin(self.atomics_alt['theta']/2)],
                                     [-sp.I*sp.sin(self.atomics_alt['theta']/2), sp.cos(self.atomics_alt['theta']/2)]])
        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_alt = [[self.matrix_alt]]
        self.matrix_numeric = np.array(self.matrix_alt.subs({self.atomics_alt[key]: val for key, val in self.parameters.items()})).astype(complex)
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric 

class RY_Gate(QuantumGateParameterized):
    '''Class of the RY(theta) gate. It is a subclass of QuantumGateParameterized. 
    Parameter must be a dict: {'theta': val}, for compatibility with gates of more parameters.  
    '''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'RY'
    def __init__(self, name: str='RY', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0, parameters: dict=None):
        super().__init__(name=name, name_short='RY', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.parameters = {'theta': parameters['theta']} if parameters is not None else {'theta': None}
        self.atomics_alt = {key: sp.symbols(key+'_' + self.atomics['00'].name[:-4]) for key in self.parameters.keys()} # {'theta': 'theta_U_qt0_qc0_s0_p00', ...}
        self.matrix_alt = sp.Matrix([[sp.cos(self.atomics_alt['theta']/2), -sp.sin(self.atomics_alt['theta']/2)],
                                     [sp.sin(self.atomics_alt['theta']/2),  sp.cos(self.atomics_alt['theta']/2)]])
        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_alt = [[self.matrix_alt]]
        self.matrix_numeric = np.array(self.matrix_alt.subs({self.atomics_alt[key]: val for key, val in self.parameters.items()})).astype(complex)
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class RZ_Gate(QuantumGateParameterized):
    '''Class of the RZ(lambda) gate. It is a subclass of QuantumGateParameterized. 
    Parameter must be a dict: {'lambda': val}, for compatibility with gates of more parameters.  
    '''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'RZ'
    def __init__(self, name: str='RZ', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0, parameters: dict=None):
        assert parameters is not None and all(key in parameters for key in ['lambda']), "RZ gate requires parameters: {'lambda': val}"
        super().__init__(name=name, name_short='RZ', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.parameters = {'lambda': parameters['lambda']} if parameters is not None else {'lambda': None}
        self.atomics_alt = {key: sp.symbols(key+'_' + self.atomics['00'].name[:-4]) for key in self.parameters.keys()} # {'lambda': 'lambda_U_qt0_qc0_s0_p00', ...}
        self.matrix_alt = sp.Matrix([[    sp.exp(-sp.I*self.atomics_alt['lambda']/2), 0],
                                     [0,  sp.exp( sp.I*self.atomics_alt['lambda']/2)]])
        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_alt = [[self.matrix_alt]]
        self.matrix_numeric = np.array(self.matrix_alt.subs({self.atomics_alt[key]: val for key, val in self.parameters.items()})).astype(complex)
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric


class P_Gate(QuantumGateParameterized):
    '''Class of the OpenQasm3 built-in P gate (phase). It is a subclass of QuantumGate. https://openqasm.com/language/gates.html#built-in-gates 
    parameters must be a dict: {'lambda': val}  '''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'P'
    def __init__(self, name: str='P', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0, parameters: dict=None):
        super().__init__(name=name, name_short='P', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.parameters = {'lambda': parameters['lambda']}
        self.atomics_alt = {key: sp.symbols(key+'_' + self.atomics['00'].name[:-4]) for key in self.parameters.keys()} # {'lambda': 'lambda_P_qt0_qc0_s0_p00'}
        self.matrix_alt = sp.Matrix([[1, 0],
                                     [0, sp.exp(sp.I*self.atomics_alt['lambda'])]])
        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_alt = [[self.matrix_alt]]

        self.matrix_numeric = np.array(self.matrix_alt.subs({self.atomics_alt[key]: val for key, val in self.parameters.items()})).astype(complex)
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class CX_Gate(QuantumGateMultiQubit):
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'CX'
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0):
        assert len(qubits_t) == 1, 'CX gate can only be applied to a single target qubit.'
        assert len(qubits_c) == 1, 'CX gate can only be applied to a single control qubit.'
        assert qubits_t[0] != qubits_c[0], 'Control and target qubit must be different.'
        super().__init__(name='CX', name_short='CX', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','X')}, control_values_c={qubits_c[0]:(0,1)})

class CCX_Gate(QuantumGateMultiQubit):
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'CCX'
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0):
        assert len(qubits_t) == 1, 'CCX gate can only be applied to a single target qubit.'
        assert len(qubits_c) == 2, 'CCX gate must be applied to two control qubits.'
        assert qubits_t[0] not in qubits_c, 'Control and target qubits must be different.'
        super().__init__(name='CCX', name_short='CCX', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','I','I','X')}, control_values_c={qubits_c[0]:(0,0,1,1), qubits_c[1]:(0,1,0,1)})

class CCXX_Gate(QuantumGateMultiQubit):
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'CCXX'
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0):
        assert len(qubits_t) == 2, 'CCXX gate can only be applied to two target qubits.'
        assert len(qubits_c) == 2, 'CCXX gate must be applied to two control qubits.'
        assert qubits_t[0] not in qubits_c, 'Control and target qubits must be different.'
        super().__init__(name='CCXX', name_short='CCXX', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','I','I','X'), qubits_t[1]: ('I','I','I','X')}, 
                         control_values_c={qubits_c[0]:(0,0,1,1), qubits_c[1]:(0,1,0,1)})

class U_for_CU_Gate(QuantumGateParameterized):
    '''Class of the U gate as used in the controlled-U gate AS DEFINED IN OPENQASM3, which includes a 4th parameter that is used as global phase.
      It is a subclass of QuantumGateParameterized.'''
    is_should_be_listed_in_gate_collection = False
    name_gate_collection = 'U_for_CU_Gate'
    def __init__(self, name: str='U', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0, parameters: dict=None):
        assert parameters is not None and all(key in parameters for key in ['theta', 'phi', 'lambda', 'gamma']), "U_for_CU_Gate requires parameters: {'theta': val, 'phi': val, 'lambda': val, 'gamma': val}"
        super().__init__(name=name, name_short='U',shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step, parameters=parameters)
        self.parameters = {'theta': parameters['theta'], 'phi': parameters['phi'], 'lambda': parameters['lambda'], 'gamma': parameters['gamma']} if parameters is not None else {'theta': None, 'phi': None, 'lambda': None, 'gamma': None}
        self.atomics_alt = {key: sp.symbols(key+'_' + self.atomics['00'].name[:-4]) for key in self.parameters.keys()} # {'theta': 'theta_U_qt0_qc0_s0_p00', ...}
        self.matrix_alt = sp.exp(sp.I*self.atomics_alt['gamma']) \
                        * sp.Matrix([[sp.cos(self.atomics_alt['theta']/2)                                     , -sp.exp(sp.I*self.atomics_alt['lambda'])*sp.sin(self.atomics_alt['theta']/2)],
                                     [sp.exp(sp.I*self.atomics_alt['phi'])*sp.sin(self.atomics_alt['theta']/2), sp.exp(sp.I*(self.atomics_alt['phi']+self.atomics_alt['lambda']))*sp.cos(self.atomics_alt['theta']/2)]])
        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_alt = [[self.matrix_alt]]
        self.matrix_numeric = np.array(self.matrix_alt.subs({self.atomics_alt[key]: val for key, val in self.parameters.items()})).astype(complex)
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric
    
class CU_Gate(QuantumGateMultiQubit):
    '''Class of the controlled-U gate AS DEFINED IN OPENQASM3, which includes a 4th parameter that is used as global phase.
      It is a subclass of QuantumGateMultiQubit.'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'CU'
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0, parameters: dict=None):
        assert len(qubits_t) == 1, 'CU gate can only be applied to a single target qubit.'
        assert len(qubits_c) == 1, 'CU gate can only be applied to a single control qubit.'
        assert qubits_t[0] != qubits_c[0], 'Control and target qubit must be different.'
        assert parameters is not None and all([key in parameters.keys() for key in ['theta', 'phi', 'lambda', 'gamma']]), "CU gate needs parameters dict with keys 'theta', 'phi', 'lambda', 'gamma'." 
        super().__init__(name='CU', name_short='CU', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','U_for_CU_Gate')}, control_values_c={qubits_c[0]:(0,1)},
                         parameters=[parameters])
        

class CUU_Gate(QuantumGateMultiQubit):
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'CUU'
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0, parameters: dict=None):
        assert len(qubits_t) == 2, 'CUU gate can only be applied to two target qubits.'
        assert len(qubits_c) == 1, 'CU gate can only be applied to a single control qubit.'
        assert qubits_c[0] not in qubits_t, 'Control and target qubit must be different.'
        assert type(parameters) is list and len(parameters)==2 and all([key in ['theta', 'phi', 'lambda'] for p in parameters for key in p.keys()]), "CUU gate needs list of parameters dicts with keys 'theta', 'phi', 'lambda'."
        super().__init__(name='CUU', name_short='CUU', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','U'), qubits_t[1]: ('I','U')}, control_values_c={qubits_c[0]:(0,1)},
                         parameters=parameters)
        print('self.atomics_alt in CUU:', self.atomics_alt)

class CY_Gate(QuantumGateMultiQubit):
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'CY'
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0):
        assert len(qubits_t) == 1, 'CY gate can only be applied to a single target qubit.'
        assert len(qubits_c) == 1, 'CY gate can only be applied to a single control qubit.'
        assert qubits_t[0] != qubits_c[0], 'Control and target qubit must be different.'
        super().__init__(name='CY', name_short='CY', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','Y')}, control_values_c={qubits_c[0]:(0,1)})

class CZ_Gate(QuantumGateMultiQubit):
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'CZ'
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0):
        assert len(qubits_t) == 1, 'CZ gate can only be applied to a single target qubit.'
        assert len(qubits_c) == 1, 'CZ gate can only be applied to a single control qubit.'
        assert qubits_t[0] != qubits_c[0], 'Control and target qubit must be different.'
        super().__init__(name='CZ', name_short='CZ', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','Z')}, control_values_c={qubits_c[0]:(0,1)})

class CP_Gate(QuantumGateMultiQubit):
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'CP'
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0, parameters: dict=None):
        assert len(qubits_t) == 1, 'CP gate can only be applied to a single target qubit.'
        assert len(qubits_c) == 1, 'CP gate can only be applied to a single control qubit.'
        assert qubits_t[0] != qubits_c[0], 'Control and target qubit must be different.'
        assert parameters is not None and all([key in parameters.keys() for key in ['lambda']]), "CP gate needs parameters dict with key 'lambda'."
        super().__init__(name='CP', name_short='CP', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','P')}, control_values_c={qubits_c[0]:(0,1)},
                        parameters=[parameters])

class CRX_Gate(QuantumGateMultiQubit):
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'CRX'
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0, parameters: dict=None):
        assert len(qubits_t) == 1, 'CRX gate can only be applied to a single target qubit.'
        assert len(qubits_c) == 1, 'CRX gate can only be applied to a single control qubit.'
        assert qubits_t[0] != qubits_c[0], 'Control and target qubit must be different.'
        assert parameters is not None and all([key in parameters.keys() for key in ['theta']]), "CRX gate needs parameters dict with key 'theta'." 
        super().__init__(name='CRX', name_short='CRX', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','RX')}, control_values_c={qubits_c[0]:(0,1)},
                         parameters=[parameters])

class CRY_Gate(QuantumGateMultiQubit):
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'CRY'
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0, parameters: dict=None):
        assert len(qubits_t) == 1, 'CRY gate can only be applied to a single target qubit.'
        assert len(qubits_c) == 1, 'CRY gate can only be applied to a single control qubit.'
        assert qubits_t[0] != qubits_c[0], 'Control and target qubit must be different.'
        assert parameters is not None and all([key in parameters.keys() for key in ['theta']]), "CRY gate needs parameters dict with key 'theta'." 
        super().__init__(name='CRY', name_short='CRY', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','RY')}, control_values_c={qubits_c[0]:(0,1)},
                         parameters=[parameters])

class CRZ_Gate(QuantumGateMultiQubit):
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'CRZ'
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0, parameters: dict=None):
        assert len(qubits_t) == 1, 'CRZ gate can only be applied to a single target qubit.'
        assert len(qubits_c) == 1, 'CRZ gate can only be applied to a single control qubit.'
        assert qubits_t[0] != qubits_c[0], 'Control and target qubit must be different.'
        assert parameters is not None and all([key in parameters.keys() for key in ['theta']]), "CRZ gate needs parameters dict with key 'theta'." 
        super().__init__(name='CRZ', name_short='CRZ', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','RZ')}, control_values_c={qubits_c[0]:(0,1)},
                         parameters=[{'lambda': parameters['theta']}]) # OPENQasm uses 'lambda' as parameter name for RZ gate, but theta for CRZ gate

class CH_Gate(QuantumGateMultiQubit):
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'CH'
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0):
        assert len(qubits_t) == 1, 'CH gate can only be applied to a single target qubit.'
        assert len(qubits_c) == 1, 'CH gate can only be applied to a single control qubit.'
        assert qubits_t[0] != qubits_c[0], 'Control and target qubit must be different.'
        super().__init__(name='CH', name_short='CH', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','H')}, control_values_c={qubits_c[0]:(0,1)})

class SWAP_Gate(QuantumGateMultiQubit):
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'SWAP'
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0):
        assert len(qubits_t) == 2, 'SWAP gate can only be applied to two target qubits.'
        assert qubits_c is None or len(qubits_c) == 0, 'SWAP gate can only be applied to target qubits.'
        assert qubits_t[0] != qubits_t[1], 'Target qubits must be different.'
        super().__init__(name='SWAP', name_short='SWAP', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qt: ('I','X','Y','Z') for qt in qubits_t}, control_values_c=None)
        self.matrix = sp.Rational(1,2) * self.matrix
        self.matrix_numeric = 1/2 * self.matrix_numeric
        self.matrix22_t = {key: [sp.Rational(1, sp.Pow(2, 1/2))*m for m in l] for key,l in self.matrix22_t.items()}
        self.matrix22_t_numeric = {key: [1/np.sqrt(2)*m for m in l] for key,l in self.matrix22_t_numeric.items()}

class CSWAP_Gate(QuantumGateMultiQubit):
    def __init__(self, *args, **kwargs):
        raise NotImplementedError('CSWAP gate is not implemented yet.')
