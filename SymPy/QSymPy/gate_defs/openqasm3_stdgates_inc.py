import numpy as np
import itertools
import sympy as sp
import sympy.physics.quantum as sp_phy_qant

#from .. import gate_bases
from ..gate_bases import QuantumGate, QuantumGateParameterized, QuantumGateMultiQubit, QuantumGateMultiQubit_new





class Identity_Gate(QuantumGate):
    '''Class of the identity gate. It is a subclass of QuantumGate.'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'I'
    def __init__(self,name: str='I', qubits_t: list[int]=[0], qubits_c: list[int]=[0], step: int=0):
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
    def __init__(self, qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name='Y', name_short='Y', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[0, -1j], [1j, 0]])

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class Pauli_Z_Gate(QuantumGate):
    '''Class of the Pauli Z gate. It is a subclass of QuantumGate.'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'Z'
    def __init__(self, qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name='Z', name_short='Z', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 0], [0, -1]])

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class Hadamard_Gate(QuantumGate):
    '''Class of the Hadamard gate. It is a subclass of QuantumGate.'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'H'
    def __init__(self, qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name='H', name_short='H', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 1], [1, -1]])/np.sqrt(2)

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class Hadamard_error_I_Gate(QuantumGate):
    '''Class of a Hadamard gate, which fails towards identity with varyaing parameter p. 
    p=0 is 0 error, p=1 means instead of Hadamard Identity is applied. Numeric gates are equal to Hadamard gate as intented way to use is to replace symbolic values.
    It is a subclass of QuantumGate.'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'H_eI'
    def __init__(self, qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name='H_eI', name_short='H_eI', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 1], [1, -1]])/np.sqrt(2)

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class U_Gate(QuantumGateParameterized):
    '''Class of the OpenQasm3 built-in U gate. It is a subclass of QuantumGate. https://openqasm.com/language/gates.html#built-in-gates 
    parameters must be a dict: {'theta': val, 'phi': val, 'lambda': val}  
    NOTE: There is a difference in the definition between OpenQasm2 and OpenQasm3, OQ3 = e^i(ϕ+λ)/2 OQ2'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'U'
    def __init__(self, qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0, parameters: dict=None):
        super().__init__(name='U', name_short='U', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.parameters = {'theta': parameters['theta'], 'phi': parameters['phi'], 'lambda': parameters['lambda']}
        self.atomics_alt = {key: sp.symbols(key+'_' + self.atomics['00'].name[:-4]) for key in self.parameters.keys()} # {'theta': 'theta_U_qt0_qc0_s0_p00', ...}
        self.matrix_alt = sp.Rational(1,2) * sp.Matrix([[1+sp.exp(sp.I*self.atomics_alt['theta'])                                        , -sp.I*sp.exp(sp.I*self.atomics_alt['lambda'])*(1-sp.exp(sp.I*self.atomics_alt['theta']))],
                                                        [sp.I*sp.exp(sp.I*self.atomics_alt['phi'])*(1-sp.exp(sp.I*self.atomics_alt['theta'])), sp.exp(sp.I*(self.atomics_alt['phi']+self.atomics_alt['lambda']))*(1+sp.exp(sp.I*self.atomics_alt['theta']))]])
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
    def __init__(self, qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0, parameters: dict=None):
        super().__init__(name='GP', name_short='GP', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.parameters = {'gamma': parameters['gamma']}
        self.atomics_alt = {key: sp.symbols(key+'_' + self.atomics['00'].name[:-4]) for key in self.parameters.keys()} # {'gamma': 'gamma_GP_qt0_qc0_s0_p00'}
        self.matrix_alt = sp.Matrix([[sp.exp(sp.I*self.atomics_alt['gamma']), 0],
                                     [0, sp.exp(sp.I*self.atomics_alt['gamma'])]])
        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_alt = [[self.matrix_alt]]

        self.matrix_numeric = np.array(self.matrix_alt.subs({self.atomics_alt[key]: val for key, val in self.parameters.items()})).astype(complex)
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class P_Gate(QuantumGateParameterized):
    '''Class of the OpenQasm3 built-in P gate (phase). It is a subclass of QuantumGate. https://openqasm.com/language/gates.html#built-in-gates 
    parameters must be a dict: {'lambda': val}  '''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'P'
    def __init__(self, qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0, parameters: dict=None):
        super().__init__(name='P', name_short='P', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.parameters = {'lambda': parameters['lambda']}
        self.atomics_alt = {key: sp.symbols(key+'_' + self.atomics['00'].name[:-4]) for key in self.parameters.keys()} # {'lambda': 'lambda_P_qt0_qc0_s0_p00'}
        self.matrix_alt = sp.Matrix([[1, 0],
                                     [0, sp.exp(sp.I*self.atomics_alt['lambda'])]])
        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_alt = [[self.matrix_alt]]

        self.matrix_numeric = np.array(self.matrix_alt.subs({self.atomics_alt[key]: val for key, val in self.parameters.items()})).astype(complex)
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class Barrier():

    def __init__(self, step: int=0):
        self.step = step
    
    def __str__(self):
        return f'Barrier at step {self.step}'
    def __repr__(self):
        return self.__str__()


class CX_Gate(QuantumGateMultiQubit):
    '''Class of the CX gate. It is a subclass of QuantumGateMultiQubit.
    For CX_c_t the matrix is |0_c⟩⟨0_c|⊗I_t+|1_c⟩⟨1_c|⊗X_t. For CX_t_c the matrix is I_t⊗|0_c⟩⟨0_c|+X_t⊗|1_c⟩⟨1_c|.
    For multiple qubits surrounding and/or between c and t qubit, Identity gates are inserted, e.g. for CX_3_c_1_t
    the matrix is I_3⊗|0_c⟩⟨0_c|⊗I_1⊗I_t+I_3⊗|1_c⟩⟨1_c|⊗I_1⊗X_1.'''
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'CX'
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0):
        super().__init__(name='CX', name_short='CX', shape=(4,4), qubits_t=qubits_t, qubits_c=qubits_c, step=step, num_summands_decomposed=2)
        assert len(qubits_t) == 1, 'CX gate can only be applied to a single target qubit.'
        assert len(qubits_c) == 1, 'CX gate can only be applied to a single control qubit.'
        assert qubits_t[0] != qubits_c[0], 'Control and target qubit must be different.'

        #self.num_summands_decomposed = 2

        iter_indices_22 = [str(sub[0])+str(sub[1]) for sub in itertools.product(range(2), range(2))] # (2,2) -> ['00', '01', '10', '11']
        iter_name_ext_orig = ['0c', '1c', 'It', 'Xt']
        iter_name_ext = [i[0]+i[1] for i in itertools.product(iter_name_ext_orig, iter_indices_22)] # ['0c00', '0c01', ... , 'Xt10', 'Xt11']
        iter_indices = [str(sub[0])+str(sub[1]) for sub in itertools.product(range(self.shape[0]), range(self.shape[1]))] # (4,4) -> ['00', '01', '02', '03', '10', '11', '12', '13', ...]
        #print(iter_name_ext)
        str_qubits_t = ''.join([str(q)+'_' for q in qubits_t])
        str_qubits_c = ''.join([str(q)+'_' for q in qubits_c]) if qubits_c is not None else ''
        self.atomics = {sub: sp.symbols(self.name+nex+'_qt'+str_qubits_t+'qc'+str_qubits_c+'s'+str(step)+'_p'+sub) for nex, sub in zip(iter_name_ext, iter_indices)}

        del iter_indices_22, iter_name_ext_orig, iter_name_ext, iter_indices, str_qubits_t, str_qubits_c


        statevec_zero = np.array([[1],[0]], dtype=np.int8)
        statevec_one  = np.array([[0],[1]], dtype=np.int8)
        densmat_zero  = np.kron(statevec_zero, statevec_zero.T)
        densmat_one   = np.kron(statevec_one, statevec_one.T)
        sym_statevec_zero = sp.Matrix([[1],[0]])
        sym_statevec_one  = sp.Matrix([[0],[1]])
        sym_densmat_zero  = sp_phy_qant.TensorProduct(sym_statevec_zero, sym_statevec_zero.transpose())
        sym_densmat_one   = sp_phy_qant.TensorProduct(sym_statevec_one, sym_statevec_one.transpose())
        #sym_eye_t = _Identity_Gate_in_MQG(name='CX_I', qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        #sym_x_t   =  _Pauli_X_Gate_in_MQG(name='CX_X', qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        sym_eye_t = Identity_Gate(name='CX_I', qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        sym_x_t   =  Pauli_X_Gate(name='CX_X', qubits_t=qubits_t, qubits_c=qubits_c, step=step)

        self.matrix22_t[qubits_t[0]][0] = sym_x_t.matrix
        self.matrix22_c[qubits_c[0]][0] = sym_densmat_one
        self.matrix22_t[qubits_t[0]][1] = sym_eye_t.matrix
        self.matrix22_c[qubits_c[0]][1] = sym_densmat_zero

        self.matrix22_t_numeric[qubits_t[0]][0] = np.array([[0, 1], [1, 0]],dtype=np.int8)
        self.matrix22_c_numeric[qubits_c[0]][0] = densmat_one
        self.matrix22_t_numeric[qubits_t[0]][1] = np.array([[1, 0], [0, 1]],dtype=np.int8)
        self.matrix22_c_numeric[qubits_c[0]][1] = densmat_zero

        if qubits_c[0] > qubits_t[0]:
            m1 = sp_phy_qant.TensorProduct(sym_densmat_one, sym_x_t.matrix)
            m2 = sp_phy_qant.TensorProduct(sym_densmat_zero, sym_eye_t.matrix)
            #display(m1)
            #display(m2)
            self.matrix = m1 + m2
            self.matrix_numeric = np.kron(densmat_one, np.array([[0, 1], [1, 0]],dtype=np.int8)) + np.kron(densmat_zero, np.eye(2,dtype=np.int8))
            # = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
        elif qubits_c[0] < qubits_t[0]:
            m1 = sp_phy_qant.TensorProduct(sym_x_t.matrix, sym_densmat_one)
            m2 = sp_phy_qant.TensorProduct(sym_eye_t.matrix, sym_densmat_zero)
            #display(m1)
            #display(m2)
            self.matrix = m1 + m2
            self.matrix_numeric = np.kron(np.array([[0, 1], [1, 0]],dtype=np.int8), densmat_one) + np.kron(np.eye(2,dtype=np.int8), densmat_zero)
            # = np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0])
        del statevec_zero, statevec_one, densmat_zero, densmat_one, sym_statevec_zero, sym_statevec_one, sym_densmat_zero, sym_densmat_one, sym_eye_t, sym_x_t

class CX_Gate_new(QuantumGateMultiQubit_new):
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'CX_new'
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0):
        assert len(qubits_t) == 1, 'CX gate can only be applied to a single target qubit.'
        assert len(qubits_c) == 1, 'CX gate can only be applied to a single control qubit.'
        assert qubits_t[0] != qubits_c[0], 'Control and target qubit must be different.'
        super().__init__(name='CX_new', name_short='CX', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','X')}, control_values_c={qubits_c[0]:(0,1)})

class CCX_Gate_new(QuantumGateMultiQubit_new):
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'CCX_new'
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0):
        assert len(qubits_t) == 1, 'CCX gate can only be applied to a single target qubit.'
        assert len(qubits_c) == 2, 'CCX gate must be applied to two control qubits.'
        assert qubits_t[0] not in qubits_c, 'Control and target qubits must be different.'
        super().__init__(name='CCX_new', name_short='CCX', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','I','I','X')}, control_values_c={qubits_c[0]:(0,0,1,1), qubits_c[1]:(0,1,0,1)})

class CCXX_Gate_new(QuantumGateMultiQubit_new):
    is_should_be_listed_in_gate_collection = True
    name_gate_collection = 'CCXX_new'
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0):
        assert len(qubits_t) == 2, 'CCXX gate can only be applied to two target qubits.'
        assert len(qubits_c) == 2, 'CCXX gate must be applied to two control qubits.'
        assert qubits_t[0] not in qubits_c, 'Control and target qubits must be different.'
        super().__init__(name='CCXX_new', name_short='CCXX', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','I','I','X'), qubits_t[1]: ('I','I','I','X')}, 
                         control_values_c={qubits_c[0]:(0,0,1,1), qubits_c[1]:(0,1,0,1)})
