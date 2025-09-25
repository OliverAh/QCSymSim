import pathlib
import sympy as sp
import itertools
import numpy as np
#import sympy.physics as sp_phy
import sympy.physics.quantum as sp_phy_qant
import openqasm3
import copy
from typing import Literal, Mapping, Iterable
import functools


class QuantumGate:
    '''Base class for quantum gates. It holds information of the operation that is applied to the qubit(s).
    It does not hold the matrix for the whole quantum circuit, but just for the qubit(s) it's supposed to change the state of.
    Lists of qubits and matrices are assumed to be encoded as most significant qubit first, i.e. left and bottom, e.g. |q2 q1 q0⟩ = |q2> kron |q1> kron |q0>, 
    |0 0 1> = (1,0) kron (1,0) kron (0,1) = (0,1,0,0,0,0,0,0).'''
    # We need information about the qubit and step here, because otherwise the atomic symbols could not be distinguished from those of other gates
    
    def __init__(self, name: str='', shape:tuple[int,int]=[2,2], qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0, num_summands_decomposed: int=1):
        self.name = name
        self.shape = shape
        iter_indices = [str(sub[0])+str(sub[1]) for sub in itertools.product(range(shape[0]), range(shape[1]))] # (2,2) -> ['00', '01', '10', '11']
        str_qubits_t = ''.join([str(q)+'' for q in qubits_t])
        str_qubits_c = ''.join([str(q)+'' for q in qubits_c]) if qubits_c is not None else ''
        self.atomics = {sub: sp.symbols(self.name+'_qt'+str_qubits_t+'_qc'+str_qubits_c+'_s'+str(step)+'_p'+str(sub[0])+str(sub[1])) for sub in iter_indices} # {'00': I_qt0_qc0_s0_p00, '01': I_qt0_qc0_s0_p01, ...}
        self.matrix = sp.Matrix([[self.atomics[col] for col in itertools.islice(iter_indices, row*shape[1], row*shape[1]+shape[1])] for row in range(shape[0])]) # [[I_..._p00, I_..._p01], [I_..._p10, I_..._p11]]
        self.ids_matrix_zeros = None
        self.matrix_numeric = None
        self.qubits = [qubits_t] if qubits_c is None else [qubits_t, qubits_c]
        self.qubits_t = qubits_t
        self.qubits_c = qubits_c
        self.step = step
        self.num_summands_decomposed = num_summands_decomposed
        self.matrix22_t = {q_t: [None]*self.num_summands_decomposed for q_t in qubits_t}
        self.matrix22_c = {q_c: [None]*self.num_summands_decomposed for q_c in qubits_c} if qubits_c is not None else None
        self.matrix22_t_numeric = {q_t: [None]*self.num_summands_decomposed for q_t in qubits_t}
        self.matrix22_c_numeric = {q_c: [None]*self.num_summands_decomposed for q_c in qubits_c} if qubits_c is not None else None
        
    #def __dict__(self):
    #    return {'name:': self.name, 'shape': self.shape, 'atomics': self.atomics, 'matrix': self.matrix, 'qubits': self.qubits, 'step':self.step}
    #def __repr__(self):
    #    return self.__dict__().__str__()
    def __iter__(self):
        for key in self.__dict__:
            yield key, getattr(self, key)
    def items(self):
        return self.__iter__()

class QuantumGateParameterized(QuantumGate):
    '''Base class for parameterized quantum gates. It holds information of the operation that is applied to the qubit(s).
      It is a subclass of QuantumGate, and uses an additional symbolic matrix.'''
    def __init__(self, name: str='', shape:tuple[int,int]=[2,2], qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0, num_summands_decomposed: int=1):
        super().__init__(name, shape, qubits_t, qubits_c, step, num_summands_decomposed)
        self.parameters = None # {name: value, ...}
        self.matrix_alt = self.matrix.copy()
        self.matrix22_t_alt = copy.deepcopy(self.matrix22_t) # deepcopy just to be sure that _alt and non-alt are not linked
        

class QuantumGateMultiQubit(QuantumGate):
    '''Base class for multi-qubit quantum gates. It holds information of the operation that is applied to the qubits. It is a subclass of QuantumGate,
      and uses a different way to initialize self.matrix'''


    def __init__(self, name: str='', shape:tuple[int,int]|None=[4,4], qubits_t: list[int]=[0], qubits_c: list[int]|None=None, step: int=0, num_summands_decomposed: int|None=2):
        super().__init__(name, shape, qubits_t, qubits_c, step, num_summands_decomposed)
        self.atomics = None
        self.matrix = None
        #self.matrix_numeric = None
        #self.num_summands_decomposed = self._num_summands_decomposed
        #self.matrix22_t = {q_t: [None]*self.num_summands_decomposed for q_t in qubits_t}
        #self.matrix22_c = {q_c: [None]*self.num_summands_decomposed for q_c in qubits_c}
        #self.matrix22_t_numeric = {q_t: [None]*self.num_summands_decomposed for q_t in qubits_t}
        #self.matrix22_c_numeric = {q_c: [None]*self.num_summands_decomposed for q_c in qubits_c}

class QuantumGateMultiQubit_new(QuantumGate):
    '''Base class for multi-qubit quantum gates. It holds information of the operation that is applied to the qubits. It is a subclass of QuantumGate,
      and uses a different way to initialize self.matrix'''


    def __init__(self, name: str='', qubits_t: Iterable[int]=[0], qubits_c: Iterable[int]|None=None, step: int=0, 
                 gates_t: Mapping[int,Iterable[str]]=None, control_values_c: Mapping[int,Iterable[int]]|None=None):
        _shape = 2**(len(qubits_t) + (len(qubits_c) if qubits_c is not None else 0))
        _shape = (_shape, _shape)
        _num_summands_decomposed = len(gates_t[qubits_t[0]])
        assert len(qubits_t) == len(gates_t)
        assert all([len(gt) == _num_summands_decomposed for gt in gates_t.values()])
        assert len(qubits_c) == len(control_values_c)
        assert set(qubits_t).isdisjoint(set(qubits_c)) if qubits_c is not None else True
        
        super().__init__(name=name, shape=_shape, qubits_t=qubits_t, qubits_c=qubits_c, step=step, num_summands_decomposed=_num_summands_decomposed)
        del _shape, _num_summands_decomposed
        ####
        ## Redefine atomics
        ####
        #_str_qubits_t = ''.join([str(q)+'_' for q in qubits_t])
        #_str_qubits_c = ''.join([str(q)+'_' for q in qubits_c]) if qubits_c is not None else ''
        #_iter_indices = [str(a)+str(b) for a,b in itertools.product(range(self.shape[0]), range(self.shape[1]))] # (4,4) -> ['00', '01', '02', '03', '10', '11', '12', '13', ...]
        #self.atomics = {sub: sp.symbols(self.name+'_qt'+_str_qubits_t+'qc'+_str_qubits_c+'s'+str(step)+'_p'+str(sub[0])+'_'+str(sub[1])) for sub in _iter_indices}
        #del _str_qubits_t, _str_qubits_c, _iter_indices

        _statevec_zero = np.array([[1],[0]], dtype=np.int8)
        _statevec_one  = np.array([[0],[1]], dtype=np.int8)
        _densmat_zero = np.kron(_statevec_zero, _statevec_zero.T)
        _densmat_one  = np.kron(_statevec_one, _statevec_one.T)
        _sym_statevec_zero = sp.Matrix([[1],[0]])
        _sym_statevec_one  = sp.Matrix([[0],[1]])
        _sym_densmat_zero  = sp_phy_qant.TensorProduct(_sym_statevec_zero, _sym_statevec_zero.transpose())
        _sym_densmat_one   = sp_phy_qant.TensorProduct(_sym_statevec_one, _sym_statevec_one.transpose())
        _dict_sym_densmats = {0: _sym_densmat_zero, 1: _sym_densmat_one}
        _dict_num_densmats = {0: _densmat_zero,     1: _densmat_one}
        ###
        # Overwrite matrix22_...
        ###
        _qubits_c = [] if qubits_c is None else qubits_c
        _sym_gates_t = {qt: [_get_gate_class_from_name(vv)(name=name+'_'+vv, qubits_t=[qt], qubits_c=_qubits_c, step=step) for vv in v] for qt, v in gates_t.items()}
        #self.matrix22_t         = {k: _get_gate_class_from_name(v)(name=name+'_'+v, qubits_t=[k], qubits_c=[], step=step) for k,v in gates_t.values()}
        self.matrix22_t         = {k: [vv.matrix         for vv in v] for k,v in _sym_gates_t.items()}
        self.matrix22_t_numeric = {k: [vv.matrix_numeric for vv in v]for k,v in _sym_gates_t.items()}
        
        _control_values_c = None
        if qubits_c is None or len(qubits_c)==0 and control_values_c is None:
            pass
        elif qubits_c is not None and control_values_c is None:
            #_control_values_c = {qc: [1]*self.num_summands_decomposed for qc in qubits_c} # if no control values are given, assume all 1s
            raise ValueError('If control qubits are given, control values MUST be given as well.')
        elif qubits_c is not None and control_values_c is not None:
            _control_values_c = control_values_c
        else:
            raise ValueError('If control values are given, control qubits MUST be given as well.')

        if qubits_c is None:
            pass
        else:
            self.matrix22_c         = {qc: [_dict_sym_densmats[qcvv] for qcvv in qcv] for qc, qcv in _control_values_c.items()}
            self.matrix22_c_numeric = {qc: [_dict_num_densmats[qcvv] for qcvv in qcv] for qc, qcv in _control_values_c.items()} 
            
            _merged_matrix22         = self.matrix22_t | self.matrix22_c
            _merged_matrix22_numeric = self.matrix22_t_numeric | self.matrix22_c_numeric

            
            self.matrix_decomposed = sp.Add(*[sp_phy_qant.TensorProduct(*[_merged_matrix22[k][i]         for k in sorted(_merged_matrix22.keys(),         reverse=True)]) for i in range(self.num_summands_decomposed)])
            self.matrix_numeric    =    sum( [ functools.reduce(np.kron, [_merged_matrix22_numeric[k][i] for k in sorted(_merged_matrix22_numeric.keys(), reverse=True)]) for i in range(self.num_summands_decomposed)])
            #self.matrix_numeric    =    sum( [                  np.kron(*[_merged_matrix22_numeric[k][i] for k in sorted(_merged_matrix22_numeric.keys(), reverse=True)]) for i in range(self.num_summands_decomposed)])
        #self.matrix_numeric = None
        #self.num_summands_decomposed = self._num_summands_decomposed
        #self.matrix22_t = {q_t: [None]*self.num_summands_decomposed for q_t in qubits_t}
        #self.matrix22_c = {q_c: [None]*self.num_summands_decomposed for q_c in qubits_c}
        #self.matrix22_t_numeric = {q_t: [None]*self.num_summands_decomposed for q_t in qubits_t}
        #self.matrix22_c_numeric = {q_c: [None]*self.num_summands_decomposed for q_c in qubits_c}



def _get_gate_class_from_name(gate_name: str) -> type[QuantumGate]:
    if gate_name == 'X':
        return Pauli_X_Gate
    elif gate_name == 'Y':
        return Pauli_Y_Gate
    elif gate_name == 'Z':
        return Pauli_Z_Gate
    elif gate_name == 'H':
        return Hadamard_Gate
    elif gate_name == 'H_eI':
        return Hadamard_error_I_Gate
    elif gate_name == 'I':
        return Identity_Gate
    elif gate_name == 'U':
        return U_Gate
    elif gate_name == 'GP':
        return GP_Gate
    elif gate_name == 'P':
        return P_Gate
    else:
        raise ValueError(f'Gate name {gate_name} not recognized.')
    

# class QuantumGateMultiQubit_new(QuantumGate):
#     '''Base class for multi-qubit quantum gates. It holds information of the operation that is applied to the qubits. It is a subclass of QuantumGate,
#       and uses a different way to initialize self.matrix.  
#       num_summands_decomposed, in many cases, can be inferred by 2^number-of-control-qubits, just to be flexible it is a input parameter for now.'''


#     def __init__(self, name: str='', shape:tuple[int,int]|None=[4,4], qubits_t: list[int]=[0], qubits_c: list[int]|None=None, 
#                                                                       gates_t: list[str]=['X'],   control_values_c: list[int]=[1],
#                                                                       step: int=0, _num_summands_decomposed: int|None=2):
#         super().__init__(name, shape, qubits_t, qubits_c, step, _num_summands_decomposed)
#         # To redefine: atomics, matrix, matrix22_t, matrix22_c,
#         assert len(qubits_t) == len(gates_t)
#         assert len(qubits_c)== len(control_values_c)

#         self.atomics = None
#         self.matrix = None
#         #self.matrix22_t = {q_t: [None]*self.num_summands_decomposed for q_t in qubits_t}
#         #self.matrix22_c = {q_c: [None]*self.num_summands_decomposed for q_c in qubits_c}
#         #self.matrix22_t_numeric = {q_t: [None]*self.num_summands_decomposed for q_t in qubits_t}
#         #self.matrix22_c_numeric = {q_c: [None]*self.num_summands_decomposed for q_c in qubits_c}
#         statevec_zero, statevec_one = np.array([[1],[0]], dtype=np.int8), np.array([[0],[1]], dtype=np.int8)
#         num_densmat_zero, num_densmat_one = np.kron(statevec_zero, statevec_zero.T), np.kron(statevec_one, statevec_one.T)

#         sym_statevec_zero, sym_statevec_one = sp.Matrix([[1],[0]]), sp.Matrix([[0],[1]])
#         sym_densmat_zero, sym_densmat_one = sp_phy_qant.TensorProduct(sym_statevec_zero, sym_statevec_zero.transpose()), sp_phy_qant.TensorProduct(sym_statevec_one, sym_statevec_one.transpose())
#         sym_eyes_t = {qt: Identity_Gate(name=name+'_I', qubits_t=[qt], qubits_c=qubits_c, step=step) for qt in qubits_t}
#         sym_gates_t = {qt: _get_gate_class_from_name(gate_name)(name=name+'_'+gate_name, qubits_t=[qt], qubits_c=qubits_c, step=step) for qt, gate_name in zip(qubits_t, gates_t)}

#         id_state_do = sum((vqc*(2**qc) for qc, vqc in zip(qubits_c, control_values_c)))
#         assert isinstance(self.matrix22_t, dict) and all(isinstance(self.matrix22_t[qt], list) for qt in qubits_t)
#         self.matrix22_t         = {qt: [sym_eyes_t[qt].matrix                  for s in range(self.num_summands_decomposed)] for qt in qubits_t}
#         self.matrix22_t_numeric = {qt: [np.array([[1,0],[0,1]], dtype=np.int8) for s in range(self.num_summands_decomposed)] for qt in qubits_t}
#         for qt in qubits_t:
#             self.matrix22_t[qt][id_state_do]         = sym_gates_t[qt].matrix
#             self.matrix22_t_numeric[qt][id_state_do] = sym_gates_t[qt].matrix_numeric
#         self.matrix22_c         = {qc: [sym_densmat_zero for s in range(self.num_summands_decomposed)] for qc in qubits_c}
#         self.matrix22_c_numeric = {qc: [num_densmat_zero for s in range(self.num_summands_decomposed)] for qc in qubits_c}
#         for i, vqc in enumerate(control_values_c):
#             if vqc == 1:
#                 self.matrix22_c        [qubits_c[i]][vqc*(2**qubits_c[i])] = sym_densmat_one
#                 self.matrix22_c_numeric[qubits_c[i]][vqc*(2**qubits_c[i])] = num_densmat_one
#             else: pass # if control value is 0, the matrix is already set to |0><0|
        
        

class Identity_Gate(QuantumGate):
    '''Class of the identity gate. It is a subclass of QuantumGate.'''
    def __init__(self,name: str='I', qubits_t: list[int]=[0], qubits_c: list[int]=[0], step: int=0):
        super().__init__(name=name, shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 0], [0, 1]])

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class _Identity_Gate_in_MQG(QuantumGate):
    '''Class of the identity gate in the context of multi-qubit gates. It is a subclass of QuantumGate, not Identity_Gate, as it exists in parallel.
    The difference is that this gate can be used in the context of multi-qubit gates, where the qubits are not necessarily the target qubits.
    It should not be used manually. A proper wrapper for all MQGs will be implemented in the future.'''
    def __init__(self, name: str='MQG_I', qubits_t: list[int]=[0], qubits_c: list[int]=[0], step: int=0):
        super().__init__(name=name, shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 0], [0, 1]])

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class Pauli_X_Gate(QuantumGate):
    '''Class of the Pauli X gate. It is a subclass of QuantumGate.'''
    def __init__(self, name: str='X', qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name=name, shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[0, 1], [1, 0]])

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class _Pauli_X_Gate_in_MQG(QuantumGate):
    '''Class of the Pauli X gate in the context of multi-qubit gates. It is a subclass of QuantumGate, not Pauli_X_Gate, as it exists in parallel.
    The difference is that this gate can be used in the context of multi-qubit gates, where the qubits are not necessarily the target qubits.
    It should not be used manually. A proper wrapper for all MQGs will be implemented in the future.'''
    def __init__(self, name: str='MQG_X', qubits_t: list[int]=[0], qubits_c: list[int]=[0], step: int=0):
        super().__init__(name=name, shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[0, 1], [1, 0]])

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric


class Pauli_Y_Gate(QuantumGate):
    '''Class of the Pauli Y gate. It is a subclass of QuantumGate.'''
    def __init__(self, qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name='Y', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[0, -1j], [1j, 0]])

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class Pauli_Z_Gate(QuantumGate):
    '''Class of the Pauli Z gate. It is a subclass of QuantumGate.'''
    def __init__(self, qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name='Z', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 0], [0, -1]])

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class Hadamard_Gate(QuantumGate):
    '''Class of the Hadamard gate. It is a subclass of QuantumGate.'''
    def __init__(self, qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name='H', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 1], [1, -1]])/np.sqrt(2)

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class Hadamard_error_I_Gate(QuantumGate):
    '''Class of a Hadamard gate, which fails towards identity with varyaing parameter p. 
    p=0 is 0 error, p=1 means instead of Hadamard Identity is applied. Numeric gates are equal to Hadamard gate as intented way to use is to replace symbolic values.
    It is a subclass of QuantumGate.'''
    def __init__(self, qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name='H_eI', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 1], [1, -1]])/np.sqrt(2)

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class U_Gate(QuantumGateParameterized):
    '''Class of the OpenQasm3 built-in U gate. It is a subclass of QuantumGate. https://openqasm.com/language/gates.html#built-in-gates 
    parameters must be a dict: {'theta': val, 'phi': val, 'lambda': val}  
    NOTE: There is a difference in the definition between OpenQasm2 and OpenQasm3, OQ3 = e^i(ϕ+λ)/2 OQ2'''
    def __init__(self, qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0, parameters: dict=None):
        super().__init__(name='U', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
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
    def __init__(self, qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0, parameters: dict=None):
        super().__init__(name='GP', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
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
    def __init__(self, qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0, parameters: dict=None):
        super().__init__(name='P', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
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
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0):
        super().__init__(name='CX', shape=(4,4), qubits_t=qubits_t, qubits_c=qubits_c, step=step, num_summands_decomposed=2)
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
        sym_eye_t = _Identity_Gate_in_MQG(name='CX_I', qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        sym_x_t   =  _Pauli_X_Gate_in_MQG(name='CX_X', qubits_t=qubits_t, qubits_c=qubits_c, step=step)

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
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0):
        assert len(qubits_t) == 1, 'CX gate can only be applied to a single target qubit.'
        assert len(qubits_c) == 1, 'CX gate can only be applied to a single control qubit.'
        assert qubits_t[0] != qubits_c[0], 'Control and target qubit must be different.'
        super().__init__(name='CX_new', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','X')}, control_values_c={qubits_c[0]:(0,1)})

class CCX_Gate_new(QuantumGateMultiQubit_new):
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0):
        assert len(qubits_t) == 1, 'CCX gate can only be applied to a single target qubit.'
        assert len(qubits_c) == 2, 'CCX gate must be applied to two control qubits.'
        assert qubits_t[0] not in qubits_c, 'Control and target qubits must be different.'
        super().__init__(name='CCX_new', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','I','I','X')}, control_values_c={qubits_c[0]:(0,0,1,1), qubits_c[1]:(0,1,0,1)})

class CCXX_Gate_new(QuantumGateMultiQubit_new):
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0):
        assert len(qubits_t) == 2, 'CCXX gate can only be applied to two target qubits.'
        assert len(qubits_c) == 2, 'CCXX gate must be applied to two control qubits.'
        assert qubits_t[0] not in qubits_c, 'Control and target qubits must be different.'
        super().__init__(name='CCXX_new', qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         gates_t={qubits_t[0]: ('I','I','I','X'), qubits_t[1]: ('I','I','I','X')}, 
                         control_values_c={qubits_c[0]:(0,0,1,1), qubits_c[1]:(0,1,0,1)})



class GateCollection():
    '''Class to hold lists of gates of the same type within a single quantum circuit. It is just a helper to make it more easy to traverse the circuit for evaluation.'''
    def __init__(self):
        self.allowed_gates = ['I', 'X', 'Y', 'Z', 'H', 'CX', 'CX_new', 'H_eI', 'U', 'GP', 'P', 'CCX_new', 'CCXX_new']
        self.collections = {id: [] for id in self.allowed_gates}

class QuantumCircuit():
    def __init__(self, num_qubits: int=1, num_clbits: int=1):
        self.qubits = list(reversed(list(range(num_qubits))))
        self.clbits = list(reversed(list(range(num_clbits))))
        self.gate_collection = GateCollection()
        self.barrier_collection = []
        self.steps = {} # will contain {step_number: [gate1, gate2, ...]}
        self.unitary = None
        self.allowed_gates = ['I', 'X', 'Y', 'Z', 'H', 'CX', 'CX_new', 'H_eI', 'U', 'GP', 'P', 'CCX_new', 'CCXX_new']


    def add_gate(self, name:None|str='I', qubits_t: None|list[int]=None, qubits_c: None|list[int]=None, step: None|int=0, gate: QuantumGate|None=None, parameters: dict|None=None):
        if gate is not None:
            raise ValueError('Providing a gate object is not yet implemented.')
        if qubits_t is None or qubits_t == []:
            raise ValueError('qubits_t: Target qubit(s) must be provided.')
        else:
            if name not in self.allowed_gates:
                raise ValueError('Unknown gate name')
            
            if name in self.allowed_gates:
                gate = None
                if name == 'I':
                    gate = Identity_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'X':
                    gate = Pauli_X_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'Y':
                    gate = Pauli_Y_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'Z':
                    gate = Pauli_Z_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'H':
                    gate = Hadamard_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'CX':
                    gate = CX_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'CX_new':
                    gate = CX_Gate_new(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'CCX_new':
                    gate = CCX_Gate_new(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'CCXX_new':
                    gate = CCXX_Gate_new(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'H_eI':
                    gate = Hadamard_error_I_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'U':
                    gate = U_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step, parameters=parameters)
                elif name == 'GP':
                    gate = GP_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step, parameters=parameters)
                elif name == 'P':
                    gate = P_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step, parameters=parameters)

                self.gate_collection.collections[name].append(gate)
                
                if step not in self.steps:
                    self.steps[step] = []
                self.steps[step].append(gate)
                #self.steps.sort()
                #print(self.steps)

    def add_barrier(self, step: None|int=None):
        if step is None:
            step = max(tuple(self.steps.keys()), default=0)
        bar = Barrier(step=step)
        self.barrier_collection.append(bar)
    
    #def _extend_CX(self, gates: tuple[QuantumGate]):
    #    _unitary = None
    #    for q in self.qubits:
    #        if q not in gates[0].qubits[0]:
    #            _unitary = sp_phy_qant.TensorProduct(_unitary, sp.eye(2)) if _unitary is not None else sp.eye(2)
    #        elif q in gates[0].qubits[0]:
    #            _unitary = sp_phy_qant.TensorProduct(_unitary, gates[0].matrix) if _unitary is not None else gates[0].matrix
    #    return _unitary
    
    #def extend_by_identities(self, gates: tuple[QuantumGate]):
    #    _unitary = None
    #    if len(gates) > 1 and len(self.qubits) == len(gates[0].qubits):
    #        raise ValueError('If single gate operates on all qubits it requires an exclusive timestep, i.e. can not be applied with other gates in the same step.')
    #    if len(gates) == 1:
    #        if len(self.qubits) == len(gates[0].qubits_c+gates[0].qubits_t):
    #            if gates[0].name != 'CX':
    #                return gates[0].matrix
    #            elif gates[0].name == 'CX':
    #                _u0 = None
    #                _u1 = None
    #                for q in self.qubits:
    #                    _tmp0 = None
    #                    _tmp1 = None
    #                    if q == gates[0].qubits_c[0]:
    #                        _tmp0 = sp.Matrix([[1, 0], [0, 0]])
    #                        _tmp1 = sp.Matrix([[0, 0], [0, 1]])
    #                    elif q == gates[0].qubits_t[0]:
    #                        _tmp0 = sp.eye(2)
    #                        _tmp1 = sp.Matrix([[0, 1], [1, 0]])
    #                        #_tmp0 = gates[0].matrix[0:2, 0:2]
    #                        #_tmp1 = gates[0].matrix[2:4, 2:4]
    #                    else:
    #                        _tmp0 = sp.eye(2)
    #                        _tmp1 = sp.eye(2)
    #                    _u0 = sp_phy_qant.TensorProduct(_u0, _tmp0) if _u0 is not None else _tmp0
    #                    _u1 = sp_phy_qant.TensorProduct(_u1, _tmp1) if _u1 is not None else _tmp1
    #                return _u0 + _u1
    #            else: 
    #                raise ValueError('Something went wrong.')
    #        elif len(self.qubits) > len(gates[0].qubits_c+gates[0].qubits_t):
    #            for q in self.qubits:
    #                if q not in gates[0].qubits:
    #                    _unitary = sp_phy_qant.TensorProduct(_unitary, sp.eye(2)) if _unitary is not None else sp.eye(2)
    #                else:
    #                    _unitary = sp_phy_qant.TensorProduct(_unitary, gates[0].matrix) if _unitary is not None else gates[0].matrix
    #            return _unitary
    #    elif len(gates) > 1:
    #        gates_sorted_by_qubits = sorted(gates, key=lambda x: x.qubits[0])
    #        print('gates_sorted_by_qubits:', gates_sorted_by_qubits)
    #        qubits_gates = [g.qubits[0][0] for g in gates_sorted_by_qubits]
    #        for q in self.qubits:
    #            if q not in qubits_gates:
    #                print(f'qubit {q} not in gates')
    #                _unitary = sp_phy_qant.TensorProduct(_unitary, sp.eye(2)) if _unitary is not None else sp.eye(2)
    #            elif q in qubits_gates:
    #                for gate in gates_sorted_by_qubits:
    #                    print('gate:', gate.name, 'qubits:', gate.qubits)
    #                    if q in gate.qubits[0]:
    #                        print('    gate:', gate.name, 'qubits:', gate.qubits)
    #                        _unitary = sp_phy_qant.TensorProduct(_unitary, gate.matrix) if _unitary is not None else gate.matrix
    #            else:
    #                raise ValueError('Something went wrong.')
    #
    #        return _unitary
    #    raise ValueError('Something went wrong.')
    
    def assemble_symbolic_unitary(self):
        '''Assemble the unitary matrix of the quantum circuit.'''
        self.unitary = sp.eye(2**len(self.qubits))
        for step in self.steps.keys():# step is step number of timesteps
            unitary_step = sp.zeros(2**len(self.qubits))
            gates_step = self.steps[step]
            #print(gates_step)
            num_sum_max = max([g.num_summands_decomposed for g in gates_step])
            for i in range(num_sum_max):
                gates = [None]*len(self.qubits)
                for gate in gates_step:
                    if gate.num_summands_decomposed >= i:
                        #q_t = gate.qubits_t[0]
                        #gates[q_t] = gate.matrix22_t[q_t][i]
                        for q_t in gate.qubits_t:
                            gates[q_t] = gate.matrix22_t[q_t][i]
                        #print(gates[q_t])
                        if gate.qubits_c is not None:
                            #q_c = gate.qubits_c[0]
                            #gates[q_c] = gate.matrix22_c[q_c][i]
                            for q_c in gate.qubits_c:
                                gates[q_c] = gate.matrix22_c[q_c][i]
                for j in range(len(gates)):
                    if gates[j] is None:
                        #gates[j] = Identity_Gate(qubits_t=[j], qubits_c=None, step=step).matrix22_t[j][0]
                        gates[j] = sp.eye(2)
                gates = list(reversed(gates)) # reverse the list to match the order of the qubits (little endian) for the tensor product
                _tmp_unitary = gates[0]
                for j in range(1, len(gates)):
                    #print(j, gates[j])
                    _tmp_unitary = sp_phy_qant.TensorProduct(_tmp_unitary, gates[j])
                #display(_tmp_unitary)
                unitary_step += _tmp_unitary
            self.unitary = unitary_step @ self.unitary
    
    def subs_symbolic_zeros_in_symbolic_unitary(self):
        '''Substitute entries in the symbolic unitary that are zero with zeros.'''
        # this function based create_numeric_unitary_from_symbolic, so if making changes to either you might want to adapt the other as well 
        for gate_type, gates in self.gate_collection.collections.items():
            for gate in gates:
                for i in range(gate.num_summands_decomposed):
                    symbs_t = gate.matrix22_t[gate.qubits_t[0]][i].flat()
                    nums_t = gate.matrix22_t_numeric[gate.qubits_t[0]][i].flatten()
                    for s, n in zip(symbs_t, nums_t):
                        if n == 0:
                            self.unitary = self.unitary.subs(s, n)
                    
                    if gate.qubits_c is not None:
                        symbs_c = gate.matrix22_c[gate.qubits_c[0]][i].flat() 
                        nums_c = gate.matrix22_c_numeric[gate.qubits_c[0]][i].flatten()
                        for s, n in zip(symbs_c, nums_c):
                            if n == 0:
                                self.unitary = self.unitary.subs(s, n)
    
    @staticmethod
    def _filter_gates_to_replace_with_alternatives(gate_collection, gate_identifications):
        generator_all_gates = (gate for collection in gate_collection.collections.values() for gate in collection)
        if gate_identifications is None:
            list_gates_to_replace = [gate for gate in generator_all_gates]
        else:
            generator_all_identifiers = ((gate_identifications['steps'][i], gate_identifications['names'][i], gate_identifications['qubits_t'][i]) for i in range(len(gate_identifications['steps'])))
            list_gates_to_replace = [it[0] for it in itertools.product(generator_all_gates, generator_all_identifiers) if (it[0].step == it[1][0] and it[0].name == it[1][1] and it[1][2] in it[0].qubits_t)]
        return [gate for gate in list_gates_to_replace if hasattr(gate, 'matrix_alt')]

    def subs_symbolic_alternatives_in_symbolic_unitary(self, gate_identifications: Mapping[str, Iterable]|None=None):
        '''Substitute entries in the symbolic unitary that have alternative representations. This is typically the case for parametric gates.  
        The mapping is optional, if not provided all alternatives will be substituted.  
        If provided, the mapping must identify each gate for which the alternative representation should be used, by step, name, and target qubits.
        Multiple gates can be identified by providing lists/tuples of the same length for each key.  
        For now the keys of the mapping MUST be 'steps', 'names', and 'qubits_t. This might be changed in the future to allow more flexible identification.
        Example: gate_identifications = {'step': (0, 1), 'name': (CX, U), 'qubits_t': (1, 1)}  
        would substitute the alternatives for a CX gate at step 0 targeting qubit 1, and a U gate at step 1 targeting qubit 1.  
        It is sufficient to provide only the target qubit (and only 1 target qubit in case of a multi-qubit gate) as a qubit can only have applied 1 gate per timestep.'''
        
        # this function based create_numeric_unitary_from_symbolic, so if making changes to either you might want to adapt the other as well 

        if gate_identifications is not None:
            num_to_replace = len(gate_identifications['steps'])
            assert all(len(v) == num_to_replace for v in gate_identifications.values()), 'All values in gate_identifications must have the same length.'

            list_gates_to_replace = self._filter_gates_to_replace_with_alternatives(self.gate_collection, gate_identifications)
            len_list_gates_to_replace = len(list_gates_to_replace)
            if len_list_gates_to_replace != num_to_replace:
                print(f'Warning: Number of gates found to replace -{len_list_gates_to_replace}- does not match number of identifications provided -{num_to_replace}-')
        else:
            list_gates_to_replace = self._filter_gates_to_replace_with_alternatives(self.gate_collection, gate_identifications)
        self.unitary = self.unitary.subs({v[0]: v[1] for gate in list_gates_to_replace for v in zip(gate.matrix.flat(), gate.matrix_alt.flat())})

        #for i in range(num_to_replace):
        #    step = gate_identifications['steps'][i]
        #    name = gate_identifications['names'][i]
        #    qubit_t = gate_identifications['qubits_t'][i]
        #    
        #    for gate_type, gates in self.gate_collection.collections.items():
        #        for gate in gates:
        #            if gate.step != step or gate.name != name or qubit_t not in gate.qubits_t:
        #                raise ValueError('Gate not found in the circuit to replace with alternative representation:', 'step', step, 'name', name, 'qubit_t', qubit_t)
        #                
        #            if hasattr(gate, 'matrix_alt'):
        #                self.unitary = self.unitary.subs(zip(gate.matrix.flat(), gate.matrix_alt.flat()))
        #            else:
        #                raise ValueError('Gate does not have an alternative representation:', 'step', step, 'name', name, 'qubit_t', qubit_t)


    def create_numeric_unitary_from_symbolic(self):
        # subs_symbolic_zeros_in_symbolic_unitary based on this function, so if making changes to either you might want to adapt the other as well 
        self.unitary_numeric = self.unitary.copy()
        for gate_type, gates in self.gate_collection.collections.items():
            for gate in gates:
                for i in range(gate.num_summands_decomposed):
                    symbs_t = [e for q_t in gate.qubits_t for e in gate.matrix22_t[q_t][i].flat()]
                    nums_t = [e for q_t in gate.qubits_t for e in gate.matrix22_t_numeric[q_t][i].flatten()]
                    for s, n in zip(symbs_t, nums_t):
                        self.unitary_numeric = self.unitary_numeric.subs(s, n)
                    
                    if gate.qubits_c is not None:
                        symbs_c = [e for q_c in gate.qubits_c for e in gate.matrix22_c[q_c][i].flat()]
                        nums_c = [e for q_c in gate.qubits_c for e in gate.matrix22_c_numeric[q_c][i].flatten()]
                        for s, n in zip(symbs_c, nums_c):
                            self.unitary_numeric = self.unitary_numeric.subs(s, n)
    
class QuantumState():
    def __init__(self, num_qubits: int=1, method: str='statevector'):
        self.qubits = list(reversed(list(range(num_qubits))))
        self.state = None
        self.method = method
        self.allowed_methods = ['statevector', 'density_matrix']
        if method not in self.allowed_methods:
            raise ValueError('Unknown method')
        if method == 'statevector':
            self.symbols = sp.symbols(f'phi:{2**len(self.qubits)}', complex=True, extended_real=False)
            self.state = sp.Matrix(self.symbols, shape=(2**len(self.qubits), 1))
            #print(self.state)
            #print(self.state.shape)
            
    def set_state(self, state: dict[int:int]):
        if self.method == 'statevector':
            if len(state) != len(self.qubits):
                raise ValueError('State vector does not match number of qubits.')
            for q in self.qubits:
                v = state[q]
                if v not in [0, 1]:
                    raise ValueError('State vector must be binary.')
                
                if q == len(self.qubits)-1:
                    self.state = sp.Matrix([1-v, v], shape=(2, 1))
                else:
                    self.state = sp_phy_qant.TensorProduct(self.state, sp.Matrix([1-v, v], shape=(2, 1)))


    
def _map_qasmgatename_to_qcgatename(qasm_gate:str=None):
    if qasm_gate == 'x':
        return 'X'
    elif qasm_gate == 'y':
        return 'Y'
    elif qasm_gate == 'z':
        return 'Z'
    elif qasm_gate == 'h':
        return 'H'
    elif qasm_gate == 'cx':
        return 'CX'
    else:
        raise ValueError('Unknown qasm gate name:', qasm_gate)
    
def _iterate_over_qasm_statements(qpf:openqasm3.ast.Program, qc: QuantumCircuit, timestep:int):
    gate_name_mapping_qasm_qc = {'x': 'X', 'y': 'Y', 'z': 'Z', 'h': 'H', 'cx': 'CX'}

    for statement in qpf.statements:
        if isinstance(statement, openqasm3.ast.Include):
            print('Included filenames:', statement.filename)
        elif isinstance(statement, openqasm3.ast.ClassicalDeclaration):
            if isinstance(statement.type, openqasm3.ast.BitType):
                cbit_type = 'bit'
            else:
                raise ValueError('Unknown type', 'statement.type:', statement.type)
            cbit_name = statement.identifier.name
            cbit_length = statement.type.size.value
            if cbit_length > len(qc.clbits):
                qc.clbits = list(reversed(list(range(cbit_length))))
            #print('Classical bit name:', cbit_name)
            #print('Classical bit type:', cbit_type)
            #print('Classical bit length:', cbit_length)
        elif isinstance(statement, openqasm3.ast.QubitDeclaration):
            qbit_name = statement.qubit.name
            qbit_length = statement.size.value
            if qbit_length > len(qc.qubits):
                qc.qubits = list(reversed(list(range(qbit_length))))
            #print('Quantum bit name:', qbit_name)
            #print('Quantum bit length:', qbit_length)
        elif isinstance(statement, openqasm3.ast.QuantumGate):
            qgate_name = statement.name.name
            qgate_qbit_names = [statement.qubits[i].name.name for i in range(len(statement.qubits))]
            qgate_qbit_indices = [statement.qubits[i].indices[0][0].value for i in range(len(statement.qubits))] # why is there doubly nested list in indices?
            #print('Quantum gate name:', qgate_name)
            #print('Quantum gate qubits names:', qgate_qbit_names)
            #print('Quantum gate qubits indices:', qgate_qbit_indices)
            #qc_gate_name = _map_qasmgatename_to_qcgatename(qgate_name)
            qc_gate_name = gate_name_mapping_qasm_qc[qgate_name]
            if len (qgate_qbit_indices) == 1: # single qubit gate
                qc.add_gate(name=qc_gate_name, qubits_t=qgate_qbit_indices, step=timestep)
            elif len (qgate_qbit_indices) == 2: # single qubit gate
                qc.add_gate(name=qc_gate_name, qubits_c=[qgate_qbit_indices[0]], qubits_t=[qgate_qbit_indices[1]], step=timestep)
        elif isinstance(statement, openqasm3.ast.QuantumBarrier):
            #print(statement)
            qbarrier_name = 'barrier' # name is not stored in statement 
            qbarrier_qbit_names = [statement.qubits[i].name.name for i in range(len(statement.qubits))]
            qbarrier_qbit_indices = [statement.qubits[i].indices[0][0].value for i in range(len(statement.qubits))] # why is there doubly nested list in indices?
            #print('Quantum gate name:', qbarrier_name)
            #print('Quantum gate qubits names:', qbarrier_qbit_names)
            #print('Quantum gate qubits indices:', qbarrier_qbit_indices)
            qc.add_barrier(step=timestep)
            timestep += 1
        elif isinstance(statement, openqasm3.ast.QuantumMeasurementStatement):
            #print(statement)
            qmeasurement_name = 'measurement' # name is not stored in statement 
            qmeasurement_qbit_name = statement.measure.qubit.name.name
            qmeasurement_qbit_index = statement.measure.qubit.indices[0][0].value
            qmeasurement_cbit_name = statement.target.name.name
            qmeasurement_cbit_index = statement.target.indices[0][0].value
            #print('Quantum measurement name:', qmeasurement_name)
            #print('Quantum measurement qubit name:', qmeasurement_qbit_name)
            #print('Quantum measurement qubit index:', qmeasurement_qbit_index)
            #print('Quantum measurement cbit name:', qmeasurement_cbit_name)
            #print('Quantum measurement cbit index:', qmeasurement_cbit_index)
    
    return qc

def openqasm3_to_qc(filepath: pathlib.Path, timestep:int|None=None, qc: QuantumCircuit|None=None, qbit_mapping:dict[int:int]|None=None, cbit_mapping:dict[int:int]|None=None):
    '''Function to read an OpenQASM 3 file and convert it to a QuantumCircuit object. If a QuantumCircuit object is provided, the circuit from the
    OpenQASM 3 file is appended to it.
    In that case a qubit and cbit mapping is required as dicts: {*bitnumber qc: *bitnumber qasm}.
    If a timestep is provided, it will overwrite the automatically determined value, which is 0 if there is no qc provided, or last timestep of the qc +1.
    '''
    def _qc_qbit(q_qasm, mapping=qbit_mapping):
        return qbit_mapping[q_qasm]
    def _qc_cbit(c_qasm, mapping=cbit_mapping):
        return cbit_mapping[c_qasm]
    
    if qc is None:
        qc = QuantumCircuit()
    if timestep is None:
        timestep: int = 0
    else:
        timestep = max(tuple(qc.steps.keys())) + 1
        qc.add_barrier(step = timestep-1)
    
    with open(filepath, 'r') as f:
        qpf = f.read()
    
    qpf = openqasm3.parse(qpf)
    
    qc = _iterate_over_qasm_statements(qpf, qc, timestep)
    return qc