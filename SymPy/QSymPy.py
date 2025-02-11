import sympy as sp
import itertools
import numpy as np
import sympy.physics
import sympy.physics.quantum


class QuantumGate:
    '''Base class for quantum gates. It holds information of the operation that is applied to the qubit(s).
    It does not hold the matrix for the whole quantum circuit, but just for the qubit(s) it's supposed to change the state of.
    Lists of qubits and matrices are assumed to be encoded as most significant qubit first, i.e. left and bottom, e.g. |q2 q1 q0⟩ = |q2> kron |q1> kron |q0>, 
    |0 0 1> = (1,0) kron (1,0) kron (0,1) = (0,1,0,0,0,0,0,0).'''
    # We need information about the qubit and step here, because otherwise the atomic symbols could not be distinguished from those of other gates
    def __init__(self, name: str='', shape:tuple[int,int]=[2,2], qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        self.name = name
        self.shape = shape
        iter_indices = [str(sub[0])+str(sub[1]) for sub in itertools.product(range(shape[0]), range(shape[1]))] # (2,2) -> ['00', '01', '10', '11']
        str_qubits_t = ''.join([str(q)+'_' for q in qubits_t])
        str_qubits_c = ''.join([str(q)+'_' for q in qubits_c]) if qubits_c is not None else ''
        self.atomics = {sub: sp.symbols(self.name+'_qt'+str_qubits_t+'qc'+str_qubits_c+'s'+str(step)+'_p'+sub) for sub in iter_indices} # {'00': I_qt0_qc0_s0_p00, '01': I_qt0_qc0_s0_p01, ...}
        self.matrix = sp.Matrix([[self.atomics[col] for col in itertools.islice(iter_indices, row*shape[1], row*shape[1]+shape[1])] for row in range(shape[0])]) # [[I_..._p00, I_..._p01], [I_..._p10, I_..._p11]]
        self.matrix_numeric = None
        self.qubits = [qubits_t] if qubits_c is None else [qubits_t, qubits_c]
        self.qubits_t = qubits_t
        self.qubits_c = qubits_c
        self.step = step
        self.num_summands_decomposed = 1
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

class QuantumGateMultiQubit(QuantumGate):
    '''Base class for multi-qubit quantum gates. It holds information of the operation that is applied to the qubits. It is a subclass of QuantumGate,
      and uses a different way to initialize self.matrix'''
    def __init__(self, name: str='', shape:tuple[int,int]=[4,4], qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name, shape, qubits_t, qubits_c, step)
        self.atomics = None
        self.matrix = None
        self.matrix_numeric = None

        self.num_summands_decomposed = 2
        self.matrix22_t = {q_t: [None]*self.num_summands_decomposed for q_t in qubits_t}
        self.matrix22_c = {q_c: [None]*self.num_summands_decomposed for q_c in qubits_c}
        self.matrix22_t_numeric = {q_t: [None]*self.num_summands_decomposed for q_t in qubits_t}
        self.matrix22_c_numeric = {q_c: [None]*self.num_summands_decomposed for q_c in qubits_c}
            



class Identity_Gate(QuantumGate):
    '''Class of the identity gate. It is a subclass of QuantumGate.'''
    def __init__(self, qubits_t: list[int]=[0], qubits_c: list[int]=[0], step: int=0):
        super().__init__(name='I', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
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
    def __init__(self,  qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name='X', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
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

class Pauli_Z_Gate(QuantumGate):
    '''Class of the Pauli Z gate. It is a subclass of QuantumGate.'''
    def __init__(self, qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name='Z', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 0], [0, -1]])

class Hadamard_Gate(QuantumGate):
    '''Class of the Hadamard gate. It is a subclass of QuantumGate.'''
    def __init__(self, qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0):
        super().__init__(name='H', shape=(2,2), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        self.matrix_numeric = np.array([[1, 1], [1, -1]])/np.sqrt(2)

        self.matrix22_t[qubits_t[0]][0] = self.matrix
        self.matrix22_t_numeric[qubits_t[0]][0] = self.matrix_numeric

class CNOT_Gate(QuantumGateMultiQubit):
    '''Class of the CNOT gate. It is a subclass of QuantumGateMultiQubit.
    For CNOT_c_t the matrix is |0_c⟩⟨0_c|⊗I_t+|1_c⟩⟨1_c|⊗X_t. For CNOT_t_c the matrix is I_t⊗|0_c⟩⟨0_c|+X_t⊗|1_c⟩⟨1_c|.
    For multiple qubits surrounding and/or between c and t qubit, Identity gates are inserted, e.g. for CNOT_3_c_1_t
    the matrix is I_3⊗|0_c⟩⟨0_c|⊗I_1⊗I_t+I_3⊗|1_c⟩⟨1_c|⊗I_1⊗X_1.'''
    def __init__(self, qubits_t: tuple|list[int]=[0], qubits_c: tuple|list[int]=None, step: int=0):
        super().__init__(name='CNOT', shape=(4,4), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        assert len(qubits_t) == 1, 'CNOT gate can only be applied to a single target qubit.'
        assert len(qubits_c) == 1, 'CNOT gate can only be applied to a single control qubit.'
        assert qubits_t[0] != qubits_c[0], 'Control and target qubit must be different.'

        self.num_summands_decomposed = 2

        iter_indices_22 = [str(sub[0])+str(sub[1]) for sub in itertools.product(range(2), range(2))] # (2,2) -> ['00', '01', '10', '11']
        iter_name_ext_orig = ['0c', '1c', 'It', 'Xt']
        iter_name_ext = [i[0]+i[1] for i in itertools.product(iter_name_ext_orig, iter_indices_22)] # ['0c00', '0c01', ... , 'Xt10', 'Xt11']
        iter_indices = [str(sub[0])+str(sub[1]) for sub in itertools.product(range(self.shape[0]), range(self.shape[1]))] # (4,4) -> ['00', '01', '02', '03', '10', '11', '12', '13', ...]
        #print(iter_name_ext)
        str_qubits_t = ''.join([str(q)+'_' for q in qubits_t])
        str_qubits_c = ''.join([str(q)+'_' for q in qubits_c]) if qubits_c is not None else ''
        self.atomics = {sub: sp.symbols(self.name+nex+'_qt'+str_qubits_t+'qc'+str_qubits_c+'s'+str(step)+'_p'+sub) for nex, sub in zip(iter_name_ext, iter_indices)}
        #self.atomics22 = {sub: sp.symbols(self.name+nex+'_qt'+str_qubits_t+'qc'+str_qubits_c+'s'+str(step)) for nex, sub in zip(iter_name_ext_orig, iter_indices_22)}
        #self.matrices22 = {}
        #atomics_list = list(self.atomics.values())
        #for i, nex in enumerate(iter_name_ext_orig):
        #    idx = i*4
        #    self.matrices22[nex] = sp.Matrix([atomics_list[idx+0:idx+2],atomics_list[idx+2:idx+4]])
        del iter_indices_22, iter_name_ext_orig, iter_name_ext, iter_indices, str_qubits_t, str_qubits_c


        statevec_zero = np.array([[1],[0]], dtype=np.int8)
        statevec_one = np.array([[0],[1]], dtype=np.int8)
        densmat_zero = np.kron(statevec_zero, statevec_zero.T)
        densmat_one = np.kron(statevec_one, statevec_one.T)
        sym_statevec_zero = sp.Matrix([[1],[0]])
        sym_statevec_one = sp.Matrix([[0],[1]])
        sym_densmat_zero = sp.physics.quantum.TensorProduct(sym_statevec_zero, sym_statevec_zero.transpose())
        sym_densmat_one = sp.physics.quantum.TensorProduct(sym_statevec_one, sym_statevec_one.transpose())
        sym_eye_t = _Identity_Gate_in_MQG(name='CNOT_I', qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        sym_x_t = _Pauli_X_Gate_in_MQG(name='CNOT_X', qubits_t=qubits_t, qubits_c=qubits_c, step=step)

        self.matrix22_t[qubits_t[0]][0] = sym_x_t.matrix
        self.matrix22_c[qubits_c[0]][0] = sym_densmat_one
        self.matrix22_t[qubits_t[0]][1] = sym_eye_t.matrix
        self.matrix22_c[qubits_c[0]][1] = sym_densmat_zero

        self.matrix22_t_numeric[qubits_t[0]][0] = np.array([[0, 1], [1, 0]],dtype=np.int8)
        self.matrix22_c_numeric[qubits_c[0]][0] = densmat_one
        self.matrix22_t_numeric[qubits_t[0]][1] = np.array([[1, 0], [0, 1]],dtype=np.int8)
        self.matrix22_c_numeric[qubits_c[0]][1] = densmat_zero

        if qubits_c[0] > qubits_t[0]:
            m1 = sp.physics.quantum.TensorProduct(sym_densmat_one, sym_x_t.matrix)
            m2 = sp.physics.quantum.TensorProduct(sym_densmat_zero, sym_eye_t.matrix)
            #display(m1)
            #display(m2)
            self.matrix = m1 + m2
            self.matrix_numeric = np.kron(densmat_one, np.array([[0, 1], [1, 0]],dtype=np.int8)) + np.kron(densmat_zero, np.eye(2,dtype=np.int8))
            # = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
        elif qubits_c[0] < qubits_t[0]:
            m1 = sp.physics.quantum.TensorProduct(sym_x_t.matrix, sym_densmat_one)
            m2 = sp.physics.quantum.TensorProduct(sym_eye_t.matrix, sym_densmat_zero)
            #display(m1)
            #display(m2)
            self.matrix = m1 + m2
            self.matrix_numeric = np.kron(np.array([[0, 1], [1, 0]],dtype=np.int8), densmat_one) + np.kron(np.eye(2,dtype=np.int8), densmat_zero)
            # = np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0])
        del statevec_zero, statevec_one, densmat_zero, densmat_one, sym_statevec_zero, sym_statevec_one, sym_densmat_zero, sym_densmat_one, sym_eye_t, sym_x_t




class GateCollection():
    '''Class to hold lists of gates of the same type within a single quantum circuit. It is just a helper to make it more easy to traverse the circuit for evaluation.'''
    def __init__(self):
        self.allowed_gates = ['I', 'X', 'Y', 'Z', 'H', 'CNOT']
        self.collections = {id: [] for id in self.allowed_gates}

class QuantumCircuit():
    def __init__(self, num_qubits: int=1, num_clbits: int=1):
        self.qubits = list(reversed(list(range(num_qubits))))
        self.clbits = list(reversed(list(range(num_clbits))))
        self.gate_collection = GateCollection()
        self.steps = {} # will contain {step_number: [gate1, gate2, ...]}
        self.unitary = None
        self.allowed_gates = ['I', 'X', 'Y', 'Z', 'H', 'CNOT']

    def add_gate(self, name:None|str='I', qubits_t: None|list[int]=[0], qubits_c: None|list[int]=None, step: None|int=0, gate: QuantumGate|None=None):
        if gate is not None:
            raise ValueError('Providing a gate object is not yet implemented.')
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
                elif name == 'CNOT':
                    gate = CNOT_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                
                self.gate_collection.collections[name].append(gate)
                
                if step not in self.steps:
                    self.steps[step] = []
                self.steps[step].append(gate)
                #self.steps.sort()
                #print(self.steps)
            
    
    def _extend_CNOT(self, gates: tuple[QuantumGate]):
        _unitary = None
        for q in self.qubits:
            if q not in gates[0].qubits[0]:
                _unitary = sp.physics.quantum.TensorProduct(_unitary, sp.eye(2)) if _unitary is not None else sp.eye(2)
            elif q in gates[0].qubits[0]:
                _unitary = sp.physics.quantum.TensorProduct(_unitary, gates[0].matrix) if _unitary is not None else gates[0].matrix
        return _unitary
    def extend_by_identities(self, gates: tuple[QuantumGate]):
        _unitary = None
        if len(gates) > 1 and len(self.qubits) == len(gates[0].qubits):
            raise ValueError('If single gate operates on all qubits it requires an exclusive timestep, i.e. can not be applied with other gates in the same step.')
        if len(gates) == 1:
            if len(self.qubits) == len(gates[0].qubits_c+gates[0].qubits_t):
                if gates[0].name != 'CNOT':
                    return gates[0].matrix
                elif gates[0].name == 'CNOT':
                    _u0 = None
                    _u1 = None
                    for q in self.qubits:
                        _tmp0 = None
                        _tmp1 = None
                        if q == gates[0].qubits_c[0]:
                            _tmp0 = sp.Matrix([[1, 0], [0, 0]])
                            _tmp1 = sp.Matrix([[0, 0], [0, 1]])
                        elif q == gates[0].qubits_t[0]:
                            _tmp0 = sp.eye(2)
                            _tmp1 = sp.Matrix([[0, 1], [1, 0]])
                            #_tmp0 = gates[0].matrix[0:2, 0:2]
                            #_tmp1 = gates[0].matrix[2:4, 2:4]
                        else:
                            _tmp0 = sp.eye(2)
                            _tmp1 = sp.eye(2)
                        _u0 = sp.physics.quantum.TensorProduct(_u0, _tmp0) if _u0 is not None else _tmp0
                        _u1 = sp.physics.quantum.TensorProduct(_u1, _tmp1) if _u1 is not None else _tmp1
                    return _u0 + _u1
                else: 
                    raise ValueError('Something went wrong.')
            elif len(self.qubits) > len(gates[0].qubits_c+gates[0].qubits_t):
                for q in self.qubits:
                    if q not in gates[0].qubits:
                        _unitary = sp.physics.quantum.TensorProduct(_unitary, sp.eye(2)) if _unitary is not None else sp.eye(2)
                    else:
                        _unitary = sp.physics.quantum.TensorProduct(_unitary, gates[0].matrix) if _unitary is not None else gates[0].matrix
                return _unitary
        elif len(gates) > 1:
            gates_sorted_by_qubits = sorted(gates, key=lambda x: x.qubits[0])
            print('gates_sorted_by_qubits:', gates_sorted_by_qubits)
            qubits_gates = [g.qubits[0][0] for g in gates_sorted_by_qubits]
            for q in self.qubits:
                if q not in qubits_gates:
                    print(f'qubit {q} not in gates')
                    _unitary = sp.physics.quantum.TensorProduct(_unitary, sp.eye(2)) if _unitary is not None else sp.eye(2)
                elif q in qubits_gates:
                    for gate in gates_sorted_by_qubits:
                        print('gate:', gate.name, 'qubits:', gate.qubits)
                        if q in gate.qubits[0]:
                            print('    gate:', gate.name, 'qubits:', gate.qubits)
                            _unitary = sp.physics.quantum.TensorProduct(_unitary, gate.matrix) if _unitary is not None else gate.matrix
                else:
                    raise ValueError('Something went wrong.')

            return _unitary
        raise ValueError('Something went wrong.')
    
    def assemble_unitary(self):
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
                        q_t = gate.qubits_t[0]
                        gates[q_t] = gate.matrix22_t[q_t][i]
                        #print(gates[q_t])
                        if gate.qubits_c is not None:
                            q_c = gate.qubits_c[0]
                            gates[q_c] = gate.matrix22_c[q_c][i]
                for j in range(len(gates)):
                    if gates[j] is None:
                        gates[j] = Identity_Gate(qubits_t=[j], qubits_c=None, step=step).matrix22_t[j][0]
                        gates[j] = sp.eye(2)
                gates = list(reversed(gates)) # reverse the list to match the order of the qubits (little endian) for the tensor product
                _tmp_unitary = gates[0]
                for j in range(1, len(gates)):
                    #print(j, gates[j])
                    _tmp_unitary = sp.physics.quantum.TensorProduct(_tmp_unitary, gates[j])
                #display(_tmp_unitary)
                unitary_step += _tmp_unitary
            self.unitary = unitary_step @ self.unitary

    def create_numeric_unitary_from_symbolic(self):
        self.unitary_numeric = self.unitary.copy()
        for gate_type, gates in self.gate_collection.collections.items():
            for gate in gates:
                for i in range(gate.num_summands_decomposed):
                    symbs_t = gate.matrix22_t[gate.qubits_t[0]][i].flat()
                    nums_t = gate.matrix22_t_numeric[gate.qubits_t[0]][i].flatten()
                    for s, n in zip(symbs_t, nums_t):
                        self.unitary_numeric = self.unitary_numeric.subs(s, n)
                    
                    if gate.qubits_c is not None:
                        symbs_c = gate.matrix22_c[gate.qubits_c[0]][i].flat() 
                        nums_c = gate.matrix22_c_numeric[gate.qubits_c[0]][i].flatten()
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
            symbols = sp.symbols(f'phi:{2**len(self.qubits)}', complex=True, extended_real=False)
            self.state = sp.Matrix(symbols, shape=(2**len(self.qubits), 1))
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
                    self.state = sp.physics.quantum.TensorProduct(self.state, sp.Matrix([1-v, v], shape=(2, 1)))
