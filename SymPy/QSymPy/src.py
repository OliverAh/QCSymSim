import pathlib
import sympy as sp
import itertools
import numpy as np
import sympy.physics
import sympy.physics.quantum
import openqasm3


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
        self.ids_matrix_zeros = None
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
        super().__init__(name='CX', shape=(4,4), qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        assert len(qubits_t) == 1, 'CX gate can only be applied to a single target qubit.'
        assert len(qubits_c) == 1, 'CX gate can only be applied to a single control qubit.'
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
        sym_eye_t = _Identity_Gate_in_MQG(name='CX_I', qubits_t=qubits_t, qubits_c=qubits_c, step=step)
        sym_x_t = _Pauli_X_Gate_in_MQG(name='CX_X', qubits_t=qubits_t, qubits_c=qubits_c, step=step)

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
        self.allowed_gates = ['I', 'X', 'Y', 'Z', 'H', 'CX']
        self.collections = {id: [] for id in self.allowed_gates}

class QuantumCircuit():
    def __init__(self, num_qubits: int=1, num_clbits: int=1):
        self.qubits = list(reversed(list(range(num_qubits))))
        self.clbits = list(reversed(list(range(num_clbits))))
        self.gate_collection = GateCollection()
        self.barrier_collection = []
        self.steps = {} # will contain {step_number: [gate1, gate2, ...]}
        self.unitary = None
        self.allowed_gates = ['I', 'X', 'Y', 'Z', 'H', 'CX']

    def add_gate(self, name:None|str='I', qubits_t: None|list[int]=None, qubits_c: None|list[int]=None, step: None|int=0, gate: QuantumGate|None=None):
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
    #            _unitary = sp.physics.quantum.TensorProduct(_unitary, sp.eye(2)) if _unitary is not None else sp.eye(2)
    #        elif q in gates[0].qubits[0]:
    #            _unitary = sp.physics.quantum.TensorProduct(_unitary, gates[0].matrix) if _unitary is not None else gates[0].matrix
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
    #                    _u0 = sp.physics.quantum.TensorProduct(_u0, _tmp0) if _u0 is not None else _tmp0
    #                    _u1 = sp.physics.quantum.TensorProduct(_u1, _tmp1) if _u1 is not None else _tmp1
    #                return _u0 + _u1
    #            else: 
    #                raise ValueError('Something went wrong.')
    #        elif len(self.qubits) > len(gates[0].qubits_c+gates[0].qubits_t):
    #            for q in self.qubits:
    #                if q not in gates[0].qubits:
    #                    _unitary = sp.physics.quantum.TensorProduct(_unitary, sp.eye(2)) if _unitary is not None else sp.eye(2)
    #                else:
    #                    _unitary = sp.physics.quantum.TensorProduct(_unitary, gates[0].matrix) if _unitary is not None else gates[0].matrix
    #            return _unitary
    #    elif len(gates) > 1:
    #        gates_sorted_by_qubits = sorted(gates, key=lambda x: x.qubits[0])
    #        print('gates_sorted_by_qubits:', gates_sorted_by_qubits)
    #        qubits_gates = [g.qubits[0][0] for g in gates_sorted_by_qubits]
    #        for q in self.qubits:
    #            if q not in qubits_gates:
    #                print(f'qubit {q} not in gates')
    #                _unitary = sp.physics.quantum.TensorProduct(_unitary, sp.eye(2)) if _unitary is not None else sp.eye(2)
    #            elif q in qubits_gates:
    #                for gate in gates_sorted_by_qubits:
    #                    print('gate:', gate.name, 'qubits:', gate.qubits)
    #                    if q in gate.qubits[0]:
    #                        print('    gate:', gate.name, 'qubits:', gate.qubits)
    #                        _unitary = sp.physics.quantum.TensorProduct(_unitary, gate.matrix) if _unitary is not None else gate.matrix
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

    def create_numeric_unitary_from_symbolic(self):
        # subs_symbolic_zeros_in_symbolic_unitary based on this function, so if making changes to either you might want to adapt the other as well 
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
            qc_gate_name = _map_qasmgatename_to_qcgatename(qgate_name)
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