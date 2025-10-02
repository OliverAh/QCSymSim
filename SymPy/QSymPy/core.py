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

#from . import gate_defs
from . import gate_bases
from .gate_bases import *
from .gate_defs import *
from .src_qubit_qreg import QBit, CBit, QReg, CReg

#KNOWN_GATES_CLASSES = {} # will be set at the end of the file (end of the import of QSymPy), when all gate classes are defined. Uses QuantumGate.__subclasses__()



class GateCollection(QuantumGate):
    '''Class to hold lists of gates of the same type within a single quantum circuit. It is just a helper to make it more easy to traverse the circuit for evaluation.'''
    def __init__(self):
        #_known_gates_classes = [cls for cls in QuantumGate.__subclasses__()]
        #_known_gates_classes.extend(scls for cls in _known_gates_classes for scls in cls.__subclasses__() if scls is not None)
        #_known_gates_classes = [cls for cls in _known_gates_classes if getattr(cls, 'is_should_be_listed_in_gate_collection', False)]
        #self.known_gates = list(dict.fromkeys([cls.name_gate_collection for cls in _known_gates_classes])) # remove possible duplicates and keep order, without order we could do list(set(...))
        self.known_gates = gate_bases.get_known_gates_classes()
        self.collections = {id: [] for id, val in self.known_gates.items() if val.is_should_be_listed_in_gate_collection}

class QuantumCircuit():
    def __init__(self, num_qubits: int=1, num_clbits: int=1):
        self.qubits = list(reversed(list(range(num_qubits))))
        self.clbits = list(reversed(list(range(num_clbits))))
        self.gate_collection = GateCollection()
        self.barrier_collection = []
        self.steps = {} # will contain {step_number: [gate1, gate2, ...]}
        self.unitary = None
        self.allowed_gates = self.gate_collection.known_gates


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
                elif name == 'H_eI':
                    gate = Hadamard_error_I_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'S':
                    gate = S_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'Sdg':
                    gate = Sdg_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'T':
                    gate = T_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'Tdg':
                    gate = Tdg_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'SX':
                    gate = SX_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'U':
                    gate = U_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step, parameters=parameters)
                elif name == 'GP':
                    gate = GP_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step, parameters=parameters)
                elif name == 'RX':
                    gate = RX_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step, parameters=parameters)
                elif name == 'RY':
                    gate = RY_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step, parameters=parameters)
                elif name == 'RZ':
                    gate = RZ_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step, parameters=parameters)
                elif name == 'P':
                    gate = P_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step, parameters=parameters)
                elif name == 'CX':
                    gate = CX_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'CCX':
                    gate = CCX_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'CCXX':
                    gate = CCXX_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'CU':
                    gate = CU_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step, parameters=parameters)
                elif name == 'CUU':
                    gate = CUU_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step, parameters=parameters)
                elif name == 'CY':
                    gate = CY_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'CZ':
                    gate = CZ_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'CP':
                    gate = CP_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step, parameters=parameters)
                elif name =='CRX':
                    gate = CRX_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step, parameters=parameters)
                elif name =='CRY':
                    gate = CRY_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step, parameters=parameters)
                elif name =='CRZ':
                    gate = CRZ_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step, parameters=parameters)
                elif name == 'CH':
                    gate = CH_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'SWAP':
                    gate = SWAP_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                elif name == 'CSWAP':
                    gate = CSWAP_Gate(qubits_t=qubits_t, qubits_c=qubits_c, step=step)
                else:
                    raise ValueError(f'No implementation registered for gate {name}.')
                                
                self.gate_collection.collections[name].append(gate)
                
                if step not in self.steps:
                    self.steps[step] = []
                for gc in self.steps[step]:
                    if any(q in gc.qubits_t  or q in gc.qubits_c for q in qubits_t if gc.qubits_c is not None):
                        raise ValueError(f'Conflict: Target qubit(s) {set(gc.qubits_t).intersection(set(qubits_t))} already have a gate applied at step {step}.')
                    elif qubits_c is not None and gc.qubits_c is not None and any(q in gc.qubits_t or q in gc.qubits_c for q in qubits_c):
                        raise ValueError(f'Conflict: Control qubit(s) {set(gc.qubits_c).intersection(set(qubits_c))} already have a gate applied at step {step}.')
                self.steps[step].append(gate)
                #self.steps.sort()
                #print(self.steps)

    def add_barrier(self, step: None|int=None):
        if step is None:
            step = max(tuple(self.steps.keys()), default=0)
        bar = Barrier(step=step)
        self.barrier_collection.append(bar)
    
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
    
    def subs_symbolic_zerosones_in_symbolic_unitary(self, zeros: bool=True, ones: bool=True):
        '''Substitute entries in the symbolic unitary that are zero with zeros.'''
        # this function based create_numeric_unitary_from_symbolic, so if making changes to either you might want to adapt the other as well 
        for gate_type, gates in self.gate_collection.collections.items():
            for gate in gates:
                for i in range(gate.num_summands_decomposed):
                    for qt in gate.qubits_t:
                        symbs_t = gate.matrix22_t[qt][i].flat()
                        nums_t = gate.matrix22_t_numeric[qt][i].flatten()
                        for s, n in zip(symbs_t, nums_t):
                            if (n == 0 and zeros) or (n == 1 and ones):
                                self.unitary = self.unitary.subs(s, n)
                    
                    if gate.qubits_c is not None:
                        for qc in gate.qubits_c:
                            symbs_c = gate.matrix22_c[qc][i].flat()
                            nums_c = gate.matrix22_c_numeric[gate.qubits_c[0]][i].flatten()
                            for s, n in zip(symbs_c, nums_c):
                                if (n == 0 and zeros) or (n == 1 and ones):
                                    self.unitary = self.unitary.subs(s, n)
    
    @staticmethod
    def _filter_gates_to_replace_with_alternatives(gate_collection, gate_identifications):
        generator_all_gates = (gate for collection in gate_collection.collections.values() for gate in collection)
        if gate_identifications is None:
            list_gates_to_replace = [gate for gate in generator_all_gates]
        else:
            generator_all_identifiers = ((gate_identifications['steps'][i], gate_identifications['names'][i], gate_identifications['qubits_t'][i]) for i in range(len(gate_identifications['steps'])))
            list_gates_to_replace = [it[0] for it in itertools.product(generator_all_gates, generator_all_identifiers) if (it[0].step == it[1][0] and it[0].name == it[1][1] and it[1][2] in it[0].qubits_t)]
        #return [gate for gate in list_gates_to_replace if hasattr(gate, 'matrix_alt')]
        return [gate for gate in list_gates_to_replace if gate.is_parametric]

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
        #self.unitary = self.unitary.subs({v[0]: v[1] for gate in list_gates_to_replace for v in zip(gate.matrix.flat(), gate.matrix_alt.flat())})
        dict_for_subs = {}
        for gate in list_gates_to_replace:
            if isinstance(gate, QuantumGateMultiQubit):
                for qt, mats in gate.gates_t.items():
                    for i,subgatestring in enumerate(mats):
                        if gate.matrix22_t_alt[qt][i] is not None:
                            gate_atomics = [v for gm22 in gate.matrix22_t[qt] for v in gm22 if subgatestring+'_qt'+str(qt) in v.name] # should be 4 elems which end in _p11, _p12, _p21, _p22
                            _di_tmp = {s: n for s, n in zip(gate_atomics, gate.matrix22_t_alt[qt][i].flat())}
                            _l = list(gate.atomics.values())
                            dict_for_subs.update(_di_tmp)
            else:
                dict_for_subs.update({s: n for s, n in zip(gate.matrix.flat(), gate.matrix_alt.flat())})
        
        self.unitary = self.unitary.subs(dict_for_subs)

    def create_numeric_unitary_from_symbolic(self):
        # subs_symbolic_zerosones_in_symbolic_unitary based on this function, so if making changes to either you might want to adapt the other as well 
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

