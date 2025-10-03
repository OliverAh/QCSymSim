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
import tqdm

#from . import gate_defs
from . import gate_bases
from .gate_bases import *
from .gate_defs import *
from .bits_regs import QBit, CBit, QReg, CReg, MapBitID


class GateCollection(QuantumGate):
    '''Class to hold lists of gates of the same type within a single quantum circuit. It is just a helper to make it more easy to traverse the circuit for evaluation.'''
    def __init__(self):
        #_known_gates_classes = [cls for cls in QuantumGate.__subclasses__()]
        #_known_gates_classes.extend(scls for cls in _known_gates_classes for scls in cls.__subclasses__() if scls is not None)
        #_known_gates_classes = [cls for cls in _known_gates_classes if getattr(cls, 'is_should_be_listed_in_gate_collection', False)]
        #self.known_gates = list(dict.fromkeys([cls.name_gate_collection for cls in _known_gates_classes])) # remove possible duplicates and keep order, without order we could do list(set(...))
        self.known_gates = gate_bases.get_known_gates_classes()
        self.collections = {id: [] for id, val in self.known_gates.items() if val.is_should_be_listed_in_gate_collection}

    def __str__(self):
        s = 'GateCollection:\n'
        for k,v in self.collections.items():
            s += f'  {k}: {len(v)} gates\n'
        return s



class QuantumCircuit():
    def __init__(self, num_qubits: int=0, num_clbits: int=0):
        
        self.context = MapBitID()
        self.qubits = list(reversed(list(range(num_qubits))))
        self.clbits = list(reversed(list(range(num_clbits))))
        self.gate_collection = GateCollection()
        self.barrier_collection = []
        self.steps = {} # will contain {step_number: [gate1, gate2, ...]}
        self.unitary = None
        self.allowed_gates = self.gate_collection.known_gates

    def add_qubit(self):
        q = QBit(context=self.context)
        self.qubits = list(reversed(list(range(len(self.qubits)+1))))
        return q
    
    def add_cbit(self):
        c = CBit(context=self.context)
        self.clbits = list(reversed(list(range(len(self.clbits)+1))))
        return c
    
    def add_qreg(self, name:str='', size:int=99):
        if self.context.has_qreg(name):
            raise ValueError('Quantum register with name already exists:', name)
        qreg = QReg(context=self.context, name=name, size=size)
        self.qubits = list(reversed(list(range(len(self.qubits)+size))))
        return qreg

    def add_creg(self, name:str='', size:int=99):
        if self.context.has_creg(name):
            raise ValueError('Classical register with name already exists:', name)
        creg = CReg(context=self.context, name=name, size=size)
        self.clbits = list(reversed(list(range(len(self.clbits)+size))))
        return creg

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
    
    def assemble_symbolic_unitary(self, use_alternative_repr:bool=False, replace_symbolic_zeros_and_ones:bool=True):
        '''Assemble the unitary matrix of the quantum circuit.'''
        self.unitary = sp.eye(2**len(self.qubits))
        for step in tqdm.tqdm(self.steps.keys()):# step is step number of timesteps
            unitary_step = sp.zeros(2**len(self.qubits))
            gates_step = self.steps[step]
            #print(gates_step)
            num_sum_max = max([g.num_summands_decomposed for g in gates_step])
            for i in range(num_sum_max):
                gates = [None]*len(self.qubits)
                for gate in gates_step:
                    if gate.num_summands_decomposed > i:
                        #q_t = gate.qubits_t[0]
                        #gates[q_t] = gate.matrix22_t[q_t][i]
                        for q_t in gate.qubits_t:
                            gates[q_t] = gate.matrix22_t[q_t][i]
                            if replace_symbolic_zeros_and_ones:
                                gates[q_t] = gates[q_t].subs({s: n for s,n in zip(gate.matrix22_t[q_t][i].flat(), 
                                                                                  gate.matrix22_t_numeric[q_t][i].flatten()) if n == 0 or n == 1})
                        #print(gates[q_t])
                        if gate.qubits_c is not None:
                            #q_c = gate.qubits_c[0]
                            #gates[q_c] = gate.matrix22_c[q_c][i]
                            for q_c in gate.qubits_c:
                                gates[q_c] = gate.matrix22_c[q_c][i]
                                if replace_symbolic_zeros_and_ones:
                                    gates[q_c] = gates[q_c].subs({s: n for s,n in zip(gate.matrix22_c[q_c][i].flat(), 
                                                                                      gate.matrix22_c_numeric[q_c][i].flatten()) if n == 0 or n == 1})
                for j in range(len(gates)):
                    if gates[j] is None:
                        #gates[j] = Identity_Gate(qubits_t=[j], qubits_c=None, step=step).matrix22_t[j][0]
                        gates[j] = sp.eye(2)
                gates = list(reversed(gates)) # reverse the list to match the order of the qubits (little endian) for the tensor product
                assert all(isinstance(g, sp.MatrixBase) for g in gates)
                assert all(all(isinstance(e, sp.Expr) for e in g) for g in gates)
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
                print(gate)
                for i in range(gate.num_summands_decomposed):
                    symbs_t = [e for q_t in gate.qubits_t for e in gate.matrix22_t[q_t][i].flat()]
                    nums_t = [e for q_t in gate.qubits_t for e in gate.matrix22_t_numeric[q_t][i].flatten()]
                    print(symbs_t)
                    print(nums_t)
                    for s, n in zip(symbs_t, nums_t):
                        print(s, n)
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
    
def _parse_binary_expression_to_single_object(be:openqasm3.ast.BinaryExpression) -> object:
    _op  = be.op
    _lhs = be.lhs
    _rhs = be.rhs

    return_objs = []
    for obj in (be.lhs, be.rhs):
        if isinstance(obj, openqasm3.ast.BinaryExpression):
            obj = _parse_binary_expression_to_single_object(obj)
        elif isinstance(obj, openqasm3.ast.UnaryExpression):
            obj = _parse_unary_expression_to_single_object(obj)
        elif isinstance(obj, openqasm3.ast.Identifier):
            if obj.name == 'pi':
                obj = sp.pi
            else:
                raise ValueError('Unknown identifier:', obj)
        elif isinstance(obj, openqasm3.ast.IntegerLiteral):
            obj = int(obj.value)
        else:
            raise ValueError('Unknown object:', obj)
        return_objs.append(obj)
    _lhs, _rhs = return_objs
    assert isinstance(_op, openqasm3.ast.BinaryOperator)
    if _op.name == '/':
        return _lhs/_rhs
    elif _op.name == '*':
        return _lhs*_rhs
    elif _op.name == '+':
        return _lhs+_rhs
    elif _op.name == '-':
        return _lhs-_rhs
    elif _op.name == '**':
        return _lhs**_rhs
    else:
        raise ValueError('Unknown binary operator:', _op)

def _parse_unary_expression_to_single_object(ue:openqasm3.ast.UnaryExpression) -> object:
    _op = ue.op
    _expr = ue.expression

    if isinstance(_expr, openqasm3.ast.BinaryExpression):
        _expr = _parse_binary_expression_to_single_object(_expr)
    elif isinstance(_expr, openqasm3.ast.UnaryExpression):
        _expr = _parse_unary_expression_to_single_object(_expr)
    elif isinstance(_expr, openqasm3.ast.Identifier):
        if _expr.name == 'pi':
            _expr = sp.pi
        else:
            raise ValueError('Unknown identifier:', _expr)
    elif isinstance(_expr, openqasm3.ast.IntegerLiteral):
        _expr = int(_expr.value)
    elif isinstance(_expr, openqasm3.ast.FloatLiteral):
        _expr = float(_expr.value)
    else:
        raise ValueError('Unknown object:', _expr)
    
    assert isinstance(_op, openqasm3.ast.UnaryOperator)
    if _op.name == '+':
        return +_expr
    elif _op.name == '-':
        return -_expr
    else:
        raise ValueError('Unknown unary operator:', _op)

def _iterate_over_qasm_statements(qpf:openqasm3.ast.Program) -> QuantumCircuit:
    qc = QuantumCircuit()
    _gate_set_step = set()
    for s in qpf.statements:
        if isinstance(s, openqasm3.ast.Include):
            print('Included filenames:', s.filename)
            if s.filename != 'stdgates.inc':
                raise ValueError('Only stdgates.inc is supported as include file, not:', s.filename)
        elif isinstance(s, openqasm3.ast.ClassicalDeclaration):
            if isinstance(s.type, openqasm3.ast.BitType):
                if not qc.context.has_creg(s.identifier.name):
                    qc.add_creg(name=s.identifier.name, size=s.type.size.value)
                else:
                    raise ValueError('Classical register with name already exists:', s.identifier.name)
            else:
                raise ValueError('Unknown type', 's.type:', s.type)
        elif isinstance(s, openqasm3.ast.QubitDeclaration):
            # in contrast to ClassicalDeclaration, QubitDeclaration does not have a type field,
            # the "identifier" is always qubit
            print(s.qubit.name)
            if not qc.context.has_qreg(s.qubit.name):
                qc.add_qreg(name=s.qubit.name, size=s.size.value)
        elif isinstance(s, openqasm3.ast.QuantumGate):
            # s.name is of type openqasm3.ast.Identifier, s.name.name is the actual name string
            _gate_name = s.name.name
            _gate_modifiers = s.modifiers
            assert _gate_modifiers==[], f'gate modiefier: {_gate_modifiers} ... Gate modifiers not yet implemented, you should implement separate gate instead.'
            _gate_arguments = s.arguments
            if _gate_arguments == []:
                pass
            elif isinstance(_gate_arguments, list):
                _list_gate_arguments = []
                for arg in _gate_arguments:
                    if isinstance(arg, openqasm3.ast.BinaryExpression):
                        _list_gate_arguments.append(_parse_binary_expression_to_single_object(arg))
                    elif isinstance(arg, openqasm3.ast.UnaryExpression):
                        _list_gate_arguments.append(_parse_unary_expression_to_single_object(arg))
                    elif isinstance(arg, openqasm3.ast.Identifier):
                        if arg.name == 'pi':
                            _list_gate_arguments.append(sp.pi)
                        else:
                            raise ValueError('Unknown identifier:', arg)
                    elif isinstance(arg, openqasm3.ast.IntegerLiteral):
                        _list_gate_arguments.append(int(arg.value))
                    elif isinstance(arg, openqasm3.ast.FloatLiteral):
                        _list_gate_arguments.append(float(arg.value))
                    else:
                        raise ValueError('Unknown gate argument:', arg)
                    print('  Gate arguments:', _list_gate_arguments)
            else:
                raise ValueError('Unknown gate arguments:', _gate_arguments)
            
            _gate_qubits = s.qubits
            assert isinstance(_gate_qubits, list), f'expected list for gate qubits: {_gate_qubits}'
            _gate_qubits_glob_ids = []
            for _q in _gate_qubits:
                # indexed identifier means that object requires indexing
                assert isinstance(_q, openqasm3.ast.IndexedIdentifier)
                assert isinstance(_q.name, openqasm3.ast.Identifier)
                if not qc.context.has_qreg(_q.name.name):
                    raise ValueError('Qubit register with name does not exist:', _q.name.name, '... tried to apply gate', s)
                else:
                    assert isinstance(_q.indices, list)
                    assert len(_q.indices)==1, f'broadcasting of gates on registers is not implemented. {s}'
                    assert isinstance(_q.indices[0], list)
                    assert len(_q.indices[0])==1, f'broadcasting of gates on registers is not implemented. {s}'
                    assert isinstance(_q.indices[0][0], openqasm3.ast.IntegerLiteral)
                    _q_name = _q.name.name
                    _q_index = _q.indices[0][0].value
                    _gate_qubits_glob_ids.append(qc.context.map_bit_id['q']['l2g'][_q_name][_q_index])
            _gate_target = OPENQASM3_TO_QSYMPY[_gate_name]
            _n_q_t = _gate_target.num_qubits_t
            _n_q_c = _gate_target.num_qubits_c
            
            _step = list(qc.steps.keys())[-1] if len(qc.steps)>0 else 0
            print('Current step for gates:', _step)
            if set() == _gate_set_step.intersection(_gate_qubits_glob_ids):
                _gate_set_step.update(_gate_qubits_glob_ids)
            else:
                _step += 1
                _gate_set_step = set(_gate_qubits_glob_ids)
            print('Adding gate at step:', _step)
            print('  Gate qubits global IDs:', _gate_qubits_glob_ids)

            _parameters = {k:v for (k,v) in zip(_gate_target.parameters, _list_gate_arguments)} if _gate_arguments != [] else None
            qc.add_gate(name=_gate_target.name_gate_collection, qubits_t=_gate_qubits_glob_ids[:_n_q_t], qubits_c=_gate_qubits_glob_ids[_n_q_t:_n_q_t+_n_q_c], step=_step,
                         parameters=_parameters)
            print('Quantum gate name:', _gate_name)
            print('Quantum gate arguments:', _gate_arguments)
            print('Quantum gate qubits names + ids:', [(_q.name.name, _q.indices[0][0].value) for _q in _gate_qubits])
            print('QuantumGate to add:', OPENQASM3_TO_QSYMPY[_gate_name])
        elif isinstance(s, openqasm3.ast.QuantumBarrier):
            _step = list(qc.steps.keys())[-1]+1 if len(qc.steps)>0 else 1
            qc.add_barrier(step=_step)
            print('Barrier added at step:', _step)
        elif isinstance(s, openqasm3.ast.QuantumMeasurementStatement):
            assert isinstance(s.measure, openqasm3.ast.QuantumMeasurement)
            assert isinstance(s.measure.qubit, openqasm3.ast.IndexedIdentifier)
            assert isinstance(s.measure.qubit.name, openqasm3.ast.Identifier)
            assert isinstance(s.measure.qubit.indices, list)
            assert len(s.measure.qubit.indices)==1
            assert isinstance(s.measure.qubit.indices[0], list)
            assert len(s.measure.qubit.indices[0])==1
            _qubit_name = s.measure.qubit.name.name
            _qubit_index = s.measure.qubit.indices[0][0].value
            assert isinstance(s.target, openqasm3.ast.IndexedIdentifier)
            assert isinstance(s.target.name, openqasm3.ast.Identifier)
            assert isinstance(s.target.indices, list)
            assert len(s.target.indices)==1
            assert isinstance(s.target.indices[0], list)
            assert len(s.target.indices[0])==1
            _cbit_name = s.target.name.name
            _cbit_index = s.target.indices[0][0].value
            #print('  Measured qubit name:', _qubit_name)
            #print('  Measured qubit index:', _qubit_index)
            #print('  Target classical bit name:', _cbit_name)
            #print('  Target classical bit index:', _cbit_index)
            print('WARNING:  Measurement not yet implemented.')
            
        else:
            raise ValueError('Unknown statement:', s)
    return qc

def openqasm3_to_qc(filepath: pathlib.Path) -> QuantumCircuit:
    '''Function to read an OpenQASM 3 file and convert it to a QuantumCircuit object.
    '''

    
    with open(filepath, 'r') as f:
        qpf = f.read()
    
    qpf = openqasm3.parse(qpf)
    
    qc = _iterate_over_qasm_statements(qpf)
    return qc

