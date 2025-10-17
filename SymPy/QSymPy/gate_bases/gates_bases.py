import numpy as np
import itertools
import sympy as sp
import sympy.physics.quantum as sp_phy_qant
from typing import Iterable, Mapping
import copy
import functools

def get_known_gates_classes() -> dict[str,type]:
    '''Returns a dictionary of gates and their respective classes, that can be used for circuits.
    They are determined by looking for subclasses of QuantumGate.
    The returned dict is in the form {'CX': type}.
    '''
    _known_gates_classes = [cls for cls in QuantumGate.__subclasses__()]
    _known_gates_classes.extend(scls for cls in _known_gates_classes for scls in cls.__subclasses__() if scls is not None)
    _known_gates_classes = [cls for cls in _known_gates_classes if getattr(cls, 'name_gate_collection', False)]
    
    _known_gates_classes = {cls.name_gate_collection: cls for cls in _known_gates_classes}
    # print('Known gates:')
    # for k,v in _known_gates_classes.items():
    #     print(' -', k, ':\t', v)
    return _known_gates_classes

def _validate_inputs_to_QuantumGate(self, name: str=None, name_short: str=None, shape:tuple[int,int]=None,
                                    qubits_t: list[int]=None, qubits_c: None|list[int]=None, step: int=None,
                                    num_summands_decomposed: int=None, parameters: None|list[dict[str, float]]|dict[str, float] = None,
                                    gates_t: Mapping[int,Iterable[str]]=None, control_values_c: Mapping[int,Iterable[int]]|None=None
                                    ) -> bool:
    '''Validates the inputs to the QuantumGate class. Also aims to distinguish between (non-)parameterized gates, instead of havin a separate class, e.g. QuantumGateParameterized.'''
    assert name is not None and isinstance(name, str) and name != '', f'name must be a non-empty string, but is {name}.'
    assert name_short is not None and isinstance(name_short, str) and name_short != '', f'name_short must be a non-empty string, but is {name_short}.'
    assert shape is not None and isinstance(shape, tuple|list) and len(shape) == 2, f'Shape must be a tuple/list of length 2, but is {shape}.'
    assert qubits_t is not None and isinstance(qubits_t, tuple|list), f'Qubits_t must be a tuple/list, but is {qubits_t}.'
    assert qubits_c is None or isinstance(qubits_c, tuple|list), f'Qubits_c must be a tuple/list or None, but is {qubits_c}.'
    assert step is not None and isinstance(step, int), f'Step must be an integer, but is {step}.'
    assert num_summands_decomposed is not None or isinstance(num_summands_decomposed, int), f'Num_summands_decomposed must be an integer, but is {num_summands_decomposed}.'
    if parameters is not None:
        assert type(self).is_parametric, f'Parameters are only supported in parameterized gates, but {type(self)} is not parameterized.'
        if isinstance(parameters, tuple|list):
            assert all([isinstance(p, dict) for p in parameters]), f'If parameters is a list/tuple, all its elements must be dicts, but is {parameters}.'
        elif isinstance(parameters, dict):
            pass
        else:
            assert False, f'Parameters must be a dict, list/tuple of dicts, or None, but is {parameters}.'
    if gates_t is not None:
        assert isinstance(gates_t, Mapping), f'Gates_t must be a Mapping, but is {gates_t}.'
        assert all([isinstance(k, int) and isinstance(v, tuple|list) for k,v in gates_t.items()]), f'Gates_t must be a dict of int keys and tuple/list values, but is {gates_t}.'
        assert all([k in qubits_t for k in gates_t.keys()]), f'All keys of gates_t must be in qubits_t, but is gates_t.keys(): {list(gates_t.keys())}, qubits_t: {qubits_t}.'
        assert all([len(v) == num_summands_decomposed for v in gates_t.values()]), f'All values of gates_t must have length num_summands_decomposed={num_summands_decomposed}, but is {gates_t}.'
    return True

class QuantumGate:
    '''Base class for quantum gates. It holds information of the operation that is applied to the qubit(s).
    It does not hold the matrix for the whole quantum circuit, but just for the qubit(s) it's supposed to change the state of.
    Lists of qubits and matrices are assumed to be encoded as most significant qubit first, i.e. left and bottom, e.g. |q2 q1 q0⟩ = |q2> kron |q1> kron |q0>, 
    |0 0 1> = (1,0) kron (1,0) kron (0,1) = (0,1,0,0,0,0,0,0).'''
    # We need information about the qubit and step here, because otherwise the atomic symbols could not be distinguished from those of other gates
    is_should_be_listed_in_gate_collection = False
    name_gate_collection = ''
    num_qubits_t = 0
    num_qubits_c = 0
    is_parametric = False
    parameters=()
    def __init__(self, name: str=None, name_short: str=None, shape:tuple[int,int]=None,
                       qubits_t: list[int]=None, qubits_c: None|list[int]=None, step: int=None,
                       num_summands_decomposed: int=None, parameters: None|list[dict[str, float]] = None):
        assert _validate_inputs_to_QuantumGate(self, name=name, name_short=name_short, shape=shape,
                                               qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                                               num_summands_decomposed=num_summands_decomposed, parameters=parameters)
        self.name = name
        self.name_short = name_short
        self.shape = shape
        iter_indices = [str(sub[0])+str(sub[1]) for sub in itertools.product(range(shape[0]), range(shape[1]))] # (2,2) -> ['00', '01', '10', '11']
        str_qubits_t = ''.join([str(q)+'' for q in qubits_t])
        str_qubits_c = ''.join([str(q)+'' for q in qubits_c]) if qubits_c is not None else ''
        self.atomics = {sub: sp.symbols(self.name+'_qt'+str_qubits_t+'_qc'+str_qubits_c+'_s'+str(step)+'_p'+str(sub[0])+str(sub[1])
                                        , complex=True, finite=True, extended_finite=False) for sub in iter_indices} # {'00': I_qt0_qc0_s0_p00, '01': I_qt0_qc0_s0_p01, ...}
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
        self.treat_numeric_only = False

        self.parameters = parameters
        self.atomics_alt = None
        self.matrix_alt = None
        self.matrix22_t_alt = None

    
    def __iter__(self):
        for key in self.__dict__:
            yield key, getattr(self, key)
    def items(self):
        return self.__iter__()

    def __str__(self):
        s = ''
        s += f'({self.name}, s: {self.step}, qt: {self.qubits_t}, q_c: {self.qubits_c}, p: {self.parameters})'
        return s
    def __repr__(self):
        return self.__str__()


#class QuantumGateParameterized(QuantumGate):
#    '''Base class for parameterized quantum gates. It holds information of the operation that is applied to the qubit(s).
#      It is a subclass of QuantumGate, and uses an additional symbolic matrix.'''
#    is_should_be_listed_in_gate_collection = False
#    def __init__(self, name: str='', name_short: str='', shape:tuple[int,int]=[2,2], qubits_t: list[int]=[0], qubits_c: None|list[int]=None, step: int=0, num_summands_decomposed: int=1,
#                 parameters: None|list[dict[str, float]] = None):
#        super().__init__(name, name_short, shape, qubits_t, qubits_c, step, num_summands_decomposed)
#        self.is_parametric = True
#        self.parameters = parameters
#        self.matrix_alt = self.matrix.copy()
#        self.matrix22_t_alt = copy.deepcopy(self.matrix22_t) # deepcopy just to be sure that _alt and non-alt are not linked

class QuantumGateMultiQubit(QuantumGate):
    '''Base class for multi-qubit quantum gates. It holds information of the operation that is applied to the qubits. It is a subclass of QuantumGate,
      ALSO IF IT USES PARAMETERIZED GATES WITHIN. In the latter case parameters can be specified as an additional kwarg in the order the parameterized gates appear,
      e.g. CUU(..., parameters=[{theta:..., phi:..., lambda:...}, {theta:..., phi:..., lambda:...}])
      or   CUX(..., parameters=[{theta:..., phi:..., lambda:...}, None|[]]) as X is not parameterized
      or  CUXU(..., parameters=[{theta:..., phi:..., lambda:...}, None|[], {theta:..., phi:..., lambda:...}]) as X is not parameterized'''
    is_should_be_listed_in_gate_collection = False
    def __init__(self, name: str=None, name_short: str=None, shape:tuple[int,int]=None, 
                 qubits_t: Iterable[int]=None, qubits_c: Iterable[int]|None=None, step: int=None, 
                 num_summands_decomposed: int=None, parameters: Iterable[Mapping[str, float]]|None=None,
                 gates_t: Mapping[int,Iterable[str]]=None, control_values_c: Mapping[int,Iterable[int]]|None=None):
        _shape = 2**(len(qubits_t) + (len(qubits_c) if qubits_c is not None else 0))
        _shape = (_shape, _shape)
        _num_summands_decomposed = len(gates_t[qubits_t[0]])
        assert len(qubits_t) == len(gates_t), f'Number of target qubits qubits_t ({len(qubits_t)}) must be equal to length of gates_t ({len(gates_t)}).'
        assert all([len(gt) == _num_summands_decomposed for gt in gates_t.values()]), f'All values of gates_t must have length num_summands_decomposed={_num_summands_decomposed}, but is {gates_t}.'
        assert (qubits_c is None and control_values_c is None) or (qubits_c is not None and control_values_c is None) or (len(qubits_c) == len(control_values_c)), f'If control qubits and control values are given, their lengths must be equal, but is len(qubits_c)={len(qubits_c)} and len(control_values_c)={len(control_values_c)}.'
        assert set(qubits_t).isdisjoint(set(qubits_c)) if qubits_c is not None else True, 'Target and control qubits must be different.'
        del _shape, _num_summands_decomposed
        
        super().__init__(name=name, name_short=name_short, shape=shape,
                         qubits_t=qubits_t, qubits_c=qubits_c, step=step,
                         num_summands_decomposed=num_summands_decomposed, parameters=parameters)
        self.gates_t = gates_t
        self.control_values_c = control_values_c
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
        #_sym_gates_t = {qt: [_get_gate_class_from_name(vv)(name=name+'_'+vv, qubits_t=[qt], qubits_c=_qubits_c, step=step) for vv in v] for qt, v in gates_t.items()}
        known_gates_classes = get_known_gates_classes()
        _parameters = parameters if parameters is not None else [None]*len(qubits_t)
        _sym_gates_t = {qt: [     known_gates_classes[vv](name=name+'_'+vv, qubits_t=[qt], qubits_c=_qubits_c, step=step, parameters=_parameters[i]) if known_gates_classes[vv].is_parametric 
                             else known_gates_classes[vv](name=name+'_'+vv, qubits_t=[qt], qubits_c=_qubits_c, step=step)
                             for vv in v] for i, (qt, v) in enumerate(gates_t.items())}
        #self.matrix22_t         = {k: _get_gate_class_from_name(v)(name=name+'_'+v, qubits_t=[k], qubits_c=[], step=step) for k,v in gates_t.values()}
        self.matrix22_t         = {k: [vv.matrix         for vv in v] for k,v in _sym_gates_t.items()}
        self.matrix22_t_numeric = {k: [vv.matrix_numeric for vv in v]for k,v in _sym_gates_t.items()}
        
        #_control_values_c = None
        #if qubits_c is None or len(qubits_c)==0 and control_values_c is None:
        #    pass
        #elif qubits_c is not None and control_values_c is None:
        #    #_control_values_c = {qc: [1]*self.num_summands_decomposed for qc in qubits_c} # if no control values are given, assume all 1s
        #    raise ValueError('If control qubits are given, control values MUST be given as well.')
        #elif qubits_c is not None and control_values_c is not None:
        #    _control_values_c = control_values_c
        #else:
        #    raise ValueError('If control values are given, control qubits MUST be given as well.')

        if qubits_c is None:
            _merged_matrix22         = self.matrix22_t
            _merged_matrix22_numeric = self.matrix22_t_numeric
        else:
            if control_values_c is None:
                control_values_c = {qc: [1]*self.num_summands_decomposed for qc in qubits_c} # if no control values are given, assume all 1s
            self.matrix22_c         = {qc: [_dict_sym_densmats[qcvv] for qcvv in qcv] for qc, qcv in control_values_c.items()}
            self.matrix22_c_numeric = {qc: [_dict_num_densmats[qcvv] for qcvv in qcv] for qc, qcv in control_values_c.items()} 

            _merged_matrix22         = self.matrix22_t | self.matrix22_c
            _merged_matrix22_numeric = self.matrix22_t_numeric | self.matrix22_c_numeric

            
        self.matrix_decomposed = sp.Add(*[sp_phy_qant.TensorProduct(*[_merged_matrix22[k][i]         for k in sorted(_merged_matrix22.keys(),         reverse=True)]) for i in range(self.num_summands_decomposed)])
        self.matrix_numeric    =    sum( [ functools.reduce(np.kron, [_merged_matrix22_numeric[k][i] for k in sorted(_merged_matrix22_numeric.keys(), reverse=True)]) for i in range(self.num_summands_decomposed)])
        #self.matrix_numeric    =    sum( [                  np.kron(*[_merged_matrix22_numeric[k][i] for k in sorted(_merged_matrix22_numeric.keys(), reverse=True)]) for i in range(self.num_summands_decomposed)])
        
        if parameters is not None:
            self.is_parametric = True
            self.parameters = _parameters
            self.atomics_alt = {k: [getattr(vv, 'atomics_alt', None) for vv in v] for k,v in _sym_gates_t.items()}
            self.atomics_alt = {key: [e for e in val if e is not None] for key, val in self.atomics_alt.items()} # keep only non-None entries
            self.atomics_alt = {key: [vval[key] for val in self.atomics_alt.values() for vval in val if key in vval] for pl in self.parameters for key in pl.keys()}
            self.matrix_alt = self.matrix.copy()
            self.matrix22_t_alt = {k: [getattr(vv, 'matrix_alt', None) for vv in v] for k,v in _sym_gates_t.items()}

