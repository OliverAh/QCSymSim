import sympy as sp
import itertools
from typing import Iterable, List, Tuple, Union


###
# mappings between global and local indices. reg_name might also be bit_name
# global: list of tuples ordered by global id (str:reg_name, int:local_index) local_index=-1 if not in a register
# local: dict of ints {str:reg_name -> int:global_index}
#map_q_bit_id = {'g2l':[], 'l2g':{'default':[]}}
#map_c_bit_id = {'g2l':[], 'l2g':{'default':[]}}

class MapBitID():
    '''Class that holds the mapping between global and local indices of
    bits (qubits and classical bits). It is used as a statefull context object
    for circuits, because one might define multiple circuits within a single runtime,
    where each circuit has its own global indexing of qbits/cbits, i.e. 0...n-1, where n
    is the number of qbits/cbits in the circuit. 
       
    The class holds a single dict - map_bit_id - with the following structure:  
    {'q': {'g2l': [(reg_name, local_index)], 'l2g': {'reg_name': [global_index]}},  
     'c': {'g2l': [(reg_name, local_index)], 'l2g': {'reg_name': [global_index]}}}  
    where 'q' and 'c' stand for quantum and classical bits, respectively.  
    'g2l' (global to local) maps the global index to a tuple (reg_name, local_index), where reg_name is the name of the 
    register the bit belongs to, or 'default' if it is a single bit not explicitly assigned to any register.  
    'l2g' (local to global) maps the register name to a list of global indices of bits that belong to that register.
    '''
    def __init__(self):
        self.map_bit_id = {'q': {'g2l':[], 'l2g':{'default':[]}},
                   'c': {'g2l':[], 'l2g':{'default':[]}}}
    
    def has_qreg(self, name:str):
        return name in self.map_bit_id['q']['l2g']
    def has_creg(self, name:str):
        return name in self.map_bit_id['c']['l2g']
    def num_qbits(self):
        return len(self.map_bit_id['q']['g2l'])
    def num_cbits(self):
        return len(self.map_bit_id['c']['g2l'])

class Bit():
    '''Base class for qubits and classical bits.'''
    def __init__(self, context:MapBitID=None, index_local:int=-1, is_quantum:bool=False, is_classical:bool=False, name:str='default'):
        assert context is not None, "A Bit must be created within a context (MapBitID)"
        assert (is_quantum != is_classical), "A Bit must be either quantum or classical"
        self.name = name
        self.index_local = None
        self.index_global = None
        self.is_quantum = is_quantum
        self.is_classical = is_classical # only for convenience

        if self.is_quantum:
            t = 'q'
        else:
            t = 'c'
        if index_local == -1:
            assert name == 'default', "If no local index is given, the name must be 'default'"
            self.index_local = len(context.map_bit_id[t]['l2g']['default'])
        else:
            self.index_local = index_local
        self.index_global = len(context.map_bit_id[t]['g2l'])
        
        context.map_bit_id[t]['g2l'].append([self.name, self.index_local])
        if self.name not in context.map_bit_id[t]['l2g']:
            context.map_bit_id[t]['l2g'][self.name] = []
        context.map_bit_id[t]['l2g'][self.name].append(self.index_global)



    def __repr__(self):
        return f'{type(self).__name__}[{self.index_local}, {self.index_global}]'

class CBit(Bit):
    '''Classical Bit (CBit). Does only exist for convenience. One should preferably
    use registers, even when only a single classical bit is needed, which would be a register
    of size 1.  
    Therefore, a name can not be provided for a CBit, as it will always be assigned
    to a register named "default".
    '''
    def __init__(self, context:MapBitID=None):
        super().__init__(context=context, is_quantum=False, is_classical=True)

class QBit(Bit):
    '''Quantum Bit (Qubit). Does only exist for convenience. One should preferably
    use registers, even when only a single qubit is needed, which would be a register
    of size 1.  
    Therefore, a name can not be provided for a QBit, as it will always be assigned
    to a register named "default".
    '''
    def __init__(self, context:MapBitID=None):
        super().__init__(context=context, is_quantum=True, is_classical=False)

class BitRegister():
    ''' Base class for registers of bits. A register is either quantum or classical, but not both.'''
    def __init__(self, context:MapBitID=None, name:str=None, bit_type:str=None, size:int=None):
        assert isinstance(name, str)
        assert bit_type in ['quantum', 'classical', 'q', 'c']
        assert isinstance(size, int) and size > 0
        assert name not in ['', 'default'], "Register name cannot be '' or 'default'"
        self.name = name
        self.size = size
        if bit_type in ['q', 'quantum']:
            self.bit_type = 'quantum'
            self.is_quantum = True
            self.is_classical = False
        elif bit_type in ['c', 'classical']:
            self.bit_type = 'classical'
            self.is_quantum = False
            self.is_classical = True
        self.bits = [Bit(context=context, index_local=i, is_quantum=self.is_quantum, is_classical=self.is_classical, name=name) for i in range(size)]

    def __getitem__(self, index):
        return self.bits[index]
    
    def __len__(self):
        return self.size
    
    def __iter__(self):
        return iter(self.bits)
    
    def __repr__(self):
        return f'Register[{self.name} {self.bit_type} {self.bits}]'

class QReg(BitRegister):
    '''Quantum Register (QReg). A register that holds quantum bits (qubits).  
    Also use this class for single qubits, i.e. registers of size 1.
    name="default" is reserved for single QBits, instantiated explicitly.'''
    def __init__(self, context:MapBitID=None, name:str='', size:int=99):
        super().__init__(context=context, name=name, bit_type='quantum', size=size)

class CReg(BitRegister):
    '''Classical Register (CReg). A register that holds classical bits (cbits).  
    Also use this class for single cbits, i.e. registers of size 1.
    name="default" is reserved for single CBits, instantiated explicitly.'''
    def __init__(self, context:MapBitID=None, name:str='', size:int=99):
        super().__init__(context=context, name=name, bit_type='classical', size=size)



    
