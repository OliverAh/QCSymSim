import sympy as sp
import itertools

class Bit():
    def __init__(self, index:None|int=None, index_global:None|int=None):
        self.index = index
        self.index_global = index_global
        self.is_quantum = None
        self.is_classical = None

    def __repr__(self):
        return None

class CBit(Bit):
    def __init__(self, index:None|int=None, index_global:None|int=None):
        super().__init__(index, index_global=index_global)
        self.is_quantum = False
        self.is_classical = not self.is_quantum # only for convenience

    def __repr__(self):
        return f'CBit[{self.index}]'

class QBit(Bit):
    def __init__(self, index:None|int=None, index_global:None|int=None):
        super().__init__(index, index_global=index_global)
        self.is_quantum = True
        self.is_classical = not self.is_quantum # only for convenience
    
    def __repr__(self):
        return f'QBit[{self.index}]'

class Register():
    def __init__(self, name, bit_type, bits):
        self.name = name
        self.size = len(bits)
        self.bit_type = bit_type
        self.bits = bits
    
    def __getitem__(self, index):
        return self.bits[index]
    
    def __len__(self):
        return self.size
    
    def __iter__(self):
        return iter(self.bits)
    
    def __repr__(self):
        return f'Register[{self.name} {self.bit_type} {self.bits}]'
    
class QReg(Register):
    def __init__(self, name:str='', bits:int|list[QBit]=99):
        if isinstance(bits, int):
            bits = [QBit(i) for i in range(bits)]
        super().__init__(name, bit_type='quantum', bits=bits)

class CReg(Register):
    def __init__(self, name:str='', bits:int|list[CBit]=99):
        if isinstance(bits, int):
            bits = [CBit(i) for i in range(bits)]
        super().__init__(name, bit_type='classical', bits=bits)





