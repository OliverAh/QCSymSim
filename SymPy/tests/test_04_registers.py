import pathlib
import sys
import sympy as sp
import numpy as np
sys.path.append(str(pathlib.Path(__file__).parent.parent))

import pytest
import QSymPy as qp

class TestAdd__Bits:
    def test_add_qbit_empty(self):
        qb = qp.QBit()
        assert qb.index is None
        assert qb.index_global is None
    
    def test_add_qbit_index(self):
        qb = qp.QBit(index=9)
        assert qb.index == 9
        assert qb.index_global == None

    def test_add_qbit_index_global(self):
        qb = qp.QBit(index_global=5)
        assert qb.index == None
        assert qb.index_global == 5

    def test_add_cbit_empty(self):
        cb = qp.CBit()
        assert cb.index is None
        assert cb.index_global is None

    def test_add_cbit_index(self):
        cb = qp.CBit(index=9)
        assert cb.index == 9
        assert cb.index_global == None

    def test_add_cbit_index_global(self):
        cb = qp.CBit(index_global=5)
        assert cb.index == None
        assert cb.index_global == 5
    
class TestAdd__Registers:
    def test_add_qreg_size_int(self):
        qr = qp.QReg(name='qreg1', bits=5)
        assert qr.name == 'qreg1'
        assert qr.size == 5
        assert qr.bit_type == 'quantum'

    def test_add_qreg_size_list(self):
        qbits = [qp.QBit(i) for i in range(3)]
        qr = qp.QReg(name='qreg2', bits=qbits)
        assert qr.name == 'qreg2'
        assert qr.size == 3
        assert qr.bit_type == 'quantum'

