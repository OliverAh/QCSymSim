import pathlib
import sys
import sympy as sp
import numpy as np

sys.path.append(str(pathlib.Path(__file__).parent.parent))

import pytest
import QSymPy as qp

def test_map_bit_id_instantiation():
    mbi = qp.MapBitID()
    assert mbi.map_bit_id == {'q': {'g2l':[], 'l2g':{'default':[]}},
                   'c': {'g2l':[], 'l2g':{'default':[]}}}
class TestAdd__Bits:
    def test_add_qbit_empty(self):
        mbi = qp.MapBitID()
        qb = qp.QBit(mbi)
        assert qb.index_local == 0
        assert qb.index_global == 0
        assert mbi.map_bit_id['q']['g2l'][0] == ['default', 0]
        assert mbi.map_bit_id['q']['l2g']['default'][0] == 0

    def test_add_cbit_empty(self):
        mbi = qp.MapBitID()
        cb = qp.CBit(mbi)
        assert cb.index_local == 0
        assert cb.index_global == 0
        assert mbi.map_bit_id['c']['g2l'][0] == ['default', 0]
        assert mbi.map_bit_id['c']['l2g']['default'][0] == 0

    def test_add_2qbits_empty(self):
        mbi = qp.MapBitID()
        qb1 = qp.QBit(mbi)
        qb2 = qp.QBit(mbi)
        assert qb1.index_local == 0
        assert qb1.index_global == 0
        assert mbi.map_bit_id['q']['g2l'][0] == ['default', 0]
        assert mbi.map_bit_id['q']['l2g']['default'][0] == 0
        assert qb2.index_local == 1
        assert qb2.index_global == 1
        assert mbi.map_bit_id['q']['g2l'][1] == ['default', 1]
        assert mbi.map_bit_id['q']['l2g']['default'][1] == 1

class TestAdd__Registers:
    def test_add_qreg(self):
        mbi = qp.MapBitID()
        qr = qp.QReg(context=mbi, name='qr1', size=3)
        assert len(qr) == 3
        assert qr.size == 3
        assert qr[0].index_local == 0
        assert qr[0].index_global == 0
        assert qr[1].index_local == 1
        assert qr[1].index_global == 1
        assert qr[2].index_local == 2
        assert qr[2].index_global == 2
        assert mbi.map_bit_id['q']['g2l'][0] == ['qr1', 0]
        assert mbi.map_bit_id['q']['g2l'][1] == ['qr1', 1]
        assert mbi.map_bit_id['q']['g2l'][2] == ['qr1', 2]
        assert mbi.map_bit_id['q']['l2g']['qr1'][0] == 0
        assert mbi.map_bit_id['q']['l2g']['qr1'][1] == 1
        assert mbi.map_bit_id['q']['l2g']['qr1'][2] == 2

    def test_add_creg(self):
        mbi = qp.MapBitID()
        cr = qp.CReg(context=mbi, name='cr1', size=3)
        assert len(cr) == 3
        assert cr.size == 3
        assert cr[0].index_local == 0
        assert cr[0].index_global == 0
        assert cr[1].index_local == 1
        assert cr[1].index_global == 1
        assert cr[2].index_local == 2
        assert cr[2].index_global == 2
        assert mbi.map_bit_id['c']['g2l'][0] == ['cr1', 0]
        assert mbi.map_bit_id['c']['g2l'][1] == ['cr1', 1]
        assert mbi.map_bit_id['c']['g2l'][2] == ['cr1', 2]
        assert mbi.map_bit_id['c']['l2g']['cr1'][0] == 0
        assert mbi.map_bit_id['c']['l2g']['cr1'][1] == 1
        assert mbi.map_bit_id['c']['l2g']['cr1'][2] == 2
        assert len(mbi.map_bit_id['q']['g2l']) == 0

    def test_add_2_qregs(self):
        mbi = qp.MapBitID()
        qr1 = qp.QReg(context=mbi, name='qr1', size=2)
        qr2 = qp.QReg(context=mbi, name='qr2', size=2)
        assert len(qr1) == 2
        assert len(qr2) == 2
        assert qr1.size == 2
        assert qr2.size == 2
        assert qr1[0].index_local == 0
        assert qr1[0].index_global == 0
        assert qr1[1].index_local == 1
        assert qr1[1].index_global == 1
        assert qr2[0].index_local == 0
        assert qr2[0].index_global == 2
        assert qr2[1].index_local == 1
        assert qr2[1].index_global == 3
        assert mbi.map_bit_id['q']['g2l'][0] == ['qr1', 0]
        assert mbi.map_bit_id['q']['g2l'][1] == ['qr1', 1]
        assert mbi.map_bit_id['q']['g2l'][2] == ['qr2', 0]
        assert mbi.map_bit_id['q']['g2l'][3] == ['qr2', 1]
        assert mbi.map_bit_id['q']['l2g']['qr1'][0] == 0
        assert mbi.map_bit_id['q']['l2g']['qr1'][1] == 1
        assert mbi.map_bit_id['q']['l2g']['qr2'][0] == 2
        assert mbi.map_bit_id['q']['l2g']['qr2'][1] == 3

    def test_add_2_cregs(self):
        mbi = qp.MapBitID()
        cr1 = qp.CReg(context=mbi, name='cr1', size=2)
        cr2 = qp.CReg(context=mbi, name='cr2', size=2)
        assert len(cr1) == 2
        assert len(cr2) == 2
        assert cr1.size == 2
        assert cr2.size == 2
        assert cr1[0].index_local == 0
        assert cr1[0].index_global == 0
        assert cr1[1].index_local == 1
        assert cr1[1].index_global == 1
        assert cr2[0].index_local == 0
        assert cr2[0].index_global == 2
        assert cr2[1].index_local == 1
        assert cr2[1].index_global == 3
        assert mbi.map_bit_id['c']['g2l'][0] == ['cr1', 0]
        assert mbi.map_bit_id['c']['g2l'][1] == ['cr1', 1]
        assert mbi.map_bit_id['c']['g2l'][2] == ['cr2', 0]
        assert mbi.map_bit_id['c']['g2l'][3] == ['cr2', 1]
        assert mbi.map_bit_id['c']['l2g']['cr1'][0] == 0
        assert mbi.map_bit_id['c']['l2g']['cr1'][1] == 1
        assert mbi.map_bit_id['c']['l2g']['cr2'][0] == 2
        assert mbi.map_bit_id['c']['l2g']['cr2'][1] == 3
        assert len(mbi.map_bit_id['q']['g2l']) == 0

    def test_add_qreg_creg(self):
        mbi = qp.MapBitID()
        qr = qp.QReg(context=mbi, name='qr1', size=3)
        cr = qp.CReg(context=mbi, name='cr1', size=2)
        assert len(qr) == 3
        assert len(cr) == 2
        assert qr.size == 3
        assert cr.size == 2
        assert qr[0].index_local == 0
        assert qr[0].index_global == 0
        assert qr[1].index_local == 1
        assert qr[1].index_global == 1
        assert qr[2].index_local == 2
        assert qr[2].index_global == 2
        assert cr[0].index_local == 0
        assert cr[0].index_global == 0
        assert cr[1].index_local == 1
        assert cr[1].index_global == 1
        assert mbi.map_bit_id['q']['g2l'][0] == ['qr1', 0]
        assert mbi.map_bit_id['q']['g2l'][1] == ['qr1', 1]
        assert mbi.map_bit_id['q']['g2l'][2] == ['qr1', 2]
        assert mbi.map_bit_id['q']['l2g']['qr1'][0] == 0
        assert mbi.map_bit_id['q']['l2g']['qr1'][1] == 1
        assert mbi.map_bit_id['q']['l2g']['qr1'][2] == 2
        assert mbi.map_bit_id['c']['g2l'][0] == ['cr1', 0]
        assert mbi.map_bit_id['c']['g2l'][1] == ['cr1', 1]
        assert mbi.map_bit_id['c']['l2g']['cr1'][0] == 0
        assert mbi.map_bit_id['c']['l2g']['cr1'][1] == 1