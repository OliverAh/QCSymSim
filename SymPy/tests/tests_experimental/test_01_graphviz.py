import pathlib
import sys
import sympy as sp
import numpy as np
sys.path.append(str(pathlib.Path(__file__).parent.parent.parent))

import pytest
import QSymPy as qp

import importlib
spec = importlib.util.find_spec("graphviz")
is_graphviz_available = spec is None
skip_reason = "Graphviz is not available"

class TestGraphvizVisualization:
    @pytest.mark.skipif(is_graphviz_available, reason=skip_reason)
    def test_create_dot_string_from_expression(self):
        qc = qp.QuantumCircuit(num_qubits=2)
        qc.add_gate(name='I', qubits_t=[0], step=0)
        qc.assemble_symbolic_unitary()

        dot_string = qp.visualization.graphviz.create_dot_string_from_expression(qc.unitary)
        assert isinstance(dot_string, str)
    
    @pytest.mark.skipif(is_graphviz_available, reason=skip_reason)
    def test_render_from_expression(self):
        qc = qp.QuantumCircuit(num_qubits=2)
        qc.add_gate(name='H', qubits_t=[0], step=0)
        qc.add_gate(name='CX', qubits_t=[1], qubits_c=[0], step=1)
        qc.assemble_symbolic_unitary()

        # Render to a temporary file
        temp_filepath = pathlib.Path(__file__).parent.joinpath("temp_qc_graph")
        tmp_format = 'svg'
        qp.visualization.graphviz.render_from_expression(qc.unitary, format=tmp_format, filenamepath=temp_filepath)

        # Check if the file was created
        assert temp_filepath.with_suffix('.' + tmp_format).exists()

        # Clean up temporary files
        temp_filepath.with_suffix('.' + tmp_format).unlink()
        temp_filepath.unlink()

    @pytest.mark.skipif(is_graphviz_available, reason=skip_reason)
    def test_ImmutableMatrix(self):
        expr = sp.ImmutableMatrix([[1, 2], [3, 4]])
        assert isinstance(expr, sp.ImmutableDenseMatrix)
        dot_string = qp.visualization.graphviz.create_dot_string_from_expression(expr)
        assert isinstance(dot_string, str)
    
    @pytest.mark.skipif(is_graphviz_available, reason=skip_reason)
    def test_MutableMatrix(self):
        expr = sp.Matrix([[1, 2], [3, 4]])
        assert isinstance(expr, sp.MutableDenseMatrix)
        with pytest.raises(AttributeError) as e:
            dot_string = qp.visualization.graphviz.create_dot_string_from_expression(expr)
        assert "'MutableDenseMatrix' object has no attribute 'args'" in str(e.value)


    