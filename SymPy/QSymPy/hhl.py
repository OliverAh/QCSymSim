import sympy

from . import src

class HHL(src.QuantumCircuit):
    """
    HHL algorithm for solving linear systems of equations.
    """

    def __init__(self, num_qubits: int=1, num_clbits: int=1):
        super().__init__(num_qubits, num_clbits)

    

