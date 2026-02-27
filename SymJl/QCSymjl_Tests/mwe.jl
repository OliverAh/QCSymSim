import Symbolics
import QCSym

#include("../QSymJl/QSymJl.jl")
#include("../QSym.jl/src/QSym.jl")

qc = QCSym.Circuits.QuantumCircuit(name="SuperDuperNiceCircuit")
println(qc)
q_reg_1 = QCSym.Circuits.add_qreg(qc, "q_reg_1", 2)
c_reg_1 = QCSym.Circuits.add_creg(qc, "c_reg_1", 2)
println(qc)

QCSym.Circuits.add_gate(qc, QCSym.Gates.H_Gate; qubits_t=[q_reg_1[1]], step=1, is_treat_numeric_only=false)
#println(qc)
 QCSym.Circuits.add_gate(qc, QCSym.Gates.H_Gate; name_prefix="H2", qubits_t=[q_reg_1[2]], step=1, is_treat_numeric_only=false)
# QCSym.Circuits.add_gate(qc, QCSym.Gates.H_Gate; name_prefix="H3", qubits_t=[q_reg_1[1]], step=2, is_treat_numeric_only=false)

QCSym.Circuits.add_gate(qc, QCSym.Gates.CX_Gate; qubits_t=[q_reg_1[2]], qubits_c=[q_reg_1[1]], step=2, is_treat_numeric_only=false)

QCSym.Circuits.add_gate(qc, QCSym.Gates.U_Gate; qubits_t=[q_reg_1[1]], step=3, is_treat_numeric_only=false, is_treat_alt_only=true)

unitary = QCSym.Circuits.assemble_symbolic_unitary(qc, true, false)
display(Symbolics.simplify(unitary))

#unitary2 = Symbolics.substitute.(unitary, (Dict(zip(qc.gatecollection.collections[QCSym.Gates.CX_Gate][1].atomics, qc.gatecollection.collections[QCSym.Gates.CX_Gate][1].matrix_numeric)),))
#display(unitary2)
#println(qc.gatecollection.collections[QCSym.Gates.CX_Gate])

