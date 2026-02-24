import Symbolics
include("../QSymJl/QSymJl.jl")

qc = QSymJl.Circuits.QuantumCircuit(name="SuperDuperNiceCircuit")
println(qc)
q_reg_1 = QSymJl.Circuits.add_qreg(qc, "q_reg_1", 2)
c_reg_1 = QSymJl.Circuits.add_creg(qc, "c_reg_1", 2)
println(qc)

QSymJl.Circuits.add_gate(qc, QSymJl.Gates.H_Gate; qubits_t=[q_reg_1[1]], step=1, is_treat_numeric_only=false)
println(qc)
# QSymJl.Circuits.add_gate(qc, QSymJl.Gates.H_Gate; name_prefix="H2", qubits_t=[q_reg_1[2]], step=1, is_treat_numeric_only=false)
# QSymJl.Circuits.add_gate(qc, QSymJl.Gates.H_Gate; name_prefix="H3", qubits_t=[q_reg_1[1]], step=2, is_treat_numeric_only=false)

QSymJl.Circuits.add_gate(qc, QSymJl.Gates.CX_Gate; qubits_t=[q_reg_1[2]], qubits_c=[q_reg_1[1]], step=2, is_treat_numeric_only=false)

unitary = QSymJl.Circuits.assemble_symbolic_unitary(qc, false, false)
display(Symbolics.simplify(unitary))

unitary2 = Symbolics.substitute.(unitary, (Dict(zip(qc.gatecollection.collections[QSymJl.Gates.CX_Gate][1].atomics, qc.gatecollection.collections[QSymJl.Gates.CX_Gate][1].matrix_numeric)),))
display(unitary2)
#println(qc.gatecollection.collections[QSymJl.Gates.CX_Gate])

