import Symbolics
include("../QSymJl/QSymJl.jl")

h = QSymJl.Gates.H_Gate(name="H", qubits_t=[1], step=UInt(1), is_treat_numeric_only=false)
cx = QSymJl.Gates.CX_Gate(name="CX", qubits_t=[1], qubits_c=[2], step=UInt(1), is_treat_numeric_only=false)
println(h)
println(cx)

qc = QSymJl.Circuits.QuantumCircuit()
QSymJl.Circuits.add_gate(qc, QSymJl.Gates.H_Gate; name="H", qubits_t=[1], step=UInt(1), is_treat_numeric_only=false)

println(qc)

