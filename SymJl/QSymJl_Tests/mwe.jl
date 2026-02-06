import Symbolics
include("../QSymJl/QSymJl.jl")

h = QSymJl.Gates.H_Gate("H", [1], UInt(1), false)
cx = QSymJl.Gates.CX_Gate("CX", [1], [2], UInt(1), false)
println(h)
println(cx)