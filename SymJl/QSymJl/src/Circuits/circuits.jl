import Symbolics
import ..BitsRegs
import ..Gates

struct GateCollection
    collections::Dict{typeof(Gates.AbstractGate), Vector{Any}}
end

struct QuantumCircuit
    context::BitsRegs.MapBitID
    qregs::Vector{BitsRegs.QReg}
    cregs::Vector{BitsRegs.CReg}
    gatecollection::GateCollection
    barriers::Vector{Any}
    
    function QuantumCircuit()
        context = BitsRegs.MapBitID()
        qregs = BitsRegs.QReg[]
        cregs = BitsRegs.CReg[]
        gatecollection = GateCollection(Dict())
        barriers = Any[]
        return new(context, qregs, cregs, gatecollection, barriers)
end
end

function add_qreg(circuit::QuantumCircuit, name::String="", size::Int=99)
    @assert !has_qreg(circuit.context["q"]["l2g"], name) "Quantum register with name $name already exists"
    qreg = BitsRegs.QReg(context=circuit.context, name=name, size=size)
    push!(circuit.qregs, qreg)
    return qreg
end

function add_creg(circuit::QuantumCircuit, name::String="", size::Int=99)
    @assert !has_creg(circuit.context["c"]["l2g"], name) "Classical register with name $name already exists"
    creg = BitsRegs.CReg(context=circuit.context, name=name, size=size)
    push!(circuit.cregs, creg)
    return creg
end

function add_gate(circuit::QuantumCircuit, gate::typeof(Gates.AbstractGate); kwargs...)
    if !haskey(circuit.gatecollection.collections, gate)
        circuit.gatecollection.collections[gate] = Any[]
    end
    push!(circuit.gatecollection.collections[gate], gate(; kwargs...))
end