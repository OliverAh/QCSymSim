using Symbolics

mutable struct MapBitID
    map_bit_id::Dict{String, Dict{String, Any}}
    
    function MapBitID()
        return new(Dict(
            "q" => Dict("g2l" => Vector{Vector{Any}}(), "l2g" => Dict("default" => Int[])),
            "c" => Dict("g2l" => Vector{Vector{Any}}(), "l2g" => Dict("default" => Int[]))
            ))
    end
end

has_qreg(ctx::MapBitID, name::String) = haskey(ctx.map_bit_id["q"]["l2g"], name)
has_creg(ctx::MapBitID, name::String) = haskey(ctx.map_bit_id["c"]["l2g"], name)
num_qbits(ctx::MapBitID) = length(ctx.map_bit_id["q"]["g2l"])
num_cbits(ctx::MapBitID) = length(ctx.map_bit_id["c"]["g2l"])

function show(ctx::MapBitID)
    s = ""
    for (t, x) in ctx.map_bit_id
        for (k, v) in x
            s *= "$t $k: $v\n"
        end
    end
    println(s)
end

abstract type AbstractBit end
abstract type AbstractClassicalBit <: AbstractBit end
abstract type AbstractQuantumBit <: AbstractBit end

abstract type AbstractRegister end
abstract type AbstractClassicalRegister <: AbstractRegister end
abstract type AbstractQuantumRegister <: AbstractRegister end


mutable struct Bit <: AbstractBit
    name::String
    index_local::UInt
    index_global::UInt
    is_quantum::Bool
end

function Bit(; context::MapBitID, index_local::Int=-1, is_quantum::Bool=false, name::String="default")

    t = is_quantum ? "q" : "c"

    if index_local == -1
        @assert name == "default" "If no local index is given, the name must be 'default'"
        index_local = length(context.map_bit_id[t]["l2g"]["default"])
    end

    index_global = length(context.map_bit_id[t]["g2l"])
    push!(context.map_bit_id[t]["g2l"], [name, index_local])

    if !haskey(context.map_bit_id[t]["l2g"], name)
        context.map_bit_id[t]["l2g"][name] = Int[]
    end
    push!(context.map_bit_id[t]["l2g"][name], index_global)

    return Bit(name, index_local, index_global, is_quantum)
end

Base.show(io::IO, b::AbstractBit) = print(io, "$(nameof(typeof(b)))[$(b.index_local), $(b.index_global)]")

struct CBit <: AbstractClassicalBit
    
    function CBit(; context::MapBitID)
        return Bit(context=context)
    end
end

struct QBit <: AbstractQuantumBit
    function QBit(; context::MapBitID)
        return Bit(context=context, is_quantum=true)
    end
end


mutable struct BitRegister <: AbstractRegister
    name::String
    size::Int
    is_quantum::Bool
    bits::Vector{AbstractBit}
    
    function BitRegister(; context::MapBitID, name::String, size::Int, is_quantum::Bool)
        @assert name in ("", "default") "Register name cannot be '' or 'default'"
        bits = [Bit(context=context, index_local=-1, is_quantum=is_quantum, name="default") for i in 0:size-1]
    return new(name, size, is_quantum, bits)
    end
end

Base.show(io::IO, r::AbstractRegister) = print(io, "Register[$(r.name) is_quantum:$(r.is_quantum) bits:$(r.bits)]")

struct QReg <: AbstractQuantumRegister
    function QReg(; context::MapBitID, name::String="", size::Int=99)
        return BitRegister(context=context, name=name, is_quantum=true, size=size)
    end
end

struct CReg <: AbstractClassicalRegister
    function CReg(; context::MapBitID, name::String="", size::Int=99)
        return BitRegister(context=context, name=name, is_quantum=false, size=size)
    end
end

