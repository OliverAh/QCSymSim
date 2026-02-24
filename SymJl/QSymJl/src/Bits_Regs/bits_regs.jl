using Symbolics


abstract type AbstractBit end
abstract type AbstractClassicalBit <: AbstractBit end
abstract type AbstractQuantumBit <: AbstractBit end

abstract type AbstractRegister end
abstract type AbstractClassicalRegister <: AbstractRegister end
abstract type AbstractQuantumRegister <: AbstractRegister end


struct MapBitID
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

function add_reg_to_context(ctx::MapBitID, name::String, size::Int, is_quantum::Bool)
    if is_quantum
        t = "q"
    else
        t = "c"
    end
    @assert !haskey(ctx.map_bit_id[t]["l2g"], name) "Register if type $t with name $name already exists"
    _num_bits = length(ctx.map_bit_id[t]["g2l"])
    ctx.map_bit_id[t]["l2g"][name] = collect(range(start=_num_bits+1, length=size))
    for i in 1:size
        push!(ctx.map_bit_id[t]["g2l"], [name, i])
    end
    return is_quantum ? QReg(context=ctx, name=name, size=size) : CReg(context=ctx, name=name, size=size)
end

add_qreg(ctx::MapBitID, name::String, size::Int) = add_reg_to_context(ctx, name, size, true)
add_creg(ctx::MapBitID, name::String, size::Int) = add_reg_to_context(ctx, name, size, false)

Base.show(io::IO, ctx::MapBitID) = begin
    s = "MapBitID:\n"
    for (t, x) in ctx.map_bit_id
         s *= "   Type $t:\n"
        for (k, v) in x
            s *= "      $k: $v\n"
        end
    end
    print(io, s)
end

function _get_global_bit_index_in_context(;ctx::MapBitID, name_reg::String, index_local::Int, is_quantum::Bool)
    t = is_quantum ? "q" : "c"
    @assert haskey(ctx.map_bit_id[t]["l2g"], name_reg) "Register of type $t with name $name_reg does not exist in context"
    @assert index_local >= 1 "Local index $index_local does not exist in register $name_reg of type $t in context"
    return ctx.map_bit_id[t]["l2g"][name_reg][index_local]
end

macro _insert_fields_AbstractBit()
    quote
        $(esc(:(context::MapBitID)))
        $(esc(:(name_reg::String)))
        $(esc(:(index_local::Int)))
        $(esc(:(index_global::Int)))
        $(esc(:(is_quantum::Bool)))
    end
end

mutable struct Bit <: AbstractBit
    @_insert_fields_AbstractBit()
    function Bit(; context::MapBitID, index_local::Int=-1, is_quantum::Bool=false, name_reg::String="default")
        if is_quantum
            return QBit(context=context, index_local=index_local, name_reg=name_reg)
        else
            return CBit(context=context, index_local=index_local, name_reg=name_reg)
        end
    end
end

Base.show(io::IO, b::AbstractBit) = print(io, "$(nameof(typeof(b)))[$(b.index_local), $(b.index_global)]")


struct CBit <: AbstractClassicalBit
    @_insert_fields_AbstractBit()
    function CBit(; context::MapBitID, name_reg::String="default", index_local::Int)
        checked_index_global = _get_global_bit_index_in_context(ctx=context, name_reg=name_reg, index_local=index_local, is_quantum=false)
        return new(context, name_reg, index_local, checked_index_global, false)
    end
end

struct QBit <: AbstractQuantumBit
    @_insert_fields_AbstractBit()
    function QBit(; context::MapBitID, name_reg::String="default", index_local::Int=-1)
        checked_index_global = _get_global_bit_index_in_context(ctx=context, name_reg=name_reg, index_local=index_local, is_quantum=true)
        return new(context, name_reg, index_local, checked_index_global, true)
    end
end

mutable struct BitRegister{T<:AbstractBit} <: AbstractRegister
    context::MapBitID
    name::String
    size::Int
    is_quantum::Bool
    bits::Vector{T}
    
    function BitRegister(; context::MapBitID, name::String, size::Int, is_quantum::Bool)
        _T = is_quantum ? AbstractQuantumBit : AbstractClassicalBit
        if is_quantum && !has_qreg(context, name)
            #add_qreg(context, name, size)
            @assert !has_qreg(context, name) "Quantum register with name $name already exists in context"
        elseif !is_quantum && !has_creg(context, name)
            #add_creg(context, name, size)
            @assert !has_creg(context, name) "Classical register with name $name already exists in context"
        end
        bits = [Bit(context=context, index_local=i, name_reg=name, is_quantum=is_quantum) for i in 1:size]
        return new{_T}(context, name, size, is_quantum, bits)
    end
end

function QReg(; context::MapBitID, name::String, size::Int)
    return BitRegister(context=context, name=name, size=size, is_quantum=true)
end
function CReg(; context::MapBitID, name::String, size::Int)
    return BitRegister(context=context, name=name, size=size, is_quantum=false)
end

function Base.getindex(br::BitRegister, k::Int)
    @assert 1 <= k <= br.size "Index $k is out of bounds for register $(br.name)"
    return br.bits[k]
end


Base.show(io::IO, r::AbstractRegister) = print(io, "Register[$(r.name) is_quantum:$(r.is_quantum) bits:$(r.bits)]")

