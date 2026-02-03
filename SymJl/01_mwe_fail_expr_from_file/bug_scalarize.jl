using Symbolics
using .SymbolicUtils

@variables H[1:2,1:2]
H_scal = Symbolics.scalarize(H)

H_var = Dict(1 => Symbolics.variables(:H_var,1:2,1:2))

H_var_subs = Dict(k => Symbolics.substitute(v, Dict(H_var[1][2,1] => 1.0)) for (k,v) in H_var)

qs = collect(1:2)
ss = collect(1:2)
Hs = Dict((q,s) => Symbolics.variables(Symbol("H_q",q,"_s",s), 1:2,1:2, T=Complex{Real}) for q in qs for s in ss)
HS_subs = Dict((q,s) => Symbolics.substitute(Hs[(q,s)], Dict(Hs[(q,s)][2,1] => 1.0+0im)) for q in qs for s in ss)

@variables (Z::Complex)[1:2,1:2]
Z = Symbolics.scalarize(Z)
Symbolics.substitute(Z, Dict(Z[2,1] => 1.0+0im))

@variables Z::Complex
typeof(Z)
Symbolics.substitute(Z, Dict(Z => 1.0+0im))

@syms Z::Complex{Real}
@syms a b
typeof(Z)
Symbolics.substitute(Z, Dict(Z => a+b))
typeof(Z)