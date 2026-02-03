import Symbolics
include("kronecker_product.jl")
#import .KroneckerProduct
#using .KroneckerProduct

qs = collect(1:7)
ss = collect(1:1)

Hs = Dict((q,s) => eval(Meta.parse("Symbolics.@syms (H_q$(q)_s$(s)::Complex{Real})[1:2,1:2]")) for q in qs for s in ss)
Hs = Dict(k => Symbolics.scalarize(v[1]) for (k,v) in Hs)

C = KroneckerProduct.kron(values(Hs)..., debug=1)
#C = kron(values(Hs)..., debug=0)

println("finished kron")
Symbolics.scalarize(C)
