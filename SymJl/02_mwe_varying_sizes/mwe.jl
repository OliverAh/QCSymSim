import Symbolics
include("kronecker_product.jl")
#import .KroneckerProduct
#using .KroneckerProduct


import Logging                                  # Required for logging 1/3
ENV["JULIA_DEBUG"] = :KroneckerProduct          # Required for logging 2/3
Logging.disable_logging(Logging.LogLevel(-100)) # Required for logging 3/3, disables all logs <=-100


qs = collect(1:7)
ss = collect(1:1)

Hs = Dict((q,s) => eval(Meta.parse("Symbolics.@syms (H_q$(q)_s$(s)::Complex{Real})[1:2,1:2]")) for q in qs for s in ss)
Hs = Dict(k => Symbolics.scalarize(v[1]) for (k,v) in Hs)

C = KroneckerProduct.kron(values(Hs)...)

println("finished kron")
Symbolics.scalarize(C)
