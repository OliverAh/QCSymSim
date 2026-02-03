using Symbolics

rows, cols, summands = 16, 16, 1024

@variables (H::Complex)[1:2, 1:2]
@variables (X::Complex)[1:2, 1:2]

H_scal = Symbolics.scalarize(H)
X_scal = Symbolics.scalarize(X)

H_s = Symbolics.substitute.(H, (Dict(H[2,1] => 1.0),))
#H_subs = Symbolics.substitute.(H_scal, (Dict(H_scal[2,1] => 1.0),))
H_subs = Symbolics.substitute(H_scal, Dict(H_scal[2,i] => 1.0 for i in 1:2))
println(H_s)
println(H_subs)
println(H_scal)
println(typeof(H_scal))

@variables (XH::Complex)[1:4,1:4]
XH = Symbolics.scalarize(XH)
for i in 1:2, j in 1:2, k in 1:2, l in 1:2
    XH[(i-1)*2+k,(j-1)*2+l] = X_scal[i,j]*H_scal[k,l]
end

#println("Multiplied variables ")
#println(XH)
#display(Symbolics.scalarize(XH))

qs = collect(1:2)
ss = collect(1:2)
#Hs = Dict((q,s) => Symbolics.scalarize(Symbolics.variables(Symbol("H_q",q,"_s",s), 1:2,1:2, T=Complex{Real})) for q in qs for s in ss)
#Hs = Dict((q,s) => Symbolics.variables(Symbol("H_q",q,"_s",s), 1:2,1:2, T=Complex{Real}) for q in qs for s in ss)
    Hs = Dict((q,s) => Symbolics.scalarize(Symbolics.@variables(Symbol("H_q",q,"_s",s), 1:2,1:2, T=Complex{Real})) for q in qs for s in ss)
println(Hs)
println(Hs[(1,1)][2,1])
println(Symbolics.substitute(Hs[(1,1)], Dict(Hs[(1,1)][2,1] => 1.0)))
#Hs = Symbolics.variables(Symbol("H_q",q,"_s",s), 1:2,1:2, T=Complex{Real})
#println(Hs)
#display(Symbolics.scalarize(Hs))

function outerproduct(A,B)
    rowsA, colsA = size(A)
    rowsB, colsB = size(B)
    C = Symbolics.zeros(Complex, rowsA*rowsB, colsA*colsB)
    for i in 1:rowsA, j in 1:colsA, k in 1:rowsB, l in 1:colsB
        C[(i-1)*rowsB + k, (j-1)*colsB + l] = A[i,j]*B[k,l]
    end
    return C
end

unitary = Symbolics.Matrix{Complex}(Symbolics.I, 2^maximum(qs), 2^maximum(qs))

for s in ss
    println("Processing s = ", s)
    tmp_unitary = Hs[(1,s)]
    for q in 2:maximum(qs)
        tmp_unitary = outerproduct(tmp_unitary, Hs[(q,s)])
    end
    global unitary = tmp_unitary * unitary
end

println("Constructed unitary")
println(Dict(h[1,1] => 1.0 for h in values(Hs)))
unitary_subs = Symbolics.substitute.(unitary, (Dict(h[1,1] => 1.0 for h in values(Hs)),))

println("Substituted unitary")

#println(unitary)
#println(unitary_subs)
