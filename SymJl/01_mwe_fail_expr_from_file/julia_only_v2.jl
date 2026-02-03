#############
#
#         a = @syms (a::Complex{Real})
#         C = Symbolics.zeros(typeof(a[1]), 4, 4)
#         C[1,1] = a[1]*a[1]
#         C
#
#############





using Symbolics
using Base.Threads

qs = collect(1:3)
ss = collect(1:1)

Hs = Dict((q,s) => eval(Meta.parse("@syms (H_q$(q)_s$(s)::Complex{Real})[1:2,1:2]")) for q in qs for s in ss)
Hs = Dict(k => Symbolics.scalarize(v[1]) for (k,v) in Hs)
Hs_subs = Dict(k => Symbolics.substitute(v, Dict(Hs[(q,s)][2,1] => 10e-35+1im for q in qs for s in ss)) for (k,v) in Hs)

function outerproduct(A,B)
    rowsA, colsA = size(A)
    rowsB, colsB = size(B)
    C = Symbolics.zeros(typeof(A[1,1]), rowsA*rowsB, colsA*colsB)
    C = Symbolics.scalarize(C)
    #for i in 1:rowsA, j in 1:colsA, k in 1:rowsB, l in 1:colsB
    @threads for index in CartesianIndices((rowsA,colsA,rowsB,colsB))
        i, j, k, l = Tuple(index)
        C[(i-1)*rowsB + k, (j-1)*colsB + l] = A[i,j]*B[k,l]
    end
    return C
end

function outerproduct_length(ars...)
    println(ars)
    shapes = [size(ar) for ar in ars]
    println(shapes)
    rows = Tuple(collect(s[1] for s in shapes))
    cols = Tuple(collect(s[2] for s in shapes))
    println(rows)
    total_rows = prod(rows)
    total_cols = prod(cols)
    println(total_rows)
    C = Symbolics.zeros(typeof(ars[1][1,1]), total_rows, total_cols)
    C = Symbolics.scalarize(C)
    for cids in Iterators.product(CartesianIndices(rows[1:end-1]), CartesianIndices(cols[1:end-1]))
        crows = Tuple(cids[1])
        ccols = Tuple(cids[2])
        println(crows)
        println(ccols)
        rowid = sum((2 .^collect((length(ars)-1):-1:1)) .* (crows .-1))
        colid = sum((2 .^collect((length(ars)-1):-1:1)) .* (ccols .-1))
        glob_row_id = rowid+1:rowid+size(ars[end],1)
        glob_col_id = colid+1:colid+size(ars[end],2)
        println("Global row id: ", glob_row_id)
        println("Global col id: ", glob_col_id)
        coeff = prod(ar[crows[i], ccols[i]] for (i,ar) in enumerate(ars[1:end-1]))
        println("Coefficient: ", coeff)
        C[glob_row_id, glob_col_id] = coeff .* ars[end][:,:]
    end
    
    display(C)
    return C
end


unitary = Symbolics.ones(typeof(Hs[(1,1)][1,1]), 2^maximum(qs), 2^maximum(qs))
#for s in ss
#    println("Processing s = ", s)
#    tmp_unitary = Hs[(1,s)]
#    for q in 2:maximum(qs)
#        tmp_unitary = outerproduct(tmp_unitary, Hs[(q,s)])
#    end
#    if s > 1
#        global unitary = tmp_unitary * unitary
#    else
#        global unitary = tmp_unitary
#    end
#end

for s in ss
    println("Processing s = ", s)
    tmp_unitary = outerproduct_length((Hs[(q,s)] for q in qs)...)
    size(tmp_unitary)
    global unitary *= tmp_unitary
end

println("Constructed unitary")

dict_subs = Dict()
for h in collect(values(Hs))[collect(1:2)]
    dict_subs[h[1,1]] = 1/sqrt(2)
    dict_subs[h[1,2]] = 1/sqrt(2)
    dict_subs[h[2,1]] = 1/sqrt(2)
    dict_subs[h[2,2]] = -1/sqrt(2)
end
#display(dict_subs)
unitary_subs = substitute(unitary, dict_subs, fold=Val(false))
println("Substituted unitary")
display(unitary_subs)


