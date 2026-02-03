module KroneckerProduct
import Symbolics
import Logging

logger = Logging.ConsoleLogger(stdout)

#public kron
export kron

function kron(ars...; debug::Int<=0)
    println("ended up here")
    Logging.disable_logging(Logging.Info)
    return outerproduct_length(ars...; debug=debug)
end

function kron(ars...; debug::Int=-1)
    println("ended up there")
    Logging.@logmsg Logging.LogLevel(-10) "Computing Kronecker product of $(length(ars)) arrays."
    return outerproduct_length(ars...; debug=debug)
end

function outerproduct_length(ars...; debug::Int=0)
    debug > 0 && println(ars)
    Logging.@logmsg Logging.LogLevel(-10) "entered outerproduct_length"
    return 0
    shapes = [size(ar) for ar in ars]
    debug > 0 && println(shapes)
    rows = Tuple(collect(s[1] for s in shapes))
    cols = Tuple(collect(s[2] for s in shapes))
    debug > 0 && println(rows)
    total_rows = prod(rows)
    total_cols = prod(cols)
    debug > 0 && println(total_rows)
    C = Symbolics.zeros(typeof(ars[1][1,1]), total_rows, total_cols)
    C = Symbolics.scalarize(C)
    for cids in Iterators.product(CartesianIndices(rows[1:end-1]), CartesianIndices(cols[1:end-1]))
        crows = Tuple(cids[1])
        ccols = Tuple(cids[2])
        debug > 0 && println(crows)
        debug > 0 && println(ccols)
        rowid = sum((2 .^collect((length(ars)-1):-1:1)) .* (crows .-1))
        colid = sum((2 .^collect((length(ars)-1):-1:1)) .* (ccols .-1))
        glob_row_id = rowid+1:rowid+size(ars[end],1)
        glob_col_id = colid+1:colid+size(ars[end],2)
        debug > 0 && println("Global row id: ", glob_row_id)
        debug > 0 && println("Global col id: ", glob_col_id)
        coeff = prod(ar[crows[i], ccols[i]] for (i,ar) in enumerate(ars[1:end-1]))
        debug > 0 && println("Coefficient: ", coeff)
        C[glob_row_id, glob_col_id] = coeff .* ars[end][:,:]
    end
    
    debug > 99 && display(C)
    return C
end

end # module KroneckerProduct
