module KroneckerProduct
import Symbolics
import Logging
logger = Base.CoreLogging.ConsoleLogger(stdout)

#public kron
export kron

function kron(ars...)
    println("ended up there")
    Logging.@logmsg Logging.LogLevel(-10) "Computing Kronecker product of $(length(ars)) arrays."
    return outerproduct_length(ars...)
end

function outerproduct_length(ars...)
    Logging.@logmsg Logging.LogLevel(-10) "entered outerproduct_length"
    Logging.@logmsg Logging.LogLevel(-10000) "Arrays: $(ars)"
    #return 0
    shapes = [size(ar) for ar in ars]
    Logging.@logmsg Logging.LogLevel(-1000) "Shapes: $(shapes)"
    rows = Tuple(collect(s[1] for s in shapes))
    cols = Tuple(collect(s[2] for s in shapes))
    Logging.@logmsg Logging.LogLevel(-1000) "Rows: $(rows)"
    total_rows = prod(rows)
    total_cols = prod(cols)
    Logging.@logmsg Logging.LogLevel(-1000) "Total rows: $(total_rows)"
    C = Symbolics.zeros(typeof(ars[1][1,1]), total_rows, total_cols)
    C = Symbolics.scalarize(C)
    for cids in Iterators.product(CartesianIndices(rows[1:end-1]), CartesianIndices(cols[1:end-1]))
        crows = Tuple(cids[1])
        ccols = Tuple(cids[2])
        Logging.@logmsg Logging.LogLevel(-1000) "crows: $(crows)"
        Logging.@logmsg Logging.LogLevel(-1000) "ccols: $(ccols)"
        rowid = sum((2 .^collect((length(ars)-1):-1:1)) .* (crows .-1))
        colid = sum((2 .^collect((length(ars)-1):-1:1)) .* (ccols .-1))
        glob_row_id = rowid+1:rowid+size(ars[end],1)
        glob_col_id = colid+1:colid+size(ars[end],2)
        Logging.@logmsg Logging.LogLevel(-1000) "glob_row_id: $(glob_row_id)"
        Logging.@logmsg Logging.LogLevel(-1000) "glob_col_id: $(glob_col_id)"
        coeff = prod(ar[crows[i], ccols[i]] for (i,ar) in enumerate(ars[1:end-1]))
        Logging.@logmsg Logging.LogLevel(-1000) "Coefficient: $(coeff)"
        C[glob_row_id, glob_col_id] = coeff .* ars[end][:,:]
    end
    
    Logging.@logmsg Logging.LogLevel(-10000) "C: $(C)"
    return C
end

end # module KroneckerProduct
