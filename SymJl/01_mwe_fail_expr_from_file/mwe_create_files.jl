using Symbolics

function create_dummy_files(rows::Integer, cols::Integer, summands::Integer)
    # Create dummy variable names
    var_names = [ "x$(i)" for i in 1:(rows*cols*summands) ]
    variables_str = join(var_names[:], ", ")
    variables_str *= "\n"
    
    expression_str = "[\n"
    id = 0
    for i in 1:rows
        expr_row = "["
        for j in 1:cols
            terms = join([ "x$(id + k)" for k in 1:summands ], " + ")
            expr_row *= " " * terms * " "
            id += summands
        end
        expr_row *=  "],\n"
        expression_str *= expr_row
    end
    expression_str *= "]\n"

    #println(variables_str)
    #println(expression_str)
    
    vars = [Symbol(var) for var in var_names]
    #println(vars)

    # Write to files
    open("fancy_string_variables_$(rows)_$(cols)_$(summands).txt", "w") do f
        write(f, variables_str)
    end
    open("fancy_string_expr_$(rows)_$(cols)_$(summands).txt", "w") do f
        write(f, expression_str)
    end
end

n = 4
N = 2^n
create_dummy_files(N,N, n^5)