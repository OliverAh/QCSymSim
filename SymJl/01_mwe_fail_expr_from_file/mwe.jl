using Symbolics

rows, cols, summands = 16, 16, 1024

str_vars = chop(read("fancy_string_variables_$(rows)_$(cols)_$(summands).txt", String), tail=1)
str_expr = chop(read(     "fancy_string_expr_$(rows)_$(cols)_$(summands).txt", String), tail=1)
println("Read variables and expression from files")

eval(Meta.parse(string("@variables ", str_vars)))
println("Parsed & evaluated variables")

metaexpr = Meta.parse(str_expr)
println("Parsed expression")
symbolics_expression = eval(metaexpr)
println("Evaluated expression")