using Symbolics

example = [
    # "small",
    "intermediate", 
    # "large"
]


sympy_expression = chop(read("inputs_complex_$(example[1])/fancy_string_complex.txt", String), head=7, tail=1)

sympy_variables = chop(read("inputs_complex_$(example[1])/fancy_string_variables.txt", String), tail=2)

eval(Meta.parse(string("@variables ", sympy_variables)))







symbolics_expression = eval(Meta.parse(sympy_expression))

simplified_symbolics_expression = simplify.(symbolics_expression)
isequal(symbolics_expression, simplified_symbolics_expression)

simplified_symbolics_expression_dthetasomething = Symbolics.derivative.(simplified_symbolics_expression, theta_RY_qt0_qc_s2)

isequal(simplified_symbolics_expression_dthetasomething, simplified_symbolics_expression)