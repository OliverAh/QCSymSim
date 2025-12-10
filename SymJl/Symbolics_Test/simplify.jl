using Symbolics

@variables x, y, z

Dx = Differential(x)
Dy = Differential(y)
Dz = Differential(z)


function replace_derivative(text::String)
    # Recursively replace Derivative calls
    # Pattern: Derivative(anything, variable) -> D{variable}(anything)
    while contains(text, r"Derivative\(")
        text = replace(text, r"Derivative\(([^,]+),\s*(\w+)\)" => s"D\2(\1)")
    end
    return text
end

size_i = 6


text = read("inputs/fancy_string_$(size_i).txt", String)

expr_with_differential = replace_derivative(text)
expr_with_derivative = replace(text, "Derivative" => "Symbolics.derivative")

symbolics_expr1 = expand_derivatives(eval(Meta.parse(expr_with_differential)))
symbolics_expr2 = eval(Meta.parse(expr_with_derivative))


symbolics_expr1 === symbolics_expr2