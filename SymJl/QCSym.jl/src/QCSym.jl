module QCSym
import Symbolics
Symbolics.@register_derivative complex(a,b) 1 Symbolics.SConst(1.0)
Symbolics.@register_derivative complex(a,b) 2 Symbolics.SConst(im)
include("./macros.jl")
include("./Bits_Regs.jl")
include("./Gates.jl")
include("./Circuits.jl")
include("./States.jl")
end