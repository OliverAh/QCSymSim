import MacroTools
"""
Copied from Lazy.jl

    @forward T.x functions...

Define methods for `functions` on type `T`, which call the relevant function
on the field `x`.

# Example

```julia
struct Wrapper
    x
end

@forward Wrapper.x  Base.sqrt                                  # now sqrt(Wrapper(4.0)) == 2.0
@forward Wrapper.x  Base.length, Base.getindex, Base.iterate   # several forwarded functions are put in a tuple
@forward Wrapper.x (Base.length, Base.getindex, Base.iterate)  # equivalent to above
```
"""

macro forward(ex, fs)
  MacroTools.@capture(ex, T_.field_) || error("Syntax: @forward T.x f, g, h")
  T = esc(T)
  fs = isexpr(fs, :tuple) ? map(esc, fs.args) : [esc(fs)]
  :($([:($f(x::$T, args...; kwargs...) = (Base.@_inline_meta; $f(x.$field, args...; kwargs...)))
       for f in fs]...);
    nothing)
end


