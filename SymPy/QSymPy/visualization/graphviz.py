import importlib
import sympy
import pathlib

spec = importlib.util.find_spec("graphviz")
if spec is None:
    pass
else:
    import graphviz

def create_dot_string_from_expression(expr) -> str:
    return sympy.dotprint(expr)

def render_from_expression(expr, format:str='svg', filenamepath:str|pathlib.Path|None=None) -> None:
    dot = graphviz.Source(create_dot_string_from_expression(expr))
    dot.render(filename=filenamepath, format=format)