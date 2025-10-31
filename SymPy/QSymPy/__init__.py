import importlib

spec = importlib.util.find_spec("dotenv")
if spec is not None:
    from dotenv import load_dotenv
    load_dotenv(override=True, verbose=True)
    print("Loaded environment variables from .env file")

from .bits_regs import *
from .gate_bases import *
from .gate_defs import *
from .core import *
spec = importlib.util.find_spec("graphviz")
if spec is not None:
    from .visualization import graphviz

#from .hhl import HHL

