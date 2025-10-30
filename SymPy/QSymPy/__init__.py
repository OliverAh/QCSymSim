from dotenv import load_dotenv
load_dotenv(override=True, verbose=True)
print("Loaded environment variables from .env file")

from .bits_regs import *
from .gate_bases import *
from .gate_defs import *
from .core import *
#from .hhl import HHL

