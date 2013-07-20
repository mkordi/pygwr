# We need to add gwstatsmodels package to our Python path
# Otherwise, gwstatsmodels won't be able to refer to himself.
import os, sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

from gwr import *
from gaussian import Gaussian

__version__ = '1.0'
