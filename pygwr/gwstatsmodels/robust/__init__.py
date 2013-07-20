"""
Robust statistical models
"""
import norms
from .scale import mad, stand_mad, Huber, HuberScale, hubers_scale

from gwstatsmodels import NoseWrapper as Tester
test = Tester().test
