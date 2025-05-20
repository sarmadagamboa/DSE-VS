import pytest
import numpy as np
from sustainability_tool import calculate_sustainability

def test_sustainability_tool():

    sustainability_impact = calculate_sustainability()

    print(sustainability_impact)