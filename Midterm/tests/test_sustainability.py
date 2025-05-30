import pytest
import numpy as np
from Sustainability.sustainability_tool import calculate_sustainability

def test_sustainability_tool():

    sustainability_impact = calculate_sustainability()

    assert np.isclose(sustainability_impact[0], 58821.52667, rtol=1e-6)
    
def test_sustainability_tool_print():

    calculate_sustainability()

    assert True