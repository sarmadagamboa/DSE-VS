import pytest
import numpy as np
import matplotlib.pyplot as plt
import Structures2 as st

def test_Polygon_element_lengths():
    ### INPUTS ###
    x = [0, 1, 1, 0] #x-coordinates of the polygon points
    y = [0, 0, 1, 1] #y-coordinates of the polygon points
    stringers = [0, 0, 0, 0] #number of stringers per element
    t = np.linspace(0.0001, 0.006, 50) #thickness (range) of the skin

    test_box = st.Polygon(height=4.5, px=x, py=y, stringers=stringers, thickness=t)

    expected_lengths = [1, 1, 1, 1]

    assert test_box.elements == expected_lengths



if __name__ == "__main__":
    pytest.main([__file__])