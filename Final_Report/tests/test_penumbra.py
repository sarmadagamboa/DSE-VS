import pytest
from Final_Report.Astrodynamics.penumbra import penumbra_one
from Final_Report.Astrodynamics.penumbra import solar_coverage


def test_penumbra_one():
    # hand calculation
    a, alpha, theta, i_star = penumbra_one()
    assert a == pytest.approx(0.762, abs=1e-3)
    assert alpha == pytest.approx(2.46, abs=1e-2)
    assert theta == pytest.approx(1.22, abs=1e-2)
    assert i_star == pytest.approx(0.482, abs=1e-3)

def test_solar_coverage():
    # hand calculation
    a_star = 0.5 
    alpha, theta, i_star = 2.457347769641925, 1.222565643747163, 0.4824648557872968
    result = solar_coverage(a_star, i_star, alpha, theta)
    assert result == pytest.approx(0.07222, abs=1e-4) 
    assert 0 <= result <= 1  # Check if the result is within expected bounds (0 to 1)



if __name__ == "__main__":
    pytest.main([__file__])