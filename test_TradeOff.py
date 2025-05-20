import pytest
import numpy as np
import matplotlib.pyplot as plt
import TradeOff as to

data = {
    "Designs":                  ["ONE", "TWO", "THREE"],
    "Dry_mass":                 [200, 300, 500],                 # higher is better
    "Power":                    [0, 100, 500],                # lower is better
}

# Specify which metrics are higher-is-better
higher_better = {
    "Dry_mass":               True,
    "Power":                  False,
}

# Weights for each criterion 
weights = {
    "Dry_mass":               1.0,
    "Power":                  1.0,
}


def test_normalize_data():
    """
    Test the normalize_data function.
    """
    normalized_data = to.normalize_data(weights, data, higher_better)

    norm_test_data = {
    "Designs":                  ["ONE", "TWO", "THREE"],
    "Dry_mass":                 [2., 3., 5.],
    "Power":                    [5., 4., 0.],
    }
    
    for key in normalized_data:
        assert np.array_equal(normalized_data[key], norm_test_data[key])




def test_sensitivity_range():
    plusminus = 2.1
    step = 0.5

    expected_sens_weights = {
        "Dry_mass": np.arange(0.0, 3.1, 0.5),  # from 0.0 to (1.0 + 2.1) with step 0.5
        "Power": np.arange(0.0, 3.1, 0.5)
    }

    expected_sens_weights_xist = {
        "Dry_mass": np.arange(-1.0, 2.1, 0.5),  # from (0.0 - 1.0) to 2.1
        "Power": np.arange(-1.0, 2.1, 0.5)
    }

    sens_weights, sens_weights_xist = to.sensitivity_range(weights, plusminus, step)

    for key in weights:
        np.testing.assert_array_almost_equal(sens_weights[key], expected_sens_weights[key])
        np.testing.assert_array_almost_equal(sens_weights_xist[key], expected_sens_weights_xist[key])


def test_compute_weighted_scores():
    """
    Test the compute_weighted_scores function.
    """
    norm_data = {
    "Designs":                  ["ONE", "TWO", "THREE"],
    "Dry_mass":                 [2., 3., 5.],
    "Power":                    [5., 4., 0.],
    }

    expected_scores = {
        "Score": [2.0 + 5.0, 3.0 + 4.0, 5.0 + 0.0]
    }
    result = to.compute_weighted_scores(norm_data, weights)
    np.testing.assert_array_almost_equal(result["Score"], expected_scores["Score"])


def test_compute_weighted_scores_sensitivity():
    """
    Test the compute_weighted_scores_sensitivity function.
    """
    
    norm_data = {
    "Designs":                  ["ONE"],
    "Dry_mass":                 [1.],
    "Power":                    [5.],
    }

    sens_weights = {
        "Dry_mass": np.arange(0.0, 3.1, 0.5),  # from 0.0 to (1.0 + 2.1) with step 0.5
        "Power": np.arange(0.0, 3.1, 0.5)
    }
    
    exp_scores = {
    'Dry_mass': [np.array([5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0])],
    'Power': [np.array([1.0, 3.5, 6.0, 8.5, 11.0, 13.5, 16.0])]
    }

    scores = to.compute_weighted_scores_sensitivity(norm_data, weights, sens_weights)
    for key in scores:
        np.testing.assert_array_almost_equal(scores[key], exp_scores[key])



if __name__ == "__main__":
    pytest.main([__file__])
