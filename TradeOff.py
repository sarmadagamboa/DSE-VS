import numpy as np
import random
import matplotlib.pyplot as plt
from beautifultable import BeautifulTable


data = {
    "Designs":                  ["LRI-ACC", "LRI-CAI", "QGG", "DT"],
    "Dry_mass":                 [529, 656, 665, 440],                 # in kg (lower is better)
    "Power":                    [763, 991, 1128, 614],                # in W (lower is better)
    "Cost":                     [527.8, 645.0, 417.1, 317.9],     # in M€ (lower is better)
    "D/O":                      [115, 160, 40, 20],                   # D/O (higher is better) 
    "Temporal sensitivity":     [5, 7, 0, 2],             # error inversed (higher is better)(10+)
    "Sustainability":           [2.2, 1, 2.8, 4],                     # (higher is better)
    "Risk":                     [0.8, 1, 0.6, 0.2],                   # (lower is better)
}

# Specify which metrics are higher-is-better
higher_better = {
    "Dry_mass":               False,
    "Power":                  False,
    "Cost":                   False,
    "D/O":                    True,   
    "Temporal sensitivity":   True,  
    "Sustainability":         True,
    "Risk":                   False,
}

# Weights for each criterion 
weights = {
    "Dry_mass":               1.0,
    "Power":                  1.0,
    "Cost":                   3.0,
    "D/O":                    5.0,
    "Temporal sensitivity":   5.0,
    "Sustainability":         2.0,
    "Risk":                   3.0,
}


def sensitivity_noise(data, weights, data_sensitivity=False, weight_sensitivity=False):
    """
    Compute the weighted scores for each design.
    
    """
    if data_sensitivity:
        for key in data:
            for i in range(len(data[key])):
                data[key][i] = random.uniform(0.8, 1.2) * data[key][i]
    
    if weight_sensitivity:
        for key in weights:
            weights[key] = random.randint(-1, 1) + weights[key]

    return data, weights


def sensitivity_range(weights):
    minus = 1.1
    plus = 1.1
    step = 0.05
    
    sens_weights = {}
    for key in weights:
        sens_weights[key] = np.arange(weights[key]-minus, weights[key]+plus, step)
    sens_weights_xist = np.arange(-minus, plus, step)
    return sens_weights, sens_weights_xist


def normalize_data(weights, data, higher_is_better, scale=5):
    """
    Min–max normalize a 1D array to the range [0, scale].
    
    """

    norm_data = {}
    norm_data["Designs"] = data["Designs"]
    for key in weights:
        arr = np.array(data[key], dtype=float)
        print(f"{key}, arr: {arr}")
        mn, mx = arr.min(), arr.max()
        print(f"{key}, mn: {mn}, mx: {mx}")
        if higher_is_better[key]:
            norm_data[key] = scale * (arr / mx)
            print(f"{key}, norm: {norm_data[key]}")
        else:
            norm_data[key] = scale * (1 - (arr / mx) + (mn / mx))
            print(f"{key}, norm: {norm_data[key]}")
    
    return norm_data


def compute_weighted_scores(norm_data, weights):
    """
    Compute the weighted scores for each design.

    """
    scores = {}
    total_scores = []
    for i in range(len(norm_data["Designs"])):
        score = 0
        for key, w in weights.items():
            score += norm_data[key][i] * w
        total_scores.append(score)
    
    #winner = scores["Designs"][total_scores.index(max(total_scores))]
    scores["Score"] = total_scores

    return scores


def compute_weighted_scores_sensitivity(norm_data, weights, sens_weights):
    """
    Compute the weighted scores for each design.

    """
    scores = {}
    
    for key in weights:
        temp_weights = weights.copy()
        for k in temp_weights:
            temp_weights[k] = temp_weights[k] * np.ones_like(sens_weights[key])

        temp_weights[key] = sens_weights[key]

        total_scores = []
        for i in range(len(norm_data["Designs"])):
            score = np.zeros_like(sens_weights[key])
            for key2, w in temp_weights.items():
                score += norm_data[key2][i] * w
            total_scores.append(score)

        #winner = scores["Designs"][total_scores.index(max(total_scores))]
        print(key)
        print(total_scores)
        scores[key] = total_scores

    return scores


def print_results(score, norm_data, weights):
    headers = [" ", "Weights", " "]
    for i in range(len(norm_data["Designs"])):
        headers.append(norm_data["Designs"][i])
    table = BeautifulTable()
    table.column_headers = headers
    #print(table)
    for key in weights:
        row = [key, weights[key], " "]
        for i in range(len(norm_data["Designs"])):
            row.append(norm_data[key][i])
        #print(row)
        table.append_row(row)
    
    runner_up_row = [" ", " ", " "]
    for i in range(len(norm_data["Designs"])):
        runner_up_row.append(" ")
    table.append_row(runner_up_row)

    bottom_row = ["Total score", " ", " "]
    for i in range(len(norm_data["Designs"])):
        bottom_row.append(score["Score"][i])
    table.append_row(bottom_row)
    print(table)


def print_sensitivity(scores, sens_weights, data, sens_axis):
    for key in scores:
        plt.plot(sens_axis, scores[key][0],label=data["Designs"][0])
        plt.plot(sens_axis, scores[key][1],label=data["Designs"][1])
        plt.plot(sens_axis, scores[key][2],label=data["Designs"][2])
        plt.plot(sens_axis, scores[key][3],label=data["Designs"][3])
        plt.axvline(0, color='r', linestyle='--', label="Original Weight")
        plt.title(key)
        plt.xlabel("Weight difference")
        plt.ylabel("Score")
        plt.legend()
        plt.show()



if __name__ == "__main__":
    sensitivity = False

    if sensitivity:
        sens_weights, sens_axis = sensitivity_range(weights)
        norm_data = normalize_data(sens_weights, data, higher_better, scale=5)
        scores = compute_weighted_scores_sensitivity(norm_data, weights, sens_weights)
        print_sensitivity(scores, sens_weights, data, sens_axis)
        
    else:
        norm_data = normalize_data(weights, data, higher_better, scale=5)
        scores = compute_weighted_scores(norm_data, weights)
        print_results(scores, norm_data, weights)

