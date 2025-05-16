import numpy as np

data = {
    "Designs":         ["LRI-ACC", "LRI-CAI", "QGG", "DT"],
    "Dry_mass":        [529, 656, 665, 440],       # in kg (lower is better)
    "Power":           [763, 991, 1128, 614],      # in W (lower is better)
    "Cost":            [561, 668.5, 474, 409],     # in k€ (lower is better)
    "Science":         [200, 340, 150, 140],       # D/O (higher is better)
    "Sustainability":  [3.5, 2, 2, 2.5],           # (higher is better)
    "Risk":            [0.03, 0.04, 0.06, 0.03],   # failure prob (lower is better)
}

# Weights for each criterion (must sum to 1.0)
weights = {
    "Dry_mass":       1.0,
    "Power":          1.0,
    "Cost":           4.0,
    "Science":        7.0,
    "Sustainability": 2.0,
    "Risk":           3.0,
}

# Specify which metrics are higher-is-better
higher_better = {
    "Dry_mass":       False,
    "Power":          False,
    "Cost":           False,
    "Science":        True,   # <-- treat Science as higher is better
    "Sustainability": True,
    "Risk":           False,
}

def normalize_minmax(values, scale=5, higher_is_better=False):
    """
    Min–max normalize a 1D array to the range [0, scale].
    
    If higher_is_better:
        score = scale * (arr - min) / (max - min)
    else:
        score = scale * (max - arr) / (max - min)
    """
    arr = np.array(values, dtype=float)
    mn, mx = arr.min(), arr.max()
    if mx == mn:
        return np.full_like(arr, scale/2)  # all the same -> mid score
    if higher_is_better:
        return scale * arr / mx
    else:
        return scale * (1 - arr / mx + mn / mx)

def main():
    # 1) Normalize each metric
    norm = {}
    for key in weights:
        norm[key] = normalize_minmax(
            data[key],
            scale=5,
            higher_is_better=higher_better[key]
        )

    # 2) Compute weighted scores
    total_scores = []
    for i, design in enumerate(data["Designs"]):
        score = 0.0
        for key, w in weights.items():
            score += norm[key][i] * w
        total_scores.append(score)

    # 3) Print out a table of results
    header = f"{'Design':8s} " + " ".join(f"{k:15s}" for k in weights) + " Total"
    print(header)
    print("-" * len(header))
    for i, design in enumerate(data["Designs"]):
        row = f"{design:8s} " + " ".join(f"{norm[k][i]:15.2f}" for k in weights)
        row += f" {total_scores[i]:.2f}"
        print(row)

if __name__ == "__main__":
    main()
