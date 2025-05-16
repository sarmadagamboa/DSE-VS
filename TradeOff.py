import numpy as np

data = {
    "Designs":         ["LRI-ACC", "LRI-CAI", "QGG", "DT"],
    "Dry_mass":        [529, 656, 665, 440],                 # in kg (lower is better)
    "Power":           [763, 991, 1128, 614],                # in W (lower is better)
    "Cost":            [561.49, 668.57, 455.25, 362.85],     # in M€ (lower is better)
    "D/O":             [135, 160, 40, 20],                 # D/O (higher is better)
    "Error":           [1e15, 1e17, 1e10, 1e12],             # errer inversed (higher is better)
    "Sustainability":  [3, 1.5, 1.5, 4],                     # (higher is better)
    "Risk":            [0.6, 1, 0.6, 0.2],                   # (lower is better)
}

# Weights for each criterion (must sum to 1.0)
weights = {
    "Dry_mass":       1.0,
    "Power":          1.0,
    "Cost":           4.0,
    "D/O":            4.0,
    "Error":          3.0,
    "Sustainability": 2.0,
    "Risk":           3.0,
}

# Specify which metrics are higher-is-better
higher_better = {
    "Dry_mass":       False,
    "Power":          False,
    "Cost":           False,
    "D/O":            True,   
    "Error":          True,  
    "Sustainability": True,
    "Risk":           False,
}

def normalize_minmax(values, scale=5, higher_is_better=False):
    """
    Min–max normalize a 1D array to the range [0, scale].
    
    """
    arr = np.array(values, dtype=float)
    mn, mx = arr.min(), arr.max()
    if mx == mn:
        return np.full_like(arr, scale / 2)  # all identical → mid score
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
