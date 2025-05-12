import pandas as pd

# criteria
criteria = ["Science return", "Eclipse time", "Solar sight dir. & variability", "ΔV requirement"]
orbit_types = ["SSO", "Polar circ.", "Polar ecc."]

# Raw scores based on trade-off matrix
scores = {
    "SSO": [3, 4, 4, 4],
    "Polar circ.": [5, 2.5, 1, 4],
    "Polar ecc.": [2.5, 3, 1, 2]
}

# Define weights for each criterion (you can change these)
weights = {
    "Science return": 4,
    "Eclipse time": 2,
    "Solar sight dir. & variability": 1,
    "ΔV requirement": 2
}

# Create Dataframe
df_scores = pd.DataFrame(scores, index=criteria)

weighted_df = df_scores.copy()
for criterion in criteria:
    weighted_df.loc[criterion] *= weights.get(criterion, 1)

# Compute total scores
total_scores = weighted_df.sum()

# Find the best orbit
best_orbit = total_scores.idxmax()
best_score = total_scores.max()

# Display results
print("Orbit Trade-Off Analysis")
print(weighted_df)
print("\nTotal Weighted Scores:")
print(total_scores)
print(f"Best Orbit: {best_orbit}, score: {best_score}")

