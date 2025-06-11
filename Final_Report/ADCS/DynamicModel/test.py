import pandas as pd

# Load CSV
df = pd.read_csv("gain_sweep_results3.csv")

# Total deviation in mrad
df['TotalDeviation'] = df['MaxDevRoll_mrad'] + df['MaxDevPitch_mrad'] + df['MaxDevYaw_mrad']

# Total torque (integral of absolute torque over time)
df['TotalTorque'] = df['TotalTorqueX_Nm_s'] + df['TotalTorqueY_Nm_s'] + df['TotalTorqueZ_Nm_s']

# Acceptable deviation threshold (in mrad)
threshold = 100

# Filter acceptable solutions
acceptable_df = df[
    (df['MaxDevRoll_mrad'] < threshold) &
    (df['MaxDevPitch_mrad'] < threshold) &
    (df['MaxDevYaw_mrad'] < threshold)
]

# Normalize all three metrics
acceptable_df['NormDeviation'] = (acceptable_df['TotalDeviation'] - acceptable_df['TotalDeviation'].min()) / (acceptable_df['TotalDeviation'].max() - acceptable_df['TotalDeviation'].min())
acceptable_df['NormTorque'] = (acceptable_df['TotalTorque'] - acceptable_df['TotalTorque'].min()) / (acceptable_df['TotalTorque'].max() - acceptable_df['TotalTorque'].min())
acceptable_df['NormActivation'] = (acceptable_df['TorqueActivation'] - acceptable_df['TorqueActivation'].min()) / (acceptable_df['TorqueActivation'].max() - acceptable_df['TorqueActivation'].min())

# Weighted cost function
# You can adjust weights depending on priorities
# Example: give 0.05 to Deviation, 0.80 to Torque, 0.15 to Activation
acceptable_df['Cost'] = (
    0.6 * acceptable_df['NormDeviation'] +
    0.4 * acceptable_df['NormTorque'] +
    0.0 * acceptable_df['NormActivation']
)

# Find best combination (minimum cost)
best_index = acceptable_df['Cost'].idxmin()
best_row = acceptable_df.loc[best_index]

print("Best combination:")
print(best_row[['Kp', 'Kd', 'MaxDevRoll_mrad', 'MaxDevPitch_mrad', 'MaxDevYaw_mrad', 'TotalTorque', 'TorqueActivation']])

# Show top 10 combinations
top_10 = acceptable_df.sort_values(by='Cost').head(10)
print(top_10[['Kp', 'Kd', 'TotalDeviation', 'TotalTorque', 'TorqueActivation', 'Cost']])
