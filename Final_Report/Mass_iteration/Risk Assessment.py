import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from collections import defaultdict


def create_risk_matrix():
    fig, ax = plt.subplots(figsize=(10, 8))

    #grid for likelihood and conseuqnece
    x = np.linspace(0, 1, 100)
    y = np.linspace(0, 1, 100)
    X, Y = np.meshgrid(x, y)

    #risk scores
    Z = X * Y

    #green -> orange -> red
    colors = ['#2E7D32', '#D88F1B', '#F79A1E', '#F25A1B',
              '#DC2B1E', '#B8181C', '#930012', '#600000']
    #creates a color map
    cmap = LinearSegmentedColormap.from_list('risk', colors, N=250)

    #displays contour plor
    contour = ax.contourf(X, Y, Z, levels=50, cmap=cmap, alpha=0.8)


    #0.14 threshold line
    risk_threshold_medium = 0.14
    x_vals_medium = np.linspace(0.01, 1, 100)
    y_vals_medium = risk_threshold_medium / x_vals_medium

    # Only show values where consequence stays within plot range
    mask = y_vals_medium <= 1
    ax.plot(x_vals_medium[mask], y_vals_medium[mask], linestyle='--', color='black', linewidth=2, alpha=0.8)

    #0.41 threshold line
    risk_threshold_high = 0.41
    x_vals_high = np.linspace(0.01, 1, 100)
    y_vals_high = risk_threshold_high / x_vals_high

    mask = y_vals_high <= 1
    ax.plot(x_vals_high[mask], y_vals_high[mask], linestyle='--', color='black', linewidth=2, alpha=0.8)

    # Labels and formatting
    ax.set_xlabel('Likelihood', fontsize=12)
    ax.set_ylabel('Consequence', fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title('Low (<0.14), Medium (0.14-0.41), High (>0.41)', fontsize=14)

    #colorbar
    cbar = plt.colorbar(contour, ax=ax, shrink=0.8)
    cbar.set_label('Risk Score', fontsize=10)

    plt.tight_layout()
    plt.show()


def add_risk_points(likelihood_list, consequence_list, labels):
    """Add risk points to the matrix"""
    fig, ax = plt.subplots(figsize=(10, 8))

    # Create the same gradient matrix
    x = np.linspace(0, 1, 100)
    y = np.linspace(0, 1, 100)
    X, Y = np.meshgrid(x, y)
    Z = X * Y

    colors = ['#2E7D32', '#D88F1B', '#F79A1E', '#F25A1B',
              '#DC2B1E', '#B8181C', '#930012', '#600000']
    # creates a color map
    cmap = LinearSegmentedColormap.from_list('risk', colors, N=250)

    # displays contour plot
    contour = ax.contourf(X, Y, Z, levels=50, cmap=cmap, alpha=0.8)

    # Add threshold lines
    # 0.14 threshold line
    risk_threshold_medium = 0.14
    x_vals_medium = np.linspace(0.01, 1, 100)
    y_vals_medium = risk_threshold_medium / x_vals_medium

    # Only show values where consequence stays within plot range
    mask = y_vals_medium <= 1
    ax.plot(x_vals_medium[mask], y_vals_medium[mask], linestyle='--', color='black', linewidth=2, alpha=0.5)

    # 0.41 threshold line
    risk_threshold_high = 0.41
    x_vals_high = np.linspace(0.01, 1, 100)
    y_vals_high = risk_threshold_high / x_vals_high

    mask = y_vals_high <= 1
    ax.plot(x_vals_high[mask], y_vals_high[mask], linestyle='--', color='black', linewidth=2, alpha=0.5)

    # Plot risk points
    offsets = [(0.02, 0.02), (-0.03, 0.02), (0.02, -0.03), (-0.03, -0.03),
               (0.01, 0.03), (0.03, -0.01), (-0.02, 0.02), (0.02, 0.01),
               (-0.01, -0.02), (0.01, -0.01), (-0.02, -0.01), (0.01, 0.01),
               (0.03, 0.03), (-0.03, 0.03), (0.03, -0.03), (-0.02, 0.03),
               (0.01, -0.02), (-0.01, 0.02), (0.03, 0), (-0.02, 0.01),
               (0.01, 0.02), (-0.01, -0.01)]

    point_counts = defaultdict(int)

    for lik, cons, label in zip(likelihood_list, consequence_list, labels):
        ax.scatter(lik, cons, color='black', s=30, edgecolor='black', marker='D')

        key = (round(lik, 3), round(cons, 3))
        count = point_counts[key]

        dx = 0.02
        # If point is near top edge, adjust downward instead of stacking upward
        if cons > 0.95:
            dy = -0.015 * count  # stack downward from top edge
            va = 'top'
        else:
            dy = -0.015 * count  # same stacking as before
            va = 'center'

        ax.text(lik + dx, cons + dy, label, fontsize=8, fontweight='bold',
                ha='left', va=va)

        point_counts[key] += 1

    ax.set_xlabel('Likelihood', fontsize=12)
    ax.set_ylabel('Consequence', fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title('Low (<0.14), Medium (0.14-0.41), High (>0.41)', fontsize=14)

    cbar = plt.colorbar(contour, ax=ax, shrink=0.8)
    cbar.set_label('Risk Score', fontsize=10)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    create_risk_matrix()

    likelihood = [
        0.7, 0.5,
        0.5, 0.5,
        0.5, 0.5, 0.4, 0.3,
        0.5, 0.3,
        0.3, 0.3, 0.5, 0.5, 0.7, 0.5,
        0.3, 0.6, 0.5,
        0.5, 0.5
    ]
    consequence = [
        0.7, 1.0,
        0.9, 0.7,
        0.8, 0.7, 0.8, 0.7,
        0.8, 0.8,
        0.9, 0.9, 0.7, 0.5, 0.7, 0.9,
        0.5, 0.6, 0.9,
        0.7, 0.6, 0.8
    ]

    labels = [
        'TR01', 'TR02',
        'TR03', 'TR04',
        'TR05/TR08', 'TR06', 'TR07', 'TR09–TR11',
        'TR12', 'TR15',
        'TR17–TR19', 'TR20', 'TR21',
        'TR39/TR66', 'TR67–TR68', 'TR24–TR34',
        'TR22/TR36', 'TR37', 'TR64',
        'TR41–TR60',
        'TR61', 'TR62'
    ]

    add_risk_points(likelihood, consequence, labels)

    create_risk_matrix()


    likelihood_2 = [
        0.2, 0.2,
        0.1, 0.1,
        0.3, 0.1, 0.2, 0.1,
        0.1, 0.1,
        0.1, 0.1, 0.3, 0.3,0.3,0.3,
        0.3, 0.2, 0.1, 0.2,
        0.1, 0.2]

    consequence_2 = [
        0.7, 0.8,
        0.5, 0.7,
        0.6, 0.5,0.6, 0.7,
        0.8, 0.6,
        0.9, 0.9, 0.5, 0.5, 0.7, 0.9,
        0.4, 0.6, 0.9, 0.7,
        0.4, 0.8]

    labels_2 = [
        'TR01', 'TR02',
        'TR03', 'TR04',
        'TR05/TR08', 'TR06', 'TR07', 'TR09–TR11',
        'TR12', 'TR15',
        'TR17–TR19', 'TR20', 'TR21',
        'TR39/TR66', 'TR67–TR68', 'TR24–TR34',
        'TR22/TR36', 'TR37', 'TR64',
        'TR41–TR60',
        'TR61', 'TR62'
    ]

    add_risk_points(likelihood_2, consequence_2, labels_2)

