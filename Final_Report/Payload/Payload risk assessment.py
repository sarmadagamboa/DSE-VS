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

    #0.14 threshold line
    risk_threshold_medium = 0.14
    x_vals_medium = np.linspace(0.01, 1, 100)
    y_vals_medium = risk_threshold_medium / x_vals_medium

    # Only show values where consequence stays within plot range
    mask = y_vals_medium <= 1
    ax.plot(x_vals_medium[mask], y_vals_medium[mask], linestyle='-', color='red', linewidth=2, alpha=0.5)

    #0.41 threshold line
    risk_threshold_high = 0.41
    x_vals_high = np.linspace(0.01, 1, 100)
    y_vals_high = risk_threshold_high / x_vals_high

    mask = y_vals_high <= 1
    ax.plot(x_vals_high[mask], y_vals_high[mask], linestyle='-', color='red', linewidth=2, alpha=0.8)

    # Labels and formatting
    ax.set_xlabel('Consequence', fontsize=12)
    ax.set_ylabel('Development Difficulty', fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title('Low (<0.14), Medium (0.14-0.41), High (>0.41)', fontsize=14)

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

    # 0.14 threshold line
    risk_threshold_medium = 0.14
    x_vals_medium = np.linspace(0.01, 1, 100)
    y_vals_medium = risk_threshold_medium / x_vals_medium

    mask = y_vals_medium <= 1
    ax.plot(x_vals_medium[mask], y_vals_medium[mask], linestyle='-', color='red', linewidth=2, alpha=0.5)

    # 0.41 threshold line
    risk_threshold_high = 0.41
    x_vals_high = np.linspace(0.01, 1, 100)
    y_vals_high = risk_threshold_high / x_vals_high

    mask = y_vals_high <= 1
    ax.plot(x_vals_high[mask], y_vals_high[mask], linestyle='-', color='red', linewidth=2, alpha=0.8)

    offsets = [(0.02, 0.02), (-0.03, 0.02), (0.02, -0.03), (-0.03, -0.03),
               (0.01, 0.03), (0.03, -0.01), (-0.02, 0.02), (0.02, 0.01),
               (-0.01, -0.02), (0.01, -0.01), (-0.02, -0.01), (0.01, 0.01),
               (0.03, 0.03), (-0.03, 0.03), (0.03, -0.03), (-0.02, 0.03),
               (0.01, -0.02), (-0.01, 0.02), (0.03, 0), (-0.02, 0.01),
               (0.01, 0.02), (-0.01, -0.01)]

    point_counts = defaultdict(int)

    for lik, cons, label in zip(likelihood_list, consequence_list, labels):
        ax.scatter(lik, cons, color='black', s=60, edgecolor='black', marker='D')

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

        ax.text(lik + dx, cons + dy, label, fontsize=10, fontweight='normal',
                ha='left', va=va)

        point_counts[key] += 1

    ax.set_xlabel('Consequence', fontsize=12)
    ax.set_ylabel('Development Difficulty', fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title('Low (<0.14), Medium (0.14-0.41), High (>0.41)', fontsize=14)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    create_risk_matrix()

    consequence = [ 0.8,0.7,0.6,0.5,0.6,0.1]
    development_difficulty = [0.11,0.36,0.27,0.27,0.06,0.86]
    labels = [
        'LRI', 'CAI',
        'Iodine Reference Unit', 'OMIS',
        'CMT', 'MiniCAS']


    add_risk_points(consequence, development_difficulty, labels)
