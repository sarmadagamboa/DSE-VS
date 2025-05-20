import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score


def linear_regression_analysis(X_vals, y_vals):
    """Perform linear regression and return model, R², RSE, and coefficients."""
    X = np.array(X_vals).reshape(-1, 1)
    y = np.array(y_vals).reshape(-1, 1)

    model = LinearRegression()
    model.fit(X, y)
    y_pred = model.predict(X)

    # R² score
    r2 = r2_score(y, y_pred)

    # Residual Standard Error (RSE)
    residuals = y - y_pred
    RSS = np.sum(residuals ** 2)
    n = len(y)
    p = 1
    RSE = np.sqrt(RSS / (n - p - 1))
    RSE_percent = (RSE / np.mean(y)) * 100

    return {
        "model": model,
        "r2": r2,
        "RSE": RSE,
        "RSE_percent": RSE_percent,
        "coef": model.coef_[0][0],
        "intercept": model.intercept_[0]
    }


def payload_vs_dry_mass_analysis(print_results=True):
    Payload_Mass = [40, 19.6, 32, 205, 209, 800, 533, 300, 45, 139, 65]
    Dry_Mass = [398, 202, 490, 872, 800, 1690, 2031, 1000, 376.3, 1031, 809]
    payload_dry_mass_results = linear_regression_analysis(Payload_Mass, Dry_Mass)
    if print_results:
        print(f"---------------------------------------------------------------- \n Payload Mass vs Dry Mass Analysis \n----------------------------------------------------------------")
        print(f"R²: {payload_dry_mass_results['r2']:.4f}")
        print(f"RSE: {payload_dry_mass_results['RSE']:.4f}")
        print(f"RSE (% of mean y): {payload_dry_mass_results['RSE_percent']:.2f}%")
        print(f"Function: Dry Mass = {payload_dry_mass_results['coef']:.2f} * Payload Mass + {payload_dry_mass_results['intercept']:.2f}")
    return payload_dry_mass_results

def payload_power_vs_total_power_analysis(print_results=True):
    Payload_Power = [180, 178, 820, 164, 53]
    Total_Power = [919, 1700, 1950, 1000, 750]
    payload_power_results = linear_regression_analysis(Payload_Power, Total_Power)
    if print_results:
        print(f"---------------------------------------------------------------- \n Payload Power vs Total Power Analysis \n----------------------------------------------------------------")
        print(f"R²: {payload_power_results['r2']:.4f}")
        print(f"RSE: {payload_power_results['RSE']:.4f}")
        print(f"RSE (% of mean y): {payload_power_results['RSE_percent']:.2f}%")
        print(f"Function: Total Power = {payload_power_results['coef']:.2f} * Payload Power + {payload_power_results['intercept']:.2f}")
    return payload_power_results


def bus_mass_vs_cost_analysis(print_results=True):
    Bus_Mass = [331.3, 553, 93, 744, 892]
    Cost = [312, 283.66, 218.6, 354.6, 495.45]
    bus_mass_cost_results = linear_regression_analysis(Bus_Mass, Cost)
    if print_results:
        print(f"---------------------------------------------------------------- \n Bus Mass vs Cost Analysis \n----------------------------------------------------------------")
        print(f"R²: {bus_mass_cost_results['r2']:.4f}")
        print(f"RSE: {bus_mass_cost_results['RSE']:.4f}")
        print(f"RSE (% of mean y): {bus_mass_cost_results['RSE_percent']:.2f}%")
        print(f"Function: Cost = {bus_mass_cost_results['coef']:.2f} * Bus Mass + {bus_mass_cost_results['intercept']:.2f}")
    return bus_mass_cost_results

payload_vs_dry_mass_analysis()
payload_power_vs_total_power_analysis()
bus_mass_vs_cost_analysis()