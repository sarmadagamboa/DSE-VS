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

def dry_mass_vs_cost_analysis(print_results=True):
    Dry_Mass = [666, 550, 809, 1031]
    Cost = [283.66, 218.6, 354.6, 495.45]
    dry_mass_cost_results = linear_regression_analysis(Dry_Mass, Cost)
    if print_results:
        print(f"---------------------------------------------------------------- \n Dry Mass vs Cost Analysis \n----------------------------------------------------------------")
        print(f"R²: {dry_mass_cost_results['r2']:.4f}")
        print(f"RSE: {dry_mass_cost_results['RSE']:.4f}")
        print(f"RSE (% of mean y): {dry_mass_cost_results['RSE_percent']:.2f}%")
        print(f"Function: Cost = {dry_mass_cost_results['coef']:.2f} * Dry Mass + {dry_mass_cost_results['intercept']:.2f}")
    return dry_mass_cost_results

def calculate_total_cost(mass, payload_cost, cost_breakdown, trl_payload=20, operational_cost=88.5, learning_curve=80, launch_cost=107.2):
    adjusted_payload_cost = payload_cost * (1 + trl_payload / 100)
    raw_cost = dry_mass_vs_cost_analysis()["coef"] * mass + dry_mass_vs_cost_analysis()["intercept"]
    bus_cost_wo_trl = (raw_cost) * (1 - 40 / 140)
    
    system_cost = np.zeros(len(cost_breakdown))
    for i, item in enumerate(cost_breakdown):
        if i == 0:
            cost = adjusted_payload_cost
        else:
            percentage = item[1]
            margin = item[2]
            cost = (percentage / 100 * bus_cost_wo_trl) * (1 + margin / 100)
        system_cost[i] = cost
    
    total_subsystem_cost = np.sum(system_cost)
    bus_cost = total_subsystem_cost - system_cost[0]
    second_unit_cost = total_subsystem_cost*2* learning_curve / 100 - total_subsystem_cost
    total_cost = total_subsystem_cost * 2 * learning_curve / 100 + launch_cost + operational_cost
    
    return total_cost, bus_cost, adjusted_payload_cost, total_subsystem_cost, second_unit_cost, system_cost


def run_cost_budget_analysis():
    dry_mass = 618.9
    wet_mass = 1060.2
    system_margin = 1.075
    payload_cost = 85
    cost_breakdown = [
        ["Payload", 40.0, 0],
        ["Structure", 18.3, 5],
        ["Thermal", 2.0, 5],
        ["EPS", 23.3, 5],
        ["TT&C", 12.6, 5],
        ["C&DH", 17.1, 5],
        ["ADCS", 18.4, 5],
        ["Propulsion", 8.6, 5]
    ]

    total_cost, bus_cost, adjusted_payload_cost, unit_cost, second_unit_cost, system_cost = calculate_total_cost(
        dry_mass, payload_cost, cost_breakdown=cost_breakdown, launch_cost=107.2
    )
    cost_breakdown2 = [
        ["Payload", 40.0, 0],
        ["Structure", 18.3, 0],
        ["Thermal", 2.0, 0],
        ["EPS", 23.3, 0],
        ["TT&C", 12.6, 0],
        ["C&DH", 17.1, 0],
        ["ADCS", 18.4, 0],
        ["Propulsion", 8.6, 0]]
    total_cost2, bus_cost2, adjusted_payload_cost2, unit_cost2, second_unit_cost2, system_cost2 = calculate_total_cost(
        dry_mass, payload_cost, cost_breakdown=cost_breakdown2, trl_payload=0, launch_cost=107.2
    )
    # Define width for left and right columns
    label_width = 40
    value_width = 12
    print("\n" + "=" * (label_width + value_width))
    print(f"{'COST SUMMARY':^{label_width + value_width}}")
    print("=" * (label_width + value_width))
    print(f"{'Adjusted Payload Cost:':<{label_width}} {adjusted_payload_cost:>{value_width}.2f} M€")
    print(f"{'Bus Cost (excluding payload):':<{label_width}} {bus_cost:>{value_width}.2f} M€")
    print(f"{'Unit Cost (Total Subsystem Cost):':<{label_width}} {unit_cost:>{value_width}.2f} M€")
    print(f"{'Second Unit Cost (w/ learning curve):':<{label_width}} {second_unit_cost:>{value_width}.2f} M€")
    print(f"{'Launch Cost:':<{label_width}} {107.20:>{value_width}.2f} M€")
    print(f"{'Operational Cost:':<{label_width}} {88.5:>{value_width}.2f} M€")
    print(f"{'System Margin Multiplier:':<{label_width}} {system_margin:>{value_width}.2f}x")
    print(f"{'Total Cost (with margin):':<{label_width}} {(total_cost * system_margin):>{value_width}.2f} M€")
    print(f"{'Effective margin:':<{label_width}} {((total_cost/total_cost2 * 100*system_margin)-100):>{value_width}.2f} %")


    print("\n" + "-" * (label_width + value_width))
    print(f"{'SYSTEM COST BREAKDOWN':^{label_width + value_width}}")
    print("-" * (label_width + value_width))
    for item, cost in zip(cost_breakdown, system_cost):
        print(f"{item[0] + ':':<{label_width}} {cost:>{value_width}.2f} M€")
    print("=" * (label_width + value_width))
    return bus_cost, adjusted_payload_cost, unit_cost, second_unit_cost, system_cost

run_cost_budget_analysis()