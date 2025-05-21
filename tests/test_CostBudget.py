from Budgets.CostBudget import calculate_total_cost  # Replace with actual file name (without .py)
from Budgets.InitialBudgets import bus_mass_vs_cost_analysis

Cost_coef = bus_mass_vs_cost_analysis()["coef"]
Cost_intercept = bus_mass_vs_cost_analysis()["intercept"]

def test_cost_tool():
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
    
    result = calculate_total_cost(529, 52, 7.5, cost_breakdown, dual_SC=True)
    
    # This expected value must match your actual printed result from your script
    expected = 561.5 
    assert result == expected