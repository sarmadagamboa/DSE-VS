from Budgets.InitialBudgets import payload_vs_dry_mass_analysis, payload_power_vs_total_power_analysis, bus_mass_vs_cost_analysis



def test_payload_vs_dry_mass_analysis():
    result = payload_vs_dry_mass_analysis()
    assert round(result["r2"], 4) == 0.7846
    assert round(result["coef"], 2) == 2.01
    assert round(result["intercept"], 2) == 446.46

def test_payload_power_vs_total_power_analysis():
    result = payload_power_vs_total_power_analysis()
    assert round(result["r2"], 4) == 0.6142
    assert round(result["coef"], 2) == 1.35
    assert round(result["intercept"], 2) == 888.01

def test_bus_mass_vs_cost_analysis():
    result = bus_mass_vs_cost_analysis()
    assert round(result["r2"], 4) == 0.7771
    assert round(result["coef"], 2) == 0.29
    assert round(result["intercept"], 2) == 183.52