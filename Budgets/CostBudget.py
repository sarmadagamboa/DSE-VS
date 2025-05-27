from InitialBudgets import dry_mass_vs_cost_analysis
import numpy as np
import matplotlib.pyplot as plt

Cost_coef = dry_mass_vs_cost_analysis()["coef"]
Cost_intercept = dry_mass_vs_cost_analysis()["intercept"]

def calculate_total_cost(mass, payload_cost, trl_margin, cost_breakdown, dual_SC, learning_curve=80, launch_cost=100):
    adjusted_payload_cost = payload_cost * (1 + trl_margin / 100)
    raw_cost = Cost_coef * mass + Cost_intercept
    bus_cost = raw_cost * (1 - 40 / 140)
    
    total_subsystem_cost = 0
    for i, item in enumerate(cost_breakdown):
        if i == 0:
            cost = adjusted_payload_cost
        else:
            percentage = item[1]
            margin = item[2]
            cost = (percentage / 100 * bus_cost) * (1 + margin / 100)
        total_subsystem_cost += cost
    if dual_SC:
        total_cost = total_subsystem_cost * 2 * learning_curve / 100 + launch_cost
    else:
        total_cost = total_subsystem_cost + launch_cost
    return round(total_cost, 1)

## dry mass
m_SST_LRI_ACC = 529
m_SST_LRI_CAI = 656
m_QGG = 665
m_DT = 440

## payload cost
C_p_SST_LRI_ACC = 52 # M€
C_p_SST_LRI_CAI = 85 # M€
C_p_QGG = 60 # M€
C_p_DT = 27 # M€

## technology readiness level cost margin payload
TRL_SST_LRI_ACC = 7.5 #%
TRL_SST_LRI_CAI = 15 #%
TRL_QGG = 20 #%
TRL_DT = 5 #%

cost_breakdown = [
    ["Payload", 40.0, 0],
    ["Structure", 18.3, 5],
    ["Thermal", 2.0, 5],
    ["EPS", 23.3, 5],
    ["TT&C", 12.6, 5],
    ["C&DH", 17.1, 5],
    ["ADCS", 18.4, 5],
    ["Propulsion", 8.6, 5]]

print("Cost of SST-LRI-ACC: ", calculate_total_cost(m_SST_LRI_ACC, C_p_SST_LRI_ACC, TRL_SST_LRI_ACC, cost_breakdown, dual_SC=True))
print("Cost of SST-LRI-CAI: ", calculate_total_cost(m_SST_LRI_CAI, C_p_SST_LRI_CAI, TRL_SST_LRI_CAI, cost_breakdown, dual_SC=True))
print("Cost of QGG: ", calculate_total_cost(m_QGG, C_p_QGG, TRL_QGG, cost_breakdown, dual_SC=False))
print("Cost of DT: ", calculate_total_cost(m_DT, C_p_DT, TRL_DT, cost_breakdown, dual_SC=False))