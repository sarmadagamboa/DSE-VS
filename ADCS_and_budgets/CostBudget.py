from InitialBudgets import Cost_coef, Cost_intercept
import numpy as np
import matplotlib.pyplot as plt

## dry bus mass
m_SST_LRI_ACC = 529
m_SST_LRI_CAI = 656
m_QGG = 665
m_DT = 440

## payload cost
C_p_SST_LRI_ACC = 52 # M€
C_p_SST_LRI_CAI = 85 # M€
C_p_QGG = 60 # M€
C_p_DT = 20 # M€

## technology readiness level cost margin payload
TRL_SST_LRI_ACC = 7.5 #%
TRL_SST_LRI_CAI = 15 #%
TRL_QGG = 20 #%
TRL_DT = 5 #%

## learning curve
Learning_Curve = 75 # percent

## cost of the payload
C_p_SST_LRI_ACC = C_p_SST_LRI_ACC * (1+TRL_SST_LRI_ACC/100)
C_p_SST_LRI_CAI = C_p_SST_LRI_CAI * (1+TRL_SST_LRI_CAI/100)
C_p_QGG = C_p_QGG * (1+TRL_QGG/100)
C_p_DT = C_p_DT * (1+TRL_DT/100)

## cost according to parametric relationship
Cost_SST_LRI_ACC = Cost_coef * m_SST_LRI_ACC + Cost_intercept
Cost_SST_LRI_CAI = Cost_coef * m_SST_LRI_CAI + Cost_intercept
Cost_QGG = Cost_coef * m_QGG + Cost_intercept
Cost_DT = Cost_coef * m_DT + Cost_intercept

## bus cost
Cost_bus_SST_LRI_ACC = Cost_SST_LRI_ACC*(1-40/140)
Cost_bus_SST_LRI_CAI = Cost_SST_LRI_CAI*(1-40/140)
Cost_bus_QGG = Cost_QGG*(1-40/140)
Cost_bus_DT = Cost_DT*(1-40/140)

## cost breakdown for ACC
spacecraft_cost_distribution_SST_ACC = [
    ["Payload", 40.0, 0],
    ["Structure", 18.3, 5],
    ["Thermal", 2.0, 5],
    ["EPS", 23.3, 5],
    ["TT&C", 12.6, 5],
    ["C&DH", 17.1, 5],
    ["ADCS", 18.4, 5],
    ["Propulsion", 8.6, 5],
]

# Calculate the cost distribution for SST LRI ACC
for i in range(len(spacecraft_cost_distribution_SST_ACC)):
    if i == 0:
        spacecraft_cost_distribution_SST_ACC[i].append(C_p_SST_LRI_ACC)
    else:
        spacecraft_cost_distribution_SST_ACC[i].append(
        (spacecraft_cost_distribution_SST_ACC[i][1] / 100 * Cost_bus_SST_LRI_ACC)*(1+spacecraft_cost_distribution_SST_ACC[i][2]/100))

# Calculate the total cost of the spacecraft
Total_SSC_SST_LRI_ACC = 0
for i in range(len(spacecraft_cost_distribution_SST_ACC)):
    Total_SSC_SST_LRI_ACC += spacecraft_cost_distribution_SST_ACC[i][3]
Final_SST_LRI_ACC = Total_SSC_SST_LRI_ACC*2*Learning_Curve/100

# Cost breakdown for CAI
Spacecraft_cost_distribution_SST_CAI = [
    ["Payload", 40.0, 0],
    ["Structure", 18.3, 5],
    ["Thermal", 2.0, 5],
    ["EPS", 23.3, 5],
    ["TT&C", 12.6, 10],
    ["C&DH", 17.1, 5],
    ["ADCS", 18.4, 5],
    ["Propulsion", 8.6, 7.5],
]

# Calculate the cost distribution for SST LRI CAI
for i in range(len(Spacecraft_cost_distribution_SST_CAI)):
    if i == 0:
        Spacecraft_cost_distribution_SST_CAI[i].append(C_p_SST_LRI_CAI)
    else:
        Spacecraft_cost_distribution_SST_CAI[i].append(
        (Spacecraft_cost_distribution_SST_CAI[i][1] / 100 * Cost_bus_SST_LRI_CAI)*(1+Spacecraft_cost_distribution_SST_CAI[i][2]/100))

# Calculate the total cost of the spacecraft
Total_SSC_SST_LRI_CAI = 0
for i in range(len(Spacecraft_cost_distribution_SST_CAI)):
    Total_SSC_SST_LRI_CAI += Spacecraft_cost_distribution_SST_CAI[i][3]
Final_SST_LRI_CAI = Total_SSC_SST_LRI_CAI*2*Learning_Curve/100

## cost breakdown for QGG
spacecraft_cost_distribution_QGG = [
    ["Payload", 40.0, 0],
    ["Structure", 18.3, 5],
    ["Thermal", 2.0, 5],
    ["EPS", 23.3, 5],
    ["TT&C", 12.6, 10],
    ["C&DH", 17.1, 5],
    ["ADCS", 18.4, 5],
    ["Propulsion", 8.6, 7.5],
]

# Calculate the cost distribution for QGG
for i in range(len(spacecraft_cost_distribution_QGG)):
    if i == 0:
        spacecraft_cost_distribution_QGG[i].append(C_p_QGG)
    else:
        spacecraft_cost_distribution_QGG[i].append(
        (spacecraft_cost_distribution_QGG[i][1] / 100 * Cost_bus_QGG)*(1+spacecraft_cost_distribution_QGG[i][2]/100))

# Calculate the total cost of the spacecraft
Total_SSC_QGG = 0
for i in range(len(spacecraft_cost_distribution_QGG)):
    Total_SSC_QGG += spacecraft_cost_distribution_QGG[i][3]

## cost breakdown for DT
spacecraft_cost_distribution_DT = [
    ["Payload", 40.0, 0],
    ["Structure", 18.3, 5],
    ["Thermal", 2.0, 5],
    ["EPS", 23.3, 5],
    ["TT&C", 12.6, 10],
    ["C&DH", 17.1, 5],
    ["ADCS", 18.4, 5],
    ["Propulsion", 8.6, 7.5],
]

# Calculate the cost distribution for DT
for i in range(len(spacecraft_cost_distribution_DT)):
    if i == 0:
        spacecraft_cost_distribution_DT[i].append(C_p_DT)
    else:
        spacecraft_cost_distribution_DT[i].append(
        (spacecraft_cost_distribution_DT[i][1] / 100 * Cost_bus_DT)*(1+spacecraft_cost_distribution_DT[i][2]/100))

# Calculate the total cost of the spacecraft
Total_SSC_DT = 0
for i in range(len(spacecraft_cost_distribution_DT)):
    Total_SSC_DT += spacecraft_cost_distribution_DT[i][3]


## cost of the spacecraft
print(f"Cost of SST LRI ACC: {Final_SST_LRI_ACC:.2f} M€")
print(f"Cost of SST LRI CAI: {Final_SST_LRI_CAI:.2f} M€")
print(f"Cost of QGG: {Cost_QGG:.2f} M€")
print(f"Cost of DT: {Cost_DT:.2f} M€")
