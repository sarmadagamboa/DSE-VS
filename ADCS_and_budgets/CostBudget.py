from InitialBudgets import Cost_coef, Cost_intercept
import numpy as np
import matplotlib.pyplot as plt

## dry bus mass
m_SST_LRI_ACC = 627.2-40 # kg
m_SST_LRI_CAI = 849.7-150
m_QGG = 0
m_DT = 0

## technology readiness level of payload
TRL_SST_LRI_ACC = 9
TRL_SST_LRI_CAI = 6
TRL_QGG = 6
TRL_DT = 9

## learning curve
Learning_Curve = 80 # percent

## cost according to parametric relationship
Cost_SST_LRI_ACC = Cost_coef * m_SST_LRI_ACC + Cost_intercept
Cost_SST_LRI_CAI = Cost_coef * m_SST_LRI_CAI + Cost_intercept
Cost_QGG = Cost_coef * m_QGG + Cost_intercept
Cost_DT = Cost_coef * m_DT + Cost_intercept

## increase cost according to TRL
# Cost_SST_LRI_ACC = Cost_SST_LRI_ACC * (1 + 0.1 * (9 - TRL_SST_LRI_ACC))
# Cost_SST_LRI_CAI = Cost_SST_LRI_CAI * (1 + 0.1 * (9 - TRL_SST_LRI_CAI))
# Cost_QGG = Cost_QGG * (1 + 0.1 * (9 - TRL_QGG))
# Cost_DT = Cost_DT * (1 + 0.1 * (9 - TRL_DT))

## factor in several units, getting full cost
B = 1-np.log(100/Learning_Curve)/np.log(2)
L = 2**B

Cost_SST_LRI_ACC = Cost_SST_LRI_ACC * L
Cost_SST_LRI_CAI = Cost_SST_LRI_CAI * L
Cost_QGG = Cost_QGG
Cost_DT = Cost_DT


## cost of the spacecraft
print(f"Cost of SST LRI ACC: {Cost_SST_LRI_ACC:.2f} M€")
print(f"Cost of SST LRI CAI: {Cost_SST_LRI_CAI:.2f} M€")
print(f"Cost of QGG: {Cost_QGG:.2f} M€")
print(f"Cost of DT: {Cost_DT:.2f} M€")