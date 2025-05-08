import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt


Gravity_Name = ['GRACE', 'GRACE-FO','GRAIL','CHAMP','GOCE']
Gravity_PM = [40, np.nan, 19.6, 32, 205]
Gravity_TM = [432,600,307,522,1077]
Gravity_DM = [398,569,202,490,872]
Gravity_PP = [75,np.nan,np.nan,50,100]
Gravity_TP = [210,355,763,150,1600]

Gravity_PM2 = [40, 19.6, 32, 205]
Gravity_DM2 = [398,202,490,872]

MClass_Name = ['Solar Orbiter','Euclid','PLATO','ARIEL']
MClass_PM = [209,800,533,300]
MClass_TM = [1800,2100,2134,1300]
MClass_DM = [800,1690,2031,1000]
MClass_PP = [180,178,820,164]
MClass_TP = [919,1700,1950,1000]

Mars_Name = ['Mars Odyssey','Mars Reconnaissance Orbiter','MAVEN']
Mars_PM = [45,139,65]
Mars_TM = [725,2180,2454]
Mars_DM = [376.3,1031,809]
Mars_PP = [53,np.nan,np.nan]
Mars_TP = [750,2000,1150]
Mars_PP2 = [53]
Mars_TP2 = [750]

# Combine Dry and payload mass into a single list for each category
Payload_Mass = Gravity_PM2 + MClass_PM + Mars_PM
Dry_Mass = Gravity_DM2 + MClass_DM + Mars_DM


X = np.array(Payload_Mass).reshape(-1, 1)  # Reshape for sklearn
y = np.array(Dry_Mass).reshape(-1, 1)  # Reshape for sklearn

# Fit linear regression model
model = LinearRegression()
model.fit(X, y)

# Predict values
y_pred = model.predict(X)

# Calculate R²
r2 = r2_score(y, y_pred)

# Calculate Residual Standard Error (RSE)
residuals = y - y_pred
RSS = np.sum(residuals ** 2)
n = len(y)
p = 1  # number of predictors
RSE = np.sqrt(RSS / (n - p - 1))
RSE_percent = (RSE / np.mean(y)) * 100
print(f"RSE (% of mean y): {RSE_percent:.2f}%")

# Output
print(f"R² Payload Mass vs Dry Mass: {r2:.4f}")
print(f"RSE Payload MAss vs Dry Mass: {RSE:.4f}")
Payload_Mass_Estimate = 75
Dry_Mass_Estimate = model.predict([[Payload_Mass_Estimate]])[0][0]
print(f"Estimated Dry Mass for Payload Mass {Payload_Mass_Estimate} kg: {Dry_Mass_Estimate:.2f} kg")
plt.plot(X, y, 'o', label='Original data')
plt.plot(
    X, y_pred, 'r',
    label=fr'$m_{{\mathrm{{Dry}}}} = {model.coef_[0][0]:.2f} \cdot m_{{\mathrm{{Payload}}}} + {model.intercept_[0]:.2f}$'
)

plt.xlabel('Payload Mass (kg)')
plt.ylabel('Dry Mass (kg)')
plt.plot(X, y_pred - RSE, 'g--', label='RSE')
plt.plot(X, y_pred + RSE, 'g--')
plt.title('Payload Mass vs Dry Mass')
plt.legend()
plt.grid()
plt.show()

# Combine Dry and payload mass into a single list for each category
Payload_Power = MClass_PP + Mars_PP2
Total_Power = MClass_TP + Mars_TP2
print(Payload_Power)
print(Total_Power)
X = np.array(Payload_Power).reshape(-1, 1)  # Reshape for sklearn
y = np.array(Total_Power).reshape(-1, 1)  # Reshape for sklearn

# fit linear regression model
model = LinearRegression()
model.fit(X, y)

# Predict values
y_pred = model.predict(X)

# Calculate R²
r2 = r2_score(y, y_pred)

# Calculate Residual Standard Error (RSE)
residuals = y - y_pred
RSS = np.sum(residuals ** 2)
n = len(y)
p = 1  # number of predictors
RSE = np.sqrt(RSS / (n - p - 1))
RSE_percent = (RSE / np.mean(y)) * 100
print(f"RSE (% of mean y) for Power: {RSE_percent:.2f}%")

# Output
print(f"R² Payload Power vs Total Power: {r2:.4f}")
print(f"RSE Payload Power vs Total Power: {RSE:.4f}")
Payload_Power_Estimate = 350
Total_Power_Estimate = model.predict([[Payload_Power_Estimate]])[0][0]
print(f"Estimated Total Power for Payload Power {Payload_Power_Estimate} W: {Total_Power_Estimate:.2f} W")
plt.plot(X, y, 'o', label='Original data')
plt.plot(
    X, y_pred, 'r',
    label=fr'$P_{{\mathrm{{Total}}}} = {model.coef_[0][0]:.2f} \cdot P_{{\mathrm{{Payload}}}} + {model.intercept_[0]:.2f}$'
)
plt.xlabel('Payload Power (W)')
plt.ylabel('Total Power (W)')
plt.plot(X, y_pred - RSE, 'g--', label='RSE')
plt.plot(X, y_pred + RSE, 'g--')
plt.title('Payload Power vs Total Power')
plt.legend()
plt.grid()
plt.show()


Cost = [298,284,247,504,662]
Bus_Mass = [331.3,553,93,744,892]

X = np.array(Bus_Mass).reshape(-1, 1)  # Reshape for sklearn
y = np.array(Cost).reshape(-1, 1)  # Reshape for sklearn

# fit linear regression model
model = LinearRegression()
model.fit(X, y)

# Predict values
y_pred = model.predict(X)

# Calculate R²
r2 = r2_score(y, y_pred)

# Calculate Residual Standard Error (RSE)
residuals = y - y_pred
RSS = np.sum(residuals ** 2)
n = len(y)
p = 1  # number of predictors
RSE = np.sqrt(RSS / (n - p - 1))
RSE_percent = (RSE / np.mean(y)) * 100
print(f"RSE (% of mean y) for Cost: {RSE_percent:.2f}%")

# Output
print(f"R² Bus Mass vs Cost: {r2:.4f}")
print(f"RSE Bus Mass vs Cost: {RSE:.4f}")
Bus_Mass_Estimate = 521.865
Cost_Estimate = model.predict([[Bus_Mass_Estimate]])[0][0]
print(f"Estimated Cost for Bus Mass {Bus_Mass_Estimate} kg: {Cost_Estimate:.2f} M$")
plt.plot(X, y, 'o', label='Original data')
plt.plot(
    X, y_pred, 'r',
    label=fr'$C = {model.coef_[0][0]:.2f} \cdot m_{{\mathrm{{Bus}}}} + {model.intercept_[0]:.2f}$'
)
plt.xlabel('Bus Mass (kg)')
plt.ylabel('Cost (M$)')
plt.plot(X, y_pred - RSE, 'g--', label='RSE')
plt.plot(X, y_pred + RSE, 'g--')
plt.title('Bus Mass vs Cost')
plt.legend()
plt.grid()
plt.show()