import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

#All data
# Graph: Power (W) vs Thrust (mN)
thrust_data = np.array([6, 8, 15, 17, 20, 24])
power_data = np.array([190, 260, 450, 490, 600, 670])
# Graph: Isp (s) vs Thrust (mN)
isp_data = np.array([1750, 2000, 2450, 3000, 3150, 3400])

target_thrust = np.array([1.0, 1.1, 1.2, 1.3, 1.5, 2.0, 2.5, 2.8, 3.0])

confidence = 0.95
# Linear regression for power vs thrust
slope_pt, intercept_pt, r_value_pt, p_value_pt, std_error_pt = stats.linregress(thrust_data, power_data)
print("Power vs Thrust Relationship")
print(f"   Linear fit: Power = {slope_pt:.2f} * Thrust + {intercept_pt:.2f}")
print(f"   R² = {r_value_pt**2:.4f}")
print(f"   Standard error = {std_error_pt:.2f}")

# Polynomial regression: Isp vs Thrust
def poly_fit_with_stats(x, y, degree=2):
    coeffs = np.polyfit(x, y, degree)
    poly = np.poly1d(coeffs)
    residuals = y - poly(x)
    r2 = 1 - np.sum(residuals**2) / np.sum((y - np.mean(y))**2)
    se = np.sqrt(np.mean(residuals**2))
    return coeffs, poly, r2, se, residuals

poly_coeffs, poly_func, r2_poly, se_poly, isp_residuals = poly_fit_with_stats(thrust_data, isp_data)
print("\nIsp vs Thrust Relationship (Polynomial)")
print(f"   Polynomial fit: Isp = {poly_coeffs[0]:.4f} * Thrust² + {poly_coeffs[1]:.4f} * Thrust + {poly_coeffs[2]:.2f}")
print(f"   R² = {r2_poly:.4f}")
print(f"   Standard error = {se_poly:.2f}")

# Extrapolation for target thrusts
print("\nExtrapolation to 1-3 mN Thrust:")
print("   Target | Power  |  Isp Poly" )
print("   Thrust | (W)    |  (s)      ")
print("   -------|--------|-----------")

for t in target_thrust:
    power = slope_pt * t + intercept_pt
    isp_poly = poly_func(t)
    print(f"   {t:.1f}  | {power:.1f} | {isp_poly:.1f}")



# Create confidence bands for linear regression
def linear_confidence(x, x_data, y_data, slope, intercept, std_err, confidence=0.90):
    n = len(x_data)
    t_value = stats.t.ppf((1 + confidence) / 2., n - 2) #critical t_value for a two-tailed test
    y_pred = slope * x + intercept
    mean_x = np.mean(x_data)
    se_fit = std_err * np.sqrt(1/n + (x - mean_x)**2 / np.sum((x_data - mean_x)**2))
    y_lower = y_pred - t_value * se_fit
    y_upper = y_pred + t_value * se_fit
    return y_pred, y_lower, y_upper

# Confidence band for polynomial
def polynomial_confidence(x, poly, residuals, confidence=0.90):
    n = len(residuals)
    dof = n - len(poly.coefficients)
    t_value = stats.t.ppf((1 + confidence) / 2., dof)
    y_pred = poly(x)
    se = np.std(residuals) #standard error
    y_lower = y_pred - t_value * se
    y_upper = y_pred + t_value * se
    return y_pred, y_lower, y_upper


# Plotting
thrust_ext = np.linspace(0,25,200)

# Power vs Thrust (linear)
power_fit, power_lower, power_upper = linear_confidence(
    thrust_ext, thrust_data, power_data, slope_pt, intercept_pt, std_error_pt, confidence
)

# Isp vs Thrust (polynomial)
isp_fit, isp_lower, isp_upper = polynomial_confidence(
    thrust_ext, poly_func, isp_residuals, confidence
)

plt.figure(figsize=(12, 6))

# Power vs Thrust
plt.subplot(1, 2, 1)
plt.fill_between(thrust_ext, power_lower, power_upper, color='blue', alpha=0.2, label='95% CI')
plt.plot(thrust_ext, power_fit, 'b-', label='Linear Fit')
plt.scatter(thrust_data, power_data, color='black', label='Data')
plt.title('Power vs Thrust with 95% Confidence Band')
plt.xlabel('Thrust (mN)')
plt.ylabel('Power (W)')
plt.legend()
plt.grid(True, alpha=0.3)

# Isp vs Thrust
plt.subplot(1, 2, 2)
plt.fill_between(thrust_ext, isp_lower, isp_upper, color='red', alpha=0.2, label='95% CI')
plt.plot(thrust_ext, isp_fit, 'r-', label='Polynomial Fit')
plt.scatter(thrust_data, isp_data, color='black', label='Data')
plt.title('Isp vs Thrust with 95% Confidence Band')
plt.xlabel('Thrust (mN)')
plt.ylabel('Isp (s)')
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()