import numpy as np
import pandas as pd

def repeat_sso(sol_range=(5, 30), tol=1e-6, max_iter=100, period_bounds_hr=(1.5, 3.0)):

    # Constants
    mu_mars = 42828.376645          # km^3 / s^2
    R_M = 3396.2                    # equatorial radius, km
    J2 = 0.001958705252674          # constant of oblateness, J2
    P_sol = 24.6230*60*60           # solar rotation period, sec
    n_sol = 1.059*10**(-7)          # mean rotation rate of Mars around the sun, rad/s

    results = []

    for martian_sols in range(sol_range[0], sol_range[1] + 1):

        lower_Q = int(np.floor((martian_sols * P_sol) / (period_bounds_hr[1] * 3600)))
        upper_Q = int(np.ceil((martian_sols * P_sol) / (period_bounds_hr[0] * 3600)))

        for Q in range(lower_Q, upper_Q + 1):

            P_omega_true = martian_sols * P_sol / Q
            a = R_M + 100               # Initial guess
            converged = False

            for _ in range(max_iter):
                cos_i = -2 * n_sol * a**(7/2) / (3 * J2 * R_M**2 * np.sqrt(mu_mars))

                if abs(cos_i) > 1:
                    break

                i_rad = np.arccos(cos_i)

                P_k = 2 * np.pi * np.sqrt(a**3 / mu_mars)

                P_omega = P_k / (1 + 3 * R_M**2 * J2 * (3 - 4 * (np.sin(i_rad))**2) / (2 * a**2))

                error = P_omega - P_omega_true

                if abs(error) < tol and (100 < a - R_M < 400) and (30 % martian_sols == 0):
                    results.append({
                        'Martian_Sols': martian_sols,
                        'Q': Q,
                        'Inclination_deg': np.degrees(i_rad),
                        'Semi_major_axis_km': a,
                        'Altitude_km': a - R_M,
                        'P_omega_hr': P_omega / 3600
                    })
                    
                    converged = True
                    
                    break

                a *= (P_omega_true / P_omega)

    df = pd.DataFrame(results)

    df['Group_ID'] = df.groupby(['Inclination_deg', 'Altitude_km']).ngroup()

    unique_orbits_df = df.sort_values(['Martian_Sols', 'Altitude_km']).drop_duplicates(subset='Group_ID', keep='first')
    unique_orbits_df.reset_index(drop=True, inplace=True)
    unique_orbits_df.head()

    print(unique_orbits_df.to_string(index=False))
        
    return results

repeat_sso()