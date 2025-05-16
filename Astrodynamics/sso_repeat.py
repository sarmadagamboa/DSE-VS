import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.cm import ScalarMappable
import matplotlib.colors as mcolors

def repeat_sso(sol_range=(5, 30), tol=1e-6, max_iter=100, period_bounds_hr=(1.5, 3.0)):

    # Constants
    mu_mars = 42828.376645          # km^3 / s^2
    R_M = 3396.2                    # equatorial radius, km
    J2 = 0.001958705252674          # constant of oblateness, J2
    P_sol = 24.6230*60*60           # solar rotation period, sec
    n_sol = 1.059*10**(-7)          # mean rotation rate of Mars around the sun, rad/s

    results = []

    for K in range(sol_range[0], sol_range[1] + 1):

        lower_Q = int(np.floor((K * P_sol) / (period_bounds_hr[1] * 3600)))
        upper_Q = int(np.ceil((K * P_sol) / (period_bounds_hr[0] * 3600)))

        for Q in range(lower_Q, upper_Q + 1):

            P_omega_true = K * P_sol / Q
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

                if abs(error) < tol and (200 < a - R_M < 350) and (30 % K == 0):
                    results.append({
                        'K': K,
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

    df['Inclination_deg'] = df['Inclination_deg'].round(4)
    df['Altitude_km'] = df['Altitude_km'].round(4)
    df['Group_ID'] = df.groupby(['Inclination_deg', 'Altitude_km']).ngroup()

    unique_orbits_df = df.sort_values(['K', 'Altitude_km']).drop_duplicates(subset=['Inclination_deg', 'Altitude_km'], keep='first').reset_index(drop=True)
    unique_orbits_df.head()

    print(unique_orbits_df.to_string(index=False))

    qk_pairs_used = unique_orbits_df[['Q', 'K']].drop_duplicates()
    qk_list = list(map(tuple, qk_pairs_used.values))
        
    return unique_orbits_df, qk_list

def generate_repeat_curves(QK_list, i_range_deg=(90, 95), tol=1e-6, max_iter=100):
    
    # Constants
    mu_mars = 42828.376645          # km^3 / s^2
    R_M = 3396.2                    # equatorial radius, km
    J2 = 0.001958705252674          # constant of oblateness, J2
    P_sol = 24.6230*60*60           # solar rotation period, sec
    n_sol = 1.059*10**(-7)          # mean rotation rate of Mars around the sun, rad/s

    results = []

    for Q, K in QK_list:
        P_omega_target = K * P_sol / Q
        inclinations = np.linspace(i_range_deg[0], i_range_deg[1], 300)

        for inc_deg in inclinations:
            i_rad = np.radians(inc_deg)

            # Initial guess for a
            a = R_M + 100
            for _ in range(max_iter):
                P_k = 2 * np.pi * np.sqrt(a**3 / mu_mars)
                P_omega = P_k / (1 + (3 * J2 * R_M**2 * (3 - 4 * np.sin(i_rad)**2)) / (2 * a**2))
                error = P_omega - P_omega_target

                if abs(error) < tol and (200 < a - R_M < 350) and (30 % K == 0):
                    results.append({
                        'Q': Q,
                        'K': K,
                        'Inclination': inc_deg,
                        'a_km': a
                    })
                    break

                a *= (P_omega_target / P_omega)

    return pd.DataFrame(results)

def generate_sso_curve(a_range_km):

    # Constants
    mu_mars = 42828.376645          # km^3 / s^2
    R_M = 3396.2                    # equatorial radius, km
    J2 = 0.001958705252674          # constant of oblateness, J2
    P_sol = 24.6230*60*60           # solar rotation period, sec
    n_sol = 1.059*10**(-7)          # mean rotation rate of Mars around the sun, rad/s

    a = np.array(a_range_km)
    cos_i = -2 * n_sol * a**(7/2) / (3 * J2 * R_M**2 * np.sqrt(mu_mars))
    cos_i = np.clip(cos_i, -1, 1)
    i_rad = np.arccos(cos_i)

    return np.degrees(i_rad)

def plot_repeat_curves(repeat_df, sso=True, sso_a_min=3500, sso_a_max=3850):

    fig, ax = plt.subplots(figsize=(12, 8))

    if sso:
        a_vals = np.linspace(sso_a_min, sso_a_max, 300)
        sso_incl = generate_sso_curve(a_vals)
        ax.plot(a_vals - 3396.2, sso_incl, 'r-', linewidth=2, label='Sun-sync')

    R_M = 3396.2
    altitudes = repeat_df['a_km'] - R_M
    norm = mcolors.Normalize(vmin=altitudes.min(), vmax=altitudes.max())
    cmap = colormaps['viridis'].resampled(256)
    cmap = mcolors.LinearSegmentedColormap.from_list("trimmed_viridis", cmap(np.linspace(0.0, 0.7, 256)))

    sorted_rows = sorted(repeat_df.itertuples(), key=lambda r: r.a_km)

    for i, row in enumerate(sorted_rows):
        alt = row.a_km - R_M
        color = cmap(norm(alt))
        plt.plot([alt, alt], [90, 95], linestyle='-', linewidth=1.5, color=color)
        
        y_offset = 95.1 if i % 2 == 0 else 95.0
        plt.text(alt, y_offset, f"({int(row.Q)},{int(row.K)})", fontsize=8, ha='center', va='bottom', color=color)

    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    ax.set_xlabel('Altitude [km]', fontsize=12)
    ax.set_ylabel('Inclination [deg]', fontsize=12)
    plt.xlim(180, 360)
    ax.set_title('Sun-Synchronous Mars Repeat Orbits', fontsize=14)
    ax.grid(True)
    ax.legend(loc='upper left', fontsize=8)
    fig.tight_layout()
    fig.savefig("repeat_orbit_colored_altitude.pdf", format='pdf', bbox_inches='tight')
    plt.show()

repeat_df, qk_pairs = repeat_sso()
repeat_curve_df = generate_repeat_curves(qk_pairs)
repeat_curve_df = repeat_curve_df.drop_duplicates(subset=['Q', 'K'])
plot_repeat_curves(repeat_curve_df)