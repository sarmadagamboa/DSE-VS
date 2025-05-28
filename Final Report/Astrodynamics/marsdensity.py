import pandas as pd
import matplotlib.pyplot as plt

def plot_mars_density(year=2):

    if year == 2:
        file_path = r"Final Report\Astrodynamics\tpdhsy21.txt"
    elif year == 1:
        file_path = r"Final Report\Astrodynamics\tpdhsy11.txt"
    
    with open(file_path, 'r') as f:
        lines = f.readlines()

    rows = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) < 16:
            continue
        try:
            ls = int(parts[0])           
            alt = float(parts[1])        
            rho = float(parts[15])       
        except ValueError:
            continue
        if ls == 270:
            rows.append((alt, rho * 1e-12))  

    df = pd.DataFrame(rows, columns=["Altitude_km", "Density_kg_per_m3"])
    df_max = df.groupby("Altitude_km").max().reset_index().sort_values("Altitude_km")

    plt.figure()
    plt.plot(df_max["Density_kg_per_m3"], df_max["Altitude_km"])
    plt.xlabel("Max Density [kg/m³]")
    plt.ylabel("Altitude [km]")
    if year == 2:
        plt.title("Worst-Case Mars Atmospheric Density vs Altitude\n(Ls = 270°, High Solar, TES Year 2)")
    elif year == 1:
        plt.title("Worst-Case Mars Atmospheric Density vs Altitude\n(Ls = 270°, High Solar, TES Year 1)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    #df_max.to_csv("mars_density_highsolar_ls270.csv", index=False)

plot_mars_density(year=2)