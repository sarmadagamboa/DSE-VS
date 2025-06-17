import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from CostBudget import run_cost_budget_analysis

bus_cost, adjusted_payload_cost, unit_cost, second_unit_cost, system_cost = run_cost_budget_analysis()
operational_cost = 88.5*1.075 
launch_cost = 107.2
development_cost = (np.sum(system_cost)*2*0.8 )*1.075 - 20
Years = np.arange(2025, 2045, 1)
Costs = np.zeros(len(Years))
print("Development Cost: ", development_cost)
# 

previous_cum_cost = 0
for i in range(len(Costs)):
    if Years[i] < 2028:
        Costs[i] = 0
    elif Years[i] < 2042:
        T = (Years[i] - 2027) / (2041 - 2027)
        cum_cost = development_cost * T**4 * (5 - 4*T)
        Costs[i] = cum_cost - previous_cum_cost
        previous_cum_cost = cum_cost
    else:
        Costs[i] = operational_cost / len(Years[Years >= 2042])

    if Years[i] == 2040:
        Costs[i] += launch_cost*1.075

rest = Costs[Years == 2040] - launch_cost*1.075
Costs[Years == 2040] = launch_cost*1.075
previous_cum_cost = 0
for i in range(len(Costs)):
    if Years[i] < 2037:
        Costs[i] += 0
    elif Years[i] < 2042:
        T = (Years[i] - 2036) / (2042 - 2036)
        cum_cost = rest[0] *(10*T**2*(1-T)**2*(T) + T**4 * (5 - 4*T))

        Costs[i] += cum_cost - previous_cum_cost
        previous_cum_cost = cum_cost
    else:
        Costs[i] += 0

Costs[Years == 2041] += 15
Costs[Years == 2039] += 5
print(np.sum(Costs[:np.where(Years == 2041)[0][0]])) # Print total cost until launch
print(np.sum(Costs)) # Print total cost until end of mission
data = {
    'Year': Years,
    'Annual Cost (Millions EUR)': Costs
}

# Create DataFrame
df = pd.DataFrame(data)
df['Cumulative Cost'] = df['Annual Cost (Millions EUR)'].cumsum()

# === Plot 1: Cost Over Time ===
plt.figure(figsize=(10, 5))
plt.plot(df['Year'], df['Annual Cost (Millions EUR)'], linestyle='-', color='#cd4631')
plt.fill_between(df['Year'], df['Annual Cost (Millions EUR)'], color='#cd4631', alpha=0.7) 
plt.title('MarsExplore Annual Cost')
plt.xlabel('Year')
plt.xticks(Years)  # Use only the exact year values
plt.ylabel('Annual Cost (Millions EUR)')
plt.grid(True)
plt.tight_layout()
plt.savefig('Final_Report\Budgets\Figures\MarsExplore_annual_cost.png')
plt.show()

# === Plot 2: Cumulative Cost Over Time ===
plt.figure(figsize=(10, 5))
plt.plot(df['Year'], df['Cumulative Cost'], linestyle='-', color='#cd4631')
plt.fill_between(df['Year'], df['Cumulative Cost'], color='#cd4631', alpha=0.7) 
plt.title('MarsExplore Cumulative Cost')
plt.xlabel('Year')
plt.xticks(Years)  # Use only the exact year values
plt.ylabel('Cumulative Cost (Millions EUR)')
plt.grid(True)
plt.tight_layout()
plt.savefig('Final_Report\Budgets\Figures\MarsExplore_cumulative_cost.png')
plt.show()
