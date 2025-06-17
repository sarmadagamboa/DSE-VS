import numpy as np
import matplotlib.pyplot as plt

plt.plot([0, 1, 2, 3, 4, 5, 6], [596.9, 635, 656,573.63,595.21,591.51,618.9], color='blue', label='Baseline Dry Mass')
plt.plot([0.5, 1, 2, 3, 4, 5, 6], [(635+596.9)/2, 635, 656, 573.63, 595.21, 591.51, 618.9], color='green', label='Midterm Dry Mass')
plt.plot([2.5, 3, 4, 5, 6], [(656+573.63)/2, 573.63,595.21,591.51,618.9], color='yellow', label='Final Dry Mass')
plt.plot([0, 6], [635*1.25, 635*1.15], color='#cd4631', linestyle='--', label='Upper Margin')
plt.plot([0, 6], [635*0.75, 635*0.85], color='#cd4631', linestyle='-.', label='Lower Margin')
plt.legend()
plt.xlabel('Iteration number')
plt.ylabel('Dry Mass (kg)')
plt.title('Dry Mass Evolution')
plt.grid(True)
plt.tight_layout()
plt.savefig('Final_Report\Budgets\Figures\MarsExplore_dry_mass_evolution.png')
plt.show()


plt.plot([0, 1, 2, 3, 4, 5, 6], [1360, 972, 991, 703, 956.57, 912, 830], color='blue', label='Baseline Power')
plt.plot([0.5, 1, 2, 3, 4, 5, 6], [(972+1360)/2, 972, 991, 703, 956.57, 912, 830], color='green', label='Midterm Power')
plt.plot([2.5, 3, 4, 5, 6], [(991+703)/2, 703, 956.57, 912, 830], color='yellow', label='Final Power')
plt.plot([0, 6], [972*1.25, 972*1.15], color='#cd4631', linestyle='--', label='Upper Margin')
plt.plot([0, 6], [972*0.75, 972*0.85], color='#cd4631', linestyle='-.', label='Lower Margin')
plt.legend()
plt.xlabel('Iteration number')
plt.ylabel('Power (W)')
plt.title('Power Evolution')
plt.grid(True)
plt.tight_layout()
plt.savefig('Final_Report\Budgets\Figures\MarsExplore_power_evolution.png')
plt.show()


plt.plot([0, 1, 2, 3, 4, 5, 6], [(400*2*0.8+88.5+107.2)*1.075, 797.9725, 790.4475, 755.4025, 764.6475, 690.15, 717.74], color='blue', label='Baseline Cost')
plt.plot([0.5, 1, 2, 3, 4, 5, 6], [((400*2*0.8+88.5+107.2)*1.075+797.9725)/2, 797.9725, 790.4475, 755.4025, 764.6475, 690.15, 717.74], color='green', label='Midterm Cost')
plt.plot([2.5, 3, 4, 5, 6], [(790.4475+755.4025)/2, 755.4025, 764.6475, 690.15, 717.74], color='yellow', label='Final Cost')
plt.plot([0, 6], [797.9725*1.25, 797.9725*1.15], color='#cd4631', linestyle='--', label='Upper Margin')
plt.plot([0, 6], [797.9725*0.75, 797.9725*0.85], color='#cd4631', linestyle='-.', label='Lower Margin')
plt.legend()
plt.xlabel('Iteration number')
plt.ylabel('Cost (Million EUR)')
plt.title('Cost Evolution')
plt.grid(True)
plt.tight_layout()
plt.savefig('Final_Report\Budgets\Figures\MarsExplore_cost_evolution.png')
plt.show()


fig, ax = plt.subplots()

# Outer pie: Breakdown of 'First unit'
labelsFirstUnit = ['Payload', 'Structure', 'Thermal', 'EPS', 'TT&C', 'C&DH', 'ADCS', 'Propulsion', '']
valuesFirstUnit = [85, 33.34, 3.64, 42.45, 22.95, 31.15, 33.52, 15.67, 448.08]
outer_colors = ['#cd4631', '#f2a900', '#00bfae', '#00a0ff', '#ff6f00', '#ff0000', '#ff00ff', '#00ff00', 'white']

ax.pie(
    valuesFirstUnit,
    labels=labelsFirstUnit,
    startangle=90,
    radius=1.1,
    wedgeprops=dict(width=0.1, edgecolor='w'),
    colors=outer_colors
)

# Inner pie: Top-level budget categories
labels = [
    'First unit', 'Second unit', 'Launch', 'Operations',
    'Subsystem level margin', 'System level margin'
]
values = [267.72, 160.63, 107.2, 88.5, 26.13619 * 2 * 0.8, 49.94]
inner_colors = ['#1e90ff', '#ff8c00', '#ff0000', '#ff00ff', '#32cd32', '#40e0d0']

inner_pie = ax.pie(
    values,
    labels=None,
    autopct='%1.1f%%',
    pctdistance=0.55,
    startangle=90,
    radius=1.0,
    wedgeprops=dict(width=0.3, edgecolor='w'),
    labeldistance=1.3,
    colors=inner_colors
)

# Add legend with correct label-color association
ax.legend(
    inner_pie[0],  # Use patch objects from inner pie
    labels
)

ax.set_title('MarsExplore Cost Breakdown', fontsize=16, fontweight='bold')
ax.axis('equal')  # Equal aspect ratio ensures the pie is circular
plt.tight_layout()
plt.savefig('Final_Report\Budgets\Figures\MarsExplore_cost_breakdown.png')
plt.show()

