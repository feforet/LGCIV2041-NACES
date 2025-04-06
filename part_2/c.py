import matplotlib.pyplot as plt
import pandas as pd

data_5mm = pd.read_csv("part_2/c/c_5mm.csv")
data_10mm = pd.read_csv("part_2/c/c_10mm.csv")
data_20mm = pd.read_csv("part_2/c/c_20mm.csv")

data_5mm["U-U2"] = -data_5mm["U-U2"]
data_10mm["U-U2"] = -data_10mm["U-U2"]
data_20mm["U-U2"] = -data_20mm["U-U2"]

funct = lambda x: int(float(x[-5:])) * 30
data_5mm['Force'] = data_5mm['Frame'].apply(funct)
data_10mm['Force'] = data_10mm['Frame'].apply(funct)
data_20mm['Force'] = data_20mm['Frame'].apply(funct)
data_5mm = data_5mm.sort_values(by='Force')
data_10mm = data_10mm.sort_values(by='Force')
data_20mm = data_20mm.sort_values(by='Force')
non_zero_ac_yield_5mm = data_5mm[data_5mm['AC YIELD'] != 0].iloc[0]
non_zero_ac_yield_10mm = data_10mm[(data_10mm['AC YIELD'] != 0) & (data_10mm['Node Label'] == 207)].sort_values(by='Force').iloc[0]
non_zero_ac_yield_20mm = data_20mm[data_20mm['AC YIELD'] != 0].iloc[0]
node_label_5mm = non_zero_ac_yield_5mm['Node Label']
node_label_10mm = non_zero_ac_yield_10mm['Node Label']
node_label_20mm = non_zero_ac_yield_20mm['Node Label']
force_5mm = non_zero_ac_yield_5mm['Force']
force_10mm = non_zero_ac_yield_10mm['Force']
force_20mm = non_zero_ac_yield_20mm['Force']
x_5mm = non_zero_ac_yield_5mm['X']
x_10mm = non_zero_ac_yield_10mm['X']
x_20mm = non_zero_ac_yield_20mm['X']
y_5mm = non_zero_ac_yield_5mm['Y']
y_10mm = non_zero_ac_yield_10mm['Y']
y_20mm = non_zero_ac_yield_20mm['Y']
z_5mm = non_zero_ac_yield_5mm['Z']
z_10mm = non_zero_ac_yield_10mm['Z']
z_20mm = non_zero_ac_yield_20mm['Z']

print("5mm")
print(f"The plastification occurs at force: {force_5mm} in node {node_label_5mm} at coordinates ({x_5mm}, {y_5mm}, {z_5mm})")
print(f"The corresponding displacement is: {non_zero_ac_yield_5mm['U-U2']} mm")
print("-" * 50)
print("10mm")
print(f"The plastification occurs at force: {force_10mm} in node {node_label_10mm} at coordinates ({x_10mm}, {y_10mm}, {z_10mm})")
print(f"The corresponding displacement is: {non_zero_ac_yield_10mm['U-U2']} mm")
print("-" * 50)
print("20mm")
print(f"The plastification occurs at force: {force_20mm} in node {node_label_20mm} at coordinates ({x_20mm}, {y_20mm}, {z_20mm})")
print(f"The corresponding displacement is: {non_zero_ac_yield_20mm['U-U2']} mm")

print("-" * 50)
print(f"The maximum force pour 5mm is: {data_5mm['Force'].max()}")


node_5mm = data_5mm[data_5mm['Node Label'] == node_label_5mm]
node_10mm = data_10mm[data_10mm['Node Label'] == node_label_10mm]
node_20mm = data_20mm[data_20mm['Node Label'] == node_label_20mm]

plt.plot(node_5mm['U-U2'], node_5mm['Force'], label='5 mm')
plt.plot(node_10mm['U-U2'], node_10mm['Force'], label='10 mm')
plt.plot(node_20mm['U-U2'], node_20mm['Force'], label='20 mm')
plt.scatter(non_zero_ac_yield_5mm['U-U2'], force_5mm, color='blue', label='Plastification 5 mm')
plt.scatter(non_zero_ac_yield_10mm['U-U2'], force_10mm, color='red', label='Plastification 10 mm')
plt.scatter(non_zero_ac_yield_20mm['U-U2'], force_20mm, color='green', label='Plastification 20 mm')
plt.xlabel('Vertical displacement [mm]')
plt.ylabel('Force [N]')
plt.legend()
plt.grid()
plt.title('Force-displacement curve considering perfect plasticity')
plt.tight_layout()
plt.savefig('part_2/plots/c_f-u.pdf')
plt.show()