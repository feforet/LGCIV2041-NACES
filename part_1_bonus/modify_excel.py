import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


data_eul = pd.read_csv("part_1_bonus/abaqus_eul.csv")
data_tim = pd.read_csv("part_1_bonus/abaqus_tim.csv")


eul_2m = data_eul[data_eul['Y'] == 10].sort_values(by='X')
print(eul_2m)
eul_20m = data_eul[data_eul['Y'] == 0].sort_values(by='X')
eul_200m = data_eul[data_eul['Y'] == -10].sort_values(by='X')

tim_2m = data_tim[data_tim['Y'] == 10].sort_values(by='X')
tim_20m = data_tim[data_tim['Y'] == 0].sort_values(by='X')
tim_20m = tim_20m[~((tim_20m['X'] % 10 == 0) & (tim_20m['UT-UT2'] == 0) & (tim_20m['X'] != 20))]
tim_200m = data_tim[data_tim['Y'] == -10].sort_values(by='X')
tim_200m = tim_200m[~((tim_200m['X'] % 10 == 0) & (tim_200m['UT-UT2'] == 0) & (tim_200m['X'] != 200))]


plt.plot(eul_2m['X'], eul_2m['UT-UT2'], label='Euler-Bernoulli')
plt.plot(tim_2m['X'], tim_2m['UT-UT2'], label='Timoshenko')
plt.xlabel('X')
plt.ylabel('vertical displacement')
plt.title('Vertical displacement for L=2m')
plt.grid()
plt.legend()
plt.show()

plt.plot(eul_20m['X'], eul_20m['UT-UT2'], label='Euler-Bernoulli')
plt.plot(tim_20m['X'], tim_20m['UT-UT2'], label='Timoshenko')
plt.xlabel('X')
plt.ylabel('vertical displacement')
plt.title('Vertical displacement for L=20m')
plt.grid()
plt.legend()
plt.show()

plt.plot(eul_200m['X'], eul_200m['UT-UT2'], label='Euler-Bernoulli')
plt.plot(tim_200m['X'], tim_200m['UT-UT2'], label='Timoshenko')
plt.xlabel('X')
plt.ylabel('vertical displacement')
plt.title('Vertical displacement for L=200m')
plt.grid()
plt.legend()
plt.show()



with pd.ExcelWriter("output_bonus.xlsx") as writer:
    # Prepare data for L=2m
    df_2m = pd.DataFrame({
        'X': eul_2m['X'],
        'Euler-Bernoulli': eul_2m['UT-UT2'],
        'Timoshenko': tim_2m['UT-UT2'],
    })
    df_2m.to_excel(writer, sheet_name='L=2m', index=False)
    
    # Prepare data for L=20m
    df_20m = pd.DataFrame({
        'X': eul_20m['X'],
        'Euler-Bernoulli': eul_20m['UT-UT2'],
        'Timoshenko': tim_20m['UT-UT2'],
    })
    df_20m.to_excel(writer, sheet_name='L=20m', index=False)

    # Prepare data for L=200m
    df_200m = pd.DataFrame({
        'X': eul_200m['X'],
        'Euler-Bernoulli': eul_200m['UT-UT2'],
        'Timoshenko': tim_200m['UT-UT2'],
    })
    df_200m.to_excel(writer, sheet_name='L=200m', index=False)
