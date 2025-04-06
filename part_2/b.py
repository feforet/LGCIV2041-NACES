import pandas as pd
import matplotlib.pyplot as plt

data_1mm_3d = pd.read_csv("b/3d_1mm.csv")
data_2mm_3d = pd.read_csv("b/3d_2mm.csv")
data_5mm_3d = pd.read_csv("b/3d_5mm.csv")
data_10mm_3d = pd.read_csv("b/3d_10mm.csv")
data_20mm_3d = pd.read_csv("b/3d_20mm.csv")

data_1mm_2d = pd.read_csv("b/2d_1mm.csv")
data_2mm_2d = pd.read_csv("b/2d_2mm.csv")
data_5mm_2d = pd.read_csv("b/2d_5mm.csv")
data_10mm_2d = pd.read_csv("b/2d_10mm.csv")
data_20mm_2d = pd.read_csv("b/2d_20mm.csv")

#extracting the data for one particular fiber
data_1mm_3d = data_1mm_3d[(data_1mm_3d['Y'] == 40) & (data_1mm_3d['Z'] == 0)].sort_values(by='X')
data_2mm_3d = data_2mm_3d[(data_2mm_3d['Y'] == 40) & (data_2mm_3d['Z'] == 0)].sort_values(by='X')
data_5mm_3d = data_5mm_3d[(data_5mm_3d['Y'] == 40) & (data_5mm_3d['Z'] == 0)].sort_values(by='X')
data_10mm_3d = data_10mm_3d[(data_10mm_3d['Y'] == 40) & (data_10mm_3d['Z'] == 0)].sort_values(by='X')
data_20mm_3d = data_20mm_3d[(data_20mm_3d['Y'] == 40) & (data_20mm_3d['Z'] == 0)].sort_values(by='X')
data_1mm_2d = data_1mm_2d[data_1mm_2d['Y'] == 40].sort_values(by='X')
data_2mm_2d = data_2mm_2d[data_2mm_2d['Y'] == 40].sort_values(by='X')
data_5mm_2d = data_5mm_2d[data_5mm_2d['Y'] == 40].sort_values(by='X')
data_10mm_2d = data_10mm_2d[data_10mm_2d['Y'] == 40].sort_values(by='X')
data_20mm_2d = data_20mm_2d[data_20mm_2d['Y'] == 40].sort_values(by='X')


def b_2d():
    var = "U-U2"
    plt.figure()
    plt.plot(data_1mm_2d['X'], data_1mm_2d[var], label='1 mm')
    plt.plot(data_2mm_2d['X'], data_2mm_2d[var], label='2 mm')
    plt.plot(data_5mm_2d['X'], data_5mm_2d[var], label='5 mm')
    plt.plot(data_10mm_2d['X'], data_10mm_2d[var], label='10 mm')
    plt.plot(data_20mm_2d['X'], data_20mm_2d[var], label='20 mm')
    plt.xlabel('X [mm]')
    plt.ylabel('Vertical displacement [mm]')
    plt.title('2D model')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/b_2d.pdf')
    plt.show()
    
    print("2D model")
    print(f"1 mm: {data_1mm_2d['U-U2'].min()}")
    print(f"2 mm: {data_2mm_2d['U-U2'].min()}")
    print(f"5 mm: {data_5mm_2d['U-U2'].min()}")
    print(f"10 mm: {data_10mm_2d['U-U2'].min()}")
    print(f"20 mm: {data_20mm_2d['U-U2'].min()}")

def b_3d():
    var = "U-U2"
    plt.figure()
    plt.plot(data_1mm_3d['X'], data_1mm_3d[var], label='1 mm')
    plt.plot(data_2mm_3d['X'], data_2mm_3d[var], label='2 mm')
    plt.plot(data_5mm_3d['X'], data_5mm_3d[var], label='5 mm')
    plt.plot(data_10mm_3d['X'], data_10mm_3d[var], label='10 mm')
    plt.plot(data_20mm_3d['X'], data_20mm_3d[var], label='20 mm')
    plt.xlabel('X [mm]')
    plt.ylabel('Vertical displacement [mm]')
    plt.title('3D model')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/b_3d.pdf')
    plt.show()
    
    print("-" * 50)
    print("3D model")
    print(f"1 mm: {data_1mm_3d['U-U2'].min()}")
    print(f"2 mm: {data_2mm_3d['U-U2'].min()}")
    print(f"5 mm: {data_5mm_3d['U-U2'].min()}")
    print(f"10 mm: {data_10mm_3d['U-U2'].min()}")
    print(f"20 mm: {data_20mm_3d['U-U2'].min()}")

b_2d()
b_3d()