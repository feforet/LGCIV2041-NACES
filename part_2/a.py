import pandas as pd
import matplotlib.pyplot as plt

data3 = pd.read_csv("a/3nodes.csv")
data6 = pd.read_csv("a/6nodes.csv")
data60 = pd.read_csv("a/60nodes.csv")

four_points3 = data3[data3['Y'] == 15].sort_values(by='X')
three_points6 = data6[data6['Y'] == 0].sort_values(by='X')
four_points6 = data6[data6['Y'] == 15].sort_values(by='X')
three_points60 = data60[data60['Y'] == 0].sort_values(by='X')
four_points60 = data60[data60['Y'] == 15].sort_values(by='X')

four_points3['X'] = four_points3['X'] - four_points3['X'].min() - 130
four_points6['X'] = four_points6['X'] - four_points6['X'].min() - 130
three_points6['X'] = three_points6['X'] - three_points6['X'].min() - 130
three_points60['X'] = three_points60['X'] - three_points60['X'].min() - 130
four_points60['X'] = four_points60['X'] - four_points60['X'].min() - 130


def a_2():
    var = 'U-U2'
    plt.plot(four_points3['X'], four_points3[var], label='4 points - 3 elements')
    plt.plot(four_points6['X'], four_points6[var], label='4 points - 6 elements')
    plt.plot(four_points60['X'], four_points60[var], label='4 points - 60 elements')
    plt.xlabel('X [mm]')
    plt.ylabel('Vertical displacement [mm]')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/a_disp.pdf')
    plt.show()
    var = 'MISESMAX'
    plt.plot(four_points3['X'], four_points3[var], label='4 points - 3 elements')
    plt.plot(four_points6['X'], four_points6[var], label='4 points - 6 elements')
    plt.plot(four_points60['X'], four_points60[var], label='4 points - 60 elements')
    plt.xlabel('X [mm]')
    plt.ylabel('Maximum von Mises stress [MPa]')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/a_vonmises.pdf')
    plt.show()
    
    print(f"displacement 3 elements: {four_points3['U-U2'].min()} at x = {four_points3['X'][four_points3['U-U2'].idxmin()]}")
    print(f"displacement 6 elements: {four_points6['U-U2'].min()} at x = {four_points6['X'][four_points6['U-U2'].idxmin()]}")
    print(f"displacement 60 elements: {four_points60['U-U2'].min()} at x = {four_points60['X'][four_points60['U-U2'].idxmin()]}")
    print(f"von mises 3 elements: {four_points3['MISESMAX'].max()} at x = {four_points3['X'][four_points3['MISESMAX'].idxmax()]}")
    print(f"von mises 6 elements: {four_points6['MISESMAX'].max()} at x = {four_points6['X'][four_points6['MISESMAX'].idxmax()]}")
    print(f"von mises 60 elements: {four_points60['MISESMAX'].max()} at x = {four_points60['X'][four_points60['MISESMAX'].idxmax()]}")


def a_3():
    var = 'SM-SM1'
    plt.plot(three_points60['X'], three_points60[var]/1000, label='3 points - 60 elements')
    plt.plot(four_points60['X'], four_points60[var]/1000, label='4 points - 60 elements')
    plt.xlabel('X [mm]')
    plt.ylabel('Bending moment (Nm)')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/a_moments.pdf')
    plt.show()
    var = 'SF-SF2'
    plt.plot(three_points60['X'], three_points60[var]/1000, label='3 points - 60 elements')
    plt.plot(four_points60['X'], four_points60[var]/1000, label='4 points - 60 elements')
    plt.xlabel('X [mm]')
    plt.ylabel('Shear force (kN)')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/a_shear.pdf')
    plt.show()

a_2()
a_3()