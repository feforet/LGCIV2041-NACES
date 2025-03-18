# This code was originally prepared for the course LGCIV1023 (Stabilité des Constructions), on Decembre 3 2019
# Igor Bouckaert, João Pacheco de Almeida, Martin Steinmetz
# Adapted for the course "Nonlinear Response Analysis", Rose, Pavia, Italy, in May 2021
# Adapted for the course LGCIV2041 (Numerical Analysis of Civil Engineering Structures) by João Pacheco de Almeida on February 1 2023
# Application of the displacement method to solve frame structures

"""### 0: INTRODUCTION

This code written for Python / Google Colaboraty / Jupyter notebook aims to go step by step through a numerical solution of a frame structure with the stiffness or displacement method.

To do this, simply run each cell of code one after the other and check what each one does.

#### IMPORTANT

Before starting, it is also necessary to run all the cells signaled with #DEFINITIONS within this section header, which define the different functions and libraries that will be used in the rest of the code.
"""

# DEFINITIONS
# Import the required libraries
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv

# DEFINITIONS
# Display of the element undeformed configuration based on the connectivity (incidence) matrix

def PlotUndeformed(coord, connect) : 

    for i in range(len(connect)) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [coord[1][connect[i][0]] , coord[1][connect[i][1]]] , 
                '#000000' ) 
    
    plt.plot(Coord[0], Coord[1], 'ro')
    
    for i in range(len(coord[1])) : 
        plt.annotate(str(i), (coord[0][i], coord[1][i]))
    
    plt.axis('equal')
    plt.grid()
    plt.show()

# DEFINITIONS
# Display of the deformed configuration of the element

def PlotDeformed(coord, connect, displ, scale) : 
    
    coord_def = coord + displ * scale
    
    for i in range(len(connect)) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [coord[1][connect[i][0]] , coord[1][connect[i][1]]] , 
                '#000000' ) 
        plt.plot([coord_def[0][connect[i][0]] , coord_def[0][connect[i][1]]] , 
                [coord_def[1][connect[i][0]] , coord_def[1][connect[i][1]]] , 
                'r-', linewidth=0.5)
    
    plt.plot(coord[0], coord[1], 'ro')
    
    for i in range(len(coord[1])) : 
        plt.annotate(str(i), (coord[0][i], coord[1][i]))
    
    plt.axis('equal')
    plt.grid()
    plt.show()

# DEFINITIONS
# Display of the shear forces

def PlotShear(coord, connect, Shear) :
    raise NotImplementedError
    
    for i in range(len(connect)) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [coord[1][connect[i][0]] , coord[1][connect[i][1]]] , 
                '#000000' ) 
    
    plt.plot(coord[0], coord[1], 'ro')
    
    for i in range(len(connect)) : 
        plt.plot( [coord[0][connect[i][0]], coord[0][connect[i][1]]], 
                [coord[1][connect[i][0]]+ Shear[i][0], coord[1][connect[i][1]]+ Shear[i][1]] , 
                'g' )
        if not (i == len(connect) or i == 0) : 
            plt.plot([coord[0][connect[i][0]], coord[0][connect[i][0]]], 
                [coord[1][connect[i][0]]+ Shear[i-1][1], coord[1][connect[i][1]]+ Shear[i][0]] , 
                'g' )
        elif i == 0 : 
            plt.plot([coord[0][connect[i][0]], coord[0][connect[i][0]]], 
                [coord[1][connect[i][0]], coord[1][connect[i][1]]+ Shear[i][0]] , 
                'g' )
        if i == len(connect) - 1 : 
            plt.plot([coord[0][connect[i][1]], coord[0][connect[i][1]]], 
                [coord[1][connect[i][1]]+ Shear[i][1], coord[1][connect[i][1]]],  
                'g' )
    plt.title('Shear forces')
    plt.xlabel('x [m]')
    plt.ylabel('Shear force [N]')
    plt.grid()
    plt.show()

def PlotBending(coord, connect, Bending, Shear) :
    raise NotImplementedError

    for i in range(len(connect)) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [coord[1][connect[i][0]] , coord[1][connect[i][1]]] , 
                '#000000' ) 
    
    plt.plot(coord[0], coord[1], 'ro')
    Starting_point = 0
    
    for i in range(len(connect)) : 
        L = abs(coord[0][connect[i][1]] - coord[0][connect[i][0]])
        x = np.linspace(0, L)
        a =  - (Shear[i][1] - Shear[i][0]) / (2*L)
        b =  - Shear[i][0]
        c = Bending[i][0]
        M = a * x ** 2 + b * x + c
        
        plt.plot(x+Starting_point,M,'r')
        
        Starting_point += L
                     
        if not (i == len(connect) or i == 0) : 
            plt.plot([coord[0][connect[i][0]], coord[0][connect[i][0]]], 
                [coord[1][connect[i][0]]+ Bending[i-1][1], coord[1][connect[i][1]]+ Bending[i][0]] , 
                'r' )
        elif i == 0 : 
            plt.plot([coord[0][connect[i][0]], coord[0][connect[i][0]]], 
                [coord[1][connect[i][0]], coord[1][connect[i][1]]+ Bending[i][0]] , 
                'r' )
        if i == len(connect) - 1 : 
            plt.plot([coord[0][connect[i][1]], coord[0][connect[i][1]]], 
                [coord[1][connect[i][1]]+ Bending[i][1], coord[1][connect[i][1]]],  
                'r' )
    plt.title('Bending moments')
    plt.xlabel('x [m]')
    plt.ylabel('Bending moment [Nm]')
    plt.grid()
    plt.show()

def plot_bending(coord, connect, xs, bending):
    for i in range(len(connect)) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [coord[1][connect[i][0]] , coord[1][connect[i][1]]] , 
                '#000000' )
    plt.plot(coord[0], coord[1], 'ro')
    plt.plot(np.ravel(xs, order='C'), np.ravel(bending, order='C'), 'r')
    plt.title('Bending moment along the beam')
    plt.xlabel('x [m]')
    plt.ylabel('Bending moment [Nm]')
    plt.grid()
    plt.show()

def plot_shear(coord, connect, xs, shear):
    for i in range(len(connect)) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [coord[1][connect[i][0]] , coord[1][connect[i][1]]] , 
                '#000000' ) 
    plt.plot(coord[0], coord[1], 'ro')
    plt.plot(np.ravel(xs, order='C'), np.ravel(shear, order='C'), 'g')
    plt.title('Shear force along the beam')
    plt.xlabel('x [m]')
    plt.ylabel('Shear force [N]')
    plt.grid()
    plt.show()

"""### 1: INPUT

All parameters that need to be set by the user are listed in this section with # USER
"""

# Parameters of the problem

# MATERIAL
EA = 4e10 # [N]
EI = 4e9 # [Nm^2]
v = 0.2 # [/]
k = 5/6 # [/]
GA = EA/(2*(1+v)) # [N]
GAc = GA*k # [N]


# USER
# Coordinates of the nodes: x in the first line and y in the second line, in [m].



def beam_mesh(n, L):
    """
    Mesh the beam with n nodes and length L
    Input:
        n: int number of nodes
        L: float length of the beam
    Output:
        Coord: np.array of shape (2,n) containing the coordinates of the nodes
        Connect: np.array of shape (n-1,2) containing the connectivity of the nodes
        Elem_Types : np.array of shape (n-1) containing the type of element (6 DoFs)
    """
    Coord = np.zeros((2,n))
    for i in range(n):
        Coord[0][i] = i*L/(n-1)
        Coord[1][i] = 0
    Connect = np.array([[i, i+1] for i in range(n-1)])
    Elem_Types = np.array([6]*(n-1))
    return Coord, Connect, Elem_Types


    
def calcul(n,L,timoshenko, SelRedInt=False):
    # Strucutre connectivity matrix and coordinates of the nodes
    # Defining the type of Element : 6 or 5 DoFs (5 - to complete in Assignment)
    Coord, Connect, Elem_Types = beam_mesh(n,L)


    #PlotUndeformed(Coord, np.array([[0,0]])) # Display of nodes

    # Total number of degrees of freedom + numbering
    No_Ddl = len(Coord[1])*3  # 3 DoF per node
    Num_Ddl = np.arange(No_Ddl) # Indexing starts at 0 in Python
    print('The structure has ' + str(No_Ddl) + ' degrees of freedom')
            
    #PlotUndeformed(Coord, Connect)

    # Total number of elements
    No_Elem = len(Connect)
    print("The structure is composed of " + str(No_Elem) + " elements.")

    # USER
    # Fixed degrees of freedom
    # Structure:
    #Fixed_DoF = np.array([0,1,7])
    Fixed_DoF = np.array([0,1,No_Ddl-2])

    # Free degrees of freedom
    Free_DoF = np.delete(Num_Ddl, Fixed_DoF)

    print("The free degrees of freedom are:")
    print(Free_DoF)

    # Nodal loads

    # Initialization
    P = np.zeros(No_Ddl)

    # USER
    # in [N] and [m]

    # Structure:
    ## Ici on prend en compte que les ddls libres, pas ceux bloqué par les réactions d'appui
    P_f = np.zeros(len(Free_DoF))
    P_f[(No_Ddl//2)-2] = -40e3 # [N]
    P_f[len(Free_DoF)-2] = -2000e3 # [N]


    # Building other vectors:
    P[Free_DoF] = P_f

    # USER

    # Bridge:
    AE_Elem = np.ones(No_Elem)*EA
    EI_Elem = np.ones(No_Elem)*EI
    GAc_Elem = np.ones(No_Elem)*GAc

    # Variable used to magnify the plot of deformed shape, if needed
    Scale = 1000

    """### 2: COMPUTATIONS

    #### Phase 1: Initialization of vectors and matrices
    """

    # Initialization of the vector of structural displacements
    U = np.zeros(No_Ddl)
    # The displacement vector corresponding to the fixed degrees of freedom is a zero vector
    U_d = U[Fixed_DoF]

    # Initialization of the local and global stiffness matrices of the elements, and of the structural stiffness matrix
    K_str = np.zeros((No_Ddl, No_Ddl))
    k_elem_loc = np.zeros((No_Elem,6,6))  # MODIFIE
    k_elem_glob = np.zeros((No_Elem,6,6))

    # Initialization of the matrices containing : 
    # 1. The length of the elements
    L_Elem = np.zeros(No_Elem)
    # 2. Rotation matrix for each element
    r_C = np.zeros((No_Elem,6,6)) #MODIFIE
    # 3. Assembly matrix
    Assemblage = np.zeros((No_Elem, 6))

    """#### Phase 2: Elements' length, rotation matrices, and assembly matrix
    Loop over the elements to calculate their respective lengths, rotation matrices, and assembly matrix
    """

    for i in range(No_Elem) : 
        
        # 1. Element's length
        L_x = Coord[0][Connect[i][1]] - Coord[0][Connect[i][0]]  # Length in x
        L_y = Coord[1][Connect[i][1]] - Coord[1][Connect[i][0]]  # Length in y 
        L_Elem[i] = np.sqrt(L_x**2 + L_y**2)
        
        # 2. Rotation matrices [[cos sin 0 0],[0 0 cos sin]]
        sin = L_y / L_Elem[i] # Sine of the rotation angle of truss element i
        cos = L_x / L_Elem[i] # Cosine of the rotation angle of truss element i
        r_C[i] = np.array([[cos, sin, 0, 0, 0, 0],
                        [-sin, cos, 0, 0, 0, 0],
                        [0, 0, 1, 0, 0, 0],
                        [0, 0, 0, cos, sin, 0],
                        [0, 0, 0, -sin, cos, 0],
                        [0, 0, 0, 0, 0, 1]])  #MODIFIE
        
        # Auxiliary matrices for the assembly: positioning of local matrices in the global matrix
        Assemblage[i] = np.array([Connect[i][0]*3, 
                                Connect[i][0]*3+1,
                                Connect[i][0]*3+2,
                                Connect[i][1]*3,
                                Connect[i][1]*3+1,
                                Connect[i][1]*3+2])
        Assemblage = Assemblage.astype(int)

    """#### Phase 3: Computation of local stiffness matrices
    Loop through the elements to calculate their respective local stiffness matrices in the local (k_loc) and global (k_glob) reference system, followed by assembly into the structural stiffness matrix 
    NOTE: the matrix product is written '@' in Python
    """

    for elem in range(No_Elem) : 
        
        # Stiffness matrices in the local reference system, 6 DoF
        if not timoshenko and Elem_Types[elem] == 6: 
            k_elem_loc[elem] = np.array([[AE_Elem[elem]/L_Elem[elem], 0, 0, -AE_Elem[elem]/L_Elem[elem], 0, 0],
                                        [0, 12*EI_Elem[elem]/L_Elem[elem]**3,  6*EI_Elem[elem]/L_Elem[elem]**2, 0,  -12*EI_Elem[elem]/L_Elem[elem]**3,  6*EI_Elem[elem]/L_Elem[elem]**2],
                                        [0, 6*EI_Elem[elem]/L_Elem[elem]**2,   4*EI_Elem[elem]/L_Elem[elem],   0,   -6*EI_Elem[elem]/L_Elem[elem]**2,   2*EI_Elem[elem]/L_Elem[elem]],
                                        [-AE_Elem[elem]/L_Elem[elem], 0, 0, AE_Elem[elem]/L_Elem[elem], 0, 0],
                                        [0, -12*EI_Elem[elem]/L_Elem[elem]**3, -6*EI_Elem[elem]/L_Elem[elem]**2, 0,  12*EI_Elem[elem]/L_Elem[elem]**3,   -6*EI_Elem[elem]/L_Elem[elem]**2],
                                        [0, 6*EI_Elem[elem]/L_Elem[elem]**2,   2*EI_Elem[elem]/L_Elem[elem],  0,    -6*EI_Elem[elem]/L_Elem[elem]**2,   4*EI_Elem[elem]/L_Elem[elem]]])  #MODIFIE
        
        # TO COMPLETE
        # Stiffness matrices in the local reference system, 5 DoF
        elif not timoshenko and Elem_Types[elem] == 5: 
            # the idea is to use a local 6x6 matrix and add 0 and 1 to cancel the effect of the DoF not considered
            k_elem_loc[elem] = np.array([[AE_Elem[elem]/L_Elem[elem], 0, 0, -AE_Elem[elem]/L_Elem[elem], 0,0],
                                        [0, 3*EI_Elem[elem]/L_Elem[elem]**3,  3*EI_Elem[elem]/L_Elem[elem]**2, 0,  -3*EI_Elem[elem]/L_Elem[elem]**3,0],
                                            [0, 3*EI_Elem[elem]/L_Elem[elem]**2,   3*EI_Elem[elem]/L_Elem[elem],   0,   -3*EI_Elem[elem]/L_Elem[elem]**2,0],
                                            [-AE_Elem[elem]/L_Elem[elem], 0, 0, AE_Elem[elem]/L_Elem[elem], 0,0],
                                            [0, -3*EI_Elem[elem]/L_Elem[elem]**3, -3*EI_Elem[elem]/L_Elem[elem]**2, 0,  3*EI_Elem[elem]/L_Elem[elem]**3,0],[0,0,0,0,0,1]])  #MODIFIE

        elif timoshenko and not SelRedInt:
            k_elem_loc[elem] = np.array([[AE_Elem[elem]/L_Elem[elem], 0, 0, -AE_Elem[elem]/L_Elem[elem], 0, 0],
                                        [0,GAc_Elem[elem]/L_Elem[elem], GAc_Elem[elem]/2, 0, -GAc_Elem[elem]/L_Elem[elem],  GAc_Elem[elem]/2],
                                        [0,GAc_Elem[elem]/2, EI_Elem[elem]/L_Elem[elem]+L_Elem[elem]*GAc_Elem[elem]/3, 0,   -GAc_Elem[elem]/2,  -EI_Elem[elem]/L_Elem[elem]+L_Elem[elem]*GAc_Elem[elem]/6],
                                        [-AE_Elem[elem]/L_Elem[elem], 0, 0, AE_Elem[elem]/L_Elem[elem], 0, 0],
                                        [0, -GAc_Elem[elem]/L_Elem[elem], -GAc_Elem[elem]/2, 0, GAc_Elem[elem]/L_Elem[elem], -GAc_Elem[elem]/2],
                                        [0,GAc_Elem[elem]/2, -EI_Elem[elem]/L_Elem[elem]+L_Elem[elem]*GAc_Elem[elem]/6, 0, -GAc_Elem[elem]/2,  EI_Elem[elem]/L_Elem[elem]+L_Elem[elem]*GAc_Elem[elem]/3]])#MODIFIE
        elif timoshenko and SelRedInt:
            k_elem_loc[elem] = np.array([[AE_Elem[elem]/L_Elem[elem], 0, 0, -AE_Elem[elem]/L_Elem[elem], 0, 0],
                                        [0,GAc_Elem[elem]/L_Elem[elem], GAc_Elem[elem]/2, 0, -GAc_Elem[elem]/L_Elem[elem],  GAc_Elem[elem]/2],
                                        [0,GAc_Elem[elem]/2, EI_Elem[elem]/L_Elem[elem]+L_Elem[elem]*GAc_Elem[elem]/4, 0,   -GAc_Elem[elem]/2,  -EI_Elem[elem]/L_Elem[elem]+L_Elem[elem]*GAc_Elem[elem]/4],
                                        [-AE_Elem[elem]/L_Elem[elem], 0, 0, AE_Elem[elem]/L_Elem[elem], 0, 0],
                                        [0, -GAc_Elem[elem]/L_Elem[elem], -GAc_Elem[elem]/2, 0, GAc_Elem[elem]/L_Elem[elem], -GAc_Elem[elem]/2],
                                        [0,GAc_Elem[elem]/2, -EI_Elem[elem]/L_Elem[elem]+L_Elem[elem]*GAc_Elem[elem]/4, 0, -GAc_Elem[elem]/2,  EI_Elem[elem]/L_Elem[elem]+L_Elem[elem]*GAc_Elem[elem]/4]])
            
        # Stiffness matrices in the global reference system
        k_elem_glob[elem] = np.transpose(r_C[elem]) @ k_elem_loc[elem] @ r_C[elem]
        
        # Assembly of the global structural stiffness matrix
        for j in range(len(k_elem_glob[elem])) : 
            for k in range(len(k_elem_glob[elem])) : 
                K_str[Assemblage[elem][j]][Assemblage[elem][k]] += k_elem_glob[elem][j][k]

    """#### Phase 4: Partitioning of the stiffness matrix"""

    # Sub-matrix for the free DoFs:
    K_ff = K_str[Free_DoF[:,None], Free_DoF[None,:]]

    #print(K_ff)

    # Sub-matrix for the fixed DoFs: 
    K_dd = K_str[Fixed_DoF[:,None], Fixed_DoF[None,:]]

    # Sub-matrices K_fd et K_df:
    K_fd = K_str[Free_DoF[:,None], Fixed_DoF[None,:]]
    K_df = np.transpose(K_fd)

    """#### Phase 5: Displacement's equation

    Solving the displacement equation to find the value of the free degrees of freedom
    """

    U_f = inv(K_ff) @ (P_f - K_fd @ U_d)

    # Completing the global displacement vector:
    U[Free_DoF] = U_f
    U[Fixed_DoF] = U_d
    #print(U)




    """#### Phase 6: Reactions' equation

    Solving the reactions' equation to find the unknown reactions in the structure
    """

    P_d = K_df @ U_f + K_dd @ U_d

    # Completing the vector of nodal forces:
    P[Fixed_DoF] = P_d

    # Computing the structural resisting forces: 
    P_r = K_str @ U

    """#### Phase 7: Computation of internal forces"""

    # Initialization of vectors u_loc et p_loc

    u_loc = np.zeros((No_Elem, 6))  #MODIFIE
    p_loc = np.zeros((No_Elem, 6))  #MODIFIE

    for i in range(No_Elem) : 
        u_loc[i] = r_C[i] @ U[Assemblage[i]]
        p_loc[i] = k_elem_loc[i] @ u_loc[i]
        
        #print(p_loc[i])
    return U, u_loc, P, P_r, p_loc, L_Elem, Scale, Coord, Connect
        
"""### 3: DISPLAY
Display of the internal forces per element


#Choice of the element to display 

Elem_ID_to_display = 4

Shear = np.zeros((1,2))
Bending = np.zeros((1,2))

Shear[0][0] = p_loc[Elem_ID_to_display][1]
Shear[0][1] = -p_loc[Elem_ID_to_display][4]
Bending[0][0] = p_loc[Elem_ID_to_display][2]
Bending[0][1] = -p_loc[Elem_ID_to_display][5]
N = p_loc[Elem_ID_to_display][3]

print(p_loc[Elem_ID_to_display])
print(Bending)

# Display of axial force
print("Axial force in element {} = {}kN".format(Elem_ID_to_display, np.around(N/1000,2)))

# Display of shear force    
PlotShear(np.array([[0, L_Elem[Elem_ID_to_display]],[0,0]]), np.array([[0,1]]), Shear)
print(Shear)

# Display of bending moment
PlotBending(np.array([[0, L_Elem[Elem_ID_to_display]],[0,0]]), np.array([[0,1]]), Bending, Shear) 
print(Bending)
"""



timoshenko = False
U, u_loc, P, P_r, p_loc, L_Elem, Scale, Coord, Connect = calcul(21,10,timoshenko)
# Print the deformed shape of the structure
Disp = np.zeros((2,len(Coord[0])))
# Pas besoin du deplacement de rotation
Disp[0] = U[np.arange(len(Coord[0]))*3]
Disp[1] = U[np.arange(len(Coord[0]))*3+1]
plotdef = False
if plotdef :
    PlotDeformed(Coord, Connect, Disp, Scale)
    print(U[31])

"""(a) (7.5 points) Derive the exact solution for the transverse displacement field 𝑢𝑦0(𝑥) and for the rotational field 𝜃(𝑥). Then,
compute the response of the beam with the Python script distributed for the first exercise session (“Python script to study
frames.py”), i.e. using Euler-Bernoulli finite elements (FEs). Use a mesh of 2 elements as well as a mesh of 20 elements.
Compare the three results and comment."""

# Calcul de la solution exacte
def exact_solution_EB(x,L,P,EI):
    """Excact solution for a beam of length L subjected to a force P in its middle Euler-Bernoulli
    x : location along the beam [m]
    L : length of the beam [m]
    P : force applied in the middle of the beam [N]
    EI : bending stiffness [Nm^2]
    """
    uy = np.zeros_like(x)
    theta = np.zeros_like(x)
    
    for i in range(len(x)):
        if x[i] <= L/2:
            uy[i] = (-P*x[i]/(48*EI))*(3*L**2 - 4*x[i]**2)
            theta[i] = (-P/(48*EI))*(3*L**2 - 12*x[i]**2)
        else:
            uy[i] = (P/(48*EI))*(L**3 - 9*x[i]*L**2 + 12*L*x[i]**2 - 4*x[i]**3)
            theta[i] = (P/(48*EI))*(-9*L**2 + 24*L*x[i] - 12*x[i]**2)
    return uy, theta

def PlotUy(Coords, Us, exactEB=None, exactT=None, lab=None,save = None):
    if lab is None:
        lab = ['' for i in range(len(Us)//3-1)]
    for i in range(len(Us)):
        plt.plot(Coords[i][0], Us[i][1::3]*1e3, label=lab[i])
    if exactEB is not None:
        plt.plot(exactEB[0], exactEB[1]*1e3, label = 'Exact solution for Euler-Bernoulli')
    if exactT is not None:
        plt.plot(exactT[0], exactT[1]*1e3, label = 'Exact solution for Timoshenko')
    plt.title('Transverse displacement field')
    plt.xlabel('x [m]')
    plt.ylabel('uy(x) [mm]')
    plt.legend()
    plt.grid()
    if save is not None:
        plt.savefig(save,bbox_inches='tight')
    plt.show()

def PlotTheta(Coords, Us, exactEB=None, exactT=None, lab=None,save=None):
    if lab == None:
        lab = ['' for i in range(len(Us)//3-1)]
    for i in range(len(Us)):
        plt.plot(Coords[i][0], Us[i][2::3], label=lab[i])
    if exactEB is not None:
        plt.plot(exactEB[0], exactEB[1], label = 'Exact solution for Euler-Bernoulli')
    if exactT is not None:
        plt.plot(exactT[0], exactT[1], label = 'Exact solution for Timoshenko')
    plt.title('Rotational field')
    plt.xlabel('x [m]')
    plt.ylabel('theta(x) [rad]')
    plt.legend()
    plt.grid()
    if save is not None:
        plt.savefig(save,bbox_inches='tight')
    plt.show()

U2, u_loc2, P2, P_r2, p_loc2, L_Elem2, Scale2, Coord2, Connect2 = calcul(3,10,False)
U20, u_loc20, P20, P_r20, p_loc20, L_Elem20, Scale20, Coord20, Connect20 = calcul(21,10,False)
x = np.linspace(0,10,100)
UEB, thetaEB = exact_solution_EB(x,10,40e3,EI)
plot_a = False
if plot_a:
    # Plot u(x) pour les deux maillages et la solution exacte
    PlotUy([Coord20, Coord2], [U20, U2], exactEB=[x, UEB], lab=['20 elements', '2 elements'], save='Uy_EB.pdf')

    # Plot theta(x) pour les deux maillages et la solution exacte
    PlotTheta([Coord20, Coord2], [U20, U2], exactEB=[x, thetaEB], lab=['20 elements', '2 elements'], save='Theta_EB.pdf')


"""(b) (7.5 points) Implement, in the same Python script, the stiffness matrix corresponding to a Timoshenko finite element,
assuming a linear approximation both for the rotations 𝜃(𝑥) and the transverse displacements 𝑢𝑦0(𝑥). Plot the transverse
displacements and rotations for a mesh of 2, 8, 20, and 200 FEs, together with the exact solution. Compare the response
obtained with the four meshes, and with respect to the exact solution. Discuss whether the FE responses satisfy the
compatibility conditions, namely: is the displacement field continuous within each element and between elements?, does it
satisfy the support conditions?"""

def exact_solution_T(x,L,P,EI,GAc):
    """Déplacement et rotation exacts pour une poutre de longueur L soumise à une force P en son milieu Timoshenko"""
    uy = np.zeros_like(x)
    theta = np.zeros_like(x)
    
    for i in range(len(x)):
        if x[i] <= L/2:
            uy[i] = (P*x[i]/(48*EI))*(-3*L**2 + 4*x[i]**2) - P*x[i]/(2*GAc)
            theta[i] = (-P/(48*EI))*(3*L**2 - 12*x[i]**2)
        else:
            uy[i] = (P/(48*EI))*(L**3 - 9*x[i]*L**2 + 12*L*x[i]**2 - 4*x[i]**3) + P*x[i]/(2*GAc)- P*L/(2*GAc)
            theta[i] = (P/(48*EI))*(-9*L**2 + 24*L*x[i] - 12*x[i]**2)
    return uy, theta

U2T, u_loc2T, P2T, P_r2T, p_loc2T, L_Elem2T, Scale2T, Coord2T, Connect2T = calcul(3,10,True)
U8T, u_loc8T, P8T, P_r8T, p_loc8T, L_Elem8T, Scale8T, Coord8T, Connect8T = calcul(9,10,True)
U20T, u_loc20T, P20T, P_r20T, p_loc20T, L_Elem20T, Scale20T, Coord20T, Connect20T = calcul(21,10,True)
U200T, u_loc200T, P200T, P_r200T, p_loc200T, L_Elem200T, Scale200T, Coord200T, Connect200T = calcul(201,10,True)
x = np.linspace(0,10,100)
UT, thetaT = exact_solution_T(x,10,40e3,EI,GAc)

plot_b = False
if plot_b:
    # Plot u(x) pour les deux maillages et la solution exacte
    PlotUy([Coord2T, Coord8T, Coord20T, Coord200T], [U2T, U8T, U20T, U200T], exactT=[x, UT], lab=['2 elements', '8 elements', '20 elements', '200 elements'], save='Uy_T.pdf')
    #PlotUy([Coord2T, Coord8T, Coord20T, Coord200T], [U2T, U8T, U20T, U200T], lab=['2 elements', '8 elements', '20 elements', '200 elements'])
    # Plot theta(x) pour les deux maillages et la solution exacte
    PlotTheta([Coord2T, Coord8T, Coord20T, Coord200T], [U2T, U8T, U20T, U200T], exactT=[x, thetaT], lab=['2 elements', '8 elements', '20 elements', '200 elements'], save='Theta_T.pdf')


"""(c) (7.5 points) Consider the response obtained in (b) with 8 FEs. Show the evolution of the bending moment and shear
force along the beam. Discuss if the FE response satisfies the equilibrium conditions, namely: locally (within each element),
between elements, and the natural boundary conditions."""

def computeBendingShear_T(u_loc, L_Elem, Connect):
    bending = np.zeros((len(Connect),2))
    shear = np.zeros((len(Connect),2))
    xs = np.zeros((len(Connect),2))
    for i in range(len(Connect)):
        L = L_Elem[i]
        u = u_loc[i]
        bending[i][0] = EI/L * (u[5] - u[2])
        bending[i][1] = EI/L * (u[5] - u[2])
        shear[i][0] = GAc * (-1/L * u[1] - u[2] + 1/L * u[4])
        shear[i][1] = GAc * (-1/L * u[1] + 1/L * U[4] - u[5])
        xs[i][0] = i*L 
        xs[i][1] = (i+1)*L
    return bending, shear, xs

def computeBendingShear_T2(u_loc, L_Elem, Connect):
    bending = np.zeros((len(Connect),2))
    shear = np.zeros((len(Connect),2))
    xs = np.zeros((len(Connect),2))
    for i in range(len(Connect)):
        L = L_Elem[i]
        Ltot = L_Elem.sum()
        print("Ltot", Ltot)
        u = u_loc[i]
        xs[i][0] = i*L 
        xs[i][1] = (i+1)*L
        bending[i][0] = EI/Ltot * (u[5] - u[2])
        bending[i][1] = EI/Ltot * (u[5] - u[2])
        shear[i][0] = GAc * (-1/Ltot * u[1] +  (xs[i][0]/Ltot -1)* u[2] + 1/Ltot * u[4]+ (-xs[i][0]/Ltot)* u[5])
        shear[i][1] = GAc * (-1/Ltot * u[1] +  (xs[i][1]/Ltot -1)* u[2] + 1/Ltot * u[4]+ (-xs[i][1]/Ltot)* u[5])
    return bending, shear, xs

bending8T, shear8T, xs8T = computeBendingShear_T2(u_loc8T, L_Elem8T, Connect8T)
#bending8T, shear8T, xs8T = computeBendingShear_T2(u_loc200T, L_Elem200T, Connect200T)
plot_c = True
if plot_c:
    plot_bending(Coord8T, Connect8T, xs8T, bending8T)
    plot_shear(Coord8T, Connect8T, xs8T, shear8T)

"""(d) (7.5 points) Consider the following alternative beam lengths: L = 2 m, 20 m, and 200 m. For each length, and considering
always a mesh of 200 FEs, plot in the same graph the transverse displacement along the beam for: Timoshenko FEs, Euler-
Bernoulli FEs, and the exact solution. Comment on the results, discussing which numerical results you would “trust” as an
engineer if you did not have access to the exact solution."""

x_2m = np.linspace(0,2,201)
x_20m = np.linspace(0,20,201)
x_200m = np.linspace(0,200,201)
U_2m_T, _, _, _, _, _, _, Coord_2m_T, _ = calcul(201,2,True)
U_2m_EB, _, _, _, _, _, _, Coord_2m_EB, _ = calcul(201,2,False)
U_2m_exact_EB, _ = exact_solution_EB(x_2m,2,40e3,EI)
U_2m_exact_T, _ = exact_solution_T(x_2m,2,40e3,EI,GAc)
U_20m_T, _, _, _, _, _, _, Coord_20m_T, _ = calcul(201,20,True)
U_20m_EB, _, _, _, _, _, _, Coord_20m_EB, _ = calcul(201,20,False)
U_20m_exact_EB, _ = exact_solution_EB(x_20m,20,40e3,EI)
U_20m_exact_T, _ = exact_solution_T(x_20m,20,40e3,EI,GAc)
U_200m_T, _, _, _, _, _, _, Coord_200m_T, _ = calcul(201,200,True)
U_200m_EB, _, _, _, _, _, _, Coord_200m_EB, _ = calcul(201,200,False)
U_200m_exact_EB, _ = exact_solution_EB(x_200m,200,40e3,EI)
U_200m_exact_T, _ = exact_solution_T(x_200m,200,40e3,EI,GAc)

plot_d = False
if plot_d:
    PlotUy([Coord_2m_T, Coord_2m_EB], [U_2m_T, U_2m_EB], exactEB=[x_2m, U_2m_exact_EB], exactT=[x_2m, U_2m_exact_T], lab=['Timoshenko FEs', 'Euler-Bernoulli FEs'], save='Uy_2m.pdf')
    PlotUy([Coord_20m_T, Coord_20m_EB], [U_20m_T, U_20m_EB], exactEB=[x_20m, U_20m_exact_EB], exactT=[x_20m, U_20m_exact_T], lab=['Timoshenko FEs', 'Euler-Bernoulli FEs'],save='Uy_20m.pdf')
    PlotUy([Coord_200m_T, Coord_200m_EB], [U_200m_T, U_200m_EB], exactEB=[x_200m, U_200m_exact_EB], exactT=[x_200m, U_200m_exact_T], lab=['Timoshenko FEs', 'Euler-Bernoulli FEs'],save='Uy_200m.pdf')
    #PlotUy([Coord_200m_T, Coord_200m_EB], [U_200m_T, U_200m_EB], exactEB=[x_200m, U_200m_exact_EB], lab=['Timoshenko', 'Euler-Bernoulli'])

"""(e) (7.5 points) Compute analytically and show the Timoshenko stiffness matrix considering selective reduced integration,
as discussed in the lecture. Implement it in the Python script and plot again the transverse displacement for the same cases
and mesh of question (d). Comment on the results obtained."""
U_2m_T_SRI, _, _, _, _, _, _, Coord_2m_T_SRI, _ = calcul(201,2,True,True)
U_20m_T_SRI, _, _, _, _, _, _, Coord_20m_T_SRI, _ = calcul(201,20,True,True)
U_200m_T_SRI, _, _, _, _, _, _, Coord_200m_T_SRI, _ = calcul(201,200,True,True)

plot_e = False
if plot_e:
    PlotUy([Coord_2m_T, Coord_2m_EB], [U_2m_T_SRI, U_2m_EB], exactEB=[x_2m, U_2m_exact_EB], exactT=[x_2m, U_2m_exact_T], lab=['Timoshenko with selective reduced integration', 'Euler-Bernoulli FEs'], save='Uy_2m_SRI.pdf')
    PlotUy([Coord_20m_T, Coord_20m_EB], [U_20m_T_SRI, U_20m_EB], exactEB=[x_20m, U_20m_exact_EB], exactT=[x_20m, U_20m_exact_T], lab=['Timoshenko with selective reduced integration', 'Euler-Bernoulli FEs'], save='Uy_20m_SRI.pdf')
    PlotUy([Coord_200m_T, Coord_200m_EB], [U_200m_T_SRI, U_200m_EB], exactEB=[x_200m, U_200m_exact_EB], exactT=[x_200m, U_200m_exact_T], lab=['Timoshenko with selective reduced integration', 'Euler-Bernoulli FEs'], save='Uy_200m_SRI.pdf')

"""(f) (22.5 points) Adapt the Python script to compute the geometrically nonlinear response of the beam using a corotational
formulation. Consider Euler-Bernoulli FEs. Show the transverse displacements for a mesh of 2 elements as well as for a
mesh of 20 elements. Compare the results with those of the geometrically linear case, in (a), and explain physically the
differences in the results. PS. You can implement a classical Newton-Raphson method and/or a displacement-control
method, as discussed in the lecture and exercise session, and shown in the shared Python code."""
