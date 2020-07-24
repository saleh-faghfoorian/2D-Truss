import numpy as np                      # This project is coded to solve each 2D truss problem
import matplotlib.pyplot as plt         # Done By Saleh, Email : saleh.faghfoorian@gmail.com , github.com/saleh-faghfoorian

N         = 0
M         = 0
Mat_types = 0

materials  = []
nodes     = []
elements  = []
ex_forces = []
FS        = []        
STR       = []
disp      = []

class material:                                                    # defining material
    def __init__(self, number, E, S_y):
        self.number = number
        self.E = E
        self.S_y = S_y

class node:                                                        # defining node
    def __init__(self,number,x,y,f_e_x,f_e_y,constraint_x,constraint_y):
        self.number                = number
        self.x                     = x
        self.y                     = y
        self.f_e_x                 = f_e_x
        self.f_e_y                 = f_e_y
        self.constraint_x          = constraint_x
        self.constraint_y          = constraint_y
        self.u                     = 0.0
        self.v                     = 0.0
        self.elements_connected    = []

class element:                                                     # defining element
    def __init__(self,number,node_i,node_j,A,Material):
        self.number   = number
        self.node_i   = node_i
        self.node_j   = node_j
        self.L        = 0.0
        self.A        = A
        self.Material = Material
        self.angle    = 0.0
        self.s        = 0.0
        self.c        = 0.0
        self.delta    = 0.0
        self.f        = 0.0
        self.stress   = 0.0
        self.FS       = 0.0
    

nodes_file = open("./nodes.txt","r+")                              # opening nodes file : properties of nodes
while nodes_file.readline():                                       # finding the number of nodes
    N += 1
nodes_file.close()
nodes_file = open("./nodes.txt","r+")



mat = open("./materials.txt","r+")                                  # opening material file : properties of the materials
while mat.readline():                                              # fining the number of materials used in the elements
    Mat_types += 1
mat.close()
mat = open("./materials.txt","r+")
for i in range(Mat_types):                                         # creating a list of materials
    M_list = mat.readline()
    M_list = M_list.split()
    M_i    = material(M_list[0], M_list[1], M_list[2])
    materials.append(M_i)
mat.close
    

for i in range(N):                                                 # getting nodes from file
    lists  = nodes_file.readline()
    lists  = lists.split()
    node_1 = node(lists[0],lists[1],lists[2],lists[3],lists[4],lists[5],lists[6])
    nodes.append(node_1)
    if int(lists[5]) == 1 :
        nodes[i].u = 0
    if int(lists[6]) == 1 :
        nodes[i].v = 0
nodes_file.close()


def plot_nodes():                                                  # drawing nodes
    x = [float(nodes[i].x) for i in range(N)]
    y = [float(nodes[j].y) for j in range(N)]
    size = 200
    offset = size / 1000
    plt.scatter(x, y, c = 'gray', s = size , zorder = 5)
    for k, location in enumerate(zip(x,y)):
        plt.annotate(k+1 ,(location[0] - offset, location[1] - offset), zorder = 1000)

def draw_elements(point_1, point_2):                               # drawing elements
    x_1 = float(nodes[point_1].x)
    y_1 = float(nodes[point_1].y)
    x_2 = float(nodes[point_2].x)
    y_2 = float(nodes[point_2].y)
    plt.plot([x_1,x_2], [y_1,y_2], color = 'black', linestyle = '-', zorder = 1)

def draw_final_elements(point_1, point_2):                         # drawing elements after deformation
    x_1 = float(nodes[point_1].x) + float(nodes[point_1].u)
    y_1 = float(nodes[point_1].y) + float(nodes[point_1].v)
    x_2 = float(nodes[point_2].x) + float(nodes[point_2].u)
    y_2 = float(nodes[point_2].y) + float(nodes[point_2].u)
    plt.plot([x_1,x_2], [y_1,y_2], color = 'gray', linestyle = '--', zorder = 1)



for i in range(N):                                                 # creating a vector of external forces
    ex_forces.append(float(nodes[i].f_e_x))
    ex_forces.append(float(nodes[i].f_e_y))

elements_file = open("./elements.txt","r+")                        # opening elements file : properties of elements
while elements_file.readline():                                    # finding the number of elements
    M += 1
elements_file.close()
elements_file = open("./elements.txt","r+")

for i in range(M):                                                 # getting elements from file
    list_2            = elements_file.readline()
    list_2            = list_2.split()
    element_1         = element(int(list_2[0]),int(list_2[1]),int(list_2[2]),float(list_2[3]),int(list_2[4]))
    elements.append(element_1)
    elements[i].L     = np.sqrt((float(nodes[elements[i].node_i].x) - float(nodes[elements[i].node_j].x))**2 \
    + (float(nodes[elements[i].node_i].y) - float(nodes[elements[i].node_j].y))**2)
    if (float(nodes[int(list_2[2])].x) - float(nodes[int(list_2[1])].x)) != 0 :
        elements[i].angle = np.arctan((float(nodes[int(list_2[2])].y) - float(nodes[int(list_2[1])].y)) \
        / (float(nodes[int(list_2[2])].x) - float(nodes[int(list_2[1])].x)))
    else :
        elements[i].angle = (np.pi / 2)
    
    elements[i].s     = np.sin(elements[i].angle)
    elements[i].c     = np.cos(elements[i].angle)

    if elements[i].s < 0 :
        elements[i].s = elements[i].s * (-1)
        elements[i].c = elements[i].c * (-1)
elements_file.close()

for i in range(N):                                                 # finding elements are connected to each node
    for j in range(M):                                             
        if int(nodes[i].number) == int(elements[j].node_i) or int(nodes[i].number) == int(elements[j].node_j) :
            nodes[i].elements_connected.append(elements[j].number) 

C = [[0]*(2*N) for i in range(2*N)]                                  


for i in range(N):                                                 # X-Axis Equations
    for j in range(2*N):  
        if i==j :
            c               = [-((elements[r].A * float(materials[int(elements[r].Material)].E) * (elements[r].c)**2 )/(elements[r].L)) for r in nodes[i].elements_connected]
            d               = [-((elements[r].A * float(materials[int(elements[r].Material)].E) * elements[r].c * elements[r].s ) / (elements[r].L)) for r in nodes[i].elements_connected]
            C[2*i][2*i]     = sum(c)
            C[2*i][2*i + 1] = sum(d)
        else :
            for r in nodes[i].elements_connected :
                if elements[r].node_i == i :
                    q = elements[r].node_j
                elif elements[r].node_j == i :
                    q = elements[r].node_i
                C[2*i][2*q]     = ((elements[r].A * float(materials[int(elements[r].Material)].E) * (elements[r].c)**2) / elements[r].L)
                C[2*i][2*q + 1] = ((elements[r].A * float(materials[int(elements[r].Material)].E) * elements[r].c * elements[r].s) / elements[r].L)

for i in range(N):                                                 # Y-Axis Equations
    for j in range(2*N):  
            if i == j :
                e  = [-((elements[r].A * float(materials[int(elements[r].Material)].E) * elements[r].c * elements[r].s )/(elements[r].L)) for r in nodes[i].elements_connected]
                f  = [-((elements[r].A * float(materials[int(elements[r].Material)].E) * (elements[r].s)**2 ) / (elements[r].L)) for r in nodes[i].elements_connected]
                C[2*i + 1][2*i]     = sum(e)
                C[2*i + 1][2*i + 1] = sum(f)
            else :
                for r in nodes[i].elements_connected :
                    if elements[r].node_i == i :
                        p = elements[r].node_j
                    elif elements[r].node_j == i :
                        p = elements[r].node_i
                    C[2*i + 1][2*p]     = ((elements[r].A * float(materials[int(elements[r].Material)].E) * elements[r].c * elements[r].s) / elements[r].L)
                    C[2*i + 1][2*p + 1] = ((elements[r].A * float(materials[int(elements[r].Material)].E) * (elements[r].s)**2) / elements[r].L)

                                                      # list of displacements

for i in range(N):                                                 # calculating new Matrix C
    if int(nodes[i].constraint_x) == 1 :
        C[2*i][2*i] = 1
        for k in range(N):
            if (k != i) and (int(nodes[k].constraint_x) == 1) :
                C[2*i][2*k] = 0
            if (int(nodes[k].constraint_y) == 1) :
                C[2*i][2*k + 1] = 0



    if int(nodes[i].constraint_y) == 1 :
        C[2*i + 1][2*i + 1] = 1
        for j in range(N):
            if (int(nodes[j].constraint_x) == 1) :
                C[2*i + 1][2*j] = 0
            if (j != i) and (int(nodes[j].constraint_y) == 1) :
                C[2*i + 1][2*j + 1] = 0



    if int(nodes[i].constraint_x) == 0 :
        for s in range(N):
            if int(nodes[s].constraint_x) == 1 :
                C[2*i][2*s] = 0
            if int(nodes[s].constraint_y) == 1 :
                C[2*i][2*s + 1] = 0



    if int(nodes[i].constraint_y) == 0 :
        for q in range(N):
            if int(nodes[q].constraint_x) == 1 :
                C[2*i + 1][2*q] = 0
            if int(nodes[q].constraint_y) == 1 :
                C[2*i + 1][2*q + 1] = 0

K = np.linalg.inv(C)                                               # inverse of matrix C

for i in range(N):                                                 # calculating new Matrix ex_forces
    ex_forces[i] = float(ex_forces[i])
    if int(nodes[i].constraint_x) == 1 :
        ex_forces[2*i] = 0
    if int(nodes[i].constraint_y) == 1 :
        ex_forces[2*i + 1] = 0

for i in range(2*N):                                               # multiplication of Matrices K and ex_forces
    d = [(-1)*(K[i][j] * ex_forces[j]) for j in range(2*N)]
    disp.append(sum(d))

for i in range(N):                                                 # determining all of variables
    if int(nodes[i].constraint_x) == 0 :
        nodes[i].u = disp[2*i]
    if int(nodes[i].constraint_y) == 0 :
        nodes[i].v = disp[2*i + 1]
    if int(nodes[i].constraint_x) == 1 :
        nodes[i].f_e_x = disp[2*i]
    if int(nodes[i].constraint_y) == 1 :
        nodes[i].f_e_y = disp[2*i + 1]

for i in range(M):                                                 # calculating forces, stresses, and Factor of Safety
    elements[i].delta = ((float(nodes[int(elements[i].node_j)].u) - float(nodes[int(elements[i].node_i)].u)) \
    * float(elements[i].c)) + ((float(nodes[int(elements[i].node_j)].v) - float(nodes[int(elements[i].node_i)].v)) \
    * float(elements[i].s))

    elements[i].f = float(elements[i].A) * float(materials[int(elements[i].Material)].E) * float(elements[i].delta) / float(elements[i].L)

    elements[i].stress = elements[i].f / float(elements[i].A)

    elements[i].FS     = np.abs(float(materials[int(elements[i].Material)].S_y) / elements[i].stress)


for k in range(M):                                                 # Creating a list for FS and Stresses
    STR.append(np.abs(float(elements[k].stress)))
    if np.abs(float(elements[k].f)) >= 0.000001 :
        FS.append(elements[k].FS)

plot_nodes()                                                       # plotting nodes

for i in range(M):                                                 # plotting elements of initial and final truss
    draw_elements(int(elements[i].node_i), int(elements[i].node_j))
    draw_final_elements(int(elements[i].node_i), int(elements[i].node_j))


outputs = open("./Report.txt","w+")                        # writing Results in the report file
outputs.write("Element Number\t\t   Force(N)\t\tForce-Type\t    Displacement\t     Stress\t\t\t FS\n")

for i in range(M):                                                 # continuing writing Results
    

    outputs.write(str(elements[i].number))                         # writing element number
    outputs.write("\t\t\t")


    if abs(float(elements[i].f)) < 0.000001 :                      # writing force
        outputs.write("      0.00")
    else :
        outputs.write("{0:10.2f}".format(np.abs(float(elements[i].f))))
    outputs.write(" N\t\t")


    if abs(float(elements[i].f)) < 0.000001 :                      # writing force-type
        outputs.write("zero member\t")
    elif float(elements[i].f) < 0 :
        outputs.write("Compression\t")
    elif float(elements[i].f) > 0 :
        outputs.write("Tension\t\t")
    outputs.write("\t")


    if abs(float(elements[i].delta)) < 0.000001 :                  # writing displacement
        outputs.write("0.00 mm")
    else :
        outputs.write("{0:2.2f}".format(np.abs(float(elements[i].delta))))
        outputs.write(" mm")
    outputs.write("\t\t")


    if abs(float(elements[i].stress)) < 0.000001 :                 # writing stress
        outputs.write("      0.00 MPa\t")
    else :
        outputs.write("{0:10.2f}".format(np.abs(float(elements[i].stress))))
        outputs.write(" MPa\t")
    outputs.write("\t\t")


    if abs(float(elements[i].f)) < 0.000001 :                      # writing factor of safety
        outputs.write("----\t\t")
    else :
        outputs.write("{0:3.2f}".format((np.abs(float(elements[i].FS)))))
    outputs.write("\n")
outputs.write("\n\n")
outputs.write("The Maximum Stress is = ")
outputs.write("{0:6.2f}".format((max(STR))))
outputs.write(" MPa\n\n")
outputs.write("Minimum Fator of Safety for the elements of this truss is = ")
outputs.write("{0:3.2f}".format((min(FS))))
outputs.write("\n\n")
outputs.write("Numner of types of materials used in the elements is = ")
outputs.write(str(Mat_types))
outputs.close()


plt.xlabel("X-Axis")                                               
plt.ylabel("Y-Axis")
plt.title("Truss")
plt.savefig("./truss.png")                                 # saving a diagram for truss

plt.clf()                                                          # clearing plot history

x_f = [int(elements[i].number) for i in range(M)]
y_f = [float(np.abs(elements[i].f)) for i in range(M)]
for i in range(M):
    plt.plot([x_f[i], x_f[i]], [0, y_f[i]], color = 'b', linestyle = '-', zorder = 5)
plt.xlabel("Element Number")
plt.ylabel("Element Force (N)")
plt.title("Element Force (N)")
plt.savefig("./Force Diagram.png")                         # saving a diagram for forces of each elements
