import numpy as np
from matplotlib import pyplot as plt

# define initial conditions

gravConst = 6.674 * pow(10, -11)

m1 = 100
m2 = 50

R1Start = np.array([0, 0, 0])
R2Start = np.array([10, 10, 0])

V1Start = np.array([10, 0, 0])
V2Start = np.array([0, 10, 0])

dt = 0.1

# create plot variables

ax = plt.axes(projection='3d')
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_zlim(0, 100)

# propagate vectors

r1_x = [R1Start[0]]
r1_y = [R1Start[1]]
r1_z = [R1Start[2]]
r2_x = [R2Start[0]]
r2_y = [R2Start[1]]
r2_z = [R2Start[2]]

R1_old = R1Start
R2_old = R2Start
V1_old = V1Start
V2_old = V2Start

for i in range(100):

    A1_old = gravConst*m2*np.add(R2_old, -R1_old) / pow(np.linalg.norm(np.add(R2_old, -R1_old)), 3)
    A2_old = gravConst*m1*np.add(-R2_old, R1_old) / pow(np.linalg.norm(np.add(R2_old, -R1_old)), 3)

    R1_new = np.add(R1_old, V1_old*dt, A1_old*0.5*dt*dt)
    R2_new = np.add(R2_old, V2_old*dt, A2_old*0.5*dt*dt)

    A1_new = gravConst*m2*np.add(R2_new, -R1_new) / pow(np.linalg.norm(np.add(R2_new, -R1_new)), 3)
    A2_new = gravConst*m1*np.add(-R2_new, R1_new) / pow(np.linalg.norm(np.add(R2_new, -R1_new)), 3)

    V1_new = np.add(V1_old, np.add(A1_old, A1_new)*0.5*dt)
    V2_new = np.add(V2_old, np.add(A2_old, A2_new)*0.5*dt)

    r1_x.append(R1_new[0])
    r1_y.append(R1_new[1])
    r1_z.append(R1_new[2])

    r2_x.append(R2_new[0])
    r2_y.append(R2_new[1])
    r2_z.append(R2_new[2])

    ax.plot3D(r1_x, r1_y, r1_z, 'green')
    ax.plot3D(r2_x, r2_y, r2_z, 'red')
    plt.pause(0.01)

    R1_old = R1_new
    R2_old = R2_new
    V1_old = V1_new
    V2_old = V2_new

plt.show()
