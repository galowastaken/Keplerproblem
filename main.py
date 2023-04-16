import numpy as np
import matplotlib.pyplot as plt

# define initial conditions

gravConst = 6.674 * pow(10, -11)

m1 = 2*pow(10, 26)
m2 = 2*pow(10, 26)

R1Start = np.array([3000000, 0, 0])
R2Start = np.array([-3000000, 0, 0])

V1Start = np.array([-7500, 15000, 0])
V2Start = np.array([7500, -15000, 0])

dt = 0.5

# create plot variables

fig = plt.figure()

#ax = plt.axes(projection='3d')
#ax.set_xlim(-1000, 1000)
#ax.set_ylim(-1000, 500000)
#ax.set_zlim(-500, 500)

# propagate vectors
# using verlet

def verlet(R1Start, R2Start, V1Start, V2Start):

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

    ax = fig.add_subplot(1, 2, 1, projection='3d')

    for i in range(500):
        
        if i==0 or i == 1 or i == 50:
            print(R1_old, R2_old)

        A1_old = gravConst * m2 * np.add(R2_old, -R1_old) / pow(np.linalg.norm(np.add(R2_old, -R1_old)), 3)
        A2_old = gravConst * m1 * np.add(-R2_old, R1_old) / pow(np.linalg.norm(np.add(R2_old, -R1_old)), 3)

        R1_new = np.add(np.add(R1_old, V1_old * dt), A1_old * 0.5 * dt * dt)
        R2_new = np.add(np.add(R2_old, V2_old * dt), A2_old * 0.5 * dt * dt)

        A1_new = gravConst * m2 * np.add(R2_new, -R1_new) / pow(np.linalg.norm(np.add(R2_new, -R1_new)), 3)
        A2_new = gravConst * m1 * np.add(-R2_new, R1_new) / pow(np.linalg.norm(np.add(R2_new, -R1_new)), 3)

        V1_new = np.add(V1_old, np.add(A1_old, A1_new) * 0.5 * dt)
        V2_new = np.add(V2_old, np.add(A2_old, A2_new) * 0.5 * dt)

        r1_x.append(R1_new[0])
        r1_y.append(R1_new[1])
        r1_z.append(R1_new[2])

        r2_x.append(R2_new[0])
        r2_y.append(R2_new[1])
        r2_z.append(R2_new[2])

        ax.plot(r1_x, r1_y, r1_z, 'green')
        ax.plot(r2_x, r2_y, r2_z, 'red')
        plt.pause(0.001)

        R1_old = R1_new
        R2_old = R2_new
        V1_old = V1_new
        V2_old = V2_new



    return True

# using runge kutta

def rungeKutta(R1Start, R2Start, V1Start, V2Start):

    r1_x = [R1Start[0]]
    r1_y = [R1Start[1]]
    r1_z = [R1Start[2]]
    r2_x = [R2Start[0]]
    r2_y = [R2Start[1]]
    r2_z = [R2Start[2]]

    R1old = R1Start
    R2old = R2Start
    V1old = V1Start
    V2old = V2Start

    ax = fig.add_subplot(1, 2, 2, projection='3d')

    for i in range(500):

        if i==0 or i == 1 or i == 50:
            print(R1old, R2old)

        F1_old = gravConst * m2 * np.add(R2old, -R1old) / pow(np.linalg.norm(np.add(R2old, -R1old)), 3)
        U1_old = np.array([R1old[0], R1old[1], V1old[0], V1old[1]])
        K1_1_new = np.array([V1old[0], V1old[1], F1_old[0], F1_old[1]])
        K1_2_new = U1_old + 0.5 * K1_1_new
        K1_3_new = U1_old + 0.5 * K1_2_new
        K1_4_new = U1_old + K1_3_new

        F2_old = gravConst * m1 * np.add(-R2old, R1old) / pow(np.linalg.norm(np.add(R2old, -R1old)), 3)
        U2_old = np.array([R2old[0], R2old[1], V2old[0], V2old[1]])
        K2_1_new = np.array([V2old[0], V2old[1], F2_old[0], F2_old[1]])
        K2_2_new = U2_old + 0.5 * K2_1_new
        K2_3_new = U2_old + 0.5 * K2_2_new
        K2_4_new = U2_old + K2_3_new

        U1_new = np.add(U1_old, dt * np.add(np.add(np.add(1 / 6 * K1_1_new, 1 / 3 * K1_2_new), 1 / 3 * K1_3_new), 1 / 6 * K1_4_new))
        U2_new = np.add(U2_old, dt * np.add(np.add(np.add(1 / 6 * K2_1_new, 1 / 3 * K2_2_new), 1 / 3 * K2_3_new), 1 / 6 * K2_4_new))

        r1_x.append(U1_new[0])
        r1_y.append(U1_new[1])
        r1_z.append(R1Start[2])

        r2_x.append(U2_new[0])
        r2_y.append(U2_new[1])
        r2_z.append(R2Start[2])

        ax.plot(r1_x, r1_y, r1_z, 'yellow')
        ax.plot(r2_x, r2_y, r2_z, 'blue')
        plt.pause(0.01)

        R1old = np.array([U1_new[0], U1_new[1], 0])
        R2old = np.array([U2_new[0], U2_new[1], 0])

        V1old = np.array([U1_new[2], U1_new[3], 0])
        V2old = np.array([U2_new[2], U2_new[3], 0])



if verlet(R1Start, R2Start, V1Start, V2Start):

    rungeKutta(R1Start, R2Start, V1Start, V2Start)

plt.show()
