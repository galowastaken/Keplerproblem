import numpy as np
from matplotlib import pyplot as plt

# initialize starting vectors

r1Start = np.array([0, 0, 0])
r2Start = np.array([0, 0, 0])

# create plot variables

ax = plt.axes(projection='3d')
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_zlim(0, 100)

# propagate vectors

r1_x = [r1Start[0]]
r1_y = [r1Start[1]]
r1_z = [r1Start[2]]
r2_x = [r2Start[0]]
r2_y = [r2Start[1]]
r2_z = [r2Start[2]]

for i in range(100):
    r1_x.append(i)
    r1_y.append(np.sin(i))
    r1_z.append(i)

    r2_x.append(i)
    r2_y.append(i)
    r2_z.append(0)

    ax.plot3D(r1_x, r1_y, r1_z, 'green')
    ax.plot3D(r2_x, r2_y, r2_z, 'red')
    plt.pause(0.01)

plt.show()
