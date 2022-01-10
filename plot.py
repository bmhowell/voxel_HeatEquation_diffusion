"""
Plot Voxel results
-Brian Howell
"""
# %%
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import matplotlib.cm as cm

# import data density
# df = pd.read_csv("voxel5.dat", sep='\t')
# df = pd.read_csv("lastTemp.dat", sep=' ')
df = pd.read_csv("density.dat", sep=' ')
data = np.asarray(df.values)
xPos = data[:, 0]
yPos = data[:, 1]
zPos = data[:, 2]
density = data[:, 3]


sideNode = int(np.round(len(xPos) ** (1 / 3), 0))
x = xPos.reshape((sideNode, sideNode, sideNode))
y = yPos.reshape((sideNode, sideNode, sideNode))
z = zPos.reshape((sideNode, sideNode, sideNode))
d = density.reshape((sideNode, sideNode, sideNode))

print('min: ', d.min().min().min())
print('max: ', d.max().max().max())
# normalize density
dNorm = (d - d.min().min().min())
dNorm = dNorm / dNorm.max().max().max()

# reshape and
# xVox = np.ones((x.shape), dtype=bool)
# yVox = np.ones((y.shape), dtype=bool)
# zVox = np.ones((z.shape), dtype=bool)
xVox = (dNorm > 0.05)
yVox = (dNorm > 0.05)
zVox = (dNorm > 0.05)

voxels = xVox | yVox | zVox
print('voxel: ', voxels.shape)

# set the colors of each object
colors = np.empty(voxels.shape, dtype=object)
ax = plt.figure().add_subplot(projection='3d')
ax.voxels(voxels, facecolors=plt.cm.jet(dNorm), edgecolor='k')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel("Z")
ax.view_init(30, 120)
plt.savefig("initParticleDensity.png")

# import data TEMPERATURE
# df = pd.read_csv("voxel5.dat", sep='\t')
df = pd.read_csv("lastTemp.dat", sep=' ')
data = np.asarray(df.values)
xPos = data[:, 0]
yPos = data[:, 1]
zPos = data[:, 2]
temp = data[:, 3]


sideNode = int(np.round(len(xPos) ** (1 / 3), 0))
x = xPos.reshape((sideNode, sideNode, sideNode))
y = yPos.reshape((sideNode, sideNode, sideNode))
z = zPos.reshape((sideNode, sideNode, sideNode))
t = temp.reshape((sideNode, sideNode, sideNode))

print('min: ', t.min().min().min())
print('max: ', t.max().max().max())
# normalize density
tNorm = (t - t.min().min().min())
tNorm = tNorm / tNorm.max().max().max()

# reshape and
# xVox = np.ones((x.shape), dtype=bool)
# yVox = np.ones((y.shape), dtype=bool)
# zVox = np.ones((z.shape), dtype=bool)
xVox = (tNorm > 0.01)
yVox = (tNorm > 0.01)
zVox = (tNorm > 0.01)

voxels = xVox | yVox | zVox
print('voxel: ', voxels.shape)

# set the colors of each object
colors = np.empty(voxels.shape, dtype=object)
ax = plt.figure().add_subplot(projection='3d')
ax.voxels(voxels, facecolors=plt.cm.hot(dNorm), edgecolor='k')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel("Z")
ax.view_init(30, 120)

# m = plt.cm.ScalarMappable(cmap=cm.hot)
# m.set_array(temp)
# m.set_clim(np.amin(temp), np.amax(temp))
#
# cbar2 = plt.colorbar(m, boundaries=np.linspace(np.amin(temp), np.amax(temp), 100))
# cbar2.set_label("Temperature (K)")

plt.savefig("initParticleTemp.png")






# %%
