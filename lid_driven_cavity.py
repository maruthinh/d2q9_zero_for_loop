import numpy as np
import matplotlib.pyplot as plt
import time

from d2q9 import D2Q9
from bc import BoundaryConditions
from sim_params import SimParams

nx = 64
ny = 64
ghost_x = 2
ghost_y = 2
tot_nx = nx + ghost_x
tot_ny = ny + ghost_y
ndim = 2

#Simulation parameters
Re=100
Mach=0.1
ref_len=1.0
res=ref_len/nx

model = D2Q9()
sim_params = SimParams(mach=Mach, reynolds_num=Re, ref_len=ref_len, res=res)
q = model.q

rho_ini = 1.0
p = rho_ini * model.T0
u_ini = 0.0
v_ini = 0.0

#initialization of arrays required
f = np.zeros((q, tot_nx, tot_ny))
moments = np.zeros((3, tot_nx, tot_ny))
moments[0] = rho_ini
moments[1] = u_ini
moments[2] = v_ini

u0 = np.zeros((tot_nx, tot_ny))
v0 = np.zeros((tot_nx, tot_ny))


f = model.feq(moments)

iter = 0
err = 1.0
start = time.time()

#main iteration loop
while(err>1e-10):
    iter+=1
    model.collision(f, sim_params.beta)
    bc = BoundaryConditions(f)
    model.advection(f)
    bc.bounce_back_left()
    bc.bounce_back_bottom()
    bc.bounce_back_right()
    bc.moving_wall_bc_top(sim_params.ref_vel)

    if iter % 1000 == 0:
        moments=model.moments(f)
        e1 = np.sum(np.sqrt((moments[1] - u0) ** 2 + (moments[2] - v0) ** 2))
        e2 = np.sum(np.sqrt(moments[1] ** 2 + moments[2] ** 2))
        err = e1 / e2
        print('error is = ', err)
    u0 = moments[1]
    v0 = moments[2]

print("total time taken: ", time.time()-start)

# to plot contours of velocity
x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, ny)
X, Y = np.meshgrid(x, y)
VelNorm = np.sqrt(moments[1, 1:-1, 1:-1] * moments[1, 1:-1, 1:-1] 
                  + moments[2, 1:-1, 1:-1] * moments[2, 1:-1, 1:-1])

plt.figure()
# cp = plt.contourf(X, Y, np.transpose(VelNorm), 25, cmap=plt.cm.viridis)
cp = plt.contour(X, Y, np.transpose(VelNorm), 15)
cbar=plt.colorbar(cp)
cbar.set_label('Velocity', fontsize=12)
plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('velocity_cont1.png', bbox_inches='tight')

