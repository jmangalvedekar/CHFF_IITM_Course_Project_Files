import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
n=10
q_ = 100.0
length = 1.0
gamma = 10.0
dx = length/n
x = np.arange(0,length,dx)
dv = dx*dx
q = q_*dx/gamma

T = np.loadtxt("2_2.txt")
Tb = np.sum(T)/n**2
Tw = np.sum(T[0][:])/n
print("Tb = ",Tb," and Tw = ",Tw)
theta = (T-Tb)/(q_*length/gamma)
plt.plot(x,theta[int(n/2),:])
plt.xlabel(r'$x$ along centreline')
plt.ylabel(r'$\theta$')
plt.title(fr'$\theta$ vs $x$ for $n = {n}$')
plt.show()

#plt.pcolormesh(theta, cmap='jet')#, vmin=30, vmax=150)
plt.contourf(x,x,np.flipud(T),levels=20,cmap='coolwarm')
plt.colorbar(label="Temperature")
plt.title(f"Temperature Distribution for n = {n}")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.show()

Nu=q_*length/(gamma*(Tw-Tb))
print("Nusselt number is",Nu)
