import numpy as np
import matplotlib.pyplot as plt
n = 20
T = np.loadtxt("3_2_CDS.txt")
x = np.arange(0,3.0,3.0/n)
plt.contourf(x,x,np.flipud(T),levels=30,cmap='coolwarm')
#plt.pcolormesh(T, cmap='jet')#, vmin=30, vmax=150)
plt.colorbar(label=fr'$\phi$')
plt.title(fr'$\phi$ Distribution Across the Domain')
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.show()

plt.plot(x,T[int(n/2),:],marker='o', linestyle='-')
plt.xlabel(r'$y$ along centreline')
plt.ylabel(fr'$\phi$')
plt.title(fr'$\phi$ vs $y$ for $n = {n}$')
plt.show()
