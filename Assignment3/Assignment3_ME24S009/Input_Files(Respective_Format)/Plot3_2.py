import numpy as np
import matplotlib.pyplot as plt
n = 50
T1 = np.loadtxt("3_2_CDS.txt")
T2 = np.loadtxt("3_2_QUICK.txt")
x = np.arange(0,3.0,3.0/n)
#plt.contourf(x,x,np.flipud(T1),levels=30,cmap='coolwarm')
#plt.pcolormesh(T, cmap='jet')#, vmin=30, vmax=150)
#plt.colorbar(label=fr'$\phi$')
#plt.title(fr'$\phi$ Distribution Across the Domain')
#plt.xlabel("X-axis")
#plt.ylabel("Y-axis")
#plt.show()

plt.plot(x,T1[int(n/2),:],marker='o', linestyle='-',label='CDS')
plt.plot(x,T2[int(n/2),:],marker='s', linestyle='-',label='QUICK')
plt.xlabel(r'$y$ along centreline')
plt.ylabel(fr'$\phi$')
plt.title(fr'$\phi$ vs $y$ for $n = {n}$')
plt.legend()
plt.show()
