import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
#plt.rcParams['figure.figsize']=10,1
f = 25
m = 80
n = 20
rows, cols = n, n
T_Exp = np.zeros((m, rows, cols))
T_Imp = np.zeros((m, rows, cols))
T_CN = np.zeros((m, rows, cols))
T_Val = np.zeros((m))

for i in range(m):
    j=f*i
    T_Exp[i] = np.loadtxt(f'output_Exp_{j}.txt')
    T_Imp[i] = np.loadtxt(f'output_Imp_{j}.txt')
    T_CN[i] = np.loadtxt(f'output_CN_{j}.txt')
# Read all numbers from the text file as flat list
data = np.loadtxt(f'Validation.txt')  # Split by whitespace
print(data)
# Convert to float (or int) and select every 25th element
#T_val = np.array([np.array(data[i]) for i in range(0, len(data), 25)])  # Index 24 = 25th number (0-based)
for i in range(m):
    T_Val[i] = data[25*i]
#print(matrices.shape)  # (n, rows, cols)
#np.meshgrid(range(n), range(n))

"""fig = plt.figure(figsize=(3.6, 3.6), dpi=80)
img = plt.pcolormesh(T_Exp[19, :, :], cmap=plt.cm.jet, vmin=296, vmax=302)
plt.title('Time = %f' %0.05 )
plt.xlabel('No. of elements along x',fontsize=10)
plt.ylabel('No. of elements along y',fontsize=10)
cbar = plt.colorbar(img)
cbar.set_label('Temperature (K)')"""
#plt.colorbar()
#plt.plot(x,T[50,:,0])
   
"""def update(i):
    plot.set_array(T[i, :, :])
    plt.title('Time step=%d' % (i))
    #plt.plot(x,T[50,:,i])"""

#ani = animation.FuncAnimation(fig, update, frames=int(m), interval = 200, repeat = True)
#plt.colorbar()
#plt.legend('Temperature (K)')
#plt.show()

t=np.arange(0,m*f,f)
#plt.plot(t,T[:,int(n/2),int(n/2)])
fig = plt.figure(figsize=(3.6, 3.6), dpi=80)
plt.plot(t,T_Exp[:,10,10],label='Explicit')
plt.plot(t,T_Imp[:,10,10],label='Implicit')
plt.plot(t,T_CN[:,10,10],label='CN')
plt.plot(t,T_Val,label='COMSOL')
plt.xlabel('Time Steps',fontsize=10)
plt.ylabel('Temperature',fontsize=10)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.title('Plot of T(t) vs t using various methods',fontsize=10)
plt.legend(fontsize=8)
plt.show()

"""fig, axs = plt.subplots(1, 3, figsize=(10.8, 3.6), dpi=80)

# Times and indices you want to plot
time_indices = [0, 19, 79]
times = [0.00, 0.05, 0.20]  # Example times for titles

for ax, idx, t in zip(axs, time_indices, times):
    img = ax.pcolormesh(T_Exp[idx, :, :], cmap=plt.cm.jet, vmin=296, vmax=302)
    ax.set_title(f'Time = {t:.2f} s')
    ax.set_xlabel('')  # Remove individual x-labels
    ax.set_ylabel('')  # Remove individual y-labels

# Shared x and y labels
fig.text(0.55, 0.2, 'No. of elements along x', ha='center', fontsize=12)
fig.text(0.2, 0.55, 'No. of elements along y', va='center', rotation='vertical', fontsize=12)

# Manually create a colorbar axis below the plots
cbar_ax = fig.add_axes([0.3, 0.15, 0.5, 0.03])  # [left, bottom, width, height]
cbar = fig.colorbar(img, cax=cbar_ax, orientation='horizontal')
cbar.set_label('Temperature (K)')

plt.show()"""
