import numpy as np
import matplotlib.pyplot as plt
T = np.loadtxt("2_1_b.txt")
plt.pcolormesh(T, cmap='jet')#, vmin=30, vmax=150)
plt.colorbar(label="Temperature")
plt.title("Temperature Distribution Across the Domain")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.show()
