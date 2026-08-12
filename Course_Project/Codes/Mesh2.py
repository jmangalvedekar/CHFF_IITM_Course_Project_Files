import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Grid settings
rows, cols = 20, 20
cell_size = 1

# Create a figure with a smaller size (3.5, 3.5)
fig, ax = plt.subplots(figsize=(3.6, 3.6))

# Draw the grid
for i in range(rows + 1):
    ax.plot([0, cols], [i, i], color='black', linewidth=0.5)
for j in range(cols + 1):
    ax.plot([j, j], [0, rows], color='black', linewidth=0.5)

# Fix axes
ax.set_xlim(0, cols)
ax.set_ylim(0, rows)
ax.set_aspect('equal')
ax.axis('off')  # Hide axes

plt.tight_layout()
plt.show()

