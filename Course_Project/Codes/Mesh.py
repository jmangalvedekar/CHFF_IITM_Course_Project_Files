import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Grid settings
rows, cols = 10, 10
cell_size = 1

fig, ax = plt.subplots(figsize=(3.6, 3.6))  # Reduced image size

# Highlight center 2x2 region
highlight_start = (rows - 2) // 2  # 4
highlight_end = highlight_start + 2  # 6

# Draw highlighted rectangles first (orange 2x2 region)
for i in range(highlight_start, highlight_end):
    for j in range(highlight_start, highlight_end):
        rect = patches.Rectangle((j, i), 1, 1, facecolor='orange', edgecolor='orange')
        ax.add_patch(rect)

# Draw the border (blue color)
border_rect = patches.Rectangle((0, 0), cols, rows, linewidth=2, edgecolor='blue', facecolor='none')
ax.add_patch(border_rect)

# Measurement line (horizontal)
x1, x2 = highlight_start, highlight_end
y_line = highlight_start - 0.5  # Just below the square
ax.annotate('', xy=(x2, y_line), xytext=(x1, y_line),
            arrowprops=dict(arrowstyle='<->', color='blue'))
ax.text((x1 + x2) / 2, y_line - 0.8, r'$x_2 - x_1 = 5\Delta x$', ha='center', va='bottom', fontsize=10, color='black')

# Measurement line (vertical)
y1, y2 = highlight_start, highlight_end
x_line = highlight_end + 0.5  # Just to the right of the square
ax.annotate('', xy=(x_line, y2), xytext=(x_line, y1),
            arrowprops=dict(arrowstyle='<->', color='blue'))
ax.text(x_line + 0.8, (y1 + y2) / 2, r'$y_2 - y_1 = 5\Delta y$', rotation=90,
        ha='right', va='center', fontsize=10, color='black')

# Add temperature text outside the blue border
temperature_text = r'$T_a = 298 \, K$'

# Top side (above the blue border)
ax.text(cols / 2, rows + 0.5, temperature_text, ha='center', va='bottom', fontsize=12, color='black')

# Bottom side (below the blue border)
ax.text(cols / 2, -0.5, temperature_text, ha='center', va='top', fontsize=12, color='black')

# Left side (to the left of the blue border)
ax.text(-0.5, rows / 2, temperature_text, ha='right', va='center', fontsize=12, color='black', rotation=90)

# Right side (to the right of the blue border)
ax.text(cols + 0.5, rows / 2, temperature_text, ha='left', va='center', fontsize=12, color='black', rotation=90)

# Add a leader line pointing to the center of the orange square with text "Localized Heat Source"
center_x = (highlight_start + highlight_end) / 2
center_y = (highlight_start + highlight_end) / 2

ax.annotate('Localized Heat Source', xy=(center_x + 0.5, center_y + 0.5), 
            xytext=(center_x - 4.5, center_y + 1.5),  # Position of the text to the left of the square
            arrowprops=dict(arrowstyle="->", color='black'),
            fontsize=10, color='black')

# Fix axes
ax.set_xlim(0, cols)
ax.set_ylim(0, rows)
ax.set_aspect('equal')
ax.axis('off')

# Save the figure as a PNG image
plt.tight_layout()
plt.savefig('figure_1.png', dpi=150)

# Display the plot (optional)
plt.show()

