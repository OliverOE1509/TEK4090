import matplotlib.pyplot as plt
import numpy as np

# Generate data for the main Nyquist path (positive frequencies)
omega = np.logspace(-2, 2, 1000)
omega = np.sort(np.concatenate((omega, [0.99, 1.0, 1.01, np.sqrt(2)-0.001, np.sqrt(2)+0.001])))
s = 1j * omega
P_jw = (s**2 + 1) / ((s**2 + 2) * (s - 2))

# Create the plot
fig, ax = plt.subplots(figsize=(8, 8))

# Plot the Nyquist plot for positive omega
ax.plot(np.real(P_jw), np.imag(P_jw), 'b-', linewidth=2, label='Nyquist Plot ($\\omega > 0$)')
# Plot the Nyquist plot for negative omega (mirror image)
ax.plot(np.real(P_jw), -np.imag(P_jw), 'b--', linewidth=1, alpha=0.5, label='Nyquist Plot ($\\omega < 0$)')

# Draw the infinite arc manually (clockwise from +inf to -inf on the real axis)
theta = np.linspace(0, -np.pi, 100)  # Clockwise angle change
arc_radius = 10  # Large radius to represent infinity
arc_x = arc_radius * np.cos(theta) + (-0.5)  # Center the arc near the critical point
arc_y = arc_radius * np.sin(theta)
ax.plot(arc_x, arc_y, 'r--', linewidth=2, label='Infinite Arc (Clockwise)')

# Mark key points
ax.plot(-0.25, 0, 'ko', markersize=5, label='Start: $P(j0) = -0.25$')
ax.plot(0, 0, 'go', markersize=5, label='Zero: $P(j1) = 0$')
ax.plot(-1, 0, 'rs', markersize=8, label='Critical Point: $-1 + 0j$')

# Draw asymptotes: as w->inf, approach origin from +90 direction
ax.arrow(0, 0.5, 0, -0.4, head_width=0.05, head_length=0.1, fc='g', ec='g', label='Asymptote to Origin')
# As w->inf from negative side, approach from -90 direction
ax.arrow(0, -0.5, 0, 0.4, head_width=0.05, head_length=0.1, fc='g', ec='g')

# Draw the real axis asymptote for the infinite arc
ax.axhline(0, color='k', linestyle=':', alpha=0.3)
ax.axvline(0, color='k', linestyle=':', alpha=0.3)

# Add labels and title
ax.set_title('Nyquist Plot of $P(s) = \\frac{s^2 + 1}{(s^2 + 2)(s - 2)}$\nWith Asymptotes and Infinite Arc')
ax.set_xlabel('Real Axis')
ax.set_ylabel('Imaginary Axis')
ax.legend(loc='upper right')
ax.grid(True, linestyle='--', alpha=0.7)
ax.axis('equal')
ax.set_xlim([-3, 2])
ax.set_ylim([-3, 3])
plt.show()