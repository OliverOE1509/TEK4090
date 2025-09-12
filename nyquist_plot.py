import numpy as np
import matplotlib.pyplot as plt

# Create the Nyquist plot manually
omega = np.linspace(0.1, 100, 1000)  # Positive frequencies
L_jw = 1 / (1j * omega)  # For K=1

plt.figure(figsize=(8, 8))
plt.plot(np.real(L_jw), np.imag(L_jw), 'b-', label='ω > 0')
plt.plot(np.real(L_jw), -np.imag(L_jw), 'b--', label='ω < 0')  # Mirror image

# Mark key points
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.plot(-1, 0, 'ro', markersize=8, label='Critical Point (-1, 0)')

plt.title('Nyquist Plot of L(s) = 1/s')
plt.xlabel('Real Axis')
plt.ylabel('Imaginary Axis')
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.axis([-1.5, 0.5, -1.5, 1.5])
plt.show()