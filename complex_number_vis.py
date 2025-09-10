import numpy as np
import matplotlib.pyplot as plt



z1 = complex(2, 3)
z2 = complex(2, -3)


# Plotting a single complex number
plt.scatter(z1.real, z1.imag, color='blue', label='z1')
plt.scatter(z2.real, z2.imag, color='red', label='z2')

# Plotting as a vector from the origin
plt.arrow(0, 0, z1.real, z1.imag, head_width=0.2, head_length=0.3, fc='red', ec='red', label='z2')
plt.arrow(0, 0, z2.real, z2.imag, head_width=0.2, head_length=0.3, fc='red', ec='red', label='z2')

# Adding labels and grid
plt.xlabel('Real Part')
plt.ylabel('Imaginary Part')
plt.axhline(0, color='gray', linewidth=0.5) # x-axis
plt.axvline(0, color='gray', linewidth=0.5) # y-axis
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.title('Argand Diagram of Complex Numbers')
plt.axis('equal') # Ensures equal scaling for x and y axes
plt.show()