import numpy as np
import matplotlib.pyplot as plt

f = lambda x: 20 + 5 * np.sin(2 * np.pi * x / 86400) 

x = np.linspace(0, 86400, 10000)


plt.plot(x, f(x))
plt.grid(True)
plt.show()
