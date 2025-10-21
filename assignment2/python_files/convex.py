import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(-5, 5, 100)

def plot_convex(a, b, f):
    if b<a:
        raise ValueError("b must be greater than a")
    
    y = f(x)
    c = np.random.choice(np.linspace(a, b, 10))
    k_1 = (f(c) - f(a)) / (c-a)
    k_2 = (f(b) - f(c)) / (b-c)
    k = (f(b) - f(a)) / (b-a)
    
    secant = f(a) + k * (x - a)
    secant2 = f(a) - k_1 * (x - a)
    secant3 = f(c) + k_2 * (x - c)

    # Plot
    plt.figure(figsize=(6,4))
    plt.plot(x, y, label=r"$f(x)$", color="blue")
    plt.plot(x, secant, "--", label="Secant line, K", color="orange")
    plt.plot(x, secant2, "--", label="Secant line 2, K_1", color="green")
    plt.plot(x, secant3, "--", label="Secant line 3, K_2", color="purple")
    plt.scatter([a, b, c], [f(a), f(b), f(c)], color="red", zorder=5)
    plt.title("Illustration of Convexity")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    #plt.ylim(min(y)-1, max(y)+1)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()


plot_convex(-2, 3, lambda x: x**2)