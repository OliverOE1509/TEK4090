import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(-5, 5, 100)

def plot_convex(a, b, f):
    if b <= a:
        raise ValueError("b must be greater than a")
    
    y = f(x)
    c = np.random.choice(np.linspace(a, b, 10))
    
    k = (f(b) - f(a)) / (b - a)  
    k_1 = (f(c) - f(a)) / (c - a)  
    k_2 = (f(b) - f(c)) / (b - c)  
    
    secant = f(a) + k * (x - a)  
    secant_left = f(a) + k_1 * (x - a)  
    secant_right = f(c) + k_2 * (x - c)  
    
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, label=r"$f(x)$", color="blue", linewidth=2)
    plt.plot(x, secant, "--", color="orange", alpha=0.8)
    plt.plot(x, secant_left, "--", color="green", alpha=0.8)
    plt.plot(x, secant_right, "--", color="purple", alpha=0.8)
    
    plt.scatter([a, b, c], [f(a), f(b), f(c)], color="red", s=50)
    plt.annotate(f"a=({a:.1f}, {f(a):.1f})", (a, f(a)))
    plt.annotate(f"b=({b:.1f}, {f(b):.1f})", (b, f(b)))
    plt.annotate(f"c=({c:.1f}, {f(c):.1f})", (c, f(c)))
    
    x_fill = np.linspace(a, b, 50)
    y_fill = f(x_fill)
    secant_fill = f(a) + k * (x_fill - a)
    plt.fill_between(x_fill, y_fill, secant_fill, alpha=0.2, color='red', 
                     label='Area showing convexity')
    
    plt.title("Illustration of Convex Function Properties")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.axis('equal')
    plt.show()

plot_convex(-2, 3, lambda x: x**2)

