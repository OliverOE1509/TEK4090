import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

def manual_d_contour(poles_on_imaginary=None, R=1e6, epsilon=1e-6, num_points=1000):
	"""
	Manual D-contour parametrization with explicit control

	Parameters:
	- poles_on_imaginary: list of poles on imaginary axis
	- R: large radius
	- epsilon: small radius around poles
	- num_points: points per segment
	"""
	s_contour = []
	print(poles_on_imaginary)
	# Segment 1: Large semicircle (Re(s) > 0)
	theta = np.linspace(-np.pi/2, np.pi/2, num_points)
	s_contour.extend(R * np.exp(1j * theta))

	# Segment 2: Positive imaginary axis (from +jR to +jε)
	omega_pos = np.linspace(R, epsilon, num_points)
	s_contour.extend(1j * omega_pos)

	# Handle poles on imaginary axis with small semicircles
	if poles_on_imaginary:
		for pole in poles_on_imaginary:
			if pole.imag > 0:
				# Small semicircle around pole
				theta_small = np.linspace(-np.pi/2, np.pi/2, num_points//4)
				s_small = pole + epsilon * np.exp(1j * theta_small)
				s_contour.extend(s_small)
	    
	# Segment 3: Negative imaginary axis (from -jε to -jR)
	omega_neg = np.linspace(-epsilon, -R, num_points)
	s_contour.extend(1j * omega_neg)

	return np.array(s_contour)

def evaluate_transfer_function(system, s_points):
    """Evaluate transfer function at complex points"""
    if hasattr(system, '__call__'):
        return system(s_points)
    else:
        # For scipy systems
        w, H = signal.freqresp(system, w=s_points.imag)
        return H

# Complete example
def complete_nyquist_example():
    # System with poles at s = 0 and s = ±j
    num = [1]
    den = [1, 0, 1]  # s³ + s² + s
    
    system = signal.TransferFunction(num, den)
    poles = signal.tf2zpk(num, den)[1]
    
    # Get imaginary axis poles
    imag_poles = [p for p in poles if abs(p.real) < 1e-10]
    
    # Create D-contour
    s_contour = manual_d_contour(poles_on_imaginary=imag_poles)
    
    # Evaluate system
    H = evaluate_transfer_function(system, s_contour)
    
    # Plot
    plt.figure(figsize=(12, 10))
    
    # Nyquist plot
    plt.subplot(2, 2, 1)
    plt.plot(H.real, H.imag, 'b-', linewidth=2)
    plt.axhline(0, color='black', linestyle='--', alpha=0.5)
    plt.axvline(0, color='black', linestyle='--', alpha=0.5)
    plt.xlabel('Real')
    plt.ylabel('Imaginary')
    plt.title('Nyquist Plot')
    plt.grid(True, alpha=0.3)
    plt.axis('equal')
    
    # D-contour
    plt.subplot(2, 2, 2)
    plt.plot(s_contour.real, s_contour.imag, 'r-', linewidth=2)
    plt.axhline(0, color='black', linestyle='--', alpha=0.5)
    plt.axvline(0, color='black', linestyle='--', alpha=0.5)
    plt.xlabel('Real(s)')
    plt.ylabel('Imag(s)')
    plt.title('D-contour in s-plane')
    plt.grid(True, alpha=0.3)
    plt.axis('equal')
    
    # Frequency response magnitude
    plt.subplot(2, 2, 3)
    plt.plot(s_contour.imag, np.abs(H), 'g-', linewidth=2)
    plt.xlabel('Frequency (rad/s)')
    plt.ylabel('|H(jω)|')
    plt.title('Magnitude Response')
    plt.grid(True, alpha=0.3)
    plt.xscale('log')
    plt.yscale('log')
    
    # Phase response
    plt.subplot(2, 2, 4)
    plt.plot(s_contour.imag, np.angle(H), 'm-', linewidth=2)
    plt.xlabel('Frequency (rad/s)')
    plt.ylabel('Phase (rad)')
    plt.title('Phase Response')
    plt.grid(True, alpha=0.3)
    plt.xscale('log')
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    complete_nyquist_example()