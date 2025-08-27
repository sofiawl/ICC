import numpy as np
import matplotlib.pyplot as plt

# Reading the data we need
data = []
with open("polinomios.dat", 'r') as infile:
    all_lines = infile.readlines()

for i in range(0, len(all_lines), 4):
    chunk = all_lines[i:i+4]
    if len(chunk) >= 3:
        line2 = chunk[1].strip()
        line3 = chunk[2].strip()
        data.append((line2, line3))

# Plotting
for polynomial in data:
    coefficients = [float(x) for x in polynomial[0].split()]

    # Getting the intervals
    a = float(polynomial[1].split()[0])
    b = float(polynomial[1].split()[1])

    # Calculating roots
    all_roots = np.roots(coefficients)
    real_roots = all_roots[np.isreal(all_roots)].real

    # The view that will be shown
    x = np.linspace(a - (b - a) * 0.1, b + (b - a) * 0.1, 400)
    y = np.polyval(coefficients, x)

    plt.figure(figsize=(8, 6))

    # The actual interval we got from the file
    plt.axvspan(a, b, color="#CBCBCB", label=f'Interval [{a}, {b}]')

    # Plotting the roots
    for root in real_roots:
        # We plot just the root that our program will find
        if a <= root <= b:
            plt.plot(root, 0, 'ro', markersize=8, label=f'Real Root')

            plt.annotate(
                f'Root ≈ {root:.2f}',
                xy=(root, 0),
                xytext=(0, 30),
                textcoords='offset points',
                ha='center',
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3")
            )

    # Label info for our plot
    plt.plot(x, y)
    plt.xlabel("x")
    plt.ylabel("y")
    title_str = 'f(x) = ' + ' + '.join([f'{c:.2g}x^{len(coefficients)-1-i}'
                                        for i, c in enumerate(coefficients)])
    plt.title(title_str, pad=20)
    plt.grid(True)
    plt.legend()

    plt.show()