import matplotlib.pyplot as plt
import numpy as np

data_t = np.genfromtxt("results/rc_filter_t.txt", delimiter='\n')
data_f = np.genfromtxt("results/rc_filter_f.txt", delimiter='\n')

plt.plot(data_t)
plt.show()

plt.plot(data_f)
plt.show()

# EYE DIAGRAM PLOT 
data_e = np.genfromtxt("results/rc_filter_eye_diagram.txt", delimiter='\n')

fs = 10
N = 100000

data_e = data_e[200:]

x = np.linspace(-0.99, 0.99, fs*2)
data_e_reshape = data_e.reshape((int(N/2), fs*2))

plt.plot(x, data_e_reshape.T, linewidth=.1, color='blue')
plt.show()

