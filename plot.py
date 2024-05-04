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

data_e = data_e[2*fs*fs:-2*fs*fs]

print(f"Max distortion: {data_e.max()}")
x = np.linspace(-1, 1, 2*fs)

data_e_reshape = data_e.reshape((int((N-2*fs)/2), 2*fs))
plt.plot(x, data_e_reshape.T, linewidth=.1, color='blue')
plt.savefig("ewe.png", dpi = 300)
