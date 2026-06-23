import numpy as np
import matplotlib.pyplot as plt

# load data
ini = np.loadtxt('freefall_initial.dat')
fin = np.loadtxt('freefall_final.dat')

# convert radius to AU
r = ini[:, 0] / 1.495978707e13

plt.figure(figsize=(10,4))

# --- Density ---
plt.subplot(1,2,1)
plt.plot(r, ini[:,1], label='initial')
plt.plot(r, fin[:,1], '--', label='final')
plt.xlabel('r (AU)')
plt.ylabel('rho (g/cm^3)')
plt.title('Density')
plt.legend()

# --- Velocity ---
plt.subplot(1,2,2)
plt.plot(r, fin[:,2]/1e5)
plt.axhline(0)
plt.xlabel('r (AU)')
plt.ylabel('v (km/s)')
plt.title('Velocity (final)')

plt.tight_layout()
plt.show()