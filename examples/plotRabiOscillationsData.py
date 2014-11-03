import numpy as np
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import matplotlib.pyplot as plt

data = np.loadtxt('rabioscillations.txt')
plt.plot(data[:, 0], data[:, 1], '-o')
plt.plot(data[:, 0], data[:, 2], '-o')
plt.xlabel(r'$t/\gamma^{-1}$')
plt.ylabel(r'$\rho_ee(t), |\rho_{ge}(t)|$')
plt.gcf().set_size_inches(4, 3)
plt.gcf().subplots_adjust(bottom = 0.15, left = 0.15, top = 0.97, right = 0.97)
plt.savefig('rabiOscillations.pdf')

