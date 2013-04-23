import numpy as np
import pylab as plt
from math import sqrt

golden_ratio = (1.0 + sqrt(5.0)) / 2.0

wl, tau, J, H, K, B = np.genfromtxt('moments.dat', unpack=True, skip_header=1)

figprops = dict(figsize=(11.0, 11.0 / golden_ratio), dpi=128)
adjustprops = dict(left=0.10, bottom=0.10, top=0.85, wspace=0.30, hspace=0.40)
legendprops = {'legend.fontsize': 10, 'legend.linewidth': 2}

fig = plt.figure(**figprops)
fig.subplots_adjust(**adjustprops)
plt.rcParams.update(**legendprops)

fig.suptitle(r'NLTE; $\epsilon= 10^{-4}$; $\lambda = 1000 \AA$')

J_plot = fig.add_subplot(321)
J_plot.plot(tau[:100], J[:100])
J_plot.set_xscale('log')
J_plot.set_xlabel(r'$\tau_\lambda$')
J_plot.set_ylabel(r'$J_\lambda$')
J_plot.tick_params(axis='both', which='major', labelsize=8)
J_plot.tick_params(axis='both', which='minor', labelsize=6)

H_plot = fig.add_subplot(322)
H_plot.plot(tau[:100], H[:100])
H_plot.set_xscale('log')
H_plot.set_xlabel(r'$\tau_\lambda$')
H_plot.set_ylabel(r'$H_\lambda$')
H_plot.tick_params(axis='both', which='major', labelsize=8)
H_plot.tick_params(axis='both', which='minor', labelsize=6)

K_plot = fig.add_subplot(323)
K_plot.plot(tau[:100], K[:100])
K_plot.set_xscale('log')
K_plot.set_xlabel(r'$\tau_\lambda$')
K_plot.set_ylabel(r'$K_\lambda$')
K_plot.tick_params(axis='both', which='major', labelsize=8)
K_plot.tick_params(axis='both', which='minor', labelsize=6)

B_plot = fig.add_subplot(324)
B_plot.plot(tau[:100], B[:100])
B_plot.set_xscale('log')
B_plot.set_xlabel(r'$\tau_\lambda$')
B_plot.set_ylabel(r'$B_\lambda$')
B_plot.tick_params(axis='both', which='major', labelsize=8)
B_plot.tick_params(axis='both', which='minor', labelsize=6)

KJ_plot = fig.add_subplot(325)
KJ_plot.plot(tau[:100], K[:100]/J[:100])
KJ_plot.set_xscale('log')
KJ_plot.set_xlabel(r'$\tau_\lambda$')
KJ_plot.set_ylabel(r'$K_\lambda / J_\lambda$')
KJ_plot.tick_params(axis='both', which='major', labelsize=8)
KJ_plot.tick_params(axis='both', which='minor', labelsize=6)

HJ_plot = fig.add_subplot(326)
HJ_plot.plot(tau[:100], H[:100]/J[:100])
HJ_plot.set_xscale('log')
HJ_plot.set_xlabel(r'$\tau_\lambda$')
HJ_plot.set_ylabel(r'$H_\lambda / J_\lambda$')
HJ_plot.tick_params(axis='both', which='major', labelsize=8)
HJ_plot.tick_params(axis='both', which='minor', labelsize=6)

plt.show()
