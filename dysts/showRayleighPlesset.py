import matplotlib.pyplot as plt
import numpy as np

from dysts.flows import RayleighPlesset
from dysts.base import make_trajectory_ensemble

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

plt.rc('text', usetex=True)

# the value dt is 100, so one cycle is determined in .json file.
model0 = RayleighPlesset()

n = int(500)
_method = "DOP853"
_resample = True
_return_times = True
_postprocess = False  # keep the raw phase, used as time in fig0

tpts, sol0 = model0.make_trajectory(n, method=_method, resample=_resample,
                                    return_times=_return_times,
                                    postprocess=_postprocess)

fig1 = plt.figure("Rayleigh-Plesset Phase")
ax1 = fig1.add_subplot(1, 1, 1)
ax1.grid()
ax1.plot(sol0[:, 0], sol0[:, 1], color='tab:blue', marker='.', markersize=10,
         linestyle='solid', linewidth=2)
ax1.set_xlabel(r'$R$ (m)')
ax1.set_ylabel(r'$\dot{R}$ (m/s)')

fig0 = plt.figure("Rayleigh-Plesset Oscillations 1")
ax0 = fig0.add_subplot(1, 1, 1)
ax0.grid()
ax0.plot(sol0[:, 2] / (2.0 * np.pi), sol0[:, 0], color='tab:green', marker='.',
         markersize=10, linestyle='solid', linewidth=2)
ax0.set_xlabel(r'$t$')
ax0.set_ylabel(r'$\dot{R}$ (m/s)')

fig2 = plt.figure("Rayleigh-Plesset Oscillations 2")
ax2 = fig2.add_subplot(1, 1, 1)
ax2.grid()
ax2.plot(tpts * model0.omega / (2.0 * np.pi), sol0[:, 0], color='tab:orange',
         marker='.', markersize=10, linestyle='solid', linewidth=2)
ax2.set_xlabel(r'$t$')
ax2.set_ylabel(r'$\dot{R}$ (m/s)')

fig3 = plt.figure("Rayleigh-Plesset Oscillations 3")
ax3 = fig3.add_subplot(1, 1, 1)
ax3.grid()
ax3.plot(tpts / model0.period, sol0[:, 0], color='tab:red', marker='.',
         markersize=10, linestyle='solid', linewidth=2)
ax3.set_xlabel(r'$t$')
ax3.set_ylabel(r'$\dot{R}$ (m/s)')

# View subset of attractors

_subset = ['RayleighPlesset', 'RikitakeDynamo', 'Rossler', 'Rucklidge']
all_trajectories = make_trajectory_ensemble(10000, subset=_subset,
                                            method="Radau", resample=True)

if (len(all_trajectories) > 10):
    num_cols = 10
else:
    num_cols = len(all_trajectories)
num_rows = int(np.ceil(len(all_trajectories) / num_cols))
fig = plt.figure(figsize=(num_cols*4, num_rows*4))

gs = plt.matplotlib.gridspec.GridSpec(num_rows, num_cols)
gs.update(wspace=0.0, hspace=0.05)

for i, attractor_name in enumerate(all_trajectories):
    sol = all_trajectories[attractor_name]
    if (attractor_name == 'RayleighPlesset'):
        plt.subplot(gs[i])
        plt.plot(sol[:, 0], sol[:, 1], 'tab:blue', linewidth=0.25)
        plt.title(attractor_name, y=-0.05)
        plt.gca().axis('off')
    else:
        plt.subplot(gs[i])
        plt.plot(sol[:, 0], sol[:, 1], 'k', linewidth=0.25)
        plt.title(attractor_name, y=-0.05)
        plt.gca().axis('off')

dpi = 150
pad = 0.0
pad_inches = 0.02
remove_border = False
if remove_border:
    plt.gca().set_axis_off()
    plt.subplots_adjust(top=1+pad, bottom=0+pad, right=1+pad, left=0+pad,
                        hspace=0, wspace=0)
    plt.margins(0+pad, 0+pad)
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())

plt.savefig('subset_attractors_RP.png', bbox_inches='tight',
            pad_inches=pad_inches, dpi=dpi)

plt.show()
