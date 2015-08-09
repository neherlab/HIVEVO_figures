import os
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import seaborn as sns
from hivevo.patients import Patient
from hivevo.samples import all_fragments
from hivevo.sequence import alpha
from filenames import get_figure_folder
from util import store_data, load_data, fig_width, fig_fontsize, patients, patient_colors, HIVEVO_colormap
plt.ion()
sns.set_style('darkgrid')

username = os.path.split(os.getenv('HOME'))[-1]
foldername = get_figure_folder(username, 'first')


cmap = HIVEVO_colormap()
p = Patient.load('p1')
fig, axs = plt.subplots(2,3, sharey=True, sharex=True)
traj = []
ti = 3
tj = ti+1
for fi, frag in enumerate(all_fragments):
    ax = axs[fi//3,fi%3]
    aft = p.get_allele_frequency_trajectories(frag)
    pos = np.linspace(0,1,aft.shape[-1])
    for pi in xrange(aft.shape[-1]):
        for ni in xrange(5):
            if aft[0,ni,pi]<0.5 and (aft[ti,ni,pi]>0.2 and aft[tj,ni,pi]>0.2):
                traj.append([frag, pi, aft[:,ni,pi]])
    try:
        for ni in xrange(5):
            ind = (aft[ti,ni,:]*(1-aft[ti,ni,:])>0.01)|(aft[tj,ni,:]*(1-aft[tj,ni,:])>0.01)
            ax.scatter(aft[ti,ni,ind], aft[tj,ni,ind], c = [cmap(x) for x in pos[ind]])
    except:
        print 'fragment didnt work'

for ax in axs[:,0]:
    ax.set_ylabel('frequency at '+str(int(p.dsi[tj]))+' days')
    ax.locator_params(nbins=5)
    ax.tick_params(axis='both', labelsize = fig_fontsize)
for ax in axs[-1,:]:
    ax.set_xlabel('frequency at '+str(int(p.dsi[ti]))+' days')
    ax.locator_params(nbins=5)
    ax.tick_params(axis='both', labelsize = fig_fontsize)
for fi, frag in enumerate(all_fragments):
    ax = axs[fi//3][fi%3]
    ax.text(0.8,0.02, frag, fontsize=1.5*fig_fontsize, transform=ax.transAxes)
ax.set_ylim([-0.02,1.02])
ax.set_xlim([-0.02,1.02])
plt.tight_layout()

for fmt in ['.pdf', '.png', '.svg']:
    plt.savefig(foldername+'nearby_freq'+fmt)


plt.figure()
dt_dt = []
for t in traj:
    if t[0]=='F5':
        continue
    plt.plot(p.ysi, t[-1], c=cm.jet(min(1,max(0,(t[-1][tj]-t[-1][ti])/1+0.5))))
    dt_dt.append([t[-1][tj]-t[-1][ti], t[-1][tj+1]-t[-1][ti-1]])

dt_dt=np.array(dt_dt)
plt.figure()
plt.scatter(dt_dt[:,0], dt_dt[:,1])
