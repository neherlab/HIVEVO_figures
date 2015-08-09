import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import seaborn as sns
from hivevo.patients import Patient
from hivevo.samples import all_fragments
from hivevo.sequence import alpha
plt.ion()
sns.set_style('darkgrid')

p = Patient.load('p3')
aft = p.get_allele_frequency_trajectories('RT1')
div = (aft*(1.0-aft)).sum(axis=1)
var_pos = div.max(axis=0)>0.1

plt.figure()
for pos in np.where(var_pos)[0]:
    for ni in range(5):
        traj = aft[:,ni,pos]
        if traj.max()>0.2 and traj[0]<0.5 : #and traj[-1]<0.2:
            plt.plot(p.ysi[~traj.mask], traj[~traj.mask], c = cm.jet(1.0*pos/aft.shape[-1]), label = str(pos+1)+alpha[ni], lw=2)
            print pos, alpha[ni], np.round(traj,2)

plt.ylabel('SNP frequency')
plt.xlabel('ETI [years]')
plt.legend(loc=2)
plt.savefig('mutatons_p3_RT1.pdf')


