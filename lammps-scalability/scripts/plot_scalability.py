import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
import numpy as np

n_bosons = []
run_times = []
errors = [] 

for fname in snakemake.input:
    # assume that the data is #bosons time(sec)
    data = np.loadtxt(fname, skiprows=1, delimiter=',')
    n_bosons.append(data[0,0])
    run_times.append(np.mean(data[:,1]))
    errors.append(np.std(data[:,1]))

n = np.array(n_bosons)
plt.bar(n, run_times, yerr=errors)
#plt.loglog(n, np.array(run_times), '.')

#plt.gca().loglog(n, (np.max(run_times)/1024**2)*n**2, 'r-')
plt.show()
plt.savefig(snakemake.output[0])
