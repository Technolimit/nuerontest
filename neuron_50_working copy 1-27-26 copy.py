#Mimic Fig 3c
import numpy as np
from neuron import n, h, load_mechanisms
from neuron.units import ms, mV, µm
import matplotlib.pyplot as plt
import pickle as pl
import random
import datetime
n.load_file("stdrun.hoc")
load_mechanisms(".")
class HHNeuron: 
    def __init__(self, gid):
        self.gid = gid 
        self.soma = n.Section(name=f'soma_{gid}')
        self.soma.L = 20 
        self.soma.diam = 20 
        self.soma.Ra = 100 
        self.soma.cm = 1

        self.soma.insert("hh")
        for seg in self.soma:
            seg.hh.gnabar = 0.12 
            seg.hh.gkbar = 0.036 
            seg.hh.gl = 0.0001
            seg.hh.el = -65

        # Inhibitory Synapse
        self.syn = n.FDSExp2Syn(self.soma(0.5))
        self.syn.tau1 = 0.5   
        self.syn.tau2 = 2.0   
        self.syn.e = -75

        self.syn.f = 1
        self.syn.d2 = 1
        self.syn.d1 = 0.4     # Fast depression
        self.syn.tau_d1 = 300

        # Excitatory Synapse
        self.syn_exc = n.Exp2Syn(self.soma(0.5))
        self.syn_exc.e = 0        
        self.syn_exc.tau1 = 0.2
        self.syn_exc.tau2 = 2.0

        
        self.v_vec = n.Vector().record(self.soma(0.5)._ref_v)
        

random.seed(67890)
numNeurons = 100;
neurons = []

def run_experiment(exc_w, inh_w, inh_tau, prob, t_stop=600):
    # 1. Reset NEURON state
    h('forall delete_section()')
    
    numNeurons = 100
    neurons = [HHNeuron(i) for i in range(numNeurons)]
    
    # Store references to prevent Python from deleting them prematurely
    all_stims = []
    all_ncs = []

    # Setup Excitation
    for neuron in neurons:
        ns = n.NetStim()
        ns.interval, ns.number, ns.noise, ns.start = 100, 1e9, 1, 0
        nc = n.NetCon(ns, neuron.syn_exc)
        nc.weight[0] = exc_w
        all_stims.append(ns)
        all_ncs.append(nc)

    # Setup Inhibitory Connections
    for pre in neurons:
        for post in neurons:
            if pre.gid != post.gid and random.random() < prob:
                nc = n.NetCon(pre.soma(0.5)._ref_v, post.syn, sec=pre.soma)
                nc.weight[0] = inh_w
                nc.delay = 1 * ms
                all_ncs.append(nc)

    # Run Simulation
    n.finitialize(-65)
    n.dt = 0.025
    n.continuerun(t_stop)

    v_all = np.array([np.array(neuron.v_vec) for neuron in neurons])

    v_avg = np.mean(v_all, axis=0)
    v_normalized = v_avg - np.mean(v_avg)
    
    N = len(v_normalized) 
    yf = np.fft.fft(v_normalized)
    xf = np.fft.fftfreq(N, d=n.dt/1000)
    magnitude = np.abs(yf) * 2 / N
    
    mask = (xf >= 1) & (xf <= 400)
    if not any(mask):
        return 0.0
    
    mag_plot = magnitude[mask]
    xf_plot = xf[mask]
    fig, axes = plt.subplots(1, 1, figsize=(12, 10), sharex=False)
    fig.suptitle(f'Mean Membrane Potential Spikes with {exc_w} excitation and {inh_w} inhibition', fontsize=16)
    axes.plot(xf_plot, mag_plot, color='m')
    axes.set_xlim(0, 400)
    axes.set_ylabel("Amplitude")
    axes.set_xlabel("Frequency (Hz)")
    name = f"figures/FFT_exc_{exc_w}_inh_{inh_w}_tau_{inh_tau}_prob_{prob}.png"
    plt.savefig(name)
    plt.close(fig)
    return xf_plot[np.argmax(mag_plot)]


inh_weights = [inh/1000 for inh in range(0, 55, 5)]
exc_values = [exc/1000 for exc in range(0, 55, 5)]
print(inh_weights)
print(exc_values)
results = np.zeros((len(inh_weights), len(exc_values)))

for i, iw in enumerate(inh_weights):
    print(f"Running Inhibitory Weight: {iw}")
    for j, exc in enumerate(exc_values):
        max_freq = run_experiment(exc_w=exc, inh_w=iw, inh_tau=1.0, prob=1.0)
        results[j, i] = max_freq

# X, Y = np.meshgrid(exc_values, inh_weights) 

# fig = plt.figure(figsize=(12, 8))
# ax = fig.add_subplot(111, projection='3d')


# surf = ax.plot_surface(X, Y, results, cmap='viridis', edgecolor='none', alpha=0.9)


# ax.set_xlabel('Inhibitory Weight')
# ax.set_ylabel('Excitatory Weight')
# ax.set_zlabel('Max Frequency (Hz)')
# ax.set_title('Inhibitory vs Excitatory Weights: Effect on Max Frequency')

# fig.colorbar(surf, shrink=0.5, aspect=10, label='Frequency (Hz)')
# ax.view_init(elev=30, azim=225) 
# timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
# title = f"figures/simulationResults_inh_vs_exc_weights_{timestamp}.png"
# pl.dump(fig, open(f"figures/interact_inh_vs_exc_weights_{timestamp}.pkl", "wb"))
# # plt.show()


fig, ax = plt.subplots(figsize=(10, 8))
im = ax.imshow(results, aspect='auto', origin='lower', cmap='plasma')
plt.colorbar(im, label='Max FFT Frequency (Hz)')

for i in range(len(inh_weights)):
    for j in range(len(exc_values)):

        text = ax.text(j, i, f'{results[i, j]:.1f}',
                       ha="center", va="center", color="w", 
                       fontsize=9, fontweight='bold')


ax.set_xticks(range(len(exc_values)))
ax.set_xticklabels([f"{x:.3f}" for x in exc_values])

ax.set_yticks(range(len(inh_weights)))
ax.set_yticklabels([f"{x:.3f}" for x in inh_weights])

ax.set_xlabel('Excitatory Weight') 
ax.set_ylabel('Inhibitory Weight') 
ax.set_title('Max Frequency per E/I Parameter Pair')

plt.tight_layout()
plt.show()

