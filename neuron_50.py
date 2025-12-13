from neuron import n
from neuron.units import ms, mV, µm
import matplotlib.pyplot as plt
import random
n.load_file("stdrun.hoc")
class HHNeuron: 
    def __init__(self, gid):
        self.gid = gid 
        
        # 1. Create a section (soma)
        self.soma = n.Section(name=f'soma_{gid}')
        self.soma.L = 20 
        self.soma.diam = 20 
        self.soma.Ra = 100 
        self.soma.cm = 1 
        
        # 2. Insert Hodgkin-Huxley biophysics ('hh')
        self.soma.insert("hh")
        # Set parameters for the hh mechanism
        for seg in self.soma:
            seg.hh.gnabar = 0.12 
            seg.hh.gkbar = 0.036 
            seg.hh.gl = 0.0003
            seg.hh.el = -54.3  # mV

        # 3. Setup recording vectors
        self.v_vec = n.Vector()
        self.t_vec = n.Vector()
        self.v_vec.record(self.soma(0.5)._ref_v)
        self.t_vec.record(n._ref_t)

    def __del__(self):
        del self.soma

numNeurons = 50;
neurons = []
for i in range(numNeurons):
    neuron = HHNeuron(gid=i)
    neurons.append(neuron)
    stim = n.IClamp(neuron.soma(0.5))
    stim.delay = random.uniform(100, 200) * ms
    stim.dur = 500 * ms
    stim.amp = random.uniform(0.05, 0.15)  #
    neuron.stim = stim
    neuron.stim_delay = stim.delay
    neuron.stim_amp = stim.amp

syn_list = []
for pre in neurons:
    for post in neurons:
        if pre.gid != post.gid:
            if random.random() < 0.1:  # 10% connection probability
                syn = n.ExpSyn(post.soma(0.5))
                syn.tau = 2 * ms
                nc = n.NetCon(pre.soma(0.5)._ref_v, syn, sec=pre.soma)
                nc.weight[0] = 0.05
                nc.delay = 5 * ms
                syn_list.append((pre.gid, post.gid, nc))

spike_times = n.Vector()
spike_gids = n.Vector()
for neuron in neurons:
    spike_detector = n.NetCon(neuron.soma(0.5)._ref_v, None, sec=neuron.soma)
    spike_detector.threshold = 0  # mV
    spike_detector.record(spike_times, spike_gids, neuron.gid)
    neuron.spike_detector = spike_detector

T_STOP = 500  # ms
n.finitialize(-65) # Initialize V_m to -65 mV
n.t = 0
n.dt = 0.025 # Smaller time step for accuracy

print(f"Running simulation for {numNeurons} neurons...")
n.continuerun(T_STOP)
print("Simulation complete.")


# --- Plotting Section ---
# Convert NEURON vectors to Python lists for plotting
t_list = list(spike_times)
gid_list = list(spike_gids)

plt.style.use('seaborn-v0_8-whitegrid')
fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
fig.suptitle(f'Network Simulation ({numNeurons} HH Neurons)', fontsize=16)

# --- Subplot 1: 6. Plot Membrane Potential of a Sample Neuron ---
target_neuron = neurons[0] 
axes[0].plot(target_neuron.t_vec, target_neuron.v_vec, color='b')
axes[0].set_ylabel("Membrane Potential (mV)")
axes[0].set_title(f"Membrane Potential of Neuron {target_neuron.gid}")
axes[0].axhline(y=-20, color='r', linestyle='--', alpha=0.7, label='Spike Threshold')
axes[0].legend()
axes[0].set_xlim(0, T_STOP)


# --- Subplot 2: 7. Plot Spike Timings (Raster Plot) ---
axes[1].scatter(t_list, gid_list, s=10, c='k', marker='|') 
axes[1].set_xlabel("Time (ms)")
axes[1].set_ylabel("Neuron ID (GID)")
axes[1].set_title("Network Raster Plot (Spike Timings)")
axes[1].set_yticks(range(0, numNeurons, 5)) # Show every 5th GID
axes[1].set_ylim(-1, numNeurons)

plt.tight_layout(rect=[0, 0, 1, 0.96]) 

plt.savefig("network_simulation_results.png")
plt.show()

print("\nResults saved to 'network_simulation_results.png'")
print("Total spikes detected across all neurons: {len(t_list)}")