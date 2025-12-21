#Trying to replicate the network simulation as shown in the Wang and Buzsaki paper
from neuron import n
from neuron import h, load_mechanisms
from neuron.units import ms, mV, µm
import matplotlib.pyplot as plt
import random
import datetime
n.load_file("stdrun.hoc")
# load_mechanisms(".")
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
        self.soma.insert("kdr")
        self.soma.insert("naf")
        # Set parameters for the hh mechanism
        for seg in self.soma:
            # seg.hh.gnabar = 0.12 
            # seg.hh.gkbar = 0.036 
            seg.naf.g = 0.035  # mS/cm^2
            seg.kdr.g = 0.009  # mS/cm^2
            seg.hh.gl = 0.0001
            seg.hh.el = -65  # mV

        # 3. Generate Synapse
        self.syn = n.ExpSyn(self.soma(0.5))
        self.syn.tau = 10 * ms
        self.syn.e = -75  # mV

        # 4. Setup recording vectors
        self.v_vec = n.Vector()
        self.t_vec = n.Vector()
        self.v_vec.record(self.soma(0.5)._ref_v)
        self.t_vec.record(n._ref_t)
        self.syn_i = n.Vector()
        self.syn_i.record(self.syn._ref_i)
        

        

    def __del__(self):
        del self.soma

random.seed(67890)
numNeurons = 100;
neurons = []
for i in range(numNeurons):
    neuron = HHNeuron(gid=i)
    neurons.append(neuron)
    stim = n.IClamp(neuron.soma(0.5))
    stim.delay = 50.0 * ms
    stim.dur = 500 * ms
    stim.amp = random.uniform(0.05, 0.15)  #
    neuron.stim = stim
    neuron.stim_delay = stim.delay
    neuron.stim_amp = stim.amp

probability = 0.1


syn_list = []
for pre in neurons:
    for post in neurons:
        if pre.gid != post.gid:
            if random.random() < probability: 
                # print("Connection from", pre.gid, "to", post.gid)
                nc = n.NetCon(pre.soma(0.5)._ref_v, post.syn, sec=pre.soma)
                nc.weight[0] = 0.00001
                nc.delay = 1 * ms

                
                syn_list.append((pre.gid, post.gid, nc))


spike_times = n.Vector()
spike_gids = n.Vector()
for neuron in neurons:
    spike_detector = n.NetCon(neuron.soma(0.5)._ref_v, None, sec=neuron.soma)
    spike_detector.threshold = 0   # mV
    spike_detector.record(spike_times, spike_gids, neuron.gid)
    neuron.spike_detector = spike_detector

T_STOP = 600  # ms
n.finitialize(-52) # Initialize V_m to -65 mV
n.t = 0
n.dt = 0.025 # Smaller time step for accuracy
t = n.Vector().record(n._ref_t)
print(f"Running simulation for {numNeurons} neurons...")
n.continuerun(T_STOP)
print("Simulation complete.")
# print(syn_record_list[0])

# --- Plotting Section ---
# Convert NEURON vectors to Python lists for plotting
t_list = list(spike_times)
gid_list = list(spike_gids)

plt.style.use('seaborn-v0_8-whitegrid')
fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
fig.suptitle(f'Network Simulation ({numNeurons} HH Neurons)', fontsize=16)

# --- Subplot 1: 6. Plot Membrane Potential of a Sample Neuron ---
target_neuron = neurons[0]
second_neuron = neurons[1] 
axes[0].plot(target_neuron.t_vec, target_neuron.v_vec, color='b')
axes[0].plot(second_neuron.t_vec, second_neuron.v_vec, color='g')
axes[0].set_ylabel("Membrane Potential (mV)")
axes[0].set_title(f"Membrane Potential of Neuron {target_neuron.gid}")
axes[0].axhline(spike_detector.threshold, color='r', linestyle='--', alpha=0.7, label='Spike Threshold')
axes[0].legend()
axes[0].set_xlim(0, T_STOP)


# --- Subplot 2: 7. Plot Spike Timings (Raster Plot) ---
axes[1].scatter(t_list, gid_list, s=10, c='k', marker='|') 
axes[1].set_xlabel("Time (ms)")
axes[1].set_ylabel("Neuron ID (GID)")
axes[1].set_title("Network Raster Plot (Spike Timings)")
axes[1].set_yticks(range(0, numNeurons, 5)) # Show every 5th GID
axes[1].set_ylim(-1, numNeurons)

axes[2].plot(target_neuron.t_vec, target_neuron.syn_i, color='m')

timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")
filename = f"figures/network_simulation_results_{probability*100}%_{timestamp}.png"

plt.tight_layout(rect=[0, 0, 1, 0.96]) 

plt.savefig(filename)
plt.show()

print("\nResults saved to 'network_simulation_results.png'")
print(f"Total spikes detected across all neurons: {len(t_list)}")