#Mimic Fig 3c
from neuron import n
import numpy as np
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
        # self.soma.insert("kdr")
        # self.soma.insert("naf")
        # Set parameters for the hh mechanism
        for seg in self.soma:
            seg.hh.gnabar = 0.12 
            seg.hh.gkbar = 0.036 
            # seg.naf.g = 0.035  # mS/cm^2
            # seg.kdr.g = 0.009  # mS/cm^2
            seg.hh.gl = 0.0001
            seg.hh.el = -65  # mV

        # 3. Generate Synapse
        self.syn = n.ExpSyn(self.soma(0.5))
        self.syn.tau = 1 * ms
        self.syn.e = -75  # mV
        #5 new synapse for excitation
        self.syn_exc = n.ExpSyn(self.soma(0.5))
        self.syn_exc.tau = 2 * ms
        self.syn_exc.e = 0
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
    # stim.delay = 50.0 * ms
    stim.dur = 1e9 * ms
    # stim.amp = 10*400*n.PI/1e5
    stim.amp = 0
    neuron.stim = stim
    neuron.stim_delay = stim.delay
    neuron.stim_amp = stim.amp

probability = 1

all_stims = []
for neuron in neurons:
    ns = n.NetStim()
    ns.interval = 100
    ns.number = 1e9     
    ns.noise = 1     
    ns.start = 0

    nc = n.NetCon(ns, neuron.syn_exc) 
    nc.weight[0] = 0.002
    all_stims.append((ns, nc))

syn_list = []
for pre in neurons:
    for post in neurons:
        if pre.gid != post.gid:
            if random.random() < probability: 
                # print("Connection from", pre.gid, "to", post.gid)
                nc = n.NetCon(pre.soma(0.5)._ref_v, post.syn, sec=pre.soma)
                nc.weight[0] = 0.005 #Inhibitory weight
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
# for neuron in neurons:
    # neuron.soma.v = random.uniform(-68, -55)
    # neuron.soma.v = -65
n.finitialize(-65) # Initialize V_m to -65 mV
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
fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=False)
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

v_all = np.array([list(neuron.v_vec) for neuron in neurons])
v_avg = np.mean(v_all, axis=0)

# Synchrony Index (k): Variance of average / Average of variances
sync_index = np.var(v_avg) / np.mean(np.var(v_all, axis=1))

v_normalized = v_avg - np.mean(v_avg)
# v_normalized = v_avg
N = len(v_normalized)
yf = np.fft.fft(v_normalized)
xf = np.fft.fftfreq(N, d=n.dt/1000)  # Convert dt from ms to s for frequency calculation
magnitude = np.abs(yf) * 2 / N

mask = (xf >= 1) & (xf <= 400) # Filters out the 0Hz spike and high noise
xf_plot = xf[mask]
# print(xf_plot)
mag_plot = magnitude[mask]
print(mag_plot)
# --- Subplot 3: 8. Plot Population Mean Membrane Potential ---
axes[2].plot(neurons[0].t_vec, v_avg, color='r')
axes[2].set_ylabel("Avg Vm (mV)")
axes[2].set_title(f"Population Mean (Sync Index: {sync_index:.3f})")
axes[2].set_xlabel("Time (ms)")
timestamp = datetime.datetime.now().strftime("%Y-%m-%d")
weight = syn_list[0][2].weight[0] if syn_list else 0
excWeight = all_stims[0][1].weight[0] if all_stims else 0
maxFreq = mag_plot[np.argmax(mag_plot)]
maxFreq = str(maxFreq)[:5]
# print(weight, excWeight)


plt.tight_layout(rect=[0, 0, 1, 0.96]) 
# plt.savefig(filename)

fig.suptitle(f'Network Simulation ({numNeurons} HH Neurons)', fontsize=16)
axes[3].plot(xf_plot, mag_plot, color='m')
axes[3].set_xlim(0, 400)
axes[3].set_ylabel("Amplitude")
axes[3].set_xlabel("Frequency (Hz)")
axes[3].set_title("FFT of Population Mean Membrane Potential")
filename = f"figures/simulation_results_{probability*100}%_{timestamp}-inhib+exc{weight}-{excWeight}maxFrequency{maxFreq}.png"
plt.savefig(filename)
# plt.grid(True)
plt.show()

print("\nResults saved to 'network_simulation_results.png'")
print(f"Total spikes detected across all neurons: {len(t_list)}")