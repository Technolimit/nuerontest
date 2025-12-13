from neuron import n
from neuron.units import ms, mV, µm
import matplotlib.pyplot as plt
n.load_file("stdrun.hoc")

soma = n.Section("soma")
soma.L = 20
soma.diam = 20
soma.Ra=100
soma.cm=1
soma.insert("hh")
for seg in soma:
            seg.hh.gnabar = 0.12  # Sodium conductance in S/cm2
            seg.hh.gkbar = 0.036  # Potassium conductance in S/cm2
            seg.hh.gl = 0.0003  # Leak conductance in S/cm2
            seg.hh.el = -54.3 * mV  # Reversal potential
stim = n.IClamp(soma(0.5))
stim.delay = 100 * ms
stim.dur = 500 * ms
stim.amp = 0.1  # nA
n.tstop = 700 * ms
v_vec = n.Vector()
t_vec = n.Vector()
v_vec.record(soma(0.5)._ref_v)
t_vec.record(n._ref_t)


def detect_spikes():
    threshold = 0
    pastThreshold = False
    spikes = 0
    for i in v_vec:
        if (pastThreshold == False and i>threshold):
            spikes += 1
            pastThreshold = True
        if (pastThreshold==True and i<threshold):
            pastThreshold = False
    print("Number of spikes:", spikes)
    return spikes

# stim2 = n.NetStim()
# syn_ = n.ExpSyn(soma(0.5))
# stim2.number = 1
# stim2.start =200 * ms
# ncstim = n.NetCon(stim2, syn_)
# ncstim.delay = 1 * ms
# ncstim.weight[0] = 0.04
# syn_.tau = 2 * ms

# nc = n.NetCon(source.soma(0.5)._ref_v, syn, sec=source.soma)
# nc.weight[0] = 0.05
# nc.delay = 5
spikeList = []
iAmp = []
for i in range(10):
    stim.amp=i*0.02
    n.run()
    spikeList.append(detect_spikes())
    iAmp.append(stim.amp)
plt.figure()
plt.plot(iAmp, spikeList)
plt.xlabel("Input Current (nA)")
plt.ylabel("Number of Spikes")
plt.title("F-I Curve")
plt.figure()
plt.plot(t_vec, v_vec)
plt.xlabel("Time (ms)")
plt.ylabel("Membrane Potential (mV)")
plt.title("Membrane Potential vs Time")
plt.show()

# class BallAndStick:
#     def __init__(self):
#         self.soma = n.Section("soma", self)
#     # def __repr__(self):
#     #     return "BallAndStick[{}]".format(self._gid)
# my_cell = BallAndStick()
# my_other_cell = BallAndStick()
n.topology()