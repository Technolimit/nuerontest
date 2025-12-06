from neuron import h, init
import matplotlib.pyplot as plt
soma = h.Section(name='soma')

soma.L = 20.0
soma.diam = 20.0
soma.nseg = 1

soma.cm = 1.0
soma.Ra = 100.0

soma.insert('hh')

soma.gnabar_hh = 0.12
soma.gkbar_hh = 0.036
soma.el_hh = -54.3
soma.gl_hh = 0.0003

stim = h.IClamp(soma(0.5))
stim.delay = 100.0
stim.dur = 50.0
stim.amp = 0.1

# Recording setup
t_vec = h.Vector()
v_vec = h.Vector()
t_vec.record(h._ref_t)
v_vec.record(soma(0.5)._ref_v)

# Simulation commands
h.tstop = 300.0
h.dt = 0.025
h.v_init = -65.0

h.finitialize(h.v_init)
h.continuerun(h.tstop)



plt.figure(figsize=(10, 4))
plt.plot(t_vec, v_vec)
plt.title('Membrane Potential vs. Time')
plt.xlabel('Time (ms)')
plt.ylabel('Membrane Potential (mV)')
plt.grid(True)
plt.show()

print(f"Simulation finished. Initial V: {h.v_init} mV. Final V: {v_vec.to_python()[-1]:.2f} mV.")