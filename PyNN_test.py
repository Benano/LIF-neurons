#%% Imports 
import numpy as np
import matplotlib.pyplot as plt
import pyNN.nest as sim

# %% PYNN Neuron

sim.setup(timestep=0.1, min_delay=0.1, max_delay=10.0)
IF_sim = sim.Population(1, sim.IF_curr_alpha(i_offset=1.0), label="IF_curr_exp")
IF_sim.record('v')

## Parameters
# {'v_rest': -65.0, 'cm': 1.0, 'tau_m': 20.0, 
# 'tau_refrac': 0.1, 'tau_syn_E': 0.5, 'tau_syn_I': 0.5, 
# 'i_offset': 1.0, 'v_reset': -65.0, 'v_thresh': -50.0}

# Running simulation in MS
sim.run(100.0)

# Data 
v_data = IF_sim.get_data()   
data = IF_sim.get_data().segments[0]
v = data.filter(name="v")[0]
plt.plot(v)
# End
sim.end()
