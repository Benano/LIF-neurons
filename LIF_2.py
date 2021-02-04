# Updated LIF Models
#
# %% Imports
import numpy as np
import matplotlib.pyplot as plt
import pyNN.nest as sim 

# %% Functions

def make_spikes(start,end,period,duration=100,step=0.1):
    
    '''
    Returns an array of length duration x step.
    The array contains 0's where there are no spikes and 1's for spikes.
    A spike train is defined by its "start", "end" and "period". 
    '''

    spikes = np.zeros(int(duration/step))
    spikes[start:end:period] = 1

    return spikes

# %% LIF CUBA

def LIF_CUBA(spikes=np.empty(1),duration=100,step=0.1,cm=1,tau_m=20,I_ext=0,
            E_rest=-65,v_reset=-65,tau_refrac=0.1,v_thresh=-50):

    ''' 
    Simulation of a LIF CUBA neuron.
    '''

    # Initializing Variables
    I_vec = []
    u_vec = []
    s_vec = []
    
    u  = E_rest
    g1 = cm/tau_m
    I_syn = 0

    for t in range(int(duration/step)):
        
        # Incoming Spikes
        di_syn = (-I_syn/tau_m)*step
        I_syn  += di_syn
        
        if np.any(spikes):
            if spikes[t]:
                I_syn += 1 
            
        # Refractory Period
        if 1 in s_vec[-(int(tau_refrac/step)):]:
            u = E_rest

        # Changing Voltage
        else:
            du = (((g1*(E_rest - u)) + I_ext + I_syn) / cm)*step
            u += du
        
        # Spiking
        if u >= v_thresh:
            s_vec.append(1)
            u = v_reset
        else:
            s_vec.append(0)
        
        u_vec.append(u)
        I_vec.append(I_syn)

    return I_vec, u_vec, s_vec
    
# %% Simulation LIF CUBA

spikes = make_spikes(40,50,3)
I_cuba, u_cuba, s_cuba = LIF_CUBA(I_ext=1)

# Plotting
plt.figure()
plt.title("Current")
plt.plot(I_cuba)

plt.figure()
plt.title("Membrane Voltage")
plt.plot(u_cuba)

plt.figure()
plt.title("Spikes")
plt.plot(s_cuba)

# %% PyNN CUBA

# Setup
sim.setup(timestep=0.1, min_delay=0.1, max_delay=10.0)
IF_sim = sim.Population(1, sim.IF_curr_alpha(i_offset=1.0), label="IF_curr_exp")
IF_sim.record('v')

# Running simulation in MS
sim.run(100.0)

# Data 
v_data = IF_sim.get_data()   
data = IF_sim.get_data().segments[0]
v = data.filter(name="v")[0]

# Plotting
plt.plot(v)

# End
sim.end()
# %%
