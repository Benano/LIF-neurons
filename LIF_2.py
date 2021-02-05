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
    start = int(start/step)
    end = int(end/step)
    period = (int(period/step))
    spikes[start:end:period] = 1

    return spikes

# %% LIF CUBA

def LIF_CUBA(spikes=np.zeros(1),duration=100,step=0.1,cm=1,tau_m=20,I_ext=0,
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
    
# %% LIF COBA

def LIF_COBA(spikes_e=np.zeros(1),spikes_i=np.zeros(1),duration=100,step=0.1,cm=1,tau_m=20,
            E_rest=-65,E_rev_e=0,E_rev_i=-70,tau_syn_e=20,tau_syn_i=20,v_reset=-65,tau_refrac=0.1,v_thresh=-50):

    ''' 
    Simulation of a LIF COBA neuron.
    '''

    # Initializing Variables
    ge_vec = []
    gi_vec = []
    u_vec = []
    s_vec = []
    
    u  = E_rest
    g1 = cm/tau_m
    I_syn = 0
    g_syn_e = 0
    g_syn_i = 0

    for t in range(int(duration/step)):
        
        # Incoming Excitatory Spikes
        dg_syn_e = (-g_syn_e/tau_syn_e)*step
        g_syn_e += dg_syn_e
        
        if np.any(spikes_e):
            if spikes_e[t]:
                g_syn_e += cm/tau_syn_e 
            
        # Incoming Inhibitory Spikes
        dg_syn_i = (-g_syn_i/tau_syn_i)*step
        g_syn_i += dg_syn_i
        
        if np.any(spikes_i):
            if spikes_i[t]:
                g_syn_i += 0.05

        # Refractory Period
        if 1 in s_vec[-(int(tau_refrac/step)):]:
            u = E_rest

        # Changing Voltage
        else:
            du = (((g1*(E_rest - u)) + g_syn_e*(E_rev_e-u) + g_syn_i*(E_rev_i-u)) / cm)*step
            u += du
        
        # Spiking
        if u >= v_thresh:
            s_vec.append(1)
            u = v_reset
        else:
            s_vec.append(0)
        
        u_vec.append(u)
        ge_vec.append(g_syn_e)
        gi_vec.append(g_syn_i)

    return ge_vec, gi_vec, u_vec, s_vec
    
# %% Simulation LIF COBA

spikes_e = make_spikes(20,22,3)
spikes_i = make_spikes(70,72,3)
ge, gi, u, sp  = LIF_COBA(spikes_e=spikes_e,tau_syn_e=5)


#%% # Plotting 
plt.figure()
plt.title("g_e")
plt.plot(ge)

plt.figure()
plt.title("g_i")
plt.plot(gi)

plt.figure()
plt.title("voltage")
plt.plot(u)

plt.figure()
plt.title("spikes")
plt.plot(sp)


# %% PyNN CUBA

# Setup
sim.setup(timestep=0.1, min_delay=0.1, max_delay=10.0)
IF_sim = sim.Population(1, sim.IF_curr_exp(i_offset=1.0), label="IF_curr_exp")
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

# %% PyNN COBA

# Setup
sim.setup(timestep=0.1, min_delay=0.1, max_delay=10.0)
IF_sim = sim.Population(1, sim.IF_cond_exp(), label="IF_curr_exp")
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
