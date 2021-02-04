# Updated LIF Models
#
# %% Imports
import numpy as np
import matplotlib.pyplot as plt
import pyNN.nest as sim 

# %% Functions

def spikes(start,end,period,duration=100,step=0.1):
    
    '''
    Returns an array of length duration x step.
    The array contains 0's where there are no spikes and 1's for spikes.
    A spike train is defined by its "start", "end" and "period". 
    '''

    spikes = np.zeros(int(duration/step))
    spikes[start,end,period] = 1

    return spikes

# %% LIF CUBA

def LIF_CUBA(spikes,duration=100,step=0.1,cm=1,tau_m=20,I_ext=0,
            v_rest=-65,v_reset=-65,tau_refrac=0.1,v_thresh=-50):

    ''' 
    Simulation of a LIF CUBA neuron.
    '''

    # Initializing Variables
    I_vec = []
    u_vec = []
    s_vec = []

    u  = v_rest
    g1 = Cm/tau_m

    for t in range(int(duration/step)):
        
        # Refractory Period
        if i in spikes[-tau_refrac:]:
            u = E_rest

        # Chaning Voltage
        else:
            du = ((g1*(E_rest - u)) + I_ext + I_syn) / cm
            u += du
        
        # Spiking
        if u >= v_thresh:
            spikes.append(1)
        else:
            spikes.append(0)
        
        uvec.append(u)

# %%











