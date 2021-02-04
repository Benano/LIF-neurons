# %% Imports 

import numpy as np
import matplotlib.pyplot as plt
import random as rd
import pyNN.nest as sim

# %% Functions 
#Injecting Current
def inject(s,a,b,n):

    Iex = np.zeros(N)
    Iex[a:b] = s

    if n:
        noise = np.random.normal(0,1,N)
        Iex = Iex + noise

    return Iex
# %% Eulers Parameters
b = 100
N = 1000
a = 0
h = (b-a)/N
f = int(N/b)

# Plotting window
w  = 100 # window in ms
wf = w*f # window in time points
xs = np.linspace(a,w,w*f,endpoint=False)

# %% Single LIF
# Current
tau = 20*f 
Cm  = 1 
Er  = -65 
v_thresh = -50
g1 = Cm/tau
I  = 1/f 
tau_ref = int(0.1*f)

spikes = []
uvec = np.zeros(N)
u = Er
for i in range(N):

    # Refractory Time
    if 1 in spikes[-tau_ref:]:
        uvec[i] = u = Er
    
    # Change in u
    else:
        cu = (g1*(Er-u) + I)/Cm
        u = u + cu
        uvec[i] = u

    # Spiking
    if u >= v_thresh:
        spikes.append(1)
    else:
        spikes.append(0)

# Plot
plt.plot(xs,uvec[:wf])
plt.xlabel("ms")
plt.ylabel("mV")
plt.title("Single LIF Neuron (resting)")

# %% Single LIF Neuron
# Resting
tau = 20*f 
Cm  = 1 
Er  = -65 
v_thresh = 0
g1 = Cm/tau
I  = 1/f 
eff_Er = Er + I/g1
tau_ref = int(0.1*f)

spikes = []
uvec = np.zeros(N)
u = Er
for i in range(N):

    # Refractory Time
    if 1 in spikes[-tau_ref:]:
        uvec[i] = u = Er
    
    # Change in u
    else:
        cu = (g1*(eff_Er-u))/Cm
        u = u + cu
        uvec[i] = u

    # Spiking
    if u >= v_thresh:
        spikes.append(1)
    else:
        spikes.append(0)

# Plot
plt.plot(xs,uvec[:wf])
plt.xlabel("ms")
plt.ylabel("mV")
plt.title("Single LIF Neuron (current)")

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


# %% Two LIF neurons

# Parameters
Cm = 1
tau = 20
E1 = -65 
g1 = Cm/tau
I1 = 1.5
I2 = 0.8
tau_ref = 3
kick = 40


# Neuron 1
uvec1 = np.zeros(N)
spikes1 = []
Iex1 = inject(0,40,100,0)

# Synapse
w0 = 0.2 
mu = 0.05

# Neuron 2
uvec2 = np.zeros(N)
spikes2 = []
Iex2 = inject(0,40,100,0)


# Loop
u1 = -80
u2 = -70
for i in range(N):

    # Neuron 1
    if 1 in spikes1[-tau_ref:]:
        uvec1[i] = u1 = E1
    
    else:
        cu1 = g1*(E1-u1) + I1 + Iex1[i]
        u1 = u1 + cu1
        uvec1[i] = u1

    if u1 >= -v_thresh:
        spikes1.append(1)
    else:
        spikes1.append(0)

    # Neuron 2
    if 1 in spikes2[-tau_ref:]:
        uvec2[i] = u2 = -65
    
    else:
        if spikes1[i]:
            sp = w0*kick
        else:
            sp = 0

        cu2 = g1*(E1-u2) + I2 + Iex2[i] + sp
        u2 = u2 + cu2
        uvec2[i] = u2

    if u2 >= -40:
        spikes2.append(1)
        if sp:
            w0 = w0 + mu
    
    else:
        spikes2.append(0)

plt.plot(uvec1[0:400])
plt.plot(uvec2[0:400])




