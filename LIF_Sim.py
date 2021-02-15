# %% Simulation of LIF
import matplotlib.pyplot as plt
import numpy as np

# %% Classes
class sim:
    def __init__(self,increment):
        self.increment = increment

    neuron_labels = []
    neurons = []
    synapses = []

    def add_neuron(self,neuron,label=None):
        self.neurons.append(neuron)

        if label:
            self.neuron_labels.append(label)
        else:
            self.neuron_labels.append("neuron_" + str(len(self.neuron)))

    def add_synapse(self,synapse):
        self.synapses.append(synapse)
    
    def current_sim(self):
        print(self.neuron_labels)
        
    def run(self,duration):
        # MS conversion
        steps = int(duration/self.increment)

        for i in range(steps):

            for neuron in sim.neurons:
                neuron.update(self.increment)

            for synapse in sim.synapses:
                if synapse.pre.s_ts[i]:
                    synapse.post.s_input =+ synapse.weight

class neuron():
    def __init__(self,tau_m=20,cm=1,tau_syn=5,E_rest=-65,v_reset=-65,tau_refrac=0.1,v_thresh=-50):
        self.tau_m = tau_m
        self.cm = cm
        self.tau_syn = tau_syn
        self.E_rest = E_rest
        self.v_reset = v_reset
        self.tau_refrac = tau_refrac
        self.v_thresh = v_thresh
        self.u = v_reset
        self.g1 = self.cm/self.tau_m

class lif_coba(neuron):
    def __init__(self,tau_m=20,cm=1,tau_syn=5,E_rest=-65,v_reset=-65,tau_refrac=0.1,v_thresh=-50,
                E_rev_e=0, E_rev_i=-70, g_syn_e=5, g_syn_i=5):
        super().__init__(tau_m,cm,tau_syn,E_rest,v_reset,tau_refrac,v_thresh)
        self.E_rev_e = E_rev_e
        self.E_rev_i = E_rev_i
        self.g_syn_e = g_syn_e
        self.g_syn_i = g_syn_i

        self.ge_ts = []
        self.gi_ts = []
        self.u_ts = []
        self.s_ts = []

    def update(self,increment):
        # Voltage
        du = (((self.g1*(self.E_rest - self.u)) + self.g_syn_e*(self.E_rev_e-self.u) + self.g_syn_i*(self.E_rev_i-self.u)) / self.cm)* increment
        self.u += du

        # Spike
        if self.u >= self.v_thresh:
            self.u = self.v_reset
            self.s_ts.append(1)
        else:
            self.s_ts.append(0)
        
        # Refractory Period
        if 1 in self.s_ts[-(int(self.tau_refrac/increment)):]:
            self.u = self.v_reset

        self.u_ts.append(self.u)

class lif_cuba(neuron):
    def __init__(self,tau_m=20,cm=1,tau_syn=5,E_rest=-65,v_reset=-65,tau_refrac=0.1,v_thresh=-50,
                E_syn_e=0, E_syn_i=-70, g_syn_e=5, g_syn_i=5, I_ext=0):
        super().__init__(tau_m,cm,tau_syn,E_rest,v_reset,tau_refrac,v_thresh)

        self.E_syn_e = E_syn_e
        self.E_syn_i = E_syn_i
        self.g_syn_e = g_syn_e
        self.g_syn_i = g_syn_i
        self.I_ext = I_ext

        self.u_ts = []
        self.s_ts = []
        self.I_syn = 0
        self.s_input = 0

    def update(self,increment):

        # Input
        di_syn = (-self.I_syn/self.tau_syn)*increment
        self.I_syn += di_syn 
        self.I_syn += self.s_input
        self.s_input = 0

        # Voltage
        du = (((self.g1*(self.E_rest - self.u)) + self.I_ext + self.I_syn) / self.cm)* increment
        self.u += du
        
        # Spike
        if self.u >= self.v_thresh:
            self.u = self.v_reset
            self.s_ts.append(1)
        else:
            self.s_ts.append(0)
        
        # Refractory Period
        if 1 in self.s_ts[-(int(self.tau_refrac/increment)):]:
            u = self.v_reset
        
        self.u_ts.append(self.u)

class spike_array():
    def __init__(self,spike_times,increment):
        self.spike_times = np.array(spike_times)/increment

        self.time = 0
        self.s_ts = []

    def update(self,increment):
        self.time += 1

        if self.time in self.spike_times:
            self.s_ts.append(1)
        else:
            self.s_ts.append(0)

class synapse():
    def __init__(self,pre,post,weight):
        self.pre = pre
        self.post = post
        self.weight = weight

# %% 
# Simulation
step = 0.1
sim1 = sim(step)

# Neuron
n1 = lif_cuba(I_ext=0)
sim1.add_neuron(n1,"Lif_Coba_1")

# Spiky Input (ms)
spike_times = [10]
s1 = spike_array(spike_times,step)
sim1.add_neuron(s1,"Spike_Input")

# Connection
sim1.add_synapse(synapse(s1,n1,1))

# Summary
print(sim1.neuron_labels)

# Running
sim1.run(100)

# Plotting
plt.plot(n1.u_ts)

# %%
