# %% Simulation of LIF
import matplotlib.pyplot as plt

# %% Classes
class sim:
    def __init__(self,increment): #input & output
        self.increment = increment

    neuron_labels = []
    neurons = []
    synapses = []

    def add_neuron(self,neuron,label=None):
        self.neurons.append(neuron)
        self.synapses.append([])

        if label:
            self.neuron_labels = label
        else:
            self.neuron_labels = "neuron_" + str(len(self.neuron))

    def add_synapse(self,pre,post):
        self.synapses[post].append(pre)
    
    def current_sim(self):
        print(self.neuron_labels)
        
    def run(self,duration):
        # MS conversion
        steps = int(duration/self.increment)
        for i in range(steps):
            for idx,neuron in enumerate(sim.neurons):
                neuron.update(self.increment,spike_input)

                if i:
                    if synapses[i]:
                        for s in range(synapses):

            

class neuron():
    def __init__(self,tau_m=20,cm=1,tau_syn=5,E_rest=-65,v_reset=-65,tau_refrac=0.1,v_thresh=-50): # Attributes shared by all
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

    def update(self,increment,spike_input=None):
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

        self.I_ext = I_ext
        self.u_ts = []
        self.s_ts = []

    def update(self,increment,spike_input=None):
        # Voltage
        du = (((self.g1*(self.E_rest - self.u)) + self.I_ext) / self.cm)* increment
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

class spike_array(sim):
    def __init__(self,spike_times):
        super().__init__(increment)

        self.spike_times = np.array(spike_times)/self.increment






# %% 
sim1 = sim(0.1)
n1 = lif_cuba(I_ext=1)
n2 = lif_cuba(I_ext=1.5)

sim1.add_neuron(n1,"Lif_Coba_1")
sim1.add_neuron(n2,"Lif_Coba_2")
sim1.add_synapse(0,1)
sim1.current_sim()

sim1.run(100)

plt.plot(n1.u_ts)

print(sim.synapses)


