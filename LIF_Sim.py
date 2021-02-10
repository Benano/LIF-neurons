# %% Simulation of LIF neuron interactions


# %% Classes

class sim:
    def __init__(self,increment,neurons=None): #input & output
        self.increment = increment

    neurons = []
    rec_var = []

    def add_neuron(self,neuron):
        self.neurons.append(neuron)

    def record(self,var):
        self.rec_var.append(var)
        
    def run(self,duration):

        # MS conversion
        steps = int(duration/self.increment)

        for i in range(steps):
            for neuron in neurons:
                neuron.update(self.increment)


class neuron():
    def __init__(self,tau_m=20,cm=1,tau_syn=5,E_rest=-65,v_reset=-65,tau_refrac=0.1,v_thresh=-50, I_ext=0): # Attributes shared by all
        self.tau_m = tau_m
        self.cm = cm
        self.tau_syn = tau_syn
        self.E_rest = E_rest
        self.v_reset = v_reset
        self.tau_refrac = tau_refrac
        self.v_thresh = v_thresh
        self.I_ext = I_ext
        self.u = E_rest
        self.g1 = self.cm/self.tau_m

    def update(self,increment):
        du = (((self.g1*(self.E_rest - self.u)) + self.I_ext) / self.cm)*self.step
        self.u += du
        
        # Spike
        if self.u >= self.v_thresh:
            self.u = self.v_reset


sim1 = sim(100,0.1)
n1 = neuron()
sim1.add_neuron(n1)
sim1.record(v)
sim1.run(10)

# %% 
class lif_cuba(neuron):
    def __init__(self,tau_m=20,cm=1,tau_syn=5,E_rest=-65,v_reset=-65,tau_refrac=0.1,v_thresh=-50):
        super().__init__(tau_m,cm,tau_syn,E_rest,v_reset,tau_refrac,v_thresh)


    def inject_current(self,strength,start,duration):

        if self.kind == "lif_cond":
            if (start + duration) > self 
            
            current_inj = np.arange()

            self.I_ext = stren

    def run(self):
        for i in range(self.steps):
            for neuron in self.neurons:
                neuron.update()

                yield neuron

    def record()