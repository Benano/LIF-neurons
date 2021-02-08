# %% Simulation of LIF neuron interactions


# %% Classes

class sim:
    def __init__(self,runtime,increment):
        self.runtime = runtime
        self.increment = increment
        self.steps = int(runtime/increment)

    def add_neuron(self,kind,tau_m=20,cm=1,tau_syn=5,E_rest=-65,v_reset=-65,tau_refrac=0.1,v_thresh=-50):
        self.kind = kind
        self.tau_m = tau_m
        self.cm = cm
        self.tau_syn = tau_syn
        self.E_rest = E_rest
        self.v_reset = v_reset
        self.tau_refrac = tau_refrac
        self.v_thresh = v_thresh

   
sim = sim(100,0.1)
sim.add_neuron("CUBA")
# %% 
        if kind == "lif_cond":
            self.I_ext = 0
    
    def update_neuron(self):
    def inject_current(self,strength,start,duration):

        if self.kind == "lif_cond":
            if (start + duration) > self 
            
            current_inj = np.arange()

            self.I_ext = stren

    def run(self):
        for i in range(self.steps):
            print("I am happy")

    def add_neuron()

    def update_neuron()

class synapse:

    def run():
    for neuron in ...:
        update ()

        yield ...

