# from mpi4py import MPI
from neuron import h 
h.load_file('stdrun.hoc')

pc = h.ParallelContext()
pcid = pc.id()
nhost =pc.nhost()
pc.timeout(0)
pc.set_maxstep(100)

print("I am " + str(pcid))

h.quit()