from neuron import h, crxd as rxd
from neuron.crxd import v
from matplotlib import pyplot
import numpy 
import os 
h.load_file('stdrun.hoc')


whichSave = 'savestate' # alternatively 'bbsavestate' 'manual'
outdir = 'saveStateDebug/mpi_restore/' # 'saveStateDebug/save_state_mod_mpi/' #
restoredir = 'saveStateDebug/save_state_mod_mpi/' # None # assign a path to restore sim from saved state 

try:
    os.makedirs(outdir)
except:
    pass

pc = h.ParallelContext()
pcid = pc.id()
nhost =pc.nhost()
pc.timeout(0)
pc.set_maxstep(100)

def saveRxd():
    for sp in rxd.species._all_species:
        s = sp()
        numpy.save(os.path.join(outdir, s.name + '_concentrations_' + str(pcid) + '.npy'), s.nodes.concentration)

def saveVs():
    all_v = []
    for sec in h.allsec():
        for seg in sec.allseg():
            all_v.append(seg.v)
    numpy.save(os.path.join(outdir, 'pcid'+str(pcid)+'_allv.npy'), all_v)
    
def runSS():
    svst = h.SaveState()
    svst.save()
    f = h.File(os.path.join(outdir,'save_test_' + str(pcid) + '.dat'))
    svst.fwrite(f)

def runBBSS():
    svst = h.BBSaveState()
    svst.save_test()

# parameters
h.celsius = 6.3 

###### uncomment to add traditional cell ######
somaB = h.Section('somaB')
somaB.pt3dclear()
somaB.pt3dadd(90,0,0,30)
somaB.pt3dadd(60,0,0,30)
somaB.nseg = 11

###### uncomment to add traditional cell ######
somaB.insert('hh')

###### uncomment to add traditional cell ######
stimB = h.IClamp(somaB(0.5))
stimB.delay = 50
stimB.amp = 1
stimB.dur = 50

# record
tvec = h.Vector().record(h._ref_t)
vvecB = h.Vector().record(somaB(0.5)._ref_v)
kvecB = h.Vector().record(somaB(0.5)._ref_ik)
navecB = h.Vector().record(somaB(0.5)._ref_ina)
mvecB = h.Vector().record(somaB(0.5).hh._ref_m)
nvecB = h.Vector().record(somaB(0.5).hh._ref_n)
hvecB = h.Vector().record(somaB(0.5).hh._ref_h)

if restoredir:
    def restoreSS():
        svst = h.SaveState()
        f = h.File(os.path.join(restoredir, 'save_test_'+str(pcid) + '.dat'))
        svst.fread(f)
        svst.restore()

    def restoreSim():
        restoreSS()
        # for sp in rxd.species._all_species:
        #     s = sp()
        #     print(s.name)
        #     temp = numpy.load(os.path.join(restoredir, s.name + '_concentrations_' + str(pcid) + '.npy'))
        #     s.nodes.concentration = list(temp)
        #     print('PCID ' + str(pcid) + ': resotred ' + s.name)
    fih = fih = h.FInitializeHandler(1, restoreSS)
    h.finitialize()
    tstop = 100
else:
    h.finitialize(-70)
    tstop = 70

while h.t < tstop:
    pc.psolve(pc.t(0)+h.dt)

if whichSave == 'savestate':
    runSS()
    # saveRxd()
elif whichSave == 'bbsavestate':
    runBBSS()
elif whichSave == 'manual':
    saveRxd()
    saveVs()

out = {'kvec' : kvecB.as_numpy(),
    'navec' : navecB.as_numpy(),
    'tvec' : tvec.as_numpy(),
    'hvec' : hvecB.as_numpy(),
    'mvec' : mvecB.as_numpy(),
    'nvec' : nvecB.as_numpy(),
    'vvec' : vvecB.as_numpy()}

from scipy.io import savemat
savemat(os.path.join(outdir, 'data_' + str(pcid) + '.mat'), out)

# v0.00 - remove traditional cell, save state at 30 ms 