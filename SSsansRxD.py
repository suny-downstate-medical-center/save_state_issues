from neuron import h, crxd as rxd
from neuron.crxd import v
from neuron.crxd.rxdmath import vtrap, exp, log
from math import pi
from matplotlib import pyplot
import numpy 
import os 
h.load_file('stdrun.hoc')


whichSave = 'savestate' # alternatively 'bbsavestate' 'manual'
outdir = 'saveStateDebug/save_state_mod/'
try:
    os.makedirs(outdir)
except:
    pass

pcid = 0
# pc = h.ParallelContext()
# pcid = pc.id()
# nhost =pc.nhost()
# pc.timeout(0)
# pc.set_maxstep(100)

tstop = 100

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
somaB.pt3dadd(-90,0,0,30)
somaB.pt3dadd(-60,0,0,30)
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

###### uncomment to add traditional cell ######
vvecB = h.Vector().record(somaB(0.5)._ref_v)
kvecB = h.Vector().record(somaB(0.5)._ref_ik)
navecB = h.Vector().record(somaB(0.5)._ref_ina)
mvecB = h.Vector().record(somaB(0.5).hh._ref_m)
nvecB = h.Vector().record(somaB(0.5).hh._ref_n)
hvecB = h.Vector().record(somaB(0.5).hh._ref_h)

h.finitialize(-70)

h.continuerun(70)

if whichSave == 'savestate':
    runSS()
    # saveRxd()
elif whichSave == 'bbsavestate':
    runBBSS()
elif whichSave == 'manual':
    saveRxd()
    saveVs()

h.continuerun(100)

# plotting 
fig = pyplot.figure()
pyplot.plot(tvec, vvecB, label="mod")
pyplot.xlabel('t (ms)')
pyplot.ylabel('V$_m$ (mV)')
pyplot.legend(frameon=False)
pyplot.savefig(os.path.join(outdir, 'fig1.png'))

fig = pyplot.figure()
pyplot.plot(tvec, hvecB, ':b', label='h')
pyplot.plot(tvec, mvecB, ':r', label='m')
pyplot.plot(tvec, nvecB, ':g', label='n')
pyplot.xlabel('t (ms)')
pyplot.ylabel('state')
pyplot.legend(frameon=False)
pyplot.savefig(os.path.join(outdir, 'fig2.png'))


fig = pyplot.figure()
pyplot.plot(tvec, kvecB, ':b', label='k')
pyplot.plot(tvec, navecB, ':r', label='na')
pyplot.xlabel('t (ms)')
pyplot.ylabel('current (mA/cm$^2$)')
pyplot.legend(frameon=False)
pyplot.savefig(os.path.join(outdir, 'fig3.png'))

out = {'kvec' : kvecB.as_numpy(),
    'navec' : navecB.as_numpy(),
    'tvec' : tvec.as_numpy(),
    'hvec' : hvecB.as_numpy(),
    'mvec' : mvecB.as_numpy(),
    'nvec' : nvecB.as_numpy(),
    'vvec' : vvecB.as_numpy()}

from scipy.io import savemat
savemat(os.path.join(outdir, 'data.mat'), out)

# v0.00 - remove traditional cell, save state at 30 ms 