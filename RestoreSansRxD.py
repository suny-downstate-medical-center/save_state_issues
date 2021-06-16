from neuron import h, crxd as rxd
from neuron.crxd import v
from neuron.crxd.rxdmath import vtrap, exp, log
from math import pi
from matplotlib import pyplot
# pyplot.ion()
import numpy 
import os 
h.load_file('stdrun.hoc')

outdir = 'saveStateDebug/restore_run70_ssupdate/'
restoredir = 'saveStateDebug/save_state_mod/'

try:
    os.makedirs(outdir)
except:
    pass

# pc = h.ParallelContext()
# pcid = pc.id()
# nhost =pc.nhost()
# pc.timeout(0)
# pc.set_maxstep(100)
pcid = 0

endt = 70
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

def restoreSim():
    # dat_v = numpy.load(os.path.join(restoredir, 'pcid'+str(pcid)+'_allv.npy'))
    # ind = 0
    # for sec in h.allsec():
    #     for seg in sec.allseg():
    #         seg.v = dat_v[ind]
    #         ind = ind + 1
    # print('Restored voltages from pcid'+str(pcid)+'_allv.npy')
    
    for sp in rxd.species._all_species:
        s = sp()
        print(s.name)
        temp = numpy.load(os.path.join(restoredir, s.name + '_concentrations_' + str(pcid) + '.npy'))
        s.nodes.concentration = list(temp)
        print('PCID ' + str(pcid) + ': resotred ' + s.name)
    restoreSS()

    # h.t = endt

def restoreSS():
    svst = h.SaveState()
    f = h.File(os.path.join(restoredir, 'save_test_'+str(pcid) + '.dat'))
    svst.fread(f)
    svst.restore()

# parameters
h.celsius = 6.3 

somaB = h.Section('somaB')
somaB.pt3dclear()
somaB.pt3dadd(60,0,0,30)
somaB.pt3dadd(90,0,0,30)
somaB.nseg = 11

somaB.insert('hh')

# stimulate
stimB = h.IClamp(somaB(0.5))
stimB.delay = 50
stimB.amp = 1
stimB.dur = 50

# record
vvecB = h.Vector().record(somaB(0.5)._ref_v)
kvecB = h.Vector().record(somaB(0.5)._ref_ik)
navecB = h.Vector().record(somaB(0.5)._ref_ina)
mvecB = h.Vector().record(somaB(0.5).hh._ref_m)
nvecB = h.Vector().record(somaB(0.5).hh._ref_n)
hvecB = h.Vector().record(somaB(0.5).hh._ref_h)
tvec = h.Vector().record(h._ref_t)

fih = h.FInitializeHandler(1, restoreSS)
h.finitialize()

# h.continuerun(100)

# while h.t < tstop:
#     # pc.psolve(h.t + h.dt)
#     h.fadvance()
#     print(h.t)

# # plotting 
# fig = pyplot.figure()
# pyplot.plot(tvec, vvecA, label="rxd")
# # pyplot.plot(tvec, vvecB, label="mod")
# pyplot.xlabel('t (ms)')
# pyplot.ylabel('V$_m$ (mV)')
# pyplot.legend(frameon=False)
# # pyplot.savefig(os.path.join(outdir, 'fig1.png'))

# fig = pyplot.figure()
# pyplot.plot(tvec, hvecA, '-b', label='h')
# pyplot.plot(tvec, mvecA, '-r', label='m')
# pyplot.plot(tvec, nvecA, '-g', label='n')
# # pyplot.plot(tvec, hvecB, ':b')
# # pyplot.plot(tvec, mvecB, ':r')
# # pyplot.plot(tvec, nvecB, ':g')
# pyplot.xlabel('t (ms)')
# pyplot.ylabel('state')
# pyplot.legend(frameon=False)
# # pyplot.savefig(os.path.join(outdir, 'fig2.png'))


# fig = pyplot.figure()
# pyplot.plot(tvec, kvecA.as_numpy(), '-b', label='k')
# pyplot.plot(tvec, navecA.as_numpy(), '-r', label='na')
# # pyplot.plot(tvec, kvecB, ':b')
# # pyplot.plot(tvec, navecB, ':r')
# pyplot.xlabel('t (ms)')
# pyplot.ylabel('current (mA/cm$^2$)')
# pyplot.legend(frameon=False)
# # pyplot.savefig(os.path.join(outdir, 'fig3.png'))

# out = {'kvec' : kvecA.as_numpy(),
#     'navec' : navecA.as_numpy(),
#     'tvec' : tvec.as_numpy(),
#     'hvec' : hvecA.as_numpy(),
#     'mvec' : mvecA.as_numpy(),
#     'nvec' : nvecA.as_numpy(),
#     'vvec' : vvecA.as_numpy()}

# from scipy.io import savemat
# savemat(os.path.join(outdir, 'data.mat'), out)