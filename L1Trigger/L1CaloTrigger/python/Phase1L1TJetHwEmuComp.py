#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import math
import numpy
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors

ptLSB = 0.25;
etaLSB = 0.0043633231;
phiLSB = 0.0043633231;
hwJets = []
emJets = []
hwData = []
emData = []
hwDataNoZ = []
emDataNoZ = []
nHw = 0
nEm = 0
nEv = 5000
nEvNoZ = 0


for hwRx in xrange(0,70):
    with open("datFull/" + str(hwRx) + "/tx_summary.txt", "r") as inFile:
        frameIt = -1
        for line in inFile:
            if('1v' in line):
                frameIt += 1
                if frameIt < 20:
                    continue
                if hwRx == 0 and frameIt < 33:
                    continue
                linkData = line.split('1v')
                for wordIt in xrange(1,25):
                    word = linkData[wordIt].replace(' ','').replace('\n','')
                    if int(word, 16) & 0xffff:
                        jet = word[8:]
                        hwJets.append([(int(jet,16)&0xffff)*ptLSB, 
                                         ((((int(jet,16)>>24)&0xff)*19)+9)*etaLSB,
                                         ((((int(jet,16)>>16)&0xff)*20)+10)*phiLSB])
                    if (int(word, 16)>>32) & 0xffff:
                        jet = word[:8]
                        hwJets.append([(int(jet,16)&0xffff)*ptLSB, 
                                         ((((int(jet,16)>>24)&0xff)*19)+9)*etaLSB,
                                         ((((int(jet,16)>>16)&0xff)*20)+10)*phiLSB])
                if (frameIt%13) == 6:
                    if(nHw>=nEv):
                        break
                    nHw+=1
                    if len(hwJets)==0:
                        hwJets.append([0,0,0])
                    hwData.append(hwJets)
                    del hwJets
                    hwJets = []


with open("emuout.txt", "r") as inFile:
    for line in inFile:
        if " " in line:
            if(nEm>=nEv):
                break
            nEm+=1
            if len(emJets)>0:
                emData.append(emJets)
            del emJets
            emJets = []
        else:
            jet = [float(line.split("\t")[0]),
                   float(line.split("\t")[1]),
                   float(line.split("\t")[2])]
            emJets.append(jet)




print("=====================================================================================")
print("\t\tFirmware Events: " + str(nHw) + "\t\t" + "Emulator Events: " + str(nEm))
print("=====================================================================================")
print("\t\tpT\t" + "eta\t" + "phi\t\t" + "pT\t" + "eta\t" + "phi\t")
print("=====================================================================================")


for evIt in xrange(0,nEv):
    if hwData[evIt][0][0] > 0:
        hwDataNoZ.append(hwData[evIt])
    if emData[evIt][0][0] > 0:
        emDataNoZ.append(emData[evIt])
    nEvNoZ+=1


for evIt in xrange(0,nEv):
    if hwData[evIt][0][0] ==0 and emData[evIt][0][0] == 0:
        continue
    jetCount=0
    jetDiff = len(hwData[evIt]) - len(emData[evIt])
    print("")
    if jetDiff==0:
        for jetIt in xrange(len(hwData[evIt])):
            print(str(evIt) + "\t\t" + str(hwData[evIt][jetIt][0]) + "\t" + str(hwData[evIt][jetIt][1])[:4] + "\t" + str(hwData[evIt][jetIt][2])[:4] + "\t\t" +
                  str(emData[evIt][jetIt][0]) + "\t" + str(emData[evIt][jetIt][1])[:4] + "\t" + str(emData[evIt][jetIt][2])[:4])
    if jetDiff>0:
        for jetIt in xrange(len(hwData[evIt])):
            jetCount+=1
            if jetCount > len(emData[evIt]):
                emData[evIt].append([0,0,0])
            print(str(evIt) + "\t\t" + str(hwData[evIt][jetIt][0]) + "\t" + str(hwData[evIt][jetIt][1])[:4] + "\t" + str(hwData[evIt][jetIt][2])[:4]  + "\t\t" +
                  str(emData[evIt][jetIt][0]) + "\t" + str(emData[evIt][jetIt][1])[:4] + "\t" + str(emData[evIt][jetIt][2])[:4])
    if jetDiff<0:
        for jetIt in xrange(len(emData[evIt])):
            jetCount+=1
            if jetCount > len(hwData[evIt]):
                hwData[evIt].append([0,0,0])
            print(str(evIt) + "\t\t" + str(hwData[evIt][jetIt][0]) + "\t" + str(hwData[evIt][jetIt][1])[:4] + "\t" + str(hwData[evIt][jetIt][2])[:4]  + "\t\t" +
                  str(emData[evIt][jetIt][0]) + "\t" + str(emData[evIt][jetIt][1])[:4] + "\t" + str(emData[evIt][jetIt][2])[:4])
        



fig, axs =   plt.subplots(2,3, figsize=(20, 10), gridspec_kw={'height_ratios': [3, 1]})

fig.patch.set_facecolor( '#ffffff')


nPtHw  = axs[0,0].hist([jet[0] for event in hwDataNoZ for jet in event], bins=50, range=(0,200), histtype='step', linewidth=1.5, label='Firmware', color='#000000')[0]
nEtaHw = axs[0,1].hist([jet[1] for event in hwDataNoZ for jet in event], bins=18, range=(0,1.5), histtype='step', linewidth=1.5, label='Firmware', color='#000000')[0]
nPhiHw = axs[0,2].hist([jet[2] for event in hwDataNoZ for jet in event], bins=8,  range=(0,0.7), histtype='step', linewidth=1.5, label='Firmware', color='#000000')[0]

nPtEm,  bPtEm  = np.histogram([jet[0] for event in emDataNoZ for jet in event], bins=50, range=(0,200))
nEtaEm, bEtaEm = np.histogram([jet[1] for event in emDataNoZ for jet in event], bins=18, range=(0,1.5))
nPhiEm, bPhiEm = np.histogram([jet[2] for event in emDataNoZ for jet in event], bins=8,  range=(0,0.7))

meansPt  = [0.5*(bPtEm[i]  + bPtEm[i+1])  for i in range(len(nPtEm))]
meansEta = [0.5*(bEtaEm[i] + bEtaEm[i+1]) for i in range(len(nEtaEm))]
meansPhi = [0.5*(bPhiEm[i] + bPhiEm[i+1]) for i in range(len(nPhiEm))]

axs[0,0].scatter(meansPt,  nPtEm,  label='Emulator', c='#380282', linewidths=0.5, s=15)
axs[0,1].scatter(meansEta, nEtaEm, label='Emulator', c='#380282', linewidths=0.5, s=15)
axs[0,2].scatter(meansPhi, nPhiEm, label='Emulator', c='#380282', linewidths=0.5, s=15)

axs[1,0].scatter(meansPt,  [(hw/em) for hw,em in zip(nPtHw,nPtEm)] , c='#380282', linewidths=0.5, s=15)
axs[1,1].scatter(meansEta, [(hw/em) for hw,em in zip(nEtaHw,nEtaEm)], c='#380282', linewidths=0.5, s=15)
axs[1,2].scatter(meansPhi, [(hw/em) for hw,em in zip(nPhiHw,nPhiEm)], c='#380282', linewidths=0.5, s=15)

axs[1,0].axhline(y=0.993, linewidth=1.5, linestyle='--', c='#000000')
axs[1,1].axhline(y=0.993, linewidth=1.5, linestyle='--', c='#000000')
axs[1,2].axhline(y=0.993, linewidth=1.5, linestyle='--', c='#000000')

axs[1,0].set(ylim=(0.5,1.5))
axs[1,1].set(ylim=(0.5,1.5))
axs[1,2].set(ylim=(0.5,1.5))

axs[0,0].set(ylabel="Events")
axs[1,0].set(ylabel="FW / EMU")

axs[0,0].legend(prop={'size': 10})
axs[0,1].legend(prop={'size': 10})
axs[0,2].legend(prop={'size': 10})

ymaxPt  = max(np.concatenate([nPtHw,nPtEm]))
ymaxEta = max(np.concatenate([nEtaHw,nEtaEm]))
ymaxPhi = max(np.concatenate([nPhiHw,nPhiEm]))

axs[0,0].set(xlim=(0,200))
axs[0,1].set(xlim=(0,1.5))
axs[0,2].set(xlim=(0,0.7))

axs[0,0].set(ylim=(0,ymaxPt +(0.05*ymaxPt)))
axs[0,1].set(ylim=(0,ymaxEta+(0.05*ymaxEta)))
axs[0,2].set(ylim=(0,ymaxPhi+(0.05*ymaxPhi)))


#axs[0,0].set_title("Histogrammed PF Jet FW vs EMU: pT, ttbar, 3900 events")  
#axs[0,1].set_title("Histogrammed PF Jet FW vs EMU: Eta, ttbar, 3900 events")  
#axs[0,2].set_title("Histogrammed PF Jet FW vs EMU: Phi, ttbar, 3900 events") 

axs[0,0].set(xlabel="Jet $p_T$ (GeV)")
axs[0,1].set(xlabel="Jet $\eta$")
axs[0,2].set(xlabel="Jet $\phi$")

plt.savefig('ttbarPU200_3900.pdf', bbox_inches='tight')
plt.show()
