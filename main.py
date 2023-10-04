'''NMR is an chemical analysis method used to identify and to quantify a given chemical
compound in
a mixture. An NMR spectrum, i.e. the data returned by the NMR measuring
instrument, is a 2-column table. The first column (x-axis) corresponds to a frequency
shift of the proton resonance frequency, expressed in ppm (part per million) with
respect to a fundamental frequency. The second column (y-axis) represents the
intensity of the signal at a given frequency shift. The signature of a molecule, used to
identify it, is a series of peaks more or less high at well-defined frequencies shift. 

'''

'''
    Kinetic of the reaction
Again, no prior knowledge on chemical kinetics is required to understand the problem.
Basically, the esterification of ethanoic acid (or acetic acid) with propanol is considered
as a non-reversible first-order reaction. Thus, the reaction rate is given by:
        v= k * [EA] * [Pro]
where [EA] and [Pro] are the concentrations of ethanoic acid and propanol, respectively.

Thus, the evolution of [EA], [Pro] and the concentration of propyl-ethanoate [PE]
during the reaction is given by the three differential equations:
        d[EA]/dt = -v;
        d[Pro]/dt = -v;
        d[PE]/dt = v;

For the sake of simplicity, the solvent (water) is never considered

'''

'''
Provided Dataset

The file you are provided with is a .dat file which can be imported on Python using
load() function from Pickle module. The Dataset consists of a dictionary which
contains the following information:
• Title of the dataset (string)
• The frequency shift corresponding to the first value in NMR spectrums (float),
the last value in NMR spectrums (float) and the step value in NMR spectrums
(float). These values can be used to calculate the x-axis of a spectrum
• A dictionary giving the NMR spectrum, i.e. y-axis values of:
o A 1 µmol/L solution of ethanoic acid (float array)
o A 1 µmol/L solution of propanol (float array)
o A 1 µmol/L solution of propyl-ethanoate (float array)
• A table which contains NMR spectrum measured during the reaction process.
Each element of the table is a dictionary that contains:
o The time in hour at which the measurement is performed (float)
o The measured NMR spectrum (float array)

'''


#1.Plot the NMR spectrum of ethanoic acid,propanol and propyl-ethanoate

import pickle
import numpy as np
import matplotlib.pyplot as plt

#load the dataset
with open('Homework_L2_2022.dat','rb') as f:
    data = pickle.load(f)
    
print(data)


#extract the data
title = data['title']
first = data['Frequency Min']
last = data['Frequency Max']
step = data['Frequency Step']
nmr = data['Pure Spectrums']
table = data['Measurements']

#calculate the x-axis
x = np.arange(first,last+step,step)

#plot the NMR spectrum
plt.figure()
plt.plot(x,nmr['EA'],'r',label='EA')
plt.plot(x,nmr['Pro'],'b',label='Pro')
plt.plot(x,nmr['PE'],'g',label='PE')
plt.xlabel('ppm')
plt.ylabel('Intensity')
plt.legend()
plt.title(title)
plt.show()


#2.Write a function that calculate the position and the intensity of each peak in a
#NMR spectrum. Apply it to ethanoic acid, propanol and propyl-ethanoate

def peak(x,y):
    '''Calculate the position and the intensity of each peak in a NMR spectrum
    x: x-axis of the spectrum
    y: y-axis of the spectrum
    '''
    #find the peaks
    peaks = np.where(y>0.5*np.max(y))[0]
    #calculate the position and the intensity of each peak
    pos = x[peaks]
    inten = y[peaks]
    return pos,inten

#apply the function to ethanoic acid, propanol and propyl-ethanoate
pos_EA,inten_EA = peak(x,nmr['EA'])
pos_Pro,inten_Pro = peak(x,nmr['Pro'])
pos_PE,inten_PE = peak(x,nmr['PE'])

#plot the peaks
plt.figure()
plt.plot(x,nmr['EA'],'r',label='EA')
plt.plot(x,nmr['Pro'],'b',label='Pro')
plt.plot(x,nmr['PE'],'g',label='PE')
plt.plot(pos_EA,inten_EA,'ro')
plt.plot(pos_Pro,inten_Pro,'bo')
plt.plot(pos_PE,inten_PE,'go')
plt.xlabel('ppm')
plt.ylabel('Intensity')
plt.legend()
plt.title(title)
plt.show()

#3.. From the initial spectrum (at t=0), estimate the initial concentration of ethanoic
#acid, propanol and propyl-ethanoate

#initial concentration of ethanoic acid
EA0 = 1e-6*inten_EA[0]
#initial concentration of propanol
Pro0 = 1e-6*inten_Pro[0]
#initial concentration of propyl-ethanoate
PE0 = 1e-6*inten_PE[0]

#4.Analyze the NMR spectrums to plot the evolution of the concentration of
#ethanoic acid, propanol and propyl-ethanoate during the reaction.

#calculate the concentration of ethanoic acid, propanol and propyl-ethanoate
#during the reaction
EA = np.zeros(len(table))
Pro = np.zeros(len(table))
PE = np.zeros(len(table))
for i in range(len(table)):
    pos,inten = peak(x,table[i]['NMR'])
    EA[i] = 1e-6*inten[0]
    Pro[i] = 1e-6*inten[1]
    PE[i] = 1e-6*inten[2]

#plot the evolution of the concentration of ethanoic acid, propanol and
#propyl-ethanoate during the reaction
plt.figure()
plt.plot(table['Time'],EA,'r',label='EA')
plt.plot(table['Time'],Pro,'b',label='Pro')
plt.plot(table['Time'],PE,'g',label='PE')
plt.xlabel('Time (h)')
plt.ylabel('Concentration (mol/L)')
plt.legend()
plt.title(title)
plt.show()

#5.From these data, estimate the value of the reaction constant k.

#calculate the reaction rate
v = (EA[1:]-EA[:-1])/(table['Time'][1:]-table['Time'][:-1])

#plot the reaction rate
plt.figure()
plt.plot(table['Time'][:-1],v,'k')
plt.xlabel('Time (h)')
plt.ylabel('Reaction rate (mol/L/h)')
plt.title(title)
plt.show()

#estimate the value of the reaction constant k
k = v[0]/(EA0*Pro0)



