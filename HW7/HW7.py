# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

data = np.load('strain.npy')
days = np.linspace(0,len(data)/60/24,len(data))
seconds = np.linspace(0,len(data)*60,len(data))

plt.figure()
plt.plot(days,data)
plt.title('Raw Data')
plt.xlabel('Days')
plt.ylabel('Strain h')

#Sampling rate and Nyquist frequency
Fs = abs(1/(seconds[1]-seconds[0]))
nyq = Fs/2
freq = np.arange(0,nyq+Fs/len(data),Fs/len(data))

f = np.fft.rfft(data - np.average(data))
f = [a*2/len(data) for a in f]
f[0] = f[0]/2
f[-1] = f[-1]/2
p = [x*np.conj(x) for x in f]

plt.figure()
plt.plot(np.log10(freq),np.log10(p))
plt.xlabel('$Log_{10}(ƒ)$ ')
plt.ylabel('$Log_{10}$ of Power Spectrum')


pts = []
for i in np.log10(p):
    if i > -42.5:
        pts.append(i)
index = np.where(np.log10(p) == pts[-1])[0]

#Log10(f_GW), Log10(Power Spectrum), Day, Strain
print(np.log10(freq)[index],np.log10(p)[index],days[index],data[index])

#2 equations, 2 unknowns

