#!/usr/bni/env python

import os
import numpy as np
import astropy.io.ascii as ascii

inpfiles = ['my_few-lines', 'my_emission', 'my_absorption', 'my_absorption_emission', 'my_narrow-range', 'my_broad-lines', 'my_extragalactic']
inpnoises = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]

for inpnoise in inpnoises:
    
    for inpfile in inpfiles:

        f = open(inpfile+'.dat', 'r')

        lflux = []
        freqs = []
        counter = 0
        for line in f:
            counter = counter+1
            freqs.append(float(line.split()[0]))
            lflux.append(float(line.split()[1]))
        flux = np.array(lflux)
        freq = np.array(freqs)
        noise = np.random.normal(loc=0.0, scale=inpnoise, size=counter)
        
        
        if (inpnoise < 10):
            ascii.write((freqs, flux+noise), output='noise/'+inpfile+'_n0'+str(inpnoise)+'.dat', format='no_header')
        if (inpnoise >= 10):
            ascii.write((freqs, flux+noise), output='noise/'+inpfile+'_n'+str(inpnoise)+'.dat', format='no_header')
