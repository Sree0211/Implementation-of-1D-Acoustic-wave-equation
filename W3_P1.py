import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec

xmax = 10000                            # number of grid points in x-direction
nx = 10000                              # physical domain (m)
dx = xmax/nx                            # grid point distance in x-direction
c0 = 334                                # wave speed in medium (m/s)
isrc = int(nx/2)                        # source location in grid in x-direction
#ir = isrc + 100
nt = 1001                               # maximum number of time steps
dt = 0.001                               # time step


f0 = 25
t0 = 4/f0

idisp = 5

#Plot for Source Time function

#Source time function

src = np.zeros(nt+1)
time = np.linspace(0,nt*dt,nt)

src = -8 * (time-t0)*f0*(np.exp(-1*(4*f0)**2*(time-t0)**2))

plt.ion()
fig1 = plt.figure(figsize=(10,6))
gs1 = gridspec.GridSpec(1,2,width_ratios=[1,1],hspace=0.3,wspace=0.3)

ax1 = plt.subplot(gs1[0])
ax1.plot(time,src)
ax1.set_title('Source Time Function')
ax1.set_xlim(time[0],time[-1])
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Amplitude')

ax2  = plt.subplot(gs1[1])
spec = np.fft.fft(src) # source time function in frequency domain
freq = np.fft.fftfreq(spec.size, d = dt ) # time domain to frequency domain
ax2.plot(np.abs(freq), np.abs(spec)) # plot frequency and amplitude
ax2.set_xlim(0, 250) # only display frequency from 0 to 250 Hz
ax2.set_title('Source Spectrum')
ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Amplitude')

ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")

plt.show()


#Plot Snapshot & Seismogram

p = np.zeros(nx)
pold = np.zeros(nx)
pnew = np.zeros(nx)
d2px = np.zeros(nx)

c = np.zeros(nx)
c = c+c0

x = np.arange(nx)
x = x*dx

plt.ion()
fig2 = plt.figure(figsize=(10,6))
gs2 = gridspec.GridSpec(1,1,width_ratios=[1],hspace=0.3,wspace=0.3)

ax3 = plt.subplot(gs2[0])
leg1,= ax3.plot(isrc,0,'r*',markersize = 11)
up31,= ax3.plot(p)
ax3.set_xlim(0,xmax)
ax3.set_ylim(-np.max(p),np.max(p))
ax3.set_title('Time Step (nt) = 0')
ax3.set_xlabel('x (m)')
ax3.set_ylabel('Pressure Amplitude')

plt.show()

#1-D Wave propagation (FD)

for it in range(nt):
    for i in range(1,nx-1):
        d2px[i] = (p[i+1]-2*p[i]+p[i-1])/dx**2

    pnew = 2*p - pold + c**2*dt**2*d2px

    pnew[isrc] = pnew[isrc] + src[it]/(dx)*dt**2

    pold,p = p,pnew

    if(it%idisp)==0:
        ax3.set_title('Time Step (nt) = %d'%it)
        ax3.set_ylim(-1.1*np.max(abs(p)),1.1*np.max(abs(p)))
        window = 100;xshift = 25
        ax3.set_xlim(isrc*dx+c0*it*dt-window*dx-xshift,isrc*dx+c0*it*dt+window*dx-xshift)
        up31.set_ydata(p)
        plt.gcf().canvas.draw()
            

