import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from ROOT_IO import readROOT
import scipy.signal
from scipy.signal import medfilt

led = 3

def edge_break(t,y,edge):
    filt = []
    filtime = []

    for j in range(0,1024):
        if y[j]>(max(y)-1):
            break
        else:
            filt.append(y[j])
            filtime.append(t[j])

    tled = np.interp(edge,filt,filtime)
    return tled

#finds baseline of waveform
def baseline(y):
    baseline1 = []
    for i in range(0,len(y)):
        base = 0
        count = 0
        for j in range(4,len(y[i])):
            if y[i][j]>-3 and count>0:
                baseline1.append(base/count)
                break
            if j==len(y[i]):
		baseline1.append(base/count)
                break
	    else:
                base = base + y[i][j]
                count = count + 1
    return baseline1

def fitfunc(p,x):
    return p[0]/(np.sqrt(2*np.pi))*np.exp(-((x-p[1])**2)/(2*p[2]**2))+p[3]
def residual(p,x,y,dy):
    return (fitfunc(p,x)-y)/dy

t = []
c1 = []
c2 = []

for i in range(0,3):
    #filename=u'/home/amcvea/Documents/ketek_data/BGOxLYSO_10_03_2017/20171003-bgo-ketek-lyso-hamamatsu-29.0mV-55.5mV-20uCi-120mV-trigger-side-setup-amplifier-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Energy_Res/20170925-bgo-ketek-33.5V-20uCi-150ns-delay-5mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Energy_Res/20171002-15mm-bgo-ketek-29.0mV-20uCi-5mV-trigger-side-setup-amplifier-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Hamamatsu_Energy_Res/20171005-bgo-hamamatsu-29.0mV-57.5mV-20uCi-5mV-trigger-side-setup-amplifier-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Hamamatsu_Energy_Res/20171005-bgo-hamamatsu-29.0mV-57.75mV-20uCi-5mV-trigger-side-setup-amplifier-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Hamamatsu_Energy_Res/20171005-bgo-hamamatsu-29.0mV-58.0mV-20uCi-5mV-trigger-side-setup-amplifier-%d.root' %(i)
    filename=u'/home/amcvea/Documents/ketek_data/BGO_Hamamatsu_Energy_Res/20171009-bgo-hamamatsu-58.5mV-20uCi-5mV-trigger-inverted-amplifier-%d.root' %(i)
    t_1,c1_1,c2_1,c3_1,c4_1 = readROOT(filename,'2.5G')
    t.extend(t_1)
    #c1.extend(c1_1)
    c2.extend(c1_1)

#c1 = scipy.signal.medfilt(c1,kernel_size=9)
c2 = scipy.signal.medfilt(c2,kernel_size=9)

#base1 = baseline(c1)
base2 = baseline(c2)

#amp1 = []
amp2 = []
#y1 = []
y2 = []
timing = []
for i in range(0,len(c2)):
        #c1[i] = c1[i] - base1[i]
        c2[i] = c2[i] - base2[i]
        #amp1.append(np.max(c1[i]))
        amp2.append(np.max(c2[i]))

	for j in range(0,len(c2[i])):
	    #if j<5 and np.abs(c1[i][j])>led:
		#continue
	    if j<5 and np.abs(c2[i][j])>led:
		continue
	    if np.abs(c2[i][j])>led: #or np.abs(c2[i][j])>led:
		break
	    else:
		#y1.append(c1[i][j])
		y2.append(c2[i][j])


bins1 = 50
bins2 = 100

#counts1, time1 = np.histogram(y1,bins1)
counts2, time2 = np.histogram(y2,bins2)

#center1 = (time1[1:] + time1[:-1])/2
center2 = (time2[1:] + time2[:-1])/2
"""
p01 = [10, 0, 1, 0]
pf1, cov1, info1, mesg1, success1 = optimize.leastsq(residual, p01, args=(center1, counts1, 0.001), full_output=1)


if success1 >= 4:
        print "Not converged"
        print mesg1
else:
        print "Converged"
        print "Number of Events = ", len(y1)
        print "Number of Bins = ", bins1
        print "Mean = ", pf1[1], "mV"
        print "Standard Deviation = ", pf1[2], "mV"

        fwhm1 = pf1[2]*2*np.sqrt(2*np.log(2))
        print "FWHM = ", fwhm1, "mV"
        print "\n"
"""

p02 = [10, 0, 1, 0]
pf2, cov2, info2, mesg2, success2 = optimize.leastsq(residual, p02, args=(center2, counts2, 0.001), full_output=1)


if success2 >= 4:
        print "Not converged"
        print mesg2
else:
        print "Converged"
        print "Number of Events = ", len(y2)
        print "Number of Bins = ", bins2
        print "Mean = ", pf2[1], "mV"
        print "Standard Deviation = ", pf2[2], "mV"

        fwhm2 = pf2[2]*2*np.sqrt(2*np.log(2))
        print "FWHM = ", fwhm2, "mV"
        print "\n"

#channel1 = np.linspace(min(center1),max(center1),1000)
channel2 = np.linspace(min(center2),max(center2),1000)
"""
plt.figure(1)
plt.hist(y1,50,histtype='step')
plt.plot(center1, fitfunc(pf1,center1),'r-')
plt.xlabel("Voltage (mV)")
plt.ylabel("Counts")
plt.title("Ketek Noise")
"""
plt.figure(2)
plt.hist(y2,50,histtype='step')
plt.plot(center2, fitfunc(pf2,center2),'r-')
plt.xlabel("Time (ns)")
plt.ylabel("Counts")
plt.title("Hamamatsu Noise")

plt.show()
