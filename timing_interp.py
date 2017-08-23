import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from ROOT_IO import readROOT
import scipy.interpolate

#Leading edge discrimination
def led(t,y1,y2,edge,left1,right1,left2,right2):
    time1 = []
    time2 = []
    filt1 = []
    filt2 = []
    filtime = []
    array_y1 = []
    array_y2 = []
    array_t1 = []
    array_t2 = []
    time_diff = []

    #Select only events in the 511keV peak and create filtered arrays
    for i in range(0,len(y1)):
	if right1>np.max(y1[i])>left1 and right2>np.max(y2[i])>left2:
	    for j in range(0,1024):
		if y1[i][j]==np.max(y1[i]):
		    filt1.append(array_y1)
		    filt2.append(array_y2)
		    filtime.append(array_t1)
		    array_y1 = []
		    array_t1 = []
		    array_y2 = []
		    break
		else:
		     array_y1.append(y1[i][j])
		     array_y2.append(y2[i][j])
		     array_t1.append(t[i][j])


    #Find time when events break threshold
    for j in range(0,len(filt1)):
        time1.append(np.interp(edge,filt1[j],filtime[j]))
        time2.append(np.interp(edge,filt2[j],filtime[j]))

    time_diff = np.subtract(time1,time2)

    return time_diff

def fitfunc(p,x):
    return p[0]/(np.sqrt(2*np.pi))*np.exp(-((x-p[1])**2)/(2*p[2]**2))
def residual(p,x,y,dy):
    return (fitfunc(p,x)-y)/dy

t = []
c1 = []
c2 = []

for i in range(0,7):
    filename=u'/home/amcvea/Documents/ketek_data/LYSO_07-21-2017/20170721-b2210-30.0V-73.5V-coin-ketex-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_08-20-2017/2017818-bgo-ketek-31.0V-56.8V-%d.root' %(i)
    t_1,c1_1,c2_1,c3_1,c4_1 = readROOT(filename,'2.5G')
    t.extend(t_1)
    c1.extend(c1_1)
    c2.extend(c2_1)

amp1 = []
amp2 = []
for i in range(0,len(c1)):
        amp1.append(np.max(c1[i]))
        amp2.append(np.max(c2[i]))

t_diff = []
y2 = []

for i in range(10,30):
    t_diff.extend(led(t,c1,c2,i,260,340,330,400))

    y2[i] = []

    for l in range(0,len(t_diffl)):
        if np.abs(t_diffl[l])<10:
            y2[i].append(t_diffl[l])

    largest2 = np.ceil(np.max(y2[i]))
    smallest2 = np.floor(np.min(y2[i]))

bins2 = 100*np.int(np.subtract(largest2,smallest2))

counts2, time2 = np.histogram(y2,bins2)
center2 = (time2[1:] + time2[:-1])/2

p02 = [10, 0, 1]

pf2, cov2, info2, mesg2, success2 = optimize.leastsq(residual, p02, args=(center2, counts2, 0.001), full_output=1)

if success2 >= 4:
    print "Not converged"
    print mesg2
else:
    print "LED Converged"
    print "Number of Events = ", len(y2)
    print "Number of Bins = ", bins2
    print "Mean = ", pf2[1], "ns"
    print "Standard Deviation = ", pf2[2], "ns"

    fwhm2 = pf2[2]*2*np.sqrt(2*np.log(2))
    print "FWHM = ", fwhm2, "ns"

    channel2 = np.linspace(min(center2),max(center2),1000)

    plt.figure(1)
    plt.hist(y2,bins2,histtype='step')
    plt.plot(center2, fitfunc(pf2,center2),'r-')
    plt.xlabel("Time Difference (ns)")
    plt.ylabel("Counts")
    plt.title("Time Resolution (LED)")

    plt.figure(2)
    #plt.subplot(211)
    plt.hist(amp1,300,histtype='step')

    plt.figure(3)
    #plt.subplot(212)
    plt.hist(amp2,300,histtype='step')

    plt.show()
