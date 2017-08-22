import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from ROOT_IO import readROOT

#Constant fraction discrimination
def cfd(t,y1,y2,frac,left1,right1,left2,right2):
    discrim1 = []
    discrim2 = []
    time1 = []
    time2 = []
    filt1 = []
    filt2 = []
    filtime = []

    for i in range(0,len(y1)):
	if right1>max(y1[i])>left1 and right2>max(y2[i])>left2:
	    filt1.append(y1[i])
	    filt2.append(y2[i])
            discrim1.append(np.divide(max(y1[i]),frac))
            discrim2.append(np.divide(max(y2[i]),frac))
	    filtime.append(t[i])
 	else:
	    continue

    for i in range(0,len(filt1)):
        for j in range(0,1024):
            #if filt1[i][j]>discrim1[i] and filt1[i][5]<discrim1[i] and filt2[i][5]<discrim2[i]:
            if filt1[i][j]>discrim1[i] and filt1[i][5]<0 and filt2[i][5]<0:
                time1.append(filt1[i][j])
		break
            else:
	        continue

	for k in range(0,1024):
	    #if filt2[i][k]>discrim2[i] and filt1[i][5]<discrim1[i] and filt2[i][5]<discrim2[i]:
            if filt2[i][k]>discrim2[i] and filt1[i][5]<0 and filt2[i][5]<0:
                time2.append(filt2[i][k])
                break
	    else:
		continue

    time_diff = np.subtract(time1,time2)
    return time_diff

#Leading edge discrimination
def led(t,y1,y2,edge,left1,right1,left2,right2):
    time1 = []
    time2 = []
    filt1 = []
    filt2 = []
    filtime = []
    time_diff = []

    #Select only events in the 511keV peak and create filtered arrays
    for i in range(0,len(y1)):
	if right1>max(y1[i])>left1 and right2>max(y2[i])>left2:
	    filt1.append(y1[i])
	    filt2.append(y2[i])
	    filtime.append(t[i])

    #Find time at which events go above threshold
    for i in range(0,len(filt1)):
        for j in range(0,1024):
            if filt1[i][j]>edge and filt1[i][5]<0 and filt2[i][5]<0:
                time1.append(filtime[i][j])
		break
	    else:
		continue

        for k in range(0,1024):
	    if filt2[i][k]>edge and filt1[i][5]<0 and filt2[i][5]<0:
                time2.append(filtime[i][k])
                break
	    else:
		continue

    time_diff = np.subtract(time1,time2)
    print time_diff
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
        amp1.append(max(c1[i]))
        amp2.append(max(c2[i]))

t_diffc = cfd(t,c1,c2,25,260,340,330,400)
t_diffl = led(t,c1,c2,15,260,340,330,400)

y1 = []
y2 = []

for k in range(0,len(t_diffc)):
    if np.abs(t_diffc[k])<3:
        y1.append(t_diffc[k])

for l in range(0,len(t_diffl)):
    if np.abs(t_diffl[l])<3:
        y2.append(t_diffl[l])

largest1 = np.ceil(np.max(y1))
smallest1 = np.floor(np.min(y1))
largest2 = np.ceil(np.max(y2))
smallest2 = np.floor(np.min(y2))

bins1 = 100*np.int(np.subtract(largest1,smallest1))
bins2 = 100*np.int(np.subtract(largest2,smallest2))

counts1, time1 = np.histogram(y1,bins1)
center1 = (time1[1:] + time1[:-1])/2
counts2, time2 = np.histogram(y2,bins2)
center2 = (time2[1:] + time2[:-1])/2

p01 = [10, 0, 1]
p02 = [10, 0, 1]

pf1, cov1, info1, mesg1, success1 = optimize.leastsq(residual, p01, args=(center1, counts1, 0.001), full_output=1)
pf2, cov2, info2, mesg2, success2 = optimize.leastsq(residual, p02, args=(center2, counts2, 0.001), full_output=1)

if success1 >= 4:
    print "Not converged"
    print mesg1
else:
    print "CFD Converged"
    print "Number of Events = ", len(y1)
    print "Number of Bins = ", bins1
    print "Mean = ", pf1[1], "ns"
    print "Standard Deviation = ", pf1[2], "ns"

    fwhm1 = pf1[2]*2*np.sqrt(2*np.log(2))
    print "FWHM = ", fwhm1, "ns"

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

    channel1 = np.linspace(min(center1),max(center1),1000)
    channel2 = np.linspace(min(center2),max(center2),1000)

    plt.figure(1)
    plt.hist(y1,bins1,histtype='step')
    plt.plot(center1, fitfunc(pf1,center1),'r-')
    plt.xlabel("Time Difference (ns)")
    plt.ylabel("Counts")
    plt.title("Time Resolution (CFD)")

    plt.figure(2)
    plt.hist(y2,bins2,histtype='step')
    plt.plot(center2, fitfunc(pf2,center2),'r-')
    plt.xlabel("Time Difference (ns)")
    plt.ylabel("Counts")
    plt.title("Time Resolution (LED)")

    plt.figure(3)
    plt.subplot(211)
    plt.hist(amp1,300,histtype='step')

    plt.subplot(212)
    plt.hist(amp2,300,histtype='step')

    plt.show()
