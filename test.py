import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from ROOT_IO import readROOT
import scipy.interpolate
from scipy.signal import medfilt #***

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

#finds baseline of waveform ***all new
def baseline(y):
    baseline1 = []
    for i in range(0,len(y)):
        base = 0
        count = 0
        for j in range(4,len(y[i])):
            if y[i][j]>-2 and count>0:
                baseline1.append(base/count)
                break
            else:
                base = base + y[i][j]
                count = count + 1
    return baseline1

#Leading edge discrimination
def led(t,y1,y2,edge1,edge2,left1,right1,left2,right2):
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

    for a in range(0,len(y1)): #***for loop
        for b in range(0,4):
            if y1[a][b]>0 or y1[a][b]<0:
                y1[a][b]=-2
            if y2[a][b]>0 or y2[a][b]<0:
                y2[a][b]=-2

    #Select only events in the 511keV peak and create filtered arrays
    for i in range(0,len(y1)):
	if y1[i][4]>20 or y2[i][4]>20: #*** if statement
            continue
        if right1>np.max(y1[i])>left1 and right2>np.max(y2[i])>left2:
            for j in range(0,1024):
		if y1[i][j]>edge1 and y2[i][j]>edge2: #***changed filling of arrays
                    	   array_y1.append(y1[i][j])
                    	   array_y2.append(y2[i][j])
                           array_t1.append(t[i][j])
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
        time1.append(np.interp(edge1,filt1[j],filtime[j]))
        time2.append(np.interp(edge2,filt2[j],filtime[j]))

    time_diff = np.subtract(time1,time2)
    return time_diff

def fitfunc(p,x):
    return p[0]/(np.sqrt(2*np.pi))*np.exp(-((x-p[1])**2)/(2*p[2]**2))+p[3]
def residual(p,x,y,dy):
    return (fitfunc(p,x)-y)/dy

t = []
c1 = []
c2 = []

for i in range(0,2):
    #filename=u'/home/amcvea/Documents/ketek_data/20170830-bgo-ketek-hamamatsu-31.0V-56.0V-20uCi-150ns-delay-study-0.root'
    #filename=u'/home/amcvea/Documents/ketek_data/BGOxLYSO_09_08_2017/20170908-bgo-ketek-lyso-hamamatsu-32.5V-54.5V-20uCi-150ns-delay-20mV-trigger-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGOxLYSO_09_12_2017/20170912-lyso-hamamatsu-bgo-ketek-34.0V-55.5V-20uCi-150ns-delay-20mV-trigger-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_09_02_2017/20170901-lyso-ketek-lyso-hamamatsu-30.0V-55.5V-20uCi-150ns-delay-20mV-trigger-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Side-Setup_09_12_2017/20170912-bgo-hamamatsu-bgo-ketek-34.5V-57.5V-20uCi-150ns-delay-20mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Side-Setup_09_14_2017/20170914-bgo-ketek-bgo-hamamatsu-33.5V-57.5V-20uCi-150ns-delay-20mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_09_12_2017/20170911-lyso-hamamatsu-lyso-ketek-30.0V-54.5V-20uCi-150ns-delay-20mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_Side-Setup_09_18_2017/20170918-lyso-ketek-lyso-hamamatsu-30.0V-55.5V-20uCi-150ns-delay-20mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_09-05-2017/20170905-bgo-ketek-hamamatsu-32.5V-57.5V-20uCi-150ns-delay-20mV-trigger-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGOxLYSO_Side-Setup_09_19_2017/20170919-bgo-ketek-lyso-hamamatsu-32.0V-55.5V-20uCi-150ns-delay-20mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_Side-Setup_09_19_2017/20170919-lyso-ketek-lyso-hamamatsu-30.0V-55.5V-20uCi-150ns-delay-200mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGOxLYSO_09_25_2017/20170925-bgo-ketek-lyso-hamamatsu-33.5V-55.5V-20uCi-150ns-delay-80mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Amplifier_09_26_2017/20179027-bgo-ketek-bgo-hamamatsu-29.0mV-57.5mV-20uCi-60mV-trigger-side-setup-amplifier-0.root'
    #filename=u'/home/amcvea/Documents/ketek_data/BGOxLYSO_10_03_2017/20171003-bgo-ketek-lyso-hamamatsu-29.0mV-55.5mV-20uCi-120mV-trigger-side-setup-amplifier-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Inverted_10_06_2017/20171006-bgo-hamamatsu-bgo-ketek-58.0mV-29.0mV-20uCi-35mV-trigger-inverted-amplifier-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Amplifier_10_09_2017/20171009-bgo-ketek-bgo-hamamatsu-29.0mV-58.0mV-20uCi-35mV-trigger-amplifier-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Amplifier_10_12_2017/20171012-bgo-ketek-bgo-hamamatsu-29.0mV-58.5mV-20uCi-30mV-trigger-amplifier-%d.root' %(i)
    filename=u'/home/amcvea/Documents/ketek_data/BGO_Amplifier_10_17_2017/20171016-bgo-ketek-bgo-hamamatsu-28.5mV-58.0mV-20uCi-30mV-trigger-amplifier-%d.root' %(i)
    t_1,c1_1,c2_1,c3_1,c4_1 = readROOT(filename,'2.5G')
    t.extend(t_1)
    c1.extend(c1_1)
    c2.extend(c2_1)

c1 = scipy.signal.medfilt(c1,kernel_size=9) #***
c2 = scipy.signal.medfilt(c2,kernel_size=9)

base1 = baseline(c1) #***
base2 = baseline(c2)

amp1 = []
amp2 = []
for i in range(0,len(c1)):
        c1[i] = c1[i] - base1[i] #***
        c2[i] = c2[i] - base2[i]
        amp1.append(np.max(c1[i]))
        amp2.append(np.max(c2[i]))

start = 2
end = 9
res = []
counter1 = 0
for i in range(start,end):
    y2 = []
    t_diff = []
    #t_diff.extend(led(t,c1,c2,i,i,60,200,40,210))
    #t_diff.extend(led(t,c1,c2,i,i,40,90,30,60))
    t_diff.extend(led(t,c1,c2,i,i,70,100,40,70))
    #t_diff.extend(led(t,c1,c2,3+(i*0.1),3+(i*0.1),100,140,50,60))
    #print 3+(i*0.1)

    for l in range(1000,len(t_diff)):
        if np.abs(t_diff[l])<20 and counter1<20:
            y2.append(t_diff[l])
            """
	    plt.figure(l)
            plt.plot(t[l],c1[l])
            counter1 = counter1 + 1
            plt.plot(t[l],c2[l])
            ted1 = edge_break(t[l],c1[l],i)
            ted2 = edge_break(t[l],c2[l],i)
            plt.plot(ted1,i,'ro')
            plt.plot(ted2,i,'go')
	    """

    largest1 = np.ceil(np.max(y2))
    smallest1 = np.floor(np.min(y2))
    bins1 = 10*np.int(np.subtract(largest1,smallest1))
    #bins1 = 100
    counts2, time2 = np.histogram(y2,bins1)
    center2 = (time2[1:] + time2[:-1])/2

    p02 = [10, 0, 1, 0]

    pf2, cov2, info2, mesg2, success2 = optimize.leastsq(residual, p02, args=(center2, counts2, 0.001), full_output=1)

    if success2 >= 4:
        print "Not converged"
        print mesg2
    else:
        print "LED Converged: %dns" %(i) #***
        print "Number of Events = ", len(y2)
        print "Number of Bins = ", bins1
        print "Mean = ", pf2[1], "ns"
        print "Standard Deviation = ", pf2[2], "ns"

        fwhm2 = pf2[2]*2*np.sqrt(2*np.log(2))
        print "FWHM = ", fwhm2, "ns"
        res.append(np.abs(fwhm2))
	print "\n"

        channel2 = np.linspace(min(center2),max(center2),1000)

        plt.figure(i+2)
        plt.hist(y2,bins1,histtype='step')
        plt.plot(center2, fitfunc(pf2,center2),'r-')
        plt.xlabel("Time Difference (ns)")
        plt.ylabel("Counts")
        plt.title("Time Resolution (LED)")
        """
        plt.figure(2)
        plt.hist(amp1,300,histtype='step')
        plt.title("Energy Resolution (Ketek)")
        plt.xlabel("Energy (mV)")
        plt.ylabel("Counts")

        plt.figure(3)
        plt.hist(amp2,300,histtype='step')
        plt.title("Energy Resolution (Hamamatsu)")
        plt.xlabel("Energy (mV)")
        plt.ylabel("Counts")
        """

discrim = np.linspace(start,end-1,num=len(res))
plt.figure(1)
plt.plot(discrim,res,'b-')
plt.ylabel("Discrimination (mV)")
plt.ylabel("Resolution (ns)")
plt.title("Time Resolution (LED)")

plt.show()
