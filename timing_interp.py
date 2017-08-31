import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from ROOT_IO import readROOT
import scipy.interpolate

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
    return p[0]/(np.sqrt(2*np.pi))*np.exp(-((x-p[1])**2)/(2*p[2]**2))+p[3]
def residual(p,x,y,dy):
    return (fitfunc(p,x)-y)/dy

t = []
c1 = []
c2 = []

for i in range(0,2):
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_07-21-2017/20170721-b2210-30.0V-73.5V-coin-ketex-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_08-20-2017/2017818-bgo-ketek-31.0V-56.8V-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_08-22-2017/20170822-bgo-ketek-31.0V-56.8V-20uCi-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_08-24-2017/20170824-lyso-ketek-hamamatsu-31.0V-56.8V-20uCi-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_08-28-2017/20170828-lyso-ketek-hamamatsu-31.0V-56.8V-20uCi-30mV-trig-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_08-29-2017/20170829-lyso-ketek-hamamatsu-29.0V-54.5V-20uCi-30mV-trig-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_08-30-2017/20170830-lyso-ketek-hamamatsu-29.0V-54.5V-20uCi-100ns-delay-study-0.root'
    filename=u'/home/amcvea/Documents/ketek_data/20170830-bgo-ketek-hamamatsu-31.0V-56.0V-20uCi-150ns-delay-study-0.root'
    t_1,c1_1,c2_1,c3_1,c4_1 = readROOT(filename,'2.5G')
    t.extend(t_1)
    c1.extend(c1_1)
    c2.extend(c2_1)

amp1 = []
amp2 = []
for i in range(0,len(c1)):
        amp1.append(np.max(c1[i]))
        amp2.append(np.max(c2[i]))

res = []
counter1 = 0
for i in range(10,20):
    y2 = []
    t_diff = []
    t_diff.extend(led(t,c1,c2,i,0,500,0,500))
    #t_diff.extend(led(t,c1,c2,i,180,240,230,280))

    for l in range(0,len(t_diff)):
        if np.abs(t_diff[l])<10 and counter1<20:
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
    #bins1 = 50*np.int(np.subtract(largest1,smallest1))
    bins1 = 100
    counts2, time2 = np.histogram(y2,bins1)
    center2 = (time2[1:] + time2[:-1])/2
    """
    p02 = [10, 0, 1, 0]

    pf2, cov2, info2, mesg2, success2 = optimize.leastsq(residual, p02, args=(center2, counts2, 0.001), full_output=1)

    if success2 >= 4:
        print "Not converged"
        print mesg2
    else:
        print "LED Converged"
        print "Number of Events = ", len(y2)
        print "Number of Bins = ", bins1
        print "Mean = ", pf2[1], "ns"
        print "Standard Deviation = ", pf2[2], "ns"
        fwhm2 = pf2[2]*2*np.sqrt(2*np.log(2))
        print "FWHM = ", fwhm2, "ns"
        res.append(np.abs(fwhm2))

        channel2 = np.linspace(min(center2),max(center2),1000)

        plt.figure(i+4)
        plt.hist(y2,bins1,histtype='step')
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
     """
discrim = np.linspace(10,20,num=len(res))
plt.figure(1)
plt.plot(discrim,res,'b-')
plt.ylabel("Discrimination (mV)")
plt.ylabel("Resolution (ns)")
plt.title("Time Resolution (LED)")

plt.hist(y2,bins1,histtype='step')
plt.plot(center2, fitfunc(pf2,center2),'r-')
plt.xlabel("Time Difference (ns)")
plt.ylabel("Counts")
plt.title("Time Resolution (LED)")


print counter1
plt.show()
