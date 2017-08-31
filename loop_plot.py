import matplotlib.pyplot as plt
from ROOT_IO import readROOT
import numpy as np
import scipy.interpolate

led = 20

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

t = []
c1 = []
c2 = []
amp1 = []
amp2 = []
for i in range(0,1):
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_07-21-2017/20170721-b2210-30.0V-73.5V-coin-ketex-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_08-22-2017/20170822-bgo-ketek-31.0V-56.8V-20uCi-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_08-24-2017/20170824-lyso-ketek-hamamatsu-31.0V-56.8V-20uCi-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_08-28-2017/20170828-lyso-ketek-hamamatsu-31.0V-56.8V-20uCi-30mV-trig-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_08-29-2017/20170829-bgo-ketek-hamamatsu-31.0V-56.8V-20uCi-30mV-trig-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_08-29-2017/20170829-lyso-ketek-hamamatsu-29.0V-54.5V-20uCi-30mV-trig-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_08-30-2017/20170830-lyso-ketek-hamamatsu-29.0V-54.5V-20uCi-100ns-delay-study-0.root'
    filename=u'/home/amcvea/Documents/ketek_data/20170830-bgo-ketek-hamamatsu-31.0V-56.0V-20uCi-150ns-delay-study-0.root'

    t_1,c1_1,c2_1,c3_1,c4_1 = readROOT(filename,'2.5G')
    t.extend(t_1)
    c1.extend(c1_1)
    c2.extend(c2_1)

for k in range(0,len(c1)):
    amp1.append(np.max(c1[k]))
    amp2.append(np.max(c2[k]))

print len(c1)

plt.figure(1)
plt.subplot(211)
plt.hist(amp1,500,histtype='step')

plt.subplot(212)
plt.hist(amp2,500,histtype='step')

counter1= 0
counter2= 0


for j in range(187,207):
#for j in range(0,len(c1)):
    for k in range(0,4):
        #if c1[j][k]>10 or c2[j][k]>10:
        #    break
        if c1[j][k]>5 or c1[j][k]<-10:
            c1[j][k]=-4
        if c2[j][k]>5 or c2[j][k]<-10:
            c2[j][k]=-4
        if k==3:
            plt.figure(j)
            plt.plot(t[j],c1[j])
            counter1 = counter1 + 1
            plt.plot(t[j],c2[j])
            ted1 = edge_break(t[j],c1[j],led)
            ted2 = edge_break(t[j],c2[j],led)
            plt.plot(ted1,led,'ro')
            plt.plot(ted2,led,'go')
            plt.title("Captured Waveform")
            plt.xlabel("Time (ns)")
            plt.ylabel("Amplitude (mV)")
                        
print counter1

plt.show()
