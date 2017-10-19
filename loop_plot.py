import matplotlib.pyplot as plt
from ROOT_IO import readROOT
import numpy as np
import scipy.interpolate
from scipy.signal import medfilt

led = 3

def edge_break(t,y,edge):
    filt = []
    filtime = []

    for j in range(0,1024):
        if y[j]>led:
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
            if y[i][j]>-2 and count>0:
                baseline1.append(base/count)
                break
            else:
                base = base + y[i][j]
                count = count + 1
    return baseline1

t = []
c1 = []
c2 = []
amp1 = []
amp2 = []
for i in range(0,2):
    #filename=u'/home/amcvea/Documents/ketek_data/BGOxLYSO_10_03_2017/20171003-bgo-ketek-lyso-hamamatsu-29.0mV-55.5mV-20uCi-120mV-trigger-side-setup$
    filename=u'/home/amcvea/Documents/ketek_data/BGO_Inverted_10_06_2017/20171006-bgo-hamamatsu-bgo-ketek-58.0mV-29.0mV-20uCi-35mV-trigger-inverted-$
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Amplifier_10_09_2017/20171009-bgo-ketek-bgo-hamamatsu-29.0mV-58.0mV-20uCi-35mV-trigger-amplifi$
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Amplifier_10_12_2017/20171012-bgo-ketek-bgo-hamamatsu-29.0mV-58.5mV-20uCi-30mV-trigger-amplifi$
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Amplifier_10_17_2017/20171016-bgo-ketek-bgo-hamamatsu-28.5mV-58.0mV-20uCi-30mV-trigger-amplifi$
    t_1,c1_1,c2_1,c3_1,c4_1 = readROOT(filename,'2.5G')
    t.extend(t_1)
    c1.extend(c1_1)
    c2.extend(c2_1)

c1 = scipy.signal.medfilt(c1,kernel_size=9)
c2 = scipy.signal.medfilt(c2,kernel_size=9)

base1 = baseline(c1)
base2 = baseline(c2)


for k in range(0,len(c1)):
    c1[k] = c1[k] - base1[k]
    c2[k] = c2[k] - base2[k]
    amp1.append(np.max(c1[k]))
    amp2.append(np.max(c2[k]))

print len(c1)

#plt.figure(1)
#plt.hist(amp1,100,histtype='step')

#plt.figure(2)
#plt.hist(amp2,100,histtype='step')

counter1= 331
counter2= 0

for j in range(0,40):
#for j in range(0,len(c1)):
    for k in range(0,4):
        if c1[j][k]>10 or c2[j][k]>10:
            break
        if c1[j][k]>2 or c1[j][k]<-10:
            c1[j][k]=0
        if c2[j][k]>2 or c2[j][k]<-10:
            c2[j][k]=0
        if k==3 and counter1<340:
            plt.figure(j)
            #plt.subplot(counter1)
            plt.axis([0,400,-10,140])
            plt.plot(t[j],c1[j])
            counter1 = counter1 + 1
            plt.plot(t[j],c2[j])
            ted1 = edge_break(t[j],c1[j],led)
            ted2 = edge_break(t[j],c2[j],led)
            plt.plot(ted1,led,'ro')
            plt.plot(ted2,led,'go')
            plt.title("Captured Waveform")
            plt.xlabel("Time (ns)")
            #plt.ylabel("Amplitude (mV)")

plt.show()
