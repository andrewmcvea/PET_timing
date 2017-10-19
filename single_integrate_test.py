import numpy as np
import matplotlib.pyplot as plt
from ROOT_IO import readROOT
from scipy.signal import medfilt

#finds baseline of waveform
def baseline(y):
    baseline1 = []
    for i in range(0,len(y)):
        base = 0
        count = 0
        for j in range(4,len(y[i])):
            if y[i][j]>0 and count>0:
                baseline1.append(base/count)
                break
            else:
                base = base + y[i][j]
                count = count + 1
    return baseline1

t = []
c1 = []
c2 = []

for i in range(0,3):
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Side-Setup_09_12_2017/20170912-bgo-hamamatsu-bgo-ketek-34.5V-57.5V-20uCi-150ns-delay-20mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_09_02_2017/20170901-lyso-ketek-lyso-hamamatsu-30.0V-55.5V-20uCi-150ns-delay-20mV-trigger-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Side-Setup_09_14_2017/20170914-bgo-ketek-bgo-hamamatsu-33.5V-57.5V-20uCi-150ns-delay-20mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_09_12_2017/20170911-lyso-hamamatsu-lyso-ketek-30.0V-54.5V-20uCi-150ns-delay-20mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_Side-Setup_09_18_2017/20170918-lyso-ketek-lyso-hamamatsu-30.0V-55.5V-20uCi-150ns-delay-20mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_09-05-2017/20170905-bgo-ketek-hamamatsu-32.5V-57.5V-20uCi-150ns-delay-20mV-trigger-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGOxLYSO_Side-Setup_09_19_2017/20170919-bgo-ketek-lyso-hamamatsu-32.0V-55.5V-20uCi-150ns-delay-20mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGOxLYSO_09_12_2017/20170912-lyso-hamamatsu-bgo-ketek-34.0V-55.5V-20uCi-150ns-delay-20mV-trigger-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/LYSO_Side-Setup_09_19_2017/20170919-lyso-ketek-lyso-hamamatsu-30.0V-55.5V-20uCi-150ns-delay-200mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Energy_Res/20170925-bgo-ketek-33.5V-20uCi-150ns-delay-5mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Energy_Res/20170925-bgo-ketek-35.0V-20uCi-150ns-delay-5mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Energy_Res/20170926-bgo-hamamatsu-57.5V-20uCi-150ns-delay-5mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Energy_Res/20170926-bgo-hamamatsu-57.5V-20uCi-150ns-delay-60mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGOxLYSO_09_25_2017/20170925-bgo-ketek-lyso-hamamatsu-33.5V-55.5V-20uCi-150ns-delay-80mV-trigger-side-setup-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Amplifier_09_26_2017/20170926-bgo-ketek--29.0V-20uCi-150ns-delay-5mV-trigger-side-setup-amplifier-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Amplifier_09_26_2017/20179027-bgo-ketek-bgo-hamamatsu-29.0mV-57.5mV-20uCi-60mV-trigger-side-setup-amplifier-0.root'
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Energy_Res/20171002-15mm-bgo-ketek-29.0mV-20uCi-5mV-trigger-side-setup-amplifier-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGOxLYSO_10_03_2017/20171003-bgo-ketek-lyso-hamamatsu-29.0mV-55.5mV-20uCi-120mV-trigger-side-setup-amplifier-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Hamamatsu_Energy_Res/20171005-bgo-hamamatsu-29.0mV-57.5mV-20uCi-5mV-trigger-side-setup-amplifier-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Hamamatsu_Energy_Res/20171005-bgo-hamamatsu-29.0mV-57.75mV-20uCi-5mV-trigger-side-setup-amplifier-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Hamamatsu_Energy_Res/20171005-bgo-hamamatsu-29.0mV-58.0mV-20uCi-5mV-trigger-side-setup-amplifier-%d.root' %(i)
    #filename=u'/home/amcvea/Documents/ketek_data/BGO_Hamamatsu_Energy_Res/20171005-bgo-hamamatsu-57.5mV-20uCi-5mV-trigger-inverted-amplifier-%d.root' %(i)
    filename=u'/home/amcvea/Documents/ketek_data/BGO_Hamamatsu_Energy_Res/20171009-bgo-hamamatsu-58.5mV-20uCi-5mV-trigger-inverted-amplifier-%d.root' %(i)
    t_1,c1_1,c2_1,c3_1,c4_1 = readROOT(filename,'2.5G')
    t.extend(t_1)
    #c1.extend(c1_1)
    c2.extend(c2_1)

#c1 = medfilt(c1,kernel_size=9)
c2 = medfilt(c2,kernel_size=9)

#base1 = baseline(c1)
base2 = baseline(c2)

#amp1 = []
amp2 = []
#area1 = []
area2 = []
counter = 0

for k in range(0,len(c2)):
    #y1 = []
    y2 = []

    #c1[k] = c1[k] - base1[k]
    c2[k] = c2[k] - base2[k]
    #amp1.append(np.max(c1[k]))
    amp2.append(np.max(c2[k]))
    for i in range(0,len(c2[k])):
    	#if c1[k][i]<0:
	#    y1.append(0)
	#else:
	#    y1.append(c1[k][i])

	if c2[k][i]<0:
	    y2.append(0)
    	else:
	    y2.append(c2[k][i])

    #int1 = np.trapz(y1,x=t[k])
    int2 = np.trapz(y2,x=t[k])
    """
    if int2<1 and counter<15:
	normalizer = np.sum(int2)
	int2 = np.divide(int1,normalizer)
	plt.figure(k+4)
	plt.plot(t[k],c1[k])
	plt.plot(t[k],c2[k])
	counter = counter + 1
    """
    #area1.append(int1)
    area2.append(int2)


plt.figure(1)
plt.title("Ketek Amplitude")
#plt.hist(amp1,100,histtype='step')

plt.figure(2)
plt.title("Hamamatsu Amplitude")
plt.hist(amp2,100,histtype='step')

plt.figure(3)
plt.title("Ketek Integral")
#plt.hist(area1,100,histtype='step')

plt.figure(4)
plt.title("Hamamatsu Integral")
plt.hist(area2,100,histtype='step')

plt.show()
