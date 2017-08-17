import matplotlib.pyplot as plt
from ROOT_IO import readROOT

t = []
c1 = []
c2 = []

for i in range(0,1):
    filename=u'/home/amcvea/Documents/ketek_data/20170721-b2210-30.0V-73.5V-coin-ketex-%d.root' %(i)
    t_1,c1_1,c2_1,c3_1,c4_1 = readROOT(filename,'2.5G')
    for k in range(0,len(t_1)):
    	if c1_1[k][0]<100 and c2_1[k][0]<100:
    	    t.extend(t_1)
    	    c1.extend(c1_1)
    	    c2.extend(c2_1)

for j in range(0,len(t)):
    plt.figure(1)
    plt.plot(t[j],c1[j])
    #plt.plot(t[j],c1[j],'ro')
    plt.figure(2)
    plt.plot(t[j],c2[j])
    #plt.plot(t[2],c2[2],'ro')

plt.show()

