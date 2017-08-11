import matplotlib.pyplot as plt
import numpy as np
from ROOT_IO import readROOT

#Constant fraction discrimination
def cfd(t,y,frac):
    amp = []
    discrim = []
    time = []
    for i in range(0,500):
        amp.append(min(y[i]))
        discrim.append(np.multiply(amp[i],frac))

        for j in range(0,500):
            if y[i][j]<discrim[i]:
                time.append(t[i][j])
                break
    return time

#Add later Leading edge discrimination
#def led(x,y):


filename=u'/home/amcvea/Documents/ketek_data/20170721-b2210-30.0V-73.5V-coin-ketex-0.root'
t,c1,c2,c3,c4 = readROOT(filename,'2.5G')

t1 = cfd(t,c1,0.1)
t2 = cfd(t,c2,0.1)

t_diff = np.subtract(t1,t2)

print len(t_diff)

largest = np.ceil(np.max(t_diff))
smallest = np.floor(np.min(t_diff))

bins = 10*np.int(np.subtract(largest,smallest))

plt.hist(t_diff,bins)
plt.axis([-5,5,0,25])

plt.show()


