import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from ROOT_IO import readROOT

#Constant fraction discrimination
def cfd(t,y,frac):
    amp = []
    discrim = []
    time = []
    for i in range(0,1000):
        amp.append(min(y[i]))
        discrim.append(np.multiply(amp[i],frac))

        for j in range(0,1024):
            if y[i][j]<discrim[i]:
                time.append(t[i][j])
                break
    return time

#Add later Leading edge discrimination
def led(t,y,edge):
    time = []

    for i in range(0,1000):
        for j in range(0,1024):
            if y[i][j]<edge:
                time.append(t[i][j])
                break
    return time

def fitfunc(p,x):
    return p[0]/(np.sqrt(2*np.pi))*np.exp(-((x-p[1])**2)/(2*p[2]**2))
def residual(p,x,y,dy):
    return (fitfunc(p,x)-y)/dy

filename=u'/home/amcvea/Documents/ketek_data/20170721-b2210-30.0V-73.5V-coin-ketex-0.root'
t,c1,c2,c3,c4 = readROOT(filename,'2.5G')

#t1 = cfd(t,c1,0.1)
#t2 = cfd(t,c2,0.1)

t1 = led(t,c1,-5)
t2 = led(t,c2,-5)

t_diff = np.subtract(t1,t2)

y = []
for k in range(0,1000):
    if np.abs(t_diff[k])<3:
        y.append(t_diff[k])

largest = np.ceil(np.max(y))
smallest = np.floor(np.min(y))

bins = 10*np.int(np.subtract(largest,smallest))

counts, time = np.histogram(y,bins)
center = (time[1:] + time[:-1])/2

p0 = [10, 0, 1]

pf, cov, info, mesg, success = optimize.leastsq(residual, p0, args=(center, counts, 0.1), full_output=1)

if success >= 4:
    print "Not converged"
    print mesg
else:
    print "Converged"
    print "Mean = ", pf[1], "ns"
    print "Standard Deviation = ", pf[2], "ns"

    fwhm = pf[2]*2*np.sqrt(2*np.log(2))
    print "FWHM = ", fwhm, "ns"

    channel = np.linspace(min(center),max(center),1000)
    plt.hist(y,bins,histtype='step')
    plt.plot(center, fitfunc(pf,center),'r-')
    plt.xlabel("Time Difference (ns)")
    plt.ylabel("Counts")
    plt.title("Time Resolution")

    plt.show()
