from matplotlib.pyplot import*
from matplotlib.patches import Rectangle
style.use("dark_background")
from numpy import *

c         = 3e5     # speed of light km/s
s11       = 0       # s11
QL        = 70_000  # Q factor
N         = 10_000  # Integration time
T         = 5.6     # (k)
kb        = 1.38e-1 # 1e-22
bandwidth = 125e-6  # MHz
scan      = 0.5     # MHz
fl        = 749     # MHz
fr        = 751     # MHz
shift     = 2e-2    # MHz

sigma     = kb*T*bandwidth*1e6/sqrt(N) # sigma

def sigle_power(x,f0):
    return 3*(f0/750)*(1-2*s11)/(1-s11) * (1/(1+4*QL*QL*pow(x/f0-1,2)))

def single_bin():
    x = linspace(fl,fr,10000)  # x axis
    grand_array = x*0          # [0,0,...,0] 
    deltax = x[1] - x[0]       # dx
    bin_index = random.randint(5000,7500) # randomly choose  bin position
    bin_pos = x[bin_index]                # bin x axis

    noise = np.random.normal(0,sigma,10000) # noise
    bin_pow = max(noise)/3#abs(noise[bin_index]) # signal power
    step = int((fr-fl-scan)/shift)  # How many step we need to scan full range
    each_len = int(scan/deltax) # How many data point in scan
    each_change = int(shift/deltax) # How many data point in a shift

    #fig, ax = subplots()
    for i in range(step):
        lo = fl + shift * i  # left boundary 
        hi = lo + scan       # right boundary 
        resance =  lo + scan/2  # resonance frequency
        noise = np.random.normal(0,sigma,10000) # noise
        noise[bin_index] += bin_pow # add the bin in noise 

        scan_x = x[i*each_change:each_len+i*each_change] # x axis in this scan
        scan_result = noise[i*each_change:each_len+i*each_change] # power get in this scan
        w = sum(sigle_power(scan_x,resance)/sigma) # weight

         # add in grand array
        grand_array[i*each_change:each_len+i*each_change] = grand_array[i*each_change:each_len+i*each_change] + scan_result*w

        ###### plot part
        rect = Rectangle((lo,0),scan,bin_pow,fill=None,edgecolor="r",linestyle="--")
        ax = subplot(121)
        cla()
        title(f"signal : {bin_pos:6.3f}(MHz)")
        plot([bin_pos,bin_pos],[0,bin_pow],label="signal")
        plot([resance,resance],[0,bin_pow],'r',label="resonace frequency")
        ylabel("$power [10^{-22}]$")
        xlabel("Frequency[MHz]")
        ax.add_patch(rect)
        legend()

        plot(x,x*0,"black")
        ax2 = subplot(122)
        cla()
        title(f"{lo} ~ {hi}")
        plot(scan_x,scan_result)
        if (i*each_change<=bin_index<each_len+i*each_change):
            ax2.arrow(x[bin_index],7,0,-1,color="r",head_width=0.04, head_length=0.4)
        ylim(-5,10)
        savefig(f"subgle_bin/{i}.png")
        pause(0.001)
                # subplot(212)
                # plot([bin_pos,bin_pos],[0,bin_pow])
    figure()
    max_pow = argmax(grand_array)
    title(f"Max : {x[max_pow]:8.5f} signal : {bin_pos:8.5f}(MHz)")
    plot(x,grand_array)
    xlabel("Frequency[MHz]")
    savefig(f"result.png")
    show()

if __name__=="__main__":
    single_bin()