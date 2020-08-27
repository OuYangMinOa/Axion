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
def six_power(x,f0):
    return 3*(f0/750)*(1-2*s11)/(1-s11) * (1/(1+4*QL*QL*pow(x/f0-1,2)))

def maxwell(x,x0):
    v = 0.013*5.7
    norm = 4*pi / (v*v*pi)**1.5
    out = []
    for each in x:
        if each<x0:
            out.append(0)
        else:
            now_X = each - x0
            out.append(norm * now_X*now_X*exp(-now_X*now_X/v/v))
    return array(out)

def single_bin():
    x = linspace(fl,fr,10000)
    grand_array = x*0
    deltax = x[1] - x[0]
    bin_index = random.randint(5000,7500)
    bin_pos = x[bin_index]

    noise = np.random.normal(0,sigma,10000)
    bin_pow = max(noise)/3#abs(noise[bin_index])
    step = int((fr-fl-scan)/shift)
    each_len = int(scan/deltax)
    each_change = int(shift/deltax)

    #fig, ax = subplots()
    for i in range(step):
        lo = fl + shift * i
        hi = lo + scan
        resance =  lo + scan/2
        noise = np.random.normal(0,sigma,10000)
        noise[bin_index] += bin_pow

        scan_x = x[i*each_change:each_len+i*each_change]
        scan_result = noise[i*each_change:each_len+i*each_change]
        w = sum(sigle_power(scan_x,resance)/sigma)

        grand_array[i*each_change:each_len+i*each_change] = grand_array[i*each_change:each_len+i*each_change] + scan_result*w

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