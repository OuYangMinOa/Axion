from matplotlib.pyplot import*
from matplotlib.patches import Rectangle
style.use("dark_background")
from numpy import *

c         = 3e5     # speed of light km/s
s11       = 0       # s11
QL        = 70_000  # Q factor
N         = 16_001  # Integration time
T         = 5.6     # (k)
kb        = 1.38e-1 # 1e-22
scan      = 5e-2    # MHz
fl        = 749     # MHz
fr        = 751     # MHz
bandwidth = 125e-6#(fr-fl)/(N-1)  # MHz
shift     = 2e-2    # MHz
sigma     = kb*T*bandwidth*1e6/sqrt(N-1) # sigma

print(bandwidth)

def sigle_power(x,f0):
    return 3.0*(f0/750.0)*(1-2*s11)/(1-s11)* (x*0+1) *lortz(x,f0)

def lortz(x,f0):
    out = 1/ (1+4*QL*QL*(x/f0-1)**2)
    return out #/ max(out)

def single_bin():
    x = linspace(fl,fr,N)  # x axis
    nowei_array = x*0          # [0,0,...,0] # no weight
    addti_array = x*0          # [0,0,...,0] 
    sigm1_array = x*0          # [0,0,...,0] 

    grand_array = x*0          # [0,0,...,0] # expect weight
    weigh_array = x*0          # [0,0,...,0] 
    sigm2_array = x*0          # [0,0,...,0] 

    deltax = x[1] - x[0]       # dx
    print(deltax)
    step = int((fr-fl-scan)/shift)  # How many step we need to scan full range
    each_len = int(scan/deltax) # How many data point in scan
    each_change = int(shift/deltax) # How many data point in a shift

    bin_index = 8150#random.randint(N-2*each_len,N-each_len) # randomly choose bin position
    bin_pos = x[bin_index]                # bin x axis

    noise = np.random.normal(0,sigma,N) # noise
    bin_pow = 5#max(noise) #abs(noise[bin_index]) # signal power
    
    print(step)
    #fig, ax = subplots()

    last_1, last_2 = [1 for i in range((each_len))], [1 for i in range((each_len))]

    for i in range(step):
        lo = fl + shift * i  # left boundary 
        hi = lo + scan       # right boundary 

        scan_x = x[i*each_change:each_len+i*each_change]

        resance =  scan_x[len(scan_x)//2]#lo + scan/2  # resonance frequency

        noise = np.random.normal(0,sigma,len(scan_x)) # noise
        if (i*each_change<bin_index<each_len+i*each_change):
            noise[bin_index-i*each_change] += bin_pow*lortz(resance,bin_pos)  # add the bin in noise 

        scan_result = noise#[i*each_change:each_len+i*each_change] # power get in this scan

        sigma_this = sqrt(mean(scan_result**2))

        w    = sigle_power(scan_x,resance)/(sigma_this**2)# weight

        sig  = (sigle_power(scan_x,resance)/sigma_this)**2# sigma# sigle_power(scan_x,resance)/sigma_this

        grand_array[i*each_change:each_len+i*each_change] = grand_array[i*each_change:each_len+i*each_change] + scan_result*w
        weigh_array[i*each_change:each_len+i*each_change] = weigh_array[i*each_change:each_len+i*each_change] + w #[w for i in range(len(scan_result))]
        sigm1_array[i*each_change:each_len+i*each_change] = sigm1_array[i*each_change:each_len+i*each_change] + sig

        nowei_array[i*each_change:each_len+i*each_change] = nowei_array[i*each_change:each_len+i*each_change] + scan_result
        addti_array[i*each_change:each_len+i*each_change] = addti_array[i*each_change:each_len+i*each_change] + [1 for i in range(len(scan_result))]
        sigm2_array[i*each_change:each_len+i*each_change] = sigm2_array[i*each_change:each_len+i*each_change] + [sigma_this**2 for i in range(len(scan_result))]

        if (0):
            rect = Rectangle((lo,0),scan,bin_pow,fill=None,edgecolor="r",linestyle="--")
            ax = subplot(121)
            cla()
            title(f"signal : {bin_pos:8.5f}(MHz)")
            plot([bin_pos,bin_pos],[0,bin_pow],label="signal")
            plot([resance,resance],[0,bin_pow],'r',label="resonace frequency")
            
            ylabel("$power [10^{-22}]$")
            xlabel("Frequency[MHz]")

            ax.add_patch(rect)

            plot(x,x*0,"black")
            legend()
            ax2 = subplot(122)
            cla()
            title(f"{lo} ~ {hi}")
            plot(scan_x,scan_result)
            if (i*each_change<=bin_index<each_len+i*each_change):
                ax2.arrow(x[bin_index],7,0,-1,color="r",head_width=scan/25, head_length=0.4)

            ylim(-8,12)
            savefig(f"subgle_bin/{i}.png")
            pause(0.01)

    figure()
    gx = subplot(121)
    max_pow = argmax(nowei_array)
    title(f"no weight Max : {x[max_pow]:8.5f}")
    
    plot([bin_pos,bin_pos],[min(nowei_array/addti_array),max(nowei_array/addti_array)],"r--")
    plot(x,nowei_array/addti_array)#/max(nowei_array/addti_array)*max(grand_array/weigh_array))

    print(addti_array[-1])
    xlabel("Frequency[MHz]")
    ylabel("$power\t[10^{-22}]$")
    subplot(122,sharex=gx,sharey=gx)
    
    max_pow = argmax((grand_array/weigh_array)[each_len:-each_len])+each_len
    title(f"weight Max : {x[max_pow]:8.5f}")
    plot([bin_pos,bin_pos],[min(grand_array/weigh_array),max(grand_array/weigh_array)],"r--")
    plot(x,grand_array/weigh_array)
    xlabel("Frequency[MHz]")

    figure()
    gx = subplot(121)
    max_pow = argmax(nowei_array)
    title(f"no weight Max : {x[max_pow]:8.5f}")
    gra_snr = (grand_array/sqrt(sigm1_array))
    now_snr = (nowei_array/sqrt(sigm2_array))#[each_len:-each_len]
    new_x = x#[each_len:-each_len]
    #plot([bin_pos,bin_pos],[min(now_snr),max(now_snr)],"r--")
    plot(new_x,now_snr)#/max(now_snr) * max(gra_snr))
    
    xlabel("Frequency[MHz]")
    ylabel("SNR")
    subplot(122,sharex=gx,sharey=gx)
    #plot([bin_pos,bin_pos],[min(gra_snr),max(gra_snr)],"r--")
    #[each_len:-each_len]
    plot(new_x,gra_snr)
    max_pow = argmax(gra_snr[each_len:-each_len])+each_len
    
    title(f"weight Max : {new_x[max_pow]:8.5f}")
    xlabel("Frequency[MHz]")
    savefig(f"result.png")

    show()

if __name__=="__main__":
    single_bin()



