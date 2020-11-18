from matplotlib.pyplot import*
from matplotlib.patches import Rectangle
style.use("dark_background")
from numpy import *

c         = 3e8     # speed of light km/s
s11       = 0       # s11
QL        = 10_000  # Q factor

N         = 662625  # Average time


T         = 0.3     # (k)
kb        = 1.38e-1 # 1e-22 
scan      = 1.8     # MHz
fl        = 4.95e3     # MHz
fr        = 5.05e3   # MHz

bin_pos = (fl + fr)/2
print(bin_pos)
bandwidth = 1125e-6  #(fr-fl)/(N-1)  # MHz
shift     = 0.25    # MHz

sigma     = kb*T*bandwidth*1e6/sqrt(N) # sigma


grid = 20_000

kb_T_v = kb*T*bandwidth*1e6

print("sigma",sigma)

SNR = 2
g_gamma = 3.27e-14 * SNR**0.5
B = 9
V = 48338e-9
cnml = 0.5
beta = 0.5
def expected_power(fe):
    f = fe * 1e6
    return 6.363e30*g_gamma**2/(2*pi*f)*B**2*QL*beta

print(expected_power(5e3))

def sigle_power(x,f0):
    return expected_power(f0)* (x*0+1) *lortz(x,f0)

def lortz(x,f0):
    out = 1/ (1+4*QL*QL*(x/f0-1)**2)
    return out #/ max(out)

def single_bin(number):
      # x axis
    half_windows = scan/2
    real_x = linspace(fl-half_windows,fr+half_windows,grid)
    show_x = where((fl<=real_x) & (real_x<=fr))

    x = real_x

    nowei_array = x*0          # [0,0,...,0] # no weight
    addti_array = x*0          # [0,0,...,0] 
    sigm1_array = x*0          # [0,0,...,0] 

    grand_array = x*0          # [0,0,...,0] # expect weight
    weigh_array = x*0          # [0,0,...,0] 
    sigm2_array = x*0          # [0,0,...,0] 

    grhay_array = x*0          # [0,0,...,0] # expect weight
    wehay_array = x*0          # [0,0,...,0] 
    wehay2_array = x*0          # [0,0,...,0] 
    sigm3_array = x*0          # [0,0,...,0] 

    
    dx = real_x[1] - real_x[0]
    step = int((fr-fl+scan)/shift)
    each_len = int(scan/dx)

    if (shift < dx):
        print("shift can't smaller than dx")
        raise ValueError("shift can't smaller than dx")
    each_change = int(shift/dx) # How many data point in a shift
    #print(step,each_len,each_change)
    bin_index = argmin(abs(real_x-bin_pos))
    #random.randint(N-2*each_len,N-each_len) # randomly choose bin position
                # bin x axis
    #print("bin position",real_x[bin_index])
    #print("step=",step)

    noise = np.random.normal(0,sigma,grid) # noise

    bin_pow = expected_power(5e3)
    
    #print(step)

    for i in range(step):
        lo = fl + shift * i - half_windows # left boundary 
        hi = lo + scan       # right boundary 

        scan_x = real_x[i*each_change:each_len+i*each_change]

        resance =  scan_x[len(scan_x)//2]#lo + scan/2  # resonance frequency


        if (0):

            noise = np.random.normal(0,sigma,len(scan_x)) # noise

            if (min(scan_x)<=bin_pos<=max(scan_x)):
                noise[bin_index-i*each_change] += bin_pow*lortz(resance,bin_pos)  # add the bin in noise 
            scan_result = noise
        else:
            noise = np.random.normal(0,sigma,len(scan_x)) # noise
            scan_result = noise*0

            for _ in range(N):
                noise = np.random.normal(0,sigma,len(scan_x)) # noise
                if (min(scan_x)<=bin_pos<=max(scan_x)):
                    noise[bin_index-i*each_change] += bin_pow*lortz(resance,bin_pos)  # add the bin in noise 
                scan_result += noise**2

            scan_result = sqrt(scan_result/N)  # root mean square

        sigma_this = sqrt(mean(scan_result**2))

        w    = sigle_power(scan_x,resance)/(sigma_this**2)# weight
        sig  = ( sigle_power(scan_x,resance)/sigma_this)**2# sigma# sigle_power(scan_x,resance)/sigma_this

        w2    = sigle_power(scan_x,resance)**2/(sigma_this**2)# weight
        sig2  = ( w2*sigma_this)**2# sigma

        grand_array[i*each_change:each_len+i*each_change] = grand_array[i*each_change:each_len+i*each_change] + scan_result*w
        weigh_array[i*each_change:each_len+i*each_change] = weigh_array[i*each_change:each_len+i*each_change] + w #[w for i in range(len(scan_result))]
        sigm1_array[i*each_change:each_len+i*each_change] = sigm1_array[i*each_change:each_len+i*each_change] + sig

        grhay_array[i*each_change:each_len+i*each_change] = grhay_array[i*each_change:each_len+i*each_change] + scan_result*w2*kb_T_v/sigle_power(scan_x,resance)
        wehay_array[i*each_change:each_len+i*each_change] = wehay_array[i*each_change:each_len+i*each_change] + w2 #[w for i in range(len(scan_result))]
        wehay2_array[i*each_change:each_len+i*each_change] = wehay2_array[i*each_change:each_len+i*each_change] + w2**2 #[w for i in range(len(scan_result))]
        sigm3_array[i*each_change:each_len+i*each_change] = sigm3_array[i*each_change:each_len+i*each_change] + sig2*kb_T_v**2
        
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

    # figure()
    # gx = subplot(121)
    # max_pow = argmax(nowei_array)
    # title(f"no weight Max : {x[max_pow]:8.5f}")
    
    # plot([bin_pos,bin_pos],[min(nowei_array/addti_array),max(nowei_array/addti_array)],"r--")
    # plot(x,nowei_array/addti_array)#/max(nowei_array/addti_array)*max(grand_array/weigh_array))

    # print(addti_array[-1])
    # xlabel("Frequency[MHz]")
    # ylabel("$power\t[10^{-22}]$")
    # subplot(122,sharex=gx,sharey=gx)
    
    # max_pow = argmax((grand_array/weigh_array)[each_len:-each_len])+each_len
    # title(f"weight Max : {x[max_pow]:8.5f}")
    # plot([bin_pos,bin_pos],[min(grand_array/weigh_array),max(grand_array/weigh_array)],"r--")
    # plot(x,grand_array/weigh_array)
    # xlabel("Frequency[MHz]")

    
    gra_snr = (grand_array)
    now_snr = (nowei_array)
    # figure()
    # gx = subplot(121)
    # max_pow = argmax(nowei_array)
    # title(f"no weight Max : {x[max_pow]:8.5f}")
    # new_x = x
    # plot(new_x,now_snr)
    
    # xlabel("Frequency[MHz]")
    # ylabel("SNR")
    # subplot(122,sharex=gx,sharey=gx)
    # plot(new_x,gra_snr)
    # max_pow = argmax(grand_array)#+each_len
    
    # title(f"weight Max : {new_x[max_pow]:8.5f}")
    # xlabel("Frequency[MHz]")
    if (1):
        gx = subplot(121)
        cla()
        max_pow = argmax(nowei_array)
        title(f"no weight Max : {x[max_pow]:8.5f}")
        gra_snr = (grand_array/sqrt(sigm1_array))[show_x]
        now_snr = (nowei_array/sqrt(sigm2_array))[show_x]
        new_x = real_x[show_x]
        plot(new_x,now_snr)
        xlabel("Frequency[MHz]")
        ylabel("SNR")
        gx2 = subplot(122,sharex=gx,sharey=gx)
        cla()
        plot(new_x,gra_snr)
        max_pow = argmax(gra_snr)#+each_len
        title(f"weight Max : {new_x[max_pow]:8.5f}")
        xlabel("Frequency[MHz]")
    else:
        gx = subplot(121)
        cla()
        max_pow = argmax(nowei_array)
        title(f"no weight Max : {x[max_pow]:8.5f}")
        gra_snr = ((grhay_array/wehay_array)/sqrt(sigm3_array/wehay2_array))[show_x]
        now_snr = (nowei_array/sqrt(sigm2_array))[show_x]
        new_x = real_x[show_x]
        plot(new_x,now_snr)
        xlabel("Frequency[MHz]")
        ylabel("SNR")
        gx2 = subplot(122,sharex=gx,sharey=gx)
        cla()
        plot(new_x,gra_snr)
        max_pow = argmax(gra_snr)#+each_len
        title(f"weight Max : {new_x[max_pow]:8.5f}")
        xlabel("Frequency[MHz]")
    gx.set_facecolor("red")
    gx2.set_facecolor("red")
    #print(abs(new_x[max_pow]-bin_pos))
    if_find = False
    second_max = sort(gra_snr)[-2]
    bigger = abs(second_max-gra_snr[max_pow])/gra_snr[max_pow]
    print(second_max,gra_snr[max_pow],"bigger",bigger)
    if (abs(new_x[max_pow]-bin_pos)<=0.1 and bigger>(0.1)):
        print(number,"find bin!!")
        if_find = True
        gx.set_facecolor("black")
        gx2.set_facecolor("black")

    show()
    savefig(f"multiply/single_bin_{number}.png")
    return if_find

if __name__=="__main__":
    num = 0
    figure()
    for i in range(1):
        if (single_bin(i)):
            num += 1
    print(num)
    print(num/768)
    print("stop")

