from matplotlib.pyplot import*
from matplotlib.patches import Rectangle
import sys
style.use("dark_background")
from numpy import *
import time
import threading

c         = 3e8     # speed of light km/s
s11       = 0       # s11
QL        = 10_000  # Q factor
N         = 662625  # Integration time
T         = 0.3     # (k)
kb        = 1.38e-1 # 1e-22 
scan      = 1.8     # MHz
fl        = 4.9875e3     # MHz
fr        = 5.0125e3   # MHz


bandwidth = 1125e-6  #(fr-fl)/(N-1)  # MHz
griding   = int((fr-fl)/bandwidth)
print(f"griding = {griding}")

shift     = 0.25    # MHz

sigma     = kb*T*bandwidth*1e6/sqrt(N) # sigma




kb_T_v = kb*T*bandwidth*1e6

print("sigma",sigma)



def single_bin(number,SNR= 2):
    g_gamma = 3.27e-14 * SNR**0.5
    B = 9
    V = 48338e-9
    cnml = 0.5
    beta = 0.5
    def expected_power(fe):
        f = fe * 1e6
        return 6.363e30*g_gamma**2/(2*pi*f)*B**2*QL*beta

    #print(expected_power(5e3))

    def sigle_power(x,f0):
        return expected_power(f0)* (x*0+1) *lortz(x,f0)

    def lortz(x,f0):
        out = 1/ (1+4*QL*QL*(x/f0-1)**2)
        return out #/ max(out)
      # x axis

    half_windows = scan/2
    real_x = linspace(fl-half_windows,fr+half_windows,griding)
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
    
    #random.randint(N-2*each_len,N-each_len) # randomly choose bin position
                # bin x axis
    #print("bin position",real_x[bin_index])
    #print("step=",step)

    noise = np.random.normal(0,sigma,griding) # noise
    bin_pos = random.choice(real_x[show_x])

    bin_index = argmin(abs(real_x-bin_pos))
    
    #print(bin_pos,bin_index)
    bin_pow = expected_power(5e3)
    
    #print(step)

    for i in range(step):
        lo = fl + shift * i - half_windows # left boundary 
        hi = lo + scan       # right boundary 

        scan_x = real_x[i*each_change:each_len+i*each_change]

        resance =  scan_x[len(scan_x)//2]#lo + scan/2  # resonance frequency


        if (1):
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

        w2    = sigle_power(scan_x,resance)**2/(kb_T_v*sigma_this**2)# weight
        sig2  = ( w2*sigma_this)**2# sigma

        grand_array[i*each_change:each_len+i*each_change] = grand_array[i*each_change:each_len+i*each_change] + scan_result*w
        weigh_array[i*each_change:each_len+i*each_change] = weigh_array[i*each_change:each_len+i*each_change] + w #[w for i in range(len(scan_result))]
        sigm1_array[i*each_change:each_len+i*each_change] = sigm1_array[i*each_change:each_len+i*each_change] + sig

        grhay_array[i*each_change:each_len+i*each_change] = grhay_array[i*each_change:each_len+i*each_change] + scan_result*w2*kb_T_v/sigle_power(scan_x,resance)
        wehay_array[i*each_change:each_len+i*each_change] = wehay_array[i*each_change:each_len+i*each_change] + w2 #[w for i in range(len(scan_result))]
        wehay2_array[i*each_change:each_len+i*each_change] = wehay2_array[i*each_change:each_len+i*each_change] + w2**2 #[w for i in range(len(scan_result))]
        sigm3_array[i*each_change:each_len+i*each_change] = sigm3_array[i*each_change:each_len+i*each_change] + sig##*kb_T_v**2
        
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
    else:
        gra_snr = (grand_array/sqrt(sigm1_array))[show_x]
        now_snr = (nowei_array/sqrt(sigm2_array))[show_x]

        # gra_snr = ((grhay_array/wehay_array)/sqrt(sigm3_array/wehay2_array))[show_x]
        # now_snr = (nowei_array/sqrt(sigm2_array))[show_x]

        new_x = real_x[show_x]
        max_pow = argmax(gra_snr)
    #print(abs(new_x[max_pow]-bin_pos))
    if_find = False
    second_max = sort(gra_snr)[-2]
    bigger = abs(second_max-gra_snr[max_pow])/gra_snr[max_pow]
    #print(second_max,gra_snr[max_pow],"bigger",bigger)
    gra_sigma =  sqrt(mean(gra_snr**2.0))
    # print(SNR, gra_snr[max_pow], gra_sigma)
    if_find_wrong = False

    area = 30
    condition = gra_sigma# gra_sigmamax(max(gra_snr[max_pow-area:max_pow]),max(gra_snr[max_pow+1:max_pow+area]))

    if (abs(new_x[max_pow]-bin_pos)<=0.1 and gra_snr[max_pow]>SNR*condition):
        #print(number,"find bin!!")
        if_find = True
        gx.set_facecolor("black")
        gx2.set_facecolor("black")
    elif (gra_snr[max_pow]>SNR*condition):
        if_find_wrong = True

    # print(bin_pos,gra_snr[max_pow],SNR*condition)
    show()
    #savefig(f"multiply/single_bin_{number}.png")
    return if_find, if_find_wrong

def six_bin(number,SNR= 2,MERGE_NUM=20,v = 0.005):
    g_gamma = 3.27e-14 * SNR**0.5
    B = 9
    V = 48338e-9
    cnml = 0.5
    beta = 0.5
    def expected_power(fe):
        f = fe * 1e6
        return 6.363e30*g_gamma**2/(2*pi*f)*B**2*QL*beta

    #print(expected_power(5e3))
    def FWHM(x,f):
        max_x = argmax(f)

        f = f - max(f)/2
        f = abs(f)

        l_min = argmin(f[:max_x])
        r_min = max_x+argmin(f[max_x:])
        print(l_min,r_min)
        return x[r_min] - x[l_min]

    def sigle_power(x,f0):
        return expected_power(f0)* (x*0+1) *lortz(x,f0)

    def lortz(x,f0):
        out = 1/ (1+4*QL*QL*(x/f0-1)**2)
        return out #/ max(out) 
      # x axis

    def coadd(array,sub_len):
        new_grand = []
        for i in range(len(array)-sub_len):
            new_grand.append(sum(array[i:i+sub_len]))
        return new_grand

    def merge_bin(x,y,sub_len):
        new_x = []
        new_y = []
        for i in range(len(x)//sub_len):
               new_x.append(mean(x[i*sub_len:(i+1)*sub_len])) 
               new_y.append(sum(y[i*sub_len:(i+1)*sub_len])) 
        return array(new_x), array(new_y)
    
    def maxwell(x,x0):
        norm = 4*pi / (v*v*pi)**1.5
        out = []
        for each in x:
            if each<x0:
                out.append(0)
            else:
                now_X = each - x0
                out.append(norm * now_X*now_X*exp(-now_X*now_X/v/v))
        return array(out)

    half_windows = scan/2
    real_x = linspace(fl-half_windows,fr+half_windows,griding)
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
    
    #random.randint(N-2*each_len,N-each_len) # randomly choose bin position
                # bin x axis
    #print("bin position",real_x[bin_index])
    #print("step=",step)

    noise = np.random.normal(0,sigma,griding) # noise
    bin_pos = np.random.choice(real_x[show_x])

    whole_signal = maxwell(real_x,bin_pos)
    #print(sum(whole_signal))
    whole_signal = whole_signal/ sum(whole_signal) * expected_power(real_x)
    
    bin_pos = real_x[argmax(maxwell(real_x,bin_pos))]

    # title(f"v={v} bin={bin_pos}")
    # xlabel("f[MHZ]")
    # ylabel("power[10e-22]")
    # noise = np.random.normal(0,sigma,len(real_x))
    # plot(real_x,whole_signal)
    # print("FWHM",FWHM(real_x,whole_signal))
    # show()
    

    bin_index = argmin(abs(real_x-bin_pos))
    
    #print(bin_pos,bin_index)
    bin_pow = max(whole_signal)
    
    #print(step)

    for i in range(step):
        lo = fl + shift * i - half_windows # left boundary 
        hi = lo + scan       # right boundary 

        scan_x = real_x[i*each_change:each_len+i*each_change]

        resance =  scan_x[len(scan_x)//2]#lo + scan/2  # resonance frequency


        if (1):
            noise = np.random.normal(0,sigma,len(scan_x)) # noise

            bin_signal = whole_signal[i*each_change:each_len+i*each_change]*sigle_power(scan_x,resance)# * lortz(scan_x,resance)

            noise = noise + bin_signal # add the bin in noise 
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
        sig2  = ( w2/kb_T_v**2)# sigma

        grand_array[i*each_change:each_len+i*each_change] = grand_array[i*each_change:each_len+i*each_change] + scan_result*w
        weigh_array[i*each_change:each_len+i*each_change] = weigh_array[i*each_change:each_len+i*each_change] + w #[w for i in range(len(scan_result))]
        sigm1_array[i*each_change:each_len+i*each_change] = sigm1_array[i*each_change:each_len+i*each_change] + sig

        grhay_array[i*each_change:each_len+i*each_change] = grhay_array[i*each_change:each_len+i*each_change] + scan_result*w2*kb_T_v/sigle_power(scan_x,resance)
        wehay_array[i*each_change:each_len+i*each_change] = wehay_array[i*each_change:each_len+i*each_change] + w2 #[w for i in range(len(scan_result))]
        wehay2_array[i*each_change:each_len+i*each_change] = wehay2_array[i*each_change:each_len+i*each_change] + w2**2 #[w for i in range(len(scan_result))]
        sigm3_array[i*each_change:each_len+i*each_change] = sigm3_array[i*each_change:each_len+i*each_change] + sig##*kb_T_v**2
        
        nowei_array[i*each_change:each_len+i*each_change] = nowei_array[i*each_change:each_len+i*each_change] + scan_result
        addti_array[i*each_change:each_len+i*each_change] = addti_array[i*each_change:each_len+i*each_change] + [1 for i in range(len(scan_result))]
        sigm2_array[i*each_change:each_len+i*each_change] = sigm2_array[i*each_change:each_len+i*each_change] + [sigma_this**2 for i in range(len(scan_result))]

        if (0):
            rect = Rectangle((lo,0),scan,bin_pow,fill=None,edgecolor="r",linestyle="--")
            ax = subplot(121)
            cla()
            title(f"signal : {bin_pos:8.5f}(MHz)")
            plot(real_x,whole_signal,label="signal")
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
    
    if (0):
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
        if_find_1 = False
        if (abs(new_x[argmax(now_snr)]-bin_pos)<=1):
            gx.set_facecolor("black")
            pass
        if (abs(new_x[argmax(gra_snr)]-bin_pos)<=1 ):
            #print(number,"find bin!!")
            if_find_1 = True
            gx2.set_facecolor("black")
    else:
        gra_snr = (grand_array/sqrt(sigm1_array))[show_x]
        now_snr = (nowei_array/sqrt(sigm2_array))[show_x]

        # gra_snr = ((grhay_array/wehay_array)/sqrt(sigm3_array/wehay2_array))[show_x]
        # now_snr = (nowei_array/sqrt(sigm2_array))[show_x]

        new_x = real_x[show_x]
        max_pow = argmax(gra_snr)

        if_find_1 = False
        if (abs(new_x[argmax(gra_snr)]-bin_pos)<=1 ):
            #print(number,"find bin!!")
            if_find_1 = True
    #print(abs(new_x[max_pow]-bin_pos))
    
    second_max = sort(gra_snr)[-2]
    bigger = abs(second_max-gra_snr[max_pow])/gra_snr[max_pow]
    #print(second_max,gra_snr[max_pow],"bigger",bigger)
    gra_sigma =  sqrt(mean(gra_snr**2.0))
    #print(SNR, gra_snr[max_pow], gra_sigma)



    #show()
    

    grand_array = coadd(gra_snr,MERGE_NUM)
    nowei_array = coadd(now_snr,MERGE_NUM)

    # figure()
    if (0):
        cla()
        gx = subplot(121)
        max_x = new_x[argmax(nowei_array)]
        title(f"no weight(coadd) {max_x}")

        plot(new_x[:-MERGE_NUM],nowei_array)
        print(max_x, bin_pos)
        if (abs(max_x-bin_pos)>1):
            gx.set_facecolor("red") 
        else:
            gx.set_facecolor("black")
        gx2 = subplot(122,sharex=gx,sharey=gx)
        max_x = new_x[argmax(grand_array)]
        title(f"weighted(coadd) {max_x}")
        print(max_x, bin_pos)
        plot(new_x[:-MERGE_NUM],grand_array)
        if_find_2 = False
        if (abs(max_x-bin_pos)>1):
            gx2.set_facecolor("red")
        else:
            gx2.set_facecolor("black")
            if_find_2 = True
    else:
        if_find_2 = False
        max_x = new_x[argmax(grand_array)]
        if (abs(max_x-bin_pos)<1):
            if_find_2 = True
    #show()


    xx,grand_array = merge_bin(new_x,gra_snr,MERGE_NUM)
    xx,nowei_array = merge_bin(new_x,now_snr,MERGE_NUM)

    #figure()
    if (0):
        cla()
        gx = subplot(121)
        max_x = xx[argmax(nowei_array)]
        title(f"no weight(merge) {max_x}")

        plot(xx,nowei_array)
        print(max_x, bin_pos)
        if (abs(max_x-bin_pos)>1):
            gx.set_facecolor("red") 
        else:
            gx.set_facecolor("black")
        gx2 = subplot(122,sharex=gx,sharey=gx)
        max_x = xx[argmax(grand_array)]
        title(f"weighted(merge) {max_x}")
        print(max_x, bin_pos)
        plot(xx,grand_array)
        if_find_3 = False
        if (abs(max_x-bin_pos)>1):
            gx2.set_facecolor("red")
        else:
            gx2.set_facecolor("black")
            if_find_3 = True
    else:
        if_find_3 = False
        max_x = xx[argmax(grand_array)]
        if (abs(max_x-bin_pos)<1):
            if_find_3 = True
    #show()
    #savefig(f"multiply/single_bin_{number}.png")
    del grand_array
    del nowei_array
    return if_find_1, if_find_2, if_find_3

def run_mutiple_time():
    if (1):
        num = 0
        #figure()
        find_prob = []
        total_time = 200
        sigma_x = np.linspace(1,5,100)
        for snr in sigma_x:
            num = 0
            for i in range(total_time):
                # print(i)
                a, b = single_bin(i,snr)
                if (a):
                    num += 1

            print(snr, num/total_time)
            find_prob.append(num/total_time)
        # print(snr, )
        np.save("sigma_prob_hay.npy",[sigma_x,find_prob])

def show_result():
    x1,y1 = np.load("sigma_prob_hay.npy")
    x2,y2 = np.load("sigma_prob.npy")

    g_gamma = lambda SNR:3.27e-14 * SNR **0.5

    fig = figure()
    ax1 = fig.add_subplot(111)
    title("single")
    ax1.set_xlabel("SNR")
    ax1.set_ylabel("%")
    xlim(1,5)
    ax1.plot(x1,100*y1,"-",label="HAYSTAC")
    ax1.plot(x2,100*y2,"-",label="ADMX")
    xticks_arr = array([1+i/2 for i in range(0,9)])
    ax1.set_xticks(xticks_arr)
    legend()
    grid()

    ax2 = ax1.twiny()
    g_x1 = g_gamma(x1)
    ax2.plot(g_x1,100*y1,alpha=0)
    ax2.set_xlim(g_gamma(array(ax1.get_xlim())))
    ax2.set_xticks(np.linspace(ax2.get_xticks()[0], ax2.get_xticks()[-1], len(ax1.get_xticks())))

    # ax2.set_xticks(g_gamma(linspace(1,5, len(ax1.get_xticks()) )))

    ax2.set_xlabel("g_gamma")
    tight_layout()
    print("stop")
    
    
    show()

def run_six_mutime( MERGE_NUM = 20,v=0.005):
    print(MERGE_NUM,v)
    if (1):
        num = 0
        num_com = 0
        num_mer = 0
        #figure()
        find = []
        find_com = []
        find_mer = []

        total_time = 200
        sigma_x = np.linspace(1,5,100)
        for snr in sigma_x:
            num = 0
            num_com = 0
            num_mer = 0
            for i in range(total_time):
                fi, fi2, fi3 = six_bin(i,snr,MERGE_NUM,v)
                if (fi):
                    num += 1
                if (fi2):
                    num_com += 1
                if (fi3):
                    num_mer += 1

            print(snr, num/total_time,num_com/total_time,num_mer/total_time)

            find.append(num/total_time)
            find_com.append(num_com/total_time)
            find_mer.append(num_mer/total_time)


        np.save(f"sixbin_data_{MERGE_NUM}_{v}.npy",[sigma_x,find,find_com,find_mer])
        print(f"sixbin_data_{MERGE_NUM}_{v}.npy saved")
        

def show_six_result(file_name):
    x1,find, find_com, find_mer = np.load(file_name)
    print(x1)
    g_gamma = lambda SNR:3.27e-14 * SNR **0.5

    fig = figure()
    ax1 = fig.add_subplot(111)
    title("single")
    ax1.set_xlabel("SNR")
    ax1.set_ylabel("%")
    xlim(1,5)
    ax1.plot(x1,100*find,"-",label="weighted")
    ax1.plot(x1,100*find_com,"-",label="combined")
    ax1.plot(x1,100*find_mer,"-",label="merge")
    xticks_arr = array([1+i/2 for i in range(0,9)])
    ax1.set_xticks(xticks_arr)
    legend()
    grid()

    ax2 = ax1.twiny()
    g_x1 = g_gamma(x1)
    ax2.plot(g_x1,100*find,alpha=0)
    ax2.set_xlim(g_gamma(array(ax1.get_xlim())))
    ax2.set_xticks(np.linspace(ax2.get_xticks()[0], ax2.get_xticks()[-1], len(ax1.get_xticks())))

    # ax2.set_xticks(g_gamma(linspace(1,5, len(ax1.get_xticks()) )))

    ax2.set_xlabel("g_gamma")
    tight_layout()
    print("stop")
    
    
    show()

if __name__=="__main__":
    before = time.time()
    # run_mutiple_time()
    # a, b = single_bin(0,1.5)
    # print(b)

    # show_six_result("sixbin_data_250_0.5.npy")

    # threading.Thread(target=run_mutiple_time ,daemon=False).start()
    # show_result()
    # six_bin(0)
    # single_bin(1,2.5)

    # v = eval(sys.argv[1])
    # print(v)
    # run_six_mutime(v,0.05)
    # six_bin(0,5,500,0.5)

    # v = 0.5#[0.05,0.005]
    # q = 250#[4,8,12,16,20]
    # run_six_mutime(q,v)
    
    # threading.Thread(target=run_six_mutime,daemon=False,args=(q,v)).start()
    # six_bin(0,5,1,0.005)
    print(time.time()-before)