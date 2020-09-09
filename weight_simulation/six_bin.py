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
bandwidth = 125e-6  # MHz
scan      = 5e-1     # MHz
fl        = 749     # MHz
fr        = 751     # MHz
shift     = 2e-2    # MHz

sigma     = kb*T*bandwidth*1e6/sqrt(N) # sigma

def sigle_power(x,f0):
    return 3*(f0/750)*(1-2*s11)/(1-s11) * (1/(1+4*QL*QL*pow(x/f0-1,2)))

def lortz(x,f0):
    out = 1/ (1+4*QL*QL*(x/f0-1)**2)
    return out #/ max(out)

def coadd(array,sub_len):
    new_grand = []
    for i in range(len(array)-sub_len):
        new_grand.append(sum(array[i:i+sub_len]))
    return new_grand
    
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

def six_bin():
    plot()
    x = linspace(fl,fr,10000)
    noise = np.random.normal(0,sigma,10000)

    grand_array = x*0          # [0,0,...,0] 
    weigh_array = x*0          # [0,0,...,0] 
    sigm1_array = x*0          # [0,0,...,0]

    nowei_array = x*0          # [0,0,...,0] # no weight
    addti_array = x*0          # [0,0,...,0] 
    sigm2_array = x*0          # [0,0,...,0] 

    deltax = x[1] - x[0]
    start_data = 5000#random.randint(5000,7500)
    x0 = x[start_data]
    bin_data = maxwell(x,x0)
    bin_data = bin_data / max(bin_data) * 3
    bin_pos = x[argmax(bin_data)]
    bin_pow = max(bin_data)#/2#abs(noise[bin_index])

    step = int((fr-fl-scan)/shift)
    each_len = int(scan/deltax)
    each_change = int(shift/deltax)
    #fig, ax = subplots()
    for i in range(step):
        lo = fl + shift * i
        hi = lo + scan
        resance =  lo + scan/2
        noise = np.random.normal(0,sigma,10000)
        noise = noise + bin_data*lortz(x,resance)

        scan_x = x[i*each_change:each_len+i*each_change]
        scan_result = noise[i*each_change:each_len+i*each_change]

        sigma_this = sqrt(mean(scan_result**2))
        w = sigle_power(scan_x,resance)/sigma_this**2
        sig  = (sigle_power(scan_x,resance)/sigma_this)**2# sigma

        grand_array[i*each_change:each_len+i*each_change] = grand_array[i*each_change:each_len+i*each_change] + scan_result*w
        weigh_array[i*each_change:each_len+i*each_change] = weigh_array[i*each_change:each_len+i*each_change] + w#[w for i in range(len(scan_result))]
        sigm1_array[i*each_change:each_len+i*each_change] = sigm1_array[i*each_change:each_len+i*each_change] + sig
        
        nowei_array[i*each_change:each_len+i*each_change] = nowei_array[i*each_change:each_len+i*each_change] + scan_result
        addti_array[i*each_change:each_len+i*each_change] = addti_array[i*each_change:each_len+i*each_change] + [1 for i in range(len(scan_result))]
        sigm2_array[i*each_change:each_len+i*each_change] = sigm2_array[i*each_change:each_len+i*each_change] + [sigma_this**2 for i in range(len(scan_result))]

        if (1):
            rect = Rectangle((lo,0),scan,bin_pow,fill=None,edgecolor="r",linestyle="--")
            ax = subplot(121)
            cla()
            title(f"Max signal : {bin_pos:6.3f}(MHz)")
            plot(x,bin_data,label="signal")
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
            if (i*each_change<=start_data<each_len+i*each_change):
                ax2.arrow(x[start_data],7,0,-1,color="r",head_width=0.04, head_length=0.4)
            ylim(-5,10)
            savefig(f"six_bin/{i}.png")
            pause(0.001)
                # subplot(212)
                # plot([bin_pos,bin_pos],[0,bin_pow])
    timing = 1
    grand_array = grand_array/sqrt(sigm1_array)
    nowei_array = nowei_array/sqrt(sigm2_array)
    judge = max(grand_array)+0.1

    figure()
    gx = subplot(121)
    title("no weight")
    plot(x,nowei_array)
    subplot(122,sharex=gx,sharey=gx)
    title("weighted")
    plot(x,grand_array)
    savefig(f"six_bin_weight.png")

    subplots()
    plot(x,bin_data)
    
    subplots()
    grand_array = coadd(grand_array,6)
    nowei_array = coadd(nowei_array,6)

    figure()
    gx = subplot(121)
    title("no weight(coadd)")
    plot(x[:-6],nowei_array)
    subplot(122,sharex=gx,sharey=gx)
    title("weighted(coadd)")
    plot(x[:-6],grand_array)
    savefig(f"six_bin_weight_coadd.png")



    figure()
    candidate = []
    for _ in range(0):
        judge = max(grand_array)+0.1
        while len(candidate)<5:
            for i in range(len(grand_array)):
                if (grand_array[i] > judge):
                    print(x[i],grand_array[i],judge)
                    candidate.append(x[i])
                    if (i-8>0 and i+8<len(x)):
                        for each in range(i-8,i+8):
                            grand_array[each] = -10
                    elif (i-8<0):
                        for each in range(0,i+8):
                            grand_array[each] = -10
                    elif (i+8>len(x)):
                        for each in range(i-8,len(x)):
                            grand_array[each] = -10
            cla()
            title(f"judge: {judge:8.5f} candidate num : {len(candidate)}/5")
            plot([0,len(grand_array)],[judge,judge],"r--")
            plot(grand_array)
            xlabel("Frequency[MHz]")
            savefig(f"six_judge\\{timing}.png")
            timing += 1
            judge -= 0.1*sigma
            pause(0.1)
    show()



#single_bin()
six_bin()