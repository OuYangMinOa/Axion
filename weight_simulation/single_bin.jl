using Distributed,Pkg
if (nprocs()<2)
    addprocs(3, exeflags="--project=.")
    Pkg.activate(".")
    Pkg.instantiate()
    current_dir=pwd()
    @everywhere current_dir=$current_dir
    @everywhere import Pkg
    @everywhere Pkg.activate(current_dir)

    using Distributions
end


# addprocs(7, exeflags="--project=.")
@everywhere using Distributed
@everywhere using PyPlot
@everywhere using Printf
@everywhere using CurveFit
@everywhere using CSV, DataFrames
@everywhere using Random
@everywhere using Distributions

@everywhere begin
    pygui(true)
    global c         = 3e8     # speed of light km/s
    global s11       = 0       # s11
    global QL        = 10_000  # Q factor
    global N         = 662625  # Integration time
    global T         = 0.3     # (k)
    global kb        = 1.38e-1 # 1e-22 
    global scan      = 1.8     # MHz
    global fl        = 4.9875e3     # MHz
    global fr        = 5.0125e3   # MHz
    global B         = 9
    global V         = 48338e-9
    global cnml      = 0.5
    global beta      = 0.5
    global bandwidth = 1125e-6  #(fr-fl)/(N-1)  # MHz
    global griding   = trunc(Int, (fr-fl+scan)/bandwidth)+1
    global shift = 0.25
    global sigma = kb*T*bandwidth*1e6 / sqrt(N)
    global kb_T_v = kb*T*bandwidth
end
println("griding = $griding")
@everywhere function expected_power(fe, g_gamma)
    f = fe * 1e6
    return 6.363e30*g_gamma^2/(2*pi*f)*B^2*QL*beta  
end

@everywhere function lortz(x, f0)
    out = 1 ./ (1. .+ 4*QL*QL*(x/f0.-1).^2)
    return out
end

@everywhere function single_bin(plotting = true,number=1, SNR=2)

    g_gamma = 3.27e-14 * SNR^0.5

    function sigle_power(x,f)
        return expected_power(f,g_gamma) .* (x.*0. .+1) .* lortz(x,f)
        
    end
    
    half_windows = scan/2
    real_x = collect(range(fl-half_windows,stop=fr+half_windows,length=griding))
    show_x = findall(x->x>fl&&x<fr,real_x)
    nowei_array = real_x.*0          # [0,0,...,0] # no weight
    addti_array = real_x.*0          # [0,0,...,0] 
    sigm1_array = real_x.*0          # [0,0,...,0] 

    grand_array = real_x.*0          # [0,0,...,0] # expect weight
    weigh_array = real_x.*0          # [0,0,...,0] 
    sigm2_array = real_x.*0          # [0,0,...,0] 

    grhay_array = real_x.*0          # [0,0,...,0] # expect weight
    wehay_array = real_x.*0          # [0,0,...,0] 
    wehay2_array =real_x.*0          # [0,0,...,0] 
    sigm3_array = real_x.*0          # [0,0,...,0] 

    dx = real_x[2] - real_x[1]

    step = trunc(Int,(fr-fl)/shift)
    each_len = trunc(Int,scan/dx)
    each_change = trunc(Int,shift/dx)

    noise = rand(Normal(0,sigma), 100)

    bin_index = rand(show_x,1)[1]
    bin_pos = real_x[bin_index][1]
    bin_pow = expected_power((fl+fr)/2, g_gamma)
    # println(bin_index,bin_pos)

    for i in 0:step - 1

        lo = fl + shift * i- half_windows
        hi = lo + scan
        scan_x = real_x[i*each_change+1:each_len+i*each_change+1]
        resance = scan_x[div(length(scan_x),2)]
        
        noise = rand(Normal(0,sigma),length(scan_x))

        

        if (minimum(scan_x)<=bin_pos<=maximum(scan_x))
            noise[bin_index-i*each_change] += bin_pow*lortz(resance,bin_pos)  # add the bin in noise 
        end # end if
        scan_result = noise    

        sigma_this = sqrt(mean(scan_result.^2))

        w    = sigle_power(scan_x,resance)/(sigma_this.^2) # weight
        sig  = ( sigle_power(scan_x,resance)/sigma_this).^2 # sigma# sigle_power(scan_x,resance)/sigma_this

        w2    = sigle_power(scan_x,resance).^2/(sigma_this).^2 # weight
        sig2  = ( w2./kb_T_v^2) # sigma

        
        grand_array[i*each_change + 1:each_len+i*each_change+ 1] = grand_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ scan_result.*w
        weigh_array[i*each_change+ 1:each_len+i*each_change+ 1] = weigh_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ w #[w for i in range(len(scan_result))]
        sigm1_array[i*each_change+ 1:each_len+i*each_change+ 1] = sigm1_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ sig

        grhay_array[i*each_change+ 1:each_len+i*each_change+ 1] = grhay_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ scan_result.*w2./sigle_power(scan_x,resance)
        wehay_array[i*each_change+ 1:each_len+i*each_change+ 1] = wehay_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ w2 #[w for i in range(len(scan_result))]
        wehay2_array[i*each_change+ 1:each_len+i*each_change+ 1] = wehay2_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ w2.^2 #[w for i in range(len(scan_result))]
        sigm3_array[i*each_change+ 1:each_len+i*each_change+ 1] = sigm3_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ sig##*kb_T_v**2
        
        nowei_array[i*each_change+ 1:each_len+i*each_change+ 1] = nowei_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ scan_result
        addti_array[i*each_change+ 1:each_len+i*each_change+ 1] = addti_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ [1 for i in 1:length(scan_result)]
        sigm2_array[i*each_change+ 1:each_len+i*each_change+ 1] = sigm2_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ [sigma_this.^2 for i in 1:length(scan_result)]

    end # end for

    now_snr = (nowei_array./sqrt.(sigm2_array))[show_x]
    now_max = findall(x->x==maximum(now_snr),now_snr)[1]

    gra_snr = (grand_array./sqrt.(sigm1_array))[show_x]
    gra_max = findall(x->x==maximum(gra_snr),gra_snr)[1]

    gra2_snr = ((grhay_array./wehay_array)./sqrt.(sigm3_array./wehay2_array))[show_x]
    gra2_max = findall(x->x==maximum(gra2_snr),gra2_snr)[1]

    new_x = real_x[show_x]

    if (plotting)
        gx = subplot(131)
        title("no weighting")
        plot(new_x,now_snr)
        gx.set_facecolor("red")
        if (abs(new_x[now_max]-bin_pos)<=1)
            gx.set_facecolor("white")
        end

        gx2 = subplot(132)
        title("ADMX")
        plot(new_x,gra_snr)
        gx2.set_facecolor("red")
        if (abs(new_x[gra_max]-bin_pos)<=1)
            gx2.set_facecolor("white")
        end

        gx3 = subplot(133)
        title("HAYSTAC")
        plot(new_x,gra2_snr)
        gx3.set_facecolor("red")
        if (abs(new_x[gra2_max]-bin_pos)<=1)
            gx3.set_facecolor("white")
        end
    end # end if
    
    find_1 = false
    find_2 = false
    find_3 = false
    
    if (abs(new_x[now_max]-bin_pos)<=1)
        find_1 = true
    end
    if (abs(new_x[gra_max]-bin_pos)<=1)
        find_2 = true
    end
    if (abs(new_x[gra2_max]-bin_pos)<=1)
        find_3 = true
    end

    return find_1, find_2, find_3
end

@everywhere function merging(x, y, sub_len)
    this_new_x = []
    this_new_y = []
    for i = 0:div(length(x),sub_len)-2
        append!(this_new_x, mean(x[i*sub_len+1:(i+1)*sub_len+1] ))
        append!(this_new_y, sum(y[i*sub_len+1:(i+1)*sub_len+1]))
    end
    return this_new_x, this_new_y
end

@everywhere function maxwell(x,x0,v)
    x_len = griding
    xx = x .- x0
    out =  xx.*xx.*exp.(-xx.*xx/v/v)
    for i in 1:x_len
        if (x[i] < x0)
            out[i] = 0.0
        else
            break
        end # end if
    end # end for
    return out/sum(out)
end

@everywhere function coadd(array, sub_len)
    out = []
    for i = 1:length(array)-sub_len
        append!(out, sum(array[i:i+sub_len]))
    end
    return out
end

@everywhere function fwHM(x,y)
    max_y = maximum(y)
    max_y_index = argmax(y)
    ny = y .- max_y/2
    ny = abs.(ny)

    l = argmin(ny[1:max_y_index])
    r = max_y_index + argmin(ny[max_y_index:end])
    return x[r] - x[l] 
end 

@everywhere function six_bin_plot(number=1, SNR=1.645, v=0.005, M_num=6)

    g_gamma = 3.27e-14 * SNR^0.5* (bandwidth/(1125e-6))^(1/4)

    function sigle_power(x,f)
        return expected_power(f,g_gamma) .* (x.*0. .+1) .* lortz(x,f)
    end

    half_windows = scan/2
    real_x = collect(range(fl-half_windows,stop=fr+half_windows,length=griding))
    show_x = findall(x->x>fl&&x<fr,real_x)

    nowei_array = real_x.*0          # [0,0,...,0] # no weight
    addti_array = copy(nowei_array)         # [0,0,...,0] 
    sigm1_array = copy(nowei_array)         # [0,0,...,0] 

    grand_array = copy(nowei_array)          # [0,0,...,0] # expect weight
    weigh_array = copy(nowei_array)          # [0,0,...,0] 
    sigm2_array = copy(nowei_array)          # [0,0,...,0] 

    grhay_array = copy(nowei_array)          # [0,0,...,0] # expect weight
    wehay_array = copy(nowei_array)         # [0,0,...,0] 
    wehay2_array =copy(nowei_array)          # [0,0,...,0] 
    sigm3_array = copy(nowei_array)          # [0,0,...,0] 

    dx = real_x[2] - real_x[1]

    step = trunc(Int,(fr-fl)/shift)
    each_len = trunc(Int,scan/dx)
    each_change = trunc(Int,shift/dx)

    noise = rand(Normal(0,sigma), griding)
    bin_pos = rand(real_x[show_x],1)[1]
    whole_signal = maxwell(real_x,bin_pos,v)
    whole_signal_max = argmax(whole_signal)

    figure()
    println("sigma = ",sigma)
    subplot(121)
    title("noise Distributions")
    singma_num = SNR*sigma
    plot([singma_num,singma_num],[0,700],"r--")
    hist(noise,100)
    subplot(122)
    title("noise")
    plot(real_x,noise)
    # plot([fl,fr],[singma_num,singma_num],"r--")

    figure()
    
    plot(real_x, whole_signal)
    println("FWHM = $(fwHM(real_x, whole_signal))")
    bin_pos = real_x[whole_signal_max][1]

    for i in 0:step - 1

        scan_x = real_x[i*each_change+1:each_len+i*each_change+1]
        resance = scan_x[div(length(scan_x),2)]
        
        noise = rand(Normal(0,sigma),length(scan_x))
        noise = noise .+ whole_signal[i*each_change+1:each_len+i*each_change+1] .* sigle_power(scan_x,resance)
        
        scan_result = noise    

        sigma_this = sqrt(mean(scan_result.^2))

        w    = sigle_power(scan_x,resance)/(sigma_this.^2) # weight
        sig  = ( sigle_power(scan_x,resance)/sigma_this).^2 # sigma# sigle_power(scan_x,resance)/sigma_this

        w2    = sigle_power(scan_x,resance).^2/(sigma_this).^2 # weight
        sig2  = ( w2/kb_T_v.^2) # sigma

        
        grand_array[i*each_change + 1:each_len+i*each_change+ 1] .= grand_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ scan_result.*w
        weigh_array[i*each_change+ 1:each_len+i*each_change+ 1] .= weigh_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ w #[w for i in range(len(scan_result))]
        sigm1_array[i*each_change+ 1:each_len+i*each_change+ 1] .= sigm1_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ sig

        grhay_array[i*each_change+ 1:each_len+i*each_change+ 1] .= grhay_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ scan_result.*sigle_power(scan_x,resance).*kb_T_v./(sigma_this^2)
        wehay_array[i*each_change+ 1:each_len+i*each_change+ 1] .= wehay_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ w2 #[w for i in range(len(scan_result))]
        wehay2_array[i*each_change+ 1:each_len+i*each_change+ 1] .= wehay2_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ w2.^2 #[w for i in range(len(scan_result))]
        sigm3_array[i*each_change+ 1:each_len+i*each_change+ 1] .= sigm3_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ w2##*kb_T_v**2
            
        nowei_array[i*each_change+ 1:each_len+i*each_change+ 1] .= nowei_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ scan_result
        addti_array[i*each_change+ 1:each_len+i*each_change+ 1] .= addti_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ [1 for i in 1:length(scan_result)]
        sigm2_array[i*each_change+ 1:each_len+i*each_change+ 1] .= sigm2_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ [sigma_this.^2 for i in 1:length(scan_result)]
        
    end # end for

    new_x = real_x[show_x]

    find_1 = false
    find_2 = false
    find_3 = false

    # NO merging
    

    gra_snr = (grand_array./sqrt.(sigm1_array))[show_x]
    gra_max = argmax(gra_snr)



    if (abs(new_x[gra_max]-bin_pos)<=1)
        find_1 = true
    end
    


    now_snr = (nowei_array./sqrt.(sigm2_array))[show_x]
    now_max = argmax(now_snr)

    gra2_snr = (grhay_array./sqrt.(sigm3_array))[show_x]#((grhay_array./wehay_array)./sqrt.(sigm3_array./wehay2_array))[show_x]
    gra2_max = argmax(gra2_snr)

    # NO merging

        figure()
        gx = subplot(131)
        title("no weighting")
        plot(new_x,now_snr)
        gx.set_facecolor("red")
        if (abs(new_x[now_max]-bin_pos)<=1)
            gx.set_facecolor("white")
        end

        gx2 = subplot(132)
        title("ADMX")
        plot(new_x,gra_snr)
        gx2.set_facecolor("red")
        if (abs(new_x[gra_max]-bin_pos)<=1)
            gx2.set_facecolor("white")
        end

        gx3 = subplot(133)
        title("HAYSTAC")
        plot(new_x,gra2_snr)
        gx3.set_facecolor("red")
        if (abs(new_x[gra2_max]-bin_pos)<=1)
            gx3.set_facecolor("white")
        end
    #

    # merging

        gra_snr_merge_x, gra_snr_merge_y = merging(new_x,gra_snr,M_num)
        gra_max_merge = argmax(gra_snr_merge_y)
        
        if (abs(gra_snr_merge_x[gra_max_merge]-bin_pos)<=1)
            find_2 = true
        end

        now_snr_merge_x, now_snr_merge_y = merging(new_x,now_snr,M_num)
        now_max_merge = argmax(now_snr_merge_y)

        gra2_snr_merge_x, gra2_snr_merge_y = merging(new_x,gra2_snr,M_num)
        gra2_max_merge = argmax(gra2_snr_merge_y)
        # merging

        figure()
        gx = subplot(131)
        title("no weighting merge")
        plot(now_snr_merge_x, now_snr_merge_y)
        gx.set_facecolor("red")
        if (abs(now_snr_merge_x[now_max_merge]-bin_pos)<=1)
            gx.set_facecolor("white")
        end

        gx2 = subplot(132)
        title("ADMX merge")
        plot( gra_snr_merge_x, gra_snr_merge_y)
        gx2.set_facecolor("red")
        if (abs(gra_snr_merge_x[gra_max_merge]-bin_pos)<=1)
            gx2.set_facecolor("white")
        end

        gx3 = subplot(133)
        title("HAYSTAC merge")
        plot(gra2_snr_merge_x, gra2_snr_merge_y)
        gx3.set_facecolor("red")
        if (abs(gra2_snr_merge_x[gra2_max_merge]-bin_pos)<=1)
            gx3.set_facecolor("white")
        end
    #

    # combined

        gra_snr_combine_y = coadd(gra_snr,M_num)
        gra_snr_combine_x = new_x[1:end-M_num]
        gra_max_combine = argmax(gra_snr_combine_y)
        
        if (abs(gra_snr_combine_x[gra_max_combine]-bin_pos)<=1)
            find_3 = true
        end

        now_snr_combine_y = coadd(now_snr,M_num)
        now_snr_combine_x = gra_snr_combine_x
        now_max_combine = argmax(now_snr_combine_y)

        gra2_snr_combine_y = coadd(gra2_snr,M_num)
        gra2_snr_combine_x = gra_snr_combine_x
        gra2_max_combine = argmax(gra2_snr_combine_y)
        # combined

        figure()
        gx = subplot(131)
        title("no weighting combined")
        plot(now_snr_combine_x, now_snr_combine_y)
        gx.set_facecolor("red")
        if (abs(now_snr_combine_x[now_max_combine]-bin_pos)<=1)
            gx.set_facecolor("white")
        end

        gx2 = subplot(132)
        title("ADMX combined")
        plot( gra_snr_combine_x, gra_snr_combine_y)
        gx2.set_facecolor("red")
        if (abs(gra_snr_combine_x[gra_max_combine]-bin_pos)<=1)
            gx2.set_facecolor("white")
        end

        gx3 = subplot(133)
        title("HAYSTAC combined")
        plot(gra2_snr_combine_x, gra2_snr_combine_y)
        gx3.set_facecolor("red")
        if (abs(gra2_snr_combine_x[gra2_max_combine]-bin_pos)<=1)
            gx3.set_facecolor("white")
        end
    #

    figure()
    gx = subplot(131)
    plot(new_x,gra_snr./gra2_snr)

    gx = subplot(132)

    plot(gra_snr_merge_x,gra_snr_merge_y./gra2_snr_merge_y)
    gx = subplot(133)
    
    plot(gra_snr_combine_x,gra_snr_combine_y./gra2_snr_combine_y)
    display([(gra_snr./gra2_snr)',(gra_snr_merge_y./gra2_snr_merge_y)',(gra_snr_combine_y./gra2_snr_combine_y)'])
    return find_1, find_2, find_3
end

@everywhere function six_bin_no_plot(number=1, SNR=1.645, v=0.005, M_num=6)

    g_gamma = 3.27e-14 * SNR^0.5 * (bandwidth/(1125e-6))^(1/4)

    function sigle_power(x,f)
        return expected_power(f,g_gamma) .* (x.*0. .+1) .* lortz(x,f)
    end

    half_windows = scan/2
    real_x = collect(range(fl-half_windows,stop=fr+half_windows,length=griding))
    show_x = findall(x->x>fl&&x<fr,real_x)

    grand_array = real_x.*0         # [0,0,...,0] # expect weight
    weigh_array = copy(grand_array)          # [0,0,...,0] 
    sigm1_array = copy(grand_array)          # [0,0,...,0] 

    dx = real_x[2] - real_x[1]

    step = trunc(Int,(fr-fl)/shift)
    each_len = trunc(Int,scan/dx)
    each_change = trunc(Int,shift/dx)

    noise = rand(Normal(0,sigma), 100)
    bin_pos = rand(real_x[show_x],1)[1]
    whole_signal = maxwell(real_x,bin_pos,v)
    whole_signal_max = argmax(whole_signal)
    
    bin_pos = real_x[whole_signal_max][1]

    for i in 0:step - 1

        scan_x = real_x[i*each_change+1:each_len+i*each_change+1]
        resance = scan_x[div(length(scan_x),2)]
        
        noise = rand(Normal(0,sigma),length(scan_x))
        noise = noise .+ whole_signal[i*each_change+1:each_len+i*each_change+1] .* sigle_power(scan_x,resance)
        
        scan_result = noise    

        sigma_this = sqrt(mean(scan_result.^2))

        w    = sigle_power(scan_x,resance)/(sigma_this.^2) # weight
        sig  = ( sigle_power(scan_x,resance)/sigma_this).^2 # sigma# sigle_power(scan_x,resance)/sigma_this

        w2    = sigle_power(scan_x,resance).^2/(kb_T_v*sigma_this).^2 # weight
        sig2  = ( w2*sigma_this).^2 # sigma

        
        grand_array[i*each_change + 1:each_len+i*each_change+ 1] .= grand_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ scan_result.*w
        weigh_array[i*each_change+ 1:each_len+i*each_change+ 1] .= weigh_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ w #[w for i in range(len(scan_result))]
        sigm1_array[i*each_change+ 1:each_len+i*each_change+ 1] .= sigm1_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ sig

    end # end for

    new_x = real_x[show_x]

    find_1 = false
    find_2 = false
    find_3 = false

    # NO merging
    
    gra_snr = (grand_array./sqrt.(sigm1_array))[show_x]
    gra_max = argmax(gra_snr)

    if (abs(new_x[gra_max]-bin_pos)<=1)
        find_1 = true
    end

    # merging

    gra_snr_merge_x, gra_snr_merge_y = merging(new_x,gra_snr,M_num)
    gra_max_merge = argmax(gra_snr_merge_y)
    
    if (abs(gra_snr_merge_x[gra_max_merge]-bin_pos)<=1)
        find_2 = true
    end

    # combined

    gra_snr_combine_y = coadd(gra_snr,M_num)
    gra_max_combine = argmax(gra_snr_combine_y)
    
    if (abs(new_x[1:end-M_num][gra_max_combine]-bin_pos)<=1)
        find_3 = true
    end

    return find_1, find_2, find_3
end

@everywhere function rescan_six_bin_no_plot(change=0, SNR=1.645, change2=-1, v=0.0041, M_num=4)

    g_gamma = 3.27e-14 * SNR^0.5 * (bandwidth/(1125e-6))^(1/4)

    function sigle_power(x,f)
        return expected_power(f,g_gamma) .* (x.*0. .+1) .* lortz(x,f)
    end

    half_windows = scan/2
    real_x = collect(range(fl-half_windows,stop=fr+half_windows,length=griding))
    show_x = findall(x->x>fl&&x<fr,real_x)

    grand_array = real_x.*0         # [0,0,...,0] # expect weight
    weigh_array = copy(grand_array)          # [0,0,...,0] 
    sigm1_array = copy(grand_array)          # [0,0,...,0] 

    dx = real_x[2] - real_x[1]

    step = trunc(Int,(fr-fl)/shift)
    each_len = trunc(Int,scan/dx)
    each_change = trunc(Int,shift/dx)

    noise = rand(Normal(0,sigma), 100)
    bin_pos = rand(real_x[show_x],1)[1]
    whole_signal = maxwell(real_x,bin_pos,v)
    whole_signal_max = argmax(whole_signal)
    
    bin_pos = real_x[whole_signal_max][1]

    for i in 0:step - 1

        scan_x = real_x[i*each_change+1:each_len+i*each_change+1]
        resance = scan_x[div(length(scan_x),2)]
        
        noise = rand(Normal(0,sigma),length(scan_x))
        noise = noise .+ whole_signal[i*each_change+1:each_len+i*each_change+1] .* sigle_power(scan_x,resance)
        
        scan_result = noise    

        sigma_this = sqrt(mean(scan_result.^2))

        w    = sigle_power(scan_x,resance)/(sigma_this.^2) # weight
        sig  = ( sigle_power(scan_x,resance)/sigma_this).^2 # sigma# sigle_power(scan_x,resance)/sigma_this

        w2    = sigle_power(scan_x,resance).^2/(kb_T_v*sigma_this).^2 # weight
        sig2  = ( w2*sigma_this).^2 # sigma

        
        grand_array[i*each_change + 1:each_len+i*each_change+ 1] .= grand_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ scan_result.*w
        weigh_array[i*each_change+ 1:each_len+i*each_change+ 1] .= weigh_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ w #[w for i in range(len(scan_result))]
        sigm1_array[i*each_change+ 1:each_len+i*each_change+ 1] .= sigm1_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ sig

    end # end for

    function rescan(candidate_array,can_y,can_w,can_sig)
        function scan_spe(left, right)
            start_l = left ÷ each_change
            step = (right-left-each_len) ÷ each_change
  
            for i in start_l:start_l+step-1
                scan_x = real_x[i*each_change+1:each_len+i*each_change+1]
                resance = scan_x[div(length(scan_x),2)]
                noise = rand(Normal(0,sigma),length(scan_x))
                noise = noise .+ whole_signal[i*each_change+1:each_len+i*each_change+1] .* sigle_power(scan_x,resance)
                
                scan_result = noise    

                sigma_this = sqrt(mean(scan_result.^2))

                w    = sigle_power(scan_x,resance)/(sigma_this.^2) # weight
                sig  = ( sigle_power(scan_x,resance)/sigma_this).^2 # sigma# sigle_power(scan_x,resance)/sigma_this

                can_y[i*each_change + 1:each_len+i*each_change+ 1] .= can_y[i*each_change+ 1:each_len+i*each_change+ 1] .+ scan_result.*w
                can_w[i*each_change+ 1:each_len+i*each_change+ 1] .= can_w[i*each_change+ 1:each_len+i*each_change+ 1] .+ w #[w for i in range(len(scan_result))]
                can_sig[i*each_change+ 1:each_len+i*each_change+ 1] .= can_sig[i*each_change+ 1:each_len+i*each_change+ 1] .+ sig
            
            end # end for
        end
        global can_index = 1
        n = length(candidate_array)
        while (can_index<=n)
            global can_index
            left  = candidate_array[can_index] - each_len÷2
            right = candidate_array[can_index] + each_len÷2
            while (can_index<n && right >= candidate_array[can_index+1])
                global can_index
                can_index = can_index + 1
                right = candidate_array[can_index] + each_len÷2
            end # end while
            scan_spe(left, right)
            can_index = can_index + 1
        end # end while
    end # end function    

    function is_signal_can(candidate_index)
        for each in candidate_index
            if (abs(real_x[each]-bin_pos)<1)
                return true
            end
        end # end for
        return false
    end

    new_x = real_x[show_x]

    find_1 = false
    find_2 = false
    find_3 = false

    find_1_rescan = false
    find_2_rescan = false
    find_3_rescan = false

    nanlee(a) = filter(!isnan, a)
    rms(x) = if (nanlee(x)==[]) 0 else sqrt(mean(nanlee(x).^2)) end

    if (change == 0)
        cut_s = SNR
    else
        cut_s = change
    end

    function cutting(cut_array,cut_sigma)
        len_show = length(cut_array)
        rms_this = rms(cut_array)
        candidate_count = 0
        condidate_arr = []
        for i in 1:len_show
            if (cut_array[i]> rms_this * cut_sigma)
                append!(condidate_arr,i+each_len÷2)
                candidate_count += 1
            end # end if
        end # end for
        return candidate_count, condidate_arr
    end

    # NO merging
    
    gra_snr = (grand_array./sqrt.(sigm1_array))[show_x]
    candidate_number_1, gra_candidate = cutting(gra_snr, cut_s)
    find_1 = is_signal_can(gra_candidate)

    gra_y   = real_x .* 0
    gra_w   = copy(gra_y)
    gra_sig = copy(gra_y)
    rescan(gra_candidate,gra_y,gra_w,gra_sig)



    gra_snr_recan = (gra_y./sqrt.(gra_sig))[show_x]
    candidate_number_1_rescan, gra_candidate_rescan = cutting(gra_snr_recan, cut_s)
    find_1_rescan = is_signal_can(gra_candidate_rescan)
    # println(candidate_number_1,"\n",candidate_number_1_rescan)

    # merging
    
    
    # gra_snr_merge_x, gra_snr_merge_y = merging(new_x,gra_snr,M_num)
    # candidate_number_2, gra_merge_candidate = cutting(gra_snr_merge_y)

    # gra_merge_y   = real_x .* 0
    # gra_merge_w   = copy(gra_merge_y)
    # gra_merge_sig = copy(gra_merge_y)
    # rescan(gra_merge_candidate,gra_merge_y,gra_merge_w,gra_merge_sig)

    # gra_merge_snr_recan = (gra_merge_y./sqrt.(gra_merge_sig))[show_x]
    # candidate_number_2_rescan, gra_merge_candidate_rescan = cutting(gra_merge_snr_recan)

    # println(candidate_number_2,"\n",candidate_number_2_rescan)

    # combined

    if (change2 == 0)
        cut_s = SNR
    elseif (change !=-1)
        cut_s = change2
    end

    gra_snr_combine_y = coadd(gra_snr,M_num)
    candidate_number_3, gra_combine_candidate = cutting(gra_snr_combine_y, cut_s)
    find_3 = is_signal_can(gra_combine_candidate)
    gra_combine_y   = real_x .* 0
    gra_combine_w   = copy(gra_combine_y)
    gra_combine_sig = copy(gra_combine_y)
    rescan(gra_combine_candidate,gra_combine_y,gra_combine_w,gra_combine_sig)

    gra_combine_snr_recan = (gra_combine_y./sqrt.(gra_combine_sig))[show_x]
    gra_combine_snr_recan = coadd(gra_combine_snr_recan,M_num)
    # gra_combine_snr_recan = coadd(gra_combine_snr_recan,M_num)
    candidate_number_3_rescan, gra_combine_candidate_rescan = cutting(gra_combine_snr_recan, cut_s)
    find_3_rescan = is_signal_can(gra_combine_candidate_rescan)
    return candidate_number_1, candidate_number_1_rescan,find_1,find_1_rescan, candidate_number_3, candidate_number_3_rescan,find_3,find_3_rescan
end

@everywhere function rescan_six_bin_plot(change=1.645, SNR=5, change2=-1, v=0.0041, M_num=4)

    g_gamma = 3.27e-14 * SNR^0.5 * (bandwidth/(1125e-6))^(1/4)

    function sigle_power(x,f)
        return expected_power(f,g_gamma) .* (x.*0. .+1) .* lortz(x,f)
    end

    half_windows = scan/2
    real_x = collect(range(fl-half_windows,stop=fr+half_windows,length=griding))
    show_x = findall(x->x>fl&&x<fr,real_x)

    grand_array = real_x.*0         # [0,0,...,0] # expect weight
    weigh_array = copy(grand_array)          # [0,0,...,0] 
    sigm1_array = copy(grand_array)          # [0,0,...,0] 

    dx = real_x[2] - real_x[1]

    step = trunc(Int,(fr-fl)/shift)
    each_len = trunc(Int,scan/dx)
    each_change = trunc(Int,shift/dx)

    noise = rand(Normal(0,sigma), 100)
    bin_pos = rand(real_x[show_x],1)[1]
    whole_signal = maxwell(real_x,bin_pos,v)
    whole_signal_max = argmax(whole_signal)
    
    bin_pos = real_x[whole_signal_max][1]

    for i in 0:step - 1

        scan_x = real_x[i*each_change+1:each_len+i*each_change+1]
        resance = scan_x[div(length(scan_x),2)]
        
        noise = rand(Normal(0,sigma),length(scan_x))
        noise = noise .+ whole_signal[i*each_change+1:each_len+i*each_change+1] .* sigle_power(scan_x,resance)
        
        scan_result = noise    

        sigma_this = sqrt(mean(scan_result.^2))

        w    = sigle_power(scan_x,resance)/(sigma_this.^2) # weight
        sig  = ( sigle_power(scan_x,resance)/sigma_this).^2 # sigma# sigle_power(scan_x,resance)/sigma_this

        w2    = sigle_power(scan_x,resance).^2/(kb_T_v*sigma_this).^2 # weight
        sig2  = ( w2*sigma_this).^2 # sigma

        
        grand_array[i*each_change + 1:each_len+i*each_change+ 1] .= grand_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ scan_result.*w
        weigh_array[i*each_change+ 1:each_len+i*each_change+ 1] .= weigh_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ w #[w for i in range(len(scan_result))]
        sigm1_array[i*each_change+ 1:each_len+i*each_change+ 1] .= sigm1_array[i*each_change+ 1:each_len+i*each_change+ 1] .+ sig

    end # end for

    function rescan(candidate_array,can_y,can_w,can_sig)
        function scan_spe(left, right)
            start_l = left ÷ each_change
            step = (right-left-each_len) ÷ each_change
  
            for i in start_l:start_l+step-1
                scan_x = real_x[i*each_change+1:each_len+i*each_change+1]
                resance = scan_x[div(length(scan_x),2)]
                noise = rand(Normal(0,sigma),length(scan_x))
                noise = noise .+ whole_signal[i*each_change+1:each_len+i*each_change+1] .* sigle_power(scan_x,resance)
                
                scan_result = noise    

                sigma_this = sqrt(mean(scan_result.^2))

                w    = sigle_power(scan_x,resance)/(sigma_this.^2) # weight
                sig  = ( sigle_power(scan_x,resance)/sigma_this).^2 # sigma# sigle_power(scan_x,resance)/sigma_this

                can_y[i*each_change + 1:each_len+i*each_change+ 1] .= can_y[i*each_change+ 1:each_len+i*each_change+ 1] .+ scan_result.*w
                can_w[i*each_change+ 1:each_len+i*each_change+ 1] .= can_w[i*each_change+ 1:each_len+i*each_change+ 1] .+ w #[w for i in range(len(scan_result))]
                can_sig[i*each_change+ 1:each_len+i*each_change+ 1] .= can_sig[i*each_change+ 1:each_len+i*each_change+ 1] .+ sig
            
            end # end for
        end
        global can_index = 1
        n = length(candidate_array)
        while (can_index<=n)
            global can_index
            left  = candidate_array[can_index] - each_len÷2
            right = candidate_array[can_index] + each_len÷2
            while (can_index<n && right >= candidate_array[can_index+1])
                global can_index
                can_index = can_index + 1
                right = candidate_array[can_index] + each_len÷2
            end # end while
            scan_spe(left, right)
            can_index = can_index + 1
        end # end while
    end # end function    

    
    new_x = real_x[show_x]

    nanlee(a) = filter(!isnan, a)
    rms(x) = if (nanlee(x)==[]) 0 else sqrt(mean(nanlee(x).^2)) end

    
    if (change == 0)
        cut_s = SNR
    else
        cut_s = change
    end

    function cutting(cut_array, cut_sigma)
        len_show = length(cut_array)
        rms_this = rms(cut_array)
        candidate_count = 0
        condidate_arr = []
        for i in 1:len_show
            if (cut_array[i]> rms_this * cut_sigma)
                append!(condidate_arr,i+each_len÷2)
                candidate_count += 1
            end # end if
        end # end for
        return candidate_count, condidate_arr
    end

    # NO merging
    
    gra_snr = (grand_array./sqrt.(sigm1_array))[show_x]
    candidate_number_1, gra_candidate = cutting(gra_snr, cut_s)

    gra_y   = real_x .* 0
    gra_w   = copy(gra_y)
    gra_sig = copy(gra_y)
    rescan(gra_candidate,gra_y,gra_w,gra_sig)


    gra_snr_recan = (gra_y./sqrt.(gra_sig))[show_x]
    candidate_number_1_rescan, gra_candidate_rescan = cutting(gra_snr_recan, cut_s)

    println(candidate_number_1,"\n",candidate_number_1_rescan)

    # merging
    
    
    # gra_snr_merge_x, gra_snr_merge_y = merging(new_x,gra_snr,M_num)
    # candidate_number_2, gra_merge_candidate = cutting(gra_snr_merge_y)

    # gra_merge_y   = real_x .* 0
    # gra_merge_w   = copy(gra_merge_y)
    # gra_merge_sig = copy(gra_merge_y)
    # rescan(gra_merge_candidate,gra_merge_y,gra_merge_w,gra_merge_sig)

    # gra_merge_snr_recan = (gra_merge_y./sqrt.(gra_merge_sig))[show_x]
    # candidate_number_2_rescan, gra_merge_candidate_rescan = cutting(gra_merge_snr_recan)

    # println(candidate_number_2,"\n",candidate_number_2_rescan)

    # combined

    if (change2 == 0)
        cut_s = SNR
    elseif (change2 != -1)
        cut_s = change2
    end

    gra_snr_combine_y = coadd(gra_snr,M_num)
    gra_snr_combine_x = new_x[1:end-M_num]
    candidate_number_3, gra_combine_candidate = cutting(gra_snr_combine_y, cut_s)
    # println(length(nanlee(gra_combine_candidate)))
    # if (length(nanlee(gra_combine_candidate))==0)
    #     println("no bins pass",length(nanlee(gra_combine_candidate)))
    # end
    gra_combine_y   = real_x .* 0
    gra_combine_w   = copy(gra_combine_y)
    gra_combine_sig = copy(gra_combine_y)
    rescan(gra_combine_candidate,gra_combine_y,gra_combine_w,gra_combine_sig)

    gra_combine_snr_recan = (gra_combine_y./sqrt.(gra_combine_sig))[show_x]
    gra_combine_snr_recan = coadd(gra_combine_snr_recan,M_num)
    gra_combine_snr_recan_x = new_x[1:length(gra_combine_snr_recan)]
    # gra_combine_snr_recan = coadd(gra_combine_snr_recan,M_num)
    candidate_number_3_rescan, gra_combine_candidate_rescan = cutting(gra_combine_snr_recan, cut_s)
    
    println(candidate_number_3,"\n",candidate_number_3_rescan)
    # display(gra_combine_snr_recan)
    figure()
    title("SNR = $SNR Combine Oringin")
    plot(new_x[1:end-M_num],gra_snr_combine_y)
    plot([minimum(gra_snr_combine_x),maximum(gra_snr_combine_x)],[cut_s.*rms(gra_snr_combine_y),cut_s.*rms(gra_snr_combine_y)],label="cut on $cut_s")
    legend()
    figure()
    title("SNR = $SNR Combine rescan")
    plot(gra_combine_snr_recan_x,gra_combine_snr_recan)
    plot([minimum(gra_combine_snr_recan_x),maximum(gra_combine_snr_recan_x)],[cut_s.*rms(gra_combine_snr_recan),cut_s.*rms(gra_combine_snr_recan)],label="cut on $cut_s")
    legend()
    return candidate_number_1, candidate_number_1_rescan, candidate_number_3, candidate_number_3_rescan
end

@everywhere function six_bin(plotting = true,number=1, SNR=2, v = 0.00375, M_num=6)
    if (plotting)
        return six_bin_plot(number,SNR,v,M_num)
    else
        return six_bin_no_plot(number,SNR,v,M_num)
    end
end

@everywhere function RUN_SINGLE()

    @time single_bin(false,1,1)

    SRN_len = 100
    total = 200
    work_finish = 0
    file_name = "C:\\Users\\$(splitdir(homedir())[end])\\OneDrive - cc.ncu.edu.tw\\研究\\for_git\\weight_simulation\\julia_sim_data\\single\\single.csv"

    SNR = range(1, stop=5, length=SRN_len)
    f1_arr = Array{Float32,1}(undef,SRN_len)
    f2_arr = Array{Float32,1}(undef,SRN_len)
    f3_arr = Array{Float32,1}(undef,SRN_len)
    # println("SNR | f1  | f2  | f3  ")
    for j in 1:SRN_len
        snr = SNR[j]
        f1_total = 0
        f2_total = 0
        f3_total = 0
        for i in 1:total
            f1, f2, f3 = single_bin(false,1,snr)
            f1_total += f1
            f2_total += f2
            f3_total += f3
        end
        f1_arr[j] = f1_total/total
        f2_arr[j] = f2_total/total
        f3_arr[j] = f3_total/total
        println(snr," | ",f1_total/total," | "
        ,f2_total/total," | "
        ,f3_total/total)

    end
    df = DataFrame(snr = SNR, f1=f1_arr, f2=f2_arr, f3=f3_arr)
    CSV.write(file_name, df)
    println("Successful save in $file_name")
end

@everywhere function plot_single()
    file = "C:\\Users\\$(splitdir(homedir())[end])\\OneDrive - cc.ncu.edu.tw\\研究\\for_git\\weight_simulation\\julia_sim_data\\single\\single.csv"
    df = CSV.read(file,DataFrame)
    snr = df.snr
    f1 = df.f1
    f2 = df.f2
    f3 = df.f3
    figure()
    plot(snr,f1,label="no weighting")
    plot(snr,f2,label="ADMX")
    plot(snr,f3,label="HAYSTAC")
    legend()
end

@everywhere function RUN_SIX(v=0.01, M_num = 5)
    @time six_bin_no_plot(1,1,v,M_num)
    println("start ($v,$M_num)")

    SRN_len = 50
    total = 200
    file_name = "C:\\Users\\$(splitdir(homedir())[end])\\OneDrive - cc.ncu.edu.tw\\研究\\for_git\\weight_simulation\\julia_sim_data\\six\\SIX_$(v)_$(M_num)_$(bandwidth).csv"

    SNR = range(1, stop=5, length=SRN_len)
    f1_arr = Array{Float32,1}(undef,SRN_len)
    f2_arr = Array{Float32,1}(undef,SRN_len)
    f3_arr = Array{Float32,1}(undef,SRN_len)
    # println("SNR | f1  | f2  | f3  ")
    for j in 1:SRN_len
        snr = SNR[j]
        f1_total = 0
        f2_total = 0
        f3_total = 0
        for i in 1:total
            f1, f2, f3 = six_bin_no_plot(1,snr,v,M_num)
            f1_total += f1
            f2_total += f2
            f3_total += f3
        end
        f1_arr[j] = f1_total/total
        f2_arr[j] = f2_total/total
        f3_arr[j] = f3_total/total
        println(trunc(Int,snr*10)/10," | ",f1_total/total," | "
        ,f2_total/total," | "
        ,f3_total/total)
    end

    df = DataFrame(snr = SNR, f1=f1_arr, f2=f2_arr, f3=f3_arr)
    CSV.write(file_name, df)
    println("Successful save ($v,$M_num) in $file_name")
end


@everywhere function plot_six(v = 0.01,M_num =8,showing=true)
    
    save_file = "C:\\Users\\$(splitdir(homedir())[end])\\OneDrive - cc.ncu.edu.tw\\研究\\for_git\\weight_simulation\\julia_sim_data\\figure\\SIX_0.0041_14_0.0002.png"
    file = "C:\\Users\\$(splitdir(homedir())[end])\\OneDrive - cc.ncu.edu.tw\\研究\\for_git\\weight_simulation\\julia_sim_data\\six\\SIX_0.0041_14_0.0002.csv"
    df = CSV.read(file,DataFrame)
    snr = df.snr
    f1 = df.f1
    f2 = df.f2
    f3 = df.f3
    figure()
    title("Merge = $M_num")
    # ylim(0,0.6)
    plot(snr,f1,label="no merging")
    plot(snr,f2,label="merge")
    plot(snr,f3,label="combine")
    ylabel("%")
    xlabel("SNR")
    legend()
    savefig(save_file)
    if (!showing)
        close()
    end
end

function main1()
    v = 0.0041
    # title("combine")
    # file = "C:\\Users\\$(splitdir(homedir())[end])\\OneDrive - cc.ncu.edu.tw\\研究\\for_git\\weight_simulation\\julia_sim_data\\six\\SIX_0.0041_14_0.0002.csv"
    # df = CSV.read(file,DataFrame)
    # snr = df.snr
    # f3 = df.f3
    # plot(snr,f3,label="combine = 14")

    for M_num in [1,5,10,15,20,25,30]
        file = "C:\\Users\\$(splitdir(homedir())[end])\\OneDrive - cc.ncu.edu.tw\\研究\\for_git\\weight_simulation\\julia_sim_data\\six\\SIX_$(v)_$(M_num)_$(bandwidth).csv"
        df = CSV.read(file,DataFrame)
        snr = df.snr
        f3 = df.f3
        plot(snr,f3,label="combine = $M_num")
    end
    legend()
    ylabel("%")
    xlabel("SNR")
    grid()
end
# plot_six(0.1,900)
# @sync @distributed for i in 1:1
#     RUN_SIX(0.01,8)
# end

# @time RUN_SIX(v=0.0037)
function plot_num_total()
    file = "C:\\Users\\$(splitdir(homedir())[end])\\OneDrive - cc.ncu.edu.tw\\研究\\for_git\\weight_simulation\\julia_sim_data\\six\\SIX_0041_25.csv"
    df = CSV.read(file,DataFrame)
    mm = df.m
    f1 = df.f1
    f2 = df.f2
    f3 = df.f3
    figure()
    title("SNR = 5")
    plot(mm,f1,label="no merging")
    plot(mm,f2,label="merge")
    plot(mm,f3,label="combine")
    ylabel("%")
    xlabel("Merge number")
    legend()

end

@everywhere f1_arr = Array{Float32,1}(undef,21)
@everywhere f2_arr = Array{Float32,1}(undef,21)
@everywhere f3_arr = Array{Float32,1}(undef,21)

function num_total()
    @time six_bin_no_plot(1,1,0.0041,1)

    SRN_len = 21
    total = 1000
    file_name = "C:\\Users\\$(splitdir(homedir())[end])\\OneDrive - cc.ncu.edu.tw\\研究\\for_git\\weight_simulation\\julia_sim_data\\six\\SIX_0041_25.csv"

    M_num_array =  range(1, stop=40, length=SRN_len)
    
    # println("SNR | f1  | f2  | f3  ")
    @sync @distributed for j in 1:SRN_len
        num = trunc(Int,M_num_array[j])
        f1_total = 0
        f2_total = 0
        f3_total = 0
        for i in 1:total
            f1, f2, f3 = six_bin_no_plot(1,3,0041,num)
            f1_total += f1
            f2_total += f2
            f3_total += f3
        end
        f1_arr[j] = f1_total/total
        f2_arr[j] = f2_total/total
        f3_arr[j] = f3_total/total
        println(num," | ",f1_total/total," | "
        ,f2_total/total," | "
        ,f3_total/total)
    end
end
# df = DataFrame(m = M_num_array, f1=f1_arr, f2=f2_arr, f3=f3_arr)
# CSV.write(file_name, df)
# plot(M_num_array,f3_arr)


@everywhere function rescan_plan(change=0)
    six_bin_no_plot(change,1)
    onetime = @elapsed six_bin_no_plot(change,1)
    SRN_len = 200
    total = 200

    println("Predict time = $(onetime*SRN_len*total)s")
    
    file = "rescan_SIX_$(change)_$(bandwidth).csv"
    file_name = "C:\\Users\\$(splitdir(homedir())[end])\\OneDrive - cc.ncu.edu.tw\\研究\\for_git\\weight_simulation\\julia_sim_data\\six\\$file"
    SNR = range(1, stop=5, length=SRN_len)
    f1_arr = Array{Float32,1}(undef,SRN_len)
    f3_arr = Array{Float32,1}(undef,SRN_len)
    f1_arr_rescan = Array{Float32,1}(undef,SRN_len)
    f3_arr_rescan = Array{Float32,1}(undef,SRN_len)

    f1_can = Array{Float32,1}(undef,SRN_len)
    f3_can = Array{Float32,1}(undef,SRN_len)
    f1_can_rescan = Array{Float32,1}(undef,SRN_len)
    f3_can_rescan = Array{Float32,1}(undef,SRN_len)
    for j in 1:SRN_len
        snr = SNR[j]
        f1_total = 0
        f3_total = 0
        f1_total_rescan = 0
        f3_total_rescan = 0

        f1_can_total = 0
        f3_can_total = 0
        f1_can_total_rescan = 0
        f3_can_total_rescan = 0
        for i in 1:total
            c1,r1,f1,fr1,c2,r2,f2,fr2 = rescan_six_bin_no_plot(change,snr)
            f1_can_total += c1
            f3_can_total += c2
            f1_can_total_rescan += r1
            f3_can_total_rescan += r2
            if f1  f1_total += 1/c1 end
            if f2  f3_total += 1/c2 end
            if fr1 f1_total_rescan += 1/r1 end
            if fr2 f3_total_rescan += 1/r2 end
        end
        f1_arr[j] = f1_total/total
        f3_arr[j] = f3_total/total
        f1_arr_rescan[j] = f1_total_rescan/total
        f3_arr_rescan[j] = f3_total_rescan/total

        f1_can[j] = f1_can_total/total
        f3_can[j] = f3_can_total/total
        f1_can_rescan[j] = f1_can_total_rescan/total
        f3_can_rescan[j] = f3_can_total_rescan/total
    end 

    df = DataFrame(snr = SNR,c1=f1_can,r1=f1_can_rescan,f1=f1_arr,fr1=f1_arr_rescan,
    c2=f3_can,r2=f3_can_rescan,f2=f3_arr,fr2=f3_arr_rescan)
    CSV.write(file_name, df)
    println("Successful save in $file_name")
end

function rescan_plot(change=0)
    file = "rescan_SIX_$(change)_$(bandwidth).csv"
    file_name = "C:\\Users\\$(splitdir(homedir())[end])\\OneDrive - cc.ncu.edu.tw\\研究\\for_git\\weight_simulation\\julia_sim_data\\six\\$file"
    df = CSV.read(file_name,DataFrame)
    snr = df.snr
    c1 = df.c1
    r1 = df.r1
    f1 = df.f1*100
    fr1 = df.fr1*100
    figure()
    subplot(121)
    title("no merge\nAVG of passed candidate")
    xlabel("SNR")
    ylabel("Num")
    plot(snr,c1,label="Oringin")
    plot(snr,r1,label="Rescan")
    legend()
    subplot(122)
    title("no merge\nAVG of posibility")
    xlabel("SNR")
    ylabel("%")
    plot(snr,f1,label="Oringin")
    plot(snr,fr1,label="Rescan")
    legend()
    c2 = df.c2
    r2 = df.r2
    f2 = df.f2*100
    fr2 = df.fr2*100
    figure()
    subplot(121)
    title("combine\nAVG of passed candidate")
    xlabel("SNR")
    ylabel("Num")
    plot(snr,c2,label="Oringin")
    plot(snr,r2,label="Rescan")
    legend()
    subplot(122)
    title("combine\nAVG of posibility")
    xlabel("SNR")
    ylabel("%")
    plot(snr,f2,label="Oringin")
    plot(snr,fr2,label="Rescan")
    legend()
end

function SNR_CUT_SERVE()
    cut = collect(range(1,5,length=200))
    cut = cut[100:end]
    SNR = collect(range(1,5,length=200))

    Z_1 = Array{Float64,2}(undef,length(cut),length(SNR))
    Z_2 = Array{Float64,2}(undef,length(cut),length(SNR))
    Z_3 = Array{Float64,2}(undef,length(cut),length(SNR))
    Z_4 = Array{Float64,2}(undef,length(cut),length(SNR))
    X, Y = repeat(SNR', length(cut), 1), repeat(cut, 1, length(SNR))

    which = "combine"

    for each in 1:length(cut)
        file = "rescan_SIX_$(cut[each])_$(bandwidth).csv"
        file_name = "C:\\Users\\$(splitdir(homedir())[end])\\OneDrive - cc.ncu.edu.tw\\研究\\for_git\\weight_simulation\\julia_sim_data\\six\\$file"
        df = CSV.read(file_name,DataFrame)
        Z_1[each,:] = df.f2 .*100
        Z_2[each,:] = df.fr2 .*100

        Z_3[each,:] = df.c2
        Z_4[each,:] = df.r2
    end #end for
    figure()
    title("$which\nBefore rescan prob")
    color_max = maximum(Z_1)
    clev = range(minimum(Z_1),step=color_max/1000,color_max)
    contourf(X,Y,Z_1,clev,cmap="jet")
    xlabel("SNR")
    ylabel("CUt")
    ticks = range(minimum(Z_1),length = 10,color_max)
    colorbar(fraction=0.03, pad=0.04,ticks = ticks)

    figure()
    title("$which\nAfter rescan prob")
    color_max = maximum(Z_2)
    clev = range(minimum(Z_2),step=color_max/1000,color_max)
    contourf(X,Y,Z_2,clev,cmap="jet")
    xlabel("SNR")
    ylabel("CUt")
    ticks = range(minimum(Z_2),length = 10,color_max)
    colorbar(fraction=0.03, pad=0.04,ticks = ticks)

    figure()
    title("$which\nBefore rescan cand")
    color_max = maximum(Z_3)
    clev = range(minimum(Z_3),step=color_max/1000,color_max)
    contourf(X,Y,Z_3,clev,cmap="jet")
    xlabel("SNR")
    ylabel("CUt")
    ticks = range(minimum(Z_3),length = 10,color_max)
    colorbar(fraction=0.03, pad=0.04,ticks = ticks)

    figure()
    title("$which\nAfter rescan prob")
    color_max = maximum(Z_4)
    clev = range(minimum(Z_4),step=color_max/1000,color_max)
    contourf(X,Y,Z_4,clev,cmap="jet")
    xlabel("SNR")
    ylabel("CUt")
    ticks = range(minimum(Z_4),length = 10,color_max)
    colorbar(fraction=0.03, pad=0.04,ticks = ticks)
end

function poly(x,c)
    ans = Array{typeof(x[1]),1}(undef,length(x))
    for j = 1:length(x)
        out = 0
        for i in 1:length(c)
            out += x[j]^(i-1) * c[i]
        end
        ans[j] =  out
    end
    return ans
end

function construct_word(c)
    out = ""
    for i in 1:length(c)
        this_c = @sprintf("%.2f", c[i])
        if (i==1)
            out = string(out,"$(this_c)")
   
        # else if (i==length(c))
        #     out = string(out,"+")
        else
            out = string("$(this_c)x^$(i-1)+",out)
  
        end
    end
    return out
end

function find_curve(n = 4)

    cut = collect(range(1,5,length=200))
    SNR = collect(range(1,5,length=200))

    Z_1 = Array{Float64,2}(undef,length(cut),length(SNR))
    Z_2 = Array{Float64,2}(undef,length(cut),length(SNR))
    Z_3 = Array{Float64,2}(undef,length(cut),length(SNR))
    Z_4 = Array{Float64,2}(undef,length(cut),length(SNR))

    f1_curve = Array{Float64,1}(undef,length(SNR))
    fr1_curve = Array{Float64,1}(undef,length(SNR))
    f2_curve = Array{Float64,1}(undef,length(SNR))
    fr2_curve = Array{Float64,1}(undef,length(SNR))
    for each in 1:length(cut)
        file = "rescan_SIX_$(cut[each])_$(bandwidth).csv"
        file_name = "C:\\Users\\$(splitdir(homedir())[end])\\OneDrive - cc.ncu.edu.tw\\研究\\for_git\\weight_simulation\\julia_sim_data\\six\\$file"
        df = CSV.read(file_name,DataFrame)
        Z_1[each,:] = Float64.(df.f1)*100
        Z_2[each,:] = df.fr1*100
        Z_3[each,:] = df.f2*100
        Z_4[each,:] = df.fr2*100
    end
    for each_snr in 1:length(SNR)
        f1_curve[each_snr] = cut[argmax(Z_1[:,each_snr])]
        fr1_curve[each_snr] = cut[argmax(Z_2[:,each_snr])]
        f2_curve[each_snr] = cut[argmax(Z_3[:,each_snr])]
        fr2_curve[each_snr] = cut[argmax(Z_4[:,each_snr])]
    end
    
    figure()
    title("No merge org")

    plot(SNR,f1_curve)

    f1_a = poly_fit(SNR,f1_curve,n)
    println(f1_a)
    plot(SNR,poly(SNR,f1_a),label="$(construct_word(f1_a))")
    xlabel("SNR")
    ylabel("cut")
    ylim(1,5)
    legend()

    figure()
    title("No merge rescan")
    plot(SNR,fr1_curve)
    fr1_a = poly_fit(SNR,fr1_curve,n)
    plot(SNR,poly(SNR,fr1_a),label="$(construct_word(fr1_a))")
    xlabel("SNR")
    ylabel("cut")
    ylim(1,5)
    legend()

    figure()
    title("combine org")
    plot(SNR,f2_curve)
    f2_a = poly_fit(SNR,f2_curve,n)
    plot(SNR,poly(SNR,f2_a),label="$(construct_word(f2_a))")
    xlabel("SNR")
    ylabel("cut")
    ylim(1,5)
    legend()

    figure()
    title("combine rescan")
    plot(SNR,fr2_curve)
    fr2_a = poly_fit(SNR,fr2_curve,n)
    plot(SNR,poly(SNR,fr2_a),label="$(construct_word(fr2_a))")
    xlabel("SNR")
    ylabel("cut")
    ylim(1,5)
    
    legend()
    # plot(SNR,poly(SNR,c))
    # fr1_a = poly_fit(SNR,fr1_curve,n)
    # f2_a = poly_fit(SNR,f2_curve,n)
    
end

@everywhere function rescan_plan_3D(change=0, change2=0)
    six_bin_no_plot(change,1,change2)
    onetime = @elapsed six_bin_no_plot(change,1,change2)
    SRN_len = 10
    total = 200

    println("Predict time = $(onetime*SRN_len*total)s")
    
    file = "rescan_SIX_$(change)_$(change2)_$(bandwidth).csv"
    file_name = "C:\\Users\\$(splitdir(homedir())[end])\\OneDrive - cc.ncu.edu.tw\\研究\\for_git\\weight_simulation\\julia_sim_data\\rescan_3D\\$file"
    SNR = range(1, stop=5, length=SRN_len)
    f1_arr = Array{Float32,1}(undef,SRN_len)
    f3_arr = Array{Float32,1}(undef,SRN_len)
    f1_arr_rescan = Array{Float32,1}(undef,SRN_len)
    f3_arr_rescan = Array{Float32,1}(undef,SRN_len)

    f1_can = Array{Float32,1}(undef,SRN_len)
    f3_can = Array{Float32,1}(undef,SRN_len)
    f1_can_rescan = Array{Float32,1}(undef,SRN_len)
    f3_can_rescan = Array{Float32,1}(undef,SRN_len)
    for j in 1:SRN_len
        snr = SNR[j]
        f1_total = 0
        f3_total = 0
        f1_total_rescan = 0
        f3_total_rescan = 0

        f1_can_total = 0
        f3_can_total = 0
        f1_can_total_rescan = 0
        f3_can_total_rescan = 0
        for i in 1:total
            c1,r1,f1,fr1,c2,r2,f2,fr2 = rescan_six_bin_no_plot(change,snr,change2)
            f1_can_total += c1
            f3_can_total += c2
            f1_can_total_rescan += r1
            f3_can_total_rescan += r2
            if f1  f1_total += 1/c1 end
            if f2  f3_total += 1/c2 end
            if fr1 f1_total_rescan += 1/r1 end
            if fr2 f3_total_rescan += 1/r2 end
        end
        f1_arr[j] = f1_total/total
        f3_arr[j] = f3_total/total
        f1_arr_rescan[j] = f1_total_rescan/total
        f3_arr_rescan[j] = f3_total_rescan/total

        f1_can[j] = f1_can_total/total
        f3_can[j] = f3_can_total/total
        f1_can_rescan[j] = f1_can_total_rescan/total
        f3_can_rescan[j] = f3_can_total_rescan/total
    end 

    df = DataFrame(snr = SNR,c1=f1_can,r1=f1_can_rescan,f1=f1_arr,fr1=f1_arr_rescan,
    c2=f3_can,r2=f3_can_rescan,f2=f3_arr,fr2=f3_arr_rescan)
    CSV.write(file_name, df)
    println("Successful save in $file_name")
end


@everywhere function rescan_plan_3D_contourf()

    cut1 = collect(range(1.0,5.0,length=50))
    cut2 = collect(range(1.0,5.0,length=50))
    SNR = collect(range(1,5,length=50))

    Z_1 = Array{Float64,3}(undef,length(cut1),length(cut2),length(SNR))
    Z_2 = Array{Float64,3}(undef,length(cut1),length(cut2),length(SNR))
    Z_3 = Array{Float64,3}(undef,length(cut1),length(cut2),length(SNR))
    Z_4 = Array{Float64,3}(undef,length(cut1),length(cut2),length(SNR))

    X, Y = repeat(SNR', length(cut1), 1), repeat(cut1, 1, length(SNR))

    
    which = "combine"

    # color_max = maximum(Z_1)
    # clev = range(minimum(Z_1),step=color_max/1000,color_max)
    # contourf(X,Y,Z_1,clev,cmap="jet")
    # xlabel("SNR")
    # ylabel("CUt")
    # ticks = range(minimum(Z_1),length = 10,color_max)
    # colorbar(fraction=0.03, pad=0.04,ticks = ticks)

    for change2 in 1:length(cut2)
        for change in 1:length(cut1)
            file = "rescan_SIX_$(cut1[change])_$(cut2[change2])_$(bandwidth).csv"
            file_name = "C:\\Users\\$(splitdir(homedir())[end])\\OneDrive - cc.ncu.edu.tw\\研究\\for_git\\weight_simulation\\julia_sim_data\\rescan_3D\\$file"
            df = CSV.read(file_name,DataFrame)
            Z_2[change,change2,:] = df.fr2*100
            Z_4[change,change2,:] = df.r2
        end #end for
    end #end for

    figure()
    
    color_max = maximum(Z_2)
    color_min = minimum(Z_2)
    clev = range(color_min,step=color_max/10000,color_max)
    ticks = range(color_min,length=10,color_max)
    # colorbar(fraction=0.03, pad=0.04,ticks = ticks)

    Z_same = Array{Float64,2}(undef,length(cut1),length(cut2))

    for sss in 1:length(cut1)
        # file = "SNR=$(SNR[sss]).png"
        # SAVE_FLODER = "C:\\Users\\$(splitdir(homedir())[end])\\OneDrive - cc.ncu.edu.tw\\研究\\for_git\\weight_simulation\\julia_sim_data\\rescan_3D_contourf\\$file"
        Z_same[sss,:] = Z_2[sss,sss,:]
    end
    color_max = maximum(Z_same)
    color_min = minimum(Z_same)
    clev = range(color_min,step=color_max/1000,color_max)
    contourf(X,Y,Z_same,clev, cmap = "jet")
    ticks = range(color_min,length=10,color_max)
    colorbar(fraction=0.03, pad=0.04,ticks = ticks)
    #     color_max = maximum(Z_2)
    # color_min = minimum(Z_2)

    
    # ticks = range(color_min,length=10,color_max)

    # for change2 in 1:length(cut2)
    #     println(cut2[change2])
    #     # this_z = Z_2[:,change2,:].+cut2[change2]
    #     # color_max = maximum(this_z)
    #     # color_min = minimum(this_z)

    #     # clev = range(color_min,step=color_max/10000,color_max)
    #     contourf(X,Y,Z_2[:,change2,:].+cut2[change2],clev.+cut2[change2], cmap = "jet")
    # end #end for


end


@everywhere function rescan_plan_3D_scatter()

    cut1 = collect(range(1.0,5.0,length=50))
    cut2 = collect(range(1.0,5.0,length=50))
    SNR = collect(range(1,5,length=50))

    Z_1 = Array{Float64,3}(undef,length(cut1),length(cut2),length(SNR))
    Z_2 = Array{Float64,3}(undef,length(cut1),length(cut2),length(SNR))
    Z_3 = Array{Float64,3}(undef,length(cut1),length(cut2),length(SNR))
    Z_4 = Array{Float64,3}(undef,length(cut1),length(cut2),length(SNR))

    X, Y = repeat(SNR', length(cut1), 1), repeat(cut1, 1, length(SNR))

    which = "combine"
    fig = figure()
    # gca(projection="3d")
    ax = fig.add_subplot(111, projection="3d")
    # color_max = maximum(Z_1)
    # clev = range(minimum(Z_1),step=color_max/1000,color_max)
    # contourf(X,Y,Z_1,clev,cmap="jet")
    # xlabel("SNR")
    # ylabel("CUt")
    # ticks = range(minimum(Z_1),length = 10,color_max)
    # colorbar(fraction=0.03, pad=0.04,ticks = ticks)

    for change2 in 1:length(cut2)
        file = "rescan_SIX_$(cut1[change])_$(cut2[change2])_$(bandwidth).csv"
        file_name = "C:\\Users\\$(splitdir(homedir())[end])\\OneDrive - cc.ncu.edu.tw\\研究\\for_git\\weight_simulation\\julia_sim_data\\rescan_3D\\$file"
        df = CSV.read(file_name,DataFrame)
        Z_2[change,change2,:] = df.fr2*100
            Z_4[change,change2,:] = df.r2
    end #end for
end


# find_curve()
# 200 [1,5,10,15,20,25,30]
# 1125 [4,5,8,10]
# @time begin
#     @sync @distributed for i in [4,5,8,10]
#         @async RUN_SIX(0.0041,i)
#     end
# end

@everywhere poss = []
for i in range(1,5,length=50)
    for j in range(1,5,length=50)
        append!(poss,[(i,j)])
    end 
end
@time begin
    @sync @distributed for (i,j) in poss
        @async rescan_plan_3D(i,j)
    end
    
end


# @sync @distributed for i in 25:25:200
#     @async RUN_SIX(0.05,i)
# end

# @time @sync @distributed for i in [1,2,3,5,6,7,9,10,11,13,14,15]
#     RUN_SIX(0.005,i)
# end

# @sync @distributed for i in 100:200:1000
#     RUN_SIX(0.5,i)
# end


# for i in [1.645,2.5,4,5]
#     rescan_six_bin_plot(3.5,i)
# end