println("### start ###")
using PyPlot
plt.style.use("dark_background")
pygui(true)

scan_range = 50_000
step_size  = 2_000
Nₛ  = scan_range/step_size
Nₐ = 10_000  # N
t  = 3600    # one hour (s)
Δv = scan_range     #Nₐ/t    # Δv = N/t
ρₐ = 0.45    # axion density (GeV/cc)
f  = 5       # 5Ghz Frequency (10⁹ Hz)
m  = 2.06e-5 # ma (ev)
B  = 9       # Magnetic (T)
V  = 1       # Volume (L)
T  = 2       # physical + noise (k) 
Q  = 50_000  # Q factor 
c  = 0.5     # Cₘₙₚ
s  = 0       # S₁₁
α  = 1/137   # fine structure const

g_KSVZ = -0.94 # g_γ KSVZ
g_DFSZ =  0.36 # g_γ DFSZ

g_target_13 = 3e-13 * π *0.62 * 1e7/ (α*m) # g_αγγ * pi * 0.62e7 / α * m
#g_target_14 = 3e-14 * π *0.62 * 1e7/ (α*m) # g_αγγ * pi * 0.62e7 / α * m
g_target_14  = 9.979108e-15 * π *0.62 * 1e7/ (α*m)

function snr_to_t(SNR)
    up  = (g_KSVZ/0.97)^2*(ρₐ/0.45)*(f/1)*(Q/10_000)*(c/0.5)*(1-2*s)/(1-s)*(B^2)*(V)* √(Nₛ)
    low = SNR * T * √(Δv)
    return ((0.0017 * up / low)^-2)
end

function snr_to_na(SNR,g)
    up  = (g/0.97)^2*(ρₐ/0.45)*(f/1)*(Q/10_000)*(c/0.5)*(1-2*s)/(1-s)*(B^2)*(V)* √(Nₛ)
    low = SNR * T * Δv
    return ((4.1e-4 * up / low)^-2)
end

function snr_to_ns(SNR,g)
    up  = (g/0.97)^2*(ρₐ/0.45)*(f/1)*(Q/10_000)*(c/0.5)*(1-2*s)/(1-s)*(B^2)*(V)* √(Nₐ * scan_range)
    low = SNR * T * Δv
    return ((0.0017 * up / low)^2)
end

function t_SNR()
    SNR  = range(0,5,length = 1000)
    KSVZ_SNR = Array{Float64,1}(undef,length(SNR))
    for i in 1:length(SNR)
        KSVZ_SNR[i] = snr_to_t(SNR[i])
    end
    figure()
    #yscale("log")
    plot(SNR, KSVZ_SNR,label="Ns = $Nₛ")
    #plot(SNR, KSVZ_SNR.*0 .+ 3600,"r",label="3600s")
    plot(SNR.*0 .+1.654, KSVZ_SNR,"r",label=string(L"1.654$\sigma$," , "t=$(trunc(Int,snr_to_t(1.654)))"))
    fill_between(SNR,KSVZ_SNR,0)
    grid(linestyle="--")
    xticks(range(0,5,length = 11))
    yticks(range(0,maximum(KSVZ_SNR),length = 10))
    xlim(0,5)
    ylim(0,)
    ylabel("t (s)")
    xlabel("SNR")
    title("KSVZ  t v.s. SNR , Δv = $Δv")
    legend()
end

function Na_SNR()
    SNR  = range(0,7,length = 100000)
    KSVZ_SNR = Array{Float64,1}(undef,length(SNR))
    DFSZ_SNR = Array{Float64,1}(undef,length(SNR))
    target_SNR_13 = Array{Float64,1}(undef,length(SNR))
    target_SNR_14 = Array{Float64,1}(undef,length(SNR))
    for i in 1:length(SNR)
        KSVZ_SNR[i] = snr_to_na(SNR[i], g_KSVZ)
        DFSZ_SNR[i] = snr_to_na(SNR[i], g_DFSZ)
        target_SNR_13[i] = snr_to_na(SNR[i], g_target_13)
        target_SNR_14[i] = snr_to_na(SNR[i], g_target_14)
        
    end
    println(snr_to_ns(4, g_KSVZ))
    grid(zorder=0)
    plot(SNR,KSVZ_SNR,color="#2b79ff")
    plot(SNR,target_SNR_13,"#54ffaa")
    plot(SNR,target_SNR_14,"#ff5454")
    #plot(SNR,DFSZ_SNR)
    title(string(L"$N_a$ v.s. SNR , Ns =","$Nₛ , Δv = $Δv"))
    xlabel("SNR")
    ylabel(L"$N_a$")
    yscale("log")
    xlim(0,7)
    fill_between(SNR,KSVZ_SNR,target_SNR_14, alpha = 0.8,color="#2b79ff", label = "KSVZ")
    
    fill_between(SNR,target_SNR_14,target_SNR_13, alpha =  0.8,color="#ff5454", label = L"target $g_{\alpha\gamma\gamma}=9.97x10^{-15}$")
    fill_between(SNR,target_SNR_13,minimum(target_SNR_13), alpha =  0.8,color="#54ffaa", label = L"target $g_{\alpha\gamma\gamma}=3x10^{-13}$")
    legend()
end

function Ns_SNR()
    SNR  = range(0,7,length = 1000)
    KSVZ_SNR = Array{Float64,1}(undef,length(SNR))
    target_SNR_13 = Array{Float64,1}(undef,length(SNR))
    target_SNR_14 = Array{Float64,1}(undef,length(SNR))
    for i in 1:length(SNR)
        KSVZ_SNR[i] = snr_to_ns(SNR[i], g_KSVZ)
        target_SNR_13[i] = snr_to_ns(SNR[i], g_target_13)
        target_SNR_14[i] = snr_to_ns(SNR[i], g_target_14)
    end
    figure()
    yscale("log")
    grid(ls="-", color="0.8")
    plot(SNR,KSVZ_SNR,"k",color="#2b79ff")
    plot(SNR,target_SNR_13,"k",color="#54ffaa")
    plot(SNR,target_SNR_14,"k",color="#ff5454")
    #plot(SNR,DFSZ_SNR)
    xlim(0,7)
    fill_between(SNR,KSVZ_SNR,minimum(KSVZ_SNR), alpha = 0.8,color="#2b79ff", label = "KSVZ")
    fill_between(SNR,target_SNR_13,target_SNR_14, alpha = 0.8,color="#54ffaa", label = L"target $g_{\alpha\gamma\gamma}=3x10^{-13}$")
    
    fill_between(SNR,target_SNR_14,KSVZ_SNR, alpha = 0.8,color="#ff5454", label = L"target $g_{\alpha\gamma\gamma}=9.97x10^{-15}$")
    
    title(string("Step size v.s. SNR , Na =","$Nₐ , Δv = $Δv"))
    xlabel("SNR")
    ylabel("step(Hz)")
    legend()
end
plt.close()
plt.close()
plt.close()
Na_SNR()
Ns_SNR()
t_SNR()