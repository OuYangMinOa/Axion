println("### start ###")
using PyPlot
using PyCall
@pyimport  matplotlib.pyplot as pltt
pltt.style.use("dark_background")
pygui(true)

t  = 3600    # one hour (s)
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
g_target_14 = 3e-14 * π *0.62 * 1e7/ (α*m) # g_αγγ * pi * 0.62e7 / α * m

function N_to_SNR_1(SNR,g)
    top = (g/0.97)^2*(ρₐ/0.45)*(f/1)*(Q/10_000)*(c/0.5)*(1-2*s)/(1-s)*(B^2)*(V)
    low = SNR * T * t
    return (4.1e-4 * top / low)^-2
end

function N_SNR()
    SNR  = range(0,10,length = 100000)
    KSVZ_SNR = Array{Float64,1}(undef,length(SNR))
    DFSZ_SNR = Array{Float64,1}(undef,length(SNR))
    target_SNR_13 = Array{Float64,1}(undef,length(SNR))
    target_SNR_14 = Array{Float64,1}(undef,length(SNR))
    for i in 1:length(SNR)
        KSVZ_SNR[i] = N_to_SNR_1(SNR[i], g_KSVZ)
        DFSZ_SNR[i] = N_to_SNR_1(SNR[i], g_DFSZ)
        target_SNR_13[i] = N_to_SNR_1(SNR[i], g_target_13)
        target_SNR_14[i] = N_to_SNR_1(SNR[i], g_target_14)
        
    end
    yscale("log")
    grid()
    plot(SNR,KSVZ_SNR)
    plot(SNR,target_SNR_13,"g")
    plot(SNR,target_SNR_14,"r")
    #plot(SNR,DFSZ_SNR)
    title("N v.s. SNR")
    xlabel("SNR")
    ylabel("N")

    fill_between(SNR,KSVZ_SNR,target_SNR_14, alpha = 1, label = "KSVZ")
    
    fill_between(SNR,target_SNR_14,target_SNR_13, alpha = 1,color="r", label = L"target $g_{\alpha\gamma\gamma}=3x10^{-14}$")
    fill_between(SNR,target_SNR_13,minimum(target_SNR_13), alpha = 1,color="g", label = L"target $g_{\alpha\gamma\gamma}=3x10^{-13}$")
    legend()
end

N_SNR()