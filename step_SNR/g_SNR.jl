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

scan = 5e4 
gaa_to_target(g_aa) = g_aa * π *0.62 * 1e7/ (α*m) # g_αγγ * pi / α

function N_to_SNR_1(SNR,g)
    top = (g/0.97)^2*(ρₐ/0.45)*(f/1)*(Q/10_000)*(c/0.5)*(1-2*s)/(1-s)*(B^2)*(V)
    low = SNR * T * t
    return scan / ((4.1e-4 * top / low)^-2)
end

function N_SNR()
    SNR_Length = 1000
    G_r_Length = 1000
    SNR  = range(1,8,length = SNR_Length)
    G_r  = 10. .^(range(-10,-15,length = G_r_Length))
    plate = Array{Float64,2}(undef,SNR_Length,G_r_Length)
    for i in 1:SNR_Length
        for j in 1:G_r_Length
            plate[i,j] = N_to_SNR_1(SNR[i],gaa_to_target(G_r[j]))
        end
    end
    
    X, Y = repeat(SNR', length(G_r), 1), repeat(G_r, 1, length(SNR))
    xlabel("SNR", labelpad=10)
    ylabel(L"$g_{\alpha\gamma\gamma}$", labelpad=10)
    #yscale("log")
    plate_log = log10.(plate)
    color_max = maximum(plate_log)
    clev = range(minimum(plate_log),step=color_max/1000,color_max)

    contourf(X,Y,plate_log,clev,cmap="jet")

    ticks = range(minimum(plate_log),length = 10,color_max)
    colorbar(fraction=0.03, pad=0.04,ticks = ticks)
    title(L"$SNR  v.s.  g_{\alpha\gamma\gamma}$")
end

N_SNR()
