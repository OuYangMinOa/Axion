println("Program start")
t  = 3600   # one hour (s)
h  = 6.626e-34 # J/k
ħ  = 4.135e-15/(2π) 
ρₐ = 0.45e15   # axion density (GeV/cc)
T  = 2      # physical + noise (k) 
Λ  = 78e-3  # Gev  78Mev
C  = 3e8    # The Speed of light
s  = 0      # S₁₁
Δv = 50e3   # 50kHz
f  = 5e9    # 5Ghz Frequency (10⁹ Hz)
B  = 9      # Magnetic (T)
V  = 1e-3   # Volume (L)
Q  = 5e4    # Q factor 

cm = 0.5
α  = 1/137  # fine structure const
×  = *
μ₀ = 4π × 1e-7 
ωᵪ = 2π*(f)
mₐ = ħ * ωᵪ # println(mₐ)
p_αγγ = 1.0
σ =  h * f *  sqrt(Δv/t)
println(σ)
shift =  ((ħ*C)^3)*ρₐ / (mₐ^2) * 0.5 * (1/μ₀) * (B^2) * V * ωᵪ* 0.5 * Q
println(shift)
g_γ = sqrt( 1.645*σ/shift)
println("\n\nf = $f    (Frequency)")
println("B = $B    (Magnetic)")
println("V = $V    (Cavity Volume)")
println("Q = $Q    (Q factor)")
println("\ng_αγγ = $(g_γ *1e9) GeV^-1")
