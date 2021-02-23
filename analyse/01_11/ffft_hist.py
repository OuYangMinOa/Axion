from numpy import *
from matplotlib.pyplot import *
from scipy.stats import norm, chisquare, poisson, chi2
x = linspace(0,5,999)
noise = random.normal(10,10,999)
dx = x[1] - x[0]
hist(noise,100)
title("gaussion noise hist")

y = fft.fftshift(fft.fft(noise))* dx / sqrt(2*pi)
y = y.real
fx =fft.fftshift(fft.fftfreq(len(x),dx)/sqrt(2*pi))
figure()
plot(fx,y)
figure()
n, bins, patches  = hist(y,100,density=True)
(mu, sigma) = norm.fit(y)
xx_middle = 0.5 * (bins[1:] + bins[:-1])
y2 = norm.pdf( xx_middle, mu, sigma)
l = plot(xx_middle, y2, 'r--', linewidth=2)
title(r'$\mathrm{Histogram\ of\ noise:}\ \mu=%e,\ \sigma=%e$' %(mu, sigma))
# title("fft gaussion noise and abs")
# print(sqrt(mean(y**2)), mean(y))
print(f"chisqure = {chisquare(n,y2)[1]}")
show()
