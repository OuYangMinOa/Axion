import pandas as pd
from matplotlib.pyplot import *
import matplotlib.mlab as mlab
from numpy import *
from scipy.optimize import curve_fit
from scipy.special import factorial
import getpass
from scipy.stats import norm, chisquare, poisson, chi2


total_time = 5
y = np.random.normal(1,1,10000)

y[5000] += 1.645
(mu, sigma) = norm.fit(y)

y = (y - mu)/sigma
print(f"(mu, sigma) : ({mu}, {sigma})" )

n, bins, patches = hist(y,100,density=True)
title(r'$\mathrm{Histogram\ of\ noise:}\ \mu=%e,\ \sigma=%e$' %(mu, sigma))
x_middle = 0.5 * (bins[1:] + bins[:-1])

BMP_y_total = np.zeros(10000) + 1
BMP_y_sum = np.zeros(10000)
figure()
BMP_y = exp(mu*y -(mu)**2/2 )

plot(BMP_y)
figure()
subplot(total_time,2,1)
plot(y)
plot(5000,y[5000],"ro")
subplot(total_time,2,2)
plot(BMP_y)
plot(5000,BMP_y[5000],"ro")
BMP_y_total *= BMP_y
BMP_y_sum += y

for i  in range(total_time-1):
	y = np.random.normal(1,1,10000)
	y[5000] += 1.645
	(mu, sigma) = norm.fit(y)
	subplot(total_time,2,(i+2)*2-1)
	BMP_this = exp(mu*y -(mu)**2/2 )
	plot(y)
	plot(5000,y[5000],"ro")
	subplot(total_time,2,(i+2)*2)
	plot(BMP_this)
	plot(5000,BMP_this[5000],"ro")
	BMP_y_total *= BMP_this
	BMP_y_sum += y

figure()
plot(BMP_y_total)
plot(5000,BMP_y_total[5000],"ro")
figure()
plot(BMP_y_sum)
plot(5000,BMP_y_sum[5000],"ro")
show()