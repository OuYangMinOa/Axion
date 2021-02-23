import pandas as pd
from matplotlib.pyplot import *
import matplotlib.mlab as mlab
from numpy import *
from scipy.optimize import curve_fit
from scipy.special import factorial
import getpass
from scipy.stats import norm, chisquare, poisson, chi2

file_name = "01_11.txt"
data = pd.read_csv(file_name, header = None)
print(data)
data = data.T


f = data[0]
y = data[1]

def fit_function(k, lamb,  scale):
    '''poisson function, parameter lamb is the fit parameter'''
    k2 = k
    return scale

def show_bis(plotting = True):
	basis = y
	inde = where(basis<-98.4)
	basis = basis.values[inde]
	basis = 10**(basis/10)
	(mu, sigma) = norm.fit(basis)
	print((mu, sigma))
	n, bins, patches  = hist(basis,100,density=True)
	x_middle = array([i for i in range(len(n))])
	xx_middle = 0.5 * (bins[1:] + bins[:-1])
	# popt, matrix = curve_fit(fit_function, x_middle, n)
	# print(*popt)
	# plot(xx_middle, fit_function(x_middle,*popt))
	# title(r"$ poisson\ distribution\ \lambda = %f, scale = %f $" % (popt[0],popt[1])))

	# print(chisquare([0,1,2],[1.2,2.5,3]))
	y2 = norm.pdf( xx_middle, mu, sigma)
	if (plotting):
		l = plot(xx_middle, y2, 'r--', linewidth=2)
		title(r'$\mathrm{Histogram\ of\ noise:}\ \mu=%e,\ \sigma=%e$' %(mu, sigma))
		# title(r"$\mathrm{Histogram\ of\ noise}$"))
		print(f"chisqure = {chisquare(n,y2)[1]}")

		xlabel("power[dbw]")
		ylabel("count")
		figure()
		plot(f.values[inde],basis)

		show()

def smaller_signal(plotting = True,SNR = 1.645):
	signal_pos = []
	start = 4
	endd = 15
	for i in range(start,endd):
		signal_pos.append(argmax(data[i]))
	signal_pos[-2] = (signal_pos[-3]+signal_pos[-1])//2
	basis = y
	inde = where(basis>-98.4)
	for i in inde[0]:
		basis.values[i]= sum(basis.values[i-10:i])/10
	new_f = f
	(mu, sigma) = norm.fit(basis)
	noise_sigma = sigma 
	np.random.shuffle(basis)
	all_data = zeros((endd-start,len(basis)))

	arg_signal_pos = len(signal_pos)//2

	popt = see_S21(plotting)
	for i in range(1,endd-start+1):
		np.random.shuffle(basis)
		all_data[i-1] = basis
		all_data[i-1][signal_pos[arg_signal_pos]] += func(f[signal_pos[i-1]],*popt)*noise_sigma
		if (plotting):
			plot(new_f,all_data[i-1]-2*(i-1))	

	no_weight = zeros(len(basis	))
	weighted = zeros(len(basis))
	if (plotting):
		show()
		figure()
	sum_w = 0
	for i in range(0,endd-start):
		no_weight = no_weight + all_data[i]
		sigma_2 = mean(all_data[i]**2)
		f_this = new_f[signal_pos[i]]
		w = func(f_this, *popt) / sigma_2
		sum_w += w
		weighted = weighted + all_data[i]*w

	
	sub = mean(weighted/sum_w) - mean(no_weight/(endd-start))
	if (plotting):
		plot(new_f,no_weight/(endd-start))
		plot(new_f,weighted/sum_w-sub)
		show()
	return argmax(weighted/sum_w-sub) == signal_pos[arg_signal_pos]

def animation():
	user = getpass.getuser()
	floder = f"C:\\Users\\{user}\\OneDrive - cc.ncu.edu.tw\\研究\\code\\analyse\\01_11\\picture"
	temp = argmax(data[9])
	signal_pos = []
	for i in range(4,15):
		signal_pos.append(argmax(data[i]))
	signal_pos[-2] = (signal_pos[-3]+signal_pos[-1])//2
	
	no_weight = zeros(len(data[0]))
	weighted = zeros(len(data[0]))
	POPP = see_S21()
	sum_w = 0
	for i in range(4,15):
		plot(f,shiftting(data[i]-2*(i-1),temp,signal_pos[i-4]))
		shifted = shiftting(data[i],temp,signal_pos[i-4])
		no_weight = no_weight + shifted
		sigma_2 = mean(shifted**2)
		f_this = f[signal_pos[i-4]]
		w = func(f_this, *POPP) / sigma_2
		sum_w += w
		weighted = weighted + shifted*w
		# pause(0.2)
		# savefig(f"{floder}\\{i}.png")
		# cla()

	show()
	figure()
	plot(f,no_weight/12,label = "no_weight")
	sub = mean(weighted/sum_w) - mean(no_weight/12)
	plot(f,weighted/sum_w-sub,label = "weighted")
	legend()
	xlabel("frequency[Hz]")
	ylabel("power[dbw]")
	show()
	# close()

def shiftting(y,max_pos,max_Y_pos=None):
	if (max_Y_pos==None):
		max_Y_pos = argmax(y)
	print(max_Y_pos,max_pos)
	output = np.zeros(len(y)) + mean(y)
	if (max_Y_pos < max_pos):
		output[(max_pos-max_Y_pos):] = y[:(max_Y_pos-max_pos)]
	elif (max_Y_pos > max_pos):
		output[:-(max_Y_pos-max_pos)] = y[(max_Y_pos-max_pos):]
	else:
		output = y
	return output

def combining():
	no_weight = zeros(len(data[1]))
	weight = zeros(len(data[1]))
	for i in range(1,18):
		no_weight = no_weight + data[i]
		plot(f,data[i])
	no_weight = no_weight / 18
	plot(f,no_weight)
	show()
	close()

def func(x,a,c):
    return 1/(1+ 4*(a*(x/c-1))**2)


def see_S21(plotting = True):
	f_max = []
	y_max = []
	signal_pos=[]
	for i in range(4,15):
		signal_pos.append(argmax(data[i]))
	signal_pos[-2] = (signal_pos[-3]+signal_pos[-1])//2

	for i in range(4-4,15-4):
		f_max.append(f[signal_pos[i]])
		if (i==5):
			y_max.append(data[i+4][signal_pos[i]]-0.4)
		else:
			y_max.append(data[i+4][signal_pos[i]])
		if (plotting):
			plot(f,data[i+4])

	
	f_max,y_max = array(f_max),array(y_max)
	bais = - max(y_max) + 1
	y_max = y_max + bais
	if (plotting):
		plot(f_max,y_max-bais,"p--")
	b_pos = argmax(y_max)
	b = y_max[b_pos]
	c = f_max[b_pos]
	a = 4000
	param_bounds=([1000,c-10],[20000,c+10])
	popt, pcov = curve_fit(func, f_max, y_max,bounds = param_bounds)
	if (plotting):
		plot(f,func(f,*popt)-bais)
		title(f"chisqure = {chisquare(y_max,func(f_max,*popt))[1]} Q={popt[0]}")
		print(a)
		print(f_max)
		xlabel("frequency[Hz]")
		ylabel("power[dbw]")
		show()
	return popt



if __name__=="__main__":
	show_bis()
	# print(smaller_signal(1,5))
	# total = 10
	# prob = []
	# snr_X = linspace(1,5,10)
	# for i in snr_X:
	# 	count = 0
	# 	print(i)
	# 	for _ in range(total):
	# 		if (smaller_signal(False,i)):
	# 			count += 1
	# 	prob.append(count/total)
	# plot(snr_X,prob)
	# show()
	# show_bis()
	# animation()
	# see_S21()
# combining()


