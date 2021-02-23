from matplotlib.pyplot import *
import matplotlib.mlab as mlab
from numpy import *
from scipy.optimize import curve_fit
from scipy.special import factorial
from scipy.stats import norm, chisquare, poisson, chi2

class Noise:
	def __init__(self,mu,sigma,length):
		self.noise = np.random.normal(mu,sigma,length)
		self.fitting()

	def fitting(self):
		(self.mu, self.sigma) = norm.fit(self.noise)

	def cut(self,value):
		self.keep  = np.where(self.noise >=(self.mu + value))
		self.aband = np.where(self.noise < (self.mu + value))

	def rescale(self):
		self.noise = exp(self.mu * self.noise - (self.mu)**2/2 )

	def __mul__(self,T1):
		return self.noise * T1

# do one time search, return two booleen
# Indicates whether the signal was found
# The first is the THRESHOLD
# The second is the BPM
# Total_time: How many time we scan
# power : the power we insert
# PLOTTING : Whether to draw it
def compare(Total_time,power=1.645,PLOTTING=True):
	N = 10000
	x = np.arange(10000)
	BMP = np.zeros(N) + 1
	THRESH = np.arange(10000)
	ABANDEND = []
	if (PLOTTING):
		figure()
	for i in range(Total_time):
		n1 = Noise(1,1,N)
		n1.noise[N//2] += power

		# Thresh
		n1.cut(1.645 * n1.mu)
		THRESH = np.intersect1d(THRESH,n1.keep)
		ABANDEND = np.union1d(ABANDEND,n1.aband)
	
		this_noise = copy(n1.noise)
		for each in ABANDEND:
			this_noise[int(each)] = 0
		# print(ABANDEND)
		if (PLOTTING):
			subplot(Total_time,2,2*i+1)
			plot(x,this_noise)
			plot(5000,this_noise[5000],"ro")
			
		# BPM method
		n1.rescale()
		if (PLOTTING):
			subplot(Total_time,2,2*i+2)
			plot(x,n1.noise)
			plot(5000,n1.noise[5000],"ro")
		BMP = n1 * BMP



	# figure()
	# plot(x[THRESH],n1.noise[THRESH])
	this_noise = copy(n1.noise)
	for each in ABANDEND:
		this_noise[int(each)] = 0
	if (PLOTTING):
		figure()
		subplot(121)
		plot(x,this_noise)
		plot(5000,this_noise[5000],"ro")
		subplot(122)
		plot(x,BMP)
		plot(5000,BMP[5000],"ro")
		show()

	return argmax(BMP)==N//2 , argmax(this_noise)==N//2

# insert power from 1 to 5
# in a given scan time
# then calculate the probability and save in a npy file
def gothrough(Total_time):
	how_many_times = 200
	power_slices = 25
	power_arr = np.linspace(1,5,power_slices)
	BMP_arr = []
	THR_arr = []
	for each_power in power_arr:
		print(each_power)
		BMP_t = 0
		THR_t = 0
		for _ in range(how_many_times):
			this_b, this_t = compare(Total_time, each_power, False)
			BMP_t += this_b
			THR_t += this_t
		BMP_arr.append(BMP_t/how_many_times)
		THR_arr.append(THR_t/how_many_times)
	np.save(f"{Total_time}.npy",[power_arr,BMP_arr,THR_arr])

# see the result from gothrough()
def plot_go_through(Total_time):
	x, BMP, THR = np.load(f"{Total_time}.npy")
	title(f"scan {Total_time} time")
	plot(x,BMP,label="BMP")
	plot(x,THR,label="THRESHOLD")
	xlabel(r"power [$\sigma$]")
	legend()
	show()

def analy_TH():
	for Total_time in range(1,6):
		x, BMP, THR = np.load(f"{Total_time}.npy")
		title(f"THRESHOLD")
		plot(x,THR,label=f"{Total_time}")
	xlabel(r"power [$\sigma$]")
	legend()
	show()

def analy_TBM():
	for Total_time in range(1,6):
		x, BMP, THR = np.load(f"{Total_time}.npy")
		title(f"BMP")
		plot(x,BMP,label=f"{Total_time}")
	xlabel(r"power [$\sigma$]")
	legend()
	show()

if __name__=="__main__":
	# print(compare(3,2,True))  
	# for i in range(1,5):
	# 	gothrough(i)
	# gothrough(5)
	# gothrough(1)
	# plot_go_through(2)
	# for i in [3,4,5,10]:
	# 	plot_go_through(i)
	# analy_TBM()