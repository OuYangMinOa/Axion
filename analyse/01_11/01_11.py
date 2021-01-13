import pandas as pd
from matplotlib.pyplot import *
from numpy import *
from scipy.optimize import curve_fit
import getpass

file_name = "01_11.txt"
data = pd.read_csv(file_name, header = None)
print(data)
data = data.T


f = data[0]
y = data[1]

def animation():
	user = getpass.getuser()
	floder = f"C:\\Users\\{user}\\OneDrive - cc.ncu.edu.tw\\研究\\code\\analyse\\01_11\\picture"
	temp = argmax(data[9])
	signal_pos = []
	for i in range(4,15):
		signal_pos.append(argmax(data[i]))
	signal_pos[-2] = (signal_pos[-3]+signal_pos[-1])//2
	
	no_weight = zeros(len(data[0]))

	for i in range(4,15):
		plot(f,shiftting(data[i]-2*(i-1),temp,signal_pos[i-4]))
		no_weight = no_weight + shiftting(data[i],temp,signal_pos[i-4])
		# pause(0.2)
		# savefig(f"{floder}\\{i}.png")
		# cla()
	show()
	figure()
	plot(f,no_weight/12)
	show()
	close()

def shiftting(y,max_pos,max_Y_pos=None):
	if (max_Y_pos==None):
		max_Y_pos = argmax(y)
	print(max_Y_pos,max_pos)
	output = [nan for i in range(len(y))]#np.zeros(len(y)) + mean(y)
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

def see_S21():
	f_max = []
	y_max = []
	for i in range(4,13):
		print(f[argmax(data[i])],max(data[i]))
		f_max.append(f[argmax(data[i])])
		y_max.append(max(data[i]))
		plot(f,data[i])
	plot(f_max,y_max)
	for i in range(1,len(f_max)-1):
		print(f_max[i+1]-f_max[i])
	print(f_max)
	show()
	figure()
	plot(f_max)
	show()

if __name__=="__main__":
	animation()
	# see_S21()

# combining()


