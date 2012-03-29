"""
1)Improvements in the constant ADC_BP
2)in the plotpeak function chek the stuff properly
3)The corner frequencies in Vivien Thesis are .25 to 30
4)Check the peak point , by taking the window around the pake..or check the peak max in three samples(the answer being in midlle)
5)the presuure decrease must be linear.
"""
import numpy
from scipy import *
from scipy import signal
from matplotlib.pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.mlab as ax
from scipy.linalg import lstsq
import scipy.special as sp
import os
name = ""
iter_value = 1
tap_val = 6
left_crop = 3000
right_crop = -2000
ADC_BP = 60.27 # 57 was the previous value
ADC_BP_CONSTANT = 161
UPPER_BP_LIMIT  = ADC_BP*200
LOWER_BP_LIMIT  = ADC_BP*40
Dia_constants = [.95,.9,.85 , .8 ,.75,.7 ]
Sys_constants = [.35 , .4 ,.45,.475,.5,.55 ]
def main():

	A = loadtxt('testdata.dat') # loading the data from the file
	print "Started \n"
	starttest(A )
	print "\nEnd"

def starttest(A ):
	global left_crop
	global right_crop
	#Calculating the max value and max index
	Amax ,indx = max2(A)
	 # filtering out data from the file , first data of raising and some data at last , they are zeros	
	A = A[indx:]
	#filtering out data that has ADC value of greater than 3000 that is approximately 40 mm Hg
	A=numpy.array(A)
	A=A[find(A>LOWER_BP_LIMIT)]
	#filtering out data that has ADC value of lesser than 11000 that is approximately 200 mm Hg
	A=A[find(A<UPPER_BP_LIMIT)]
#	print "No pf samples having values between upper limit and l;ower limit are " , len(A)
	#Calculating moving avg to remove noise and remove unwnated peaks
	#Calculating moving average of window_size =50
	A = calcMovingAvg(A , 50)
#	plot_fft(A)
	heart_beat = difference_equation(A)	
	#Croping the values

	A = A[left_crop:right_crop]
	A = A + ADC_BP_CONSTANT*ones(len(A))

	avg_diff = calPeakValley(A) # gives the average diffrence between two heart beats ;calculates both peak diff and valley diff

	minmax_data_index , minmax_data , lenmm = removePeak(heart_beat , avg_diff /2)

	negative_data_index = 	minmax_data_index[1::2]
	negative_data = 	minmax_data[1::2]	
	positive_data_index = 	minmax_data_index[::2]
	positive_data = 	minmax_data[0::2]	

	shift_required = (sum(positive_data) + sum(negative_data)) / ( len(negative_data) +len(positive_data))
	positive_data = positive_data - shift_required*ones(len(positive_data))
	negative_data = negative_data - shift_required*ones(len(negative_data))

#	Converting the values to DC
	heart_beat = heart_beat - shift_required*ones(len(heart_beat))
	print "DC shift required " ,shift_required
	avgfile = open("testdatafrompythonheart.dat" , 'w')
	for i in heart_beat:
		avgfile.write(str(i)+" \n")	
		
	avgfile = open("testdatafrompython.dat" , 'w')
	for i in A:
		avgfile.write(str(i)+" \n")			
	#Here the positive_envelope is is the max between the postive data given
	#similar things to do in negative envelope side
	positive_envelope = zeros(len(positive_data))
	for i in range( 1 , len(positive_data_index)):
		pos1 =  positive_data_index[i-1][0]
		pos2 =  positive_data_index[i][0]
		arr_con = heart_beat[pos1:pos2] 
		positive_envelope[i] = max(arr_con ) #- min(arr_con)
	positive_envelope = calcMovingAvg(positive_envelope , tap_val)#Removing the noise from the positive envelope
#	calcTap(positive_envelope , 7)

	maxval , maxid = max2(heart_beat)
	figure()
	#ADC_BP is a number obtained by calculations
	#150 is just a scaling factor used
	plot(A/ADC_BP)
	plot(heart_beat*150/ADC_BP)
	plot(positive_data_index , positive_envelope*150/ADC_BP ) #forms the evelope required
	length_A = len(A)

#################################################################################CODE DONE TILL HERE##################
	#code checked till here
	plotpeak(positive_data_index ,positive_envelope*150/ADC_BP ,length_A ,A)

def calPeakValley(A):

	peak = []
	valley = []
	peak_diff=[]
	valley_diff=[]
	deviation = 15
	i= deviation
	peak_value = []
	valley_value = []
	while True:
		#For peaks
		if(i > len(A) - deviation -1):
			break
		if   (sum(A[i-deviation:i+deviation]) > sum(A[i+deviation:i+3*deviation])) and ( sum(A[i-deviation:i+deviation]) > sum(A[i-3*deviation:i-1*deviation])):
			peak.append([ i] )
			peak_value.append( A[i] )
			#print "PEAK " ,A[i] , " at" ,  i
			i = i +deviation
		elif   (sum(A[i-deviation:i+deviation]) < sum(A[i+deviation:i+3*deviation])) and ( sum(A[i-deviation:i+deviation]) < sum(A[i-3*deviation:i-1*deviation])):
			valley.append([ i] )
			i = i +deviation			
			valley_value.append([A[i]])
			#print "VALLEY ", A[i] , " at" ,  i			
		i +=1

	j=0
	for i in peak:
		diff = i[0] - j	
		peak_diff .append([diff])
		j=i[0]
	j=0
	for i in valley:
		diff = i[0] - j
		valley_diff .append([diff])		
		j=i[0]	
		
	peak_diff = numpy.array(peak_diff)
	valley_diff = numpy.array(valley_diff)
	peak_value = numpy.array(peak_value)
	valley_value = numpy.array(valley_value)
	#Used 100 by experiment	
	peak_diff=peak_diff[ find(peak_diff >100)]
	valley_diff=valley_diff[ find(valley_diff >100)]
	avg_diff = sum(valley_diff)/len(valley_diff) + sum(peak_diff)/len(peak_diff)
	avg_diff /=2
	print "Heart Beat Rate = " , 480*60/avg_diff
	return avg_diff
	

#This function removes the erraneous peaks
#Algorithm:yet to complete
"""

"""
def removePeak(A , allowed_gap):
	deviation = 20

	min_max = 1 # 1 -to find maxima ... 0 -to find minima
	minmax_data = []	 #starts with positive 
	minmax_data_index = []	
	prev_min_data = 0
	prev_max_data = 0
	first_positive = 0
	A = numpy.array(A)
	length = len(A)
	for i in range( deviation , length - deviation -1):
		#Maxima
		if   (sum(A[i-deviation:i+deviation]) > sum(A[i+deviation:i+3*deviation])) and ( sum(A[i-deviation:i+deviation]) > sum(A[i-3*deviation:i-1*deviation]) ):
			data =  max2(A[  i-deviation:i+deviation  ] , offset = i-deviation) #gives the max of the data with a offset
			if  data[0] < 0:#Bad peak
				continue
			if min_max == 1 :
				minmax_data.append(  [data[0]])
				minmax_data_index.append(  [data[1]])
				min_max = 0		
				first_positive = 1#So that we can start from min data also
			if prev_max_data < data[0] and min_max == 0: #the previous peak was a fase peak.so a new peak is replaces the old peak
				minmax_data[-1][0] = data[0]
				minmax_data_index[-1][0] = data[1]
				first_positive = 1#So that we can start from min data also
			prev_max_data = data[0]			
		#Minima
		elif   (sum(A[i-deviation:i+deviation]) < sum(A[i+deviation:i+3*deviation])) and ( sum(A[i-deviation:i+deviation]) < sum(A[i-3*deviation:i-1*deviation])):
			data = min2(  A[i-deviation:i+deviation],
						 offset = i-deviation)
			if  data[0] > 0:#Bad valley
				continue			
			if min_max ==0:
				minmax_data.append(  [data[0]])
				minmax_data_index.append(  [data[1]])
				min_max =1
		
			elif prev_min_data < data[0] and min_max ==1 and first_positive == 1: #first positive is to account for first negative data
				minmax_data[-1][0] = data[0]
				minmax_data_index[-1][0] = data[1]
				
			prev_min_data = data[0]
	minmax_data = numpy.array(minmax_data)	
	return minmax_data_index , minmax_data ,len(minmax_data)
	
def plotpeak(data_index ,positive_envelope,length_A ,A):
	global name
	global iter_value
	global ADC_BP
	global Sys_constants
	global Dia_constants
	maxpeak , maxpeak_index = max2(positive_envelope)
	
	systolic  = [ ]
	diastolic  = []
	n = 0
	print "MAP is at pressure = " , A[data_index[maxpeak_index]][0]/ADC_BP
	for i in arange(0 , 1 , 0.005):
		n +=1
		positive_envelope = numpy.array(positive_envelope)
		great_array = find(positive_envelope > maxpeak*i)
		#Chance for improvement here
		#For calculating systolic
		left_peak_index = data_index[(great_array[0])] 		
#		print left_peak_index , "   at "	,i , " adc value is " ,A[left_peak_index]
		sys_arr = A[left_peak_index] /ADC_BP
		systolic.append([sys_arr[0] , i])
		
		#For calculating Diastolic
		right_peak_index = data_index[(great_array[-1])] 					
		dia_arr = A[right_peak_index] /ADC_BP
		diastolic.append([dia_arr[0] , i])
	
		#Plot only 10 data
		if(n%20 == 10):
			plotarr = ones(length_A)*maxpeak*i
			plotarr[left_peak_index[0]:right_peak_index[0]] =0
			plot(range(length_A) , plotarr)
			
	for values in Sys_constants:
		print "Systolic Pressure for " , values, " is " , systolic[int(values/.005)][0]
	print "\n"
	for values in Dia_constants:
		print "Diastolic Pressure for " , values, " is " , diastolic[int(values/.005)][0]
	range1  = arange(0 ,1 ,.005)
	title("The value of ADC to BP is " +str(ADC_BP))
	figure()
	plot(range1 , systolic  )
	plot(range1 , diastolic )
	title("Dia-sys :The value of ADC to BP is " +str(ADC_BP))

#Calculates the moving WINDOW average of a particular array.The initial and final values are left out as it is
def calcMovingAvg(A , win_size):
	A_mavg=zeros(len(A))		
	#limit will be rounded off to nearest integer
	limit  =  win_size / 2 
	if( win_size % 2 == 0):
		for i , value in enumerate(A[limit:-limit]):
			j=i+limit
			A_mavg[j] = sum(A[j-limit: j+limit]) /win_size
		
		A_mavg[0:limit] = A[0:limit]
		A_mavg[-limit:] = A[-limit:]
	else:

		for i , value in enumerate(A[(1+limit):-limit]):
			j=i+limit
			A_mavg[j] = sum(A[j-(1+limit): j+limit]) /win_size
		
		A_mavg[0:(limit +1 )] = A[0:(limit+1)]
		A_mavg[-limit:] = A[-limit:]	
	return A_mavg

#Calculates the moving TAP average of a particular array.The initial and final values are left out as it is
def calcTap(A ,win_size = 6):
	arr_size = len(A)
	tapArray = zeros(arr_size)
	limit       =	 win_size /2
	if(win_size % 2 == 0): #Even win_size
		for i in range(limit,arr_size-limit):
			sum1 = 0
			baseSum = 0
			for j in range(1,limit+1):
				sum1 += j*A[i-limit +j]
				baseSum +=j+1
			for j in range(1,limit+1 ):
				sum1 += (limit+1-j)*A[i +j]
				baseSum +=j+1
			tapArray[i] = sum1/baseSum
	else:
		for i in range(limit+1,arr_size-limit):
				sum1 = 0
				baseSum = 0
				for j in range(0,limit):
					sum1 += (j+1)*A[i-limit +j]
					baseSum +=j+1
				for j in range(0,limit+1 ):
					sum1 += (limit+1-j)*A[i +j]
					baseSum +=j+1
				tapArray[i] = sum1/baseSum
	tapArray[:limit+win_size%2] = A[:limit+win_size%2]
	tapArray[-limit:] = A[-limit:]		
	subplot(2,1,1)
	plot(A)
	subplot(2,1,2)
	plot(tapArray)
	show()


def max2(A , offset = 0):
	Amax = max(A)
	indx = 0
	for i in xrange(len(A)):
		if(A[i] == Amax):
			indx = i
			break

	return [Amax , indx+offset]

def min2(A , offset = 0):
	Amin = min(A)
	indx = 0
	for i in xrange(len(A)):
		if(A[i] == Amin):
			indx = i
			break
	return [Amin , indx + offset]
	
def avg2(A , offset=0):
	A = numpy.array(A)
	avg = sum(abs(A))/len(A)
	return avg+offset
	
#This function gives the value from the diffrence equation of A matrix
def difference_equation(array):
	global left_crop
	global right_crop
	A = [1.000000000000000e+000, -7.354587627363298e+000, 2.368737561318704e+001, -4.364095783862433e+001, 
   	  5.030800039631103e+001, -3.715931887201098e+001, 1.717523899762785e+001, -4.541874057040570e+000, 
	     5.261233879178731e-001]

	B = [1.667728318461477e-004, 0, -6.670913273845909e-004, 0, 
   	  1.000636991076886e-003, 0, -6.670913273845909e-004, 0, 
   	  1.667728318461477e-004]
   	
#   	B, A = signal.butter(4,[.5/470,40/470])
   	filtered_array = signal.lfilter(B,A,array)
   	return filtered_array  [left_crop:right_crop]
