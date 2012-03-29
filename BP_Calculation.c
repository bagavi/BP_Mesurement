/*
READ - The name of the file taken as input is "testdata.dat" - in line 587.Change that accordingly to make code work
Written by Vivek Kumar Bagaria <vivekee047@gmail.com> , 9176079646 ,IIT Madras
*/
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

/*Header File*/
#include "bp_header.h"

//These value have yet not been used
float SystolicPressure, DiastolicPressure, PulseRate ;


/*
	The below function is a filter which has the co-efficients defined and will give the output using the co-efficients
	This is used to extract heart beat from the signal from the pressure sensor.
	

	Argument descrioption
	FilterInputArray - Input array
	FilterOutputArray - Outputput array
	length = length of both the above arrays
*/
void Filter (float *FilterInputArray, float *FilterOutputArray, int length)
{
	int i ; // Loop Variables
	int n ; // Loop Variable
	double Y[9] ; // Used to calulate values of output in difference equation
	float *X ;
	/*
		The below co-effcients are calculated from matlab
		They are butterWorth Filter co-effcients
			Corner frequency	=>.5/480 and 30/480
			Order 		 	=> 4
	*/
	double OutPutCoefficients[9]= {1.000000000000000e+000, -7.354587627363298e+000, 2.368737561318704e+001, -4.364095783862433e+001, 5.030800039631103e+001, -3.715931887201098e+001, 1.717523899762785e+001, -4.541874057040570e+000, 5.261233879178731e-001} ;
	
	double InPutCoefficients[9]= {1.667728318461477e-004, 0, -6.670913273845909e-004, 0,1.000636991076886e-003, 0, -6.670913273845909e-004, 0, 1.667728318461477e-004} ;
	X = (float*) malloc(sizeof(float)*(length+8)) ;

	if(X == NULL)
		return ;
	
	//Basic initialization
	for(i = 0 ; i < 9 ; i++)
	{
		Y[i] =0 ;
		X[i] =0 ;
		FilterOutputArray[i] = 0 ;
	}

	for(n = 8 ; n <= length + 7 ; n++)
	{
		/*
			First part of difference equation.it is out of loop just for code convinience
		*/
		X[n] = FilterInputArray[n - 8] ;
		Y[8] = InPutCoefficients[0] * (double) X[n] ;
			
		//Loop for each value of output
		for(i = 1 ; i <= 8 ; i++)
			{
				/*
					Applying the difference equation
				*/
				Y[8] += ( (InPutCoefficients[i] * ((double) X[n - i]) ) - ((OutPutCoefficients[i] * (double) Y[8 - i])) ) ;
			}

		FilterOutputArray[n] = (float) Y[8] ;
		
		/*
			This below part is a bit tricky in logic(easy in terms of code) 
		*/
		for(i = 0 ; i < 8 ; i++)
			Y[i] = Y[i+1] ;
	}
	//Free arrays
	free(X) ;
}



/*
	The below function calcuates the moving avg of an array.Assuming window size is lesser than array length
	Assuming that thw windowSize is even

	ArgumentDescription
	ArrayLength = length of the arrays to be operated
	windowSize = Size of the avg window
	Array = InputArray
	WindowArray = Result of windowing of InputArray
*/
void ArrayMovingAverage(int ArrayLength, int windowSize, float *Array, float *WindowArray)
	{
		int i , j ;//Loop Variables
		int sum =0 ;
		int sumlength = 0 ;
		for(i= 0 ; i< windowSize/2 ; i++)
		{
			sum = 0 ;
			sumlength = 0 ;
			for(j = 0 ; j<i+windowSize/2 ; j++ , sumlength++)
			{
				sum += Array[j] ;
			}
			WindowArray[i] = sum/sumlength ;
		}
		
	for( ; i<ArrayLength-windowSize/2 ; i++)
		{
			sum = 0 ;
			sumlength = 0 ;
			for(j = i-windowSize/2 ; j<i+windowSize/2 ; j++ , sumlength++)
			{
				sum += Array[j] ;
			}
			WindowArray[i] = sum/sumlength ;
		}
		
	for ( ; i<ArrayLength ; i++)
		{
			sum = 0 ;
			sumlength = 0 ;
			for(j = i-windowSize/2+1 ; j<ArrayLength ; j++ , sumlength++)
			{
				sum += Array[j] ;
			}
			WindowArray[i] = sum/sumlength ;
		}
		
	}



/*
	Add constant value to each element in the array
	
	Argument Description
	Array - Input Array
	Length - length of the input array
	constant - the constant value to be added to every element
	
*/
void Add_Constant(float *Array, int length, float constant)
	{
		int i ; //Loop Variable
		for( i = 0 ; i < length ; i++)
			{
				Array[i] += constant ;
			}	
	}

	
	
/*
	Initializes all the values to zero before the given index
	
	Argument Description
	Array - Input Array
	index = the Index till which all the elements of the array will be Zero
*/
void IndexCutArray(float *Array, int index)
	{
		int i ; // Loop Variable
		for(i = 0 ; i<index ; i++)
		{
			Array[i] = 0;
		}
	}
	
	
	
/*
	Cuts the array from left till it reaches a value which is lesser than the argument 2.
	Removes all the ADC input above the UPPER BP LIMIT
	
	Argument Description
	ArrayLength = length of the given array
	Array = Input array
*/
void CutArrayLeft( int ArrayLength, float *Array )
	{
	float value =UPPER_BP_LIMIT ;
	int i ;
		for(i = 0 ; i< ArrayLength ; i++)
		{
			if(( Array[i] > value ) || (Array[i] == 0))
			{
				Array[i] = 0 ;
			}
			else
			{
				break ;
			}
		}
	}		



/*
	Cuts the array from right till it reaches a value which is greater than the argument 2
	Removes all the ADC input below the LOWER BP LIMIT
	
	Argument Description
	ArrayLength = length of the given array
	Array = Input array
*/
void CutArrayRight(int ArrayLength , float *Array)
	{
	float value =LOWER_BP_LIMIT ;
	int i ;
		for(i = ArrayLength ; i>0 ; i--)
		{
			if( (Array[i] < value) ||(Array[i] == 0))
			{
				Array[i] = 0
				 ;
			}
			else
			{
				break ;
			}
		}
	}

	


/*
	Below function returns the maximum value of an array
	
	
	Argument Description
	Array = Input array
	length = Array length
*/
int MaxOfArray(float *Array, int length)
	{
		int i ; //loop variable
		float MaxVal = 0 ;
		for(i = 0 ; i< length ; i++)
		{
			if(Array[i] > MaxVal)
			{
				MaxVal = Array[i] ;
			}	
		}
		return MaxVal ;
	}




/*
	Returns the index of the array where the maximum value occurs
	
	Argument Description
	Array = Input array
	start = the index to start considering for max value
	length = Array length
*/
int IndexOfMaxArray(float *Array, int start, int length)
	{
		int i ; //loop variable
		float MaxVal = Array[start] ;
		int IndexForMax = start ;
		for(i = start ; i< length ; i++)
		{
			if(Array[i] > MaxVal)
			{
				MaxVal = Array[i] ;
				IndexForMax = i ;
			}
		}
		return IndexForMax ;
	}




/*
	Below functon gives the index where the maximum value occurs
	
	
	Argument Description
	Array = Input array
	start = the index to start considering for max value
	length = Array length
*/
int IndexOfMinArray(float *Array, int start, int length)
	{
		int i ; //loop variable
		float MinVal = Array[start] ;
		int IndexForMin = start ;
		for(i = start ; i< length ; i++)
		{
			if(Array[i] > MinVal)
			{
				MinVal = Array[i] ;
				IndexForMin = i ;
			}
		
		}
		return IndexForMin ;
	}




/*
	Count number of non-zero samples in an array and returns that number (the zeros are counted only at the extremes)
	
	
	Argument Description
	Array = Input array
	length = Array length
*/
int CountNonZeroNumber(float *Array, int length)
{
	int initialZeros = 0 , finalZeros = 0 ;
	int i ; //loop variable
	for ( i =0 ; i<length ; i++)
	{
		if(Array[i]!=0)
		{
			initialZeros = i ;
			break ;
		}
	}

	for ( i=length ; i>0 ; i--)
	{
		if(Array[i]!=0)
		{
			finalZeros = i ;
			break ;
		}
	}
	return finalZeros-initialZeros+1 ;
}


/*
	The below functions removes some values from the initial and final part of array
	startNum - number of initial values to be made zero
	endNum - number of final values to be made zero	
	*Array - The array to which this operaton will be operated
	
	Assumption startNum and endNum are lesser than length of Array
*/

void RemoveExtremeValues(int startNum, int endNum, int length, float *Array)
	{	
		int i = 0 ;
		for(i=0 ; i<startNum ; i++)
		{
			Array[i] = 0 ;
		}
		
		for(i=0 ; i<endNum ; i++)
		{
			Array[length-(i+1)] = 0 ;
		}
	}

/*

	This function is a bit complicated.
	
	
	Argument Description
	Initial - Left croping of the input array
	final - Right croping of the input array
	Array - InputArray
	HeartBeatPeak - Stores the index of the points where a peak occurs
	HeartValleyPeak - Stores the index of the points where a valley occurs
	HeartBeatPeakValue - Stores the values of the peak
	HeartBeatValleyValue - Stores the values of the valley
	
	This function calcuates the the avg diff between two beats and thus by knowing the sampling rate we can find out the heart beat rate

*/
void CalcHeartBeatRate(int initial, int final, float *Array, int *HeartBeatPeak, int *HeartBeatValley, float *HeartBeatPeakValue, float *HeartBeatValleyValue, int *NoPeaks, int *NoValleys)
	{	
		int i ,j ; //Loop Variables
		int peakNumber =0, valleyNumber =0 ;
		int deviation = 8 ;						// determines the size of the sample used to compare ,to find greatest and least values
		float FirstSum , MiddleSum , LastSum ;
		int index ;							// Stores the value of the index corresponding to peaks and valleys
		int minmax = 1 ; 						//Takes value 0 and 1 only , it it the variable which determines wheter a peak or valley must come
											// 1 stands for maxima
		int positive_first = 0 ;					 // This ensures that valley can be registered before peak

		for (i = 3*deviation + initial; i < final - 3*deviation ; i++)
			{
				FirstSum 	= 0 ;
				MiddleSum 	= 0 ;
				LastSum 		= 0 ;
				//Intial 2*deviation samples
				for(j = i - 3*deviation ; j < i - deviation ; j++)
				{
					FirstSum += Array[j] ;
				}
				//Middle 2*deviation samples
				for(j = i - deviation ; j < i + deviation ; j++)
				{
					MiddleSum += Array[j] ;
				}
				//Last 2*deviation samples
				for(j = i + deviation ; j < i + 3*deviation ;j++)
				{
					LastSum += Array[j] ;
				}
				
				/*
					Maxima(peak)
				*/
				if( (MiddleSum > FirstSum) && (MiddleSum > LastSum) )
				{

					index = IndexOfMaxArray( Array , i - deviation , i + deviation) ;	
					if(Array[index] < 4) 	// Bad peak..cannot be a peak, i have choosen 4 after seeing many graphs
										//This is just so that we ignore all the peaks with very lees amplitude
					{
						continue ;
					}
					
					else if(minmax == 1) // Expecting a peak
					{
						HeartBeatPeak[peakNumber] = index ;
						HeartBeatPeakValue[peakNumber] = Array[index] ;
						peakNumber++ ;
						i += deviation ;
						minmax = 0 ;
					}
					
					//The below case is that we are expecting a valley but we found another peak ,so this peak must be a false peak or the previous peak must be a false peak
					else 
					{

						if(Array[index] > HeartBeatPeakValue[peakNumber - 1] )
						{
							HeartBeatPeak[peakNumber-1] = index ;
							HeartBeatPeakValue[peakNumber-1] = Array[index] ;
							i += 1*deviation ;
						}
					}
				}
				
				/*
					Minima
				*/
				else if( (MiddleSum < FirstSum) && (MiddleSum < LastSum) )
				{
				
					index = IndexOfMinArray( Array , i - deviation , i + deviation) ;
					
					if(Array[index] > 0) //Bad valley...cannot be a valley
					{
						continue ;
					}

					else if(minmax == 0 || positive_first ==0) // Expecting a valley
					{
						positive_first = 1 ;
						HeartBeatValley[ valleyNumber ] = index ; 
						HeartBeatValleyValue[ valleyNumber ] = Array[index] ; 
						valleyNumber++ ;
						i += 2*deviation ;
						minmax = 1 ;
					}
					
					else if(Array[index] < HeartBeatValleyValue[ valleyNumber - 1 ] )
					{
						HeartBeatValley[ valleyNumber - 1 ] = index ; 
						HeartBeatValleyValue[ valleyNumber -1 ] = Array[index] ; 
						i += 2*deviation ;
					}
				}
			}
			
			*NoPeaks = peakNumber ;
			*NoValleys = valleyNumber ;
	}




/*
	Function which calculates Systolic and Diastolic pressure for various Constants
	
*/

	void SysDiaCalc(int length, float *HeartEnvelope, int *HeartEnvelopeIndex, float *DC_Array)
	{
		float MAP ;
		int MAP_Index = 0 ;
		float MapSys = 0 ;
		float MapDia = 0 ;
		int i , j ;//Loop Variables
		MAP = MaxOfArray(HeartEnvelope ,length) ;
		MAP_Index = IndexOfMaxArray(HeartEnvelope, 0, length) ;
		printf("MAP Value %f\n", DC_Array[ HeartEnvelopeIndex[MAP_Index] ] / ADC_BP) ;
		float Dia_constants[12]	= {.95,.9,.85 , .8 ,.75,.7,.65,.6 ,.55,.5 ,.45,.4} ;
		float Sys_constants[8]	= {.25,.3,.35 , .4 ,.45,.475,.5,.55 } ;
		
		//Systolic pressure
		for(i = 0 ; i< 8 ; i++)
		{
			MapSys = Sys_constants[i]*MAP ;
			for( j = 0 ; j < MAP_Index ; j++)
			{
				if(MapSys <= HeartEnvelope[j])
				{
					printf("Systolic : Constant %f value is %f\n" , Sys_constants[i] , DC_Array[ HeartEnvelopeIndex[j] ] / ADC_BP) ;
					break ;
				}
			}
		}

		printf("\n\n") ;	
		//Distolic pressure
		for(i = 0 ; i< 12 ; i++)
		{
			MapDia = Dia_constants[i]*MAP ;
			for( j = length -1 ; j > MAP_Index ; j--)
			{
				if ( MapDia <= HeartEnvelope[j] )
				{	
					printf("Diastolic : Constant %f value is %f\n" , Dia_constants[i] , DC_Array[ HeartEnvelopeIndex[j] ] / ADC_BP) ;
					break ;
				}
			}
		}
	}




/*
	The main function 
*/
void main()
{
	
	FILE *InputFile, *OutputFile1, *OutputFile2 , *OutputFile3, *OutputFile4 ; //File varibles used to write into files for reference
	float *AC_Data, *DC_Data, *Trunc_DC_Data, *Avg_DC_Data ;
	long NoOfSamples = 0, NumberNonZeroSamplesAfterTrunc ;
	int i,j ; //Loop variable
	int li ; //test loop variable
	int MaxPressureValue = 0 , IndexOfMaxValue = 0; //The below two variables store the maxvalue and the corresponding index for that max value
	
	int HeartBeatPeak[100] ;			// Stores the index of all the peaks in the heatrbeat diagram
	float HeartBeatPeakValue[100] ;		// Stores the value of the heartbeat at the peak(index i stored in the above array)
	float AvgHeartBeatPeakValue[100] ;	// Stores the value of the Averaged heartbeat at the peak(index i stored in the above array)
	int NoPeaks ;						//Stores the number of peaks
	float PeakValueSum = 0 ;
	
	int HeartBeatValley[100 ] ;			// Stores the index of all the valley in the heatrbeat diagram
	float HeartBeatValleyValue[100] ;		// Stores the value of the heartbeat at the valley(index i stored in the above array)
	float AvgHeartBeatValleyValue[100] ;	// Stores the value of the Averaged heartbeat at the valley(index i stored in the above array)
	int NoValleys ;					//Stores the number of valleys
	float ValleyValueSum = 0 ;
	
	float Shift = 0 ;

	printf("#################################################################################\n\n") ;
	//Allocate memory for AC_Data and DC_Data
	
	AC_Data = (float *) calloc(MAXCOUNT, sizeof(float)) ;
	if (AC_Data == NULL)
	{
		printf("Could not allocate memory\n") ;
		return;
	}

	DC_Data = (float *) calloc(MAXCOUNT, sizeof(float)) ;
	if (DC_Data == NULL)
	{
		printf("Could not allocate memory\n") ;
		free (DC_Data) ;
		return ;
	}

	InputFile = fopen("testdata.dat","r+") ;
	if(InputFile == NULL)
	{
		printf("File not present") ;
		return ;
	}


	/*
		Read the data from the file and store it in DC_Data
	*/
	do
	{
		fscanf(InputFile, "%f", &DC_Data[NoOfSamples]) ;
		NoOfSamples++ ;
	}while (!feof(InputFile));
	
	NoOfSamples -= 1 ;
	printf("The number of samples are %ld \n" , NoOfSamples) ;
	fclose(InputFile) ;
	
	
	/*
		Adding the constant to ADC Input
	*/
	Add_Constant(DC_Data, NoOfSamples ,ADC_BP_CONSTANT) ;
	
	
	/*
		Calculating the max pressure of the DC_Data to remove the data before the max pressure
	*/
	MaxPressureValue = MaxOfArray(DC_Data, NoOfSamples) ;
	
	
	/*
		Calculate the index of the max pressure
	*/
	IndexOfMaxValue = IndexOfMaxArray(DC_Data,0 , NoOfSamples) ;
	
	
	/*
		Removed all the values before the pressure is released(data during pump is not cut)
	*/
	IndexCutArray(DC_Data , IndexOfMaxValue) ;
	
	
	/*
		Removed all the BP values which correspond to pressure higher than 200 mm Hg
	*/
	CutArrayLeft( NoOfSamples ,DC_Data ) ;
	
	
	/*
		Removed all the BP values which correspond to pressure lower than 40 mm Hg
	*/
	CutArrayRight( NoOfSamples , DC_Data) ;
	
	
	/*
	 NumberNonZeroSamplesAfterTrunc stores the number of non-zero valuea after removing the initial data and 
	 Removing data out of the range of Bp,specified by LOWER_BP_LIMIT and UPPER_BP_LIMIT
	*/
	NumberNonZeroSamplesAfterTrunc = CountNonZeroNumber(DC_Data , NoOfSamples) ;
	
	
	/*
		Allocating space for the new data
	*/
	Trunc_DC_Data = (float *) calloc(NumberNonZeroSamplesAfterTrunc, sizeof(float)) ;
	if (Trunc_DC_Data == NULL)
	{
		printf("Could not allocate memory\n") ;
		free (Trunc_DC_Data) ;
		return ;
	}


	/*
		The below Part is to data non-zero data to trunc data,this is a very simple stuff,can be left if you are reading this code for the first time
	*/
		{
			j = 0;
			for(i=0; i<NoOfSamples;i++)
			{	
				if(DC_Data[i]!=0)
				{
					break;
				}
			}

			for( ; i<NoOfSamples ; i++)
			{
				if(DC_Data[i]!=0)
				{
			
					Trunc_DC_Data[j] = DC_Data[i] ;
					j++ ;
				}
				else
				{
					break;
				}
			}
		}



	/*
		Avg_DC_Data Store values after averaging the Trunc_DC_Data
	*/
	Avg_DC_Data = (float *) calloc(NumberNonZeroSamplesAfterTrunc, sizeof(float)) ;
	// Window Size is 30
	ArrayMovingAverage(NumberNonZeroSamplesAfterTrunc , 30 , Trunc_DC_Data , Avg_DC_Data) ; 
	
	
	/*
		AC_Data stores the data after passing through the filter
		AC_Data contains heart beats
	*/
	Filter(Avg_DC_Data , AC_Data , NumberNonZeroSamplesAfterTrunc) ;
	
	
	/*
		To remove some initial samples from AC_Data as they are noise introduced by the filter
	*/
	RemoveExtremeValues( LeftCrop, RightCrop , NumberNonZeroSamplesAfterTrunc ,AC_Data) ;
	
	
	/*
		Below we pass 5 arrays and calculate the heart beat rate
	*/
	 CalcHeartBeatRate(LeftCrop, NumberNonZeroSamplesAfterTrunc - RightCrop, AC_Data, HeartBeatPeak, HeartBeatValley, HeartBeatPeakValue,HeartBeatValleyValue, &NoPeaks, &NoValleys ) ;
	
	
	for(i = 0 ; i < NoPeaks ; i++)
	{
		PeakValueSum += HeartBeatPeakValue[i] ;
	}
	
		for(i = 0 ; i < NoValleys ; i++)
	{
		ValleyValueSum += HeartBeatValleyValue[i] ;
	}
	
	
	/*
		Calculation moving average of the lower and the higher envelope
		TapVal- is the window size
	*/
	
	ArrayMovingAverage(NoPeaks, tapval, HeartBeatPeakValue, AvgHeartBeatPeakValue) ;
	ArrayMovingAverage(NoValleys , tapval , HeartBeatValleyValue , AvgHeartBeatValleyValue) ;


	/*
		Calculating the required DC Shift in the HeartBeat array and then shifting fotnt he envelopes by thw shift
	*/
	Shift = (PeakValueSum + ValleyValueSum )/(NoPeaks+NoValleys) ; 
	printf("Shfit = %f " ,Shift) ;
	for(i = 0 ; i < NoPeaks ; i++)
	{
		HeartBeatPeakValue[i] -= Shift ;
		AvgHeartBeatPeakValue[i] -= Shift ;
	}
	
		for(i = 0 ; i < NoValleys ; i++)
	{
		HeartBeatValleyValue[i] -= Shift ;
		AvgHeartBeatValleyValue[i] -= Shift ;
	}


	/*
		Printing the heartbeat
	*/
	printf("peaks = %d valleys = %d ,number of samples = %ld \n" , NoPeaks , NoValleys , NumberNonZeroSamplesAfterTrunc) ;
	printf("heart beat = %ld \n" , ((NoPeaks + NoValleys)/2)*60*470/(NumberNonZeroSamplesAfterTrunc) ) ;
	
	
	
	/*
		Calculating Systolic and Diastolic for various constants
	*/
	SysDiaCalc(NoPeaks, AvgHeartBeatPeakValue , HeartBeatPeak ,DC_Data) ;
	
	
	/*
	These files are used only for reference
	*/
	
	OutputFile1 = fopen("outputtest.dat","w+") ;
	OutputFile2= fopen("outputtest2.dat","w+") ;
	OutputFile3= fopen("outputtest3.dat","w+") ;
	OutputFile4= fopen("outputtest4.dat","w+") ;
	
	/*
		Writing data in file for reference
	*/
	
	for(	li = 0 ; li <NoPeaks; li++)
	{
		fprintf(OutputFile1 ,"%d\n",	 HeartBeatValley[li]-LeftCrop) ;
		fprintf(OutputFile2 ,"%f\n",	 HeartBeatValleyValue[li]) ;
		fprintf(OutputFile4 ,"%f\n",	 AvgHeartBeatValleyValue[li]) ;
		fprintf(OutputFile1 ,"%d\n",	 HeartBeatPeak[li]-LeftCrop) ;
		fprintf(OutputFile2 ,"%f\n",	 HeartBeatPeakValue[li]) ;
		fprintf(OutputFile4 ,"%f\n",	 AvgHeartBeatPeakValue[li]) ;
	}
	
	for(	li = LeftCrop ; li<NumberNonZeroSamplesAfterTrunc - RightCrop ; li++)
	{
		fprintf(OutputFile3 ,"%f\n",	 AC_Data[li]) ;
	}
	

	/*
		Closing the file variables
	*/
	
	fclose(OutputFile1) ;
	fclose(OutputFile2) ;
	fclose(OutputFile3) ;
	fclose(OutputFile4) ;
	
	printf("\n\n#################################################################################\n\n") ;


	/*
		To trigger python code for displaying values
	*/
		system("python plotgrpah.py") ;
	
	/*
		Free Arrays
	*/
	free(DC_Data);
	free(Trunc_DC_Data);
	free(Avg_DC_Data);
	free(AC_Data);
}
