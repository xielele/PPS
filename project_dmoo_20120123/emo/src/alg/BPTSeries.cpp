#include "alg/BPTSeries.h"
#include "alg/Random.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>

namespace bpts
{
	/******************************************************************************

	====================================================
	Network:      Backpropagation Network with Bias Terms and Momentum
	====================================================

	Application:  Time-Series Forecasting
	Prediction of the Annual Number of Sunspots

	Author:       Karsten Kutza
	Date:         17.4.96

	Reference:    D.E. Rumelhart, G.E. Hinton, R.J. Williams
	Learning Internal Representations by Error Propagation
	in:
	D.E. Rumelhart, J.L. McClelland (Eds.)
	Parallel Distributed Processing, Volume 1
	MIT Press, Cambridge, MA, pp. 318-362, 1986

	******************************************************************************/

	/******************************************************************************
	D E C L A R A T I O N S
	******************************************************************************/
	typedef int           BOOL;
	typedef int           INT;
	typedef double        REAL;

#define FALSE         0
#define TRUE          1
#define NOT           !
#define AND           &&
#define OR            ||

#define MIN_REAL      -HUGE_VAL
#define MAX_REAL      +HUGE_VAL
#define MIN(x,y)      ((x)<(y) ? (x) : (y))
#define MAX(x,y)      ((x)>(y) ? (x) : (y))

#define LO            0.1
#define HI            0.9
#define BIAS          1

#define sqr(x)        ((x)*(x))

	typedef struct {                     /* A LAYER OF A NET:                     */
		INT           Units;         /* - number of units in this layer       */
		REAL*         Output;        /* - output of ith unit                  */
		REAL*         Error;         /* - error term of ith unit              */
		REAL**        Weight;        /* - connection weights to ith unit      */
		REAL**        WeightSave;    /* - saved weights for stopped training  */
		REAL**        dWeight;       /* - last weight deltas for momentum     */
	} LAYER;

	typedef struct {                     /* A NET:                                */
		LAYER**       Layer;         /* - layers of this net                  */
		LAYER*        InputLayer;    /* - input layer                         */
		LAYER*        OutputLayer;   /* - output layer                        */
		REAL          Alpha;         /* - momentum factor                     */
		REAL          Eta;           /* - learning rate                       */
		REAL          Gain;          /* - gain of sigmoid function            */
		REAL          Error;         /* - total net error                     */
	} NET;


	/******************************************************************************
	R A N D O M S   D R A W N   F R O M   D I S T R I B U T I O N S
	******************************************************************************/
	void InitializeRandoms()
	{
//		srand(4711);
	}


	INT RandomEqualINT(INT Low, INT High)
	{
		return az::rnd::rand(Low, High+1);
		//return rand() % (High-Low+1) + Low;
	}


	REAL RandomEqualREAL(REAL Low, REAL High)
	{
		return az::rnd::rand(Low, High);
		//return ((REAL) rand() / RAND_MAX) * (High-Low) + Low;
	}


	/******************************************************************************
	A P P L I C A T I O N - S P E C I F I C   C O D E
	******************************************************************************/
#define NUM_LAYERS    3
#define N             5
#define M             1
INT                   Units[NUM_LAYERS] = {N, 10, M};

#define FIRST_YEAR    1700

	int NUM_YEARS,
		TRAIN_LWB,
		TRAIN_UPB,
		TRAIN_YEARS,
		TEST_LWB,
		TEST_UPB,
		TEST_YEARS,
		EVAL_LWB,
		EVAL_UPB,
		EVAL_YEARS;
	REAL *Sunspots_,*Sunspots;

	REAL                  Mean;
	REAL                  TrainError;
	REAL                  TrainErrorPredictingMean;
	REAL                  TestError;
	REAL                  TestErrorPredictingMean;
//	FILE*                 f;
	REAL				  MaxD, MinD;

	void NormalizeSunspots()
	{
		INT  Year;

		MinD = MAX_REAL;
		MaxD = MIN_REAL;
		for (Year=0; Year<NUM_YEARS; Year++) {
			MinD = MIN(MinD, Sunspots[Year]);
			MaxD = MAX(MaxD, Sunspots[Year]);
		}
		Mean = 0;
		for (Year=0; Year<NUM_YEARS; Year++) {
			Sunspots_[Year] =
				Sunspots [Year] = ((Sunspots[Year]-MinD) / (MaxD-MinD)) * (HI-LO) + LO;
			Mean += Sunspots[Year] / NUM_YEARS;
		}
	}


	void InitializeApplication(NET* Net)
	{
		INT  Year, i;
		REAL Out, Err;

		Net->Alpha = 0.5;
		Net->Eta   = 0.05;
		Net->Gain  = 1;

		NormalizeSunspots();
		TrainErrorPredictingMean = 0;
		for (Year=TRAIN_LWB; Year<=TRAIN_UPB; Year++) {
			for (i=0; i<M; i++) {
				Out = Sunspots[Year+i];
				Err = Mean - Out;
				TrainErrorPredictingMean += 0.5 * sqr(Err);
			}
		}
		TestErrorPredictingMean = 0;
		for (Year=TEST_LWB; Year<=TEST_UPB; Year++) {
			for (i=0; i<M; i++) {
				Out = Sunspots[Year+i];
				Err = Mean - Out;
				TestErrorPredictingMean += 0.5 * sqr(Err);
			}
		}
//		f = fopen("BPN.txt", "w");
	}

	void FinalizeApplication(NET* Net)
	{
		INT l,i;
		for (l=0; l<NUM_LAYERS; l++)
		{
			if (l != 0) 
			{
				for (i=1; i<=Units[l]; i++) 
				{
					free(Net->Layer[l]->Weight[i]);
					free(Net->Layer[l]->WeightSave[i]);
					free(Net->Layer[l]->dWeight[i]);
				}
			}
			free(Net->Layer[l]->Output);
			free(Net->Layer[l]->Error);
			free(Net->Layer[l]->Weight);
			free(Net->Layer[l]->WeightSave);
			free(Net->Layer[l]->dWeight);
			free(Net->Layer[l]);
		}
		free(Net->Layer);

		free(Sunspots_);
		free(Sunspots);

//		fclose(f);
	}


	/******************************************************************************
	I N I T I A L I Z A T I O N
	******************************************************************************/


	void GenerateNetwork(NET* Net)
	{
		INT l,i;

		Net->Layer = (LAYER**) calloc(NUM_LAYERS, sizeof(LAYER*));

		for (l=0; l<NUM_LAYERS; l++) {
			Net->Layer[l] = (LAYER*) malloc(sizeof(LAYER));

			Net->Layer[l]->Units      = Units[l];
			Net->Layer[l]->Output     = (REAL*)  calloc(Units[l]+1, sizeof(REAL));
			Net->Layer[l]->Error      = (REAL*)  calloc(Units[l]+1, sizeof(REAL));
			Net->Layer[l]->Weight     = (REAL**) calloc(Units[l]+1, sizeof(REAL*));
			Net->Layer[l]->WeightSave = (REAL**) calloc(Units[l]+1, sizeof(REAL*));
			Net->Layer[l]->dWeight    = (REAL**) calloc(Units[l]+1, sizeof(REAL*));
			Net->Layer[l]->Output[0]  = BIAS;

			if (l != 0) {
				for (i=1; i<=Units[l]; i++) {
					Net->Layer[l]->Weight[i]     = (REAL*) calloc(Units[l-1]+1, sizeof(REAL));
					Net->Layer[l]->WeightSave[i] = (REAL*) calloc(Units[l-1]+1, sizeof(REAL));
					Net->Layer[l]->dWeight[i]    = (REAL*) calloc(Units[l-1]+1, sizeof(REAL));
				}
			}
		}
		Net->InputLayer  = Net->Layer[0];
		Net->OutputLayer = Net->Layer[NUM_LAYERS - 1];
		Net->Alpha       = 0.9;
		Net->Eta         = 0.25;
		Net->Gain        = 1;
	}


	void RandomWeights(NET* Net)
	{
		INT l,i,j;

		for (l=1; l<NUM_LAYERS; l++) {
			for (i=1; i<=Net->Layer[l]->Units; i++) {
				for (j=0; j<=Net->Layer[l-1]->Units; j++) {
					Net->Layer[l]->Weight[i][j] = RandomEqualREAL(-0.5, 0.5);
				}
			}
		}
	}


	void SetInput(NET* Net, REAL* Input)
	{
		INT i;

		for (i=1; i<=Net->InputLayer->Units; i++) {
			Net->InputLayer->Output[i] = Input[i-1];
		}
	}


	void GetOutput(NET* Net, REAL* Output)
	{
		INT i;

		for (i=1; i<=Net->OutputLayer->Units; i++) {
			Output[i-1] = Net->OutputLayer->Output[i];
		}
	}


	/******************************************************************************
	S U P P O R T   F O R   S T O P P E D   T R A I N I N G
	******************************************************************************/


	void SaveWeights(NET* Net)
	{
		INT l,i,j;

		for (l=1; l<NUM_LAYERS; l++) {
			for (i=1; i<=Net->Layer[l]->Units; i++) {
				for (j=0; j<=Net->Layer[l-1]->Units; j++) {
					Net->Layer[l]->WeightSave[i][j] = Net->Layer[l]->Weight[i][j];
				}
			}
		}
	}


	void RestoreWeights(NET* Net)
	{
		INT l,i,j;

		for (l=1; l<NUM_LAYERS; l++) {
			for (i=1; i<=Net->Layer[l]->Units; i++) {
				for (j=0; j<=Net->Layer[l-1]->Units; j++) {
					Net->Layer[l]->Weight[i][j] = Net->Layer[l]->WeightSave[i][j];
				}
			}
		}
	}


	/******************************************************************************
	P R O P A G A T I N G   S I G N A L S
	******************************************************************************/


	void PropagateLayer(NET* Net, LAYER* Lower, LAYER* Upper)
	{
		INT  i,j;
		REAL Sum;

		for (i=1; i<=Upper->Units; i++) {
			Sum = 0;
			for (j=0; j<=Lower->Units; j++) {
				Sum += Upper->Weight[i][j] * Lower->Output[j];
			}
			Upper->Output[i] = 1 / (1 + exp(-Net->Gain * Sum));
		}
	}

	void PropagateNet(NET* Net)
	{
		INT l;

		for (l=0; l<NUM_LAYERS-1; l++) {
			PropagateLayer(Net, Net->Layer[l], Net->Layer[l+1]);
		}
	}


	/******************************************************************************
	B A C K P R O P A G A T I N G   E R R O R S
	******************************************************************************/


	void ComputeOutputError(NET* Net, REAL* Target)
	{
		INT  i;
		REAL Out, Err;

		Net->Error = 0;
		for (i=1; i<=Net->OutputLayer->Units; i++) {
			Out = Net->OutputLayer->Output[i];
			Err = Target[i-1]-Out;
			Net->OutputLayer->Error[i] = Net->Gain * Out * (1-Out) * Err;
			Net->Error += 0.5 * sqr(Err);
		}
	}


	void BackpropagateLayer(NET* Net, LAYER* Upper, LAYER* Lower)
	{
		INT  i,j;
		REAL Out, Err;

		for (i=1; i<=Lower->Units; i++) {
			Out = Lower->Output[i];
			Err = 0;
			for (j=1; j<=Upper->Units; j++) {
				Err += Upper->Weight[j][i] * Upper->Error[j];
			}
			Lower->Error[i] = Net->Gain * Out * (1-Out) * Err;
		}
	}


	void BackpropagateNet(NET* Net)
	{
		INT l;

		for (l=NUM_LAYERS-1; l>1; l--) {
			BackpropagateLayer(Net, Net->Layer[l], Net->Layer[l-1]);
		}
	}


	void AdjustWeights(NET* Net)
	{
		INT  l,i,j;
		REAL Out, Err, dWeight;

		for (l=1; l<NUM_LAYERS; l++) {
			for (i=1; i<=Net->Layer[l]->Units; i++) {
				for (j=0; j<=Net->Layer[l-1]->Units; j++) {
					Out = Net->Layer[l-1]->Output[j];
					Err = Net->Layer[l]->Error[i];
					dWeight = Net->Layer[l]->dWeight[i][j];
					Net->Layer[l]->Weight[i][j] += Net->Eta * Err * Out + Net->Alpha * dWeight;
					Net->Layer[l]->dWeight[i][j] = Net->Eta * Err * Out;
				}
			}
		}
	}


	/******************************************************************************
	S I M U L A T I N G   T H E   N E T
	******************************************************************************/


	void SimulateNet(NET* Net, REAL* Input, REAL* Output, REAL* Target, BOOL Training)
	{
		SetInput(Net, Input);
		PropagateNet(Net);
		GetOutput(Net, Output);

		ComputeOutputError(Net, Target);
		if (Training) {
			BackpropagateNet(Net);
			AdjustWeights(Net);
		}
	}


	void TrainNet(NET* Net, INT Epochs)
	{
		INT  Year, n;
		REAL Output[M];

		for (n=0; n<Epochs*TRAIN_YEARS; n++) {
			Year = RandomEqualINT(TRAIN_LWB, TRAIN_UPB);
			SimulateNet(Net, &(Sunspots[Year-N]), Output, &(Sunspots[Year]), TRUE);
		}
	}


	void TestNet(NET* Net)
	{
		INT  Year;
		REAL Output[M];

		TrainError = 0;
		for (Year=TRAIN_LWB; Year<=TRAIN_UPB; Year++) {
			SimulateNet(Net, &(Sunspots[Year-N]), Output, &(Sunspots[Year]), FALSE);
			TrainError += Net->Error;
		}
		TestError = 0;
		for (Year=TEST_LWB; Year<=TEST_UPB; Year++) {
			SimulateNet(Net, &(Sunspots[Year-N]), Output, &(Sunspots[Year]), FALSE);
			TestError += Net->Error;
		}
//		fprintf(f, "\nNMSE is %0.3f on Training Set and %0.3f on Test Set",
//			TrainError / TrainErrorPredictingMean,
//			TestError / TestErrorPredictingMean);
	}


	void EvaluateNet(NET* Net)
	{
		INT  Year;
		REAL Output [M];
		REAL Output_[M];

//		fprintf(f, "\n\n\n");
//		fprintf(f, "Year    Sunspots    Open-Loop Prediction    Closed-Loop Prediction\n");
//		fprintf(f, "\n");
		for (Year=EVAL_LWB; Year<=EVAL_UPB; Year++) {
			SimulateNet(Net, &(Sunspots [Year-N]), Output,  &(Sunspots [Year]), FALSE);
			SimulateNet(Net, &(Sunspots_[Year-N]), Output_, &(Sunspots_[Year]), FALSE);
			Sunspots_[Year] = Output_[0];
//			fprintf(f, "%d       %0.3f                   %0.3f                     %0.3f\n",
//				FIRST_YEAR + Year,//
//				Sunspots[Year],
//				Output [0],
//				Output_[0]);
		}
	}


	/******************************************************************************
	M A I N
	******************************************************************************/
	double BPPredict(double *pdata, int num, double& std)
	{
		NUM_YEARS = num;
		TRAIN_LWB = N;
		TRAIN_UPB = num-1;
		TRAIN_YEARS = (TRAIN_UPB - TRAIN_LWB + 1);
		TEST_LWB    = N;
		TEST_UPB    = num-1;
		TEST_YEARS  = (TEST_UPB - TEST_LWB + 1);
		EVAL_LWB    = N;
		EVAL_UPB    = num-1;
		EVAL_YEARS  = (EVAL_UPB - EVAL_LWB + 1);
		Sunspots	= (REAL*)  calloc(NUM_YEARS, sizeof(REAL));
		Sunspots_	= (REAL*)  calloc(NUM_YEARS, sizeof(REAL));
		for(int i=0; i<num; i++)
			Sunspots[i] = pdata[i];

		NET  Net;
		BOOL Stop;
		REAL MinTestError;

		InitializeRandoms();
		GenerateNetwork(&Net);
		RandomWeights(&Net);
		InitializeApplication(&Net);

		Stop = FALSE;
		MinTestError = MAX_REAL;
		long ss = 0;
		do 
		{
			TrainNet(&Net, 1);
			TestNet(&Net);
			if (TestError < MinTestError) 
			{
//				fprintf(f, " - saving Weights ...");
				MinTestError = TestError;
				SaveWeights(&Net);
			}
			else if (TestError > 1.2 * MinTestError) 
			{
//				fprintf(f, " - stopping Training and restoring Weights ...");
				//Stop = TRUE;
				RestoreWeights(&Net);
			}
			if(++ss > 1000) Stop = TRUE;
		} while (NOT Stop);

//		TestNet(&Net);
//		EvaluateNet(&Net);

		INT  Year;
		REAL Output[M];
		std = 0;
		for (Year=TRAIN_LWB; Year<=TRAIN_UPB; Year++) {
			SimulateNet(&Net, &(Sunspots[Year-N]), Output, &(Sunspots[Year]), FALSE);
			std += sqr(Output[0] - Sunspots[Year]);
		}
		std = sqrt(std/(NUM_YEARS-N+0.0));

		SimulateNet(&Net, &(Sunspots[NUM_YEARS-N]), Output, &(Sunspots[NUM_YEARS-1]), FALSE);
		Output[0] = (Output[0]-LO)/(HI-LO) * (MaxD-MinD) + MinD;

		FinalizeApplication(&Net);

		return Output[0];
	}
};