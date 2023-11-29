
//////////////////////////////////////////////
//	Fully connected multilayered feed		//
//	forward	artificial neural network using	//
//	Backpropogation	algorithm for training.	//
//////////////////////////////////////////////


#ifndef backprop_h
#define backprop_h

#include<cmath>

namespace az
{

namespace alg
{
class CBackProp{

//	output of each neuron
	double **out;

//	delta error value for each neuron
	double **delta;

//	vector of weights for each neuron
	double ***weight;

//	no of layers in net
//	including input layer
	int numl;

//	vector of numl elements for size 
//	of each layer
	int *lsize;

//	learning rate
	double beta;

//	momentum parameter
	double alpha;

//	storage for weight-change made
//	in previous epoch
	double ***prevDwt;

//	squashing function
	double sigmoid(double in);

public:

	~CBackProp();

//	initializes and allocates memory
	CBackProp(int nl,int *sz,double b,double a);

//	backpropogates error for one set of input
	void bpgt(double *in,double *tgt);

//	feed forwards activations for one set of inputs
	void ffwd(double *in);

//	returns mean square error of the net
	double mse(double *tgt) const;	
	
//	returns i'th output of the net
	double Out(int i) const;
};
}
}
#endif