/*! \file	Algorithm_DMO.h
	
	\brief	Framwork for Dynamic MOEA
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
		
	\date	Jul.20 2006 create
*/


#ifndef	AZ_ALGORITHM_DMO_H
#define	AZ_ALGORITHM_DMO_H

#include <ctime>
#include <vector>
#include "Fitting.h"
#include "LogFile.h"
#include "Parameter.h"
#include "PopulationMO.h"

//!\brief namespace of evolutionary algoirhtm 
namespace EA
{
//!\brief namespace of dynamic evolutionary algoirhtm 
namespace DEA
{

	const double PI = 3.141592653589793;
	double		 T	= 0.0;	
	
	// modified FDA problem
	void ZZJ(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
	{
		F[0] = X[0];
		double Gt = sin(0.5*PI*(T+1.0));
		//double Ht = Gt < -0.001 ? 0.5 : (Gt > 0.001 ? 2.0 : 1.0);
		double Ht = 1.5 + sin(0.5*PI*(T+1.0));
		double gx = 1.0;
		for(unsigned int i=1; i<X.size(); i++)
			gx +=  pow(X[i] + Gt - pow(X[0],Ht),2.0);
		F[1] = gx*(1.0 -  pow(F[0]/gx, Ht));
	}

	//	FDA1 test problem
	void FDA1(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
	{
		double Gt= sin(0.5*PI*(T+1.0));
		F[0] = X[0];
		double gx= 1.0;
		for(unsigned int i=1; i<X.size(); i++)
			gx += pow(X[i]-Gt,2.0);
		F[1] = gx*(1.0 -  pow(F[0]/gx, 0.5));
	}

	//!\brief	history node
	class HistoryNode
	{
	public:
		std::vector< unsigned int > parent, create;
		std::vector< double >		distance, alpha;
		CPopulationMO				pop;
	public:
		HistoryNode(CParameter&	par)
			: pop(par)
		{}
	};

	//!\brief	EA framework for DMO
	//!\param	COG	offspring generator
	//!\param	CES environmental selection
	template<class COG, class CES>
		class DMO : public COG, public CES
	{
	protected:
		unsigned int
			mPopSize,	//!< population size
			mStep	,	//!< current step
			mInStep,	//!< inner step of a time window
			mChange,	//!< current change
			mChangeMax,	//!< maximum changes
			mEvas	,	//!< the calculation number of the objectives
			mTaoT,		//!< the length of each state
			mNT,		//!< the number of distinct time steps
			mType,		//!< prediction types
			mHIndex;	//!< index to history record
		double 
			mTM,		//!< time used to model
			mTS;		//!< time used to selection
		std::vector< HistoryNode* >	pHistory;	//!< history record
		CPopulationMO	mPop,	//!< current population
						mTrial;	//!< trial population used by feedback strategy
		CParameter*	pPar;		//!< pointer to the parameter object
	public:
		//!\brief	constractor
		//!\param	stepmax		maximum steps
		//!\param	interval	the length of each state
		//!\param	timestep	the number of distinct time steps
		//!\param	type		prediction model type
		//!\param	par			parameter object
		//!\param	pop			initial population
		//!\param	problem		name of test problem
		//!\return	void
		DMO(
			unsigned int	stepmax	,
			unsigned int	interval,
			unsigned int	timestep,
			unsigned int	type,
			CParameter&		par		,
			CPopulationMO&	pop		,
			std::string		problem )
			: mPop(pop), mTrial(par)
		{
			mChangeMax 	= stepmax;
			pPar		= &par;
			mStep		= 0;
			mInStep		= 0;
			mChange		= 0;
			mEvas		= 0;
			mPopSize 	= mPop.Size();
			mTM			= mTS	= 0.0;
			mTaoT		= interval;
			mNT			= timestep;
			mType		= type;
			
			switch(mType)
			{
			case 4:
				pHistory.resize(0);
				break;
			case 5:
				pHistory.resize(3);
				break;
			default:
				pHistory.resize(2);
				break;
			}
			for(unsigned int i=0; i<pHistory.size(); i++)	pHistory[i] = 0;
			mHIndex		= 0;

			if(problem == std::string("ZZJ"))
				mPar.Evaluator( ZZJ );
			else
				mPar.Evaluator( FDA1 );
		}

		//!\brief	destructor
		//!\return	void
		~DMO() 
		{
			for(unsigned int i=0; i<pHistory.size(); i++) delete pHistory[i];
			pHistory.clear();
		}

		//!\brief	get the pointer to the parameter object
		//!\return	pointer to the parameter object
		inline CParameter& P() {return *pPar;}

		//!\brief	get the population
		//1\return	reference to population
		inline CPopulationMO& Population() {return mPop;}

		//!\brief	get the objective evaluation times
		//!\return	objective evaluation times
		inline unsigned int& EvaTimes() {return mEvas;}

		//!\brief	check to see whether the terminate condition is satified
		//!\return	true if terminate
		inline bool IsTerminate() {return mChange >= mChangeMax && mInStep>=mTaoT;}
		
		//!\brief	get the current step
		//!\return	current step
		inline unsigned int CurStep() {return mStep;}

		//!\brief	get the run time
		//!\return	void
		inline void Time(double& tm, double& ts) {tm = mTM/CLOCKS_PER_SEC; ts = mTS/CLOCKS_PER_SEC;}

		//!\brief	one step 
		//!\return	current step
		unsigned int Step()
		{
			clock_t t1=0,t2=0,t3=0;
			CPopulationMO popnew(P());

			if(ChangeState())
			{
				mEvas += Population().Size();
				mChange++;
				mInStep = 1;
			}
			else
			{
				t1 = clock();
				//sample new solutions
				COG::Generate(popnew, Population());

				//check the range of new solutions
				Check(popnew);

				//remove these repeated solutions
				Check(popnew, Population());

				//evaluate new solutions	
				popnew.Evaluate();
				mEvas += popnew.Size();	

				//popnew.Write("new.txt");

				t2 = clock();

				//environmental select 
				Population().Combine(popnew);

				mInStep++;
			}

			CES::Select(Population(), mPopSize);
			
			t3 = clock();

			mTM += double(t2-t1);
			mTS += double(t3-t2);

			return ++mStep;
		}


	protected:
		//!\brief	make all new solutions in feasible region
		//!\param	pop offspring population
		//!\return	void
		void Check(CPopulationMO& pop)
		{
			for(unsigned int i=0; i<pop.Size(); i++)
				pop[i].Check();
		}

		//!\brief	make all new solutions are different from old solutions
		//!\param	popnew offspring population
		//!\param	pop old population
		//!\return	void
		void Check(CPopulationMO& popnew, CPopulationMO& pop)
		{
			CPopulationMO tmp(P());
			for(unsigned int i=0; i<popnew.Size(); i++)
				if(!pop.IsContain(popnew[i])) tmp.Combine(popnew.At(i));
			popnew.Clear();
			popnew.Combine(tmp);
		}

		//!\brief	check whether it need to change the state
		//!\return	change the state or not
		bool ChangeState()
		{
			if(mStep % mTaoT == 0) 
			{
				T = floor(double(mStep)/double(mTaoT))/double(mNT);
				
				if(mStep > 0) InitTimeWindow();

				Population().Evaluate();

				return true;
			}
			else
			{
				return false;
			}
		}

		//!\biref	initialize time windows
		//!\return	void
		void InitTimeWindow()
		{
			// Step 1: record history
			RecordHistory();

			// Step 2: analyze history
			AnalyzeHistory();
			
			// Step 3: create new trial solutions
			CreateTrial();
		}

		//!\brief	record history information
		//!\return	void
		void RecordHistory()
		{
			if(pHistory.size()<1) return;

			unsigned int i,j,k; double dis;
			// history pool is not full
			if(mHIndex < pHistory.size())
			{
				if(pHistory[mHIndex] != 0 ) delete pHistory[mHIndex];
				pHistory[mHIndex] = new HistoryNode(Population().P());
				(*pHistory[mHIndex]).pop = Population();
				(*pHistory[mHIndex]).parent.resize(Population().Size());
				(*pHistory[mHIndex]).distance.resize(Population().Size());
				(*pHistory[mHIndex]).create.resize(Population().Size());
				(*pHistory[mHIndex]).alpha.resize(Population().Size());
				mHIndex++;
			}
			// history pool is full
			else
			{
				HistoryNode *th = pHistory[0];
				for(i=0; i<pHistory.size()-1; i++)
				{
					pHistory[i] = pHistory[i+1];
				}
				  pHistory[mHIndex-1]		= th;
				(*pHistory[mHIndex-1]).pop	= Population();
				(*pHistory[mHIndex-1]).parent.resize(Population().Size());
				(*pHistory[mHIndex-1]).distance.resize(Population().Size());
				(*pHistory[mHIndex-1]).create.resize(Population().Size());
				(*pHistory[mHIndex-1]).alpha.resize(Population().Size());
			}

			// initialize alpha
			for(i=0; i<(*pHistory[mHIndex-1]).pop.Size(); i++) (*pHistory[mHIndex-1]).alpha[i] = 1.0;

			// find parents
			if(mHIndex>1)
			{
				for(i=0; i<(*pHistory[mHIndex-1]).pop.Size(); i++)
				{
					(*pHistory[mHIndex-1]).distance[i]	= 1.0E200;
					(*pHistory[mHIndex-1]).parent[i]	= 0;

					for(j=0; j<(*pHistory[mHIndex-2]).pop.Size(); j++)
					{
						dis = 0.0;
						for(k=0; k<Population().P().XSize(); k++)
							dis += ((*pHistory[mHIndex-1]).pop[i][k] - (*pHistory[mHIndex-2]).pop[j][k])
								  *((*pHistory[mHIndex-1]).pop[i][k] - (*pHistory[mHIndex-2]).pop[j][k]);
						if(dis<(*pHistory[mHIndex-1]).distance[i])
						{
							(*pHistory[mHIndex-1]).distance[i]	= dis; 
							(*pHistory[mHIndex-1]).parent[i]	= j;
						}
					}
					(*pHistory[mHIndex-1]).distance[i] = sqrt((*pHistory[mHIndex-1]).distance[i]);
				}
			}
		}

		//!\brief	analyze the history pool
		//!\return	void
		void AnalyzeHistory()
		{
			// do not need analysis
			if(mType != 6 && mType != 7) return;

			// history pool is not full
			if(mHIndex < pHistory.size()) return;

			unsigned int success = 0, model = 0; 

			unsigned int i, j, p;
			double dotm, dotp, dotr;

			for(i=0; i<(*pHistory[mHIndex-1]).pop.Size(); i++)
			{
				p = (*pHistory[mHIndex-1]).parent[i];
				
				if((*pHistory[mHIndex-2]).create[p] != 1) continue;

				model++;

				dotm = dotp = dotr = 0.0;
				for(j=0; j<P().XSize(); j++) 
				{
					dotm += ((*pHistory[mHIndex-1]).pop[i][j] - (*pHistory[mHIndex-2]).pop[p][j])*(mTrial[p][j] - (*pHistory[mHIndex-2]).pop[p][j]);
					dotr += ((*pHistory[mHIndex-1]).pop[i][j] - (*pHistory[mHIndex-2]).pop[p][j])*((*pHistory[mHIndex-1]).pop[i][j] - (*pHistory[mHIndex-2]).pop[p][j]);
					dotp += (mTrial[p][j] - (*pHistory[mHIndex-2]).pop[p][j])*(mTrial[p][j] - (*pHistory[mHIndex-2]).pop[p][j]);
				}
				dotr = sqrt(dotr); dotp = sqrt(dotp);

				// check the prediction is right or wrong
				if(dotm < dotr*dotp*0.866 && dotm > 0)
				{
					success++;
					
					// adjust alpha
					(*pHistory[mHIndex-1]).alpha[i] = dotr/dotm;
				}
				else
				{
					// reset alpha
					(*pHistory[mHIndex-1]).alpha[i] = 1.0;
				}
			}
			
			// failed to predict new trial solutions
			// the environment has a great change
			if(success < 0.5*model)
			{
				std::swap(pHistory[mHIndex-1], pHistory[mHIndex-2]);
				mHIndex = 1;
			}
			//LOG::LogFile log;
			//log<<model<<"\t"<<success<<std::endl;
		}

		//!\brief	Create new trial locations
		//!\return	void
		void CreateTrial()
		{
			unsigned int i,k;
			double std;
			std::vector<double> T(3), X(3), C(3);
			T[0] = 0.0; T[1] = 1.0; T[2] = 2.0;

			// history pool is not full
			if(mHIndex < pHistory.size())
			{				
				for(i=0; i<Population().Size(); i++)
				{
					if(mHIndex>0)	std = (*pHistory[mHIndex-1]).distance[i] * sqrt(0.25/double(Population().P().XSize()));
					
					for(k=0; k<Population().P().XSize(); k++)
					{
						if(mHIndex<1)
							Population()[i][k] = rnd::rand(Population().P().XLow(k),Population().P().XUpp(k));
						else
							Population()[i][k] = (*pHistory[mHIndex-1]).pop[i][k] + std*rnd::gaussian();
						// boundary check strategy
						if(Population()[i][k]      < Population().P().XLow(k))	Population()[i][k] = 0.5*(Population().P().XLow(k)+(*pHistory[mHIndex-1]).pop[i][k]);
						else if(Population()[i][k] > Population().P().XUpp(k))	Population()[i][k] = 0.5*(Population().P().XUpp(k)+(*pHistory[mHIndex-1]).pop[i][k]);
					}
					if(mHIndex>0) (*pHistory[mHIndex-1]).create[i] = 0;
				}
				return;
			}

			// history pool is full
			switch(mType)
			{
			//half prediction + half variation
			case 1:
			case 6:
				mTrial.Resize(Population().Size());
				for(i=0; i<Population().Size(); i++)
				{
					if(mHIndex>0)	std = (*pHistory[mHIndex-1]).distance[i] * sqrt(0.25/double(Population().P().XSize()));
					else			std = 0.0;

					if(rnd::rand(0.0,1.0)<0.5)
					{
						//predict the next location
						for(k=0; k<Population().P().XSize(); k++)
						{
							mTrial[i][k] =  (*pHistory[mHIndex-1]).pop[i][k] 
											+(*pHistory[mHIndex-1]).alpha[i] 
											*(  (*pHistory[mHIndex-1]).pop[i][k] 
											  - (*pHistory[mHIndex-2]).pop[(*pHistory[mHIndex-1]).parent[i]][k]);
							Population()[i][k] = mTrial[i][k] + std*rnd::gaussian();
							
							// boundary check strategy
							if(Population()[i][k]      < Population().P().XLow(k))	Population()[i][k] = 0.5*(Population().P().XLow(k)+(*pHistory[mHIndex-1]).pop[i][k]);
							else if(Population()[i][k] > Population().P().XUpp(k))	Population()[i][k] = 0.5*(Population().P().XUpp(k)+(*pHistory[mHIndex-1]).pop[i][k]);
						}
						(*pHistory[mHIndex-1]).create[i] = 1;
					}
					else
					{
						//variation
						for(k=0; k<Population().P().XSize(); k++)
						{
							mTrial[i][k] = (*pHistory[mHIndex-1]).pop[i][k];
							Population()[i][k] = mTrial[i][k] + std*rnd::gaussian();
							
							// boundary check strategy
							if(Population()[i][k]      < Population().P().XLow(k))	Population()[i][k] = 0.5*(Population().P().XLow(k)+(*pHistory[mHIndex-1]).pop[i][k]);
							else if(Population()[i][k] > Population().P().XUpp(k))	Population()[i][k] = 0.5*(Population().P().XUpp(k)+(*pHistory[mHIndex-1]).pop[i][k]);
						}
						(*pHistory[mHIndex-1]).create[i] = 2;
					}
				}
				break;
			//only prediction
			//predict the next location
			case 2:
			case 7:
				mTrial.Resize(Population().Size());
				for(i=0; i<Population().Size(); i++)
				{
					if(mHIndex>0)	std = (*pHistory[mHIndex-1]).distance[i] * sqrt(0.25/double(Population().P().XSize()));
					else			std = 0.0;
					for(k=0; k<Population().P().XSize(); k++)
					{
						mTrial[i][k] =  (*pHistory[mHIndex-1]).pop[i][k] 
										+(*pHistory[mHIndex-1]).alpha[i] 
										*(  (*pHistory[mHIndex-1]).pop[i][k] 
											- (*pHistory[mHIndex-2]).pop[(*pHistory[mHIndex-1]).parent[i]][k]);
						Population()[i][k] = mTrial[i][k] + std*rnd::gaussian();
						
						// boundary check strategy
						if(Population()[i][k]      < Population().P().XLow(k))	Population()[i][k] = 0.5*(Population().P().XLow(k)+(*pHistory[mHIndex-1]).pop[i][k]);
						else if(Population()[i][k] > Population().P().XUpp(k))	Population()[i][k] = 0.5*(Population().P().XUpp(k)+(*pHistory[mHIndex-1]).pop[i][k]);
					}
					(*pHistory[mHIndex-1]).create[i] = 1;
				}
				break;
			// only variation
			case 3:
				for(i=0; i<Population().Size(); i++)
				{
					if(mHIndex>0)	std = (*pHistory[mHIndex-1]).distance[i] * sqrt(0.25/double(Population().P().XSize()));
					else			std = 0.0;
					for(k=0; k<Population().P().XSize(); k++)
					{
						Population()[i][k] = (*pHistory[mHIndex-1]).pop[i][k] + std*rnd::gaussian();
						
						// boundary check strategy
						if(Population()[i][k]      < Population().P().XLow(k))	Population()[i][k] = 0.5*(Population().P().XLow(k)+(*pHistory[mHIndex-1]).pop[i][k]);
						else if(Population()[i][k] > Population().P().XUpp(k))	Population()[i][k] = 0.5*(Population().P().XUpp(k)+(*pHistory[mHIndex-1]).pop[i][k]);
					}
					(*pHistory[mHIndex-1]).create[i] = 2;
				}
				break;
			// radndom generation 
			case 4:
			default:
				for(i=0; i<Population().Size(); i++)
				{
					for(k=0; k<Population().P().XSize(); k++) 
						Population()[i][k] = rnd::rand(Population().P().XLow(k), Population().P().XUpp(k));
					//(*pHistory[mHIndex-1]).create[i] = 0;
				}
				break;
			// quadratic model prediction
			case 5:
				for(i=0; i<Population().Size(); i++)
				{
					if(mHIndex>0)	std = (*pHistory[mHIndex-1]).distance[i] * sqrt(0.25/double(Population().P().XSize()));
					else			std = 0.0;
					for(k=0; k<Population().P().XSize(); k++) 
					{
						unsigned int h1,h2;
						h1 = (*pHistory[mHIndex-1]).parent[i];
						h2 = (*pHistory[mHIndex-2]).parent[h1];
						X[2] = (*pHistory[mHIndex-1]).pop[i][k];
						X[1] = (*pHistory[mHIndex-2]).pop[h1][k];
						X[0] = (*pHistory[mHIndex-3]).pop[h2][k];
						
						LEARN::poly_fit( T, X, 2, C );

						Population()[i][k] = C[0] + 3.0*C[1] + 9.0*C[2];
												+ std*rnd::gaussian();
						
						// boundary check strategy
						if(Population()[i][k]      < Population().P().XLow(k))	Population()[i][k] = 0.5*(Population().P().XLow(k)+(*pHistory[mHIndex-1]).pop[i][k]);
						else if(Population()[i][k] > Population().P().XUpp(k))	Population()[i][k] = 0.5*(Population().P().XUpp(k)+(*pHistory[mHIndex-1]).pop[i][k]);
					}
					(*pHistory[mHIndex-1]).create[i] = 1;
				}
				break;
			}
		}
	}; //class DMO

} //namespace DEA
} //namespace EA

#endif //AZ_ALGORITHM_DMO_H
