# include <stdio.h>
# include <problem.h>
# include <doublepop.h>
# include <math.h>
# include <unistd.h>
# include <sys/times.h>
# include <time.h>
# include <getoptions.h>

extern "C"
{
extern int getdimension();
extern int geteq();
extern int getineq();
extern void getleftmargin(double *x);
extern void getrightmargin(double *x);
extern double funmin(double *x);
extern void   granal(double *x,double *g);
extern void	feq(double *x,double *k);
extern	void	fineq(double *x,double *k);
extern void	done(double *x);
}
int problem_dimension;

double	*tempx,*tempg;
double	*temp_eq;
double	*temp_ineq;

class MyProblem :public Problem
{
	public:
		MyProblem(int N,int E,int NE);
		virtual double funmin(Data x);
		virtual void   granal(Data x,Data &g);
		virtual void	feq(Data x,Data &k);
		virtual void	fineq(Data x,Data &k);
		virtual void	done(Data x);
};

MyProblem::MyProblem(int N,int E,int NE)
	:Problem(N,E,NE)
{
}


double	MyProblem::funmin(Data x)
{
	++fevals;
	for(int i=0;i<problem_dimension;i++) tempx[i]=x[i];
	return ::funmin(tempx);
}

void	MyProblem::granal(Data x,Data &g)
{
	++gevals;
	for(int i=0;i<problem_dimension;i++) tempx[i]=x[i];
	::granal(tempx,tempg);
	for(int i=0;i<problem_dimension;i++) g[i]=tempg[i];
}

void	MyProblem::feq(Data x,Data &k)
{
	++eqevals;
	for(int i=0;i<problem_dimension;i++) tempx[i]=x[i];
	::feq(tempx,temp_eq);
	for(int i=0;i<neq;i++) k[i]=temp_eq[i];
}

void	MyProblem::fineq(Data x,Data &k)
{
	++neqevals;
	for(int i=0;i<problem_dimension;i++) tempx[i]=x[i];
	::fineq(tempx,temp_ineq);
	for(int i=0;i<nineq;i++) k[i]=temp_ineq[i];
}


void	MyProblem::done(Data x)
{
	for(int i=0;i<problem_dimension;i++) tempx[i]=x[i];
	::done(tempx);
}

int main(int argc, char **argv)
{
	parse_cmd_line(argc,argv);
	srand(random_seed);
	srand48(random_seed);

#ifndef BATCH
	problem_dimension=getdimension();
	int problem_eq=geteq();
	int problem_ineq=getineq();
	tempx=    new double[problem_dimension];
	tempg=    new double[problem_dimension];
	temp_eq=  new double[problem_eq];
	temp_ineq=new double[problem_ineq];
	MyProblem myproblem(problem_dimension,problem_eq,problem_ineq);
	double *l,*r;
	l=new double[problem_dimension];
	r=new double[problem_dimension];
	getleftmargin(l);
	getrightmargin(r);
	Data L,R;
	L.resize(problem_dimension);
	R.resize(problem_dimension);
	for(int i=0;i<problem_dimension;i++)
	{
		L[i]=l[i];
		R[i]=r[i];
	}
	myproblem.setLeftMargin(L);
	myproblem.setRightMargin(R);
	delete[] l;
	delete[] r;
	DoublePop pop(chromosomes,&myproblem);
	pop.setSelectionRate(selection_rate);
	pop.setMutationRate(mutation_rate);
	pop.Solve();
		
	printf("TOTAL FUNCTION CALLS = %8d\n",myproblem.getfevals());
	delete[] tempx;
	delete[] tempg;
	delete[] temp_eq;
	delete[] temp_ineq;
#else
	const int max_runs=30;
	double avg_fevals=0.0;
	double avg_eqevals=0.0;
	double avg_neqevals=0.0;
	double min_value=1e+100;
	double max_value=-1e+100;
	double xx1=0.0;
	double xx2=0.0;
	for(random_seed=1;random_seed<=max_runs;random_seed++)
	{
		srand(random_seed);
		srand48(random_seed);
		problem_dimension=getdimension();
		int problem_eq=geteq();
		int problem_ineq=getineq();
		tempx=    new double[problem_dimension];
		tempg=    new double[problem_dimension];
		temp_eq=  new double[problem_eq];
		temp_ineq=new double[problem_ineq];
		MyProblem myproblem(problem_dimension,problem_eq,problem_ineq);
		double *l,*r;
		l=new double[problem_dimension];
		r=new double[problem_dimension];
		getleftmargin(l);
		getrightmargin(r);
		Data L,R;
		L.resize(problem_dimension);
		R.resize(problem_dimension);
		Data xbest;
		xbest.resize(problem_dimension);
		double ybest;
		for(int i=0;i<problem_dimension;i++)
		{
			L[i]=l[i];
			R[i]=r[i];
		}
		myproblem.setLeftMargin(L);
		myproblem.setRightMargin(R);
		delete[] l;
		delete[] r;
		DoublePop pop(chromosomes,&myproblem);
		pop.setSelectionRate(selection_rate);
		pop.setMutationRate(mutation_rate);
		pop.Solve();
		pop.getMinimum(xbest,ybest);
		xx1+=ybest;
		xx2+=ybest * ybest;
		if(ybest<min_value) min_value = ybest;
		if(ybest>max_value) max_value = ybest;
		delete[] tempx;
		delete[] tempg;
		delete[] temp_eq;
		delete[] temp_ineq;
		avg_fevals+=myproblem.getfevals();
		avg_eqevals+=myproblem.geteqevals();
		avg_neqevals+=myproblem.getneqevals();
		printf("RUN:%4d\tFEVALS=%10.2lf\tEQEVALS=%10.2lf\tNEQEVALS=%10.2lf\n",
				random_seed,avg_fevals/random_seed,
				avg_eqevals/random_seed,avg_neqevals/random_seed);
		fflush(stdout);
	}
	printf("MINIMUM VALUE=%.20lf\n",min_value);
	printf("MAXIMUM VALUE=%.20lf\n",max_value);
	printf("MEAN    VALUE=%.20lf\n",xx1/max_runs);
	printf("STD     VALUE=%.20lg\n",sqrt(fabs(xx2/max_runs-xx1/max_runs*xx1/max_runs)));
#endif
	return 0;
}
