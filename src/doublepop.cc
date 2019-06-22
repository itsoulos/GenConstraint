# include <doublepop.h>
# include <iostream>
# include <math.h>
# include <stdio.h>

double REND=1e-5;

static double dmax(double a,double b) {return a>b?a:b;}
static double diff(double *x,double *y,int n)
{
	double s=0.0;
	for(int i=0;i<n;i++) s+=pow(x[i]-y[i],2.0);
	return sqrt(s);
}


static double getDistance(Data x1,Data x2)
{
	double s=0.0;
	for(int i=0;i<x1.size();i++) s+=(x1[i]-x2[i])*(x1[i]-x2[i]);
	return s;
}


DoublePop::DoublePop(int gcount,Problem *p)
	:Base(p,gcount)
{
	genome_size=p->getDimension();
	elitism=1;
	selection_rate = 0.1;
	mutation_rate  = 0.05;
	genome_count   = gcount;
	generation     = 0;

	double f;
	children.resize(genome_count);
	parent0.resize(genome_size);
	parent1.resize(genome_size);
	for(int i=0;i<genome_count;i++) children[i].resize(genome_size);
	lmargin.resize(p->getDimension());
	rmargin.resize(p->getDimension());
	lmargin=problem->getLeftMargin();
	rmargin=problem->getRightMargin();
	isFeasible.resize(genome_count);
	countFeasible=0;
	for(int i=0;i<genome_count;i++) 
	{
		isFeasible[i]=problem->isFeasible(genome[i]);
		countFeasible+=isFeasible[i];
	}
}

void	DoublePop::reset()
{
}

static double G(double x)
{
	return x>=0?0:1;
}

void	DoublePop::select()
{
	Data itemp;
	itemp.resize(genome_size);
	for(int i=0;i<genome_count;i++)
	{
		for(int j=0;j<genome_count-1;j++)
		{
			if((fitness[j+1]<fitness[j])
			 || (isFeasible[j+1] && !isFeasible[j]))
			{
				//if(isFeasible[j+1]==0 && isFeasible[j]) continue;
				double dtemp;
				dtemp=fitness[j];
				fitness[j]=fitness[j+1];
				fitness[j+1]=dtemp;

				itemp=genome[j];
				genome[j]=genome[j+1];
				genome[j+1]=itemp;
			
				int ik;
				ik=isFeasible[j];
				isFeasible[j]=isFeasible[j+1];
				isFeasible[j+1]=ik;
			}
		}
	}
}

void	DoublePop::getTournamentElement(Data &x)
{
        const int tournament_size =(genome_count<=100)?4:10;
	double max_fitness=1e+100;
	int    max_index=-1;
	int    status=0;
		
	if(countFeasible && rand() % 2==1)
	{
		vector<int> Index;
		int counter=0;
		do
		{
			int r=rand() % genome_count;
			if(isFeasible[r]) {Index.push_back(r);counter++;}
		}while(counter<tournament_size);
		for(int i=0;i<Index.size();i++)
		{
			int r=Index[i];
			if(fitness[r]<max_fitness) {max_index=r;max_fitness=fitness[r];}
		}
	}
	
	if(max_index==-1) 
	{
		for(int i=0;i<tournament_size;i++)
		{
		int r=rand() % genome_count;
		if(fitness[r]<max_fitness) {max_index=r;max_fitness=fitness[r];}
		}
	}
	x=genome[max_index];
}

void	DoublePop::tournament(Data &p1,Data &p2)
{
	getTournamentElement(p1);
	getTournamentElement(p2);
}

void	DoublePop::crossover()
{
        int nchildren=(int)((1.0 - selection_rate) * genome_count);
	if(!(nchildren%2==0)) nchildren++;
	
        int count_children=0;
        while(1)
        {
		tournament(parent0,parent1);
		for(int i=0;i<genome_size;i++)
		{
			double alfa,b,u,g1,g2;
			double x1,x2;
			int p1,p2;
			x1=parent0[i];
			x2=parent1[i];
			alfa=-0.5+2.0*drand48();
			g1=alfa*x1+(1.0-alfa)*x2;
			g2=alfa*x2+(1.0-alfa)*x1;
			if(g1>rmargin[i] || g1<lmargin[i])  g1=x1;
			if(g2>rmargin[i] || g2<lmargin[i])  g2=x2;

			children[count_children][i]=g1;
			children[count_children+1][i]=g2;			
		}
		count_children+=2;
		if(count_children>=nchildren) break;
	}
	for(int i=0;i<nchildren;i++)
	{
		int status1=isFeasible[genome_count-i-1];
		int status2=problem->isFeasible(children[i]);
		if(!status2 && status1) continue;
		isFeasible[genome_count-i-1]=status2;
		genome[genome_count-i-1]=children[i];
	}
}

void	DoublePop::setElitism(int s)
{
	elitism = s;
}

void	DoublePop::mutate()
{
	int start = elitism * (int)(genome_count*selection_rate);
	start = elitism;
	for(int i=start;i<genome_count;i++)
	{
		if(isFeasible[i]) mutateFeasible(genome[i]);
		else
		for(int j=0;j<genome_size;j++)
		{
			double r=drand48();
			double old=genome[i][j];
			if(r<mutation_rate)
			{
				genome[i][j]=lmargin[j]+drand48()*(rmargin[j]-lmargin[j]);
			}
		}
	}
}

void	DoublePop::calcFitnessArray()
{
	countFeasible = 0;
	for(int i=0;i<genome_count;i++)
	{
		isFeasible[i]=problem->isFeasible(genome[i]);
		fitness[i]=evalFitness(genome[i]);
		
	/*	if(rand()%genome_count<=1) 
		{
			LocalSearch(genome[i],fitness[i]);
			fitness[i]=evalFitness(genome[i]);	
		}
	*/	
		countFeasible+=isFeasible[i];
	}
}

int	DoublePop::getGeneration() const
{
	return generation;
}

int	DoublePop::getCount() const
{
	return genome_count;
}

int	DoublePop::getSize() const
{
	return genome_size;
}

void	DoublePop::nextGeneration()
{
    extern int      localsearch_generations;
    extern int      localsearch_chromosomes;
	if(generation) mutate();
	calcFitnessArray();
        if((generation+1)%localsearch_generations==0)
    	{
            for(int i=0;i<localsearch_chromosomes;i++)
			localSearch(rand() % genome_count);
    }
	select();
    	crossover();
	++generation;
}

void	DoublePop::setMutationRate(double r)
{
	if(r>=0 && r<=1) mutation_rate = r;
}

void	DoublePop::setSelectionRate(double r)
{
	if(r>=0 && r<=1) selection_rate = r;
}

double	DoublePop::getSelectionRate() const
{
	return selection_rate;
}

double	DoublePop::getMutationRate() const
{
	return mutation_rate;
}

double	DoublePop::getBestFitness() const
{
	return  fitness[0];
}

Data	DoublePop::getBestGenome() const
{
	vector<double> g;g.resize(genome_size);
	for(int i=0;i<genome_size;i++) g[i]=genome[0][i];
	return g;
}
double	DoublePop::evaluateBestFitness() 
{
	vector<double> g;g.resize(genome_size);
	for(int i=0;i<genome_size;i++) g[i]=genome[0][i];	
	return evalFitness(g);
}

void	DoublePop::setBest(Data x,double y)
{
	for(int i=0;i<genome_size;i++)
		genome[0][i]=x[i];
	fitness[0]=y;
}

double	DoublePop::getFitness(int pos) const
{
	return fitness[pos];
}

double	DoublePop::estimateVariance()
{
	double x1=0.0;
	double x2=0.0;
	for(int i=0;i<genome_count;i++)
	{
		x1+=fitness[i];
		x2+=fitness[i]*fitness[i];
	}
	return sqrt(x2/genome_count-(x1/genome_count)*(x1/genome_count));
}

void	DoublePop::localSearch(int gpos)
{
	int pos=gpos;
	Data itemp;
	itemp.resize(genome_size);
	for(int iters=1;iters<=100;iters++)
	{
        int randgenome=rand() % genome_count;
        int feas1=problem->isFeasible(genome[pos]);
        int feas2=problem->isFeasible(genome[randgenome]);
        if(feas1 && !feas2) continue;

		int randpos=rand() % genome_size;
		for(int i=0;i<randpos;i++) itemp[i]=genome[pos][i];
		for(int i=randpos;i<genome_size;i++) itemp[i]=genome[randgenome][i];
		double f=evalFitness(itemp);
		if(fabs(f)<fabs(fitness[pos]))
		{
			for(int i=0;i<genome_size;i++) genome[pos][i]=itemp[i];
			fitness[pos]=f;
		}
		else
		{
			for(int i=0;i<randpos;i++) itemp[i]=genome[randgenome][i];
			for(int i=randpos;i<genome_size;i++) itemp[i]=genome[pos][i];
			f=evalFitness(itemp);
			if(fabs(f)<fabs(fitness[pos]))
			{
				for(int i=0;i<genome_size;i++) genome[pos][i]=itemp[i];
				fitness[pos]=f;
			}
		}
	}
}

void	DoublePop::Solve()
{
	extern int generations;
	const int maxgenerations=generations;
	double x1=0.0;
	double x2=0.0;
	double stopat=0.0;
	double old_best=1e+100;
	extern int printLevel;
	int iprint=printLevel;
	int it=0;
	Data xbest;
	xbest.resize(problem->getDimension());
	double ybest;
	for(int i=2;i<=maxgenerations+1;i++)
	{
		nextGeneration();
		
		double fmin,fmax;
		fmin=fitness[0];
		double v=(1.0+fabs(fmin));
		x1+=v;
		x2+=v * v;
		double variance=x2/i-(x1/i)*(x1/i);
		if(fabs(fmin-old_best)>1e-5) 
		{
			it++;
			old_best = fmin;
			stopat=variance/2.0;
		}
#ifndef BATCH
		char ss[10];
		if(problem->isFeasible(genome[0]))
			strcpy(ss,"YES");
		else
			strcpy(ss,"NO");
		if(iprint)printf(" GENERATION:%4d\tVALUE=%15.5lg\tFEASIBLE=%s\n",
				i,fitness[0],ss);
#endif
		int status=isFeasible[0];
		
		if(!status &&  i<maxgenerations) continue;
		if(variance<stopat && i>=10) break;
    		extern int      localsearch_generations;
		if(i%localsearch_generations==0)
		{
			REND = 1e-7;
			LocalSearch(genome[0],fitness[0]);
			REND = 1e-5;
			fitness[0]=problem->funmin(genome[0]);
		}
	}

	REND = 1e-7;
	LocalSearch(genome[0],fitness[0]);
	REND = 1e-5;
	fitness[0]=problem->funmin(genome[0]);

	for(int i=0;i<genome[0].size();i++)
		printf("%lf ",genome[0][i]);
	printf(" (%20.10lf)\n",fitness[0]);

	problem->done(genome[0]);

}

DoublePop::~DoublePop()
{
}
