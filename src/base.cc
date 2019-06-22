# include <base.h>
# include <math.h>
# include <cobyla.h>
# define LAMBDA	1e+3

static double dmax(double a,double b)
{
	return a>b?a:b;
}

static double getDistance(Data x1,Data x2)
{
	double s=0.0;
	for(int i=0;i<x1.size();i++) s+=(x1[i]-x2[i])*(x1[i]-x2[i]);
	return s;
}

Base::Base(Problem *p,int count)
{
	problem = p;
	genome.resize(count);
	fitness.resize(count);
	temp_eq.resize(problem->geteq());
	temp_ineq.resize(problem->getineq());
	for(int i=0;i<count;i++)
	{
		genome[i].resize(problem->getDimension());
		/*
		Data gg;
		gg.resize(problem->getDimension());
		double maxDist=0;
		for(int j=0;j<10;j++)
		{
			problem->getSample(genome[i]);
			double d=0.0;
			for(int k=0;k<i;k++) d+=getDistance(genome[i],genome[k]);
			d/=(i+1);
			if(d>maxDist)       {gg=genome[i];maxDist=d;}
		}
		genome[i]=gg;
		*/
		problem->getSample(genome[i]);
		fitness[i]=evalFitness(genome[i]);
	}
}
		

double	Base::evalFitness(Data &x)
{
	problem->feq(x,temp_eq);
	problem->fineq(x,temp_ineq);
	double v1=0.0,v2=0.0,v3=0.0;
	
	for(int i=0;i<temp_eq.size();i++)  
	{
		v2=v2+pow(temp_eq[i],2.0);
	}
	
	for(int i=0;i<temp_ineq.size();i++)
	{
		if(temp_ineq[i]>0) v3=v3+pow(temp_ineq[i],2.0);
	}
	v1=problem->funmin(x);
	double l=dmax(LAMBDA,1000 * fabs(v1));

	return (v1+l*(v2+v3));
}

void	Base::LocalSearch(Data &x,double &y)
{
	double v;
	Data gg;
	gg.resize(x.size());
	gg=x;
	ConstraintInfo  Info;
	Info.problem=problem;
	Info.xpoint.resize(gg.size());
	Info.xpoint = gg;
	Info.ypoint = y;
	//printf("ENTER %lf ",Info.ypoint);
	Info.lmargin.resize(gg.size());
	Info.rmargin.resize(gg.size());
	Info.lmargin=problem->getLeftMargin();
	Info.rmargin=problem->getRightMargin();
	Info.eq.resize(temp_eq.size());
	Info.eq=temp_eq;
	Info.ineq.resize(temp_ineq.size());
	Info.ineq = temp_ineq;
	cobyla(Info);
	gg = Info.xpoint;
	v  = Info.ypoint;
	//printf("LEAVE %lf\n",Info.ypoint);
	if(v<y)
	{
		x=gg;
		y=v;
	}
}

void	Base::mutateFeasible(Data &x)
{
	Data c;
	c.resize(x.size());
	c=x;
	Data lmargin,rmargin;
	lmargin.resize(x.size());
	rmargin.resize(x.size());
	lmargin=problem->getLeftMargin();
	rmargin=problem->getRightMargin();
	for(int i=0;i<x.size();i++)
	{
		if(drand48()<=0.05)
		{
			for(int k=0;k<20;k++)
			{
				double old=x[i];
				x[i]=lmargin[i]+(rmargin[i]-lmargin[i])*drand48();
				int status=problem->isFeasible(x);
				if(status) break; else x[i]=old;
			}
		}
		/*
		for(double p=0.05;p<=0.5;p+=0.05)
		{
			double old=x[i];
			int iters=0;
			do{
				x[i]=x[i]*(1+p*(2*drand48()-1));
				iters++;
				if(iters==10) break;
			}while(x[i]<lmargin[i] || x[i]>rmargin[i]);
			int status=problem->isFeasible(x);
			if(!status) x[i]=old; else break;
		}
		*/
	}
	int status=problem->isFeasible(x);
	if(!status) x=c;
}

void	Base::getMinimum(Data &x,double &y)
{
	x=genome[0];
	y=fitness[0];
}

Base::~Base()
{
}
