# include <problem.h>		
# include <stdlib.h>
# include <math.h>


Problem::Problem(int N,int E,int NE)
{
	fevals=gevals=eqevals=neqevals=0;
	neq=E;
	nineq=NE;
	dimension = N;
	lmargin.resize(0);
	rmargin.resize(0);
}

Problem::Problem(int N,int E,int NE,Data l,Data r)
{
	fevals=gevals=eqevals=neqevals=0;
	neq=E;
	nineq=NE;
	dimension = N;
	lmargin.resize(N);
	lmargin = l;
	rmargin.resize(N);
	rmargin = r;
}

int	Problem::geteqevals()
{
	return eqevals;
}

int	Problem::getneqevals()
{
	return neqevals;
}

int	Problem::getfevals()
{
	return fevals;
}

int	Problem::getgevals()
{
	return gevals;
}

void	Problem::setLeftMargin(Data l)
{
	lmargin = l;
}

void	Problem::setRightMargin(Data r)
{
	rmargin = r;
}

int	Problem::getDimension()		const
{
	return dimension;
}

Data 	Problem::getLeftMargin() 	const
{
	return lmargin;
}

Data 	Problem::getRightMargin() 	const
{
	return rmargin;
}

double 	Problem::funmin(Data x)
{
	return 0.0;
}

void	Problem::done(Data x)
{
}

void   	Problem::granal(Data x,Data &g)
{
}

void	Problem::feq(Data x,Data &eq)
{
}

void	Problem::fineq(Data x,Data &ineq)
{
}

int	Problem::geteq()
{
	return neq;
}

int	Problem::getineq()
{
	return nineq;
}

int	Problem::getGradientCriterion(Data x1,Data x2,Data g1,Data g2)
{
	double s=0.0;
	for(int i=0;i<x1.size();i++)
		s+=(x1[i]-x2[i])*(g1[i]-g2[i]);
	return s>=0;
}

int	Problem::getGradientCriterion(Data x1,Data x2)
{
	Data g1,g2;
	g1.resize(x1.size());
	g2.resize(x2.size());
	granal(x1,g1);
	granal(x2,g2);
	double s=0.0;
	for(int i=0;i<x1.size();i++)
		s+=(x1[i]-x2[i])*(g1[i]-g2[i]);
	return s>=0;
}

void	Problem::getSample(Data &x)
{
	if(lmargin.size())
	{
		for(int i=0;i<dimension;i++)
		{
		        x[i]=lmargin[i]+(rmargin[i]-lmargin[i])*(drand48());
		}

	}
	else
	{
		for(int i=0;i<dimension;i++)
			x[i]=2.0*drand48()-1.0;
	}
}

int	Problem::isFeasible(Data x)
{
	for(int i=0;i<x.size();i++)
	{
		//if(x[i]-rmargin[i]>MARGIN_EPS || lmargin[i]-x[i]>MARGIN_EPS) return 0;
		if(x[i]>rmargin[i] || lmargin[i]>x[i]) return 0;
	}
	if(neq)
	{
		Data k;
		k.resize(neq);
		feq(x,k);
		double eq_max=-1;
		double sum=0.0;
		for(int i=0;i<k.size();i++) 
		{
			if(fabs(k[i])>eq_max) eq_max=fabs(k[i]);
			sum+=k[i]*k[i];
		}
		if(sum>1e-5) return 0;
		//if(eq_max>MARGIN_EPS*10) return 0;
	}
	Data t;
	t.resize(nineq);
	fineq(x,t);
	if(t.size()) for(int i=0;i<t.size();i++) if(t[i]>MARGIN_EPS) {return 0;}
	return 1;
}


Problem::~Problem()
{
}
