# include <math.h>

extern "C"
{

/*	Return the dimension of the objective function.
 * */
int	getdimension()
{
	return 2;
}

/*
 * Return the number of equalities of the objective problem.
 */
int	geteq()
{
	return 0;
}

/*	Return the number of inequalities of the objective problem.
 * */
int	getineq()
{
	return 2;
}

/*	Return the left bounds of the objective function.
 * */
void	getleftmargin(double *x)
{
	x[0]=x[1]=0;
}

/*	Return the right bounds of the objective function.
 * */
void	getrightmargin(double *x)
{
	x[0]=x[1]=10;
}

/*	Return the objective function evaluated at a feasible point.
 * */
double	funmin(double *x)
{
	double x1=x[0],x2=x[1];
	double v=sin(2.0*M_PI*x1)*sin(2.0*M_PI*x1)*sin(2.0*M_PI*x1)*sin(2.0*M_PI*x2)/(x1*x1*x1*(x1+x2));
	return -v;
}

/*	Store in the array eq the values 0 and 1 for the equalities of the problem.
 * */
void	feq(double *x,double *eq)
{
}

/*	Store in the array eq the values 0 and 1 for the inequalities of the problem.
 * */
void	fineq(double *x,double *ineq)
{
	double x1=x[0],x2=x[1];
	ineq[0]=-(-x1*x1+x2-1.0);
	ineq[1]=-(-1.0+x1-(x2-4)*(x2-4));
}

void	done(double *x)
{
}
}
