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
	return 1;
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
	x[0]=x[1]=1;
}

/*	Return the objective function evaluated at a feasible point.
 * */
double	funmin(double *x)
{
	return -x[0]-x[1];
}


void done(double *x)
{
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
	double a=2.0;
	double b=0.25;
	ineq[0]=-(((x[0]-1)*(x[0]-1)+(x[1]-1))*(1/(2.0*a*a)-1.0/(2*b*b))+(x[0]-1)*(x[1]-1)*(1.0/(a*a)-1.0/(b*b))-1.0);
}

}
