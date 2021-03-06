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
	x[0]=x[1]=6;
}

/*	Return the objective function evaluated at a feasible point.
 * */
double	funmin(double *x)
{
	double x1=x[0],x2=x[1];
	return (x1*x1+x2-11)*(x1*x1+x2-11)+(x1+x2*x2-7)*(x1+x2*x2-7);
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
	double x1=x[0],x2=x[1];
	ineq[0]=-(4.84-(x1-0.05)*(x1-0.05)-(x2-2.5)*(x2-2.5));
	ineq[1]=-(x1*x1+(x2-2.5)*(x2-2.5)-4.84);
}

}
