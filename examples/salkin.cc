# include <math.h>

extern "C"
{

/*	Return the dimension of the objective function.
 * */
int	getdimension()
{
	return 5;
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
	return 4;
}

/*	Return the left bounds of the objective function.
 * */
void	getleftmargin(double *x)
{
	x[0]=1;x[1]=80;x[2]=30;x[3]=145;x[4]=0;
}

/*	Return the right bounds of the objective function.
 * */
void	getrightmargin(double *x)
{
	x[0]=4;x[1]=88;x[2]=35;x[3]=150;x[4]=2;
}

/*	Return the objective function evaluated at a feasible point.
 * */
double	funmin(double *x)
{
	return -(3*x[0]+x[1]+2*x[2]+x[3]-x[4]);
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
	ineq[0]=25*x[0]-40*x[1]+16*x[2]+21*x[3]+x[4]-300.0;
	ineq[1]=x[0]+20*x[1]-50*x[2]+x[3]-x[4]-200.0;
	ineq[2]=60*x[0]+x[1]-x[2]+2*x[3]+x[4]-600.0;
	ineq[3]=-7*x[0]+4*x[1]+15*x[2]-x[3]+65*x[4]-700.0;
}

void	done(double *x)
{
}
}
