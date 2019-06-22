# include <math.h>

extern "C"
{

/*	Return the dimension of the objective function.
 * */
int	getdimension()
{
	return 6;
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
	return 6;
}

/*	Return the left bounds of the objective function.
 * */
void	getleftmargin(double *x)
{
	x[0]=0;
	x[1]=0;
	x[2]=1;
	x[3]=0;
	x[4]=1;
	x[5]=0;
}

/*	Return the right bounds of the objective function.
 * */
void	getrightmargin(double *x)
{
	x[0]=10;
	x[1]=5;
	x[2]=5;
	x[3]=6;
	x[4]=5;
	x[5]=10;
}

/*	Return the objective function evaluated at a feasible point.
 * */
double	funmin(double *x)
{
	return -25*(x[0]-2)*(x[0]-2)-(x[1]-2)*(x[1]-2)-(x[2]-1)*(x[2]-1)-
		(x[3]-4)*(x[3]-4)-(x[4]-1)*(x[4]-1)-(x[5]-4)*(x[5]-4);
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
	ineq[0]=-((x[2]-3)*(x[2]-3)+x[3]-4);
	ineq[1]=-((x[4]-3)*(x[4]-3)+x[5]-4);
	ineq[2]=x[0]-3*x[1]-2;
	ineq[3]=-x[0]+x[1]-2;
	ineq[4]=x[0]+x[1]-6;
	ineq[5]=-(x[0]+x[1]-2);
}

}
