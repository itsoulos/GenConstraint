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
	return 4;
}

/*	Return the left bounds of the objective function.
 * */
void	getleftmargin(double *x)
{
	x[0]=x[1]=x[2]=x[3]=x[4]=x[5]=0.0;
}

/*	Return the right bounds of the objective function.
 * */
void	getrightmargin(double *x)
{
	x[0]=0.31;x[1]=0.046;x[2]=0.068;x[3]=0.042;x[4]=0.028;x[5]=0.0134;
}

/*	Return the objective function evaluated at a feasible point.
 * */
double	funmin(double *x)
{
	return 4.3*x[0]+31.8*x[1]+63.3*x[2]+15.8*x[3]+68.5*x[4]+4.7*x[5];
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
	double x1=x[0],x2=x[1],x3=x[2],x4=x[3],x5=x[4],x6=x[5];
	ineq[0]=-(17.1*x1+38.2*x2+204.2*x3+212.3*x4+623.4*x5+1495.5*x6-169.0*x1*x3-3580*x3*x5
		-3810*x4*x5-18500*x4*x6-24300*x5*x6-4.97);
	ineq[1]=-(1.88+17.9*x1+36.8*x2+113.9*x3+169.7*x4+337.8*x5+1385.2*x6-139*x1*x3-2450*x4*x5
		-600*x4*x6-17200*x5*x6);
	ineq[2]=-(429.08-273*x2-70*x4-819*x5+26000*x4*x5);
	ineq[3]=-(159.9*x1-311*x2+587*x4+391*x5+2198*x6-14000*x1*x6+78.02);
}

void	done(double *x)
{
}
}
