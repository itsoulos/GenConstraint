# ifndef __PROBLEM__H
# define __PROBLEM__H
# include <stdlib.h>
# include <math.h>
# include <stdio.h>
# include <string.h>
# include <vector>
using namespace std;

# define	MARGIN_EPS	1e-5

typedef vector<double> Data;

/*	06-11-2004
 *	=================================================
 *	This class implements a generic problem for 
 *	optimization purposes. 
 *	=================================================
 * */
class Problem
{
	protected:
		/*	FIELDS DESCRIPTION:
		 *	=================================================
		 *	lmargin:	The left margins of the problem.
		 *	rmargin:	The right margins of the problem.
		 *	dimension:	The dimension of the problem.
		 *	neq:		The number of equalities.
		 *	nineq:		The number of inequalities.
		 *	=================================================
		 * */
		Data	lmargin;
		Data	rmargin;
		int     dimension;
		int	fevals,gevals,eqevals,neqevals;
		int	neq,nineq;
	public:
		/*	METHODS DESCRIPTION:
		 *	=============================================================
		 *	Problem(N,E,NE):	Sets only the dimension for the 
		 *			objective problem.
		 *	Problem(N,E,NE,l,r):	Set the dimension and the margins of
		 *			the objective problem.
		 *	setLeftMargin(l): Sets the left margin of the problem to l.
		 *	setRightMargin(r):Sets the right margin of the problem to r.
		 *	getLeftMargin():  Returns the left margin of the problem.
		 *	getRightMargin(): Returns the right margin of the problem.
		 *	getDimension()    Returns the dimension of the problem.
		 *	funmin(x):	Returns the function to be minimized.
		 *	granal(x,g):	Returns to g the granal of the function to 
		 *			be minimized.
		 *	feq(x,eq):	Store in the array eq the equalities of 
		 *			the problem (0,1).
		 *	fineq(x,ineq):	Store in the array ineq the inequalities 
		 *			of the problem (0,1).
		 *	feasible(x):	Return 1 if the point x is feasible.
		 *	~Problem():	The destructor of the proposed problem.
		 *	=============================================================
		 * */
		Problem(int N,int E,int NE);
		Problem(int N,int E,int NE,Data l,Data r);
		void	setLeftMargin(Data l);
		void	setRightMargin(Data r);
		int	getDimension()  const;
		Data 	getLeftMargin() 	const;
		Data 	getRightMargin() 	const;
		int	getGradientCriterion(Data x1,Data x2);
		int	getGradientCriterion(Data x1,Data x2,Data g1,Data g2);
		virtual double 	funmin(Data x);
		virtual void   	granal(Data x,Data &g);
		virtual void   	feq(Data x,Data &eq);
		virtual void	fineq(Data x,Data &ineq);
		virtual void	done(Data x);
		int	isFeasible(Data x);
		void	getSample(Data &x);
		int	getfevals();
		int	getgevals();
		int	geteqevals();
		int	getneqevals();
		int	geteq();
		int	getineq();
		~Problem();
};

# endif
