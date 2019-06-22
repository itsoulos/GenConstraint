# ifndef __BASE__H
# include <problem.h>
# include <vector>
using namespace std;

class	Base
{
	protected:
		vector<Data>	genome;
		Data		fitness;
		Problem		*problem;
		Data		temp_eq;
		Data		temp_ineq;
	public:
		Base(Problem *p,int count);
		double	evalFitness(Data &x);
		virtual void Solve()=0;
		void	LocalSearch(Data &x,double &y);
		void	getMinimum(Data &x,double &y);
		void	mutateFeasible(Data &x);
		~Base();	
};
# define __BASE__H
# endif
