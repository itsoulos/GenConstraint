# ifndef __DOUBLEPOP__H
# include <base.h>
using namespace std;
//# define PSO_MUTATION
typedef vector<double> Data;

class DoublePop: public Base
{
	private:
		vector<Data> 	children;
#ifdef PSO_MUTATION
		vector<Data>	best_genome;
		Data		best_fitness;
#endif
		Data		parent0,parent1;
		Data 		lmargin,rmargin;
		double		mutation_rate,selection_rate;
		int		genome_count;
		int		genome_size;
		int		generation;
		vector<int>	isFeasible;
		int		countFeasible;

		void	select();
		void	crossover();
		void	mutate();
		void	calcFitnessArray();
		int	elitism;
		double	estimateVariance();
		void	getTournamentElement(Data &x);
		void	tournament(Data &p1,Data &p2);
	public:
		DoublePop(int gcount,Problem *p);
		void	setElitism(int s);
		int	getGeneration() const;
		int	getCount() const;
		int	getSize() const;
		double	getFitness(int pos) const;
		void	nextGeneration();
		void	setMutationRate(double r);
		void	setSelectionRate(double r);
		double	getSelectionRate() const;
		double	getMutationRate() const;
		double	getBestFitness() const;
		double	evaluateBestFitness();
		Data	getBestGenome() const;
		void	reset();
		void	setBest(Data x,double y);
		void	localSearch(int x);
		virtual void	Solve();
		~DoublePop();
		
};
# define __DOUBLEPOP__H
# endif
