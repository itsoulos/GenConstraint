# ifndef __GETOPTIONS__H
# define __GETOPTIONS__H

extern 	int	chromosomes;
extern 	int	generations;
extern  double	selection_rate;
extern	double	mutation_rate;
extern  int 	localsearch_generations;
extern	int 	localsearch_chromosomes;
extern 	int 	printLevel;
extern 	int 	random_seed;
void		print_usage();
void		parse_cmd_line(int argc,char **argv);
# endif
