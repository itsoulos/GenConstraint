# include <stdio.h>
# include <unistd.h>
# include <stdlib.h>
# include <string.h>

int 	chromosomes=200;
int 	printLevel=1;
int 	random_seed=1;
int	generations=200;
double	selection_rate=0.10;
double	mutation_rate=0.05;
int 	localsearch_generations=50;
int 	localsearch_chromosomes=10;


void	print_usage()
{
	printf("\t-h	Display this help information.\n"
	       "\t-c integer	 Specify the size of the sample (default 200).\n"
	       "\t-m float	 Specify mutation rate (default 0.05).\n"
	       "\t-s float	 Specify selection rate (default 0.10).\n"
	       "\t-p {0|1}	 Specify the print level (0 or 1, default 1).\n"
	       "\t-r integer	 Specify the seed for the random numbers (default 1).\n"
	       "\t-d integer	 Specify number of local search chromosomes (default 10).\n"
	       "\t-n integer	 Specify number of local search generations. (Default 50).\n"
	       "\t-g integer	Specify the number of generations (default 200).\n");

}

void	parse_cmd_line(int argc,char **argv)
{
	const char *short_options="hc:p:r:g:n:d:s:m:";

	int next_option;
	do
	{
		next_option=getopt(argc,argv,short_options);
		switch(next_option)
		{
			case 'h':
				print_usage();
				exit(EXIT_SUCCESS);
				break;
			case 'g':
				generations=atoi(optarg);
				break;
			case 'c':
				chromosomes=atoi(optarg);
				break;
			case 's':
				selection_rate=atof(optarg);
				break;
			case 'm':
				mutation_rate=atof(optarg);
				break;
			case 'd':
				localsearch_chromosomes=atoi(optarg);
				break;
			case 'n':
				localsearch_generations=atoi(optarg);
				break;
			case 'p':
				printLevel=atoi(optarg);
				break;
			case 'r':
				random_seed=atoi(optarg);
				break;
			case -1:
				break;
			case '?':
				print_usage();
				exit(EXIT_FAILURE);
				break;
			default:
				print_usage();
				exit(EXIT_FAILURE);
				break;
		}
	}while(next_option!=-1);
}
