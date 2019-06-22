# ifndef __TOLMIN__H
# define __TOLMIN__H
# include <problem.h>

#ifdef __cplusplus
extern "C"{
# endif

typedef struct MinInfo
{
	Problem		*p;
	int		iters;
	int		fevals;
	int		gevals;
};
extern double	tolmin(Data &x,MinInfo &Info);

#ifdef __cplusplus
}
#endif

# endif
