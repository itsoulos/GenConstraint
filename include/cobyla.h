# ifndef __COBYLA__H
# include <problem.h>
#ifdef __cplusplus
extern "C" {
#endif
using namespace std;
typedef struct
{
	Problem *problem;
	Data    xpoint;
	Data    ineq;	
	Data	eq;
	Data	lmargin;
	Data 	rmargin;
	double  ypoint;
}ConstraintInfo;

void	cobyla(ConstraintInfo &Info);
#ifdef __cplusplus
	}
#endif

# define __COBYLA__H
# endif
