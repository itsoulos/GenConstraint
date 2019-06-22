/* COBYLA.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/
# define MAX480 100
# include <cobyla.h>
# include <math.h>
# include <tolmin.h>
#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Common Block Declarations */

struct {
    integer nprob;
} _BLNK__;

#define _BLNK__1 _BLNK__

/* Table of constant values */

static integer c__1 = 1;

int iter=0;

int cobyla_(integer *n, integer *m, doublereal *x, 
	doublereal *rhobeg, doublereal *rhoend, integer *iprint, integer *
	maxfun, doublereal *w, integer *iact,ConstraintInfo &Info);

void	cobyla(ConstraintInfo &Info)
{
	//initialize the arguments
	integer n=Info.xpoint.size();
	integer m=Info.ineq.size()+2*Info.eq.size()+Info.lmargin.size()+Info.rmargin.size();
	double *x=new double[n];
	Data xold;
	xold.resize(Info.xpoint.size());
	xold=Info.xpoint;
	for(int i=0;i<n;i++) x[i]=Info.xpoint[i];
	double rhobeg=0.5;
	extern double REND;
	double rhoend=0.0001;
	integer iprint=0;
	integer maxfun=2001;
	int memsize=n*(3*n+2*m+11)+4*m+6;
	double *w=new double[n*(3*n+2*m+11)+4*m+6];
	integer *iacct=new long[m+1];
	//call the minimizer
	iter=0;
	cobyla_(&n,&m,x,&rhobeg,&rhoend,&iprint,&maxfun,w,iacct,Info);	
	for(int i=0;i<n;i++) Info.xpoint[i]=x[i];
	
		
	if(!Info.problem->isFeasible(Info.xpoint)) 
	{
		Info.xpoint = xold;
	}
	
	Info.ypoint=Info.problem->funmin(Info.xpoint);
	//free the memory
	delete[] w;
	delete[] iacct;
	delete[] x;
}


/* Subroutine */ int cobyla_(integer *n, integer *m, doublereal *x, 
	doublereal *rhobeg, doublereal *rhoend, integer *iprint, integer *
	maxfun, doublereal *w, integer *iact,ConstraintInfo &Info)
{
    integer ia, idx, mpp, icon, isim, isigb, idatm, iveta, isimi, ivsig, 
	    iwork;
    extern /* Subroutine */ int cobylb_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *,ConstraintInfo &);


/*     This subroutine minimizes an objective function F(X) subject to M */
/*     inequality constraints on X, where X is a vector of variables that has */
/*     N components. The algorithm employs linear approximations to the */
/*     objective and constraint functions, the approximations being formed by */
/*     linear interpolation at N+1 points in the space of the variables. */
/*     We regard these interpolation points as vertices of a simplex. The */
/*     parameter RHO controls the size of the simplex and it is reduced */
/*     automatically from RHOBEG to RHOEND. For each RHO the subroutine tries */
/*     to achieve a good vector of variables for the current size, and then */
/*     RHO is reduced until the value RHOEND is reached. Therefore RHOBEG and */
/*     RHOEND should be set to reasonable initial changes to and the required */
/*     accuracy in the variables respectively, but this accuracy should be */
/*     viewed as a subject for experimentation because it is not guaranteed. */
/*     The subroutine has an advantage over many of its competitors, however, */
/*     which is that it treats each constraint individually when calculating */
/*     a change to the variables, instead of lumping the constraints together */
/*     into a single penalty function. The name of the subroutine is derived */
/*     from the phrase Constrained Optimization BY Linear Approximations. */

/*     The user must set the values of N, M, RHOBEG and RHOEND, and must */
/*     provide an initial vector of variables in X. Further, the value of */
/*     IPRINT should be set to 0, 1, 2 or 3, which controls the amount of */
/*     printing during the calculation. Specifically, there is no output if */
/*     IPRINT=0 and there is output only at the end of the calculation if */
/*     IPRINT=1. Otherwise each new value of RHO and SIGMA is printed. */
/*     Further, the vector of variables and some function information are */
/*     given either when RHO is reduced or when each new value of F(X) is */
/*     computed in the cases IPRINT=2 or IPRINT=3 respectively. Here SIGMA */
/*     is a penalty parameter, it being assumed that a change to X is an */
/*     improvement if it reduces the merit function */
/*                F(X)+SIGMA*MAX(0.0,-C1(X),-C2(X),...,-CM(X)), */
/*     where C1,C2,...,CM denote the constraint functions that should become */
/*     nonnegative eventually, at least to the precision of RHOEND. In the */
/*     printed output the displayed term that is multiplied by SIGMA is */
/*     called MAXCV, which stands for 'MAXimum Constraint Violation'. The */
/*     argument MAXFUN is an integer variable that must be set by the user to a */
/*     limit on the number of calls of CALCFC, the purpose of this routine being */
/*     given below. The value of MAXFUN will be altered to the number of calls */
/*     of CALCFC that are made. The arguments W and IACT provide real and */
/*     integer arrays that are used as working space. Their lengths must be at */
/*     least N*(3*N+2*M+11)+4*M+6 and M+1 respectively. */

/*     In order to define the objective and constraint functions, we require */
/*     a subroutine that has the name and arguments */
/*                SUBROUTINE CALCFC (N,M,X,F,CON) */
/*                DIMENSION X(*),CON(*)  . */
/*     The values of N and M are fixed and have been defined already, while */
/*     X is now the current vector of variables. The subroutine should return */
/*     the objective and constraint functions at X in F and CON(1),CON(2), */
/*     ...,CON(M). Note that we are trying to adjust X so that F(X) is as */
/*     small as possible subject to the constraint functions being nonnegative. */

/*     Partition the working space array W to provide the storage that is needed */
/*     for the main calculation. */

    /* Parameter adjustments */
    --iact;
    --w;
    --x;

    /* Function Body */
    mpp = *m + 2;
    icon = 1;
    isim = icon + mpp;
    isimi = isim + *n * *n + *n;
    idatm = isimi + *n * *n;
    ia = idatm + *n * mpp + mpp;
    ivsig = ia + *m * *n + *n;
    iveta = ivsig + *n;
    isigb = iveta + *n;
    idx = isigb + *n;
    iwork = idx + *n;
    cobylb_(n, m, &mpp, &x[1], rhobeg, rhoend, iprint, maxfun, &w[icon], &w[
	    isim], &w[isimi], &w[idatm], &w[ia], &w[ivsig], &w[iveta], &w[
	    isigb], &w[idx], &w[iwork], &iact[1],Info);
    return 0;
} /* cobyla_ */

/* ------------------------------------------------------------------------------ */
/* Subroutine */ int cobylb_(integer *n, integer *m, integer *mpp, doublereal 
	*x, doublereal *rhobeg, doublereal *rhoend, integer *iprint, integer *
	maxfun, doublereal *con, doublereal *sim, doublereal *simi, 
	doublereal *datmat, doublereal *a, doublereal *vsig, doublereal *veta,
	 doublereal *sigbar, doublereal *dx, doublereal *w, integer *iact,ConstraintInfo &Info)
{
    /* Format strings */
    static char fmt_10[] = "(/3x,\002The initial value of RHO is\002,1pe13.6\
,2x,\002and PARMU is set to zero.\002)";
    static char fmt_50[] = "(/3x,\002Return from subroutine COBYLA because t\
he \002,\002MAXFUN limit has been reached.\002)";
    static char fmt_70[] = "(/3x,\002NFVALS =\002,i5,3x,\002F =\002,1pe13.6,\
4x,\002MAXCV =\002,1pe13.6/3x,\002X =\002,1pe13.6,1p4e15.6)";
    static char fmt_80[] = "(1pe19.6,1p4e15.6)";
    static char fmt_210[] = "(/3x,\002Return from subroutine COBYLA because\
 \002,\002rounding errors are becoming damaging.\002)";
    static char fmt_410[] = "(/3x,\002Increase in PARMU to\002,1pe13.6)";
    static char fmt_580[] = "(/3x,\002Reduction in RHO to\002,1pe13.6,\002  \
and PARMU =\002,1pe13.6)";
    static char fmt_590[] = "(/3x,\002Normal return from subroutine COBYL\
A\002)";

    /* System generated locals */
    integer sim_dim1, sim_offset, simi_dim1, simi_offset, datmat_dim1, 
	    datmat_offset, a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal f;
    integer i__, j, k, l, mp, np, iz;
    doublereal phi, rho, sum, beta, cmin, cmax;
    integer ivmc;
    doublereal weta;
    integer ivmd;
    doublereal temp, wsig, gamma;
    integer iflag;
    doublereal alpha, delta, denom, tempa, barmu;
    integer nbest, ifull, iptem, jdrop;
    doublereal ratio, vmold, parmu, error, vmnew;
    extern /* Subroutine */ int calcfc_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *,ConstraintInfo &);
    integer ibrnch;
    doublereal edgmax, pareta, prerec, phimin, parsig;
    integer isdirn, nfvals, izdota;
    doublereal cvmaxm, dxsign, prerem;
    integer iptemp;
    doublereal resmax, cvmaxp;
    integer idxnew;
    doublereal resnew, trured;
    extern /* Subroutine */ int trstlp_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *,ConstraintInfo &);

    /* Fortran I/O blocks */
    static cilist io___22 = { 0, 6, 0, fmt_10, 0 };
    static cilist io___29 = { 0, 6, 0, fmt_50, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_70, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_80, 0 };
    static cilist io___39 = { 0, 6, 0, fmt_210, 0 };
    static cilist io___59 = { 0, 6, 0, fmt_410, 0 };
    static cilist io___71 = { 0, 6, 0, fmt_580, 0 };
    static cilist io___72 = { 0, 6, 0, fmt_70, 0 };
    static cilist io___73 = { 0, 6, 0, fmt_80, 0 };
    static cilist io___74 = { 0, 6, 0, fmt_590, 0 };
    static cilist io___75 = { 0, 6, 0, fmt_70, 0 };
    static cilist io___76 = { 0, 6, 0, fmt_80, 0 };



/*     Set the initial values of some parameters. The last column of SIM holds */
/*     the optimal vertex of the current simplex, and the preceding N columns */
/*     hold the displacements from the optimal vertex to the other vertices. */
/*     Further, SIMI holds the inverse of the matrix that is contained in the */
/*     first N columns of SIM. */

    /* Parameter adjustments */
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    simi_dim1 = *n;
    simi_offset = 1 + simi_dim1;
    simi -= simi_offset;
    sim_dim1 = *n;
    sim_offset = 1 + sim_dim1;
    sim -= sim_offset;
    datmat_dim1 = *mpp;
    datmat_offset = 1 + datmat_dim1;
    datmat -= datmat_offset;
    --x;
    --con;
    --vsig;
    --veta;
    --sigbar;
    --dx;
    --w;
    --iact;

    /* Function Body */
    iptem = min(*n,5);
    iptemp = iptem + 1;
    np = *n + 1;
    mp = *m + 1;
    alpha = .25;
    beta = 2.1;
    gamma = .5;
    delta = 1.1;
    rho = *rhobeg;
    parmu = 0.;
    if (*iprint >= 2) {
    }
    nfvals = 0;
    temp = 1. / rho;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sim[i__ + np * sim_dim1] = x[i__];
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    sim[i__ + j * sim_dim1] = 0.;
/* L20: */
	    simi[i__ + j * simi_dim1] = 0.;
	}
	sim[i__ + i__ * sim_dim1] = rho;
/* L30: */
	simi[i__ + i__ * simi_dim1] = temp;
    }
    jdrop = np;
    ibrnch = 0;

/*     Make the next call of the user-supplied subroutine CALCFC. These */
/*     instructions are also used for calling CALCFC during the iterations of */
/*     the algorithm. */

L40:
    if (nfvals >= *maxfun && nfvals > 0) {
	if (*iprint >= 1) {
	}
	goto L600;
    }
    ++nfvals;
    calcfc_(n, m, &x[1], &f, &con[1],Info);
    resmax = 0.;
    if (*m > 0) {
	i__1 = *m;
	for (k = 1; k <= i__1; ++k) {
/* L60: */
/* Computing MAX */
	    d__1 = resmax, d__2 = -con[k];
	    resmax = max(d__1,d__2);
	}
    }
    if (nfvals == *iprint - 1 || *iprint == 3) {
	i__1 = iptem;
	for (i__ = 1; i__ <= i__1; ++i__) {
	}
	if (iptem < *n) {
	    i__1 = *n;
	    for (i__ = iptemp; i__ <= i__1; ++i__) {
	    }
	}
    }
    con[mp] = f;
    con[*mpp] = resmax;
    if (ibrnch == 1) {
	goto L440;
    }

/*     Set the recently calculated function values in a column of DATMAT. This */
/*     array has a column for each vertex of the current simplex, the entries of */
/*     each column being the values of the constraint functions (if any) */
/*     followed by the objective function and the greatest constraint violation */
/*     at the vertex. */

    i__1 = *mpp;
    for (k = 1; k <= i__1; ++k) {
/* L90: */
	datmat[k + jdrop * datmat_dim1] = con[k];
    }
    if (nfvals > np) {
	goto L130;
    }

/*     Exchange the new vertex of the initial simplex with the optimal vertex if */
/*     necessary. Then, if the initial simplex is not complete, pick its next */
/*     vertex and calculate the function values there. */

    if (jdrop <= *n) {
	if (datmat[mp + np * datmat_dim1] <= f) {
	    x[jdrop] = sim[jdrop + np * sim_dim1];
	} else {
	    sim[jdrop + np * sim_dim1] = x[jdrop];
	    i__1 = *mpp;
	    for (k = 1; k <= i__1; ++k) {
		datmat[k + jdrop * datmat_dim1] = datmat[k + np * datmat_dim1]
			;
/* L100: */
		datmat[k + np * datmat_dim1] = con[k];
	    }
	    i__1 = jdrop;
	    for (k = 1; k <= i__1; ++k) {
		sim[jdrop + k * sim_dim1] = -rho;
		temp = 0.;
		i__2 = jdrop;
		for (i__ = k; i__ <= i__2; ++i__) {
/* L110: */
		    temp -= simi[i__ + k * simi_dim1];
		}
/* L120: */
		simi[jdrop + k * simi_dim1] = temp;
	    }
	}
    }
    if (nfvals <= *n) {
	jdrop = nfvals;
	x[jdrop] += rho;
	goto L40;
    }
L130:
    ibrnch = 1;

/*     Identify the optimal vertex of the current simplex. */

L140:
    phimin = datmat[mp + np * datmat_dim1] + parmu * datmat[*mpp + np * 
	    datmat_dim1];
    nbest = np;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	temp = datmat[mp + j * datmat_dim1] + parmu * datmat[*mpp + j * 
		datmat_dim1];
	if (temp < phimin) {
	    nbest = j;
	    phimin = temp;
	} else if (temp == phimin && parmu == 0.) {
	    if (datmat[*mpp + j * datmat_dim1] < datmat[*mpp + nbest * 
		    datmat_dim1]) {
		nbest = j;
	    }
	}
/* L150: */
    }

/*     Switch the best vertex into pole position if it is not there already, */
/*     and also update SIM, SIMI and DATMAT. */

    if (nbest <= *n) {
	i__1 = *mpp;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    temp = datmat[i__ + np * datmat_dim1];
	    datmat[i__ + np * datmat_dim1] = datmat[i__ + nbest * datmat_dim1]
		    ;
/* L160: */
	    datmat[i__ + nbest * datmat_dim1] = temp;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    temp = sim[i__ + nbest * sim_dim1];
	    sim[i__ + nbest * sim_dim1] = 0.;
	    sim[i__ + np * sim_dim1] += temp;
	    tempa = 0.;
	    i__2 = *n;
	    for (k = 1; k <= i__2; ++k) {
		sim[i__ + k * sim_dim1] -= temp;
/* L170: */
		tempa -= simi[k + i__ * simi_dim1];
	    }
/* L180: */
	    simi[nbest + i__ * simi_dim1] = tempa;
	}
    }

/*     Make an error return if SIGI is a poor approximation to the inverse of */
/*     the leading N by N submatrix of SIG. */

    error = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    temp = 0.;
	    if (i__ == j) {
		temp += -1.;
	    }
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
/* L190: */
		temp += simi[i__ + k * simi_dim1] * sim[k + j * sim_dim1];
	    }
/* L200: */
/* Computing MAX */
	    d__1 = error, d__2 = fabs(temp);
	    error = max(d__1,d__2);
	}
    }
    if (error > .1) {
	if (*iprint >= 1) {
	}
	goto L600;
    }

/*     Calculate the coefficients of the linear approximations to the objective */
/*     and constraint functions, placing minus the objective function gradient */
/*     after the constraint gradients in the array A. The vector W is used for */
/*     working space. */

    i__2 = mp;
    for (k = 1; k <= i__2; ++k) {
	con[k] = -datmat[k + np * datmat_dim1];
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* L220: */
	    w[j] = datmat[k + j * datmat_dim1] + con[k];
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    temp = 0.;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
/* L230: */
		temp += w[j] * simi[j + i__ * simi_dim1];
	    }
	    if (k == mp) {
		temp = -temp;
	    }
/* L240: */
	    a[i__ + k * a_dim1] = temp;
	}
    }

/*     Calculate the values of sigma and eta, and set IFLAG=0 if the current */
/*     simplex is not acceptable. */

    iflag = 1;
    parsig = alpha * rho;
    pareta = beta * rho;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wsig = 0.;
	weta = 0.;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
	    d__1 = simi[j + i__ * simi_dim1];
	    wsig += d__1 * d__1;
/* L250: */
/* Computing 2nd power */
	    d__1 = sim[i__ + j * sim_dim1];
	    weta += d__1 * d__1;
	}
	vsig[j] = 1. / sqrt(wsig);
	veta[j] = sqrt(weta);
	if (vsig[j] < parsig || veta[j] > pareta) {
	    iflag = 0;
	}
/* L260: */
    }

/*     If a new vertex is needed to improve acceptability, then decide which */
/*     vertex to drop from the simplex. */

    if (ibrnch == 1 || iflag == 1) {
	goto L370;
    }
    jdrop = 0;
    temp = pareta;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (veta[j] > temp) {
	    jdrop = j;
	    temp = veta[j];
	}
/* L270: */
    }
    if (jdrop == 0) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (vsig[j] < temp) {
		jdrop = j;
		temp = vsig[j];
	    }
/* L280: */
	}
    }

/*     Calculate the step to the new vertex and its sign. */

    temp = gamma * rho * vsig[jdrop];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L290: */
	dx[i__] = temp * simi[jdrop + i__ * simi_dim1];
    }
    cvmaxp = 0.;
    cvmaxm = 0.;
    i__1 = mp;
    for (k = 1; k <= i__1; ++k) {
	sum = 0.;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L300: */
	    sum += a[i__ + k * a_dim1] * dx[i__];
	}
	if (k < mp) {
	    temp = datmat[k + np * datmat_dim1];
/* Computing MAX */
	    d__1 = cvmaxp, d__2 = -sum - temp;
	    cvmaxp = max(d__1,d__2);
/* Computing MAX */
	    d__1 = cvmaxm, d__2 = sum - temp;
	    cvmaxm = max(d__1,d__2);
	}
/* L310: */
    }
    dxsign = 1.;
    if (parmu * (cvmaxp - cvmaxm) > sum + sum) {
	dxsign = -1.;
    }

/*     Update the elements of SIM and SIMI, and set the next X. */

    temp = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dx[i__] = dxsign * dx[i__];
	sim[i__ + jdrop * sim_dim1] = dx[i__];
/* L320: */
	temp += simi[jdrop + i__ * simi_dim1] * dx[i__];
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L330: */
	simi[jdrop + i__ * simi_dim1] /= temp;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (j != jdrop) {
	    temp = 0.;
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L340: */
		temp += simi[j + i__ * simi_dim1] * dx[i__];
	    }
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L350: */
		simi[j + i__ * simi_dim1] -= temp * simi[jdrop + i__ * 
			simi_dim1];
	    }
	}
/* L360: */
	x[j] = sim[j + np * sim_dim1] + dx[j];
    }
    goto L40;

/*     Calculate DX=x(*)-x(0). Branch if the length of DX is less than 0.5*RHO. */

L370:
    iz = 1;
    izdota = iz + *n * *n;
    ivmc = izdota + *n;
    isdirn = ivmc + mp;
    idxnew = isdirn + *n;
    ivmd = idxnew + *n;
    trstlp_(n, m, &a[a_offset], &con[1], &rho, &dx[1], &ifull, &iact[1], &w[
	    iz], &w[izdota], &w[ivmc], &w[isdirn], &w[idxnew], &w[ivmd],Info);
    if (ifull == 0) {
	temp = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L380: */
/* Computing 2nd power */
	    d__1 = dx[i__];
	    temp += d__1 * d__1;
	}
	if (temp < rho * .25 * rho) {
	    ibrnch = 1;
	    goto L550;
	}
    }

/*     Predict the change to F and the new maximum constraint violation if the */
/*     variables are altered from x(0) to x(0)+DX. */

    resnew = 0.;
    con[mp] = 0.;
    i__1 = mp;
    for (k = 1; k <= i__1; ++k) {
	sum = con[k];
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L390: */
	    sum -= a[i__ + k * a_dim1] * dx[i__];
	}
	if (k < mp) {
	    resnew = max(resnew,sum);
	}
/* L400: */
    }

/*     Increase PARMU if necessary and branch back if this change alters the */
/*     optimal vertex. Otherwise PREREM and PREREC will be set to the predicted */
/*     reductions in the merit function and the maximum constraint violation */
/*     respectively. */

    barmu = 0.;
    prerec = datmat[*mpp + np * datmat_dim1] - resnew;
    if (prerec > 0.) {
	barmu = sum / prerec;
    }
    if (parmu < barmu * 1.5) {
	parmu = barmu * 2.;
	if (*iprint >= 2) {
	}
	phi = datmat[mp + np * datmat_dim1] + parmu * datmat[*mpp + np * 
		datmat_dim1];
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    temp = datmat[mp + j * datmat_dim1] + parmu * datmat[*mpp + j * 
		    datmat_dim1];
	    if (temp < phi) {
		goto L140;
	    }
	    if (temp == phi && parmu == 0.) {
		if (datmat[*mpp + j * datmat_dim1] < datmat[*mpp + np * 
			datmat_dim1]) {
		    goto L140;
		}
	    }
/* L420: */
	}
    }
    prerem = parmu * prerec - sum;

/*     Calculate the constraint and objective functions at x(*). Then find the */
/*     actual reduction in the merit function. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L430: */
	x[i__] = sim[i__ + np * sim_dim1] + dx[i__];
    }
    ibrnch = 1;
    goto L40;
L440:
    vmold = datmat[mp + np * datmat_dim1] + parmu * datmat[*mpp + np * 
	    datmat_dim1];
    vmnew = f + parmu * resmax;
    trured = vmold - vmnew;
    if (parmu == 0. && f == datmat[mp + np * datmat_dim1]) {
	prerem = prerec;
	trured = datmat[*mpp + np * datmat_dim1] - resmax;
    }

/*     Begin the operations that decide whether x(*) should replace one of the */
/*     vertices of the current simplex, the change being mandatory if TRURED is */
/*     positive. Firstly, JDROP is set to the index of the vertex that is to be */
/*     replaced. */

    ratio = 0.;
    if (trured <= 0.) {
	ratio = 1.;
    }
    jdrop = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	temp = 0.;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L450: */
	    temp += simi[j + i__ * simi_dim1] * dx[i__];
	}
	temp = fabs(temp);
	if (temp > ratio) {
	    jdrop = j;
	    ratio = temp;
	}
/* L460: */
	sigbar[j] = temp * vsig[j];
    }

/*     Calculate the value of ell. */

    edgmax = delta * rho;
    l = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (sigbar[j] >= parsig || sigbar[j] >= vsig[j]) {
	    temp = veta[j];
	    if (trured > 0.) {
		temp = 0.;
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
/* L470: */
/* Computing 2nd power */
		    d__1 = dx[i__] - sim[i__ + j * sim_dim1];
		    temp += d__1 * d__1;
		}
		temp = sqrt(temp);
	    }
	    if (temp > edgmax) {
		l = j;
		edgmax = temp;
	    }
	}
/* L480: */
    }
    if (l > 0) {
	jdrop = l;
    }
    if (jdrop == 0) {
	goto L550;
    }

/*     Revise the simplex by updating the elements of SIM, SIMI and DATMAT. */

    temp = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sim[i__ + jdrop * sim_dim1] = dx[i__];
/* L490: */
	temp += simi[jdrop + i__ * simi_dim1] * dx[i__];
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L500: */
	simi[jdrop + i__ * simi_dim1] /= temp;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (j != jdrop) {
	    temp = 0.;
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L510: */
		temp += simi[j + i__ * simi_dim1] * dx[i__];
	    }
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L520: */
		simi[j + i__ * simi_dim1] -= temp * simi[jdrop + i__ * 
			simi_dim1];
	    }
	}
/* L530: */
    }
    i__1 = *mpp;
    for (k = 1; k <= i__1; ++k) {
/* L540: */
	datmat[k + jdrop * datmat_dim1] = con[k];
    }

/*     Branch back for further iterations with the current RHO. */

    if (trured > 0. && trured >= prerem * .1) {
	goto L140;
    }
L550:
    if (iflag == 0) {
	ibrnch = 0;
	goto L140;
    }

/*     Otherwise reduce RHO if it is not at its least value and reset PARMU. */

    if (rho > *rhoend) {
	rho *= .5;
	if (rho <= *rhoend * 1.5) {
	    rho = *rhoend;
	}
	if (parmu > 0.) {
	    denom = 0.;
	    i__1 = mp;
	    for (k = 1; k <= i__1; ++k) {
		cmin = datmat[k + np * datmat_dim1];
		cmax = cmin;
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MIN */
		    d__1 = cmin, d__2 = datmat[k + i__ * datmat_dim1];
		    cmin = min(d__1,d__2);
/* L560: */
/* Computing MAX */
		    d__1 = cmax, d__2 = datmat[k + i__ * datmat_dim1];
		    cmax = max(d__1,d__2);
		}
		if (k <= *m && cmin < cmax * .5) {
		    temp = max(cmax,0.) - cmin;
		    if (denom <= 0.) {
			denom = temp;
		    } else {
			denom = min(denom,temp);
		    }
		}
/* L570: */
	    }
	    if (denom == 0.) {
		parmu = 0.;
	    } else if (cmax - cmin < parmu * denom) {
		parmu = (cmax - cmin) / denom;
	    }
	}
	if (*iprint >= 2) {
	}
	if (*iprint == 2) {
	    i__1 = iptem;
	    for (i__ = 1; i__ <= i__1; ++i__) {
	    }
	    if (iptem < *n) {
		i__1 = *n;
		for (i__ = iptemp; i__ <= i__1; ++i__) {
			    ;
		}
	    }
	}
	goto L140;
    }

/*     Return the best calculated values of the variables. */

    if (*iprint >= 1) {
    }
    if (ifull == 1) {
	goto L620;
    }
L600:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L610: */
	x[i__] = sim[i__ + np * sim_dim1];
    }
    f = datmat[mp + np * datmat_dim1];
    resmax = datmat[*mpp + np * datmat_dim1];
L620:
    if (*iprint >= 1) {
	i__1 = iptem;
	for (i__ = 1; i__ <= i__1; ++i__) {
	}
	if (iptem < *n) {
	    i__1 = *n;
	    for (i__ = iptemp; i__ <= i__1; ++i__) {
	    }
	}
    }
    *maxfun = nfvals;
    return 0;
} /* cobylb_ */

/* ------------------------------------------------------------------------------ */
/* Subroutine */ int trstlp_(integer *n, integer *m, doublereal *a, 
	doublereal *b, doublereal *rho, doublereal *dx, integer *ifull, 
	integer *iact, doublereal *z__, doublereal *zdota, doublereal *vmultc,
	 doublereal *sdirn, doublereal *dxnew, doublereal *vmultd,ConstraintInfo &Info)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, j, k;
    doublereal dd;
    integer kk;
    doublereal sd;
    integer kl, kp, kw;
    doublereal sp, ss, sum, tot, acca, accb, beta;
    integer nact, icon, mcon;
    doublereal temp, step;
    integer iout;
    doublereal alpha, tempa;
    integer isave;
    doublereal spabs;
    integer nactx;
    doublereal ratio, vsave, zdotv, zdotw, resold, zdvabs, zdwabs, sumabs, 
	    resmax, optold;
    integer icount;
    doublereal optnew, stpful;


/*     This subroutine calculates an N-component vector DX by applying the */
/*     following two stages. In the first stage, DX is set to the shortest */
/*     vector that minimizes the greatest violation of the constraints */
/*       A(1,K)*DX(1)+A(2,K)*DX(2)+...+A(N,K)*DX(N) .GE. B(K), K=2,3,...,M, */
/*     subject to the Euclidean length of DX being at most RHO. If its length is */
/*     strictly less than RHO, then we use the resultant freedom in DX to */
/*     minimize the objective function */
/*              -A(1,M+1)*DX(1)-A(2,M+1)*DX(2)-...-A(N,M+1)*DX(N) */
/*     subject to no increase in any greatest constraint violation. This */
/*     notation allows the gradient of the objective function to be regarded as */
/*     the gradient of a constraint. Therefore the two stages are distinguished */
/*     by MCON .EQ. M and MCON .GT. M respectively. It is possible that a */
/*     degeneracy may prevent DX from attaining the target length RHO. Then the */
/*     value IFULL=0 would be set, but usually IFULL=1 on return. */

/*     In general NACT is the number of constraints in the active set and */
/*     IACT(1),...,IACT(NACT) are their indices, while the remainder of IACT */
/*     contains a permutation of the remaining constraint indices. Further, Z is */
/*     an orthogonal matrix whose first NACT columns can be regarded as the */
/*     result of Gram-Schmidt applied to the active constraint gradients. For */
/*     J=1,2,...,NACT, the number ZDOTA(J) is the scalar product of the J-th */
/*     column of Z with the gradient of the J-th active constraint. DX is the */
/*     current vector of variables and here the residuals of the active */
/*     constraints should be zero. Further, the active constraints have */
/*     nonnegative Lagrange multipliers that are held at the beginning of */
/*     VMULTC. The remainder of this vector holds the residuals of the inactive */
/*     constraints at DX, the ordering of the components of VMULTC being in */
/*     agreement with the permutation of the indices of the constraints that is */
/*     in IACT. All these residuals are nonnegative, which is achieved by the */
/*     shift RESMAX that makes the least residual zero. */

/*     Initialize Z and some other variables. The value of RESMAX will be */
/*     appropriate to DX=0, while ICON will be the index of a most violated */
/*     constraint if RESMAX is positive. Usually during the first stage the */
/*     vector SDIRN gives a search direction that reduces all the active */
/*     constraint violations by one simultaneously. */

    /* Parameter adjustments */
	int count480 = 0;
    z_dim1 = *n;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --b;
    --dx;
    --iact;
    --zdota;
    --vmultc;
    --sdirn;
    --dxnew;
    --vmultd;

    /* Function Body */
    *ifull = 1;
    mcon = *m;
    nact = 0;
    resmax = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* L10: */
	    z__[i__ + j * z_dim1] = 0.;
	}
	z__[i__ + i__ * z_dim1] = 1.;
/* L20: */
	dx[i__] = 0.;
    }
    if (*m >= 1) {
	i__1 = *m;
	for (k = 1; k <= i__1; ++k) {
	    if (b[k] > resmax) {
		resmax = b[k];
		icon = k;
	    }
/* L30: */
	}
	i__1 = *m;
	for (k = 1; k <= i__1; ++k) {
	    iact[k] = k;
/* L40: */
	    vmultc[k] = resmax - b[k];
	}
    }
    if (resmax == 0.) {
	goto L480;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L50: */
	sdirn[i__] = 0.;
    }

/*     End the current stage of the calculation if 3 consecutive iterations */
/*     have either failed to reduce the best calculated value of the objective */
/*     function or to increase the number of active constraints since the best */
/*     value was calculated. This strategy prevents cycling, but there is a */
/*     remote possibility that it will cause premature termination. */

L60:
    optold = 0.;
    icount = 0;
L70:
    if (mcon == *m) {
	optnew = resmax;
    } else {
	optnew = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L80: */
	    optnew -= dx[i__] * a[i__ + mcon * a_dim1];
	}
    }
    if (icount == 0 || optnew < optold) {
	optold = optnew;
	nactx = nact;
	icount = 3;
    } else if (nact > nactx) {
	nactx = nact;
	icount = 3;
    } else {
	--icount;
	if (icount == 0) {
	    goto L490;
	}
    }

/*     If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to */
/*     the active set. Apply Givens rotations so that the last N-NACT-1 columns */
/*     of Z are orthogonal to the gradient of the new constraint, a scalar */
/*     product being set to zero if its nonzero value could be due to computer */
/*     rounding errors. The array DXNEW is used for working space. */

    if (icon <= nact) {
	goto L260;
    }
    kk = iact[icon];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L90: */
	dxnew[i__] = a[i__ + kk * a_dim1];
    }
    tot = 0.;
    k = *n;
L100:
    if (k > nact) {
	sp = 0.;
	spabs = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    temp = z__[i__ + k * z_dim1] * dxnew[i__];
	    sp += temp;
/* L110: */
	    spabs += fabs(temp);
	}
	acca = spabs + fabs(sp) * .1;
	accb = spabs + fabs(sp) * .2;
	if (spabs >= acca || acca >= accb) {
	    sp = 0.;
	}
	if (tot == 0.) {
	    tot = sp;
	} else {
	    kp = k + 1;
	    temp = sqrt(sp * sp + tot * tot);
	    alpha = sp / temp;
	    beta = tot / temp;
	    tot = temp;
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		temp = alpha * z__[i__ + k * z_dim1] + beta * z__[i__ + kp * 
			z_dim1];
		z__[i__ + kp * z_dim1] = alpha * z__[i__ + kp * z_dim1] - 
			beta * z__[i__ + k * z_dim1];
/* L120: */
		z__[i__ + k * z_dim1] = temp;
	    }
	}
	--k;
	goto L100;
    }

/*     Add the new constraint if this can be done without a deletion from the */
/*     active set. */

    if (tot != 0.) {
	++nact;
	zdota[nact] = tot;
	vmultc[icon] = vmultc[nact];
	vmultc[nact] = 0.;
	goto L210;
    }

/*     The next instruction is reached if a deletion has to be made from the */
/*     active set in order to make room for the new active constraint, because */
/*     the new constraint gradient is a linear combination of the gradients of */
/*     the old active constraints. Set the elements of VMULTD to the multipliers */
/*     of the linear combination. Further, set IOUT to the index of the */
/*     constraint to be deleted, but branch if no suitable index can be found. */

    ratio = -1.;
    k = nact;
L130:
    zdotv = 0.;
    zdvabs = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp = z__[i__ + k * z_dim1] * dxnew[i__];
	zdotv += temp;
/* L140: */
	zdvabs += fabs(temp);
    }
    acca = zdvabs + fabs(zdotv) * .1;
    accb = zdvabs + fabs(zdotv) * .2;
    if (zdvabs < acca && acca < accb) {
	temp = zdotv / zdota[k];
	if (temp > 0. && iact[k] <= *m) {
	    tempa = vmultc[k] / temp;
	    if (ratio < 0. || tempa < ratio) {
		ratio = tempa;
		iout = k;
	    }
	}
	if (k >= 2) {
	    kw = iact[k];
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L150: */
		dxnew[i__] -= temp * a[i__ + kw * a_dim1];
	    }
	}
	vmultd[k] = temp;
    } else {
	vmultd[k] = 0.;
    }
    --k;
    if (k > 0) {
	goto L130;
    }
    if (ratio < 0.) {
	goto L490;
    }

/*     Revise the Lagrange multipliers and reorder the active constraints so */
/*     that the one to be replaced is at the end of the list. Also calculate the */
/*     new value of ZDOTA(NACT) and branch if it is not acceptable. */

    i__1 = nact;
    for (k = 1; k <= i__1; ++k) {
/* L160: */
/* Computing MAX */
	d__1 = 0., d__2 = vmultc[k] - ratio * vmultd[k];
	vmultc[k] = max(d__1,d__2);
    }
    if (icon < nact) {
	isave = iact[icon];
	vsave = vmultc[icon];
	k = icon;
L170:
	kp = k + 1;
	kw = iact[kp];
	sp = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L180: */
	    sp += z__[i__ + k * z_dim1] * a[i__ + kw * a_dim1];
	}
/* Computing 2nd power */
	d__1 = zdota[kp];
	temp = sqrt(sp * sp + d__1 * d__1);
	alpha = zdota[kp] / temp;
	beta = sp / temp;
	zdota[kp] = alpha * zdota[k];
	zdota[k] = temp;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    temp = alpha * z__[i__ + kp * z_dim1] + beta * z__[i__ + k * 
		    z_dim1];
	    z__[i__ + kp * z_dim1] = alpha * z__[i__ + k * z_dim1] - beta * 
		    z__[i__ + kp * z_dim1];
/* L190: */
	    z__[i__ + k * z_dim1] = temp;
	}
	iact[k] = kw;
	vmultc[k] = vmultc[kp];
	k = kp;
	if (k < nact) {
	    goto L170;
	}
	iact[k] = isave;
	vmultc[k] = vsave;
    }
    temp = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L200: */
	temp += z__[i__ + nact * z_dim1] * a[i__ + kk * a_dim1];
    }
    if (temp == 0.) {
	goto L490;
    }
    zdota[nact] = temp;
    vmultc[icon] = 0.;
    vmultc[nact] = ratio;

/*     Update IACT and ensure that the objective function continues to be */
/*     treated as the last active constraint when MCON>M. */

L210:
    iact[icon] = iact[nact];
    iact[nact] = kk;
    if (mcon > *m && kk != mcon) {
	k = nact - 1;
	sp = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L220: */
	    sp += z__[i__ + k * z_dim1] * a[i__ + kk * a_dim1];
	}
/* Computing 2nd power */
	d__1 = zdota[nact];
	temp = sqrt(sp * sp + d__1 * d__1);
	alpha = zdota[nact] / temp;
	beta = sp / temp;
	zdota[nact] = alpha * zdota[k];
	zdota[k] = temp;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    temp = alpha * z__[i__ + nact * z_dim1] + beta * z__[i__ + k * 
		    z_dim1];
	    z__[i__ + nact * z_dim1] = alpha * z__[i__ + k * z_dim1] - beta * 
		    z__[i__ + nact * z_dim1];
/* L230: */
	    z__[i__ + k * z_dim1] = temp;
	}
	iact[nact] = iact[k];
	iact[k] = kk;
	temp = vmultc[k];
	vmultc[k] = vmultc[nact];
	vmultc[nact] = temp;
    }

/*     If stage one is in progress, then set SDIRN to the direction of the next */
/*     change to the current vector of variables. */

    if (mcon > *m) {
	goto L320;
    }
    kk = iact[nact];
    temp = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L240: */
	temp += sdirn[i__] * a[i__ + kk * a_dim1];
    }
    temp += -1.;
    temp /= zdota[nact];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L250: */
	sdirn[i__] -= temp * z__[i__ + nact * z_dim1];
    }
    goto L340;

/*     Delete the constraint that has the index IACT(ICON) from the active set. */

L260:
    if (icon < nact) {
	isave = iact[icon];
	vsave = vmultc[icon];
	k = icon;
L270:
	kp = k + 1;
	kk = iact[kp];
	sp = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L280: */
	    sp += z__[i__ + k * z_dim1] * a[i__ + kk * a_dim1];
	}
/* Computing 2nd power */
	d__1 = zdota[kp];
	temp = sqrt(sp * sp + d__1 * d__1);
	alpha = zdota[kp] / temp;
	beta = sp / temp;
	zdota[kp] = alpha * zdota[k];
	zdota[k] = temp;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    temp = alpha * z__[i__ + kp * z_dim1] + beta * z__[i__ + k * 
		    z_dim1];
	    z__[i__ + kp * z_dim1] = alpha * z__[i__ + k * z_dim1] - beta * 
		    z__[i__ + kp * z_dim1];
/* L290: */
	    z__[i__ + k * z_dim1] = temp;
	}
	iact[k] = kk;
	vmultc[k] = vmultc[kp];
	k = kp;
	if (k < nact) {
	    goto L270;
	}
	iact[k] = isave;
	vmultc[k] = vsave;
    }
    --nact;

/*     If stage one is in progress, then set SDIRN to the direction of the next */
/*     change to the current vector of variables. */

    if (mcon > *m) {
	goto L320;
    }
    temp = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L300: */
	temp += sdirn[i__] * z__[i__ + (nact + 1) * z_dim1];
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L310: */
	sdirn[i__] -= temp * z__[i__ + (nact + 1) * z_dim1];
    }
    goto L340;

/*     Pick the next search direction of stage two. */

L320:
    temp = 1. / zdota[nact];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L330: */
	sdirn[i__] = temp * z__[i__ + nact * z_dim1];
    }

/*     Calculate the step to the boundary of the trust region or take the step */
/*     that reduces RESMAX to zero. The two statements below that include the */
/*     factor 1.0E-6 prevent some harmless underflows that occurred in a test */
/*     calculation. Further, we skip the step if it could be zero within a */
/*     reasonable tolerance for computer rounding errors. */

L340:
    dd = *rho * *rho;
    sd = 0.;
    ss = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((d__1 = dx[i__], fabs(d__1)) >= *rho * 1e-6) {
/* Computing 2nd power */
	    d__2 = dx[i__];
	    dd -= d__2 * d__2;
	}
	sd += dx[i__] * sdirn[i__];
/* L350: */
/* Computing 2nd power */
	d__1 = sdirn[i__];
	ss += d__1 * d__1;
    }
    if (dd <= 0.) {
	goto L490;
    }
    temp = sqrt(ss * dd);
    if (fabs(sd) >= temp * 1e-6) {
	temp = sqrt(ss * dd + sd * sd);
    }
    stpful = dd / (temp + sd);
    step = stpful;
    if (mcon == *m) {
	acca = step + resmax * .1;
	accb = step + resmax * .2;
	if (step >= acca || acca >= accb) {
	    goto L480;
	}
	step = min(step,resmax);
    }

/*     Set DXNEW to the new variables if STEP is the steplength, and reduce */
/*     RESMAX to the corresponding maximum residual if stage one is being done. */
/*     Because DXNEW will be changed during the calculation of some Lagrange */
/*     multipliers, it will be restored to the following value later. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L360: */
	dxnew[i__] = dx[i__] + step * sdirn[i__];
    }
    if (mcon == *m) {
	resold = resmax;
	resmax = 0.;
	i__1 = nact;
	for (k = 1; k <= i__1; ++k) {
	    kk = iact[k];
	    temp = b[kk];
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L370: */
		temp -= a[i__ + kk * a_dim1] * dxnew[i__];
	    }
	    resmax = max(resmax,temp);
/* L380: */
	}
    }

/*     Set VMULTD to the VMULTC vector that would occur if DX became DXNEW. A */
/*     device is included to force VMULTD(K)=0.0 if deviations from this value */
/*     can be attributed to computer rounding errors. First calculate the new */
/*     Lagrange multipliers. */

    k = nact;
L390:
    zdotw = 0.;
    zdwabs = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp = z__[i__ + k * z_dim1] * dxnew[i__];
	zdotw += temp;
/* L400: */
	zdwabs += fabs(temp);
    }
    acca = zdwabs + fabs(zdotw) * .1;
    accb = zdwabs + fabs(zdotw) * .2;
    if (zdwabs >= acca || acca >= accb) {
	zdotw = 0.;
    }
    vmultd[k] = zdotw / zdota[k];
    if (k >= 2) {
	kk = iact[k];
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L410: */
	    dxnew[i__] -= vmultd[k] * a[i__ + kk * a_dim1];
	}
	--k;
	goto L390;
    }
    if (mcon > *m) {
/* Computing MAX */
	d__1 = 0., d__2 = vmultd[nact];
	vmultd[nact] = max(d__1,d__2);
    }

/*     Complete VMULTC by finding the new constraint residuals. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L420: */
	dxnew[i__] = dx[i__] + step * sdirn[i__];
    }
    if (mcon > nact) {
	kl = nact + 1;
	i__1 = mcon;
	for (k = kl; k <= i__1; ++k) {
	    kk = iact[k];
	    sum = resmax - b[kk];
	    sumabs = resmax + (d__1 = b[kk], fabs(d__1));
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		temp = a[i__ + kk * a_dim1] * dxnew[i__];
		sum += temp;
/* L430: */
		sumabs += fabs(temp);
	    }
	    acca = sumabs + fabs(sum) * .1;
	    accb = sumabs + fabs(sum) * .2;
	    if (sumabs >= acca || acca >= accb) {
		sum = 0.;
	    }
/* L440: */
	    vmultd[k] = sum;
	}
    }

/*     Calculate the fraction of the step from DX to DXNEW that will be taken. */

    ratio = 1.;
    icon = 0;
    i__1 = mcon;
    for (k = 1; k <= i__1; ++k) {
	if (vmultd[k] < 0.) {
	    temp = vmultc[k] / (vmultc[k] - vmultd[k]);
	    if (temp < ratio) {
		ratio = temp;
		icon = k;
	    }
	}
/* L450: */
    }

/*     Update DX, VMULTC and RESMAX. */

    temp = 1. - ratio;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L460: */
	dx[i__] = temp * dx[i__] + ratio * dxnew[i__];
    }
    i__1 = mcon;
    for (k = 1; k <= i__1; ++k) {
/* L470: */
/* Computing MAX */
	d__1 = 0., d__2 = temp * vmultc[k] + ratio * vmultd[k];
	vmultc[k] = max(d__1,d__2);
    }
    if (mcon == *m) {
	resmax = resold + ratio * (resmax - resold);
    }

/*     If the full step is not acceptable then begin another iteration. */
/*     Otherwise switch to stage two or end the calculation. */

    if (icon > 0) {
	goto L70;
    }
    if (step == stpful) {
	goto L500;
    }
L480:
	count480++;
	if(count480>=MAX480) return 0;
    mcon = *m + 1;
    icon = mcon;
    iact[mcon] = mcon;
    vmultc[mcon] = 0.;
    goto L60;

/*     We employ any freedom that may be available to reduce the objective */
/*     function before returning a DX whose length is less than RHO. */

L490:
    if (mcon == *m) {
	goto L480;
    }
    *ifull = 0;
L500:
    return 0;
} /* trstlp_ */


 int calcfc_(integer *n, integer *m, doublereal *x, 
	doublereal *f, doublereal *con,ConstraintInfo &Info)
{
	for(int i=0;i<*n;i++) Info.xpoint[i]=x[i];

	iter++;
	*f=Info.ypoint=Info.problem->funmin(Info.xpoint);
	
	MinInfo myinfo;
	myinfo.p=Info.problem;
	myinfo.iters=3;
	*f=tolmin(Info.xpoint,myinfo);
	for(int i=0;i<*n;i++) x[i]=Info.xpoint[i];

	for(int i=0;i<*n;i++) x[i]=Info.xpoint[i];
	Info.problem->fineq(Info.xpoint,Info.ineq);
	Info.problem->feq(Info.xpoint,Info.eq);

	int start=0;
	int end=Info.ineq.size();
	for(int i=start;i<end;i++) 
		con[i]=-Info.ineq[i-start];

	double sum=0.0;
	start=end;
	end=start+Info.eq.size();
	for(int i=start;i<end;i++) 
	{
		con[i]=Info.eq[i-start];
		sum+=con[i]*con[i];
	}

	start=end;
	end=start+Info.eq.size();
	for(int i=start;i<end;i++) 
	{
		con[i]=-Info.eq[i-start];
		sum+=con[i]*con[i];
	}


	start=end;
	end=start+Info.lmargin.size();
	for(int i=start;i<end;i++) con[i]=Info.xpoint[i-start]-Info.lmargin[i-start];//-MARGIN_EPS;

	start=end;
	end=start+Info.rmargin.size();
	for(int i=start;i<end;i++) con[i]=Info.rmargin[i-start]-Info.xpoint[i-start];//-MARGIN_EPS;

	for(int i=0;i<*m;i++) 
	{
		if(fabs(con[i])<1e-5) con[i]=1e-7;
	}

    return 0;
} 
#ifdef __cplusplus
	}
#endif
