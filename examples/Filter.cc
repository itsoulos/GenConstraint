# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

extern "C"
{

class ConstraintProblem
{
protected:
    int N,M;
    double *Con;
public:
    ConstraintProblem();
    ConstraintProblem(int n,int m);
    virtual void getleftmargin(double *x)=0;
    virtual void getrightmargin(double *x)=0;
    virtual double funmin(double *x);
    virtual double granal(double *x,double *g);
    virtual double objfun(int n,int m, double *x, double *con)=0;
    int getN() const;
    int getM() const;
    void    setN(int n);
    void    setM(int m);
    virtual ~ConstraintProblem();
    bool isFeasible(double *x);
};

ConstraintProblem *thisProblem;

ConstraintProblem::ConstraintProblem()
{
    N=0;
    M=0;
    thisProblem=this;
    Con=NULL;
}



ConstraintProblem::ConstraintProblem(int n, int m)
{
    N=n;
    M=m;
    thisProblem=this;
    Con=new double[M];
}

void    ConstraintProblem::setN(int n)
{
    N=n;
}

void    ConstraintProblem::setM(int m)
{
    M=m;
    if(Con!=NULL) delete[] Con;
    Con=new double[M];
}

int ConstraintProblem::getN() const
{
    return N;
}

int ConstraintProblem::getM() const
{
    return M;
}


double  ConstraintProblem::funmin(double *x)
{
    double f=objfun(N,M,x,Con);
    double sum=0.0;
    for(int i=0;i<M;i++)
        sum+=Con[i]*Con[i];
    return f+100.0 * sum;
}

bool ConstraintProblem::isFeasible(double *x)
{
    objfun(N,M,x,Con);
    double rms=0.0;
    for(int i=0;i<M;i++)
        rms+=Con[i]*Con[i];
    rms=sqrt(rms);
    if(rms<1e-2) return true;
    return false;
}

double ConstraintProblem::granal(double *x,double *g)
{
    for(int i=0;i<N;i++)
        {
            double eps=pow(1e-18,1.0/3.0)*fmax(1.0,fabs(x[i]));
            x[i]+=eps;
            double v1=funmin(x);
            x[i]-=2.0 *eps;
            double v2=funmin(x);
            g[i]=(v1-v2)/(2.0 * eps);
            x[i]+=eps;
        }
}

ConstraintProblem::~ConstraintProblem()
{
    delete[] Con;
}

class Filter2D:public ConstraintProblem
{
private:
    double **A;
    double **Ccon,**Scon;
    int K,N1,N2;
    double H0,AR,AI,B1R,B1I,B2R,B2I;
    double *B,*C,*D;
    double *Omega1,*Omega2,*MValue,*Con,*MDValue;
public:
    Filter2D(int n,int m,int n1,int n2);
    void splitX(double *x);
    virtual void getleftmargin(double *x);
    virtual void getrightmargin(double *x);
    virtual double objfun(int n,int m, double *x, double *con);
    void savePlot(char *filename);
    int	getK();
    ~Filter2D();
};

Filter2D::Filter2D(int n,int m,int n1,int n2)
 :ConstraintProblem(n,m)
{
    N1=n1;
    N2=n2;

    n=n-1;//for H0
    K=(sqrt(16.0-4.0*(1.0-n))-4.0)/2.0;
    n=n-3*K;
    A=new double*[K+1];
    for(int i=0;i<K+1;i++)
        A[i]=new double[K+1];
    B=new double[K];
    C=new double[K];
    D=new double[K];
    Ccon=new double*[K+1];
    Scon=new double*[K+1];
    for(int i=0;i<K+1;i++)
    {
        Ccon[i]=new double[K+1];
        Scon[i]=new double[K+1];
    }
    Omega1=new double[2*(N1+1)*(N2+1)];
    Omega2=new double[2*(N1+1)*(N2+1)];
    MValue=new double[2*(N1+1)*(N2+1)];
    MDValue=new double[2*(N1+1)*(N2+1)];
    Con=new double[2*K];
}

int	Filter2D::getK()
{
	return K;
}

void Filter2D::getleftmargin(double *x)
{
    for(int i=0;i<N;i++) x[i]=-5.0;
    x[N-1]=0.001;
}

void Filter2D::getrightmargin(double *x)
{
    for(int i=0;i<N;i++) x[i]= 5.0;
    x[N-1]=0.1;
}

void    Filter2D::splitX(double *x)
{
    A[0][0]=1.0;
    H0=x[N-1];
    int icount=1;
    int rowCount=0;
    int colCount=1;
    int n=N-1-3*K;
    for(int i=1;i<n;i++)
    {       
       A[rowCount][colCount]=x[icount++];
       colCount++;
       if(colCount>K) {rowCount++;colCount=0;}

    }
    for(int i=0;i<K;i++)
    {
        B[i]=x[n+i];
        C[i]=x[n+K+i];
        D[i]=x[n+2*K+i];
    }
}

double  Filter2D::objfun(int n, int m, double *x, double *con)
{
   splitX(x);

    double sum=0.0;
    int p,q;
    int icount=0;
    n=n-3*K-1;


    for(int n1=0;n1<=N1;n1++)
    {
        for(int n2=0;n2<=N2;n2++)
        {
           double omega1=(n1 *1.0)/N1*M_PI;
           double omega2=(n2 *1.0)/N2*M_PI;
           omega1=-M_PI+(2*M_PI)*n1*1.0/N1;
           omega2=-M_PI+(2*M_PI)*n2*1.0/N2;
           Omega1[icount]=omega1;
           Omega2[icount]=omega2;
           for(p=0;p<=K;p++)
           {
               for(q=0;q<=K;q++)
               {
                   Ccon[p][q]=cos(p*omega1+q*omega2);
                   Scon[p][q]=sin(p*omega1+q*omega2);
               }
           }
           AR=A[0][0];
           AI=0.0;
           B1R=1.0;
           B1I=0.0;
           B2R=1.0;
           B2I=0.0;

           for(int l1=0;l1<=K;l1++)
           {

               for(int l2=0;l2<=K;l2++)
               {
                   if(l1>=1 && l2<K)
                   {
                       B2R+=B[l1]*Ccon[l1-1][l2];
                       B2I+=B[l1]*Scon[l1-1][l2]+C[l1]*Scon[l1-1][l2]+D[l1]*Scon[l2][l2];
                   }
                   if(l2<K && l1<K)
                   {
                       B1R+=B[l1]*Ccon[l1+1][l2];
                       B1I+=B[l1]*Scon[l1+1][l2];
                   }
                  if(l2>=1)
                  {
                      AI+=A[l1][l2]*Scon[l1][l2];
                      AR+=A[l1][l2]*Ccon[l1][l2];
                  }  
               }
           }
           double Mvalue=0.0,MDvalue=0.0;
           Mvalue=H0 * (sqrt(AR * AR+AI*AI))/
                   ((B1R*B1R+B1I*B1I)*(B2R*B2R+B2I*B2I));
           if(sqrt(omega1*omega1+omega2*omega2)<=0.12)
           {
               MDvalue=1.0;

           }
           else
           if(0.08*M_PI<=sqrt(omega1*omega1+omega2*omega2) &&
                   sqrt(omega1*omega1+omega2*omega2)<=0.12)
           {

               MDvalue=0.5;
           }
           else
               MDvalue=0.0;
         MDValue[icount]=MDvalue;
         MValue[icount++]=fabs(Mvalue);
           sum+=pow(fabs(Mvalue)-MDvalue,2.0);
        }
    }
    for(int i=0;i<K;i++)
    {
       con[i]=D[i]-fabs(B[i]+C[i])+1.0;
       con[K+i]=1.0-fabs(B[i]-C[i])-D[i];
    }

    return sum;
}


void    Filter2D::savePlot(char *filename)
{
    FILE *fp=fopen(filename,"w");
    int icount=0;
    for(int n1=0;n1<=N1;n1++)
    {
        for(int n2=0;n2<=N2;n2++)
        {
              fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",Omega1[icount],Omega2[icount],MValue[icount],MDValue[icount]);
              icount++;
        }
        fprintf(fp,"\n");
    }
}

Filter2D::~Filter2D()
{
    for(int i=0;i<K+1;i++)
    {
     delete[] A[i];
     delete[] Ccon[i];
     delete[] Scon[i];
    }
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] D;
    delete[] Ccon;
    delete[] Scon;
    delete[] Omega1;
    delete[] Omega2;
    delete[] MValue;
}

    Filter2D *filter=NULL;
    double *Con;
    //=new Filter2D(15,6,50,50);*/

/*	Return the dimension of the objective function.
 * */
int	getdimension()
{
	if(filter==NULL) 
	{
		filter=new Filter2D(15,6,50,50);
		Con=new double[2*filter->getK()];
	}
	return filter->getN();
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
	return 2*filter->getK();
}

/*	Return the left bounds of the objective function.
 * */
void	getleftmargin(double *x)
{
	filter->getleftmargin(x);
}

/*	Return the right bounds of the objective function.
 * */
void	getrightmargin(double *x)
{
	filter->getrightmargin(x);
}

/*	Return the objective function evaluated at a feasible point.
 * */
double	funmin(double *x)
{
	double f= filter->objfun(filter->getN(),filter->getM(),x,Con);
	return f;
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
	double f= filter->objfun(filter->getN(),filter->getM(),x,Con);
	for(int i=0;i<getineq();i++)
		ineq[i]=Con[i];

}

void	done(double *x)
{
	char filename[50];
	strcpy(filename,"cobyla.out");
    	filter->savePlot(filename);
	delete filter;
	delete[] Con;
}

}

