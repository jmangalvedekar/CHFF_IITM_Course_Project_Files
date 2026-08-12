#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

//Temperature Boundary Conditions
const double Tx1 = 0; //left edge
const double Ty1 = 1; //bottom edge

//Defining the Domain
const int n = 20; //No. of elements/nodes
const double length = 3.0; //length of an edge of the square plate
const double gamma = 3.0; //Diffusion coefficient value
const double dx = length/n; //length and breadth of an element
const double dv = dx*dx; //volume of an element

//Function Declarations
void print_array(int m, double* arr); //Prints an array
void print_matrix(int p, int q, double** mat); //Prints a matrix
void init_d(int p, double** mat); //Initializes the diagonal elements
void init_b(int p, double** mat_T, double** mat_b, double** D, double** U, double** L); //Initializes the b-vector
void setup_guess(double guess, int p, double** mat_T); //Sets up the guess values for temperature
void tdma(int p, double* arr_b, double* arr_T, double** D, double** U, double** L); //Solves a tridiagonal matrix system
void print_Temp(int p, double** mat_T); //Prints the temperature matrix in an understandable form
void solve(int p, double** mat_T, double** mat_T_old, double** mat_b, double** D, double** U, double** L);
double calc_res(int p, double** T_old, double** T); //Calculates the residue of T wrt T_old
double u(double x);
double v(double y);
double Sc(int i, int j);

int main()
{	
	double** T = new double*[n];
	double** b = new double*[n];
	double** T_old = new double*[n];
	double** D = new double*[n];
	double** L = new double*[n];
	double** U = new double*[n];
	for(int i=0;i<n;i++) 
	{
		T[i] = new double[n]();
		b[i] = new double[n]();
		T_old[i] = new double[n]();
		D[i] = new double[n]();
		L[i] = new double[n]();
		U[i] = new double[n]();
	}
	double g;
	double residue=0;
	cout<<"Initialized the Variables and Allocated Memory to Matrices\n";
	//cout<<"Initialized the a-matrix\n";
	//cout<<"Enter the guess value of temperature you want to start with: ";
	//cin>>g;
	g=0;
	setup_guess(g,n,T);
	int m = 1000;
	for(int i=0;i<m;i++)
	{
		solve(n,T,T_old,b,D,U,L);
		residue = calc_res(n,T_old,T);
		cout<<"\n"<<i+1<<" iterations over...";
		cout<<"\nResidue is "<<residue;
		if(residue<0.000001)
		{
			cout<<"\nSolution converged after "<<i+1<<" iterations\n";
			ofstream file("3_2_QUICK.txt");
			for (int j=n-1;j>=0;j--)
			{  
        		for (int i=0;i<n;i++)
				{
            		file<<T[i][j]<<" "; 
        		}
        		file<<"\n";  
    		}
    		//print_Temp(n,T);
    		file.close();
			
			break;
		}
		//print_Temp(n,T);
	}
	
	for(int i=0;i<n;i++)
    {
        delete[] T[i];
		delete[] b[i];
		delete[] D[i];
		delete[] U[i];
		delete[] L[i];
    }
    delete[] T;
	delete[] b;
	delete[] D;
	delete[] U;
	delete[] L;
	return 0;
}

void init_b(int p, double** mat_T, double** mat_b, double** D, double** U, double** L)
{
	double aE=0.0,aW=0.0,aN=0.0,aS=0.0,aP=0.0,Fe=0.0,Fw=0.0,Fn=0.0,Fs=0.0,dfe=0.0,dfw=0.0,dfn=0.0,dfs=0.0;
	
	//For bottom left corner cell (i=0,j=0)
	Fe = 2.0*u(dx)*dx;
	Fw = 2.0*u(0)*dx;
	Fn = 2.0*v(dx)*dx;
	Fs = 2.0*v(0)*dx;
	aE = gamma + max(-Fe,0.0);
	aW = gamma + max(Fw,0.0);
	aN = gamma + max(-Fn,0.0);
	aS = gamma + max(Fs,0.0);
	mat_b[0][0] 	= -(aS*Ty1+aW*Tx1)-aE*mat_T[1][0]-Sc(0,0);
	aP = aE+aW+aN+aS+Fe-Fw+Fn-Fs;
	U[0][0] = aN;
	D[0][0] = -aP;
	L[0][0] = aS;
	
	//For top left corner cell (i=0,j=p-1)
	Fe = 2.0*u(dx)*dx;
	Fw = 2.0*u(0)*dx;
	Fn = 2.0*v(p*dx)*dx;
	Fs = 2.0*v((p-1)*dx)*dx;
	aE = 0.0 + max(-Fe,0.0);
	aW = 0.0 + max(Fw,0.0);
	aN = 0.0 + max(-Fn,0.0); //High Peclet number approximation along North
	aS = 0.0 + max(Fs,0.0); //
	mat_b[0][p-1] 	= -aW*Tx1-aE*mat_T[1][p-1]-Sc(0,p-1);
	aP = aE+aW+aN+aS+Fe-Fw+Fn-Fs;
	aP = aP-aN-aE;
	U[0][p-1] = aN;
	D[0][p-1] = -aP;
	L[0][p-1] = aS;
	
	//For bottom right corner cell (i=p-1,j=0)
	Fe = 2.0*u(p*dx)*dx;
	Fw = 2.0*u((p-1)*dx)*dx;
	Fn = 2.0*v(dx)*dx;
	Fs = 2.0*v(0)*dx;
	aE = 0.0 + max(-Fe,0.0); //High Peclet number approximation along East
	aW = 0.0 + max(Fw,0.0); //
	aN = 0.0 + max(-Fn,0.0);
	aS = 0.0 + max(Fs,0.0);
	mat_b[p-1][0] 	= -aS*Ty1-aW*mat_T[p-2][0]-Sc(p-1,0);
	aP = aE+aW+aN+aS+Fe-Fw+Fn-Fs;
	aP = aP-aE-aN;
	U[p-1][0] = aN;
	D[p-1][0] = -aP;
	L[p-1][0] = aS;
	
	//For top right corner cell (i=p-1,j=p-1)
	Fe = 2.0*u(p*dx)*dx;
	Fw = 2.0*u((p-1)*dx)*dx;
	Fn = 2.0*v(p*dx)*dx;
	Fs = 2.0*v((p-1)*dx)*dx;
	aE = 0.0 + max(-Fe,0.0); //High Peclet number approximation along East
	aW = 0.0 + max(Fw,0.0); //
	aN = 0.0 + max(-Fn,0.0); //High Peclet number approximation along North
	aS = 0.0 + max(Fs,0.0); //	
	mat_b[p-1][p-1] = -aW*mat_T[p-2][p-1]-Sc(p-1,p-1);
	aP = aE+aW+aN+aS+Fe-Fw+Fn-Fs;
	aP = aP-aN-aE;
	U[p-1][p-1] = aN;
	D[p-1][p-1] = -aP;
	L[p-1][p-1] = aS;	
	
	for(int j=1;j<p-1;j++)
	{
		//For left boundary
		Fe = 2.0*u(dx)*dx;
		Fw = 2.0*u(0)*dx;
		Fn = 2.0*v((j+1)*dx)*dx;
		Fs = 2.0*v(j*dx)*dx;
		aE = gamma + max(-Fe,0.0);
		aW = gamma + max(Fw,0.0);
		aN = gamma + max(-Fn,0.0);
		aS = gamma + max(Fs,0.0);
		mat_b[0][j]	= -aW*Tx1-aE*mat_T[1][j]-Sc(0,j)*dv;
		aP = aE+aW+aN+aS+Fe-Fw+Fn-Fs;
		U[0][j] = aN;
		D[0][j] = -aP;
		L[0][j] = aS;
		
		//For right boundary(upwind)
		Fe = 2.0*u(p*dx)*dx; //They are not actually zero, but max(-Fe,0) is 0 since Fe and Fn are positive
		Fw = 2.0*u((p-1)*dx)*dx;
		Fn = 2.0*v((j+1)*dx)*dx;
		Fs = 2.0*v(j*dx)*dx;
		aE = 0.0 + max(-Fe,0.0); //High Peclet number approximation along East
		aW = 0.0 + max(Fw,0.0); //
		aN = 0.0 + max(-Fn,0.0);
		aS = 0.0 + max(Fs,0.0);
		mat_b[p-1][j] = -aW*mat_T[p-2][j]-Sc(p-1,j)*dv;
		aP = aE+aW+aN+aS+Fe-Fw+Fn-Fs;
		aP = aP-aE-aN;
		U[p-1][j] = aN;
		D[p-1][j] = -aP;
		L[p-1][j] = aS;
	}
	
	for(int i=1;i<p-1;i++)
	{ 
		//For bottom boundary
		Fe = 2.0*u((i+1)*dx)*dx;
		Fw = 2.0*u(i*dx)*dx;
		Fn = 2.0*v(dx)*dx;
		Fs = 2.0*v(0)*dx;
		aE = gamma + max(-Fe,0.0);
		aW = gamma + max(Fw,0.0);
		aN = gamma + max(-Fn,0.0);
		aS = gamma + max(Fs,0.0);
		mat_b[i][0]		= -(aW*mat_T[i-1][0]+aE*mat_T[i+1][0])-aS*Ty1-Sc(i,0)*dv;
		aP = aE+aW+aN+aS+Fe-Fw+Fn-Fs;
		U[i][0] = aN;
		D[i][0] = -aP;
		L[i][0] = aS;
		
		//For top boundary(upwind)
		Fe = 2.0*u((i+1)*dx)*dx;
		Fw = 2.0*u(i*dx)*dx;
		Fn = 2.0*v(p*dx)*dx;
		Fs = 2.0*v((p-1)*dx)*dx;
		aE = 0.0 + max(-Fe,0.0);
		aW = 0.0 + max(Fw,0.0);
		aN = 0.0 + max(-Fn,0.0); //High Peclet number approximation along North
		aS = 0.0 + max(Fs,0.0); //
		mat_b[i][p-1] 	= -(aW*mat_T[i-1][p-1]+aE*mat_T[i+1][p-1])-Sc(i,p-1)*dv;
		aP = aE+aW+aN+aS+Fe-Fw+Fn-Fs;
		aP = aP-aN-aE;
		U[i][p-1] = aN;
		D[i][p-1] = -aP;
		L[i][p-1] = aS;		
	}
	
	//For all interior cells
	for(int i=1;i<p-1;i++)
	{
		for(int j=1;j<p-1;j++)
		{
			Fe = 2.0*u((i+1)*dx)*dx;
			Fw = 2.0*u(i*dx)*dx;
			Fn = 2.0*v((j+1)*dx)*dx;
			Fs = 2.0*v(j*dx)*dx;
			aE = gamma + max(-Fe,0.0);
			aW = gamma + max(Fw,0.0);
			aN = gamma + max(-Fn,0.0);
			aS = gamma + max(Fs,0.0);
			if(i>=2 && j>=2)
			{
				dfe = (mat_T[i+1][j]-mat_T[i-1][j])/4.0 + (mat_T[i+1][j]+mat_T[i-1][j]-2*mat_T[i][j])/8.0;
				dfw = (mat_T[i][j]-mat_T[i-2][j])/4.0 + (mat_T[i][j]+mat_T[i-2][j]-2*mat_T[i-1][j])/8.0;
				dfn = (mat_T[i][j+1]-mat_T[i][j-1])/4.0 + (mat_T[i][j+1]+mat_T[i][j-1]-2*mat_T[i][j])/8.0;
				dfs = (mat_T[i][j]-mat_T[i][j-2])/4.0 + (mat_T[i][j]+mat_T[i][j-2]-2*mat_T[i][j-1])/8.0;
			}
			else
			{
				dfe = (mat_T[i+1][j]-mat_T[i-1][j])/4.0 + (mat_T[i+1][j]+mat_T[i-1][j]-2*mat_T[i][j])/8.0;
				dfw = 0.0;
				dfn = (mat_T[i][j+1]-mat_T[i][j-1])/4.0 + (mat_T[i][j+1]+mat_T[i][j-1]-2*mat_T[i][j])/8.0;
				dfs = 0.0;
			}
			mat_b[i][j] = -(aW*mat_T[i-1][j]+aE*mat_T[i+1][j])-Sc(i,j)*dv+Fe*dfe-Fw*dfw+Fn*dfn-Fs*dfs;
			aP = aE+aW+aN+aS+Fe-Fw+Fn-Fs;
			U[i][j] = aN;
			D[i][j] = -aP;
			L[i][j] = aS;
		}
	}
}

void setup_guess(double guess, int p, double** mat_T)
{
	for(int i=0;i<p;i++)
	{
		for(int j=0;j<p;j++)
		{
			mat_T[i][j] = guess;
		}
	}
}

void tdma(int p, double* arr_b, double* arr_T, double* D, double* U, double* L)
{
	double* d_copy = new double[p];
	double* b_copy = new double[p];
	for(int i=0; i<p; i++) 
	{
		d_copy[i] = D[i];
		b_copy[i] = arr_b[i];
	}
	
	for(int i=1;i<p;i++)
	{
		d_copy[i]=d_copy[i]-(L[i]/d_copy[i-1])*U[i-1];
		b_copy[i]=b_copy[i]-(L[i]/d_copy[i-1])*b_copy[i-1];
	}
	
	//print_array(n,b);
	
	arr_T[p-1]=b_copy[p-1]/d_copy[p-1];
	for(int i=p-2;i>=0;i--)
	{
		arr_T[i]=(b_copy[i]-U[i]*arr_T[i+1])/d_copy[i];
	}
	delete[] d_copy;
	delete[] b_copy;
}

void solve(int p, double** mat_T, double** mat_T_old, double** mat_b, double** D, double** U, double** L)
{
	//Creating a copy of the Temperatures to calculate the residuals later
	for(int i=0;i<p;i++)
	{
		for(int j=0;j<p;j++)
		{
			mat_T_old[i][j]=mat_T[i][j];
		}
	}
	//The sweeps for each line start here
	for(int i=0;i<p;i++)
	{
		init_b(p,mat_T,mat_b,D,U,L);
		tdma(p,mat_b[i],mat_T[i],D[i],U[i],L[i]);
	}
}

double calc_res(int p, double** T_old, double** T)
{
	double res=0.0;
	for(int i=0;i<p;i++)
	{
		for(int j=0;j<p;j++)
		{
			if(res<fabs(T[i][j]-T_old[i][j]))
			{
				res = fabs(T[i][j]-T_old[i][j]);
			}
		}
	}
	return res;
}

void print_matrix(int p, int q, double** mat)
{
	cout<<"\n";
	for(int i=0;i<p;i++)
	{
		for(int j=0;j<q;j++)
		{
			cout<<mat[i][j]<<"\t";
		}
		cout<<"\n";
	}	
}

void print_Temp(int p, double** mat_T)
{
	cout<<"\n";
	for (int j=p-1;j>=0;j--) 
	{ // Start from last column
        for (int i=0;i<p;i++) 
		{  // Print top to bottom
            cout<<mat_T[i][j]<<"\t";
        }
        cout<<"\n"; // Newline after each column
    }
}

void print_array(int p, double* arr)
{
	cout<<"\n";
	for(int i=0;i<p;i++)
	{
		cout<<arr[i]<<"\t";
	}
	cout<<"\n";
}

double u(double x)
{
	return x*x+1.0;
}

double v(double y)
{
	return y*y+1.0;
}

double Sc(int i, int j)
{
	return 2*(i+0.5+j+0.5)*dx;
}

double max(double a, double b)
{
	if(a>b) return a;
	else return b;
}
