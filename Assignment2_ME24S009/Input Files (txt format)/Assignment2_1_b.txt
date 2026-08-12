#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

//Temperature Boundary Conditions
const float Tx1 = 50; //left edge
const float Tx2 = 150; //right edge
const float Ty1 = 200; //bottom edge
const float Ty2 = 100; //top edge

//Defining the Domain
const int n = 20; //No. of elements/nodes
const float length = 1.0; //length of an edge of the square plate
const float gamma = 0.01; //Diffusion coefficient value
const float dx = length/n; //length and breadth of an element
const float dv = dx*dx; //volume of an element
const float S = 0; //source term
const float Sc = S*dv/gamma; //S=Sc+Sp*Phip, Sp=0 and Sc is this

//Function Declarations
void print_array(int m, float* arr); //Prints an array
void print_matrix(int p, int q, float** mat); //Prints a matrix
void init_d(int p, float** mat); //Initializes the diagonal elements
void init_b(int p, float** mat_T, float** mat_b); //Initializes the b-vector
void setup_guess(float guess, int p, float** mat_T); //Sets up the guess values for temperature
void tdma(int p, float* arr_d, float* arr_b, float* arr_T); //Solves a tridiagonal matrix system
void print_Temp(int p, float** mat_T); //Prints the temperature matrix in an understandable form
void solve(int p, float** mat_T, float** mat_T_old, float** mat_b, float** d);
float calc_res(int p, float** T_old, float** T); //Calculates the residue of T wrt T_old

int main()
{	
	float** T = new float*[n];
	float** b = new float*[n];
	float** T_old = new float*[n];
	float** d = new float*[2];
	for(int i=0;i<n;i++) 
	{
		T[i] = new float[n]();
		b[i] = new float[n]();
		T_old[i] = new float[n]();
		if(i<2) d[i] = new float[n]();
	}
	float g;
	float residue=0;
	cout<<"Initialized the Variables and Allocated Memory to Matrices\n";
	init_d(n,d);
	cout<<"Initialized the a-matrix\n";
	//cout<<"Enter the guess value of temperature you want to start with: ";
	//cin>>g;
	g=190;
	setup_guess(g,n,T);
	int m = 1000;
	for(int i=0;i<m;i++)
	{
		solve(n,T,T_old,b,d);
		residue = calc_res(n,T_old,T);
		cout<<"\n"<<i+1<<" iterations over...";
		cout<<"\nResidue is "<<residue;
		if(residue<0.000001)
		{
			cout<<"\nSolution converged after "<<i+1<<" iterations\n";
			ofstream file("2_1_b.txt");
			for (int j=n-1;j>=0;j--)
			{  
        		for (int i=0;i<n;i++)
				{
            		file<<T[i][j]<<" "; 
        		}
        		file<<"\n";  
    		}
    		file.close();
			
			break;
		}
	}
	
	for(int i=0;i<n;i++)
    {
        delete[] T[i];
		delete[] b[i];
		if(i<2) delete[] d[i];
    }
    delete[] T;
	delete[] b;
	delete[] d;
	return 0;
}

void init_d(int p, float** mat)
{
	//d[0] contains the matrix diagonal for the boundary indices
	//d[1] contains the matrix diagonal for the remaining indices
	mat[0][0]   = -6;
	mat[0][p-1] = -6;
	mat[1][0]   = -5;
	mat[1][p-1] = -5;
	
	for(int i=1;i<p-1;i++)
	{
		mat[0][i] = -5;
		mat[1][i] = -4;
	}
}

void init_b(int p, float** mat_T, float** mat_b)
{
	mat_b[0][0] 	= -2*(Ty1+Tx1)-mat_T[1][0]-Sc;
	mat_b[0][p-1] 	= -2*(Ty2+Tx1)-mat_T[1][p-1]-Sc;
	mat_b[p-1][0] 	= -2*(Ty1+Tx2)-mat_T[p-2][0]-Sc;
	mat_b[p-1][p-1] = -2*(Ty2+Tx2)-mat_T[p-2][p-1]-Sc;
	
	for(int j=1;j<p-1;j++)
	{
		mat_b[0][j]		= -2*Tx1-mat_T[1][j];
		mat_b[p-1][j]	= -2*Tx2-mat_T[p-2][j];
	}
	for(int i=1;i<p-1;i++)
	{
		mat_b[i][0] 	= -(mat_T[i-1][0]+mat_T[i+1][0])-2*Ty1-Sc;
		mat_b[i][p-1] 	= -(mat_T[i-1][p-1]+mat_T[i+1][p-1])-2*Ty2-Sc;
	}
	
	for(int i=1;i<p-1;i++)
	{
		for(int j=1; j<p-1;j++)
		{
			mat_b[i][j] = -(mat_T[i-1][j]+mat_T[i+1][j])-Sc;
		}
	}
}

void setup_guess(float guess, int p, float** mat_T)
{
	for(int i=0;i<p;i++)
	{
		for(int j=0;j<p;j++)
		{
			mat_T[i][j] = guess;
		}
	}
}

void tdma(int p, float* arr_d, float* arr_b, float* arr_T)
{
	float l[p]={0},u[p]={0};
	float* d_copy = new float[p];
	float* b_copy = new float[p];
	for(int i=0; i<p; i++) 
	{
		d_copy[i] = arr_d[i];
		b_copy[i] = arr_b[i];
	}
	for(int i=0;i<p-1;i++)
	{
		l[i+1]=1;
		u[i]=1;
	}
	u[0]=1;
	
	for(int i=1;i<p;i++)
	{
		d_copy[i]=d_copy[i]-(l[i]/d_copy[i-1])*u[i-1];
		b_copy[i]=b_copy[i]-(l[i]/d_copy[i-1])*b_copy[i-1];
	}
	
	//print_array(n,b);
	
	arr_T[p-1]=b_copy[p-1]/d_copy[p-1];
	for(int i=p-2;i>=0;i--)
	{
		arr_T[i]=(b_copy[i]-u[i]*arr_T[i+1])/d_copy[i];
	}
	delete[] d_copy;
	delete[] b_copy;
}

void solve(int p, float** T, float** T_old, float** b, float** d)
{
	//Creating a copy of the Temperatures to calculate the residuals later
	for(int i=0;i<p;i++)
	{
		for(int j=0;j<p;j++)
		{
			T_old[i][j]=T[i][j];
		}
	}
	//The sweeps for each line start here
	init_b(p,T,b);
	tdma(p,d[0],b[0],T[0]);
	for(int i=1;i<p-1;i++)
	{
		init_b(p,T,b);
		tdma(p,d[1],b[i],T[i]);
	}
	init_b(p,T,b);
	tdma(p,d[0],b[p-1],T[p-1]);
}

float calc_res(int p, float** T_old, float** T)
{
	float res=0.0;
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

void print_matrix(int p, int q, float** mat)
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

void print_Temp(int p, float** mat_T)
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

void print_array(int p, float* arr)
{
	cout<<"\n";
	for(int i=0;i<p;i++)
	{
		cout<<arr[i]<<"\t";
	}
	cout<<"\n";
}
