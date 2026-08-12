#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

//Flux Boundary Conditions
const float q_ = -100.0f; //All edges; its negative because it is exactly opposite to the area

//Defining the Domain
const int n = 10; //No. of elements/nodes
const float length = 1.0f; //length of an edge of the square plate
const float gamma_ = 10.0f; //Diffusion coefficient value
const float dx = length/n; //length and breadth of an element
const float dv = dx*dx; //volume of an element
const float S = 4*q_/length; //source term
const float Sc = S*dv/gamma_; //S=Sc+Sp*Phip, Sp=0 and Sc is this
const float q = q_*dx/gamma_; //modified heat flux
const float alpha = 1; //relaxation parameter
const float beta_ = (1-alpha); //another relaxation parameter

//Function Declarations
void print_array(int m, float* arr); //Prints an array
void print_matrix(int p, int q, float** mat); //Prints a matrix
void init_d(int p, float** mat); //Initializes the diagonal elements
void init_b(int p, float** mat_T, float** mat_T_old, float** mat_b); //Initializes the b-vector
void setup_guess(float guess, int p, float** mat_T); //Sets up the guess values for temperature
void tdma(int p, float* arr_d, float* arr_b, float* arr_T); //Solves a tridiagonal matrix system
void print_Temp(int p, float** mat_T); //Prints the temperature matrix in an understandable form
void solve(int p, float** mat_T, float** mat_T_old, float** mat_b, float** d); //Solves for one iteration
float calc_res(int p, float** T_old, float** T); //Calculates the residue of T wrt T_old
float avg_temp(float** temp); //Calculates and returns the bulk temperature 
//float integrateRow(int x, float** temp); //Part of calculations of Simpson's rule

int main()
{	
	float** T = new float*[n]; //Definition of the temperature matrix
	float** b = new float*[n]; //Definition of the b-matrix
	float** T_old = new float*[n]; //Definition of previous iteration temperature matrix
	float** d = new float*[2]; //Definition of the diagonal values of the the a-matrix
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
	init_d(n,d); //Initialize the a-matrix
	cout<<"Initialized the a-matrix\n";
	g=100; //Enter the guess value of temperature you want to start with
	setup_guess(g,n,T); //Sets up a guess value to start with
	int m = 10000; //Number of iterations
	for(int i=0;i<m;i++)
	{
		solve(n,T,T_old,b,d);
		residue = calc_res(n,T_old,T);
		cout<<"\n"<<i+1<<" iterations over...";
		cout<<"\nResidue is "<<residue;
		if(residue<0.001)
		{
			cout<<"\nSolution converged after "<<i+1<<" iterations\n";
			ofstream file("2_2.txt");
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

//Function Definitions
void init_d(int p, float** mat)
{
	//d[0] contains the matrix diagonal for the boundary indices
	//d[1] contains the matrix diagonal for the remaining indices
	mat[0][0]   = -2.0/alpha;
	mat[0][p-1] = -2.0/alpha;
	mat[1][0]   = -3.0/alpha;
	mat[1][p-1] = -3.0/alpha;
	
	for(int i=1;i<p-1;i++)
	{
		mat[0][i] = -3.0/alpha;
		mat[1][i] = -4.0/alpha;
	}
}

void init_b(int p, float** mat_T, float** T_old, float** mat_b)
{
	mat_b[0][0] 	= 2*q-mat_T[1][0]-Sc-2*T_old[0][0]*beta_/alpha;
	mat_b[0][p-1] 	= 2*q-mat_T[1][p-1]-Sc-2*T_old[0][p-1]*beta_/alpha;
	mat_b[p-1][0] 	= 2*q-mat_T[p-2][0]-Sc-2*T_old[p-1][0]*beta_/alpha;
	mat_b[p-1][p-1] = 2*q-mat_T[p-2][p-1]-Sc-2*T_old[p-1][p-1]*beta_/alpha;
	
	for(int j=1;j<p-1;j++)
	{
		mat_b[0][j]		= q-mat_T[1][j]-Sc-3*T_old[0][j]*beta_/alpha;
		mat_b[p-1][j]	= q-mat_T[p-2][j]-Sc-3*T_old[p-1][j]*beta_/alpha;
	}
	for(int i=1;i<p-1;i++)
	{
		mat_b[i][0] 	= q-(mat_T[i-1][0]+mat_T[i+1][0])-Sc-3*T_old[i][0]*beta_/alpha;
		mat_b[i][p-1] 	= q-(mat_T[i-1][p-1]+mat_T[i+1][p-1])-Sc-3*T_old[i][p-1]*beta_/alpha;
	}
	
	for(int i=1;i<p-1;i++)
	{
		for(int j=1; j<p-1;j++)
		{
			mat_b[i][j] = -(mat_T[i-1][j]+mat_T[i+1][j])-Sc-4*T_old[i][j]*beta_/alpha;
		}
	}
	//Creating a copy of the Temperatures to calculate the residuals later
	for(int i=0;i<p;i++)
	{
		for(int j=0;j<p;j++)
		{
			T_old[i][j]=mat_T[i][j];
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
	//The sweeps for each line start here
	init_b(p,T,T_old,b);
	tdma(p,d[0],b[0],T[0]);
	for(int i=1;i<p-1;i++)
	{
		init_b(p,T,T_old,b);
		tdma(p,d[1],b[i],T[i]);
	}
	init_b(p,T,T_old,b);
	tdma(p,d[0],b[p-1],T[p-1]);
}

float calc_res(int p, float** T_old, float** T)
{
	float res=0.0;
	float Tb = avg_temp(T);
	float Tb_old = avg_temp(T_old);
	for(int i=0;i<p;i++)
	{
		for(int j=0;j<p;j++)
		{
			if(res<fabs(T[i][j]-Tb-T_old[i][j]+Tb_old))
			{
				res = fabs(T[i][j]-Tb-T_old[i][j]+Tb_old);
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
        cout<<"\n";
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

/*float integrateRow(int x, float** temp) 
{
    float sum = 0.0f;
    
    // Apply Simpson's rule
    sum += temp[x][0] + temp[x][n-1];  // First and last term

    for (int j=1;j<n-1;j++) 
	{
        if (j%2== 0) 
		{
            sum+=2*temp[x][j];  // Even index terms
        } 
		else 
		{
            sum+=4*temp[x][j];  // Odd index terms
        }
    }

    return sum*dx/3.0f;
}

// Function to compute Simpson's rule for the entire 2D region (double integration)
float avg_temp(float** temp) 
{
    float totalSum=0.0f;
    
    // Apply Simpson's rule over all rows
    totalSum+=integrateRow(0,temp)+integrateRow(n-1,temp);  // First and last row

    for (int i=1;i<n-1;i++) 
	{
        if (i%2==0) 
		{
            totalSum+=2*integrateRow(i,temp);  // Even index rows
        } 
		else 
		{
            totalSum+=4*integrateRow(i,temp);  // Odd index rows
        }
    }

    return (totalSum*dx/3.0f)/(n*n);  // Normalize to get the average
}*/
float avg_temp(float** temp)
{
	float sum=0.0;
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			sum+=temp[i][j];
		}
	}
	return sum/n*n;
}
