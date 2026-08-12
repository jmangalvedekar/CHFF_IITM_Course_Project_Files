#include<iostream>
#include<cmath>
#include<fstream>
#include<sstream>
using namespace std;

//Temperature Boundary Conditions
const float Tx1 = 298.0; //left edge
const float Tx2 = 298.0; //right edge
const float Ty1 = 298.0; //bottom edge
const float Ty2 = 298.0; //top edge

//Defining the Domain
const int n = 20; //No. of elements/nodes
const float length = 1.0; //length of an edge of the square plate
const float gamma = 1.0; //Diffusion coefficient value(same as thermal conductivity)
const float dx = length/n; //length and breadth of an element
const float dv = dx*dx; //volume of an element
//const float S = 0; //source term
//const float Sc = S*dv/gamma; //S=Sc+Sp*Phip, Sp=0 and Sc is this
const float time = 1.0; //Time for beam movement
const float dt = 0.0001;
const float ux = 0.0; //m/s
const float P = 2.0; //Power(W)
const float A = 0.01; //Area on which beam is concentrated
const float rho = 1.0;
const float Cp = 1.0;
const float Alpha = gamma/(rho*Cp);
const float Fo = Alpha*dt/dv;

//Function Declarations
void print_array(int m, float* arr); //Prints an array
void print_matrix(int p, int q, float** mat); //Prints a matrix
void init_d(int p, float** mat); //Initializes the diagonal elements
void init_b(int p, float** mat_T, float** T_old, float** mat_b, float** q); //Initializes the b-vector
void init_q(int k, float** q); //Initialize the power matrix(moving source)
void setup_guess(float guess, int p, float** mat_T); //Sets up the guess values for temperature
void tdma(int p, float* arr_d, float* arr_b, float* arr_T); //Solves a tridiagonal matrix system
void print_Temp(int p, float** mat_T); //Prints the temperature matrix in an understandable form
void solve(int p, float** mat_T, float** mat_T_old, float** T_old_t, float** mat_b, float** d, float** q);
void solve_time(float** T, float** T_old, float** q);
float calc_res(int p, float** T_old, float** T); //Calculates the residue of T wrt T_old

int main()
{	
	//int m = int(time/dt);
	int m = 2000;
	float** T = new float*[n];
	float** b = new float*[n];
	float** T_old = new float*[n];
	float** T_old_t = new float*[n];
	float** d = new float*[2];
	float** q = new float*[n];
	for(int j=0;j<n;j++)
	{
		q[j] = new float[n]();
	}
	for(int i=0;i<n;i++) 
	{
		T[i] = new float[n]();
		b[i] = new float[n]();
		T_old[i] = new float[n]();
		T_old_t[i] = new float[n]();
		if(i<2) d[i] = new float[n]();
	}
	float g;
	cout<<"Initialized the Variables and Allocated Memory to Matrices\n";
	init_d(n,d);
	cout<<"Initialized the a-matrix\n";
	//cout<<"Enter the guess value of temperature you want to start with: ";
	//cin>>g;
	g=298.0;
	setup_guess(g,n,T);
	cout<<"Initialized the source terms\n";
	for(int k=0;k<m-1;k++)
	{
		float residue=10.0;
		int count = 0;
		init_q(k,q);
		for (int i = 0; i < n; i++) 
		{
    		for(int j = 0; j < n; j++) 
			{
        		T_old_t[i][j] = T[i][j];
    		}
		}
		while(residue>1e-6)
		{
			/*for (int i = 0; i < n; i++) 
			{
    			for (int j = 0; j < n; j++) 
				{
        			T_old[i][j] = T[i][j];
    			}
			}*/
			solve(n,T,T_old,T_old_t,b,d,q);
			residue = calc_res(n,T_old,T);
			//cout<<"\n"<<count+1<<" iterations over...";
			//cout<<"\nResidue is "<<residue;
			count+=1;
		}
		//cout<<"\nSpatial solution converged after "<<count<<" iterations in space";
		//solve_time(T,T_old,q);
		//cout<<"\n"<<k+1<<" iterations over\n";
		ostringstream filename;
		if(k%25==0)
		{
	        filename << "output_Imp_" << k << ".txt";
	        ofstream file(filename.str());
			for (int j=n-1;j>=0;j--)
			{  
		        for (int i=0;i<n;i++)
				{
		            file<<T[i][j]<<" "; 
		        }
		        file<<"\n";  
		    }
		    file.close();
		}
		residue = calc_res(n,T_old_t,T);
		if(residue<1e-4)
		{
			//cout<<"\nThe solution converged in "<<k+1<<" time steps\n";
			//break;
			//cout<<"\n%Error is "<<100.0*residue;
		}
	}
	cout<<"\nTemperature at steady state="<<T[10][10];
	
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
	mat[0][0]   = 1.0+6.0*Fo;
	mat[0][p-1] = 1.0+6.0*Fo;
	mat[1][0]   = 1.0+5.0*Fo;
	mat[1][p-1] = 1.0+5.0*Fo;
	
	for(int i=1;i<p-1;i++)
	{
		mat[0][i] = 1.0+5.0*Fo;
		mat[1][i] = 1.0+4.0*Fo;
	}
}

void init_b(int p, float** mat_T, float** T_old_t, float** mat_b, float** q)
{
	mat_b[0][0] 	= 2.0*Fo*(Ty1+Tx1)+Fo*mat_T[1][0]+q[0][0]*dt/(rho*Cp)+T_old_t[0][0];
	mat_b[0][p-1] 	= 2.0*Fo*(Ty2+Tx1)+Fo*mat_T[1][p-1]+q[0][p-1]*dt/(rho*Cp)+T_old_t[0][p-1];
	mat_b[p-1][0] 	= 2.0*Fo*(Ty1+Tx2)+Fo*mat_T[p-2][0]+q[p-1][0]*dt/(rho*Cp)+T_old_t[p-1][0];
	mat_b[p-1][p-1] = 2.0*Fo*(Ty2+Tx2)+Fo*mat_T[p-2][p-1]+q[p-1][p-1]*dt/(rho*Cp)+T_old_t[p-1][p-1];
	
	for(int j=1;j<p-1;j++)
	{
		mat_b[0][j]		= 2.0*Fo*Tx1+Fo*mat_T[1][j]+q[0][j]*dt/(rho*Cp)+T_old_t[0][j];
		mat_b[p-1][j]	= 2.0*Fo*Tx2+Fo*mat_T[p-2][j]+q[p-1][j]*dt/(rho*Cp)+T_old_t[p-1][j];
	}
	for(int i=1;i<p-1;i++)
	{
		mat_b[i][0] 	= Fo*(mat_T[i-1][0]+mat_T[i+1][0])+2.0*Fo*Ty1+q[i][0]*dt/(rho*Cp)+T_old_t[i][0];
		mat_b[i][p-1] 	= Fo*(mat_T[i-1][p-1]+mat_T[i+1][p-1])+2.0*Fo*Ty2+q[i][p-1]*dt/(rho*Cp)+T_old_t[i][p-1];
	}
	
	for(int i=1;i<p-1;i++)
	{
		for(int j=1; j<p-1;j++)
		{
			mat_b[i][j] = Fo*(mat_T[i-1][j]+mat_T[i+1][j])+q[i][j]*dt/(rho*Cp)+T_old_t[i][j];
		}
	}
}

void init_q(int k, float** q)
{
	//for (int i = 0; i < n; i++)
        //for (int j = 0; j < n; j++)
            //q[i][j] = 0.0;
	//int i = 3+int((ux*k*dt)/dx);
	int i = int(0.5/dx);
	if(i<n-3)
	{
		for(int ix=-2;ix<3;ix++)
		{
			for(int iy=-2;iy<3;iy++)
			{
				q[i+ix][i+iy] = P/A;
			}
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
		l[i+1]=-Fo;
		u[i]=-Fo;
	}
	u[0]=-Fo;
	
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

void solve(int p, float** T, float** T_old, float** T_old_t, float** b, float** d, float** q)
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
	init_b(p,T,T_old_t,b,q);
	tdma(p,d[0],b[0],T[0]);
	for(int i=1;i<p-1;i++)
	{
		init_b(p,T,T_old_t,b,q);
		tdma(p,d[1],b[i],T[i]);
	}
	init_b(p,T,T_old_t,b,q);
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
