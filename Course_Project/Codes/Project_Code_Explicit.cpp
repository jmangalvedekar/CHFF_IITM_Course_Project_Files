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

//Function Declarations
void print_array(int m, float* arr); //Prints an array
void print_matrix(int p, int q, float** mat); //Prints a matrix
void init_q(int k, float** q); //Initialize the power matrix(moving source)
void setup_guess(float guess, int p, float** mat_T); //Sets up the guess values for temperature
void print_Temp(int p, float** mat_T); //Prints the temperature matrix in an understandable form
void solve_time(float** T, float** T_old, float** q);
float calc_res(int p, float** T_old, float** T); //Calculates the residue of T wrt T_old

int main()
{	
	//int m = int(time/dt);
	int m = 2000;
	float** T = new float*[n];
	float** b = new float*[n];
	float** T_old = new float*[n];
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
		if(i<2) d[i] = new float[n]();
	}
	float g;
	cout<<"Initialized the Variables and Allocated Memory to Matrices\n";
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
		solve_time(T,T_old,q);
		ostringstream filename;
		if(k%25==0)
		{
	        filename << "output_Exp_" << k << ".txt";
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
		residue = calc_res(n,T_old,T);
		if(residue<1e-4)
		{
			//cout<<"\nThe solution converged in "<<k+1<<" time steps\n";
			//cout<<"\n%Error is "<<100.0*residue;
			//break;
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
	/*for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			q[i][j] = P/A;
		}
	}*/
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

void solve_time(float** T, float** T_old, float** q)
{
	for (int i = 0; i < n; i++) 
	{
    	for (int j = 0; j < n; j++) 
		{
        	T_old[i][j] = T[i][j];
    	}
	}
	float Fo = Alpha*dt/dv;
	for(int i=1;i<n-1;i++)
	{
		for(int j=1;j<n-1;j++)
		{
			T[i][j] = Fo*(T_old[i+1][j]+T_old[i-1][j]+T_old[i][j+1]+T_old[i][j-1])+(1.0-4.0*Fo)*T_old[i][j]+q[i][j]*dt/(rho*Cp);
		}
	}
	
	T[0][0] = Fo*(2.0*Tx1+2.0*Ty1+T_old[1][0]+T_old[0][1])+(1.0-6.0*Fo)*T_old[0][0]+q[0][0]*dt/(rho*Cp); //Bottom left corner
	T[0][n-1] = Fo*(2.0*Tx1+2.0*Ty2+T_old[1][n-1]+T_old[0][n-2])+(1.0-6.0*Fo)*T_old[0][n-1]+q[0][n-1]*dt/(rho*Cp); //Top left corner
	T[n-1][0] = Fo*(2.0*Tx2+2.0*Ty1+T_old[n-2][0]+T_old[n-1][1])+(1.0-6.0*Fo)*T_old[n-1][0]+q[n-1][0]*dt/(rho*Cp); //Bottom right corner
	T[n-1][n-1] = Fo*(2.0*Tx2+2.0*Ty2+T_old[n-2][n-1]+T_old[n-1][n-2])+(1.0-6.0*Fo)*T_old[n-1][n-1]+q[n-1][n-1]*dt/(rho*Cp); //Top right corner
	
	//For left and right boundaries
	for(int j=1;j<n-1;j++)
	{
		T[0][j] = Fo*(2.0*Tx1+T_old[1][j]+T_old[0][j+1]+T_old[0][j-1])+(1.0-5.0*Fo)*T_old[0][j]+q[0][j]*dt/(rho*Cp);
		T[n-1][j] = Fo*(2.0*Tx2+T_old[n-2][j]+T_old[n-1][j+1]+T_old[n-1][j-1])+(1.0-5.0*Fo)*T_old[n-1][j]+q[n-1][j]*dt/(rho*Cp);
	}
	
	//For top and bottom boundaries
	for(int i=1;i<n-1;i++)
	{
		T[i][0] = Fo*(2.0*Ty1+T_old[i+1][0]+T_old[i-1][0]+T_old[i][1])+(1.0-5.0*Fo)*T_old[i][0]+q[i][0]*dt/(rho*Cp);
		T[i][n-1] = Fo*(2.0*Ty2+T_old[i+1][0]+T_old[i-1][0]+T_old[i][n-2])+(1.0-5.0*Fo)*T_old[i][n-1]+q[i][n-1]*dt/(rho*Cp);
	}
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
