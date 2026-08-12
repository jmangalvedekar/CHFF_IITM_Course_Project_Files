#include<iostream>
#include<cmath>
using namespace std;

void print_matrix(int p, int** mat);
void print_array(int p, float* arr);
int main()
{
	int n=5,m=100,converged_index=0;
	float b[n]={-204,-12,-20,-28,-1236},c[n]={0};
	int** a = new int*[n];
	float** phi = new float*[n];
	for(int i=0;i<n;i++) 
	{
	a[i] = new int[n]();
	phi[i] = new float[m]();
	}
	for(int i=0;i<n;i++)
	{
		/*for(int j=0;j<n;j++)
		{
			a[i][j]=0;
		}	
		for(int j=0;j<m;j++)
		{
			phi[i][j]=0;
		}*/
		phi[i][0]=200;
		cout<<phi[i][0]<<"\t";
	}
	for(int i=0;i<n-1;i++)
	{
		a[i][i]=-2;
		a[i][i+1]=1;
		a[i+1][i]=1;
	}
	a[0][0]=-3;
	a[n-1][n-1]=-3;
	
	print_array(n,b);
	print_matrix(n,a);
	cout<<"\n";
	
	for(int k=0;k<m-1;k++)
	{
		float max_diff = 0.0;
		for(int i=0;i<n;i++)
		{
			c[i]=0;
			for(int j=0;j<n;j++)
			{
				//if(j<i) c[i]+=a[i][j]*phi[j][k+1];
				//if(j>i) c[i]+=a[i][j]*phi[j][k+1];
				if (j != i) c[i] += a[i][j] * phi[j][k];
			}
			phi[i][k+1] = (b[i]-c[i])/a[i][i];
			float diff=fabs(phi[i][k+1]-phi[i][k]);
			if(diff>max_diff) max_diff=diff;
			cout<<phi[i][k+1]<<"\t";
		}
		cout << "\nMax change: " << max_diff << "\n";
		if(max_diff<0.0001)
		{
			cout<<"\nSolution converged at iteration "<<k+1;
			converged_index=k+1;
			break;
		}
		cout<<"\n";
		
	}
	
	float phi_exact[5]={166.5,295.5,412.5,509.5,578.5},phi_diff[n];
	cout<<"\n\nThe exact solution is ";
	print_array(5,phi_exact);
	cout<<"\nThe calculated solution is \n";
	for(int i=0;i<n;i++)
	{
		cout<<phi[i][converged_index]<<"\t";
		phi_diff[i]=100*(fabs(phi_exact[i]-phi[i][converged_index]))/phi_exact[i];
	}
	cout<<"\n\nThe percentage Error in the Computed Solution at Each Cell Centroid is \n";
	print_array(n,phi_diff);
	
	for(int i = 0; i < n; i++)
    {
        delete[] a[i];
		delete[] phi[i];
    }
    delete[] a;
	delete[] phi;
	return 0;
}

void print_matrix(int p, int** mat)
{
	cout<<"\n";
	for(int i=0;i<p;i++)
	{
		for(int j=0;j<p;j++)
		{
			cout<<mat[i][j]<<"\t";
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
