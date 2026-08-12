#include<iostream>
using namespace std;

void print_array(int m, float* arr);
void print_matrix(int p, float** mat);

int main()
{	
	int n=5;
	float b[n]={-204,-12,-20,-28,-1236};
	float x[n]={0},y=0;
	float** a = new float*[n];
	for(int i=0;i<n;i++) 
	{a[i] = new float[n]();}
	for(int i=0;i<n-1;i++)
	{
		a[i][i]=-2;
		a[i][i+1]=1;
		a[i+1][i]=1;
	}
	a[0][0]=-3;
	a[n-1][n-1]=-3;
	
	float l[n]={0},u[n]={0},d[n]={0};//,d_[n]={0},b_[n]={0};
	
	for(int i=0;i<n;i++)
	{
		d[i]=a[i][i];
	}
	for(int i=0;i<n-1;i++)
	{
		u[i]=a[i][i+1];
	}
	u[0]=1;
	
	for(int i=0;i<n-1;i++)
	{
		l[i+1]=a[i+1][i];
	}
	
	//print_matrix(n,a);
	//print_array(n,l);
	//print_array(n,u);
	//print_array(n,d);
	//print_array(n,x);
	
	for(int i=1;i<n;i++)
	{
		d[i]=d[i]-(l[i]/d[i-1])*u[i-1];
		b[i]=b[i]-(l[i]/d[i-1])*b[i-1];
	}
	
	//print_array(n,b);
	
	x[n-1]=b[n-1]/d[n-1];
	for(int i=n-2;i>=0;i--)
	{
		x[i]=(b[i]-u[i]*x[i+1])/d[i];
	}

	print_array(n,x);
	return 0;
}

void print_matrix(int p, float** mat)
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
