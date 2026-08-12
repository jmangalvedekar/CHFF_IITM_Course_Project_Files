#include<iostream>
using namespace std;

void print_matrix(int m, int** mat);
void print_array(int m, float* arr);
int main()
{
	int n=5;
	float b[n]={-204,-12,-20,-28,-1236},x[n],c[n],d[n],e[n];
	int** a = new int*[n];
	int** l = new int*[n];
	int** u = new int*[n];
	for(int i=0;i<n;i++) 
	{
	a[i] = new int[n];
	l[i] = new int[n];
	u[i] = new int[n];
	}
	//float** c = a;
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			a[i][j]=0;
		}	
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
	
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			if(j<=i)
			{l[i][j]=a[i][j];}
			else{l[i][j]=0;}
			u[i][j]=a[i][j]-l[i][j];
		}	
	}
	
	print_matrix(n,l);
	print_matrix(n,u);
	
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			c[i] += u[i][j]*x[j];
		}
		d[i]=b[i]-c[i];
	}
	
	x[0]=d[0]/l[0][0];
	for(int i=1;i<n;i++)
	{
		for(int j=0;j<i;j++)
		{
			e[i]+=l[i][j]*x[j];
		}
		x[i]=(d[i]-e[i])/l[i][i];
		if(x[i]-x[i-1]<0.001)
		{
			break;
		}
	}
	
	for(int i = 0; i < n; i++)
    {
        delete[] a[i];
		//delete[] c[i];// Deallocate each row
    }
    delete[] a;
	//delete[] c;// Deallocate the array of pointers
	return 0;
}

void print_matrix(int m, int** mat)
{
	cout<<"\n";
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<m;j++)
		{
			cout<<mat[i][j]<<"\t";
		}
		cout<<"\n";
	}
}

void print_array(int m, float* arr)
{
	cout<<"\n";
	for(int i=0;i<m;i++)
	{
		cout<<arr[i]<<"\t";
	}
	cout<<"\n";
}
