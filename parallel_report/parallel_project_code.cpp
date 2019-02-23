#include <bits/stdc++.h>
#include <omp.h>
#include <time.h>
using namespace std;

struct node
{
	long long  data;
	struct node* left;
	struct node* right;
}typedef node;

node* root;

void dfs(vector<vector<long long > > adj,vector<long long > &vis,node* par,long long  t,vector<node*> &ptr)
{
	vis[t]=1;
	node* temp=new node;
	temp->data=t;
	temp->left=NULL;
	temp->right=NULL;
	
	ptr[t]=temp;
	if(par->left==NULL)
		par->left=temp;
	else if(par->right==NULL)
		par->right=temp;
	
	for(long long i=0;i<adj[t].size();i++)
	{
		if(vis[adj[t][i]]==0)
			dfs(adj,vis,temp,adj[t][i],ptr);
	}
}


void elimination(vector<vector<long long > > adj,vector<node*> &ptr)
{
	long long  i,j,n=adj.size();
	i=n-1;
	while(adj[i].size()==0){
		i--;
	}
		
		
	root=new node;
	root->data=i;
	root->left=NULL;
	root->right=NULL;
	
	ptr[i]=root;
	
	vector<long long > vis(n,0);
	vis[i]=1;
	for(long long  j=0;j<adj[i].size();j++)
	{
		if(vis[adj[i][j]]==0)
			dfs(adj,vis,root,adj[i][j],ptr);
	}
}

void get_parent(vector<vector<double> > m ,vector<long long > &parent)
{

	#pragma omp parallel for
	for(long long i=0;i<parent.size();i++)
	{
		for(long long j=i+1;j<parent.size();j++)
		{
			if(m[j][i]!=0){
				parent[i]=j;
				break;
			}
		}
	}
}

void adjacent(vector<vector<long long > > &adj,vector<long long > parent)
{
	
	for(long long i=0;i<parent.size();i++)
	{
		if(parent[i]!=-1)
		{
			adj[i].push_back(parent[i]);
			adj[parent[i]].push_back(i);
		}
	}
}


void divide_seq(long long  i,vector<vector<double> > &a)
{
		
	a[i][i] = sqrt(a[i][i]);
	
	for(long long  k=i+1;k<a.size();k++)
	{
		a[k][i] = a[k][i]/a[i][i];
	}
}

void modify_seq(long long  i,long long  j,vector<vector<double> > &a)
{
	long long  n=a.size();
	
	for(long long  k=i;k<a.size();k++)
	{
		a[k][i] = a[k][i] -a[k][j]*a[i][j];
	}
}





void sequential_cholesky(vector<vector<double> > &a)
{
	
	long long  n=a.size();
	
	for(long long  i=0;i<n;i++)
	{
		for(long long  j=0;j<i;j++)
			modify_seq(i,j,a);

		divide_seq(i,a);
	 }	
	 
	for(long long  i=0;i<n;i++)
	{
		for(long long  j=i+1;j<n;j++)
		{
			a[i][j] = 0;
		}
	}
}


void divide_para(long long  i,vector<vector<double> > &a)
{

	a[i][i] = sqrt(a[i][i]);
	
	#pragma omp parallel for
	for(long long  k=i+1;k<a.size();k++)
	{
		a[k][i] = a[k][i]/a[i][i];
	}
}


void modify_para(long long  i,long long  j,vector<vector<double> > &a)
{
	long long  n=a.size();
	
	
	#pragma omp parallel for
	for(long long  k=i;k<a.size();k++)
	{
		a[k][i] = a[k][i] -a[k][j]*a[i][j];
	}
}






void parallel_cholesky(vector<node*> &ptr,vector<vector<double> > &a,vector<long long > parent)
{
	long long  n,sum=0;
	n=ptr.size();
	#pragma omp parallel for reduction(+ : sum)
	for(long long  i=0;i<ptr.size();i++)
	{
		if(ptr[i]!=NULL)
			sum++;
	}
	while(sum>0)
	{
		#pragma omp parallel for
		for(long long  i=0;i<ptr.size();i++)
		{
			if(ptr[i]!=NULL && ptr[i]->left==NULL && ptr[i]->right==NULL)
			{
				divide_para(i,a);
				node* temp2=ptr[i];
				ptr[i]=NULL;
				if(parent[i]!=-1){
					node* temp1=ptr[parent[i]];
				if(temp1->left==temp2)
					temp1->left=NULL;
				else if(temp1->right==temp2)
					temp1->right=NULL;
				}
			}
		}
		
		#pragma omp parallel for
		for(long long  i=0;i<ptr.size();i++)
		{
			if(ptr[i]!=NULL && ptr[i]->left==NULL && ptr[i]->right==NULL)
			{
				long long  j;
				#pragma omp parallel for
				for(long long  j=0;j<i;j++)
				{
					if(a[i][j]!=0)
						modify_para(i,j,a);
				} 
			}
		}
		sum=0;
		#pragma omp parallel for reduction(+ : sum)
		for(long long  i=0;i<ptr.size();i++)
		{
			if(ptr[i]!=NULL)
				sum++;
		}
	}
}

void print_matrix(vector<vector<double> > v)
{

	cout<<setprecision(8)<<endl;
	for(long long  i=0;i<v.size();i++)
	{
		for(long long  j=0;j<v[i].size();j++)
		{
			if(i<j)
				cout<<"0.000000"<<" ";
				else
			printf("%f ",v[i][j]);
		}
		cout<<endl;
	}
	cout<<endl;
}


int main()
{
    clock_t t0,t1;
	long long  n;
	vector< vector<double> > m1{ {32,19,0,0,0},
                   	             {19,44,25,0,0},
                   	             {0,25,74,34,0},
                   	             {0, 0, 34,143,47},
									{0, 0, 0,47,247}};
                   	
	vector<vector<double> > m2{{32,19,0,0,0},
                   	             {19,44,25,0,0},
                   	             {0,25,74,34,0},
                   	             {0, 0, 34,143,47},
									{0, 0, 0,47,247}};            	
	n=m1.size();
	
	vector<long long > parent(n,-1);
   	
   	get_parent(m1,parent);
   	
	vector<node*> ptr(n,NULL);
   	
	vector<vector<long long > > adj(n,vector<long long >());
   	
	adjacent(adj,parent);
   	
	elimination(adj,ptr);


    t0 = clock();
	sequential_cholesky(m2);
	t1 = clock();

   	printf("Lower triangular matrix for sequential time : %f seconds\n",(double)(t1-t0)/CLOCKS_PER_SEC);
	cout<<setprecision(8);
	print_matrix(m2);
	
	cout<<endl<<endl;   	
   	
	t0 = clock();  	
	parallel_cholesky(ptr,m1,parent);
    t1 = clock();
	
	printf("Lower triangular matrix for parallel method time: %f seconds\n",(double)(t1-t0)/CLOCKS_PER_SEC);
	
	cout<<setprecision(8);
   	print_matrix(m1);
   	
    return 0;
}

