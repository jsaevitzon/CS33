#include "func.h"

void Func1(int c[][TSIZE], int a[][TSIZE], int b[][TSIZE])
{
	int i, j, k;

	int B=16;
	int i1, j1, k1;

#pragma omp parallel for private (i,j,k,i1,j1,k1) shared(c)
	for(i=0; i<TSIZE; i+=B)
	  for(j=0; j<TSIZE; j+=B)
	    for(k=0; k<TSIZE; k+=B)
	      for(i1=i; i1<i+B; i1++)
		{
		for(j1=j;j1<j+B; j1++)
		  {
		    int sum = 0;
		    for(k1=k; k1<k+B; k1++)
		      sum+=a[i1][k1]*b[k1][j1];
		    c[i1][j1]+=sum;
		  }
		}
	
}

void Func2(int d[][MSIZE], int c[][MSIZE])
{
	int i,j;
	int B=16;
	int i1, j1;

#pragma omp parallel for private(i,j,i1,j1) shared(c,d)
	for(i=0; i<MSIZE; i+=B) 
	  for(j=0; j<MSIZE; j+=B)
	    for(i1=i;i1<i+B; i1++)
	      for(j1=j; j1<j+B; j1++)
		d[i1][j1]=c[j1][i1];
}

void Func3(int z[][MSIZE], int d[][MSIZE])
{
	int y, x;
	int near = 2;  		// The weight of neighbor
	int itself = 84; 	// The weight of center


#pragma omp parallel sections
	{
#pragma omp section private(x,y) shared(d)
	  {
	    /*FOR X=0, Y=0*/
	    y= 0; x = 0;
	    z[0][0] = 
	      ( d[y][x] +
		d[y][x+1] + 
		d[y+1][x+1] +
		d[y][x] +
		d[y][x+1] +
		d[y][x] +
		d[y+1][x])<<1 +
	      itself * d[y][x];
	    z[y][x]/=100;
	  }
	  /*FOR Y=0*/
#pragma omp section private(x,y) shared(d,z)
	  {
#pragma omp parallel for private(x,y) shared(d,z)
	  for(x=1; x<MSIZE-1; x++)
	    {
	      z[y][x] = 
		(d[y][x-1] +
		 d[y][x+1] +
		 d[y+1][x-1] +
		 d[y+1][x+1] +
		 d[y][x-1] +
		 d[y][x+1] +
		 d[y][x] +
		 d[y+1][x])<<1 +
		itself * d[y][x];
	      z[y][x]/=100;
	    }
	  }
	  
	  /*For Y=0, X=MSIZE-1*/
#pragma omp section private(x,y) shared(d,z)
	  {
	    x = MSIZE-1;
	    z[y][x] = 
	      (d[y][x-1] +
	       d[y][x] +
	       d[y+1][x-1] +
	       d[y+1][x] +
	       d[y][x-1] +
	       d[y][x] +
	       d[y][x] +
	       d[y+1][x])<<1 +
	      itself * d[y][x];
	    z[y][x]/=100;
	  }

	/*FOR Y=MSIZE-1, X = 0*/
#pragma omp section private(x,y) shared(d,z)
	  {
	    x = 0; y = MSIZE-1;
	    z[y][x]= 
	      (d[y-1][x] +
	       d[y-1][x+1] +
	       d[y][x] +
	       d[y][x+1] +
	       d[y][x] +
	       d[y][x+1] +
	       d[y-1][x] +
	       d[y][x])<<1 +
	      itself * d[y][x];
	    z[y][x]/=100;
	  }
#pragma omp section private(x,y) shared(d,z)
	  {
#pragma omp parallel for private(x,y) shared (d,z)
	    for(x = 1; x<MSIZE-1; x++)
	      {
		z[y][x] =       near * d[y-1][x-1] +
		  near * d[y-1][x+1] +
		  near * d[y][x-1] +
		  near * d[y][x+1] +
		  near * d[y][x-1] +
		  near * d[y][x+1] +
		  near * d[y-1][x] +
		  near * d[y][x] +
		  itself * d[y][x];
		z[y][x]/=100;
	      }
	  }
	
	/*X=MSIZE-1*/
#pragma omp section private(x,y) shared(d,z)
	  {
	    x = MSIZE-1;
	    z[y][x] =       
	      near * d[y-1][x-1] +
	      near * d[y-1][x] +
	      near * d[y][x-1] +
	      near * d[y][x] +
	      near * d[y][x-1] +
	      near * d[y][x] +
	      near * d[y-1][x] +
	      near * d[y][x] +
	      itself * d[y][x];
	    z[y][x]/=100;
	  }
	}
	  
#pragma omp parallel for private(y,x) shared(z)
	for (y=1; y<MSIZE-1; y++) {
		for (x=0; x<MSIZE; x++) {      
		    
		       	if (x==0) {
			  z[y][x] = 	near * d[y-1][x] +
			    near * d[y-1][x+1] +
			    near * d[y+1][x] +
			    near * d[y+1][x+1] +
			    near * d[y][x] +
			    near * d[y][x+1] +
			    near * d[y-1][x] +
			    near * d[y+1][x] +
			    itself * d[y][x];
			} else if (x==MSIZE-1) {
			  z[y][x] = 	near * d[y-1][x-1] +
			    near * d[y-1][x] +
			    near * d[y+1][x-1] +
			    near * d[y+1][x] +
			    near * d[y][x-1] +
			    near * d[y][x] +
			    near * d[y-1][x] +
			    near * d[y+1][x] +
			    itself * d[y][x];
			} else {
			  z[y][x] = 	near * d[y-1][x-1] +
			    near * d[y-1][x+1] +
			    near * d[y+1][x-1] +
			    near * d[y+1][x+1] +
			    near * d[y][x-1] +
			    near * d[y][x+1] +
			    near * d[y-1][x] +
			    near * d[y+1][x] +
			    itself * d[y][x];
				}
		}
			z[y][x]/=100;
		}
	}
}
						
void Func4(int b[], int a[])
{
	int chuck_size=MSIZE;	 
	int array_size=VSIZE/chuck_size;
	int chuck[chuck_size];
    int i, j;
	
	for(j=0; j<chuck_size; j++) {
		b[j*array_size]=a[j*array_size];
		for (i=1; i<VSIZE/chuck_size; i++) {
			b[j*array_size+i]=b[j*array_size+i-1]+a[j*array_size+i];
		}
		chuck[j]=b[(j+1)*array_size-1];
	}

	//loop unroll
	
	for(j=1; j<chuck_size; j+=3) {
		chuck[j]+=chuck[j-1];
		chuck[j+1]+=chuck[j];
		chuck[j+2]+=chuck[j+1];
	}

#pragma omp parallel for private(j, i) shared(b)
	for(j=1; j<chuck_size; j++) {
	  int size = j*array_size;
	  int cuck = chuck[j-1];
	  int div = j+1;
		for (i=0; i<VSIZE/chuck_size; i++) {
		  b[size+i]+=cuck/(div);
		  //b[size+i+1]+=cuck/(div);
		}
	}
}

void Func5(int b[], int a[])
{
	int i=0, j,  stride, stride2, step;
    int temp;
	long log_N=log2(VSIZE);
	int loop_log_N=(int)log_N;
#pragma omp parallel for shared(b) private(j)
	for(j=0; j<VSIZE; j+=2) {
		b[j]=a[j];
		b[j+1] = a[j] + a[j+1];
	}

	//	#pragma omp parallel for private(i,j) shared(b)
	for(i=4; i<log_N; i=i+i) {
	  int half = i/2-1;
	  int sub = i-1;
	  for(j=0; j<log_N; j+=i) {
		  b[j+sub]+= b[j+half];
		}
	}
	
	b[VSIZE-1]=0;
	for(i=(log_N-1); i>=0; i--) {
		stride2=(2<<i)-1;
		stride=(1<<i)-1;
		step=stride2+1;
		int power = (int)pow(i,2)-1;
		int pow1  = (int)pow(i+1,2)-1;
		for(j=0; j<VSIZE; j+=pow1-1) {
                        temp=b[j+power];
			b[j+power] = b[j+pow1];
			b[j+pow1] = temp+b[j+pow1];
		}
	}
}
