#include <stdio.h> 
#include <string.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include <grace_np.h>
#include <unistd.h>
#include <float.h>
#include <limits.h>
#include <signal.h>
#include <cash2003.h>
#include <cash2.h>
#include <mersenne.h>
#include <cash2-s.h>


static TYPE2** Life;

void Initial(void)
{
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 200; /* # of row (default=100)*/
  ncol = 200; /* # of column (default=100)*/
  nplane = 1; /* # of planes (default=0)*/
  scale = 2; /* size of the window (default=2)*/
  boundary = WRAP; /* the type of boundary: FIXED, WRAP, ECHO (default=WRAP). Note that
		      Margolus diffusion is not supported for ECHO. */
  ulseedG = 56; /* random seed (default=56)*/

  /* useally, one does not have to change the followings */
  /* the value of boundary (default=(TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.})*/
  boundaryvalue2 = (TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.}; 
}

void InitialPlane(void)
{
  MakePlane(&Life);

  /* InitialSet(1,2,3,4,5)
    1: name of plane
    2: number of the state other than the background state
    3: background state or empty state
    4: state I want to put
    5: fraction of cells that get S1 state (0 to 1)
  */
  InitialSet(Life,3,0,1,0.01,2,0.01,3,0.01);
  //ReadSavedData("glidergun.sav",1,Life);
  //printf("\n\nLife is ready. Click or scroll your mouse to update the CA.\n\n");
  Boundaries2(Life);
}

double q=0.5;
double p1=0.1;
double p2=100.;
double p3=1.;

double a1=1.;
double a2=0.1;
double a3=1.;
double death=0.1;
double NONE=5;

void NextState(int row,int col)
{
  int sum;
  
  if(Life[row][col].val==0){
    //replication
    double fsum1 = a1*CountMoore8(Life,1,row,col);
    double fsum2 = a2*CountMoore8(Life,2,row,col);
    double fsum3 = a3*CountMoore8(Life,3,row,col);
    if(fsum1+fsum2+fsum3 ==0 ) return;
    
    double rn = (genrand_real1() * (fsum1+fsum2+fsum3+NONE));
    
    if(rn<fsum1){
      //val=1 replicates, possibly mutates into 2
      if(genrand_real2() < q)
        Life[row][col].val=1;
      else Life[row][col].val=2;
      // there should be a possibility of making garbage, or dead stuff, or "almost" dead stuf
    }else if(rn<fsum1+fsum2){
      Life[row][col].val=2;
    }else if(rn<fsum1+fsum2+fsum3){
      Life[row][col].val=3;
    }
  }else{
    //competition - you are killed with a certain prob
    // proportional to how many of the other type are there
    
    if(Life[row][col].val==3){
      double fsum1 = p1*CountMoore8(Life,1,row,col);
      double fsum2 = p2*CountMoore8(Life,2,row,col);
      if(fsum1+fsum2 ==0 ) return;
      double rn = (genrand_real1() * (fsum1+fsum2+NONE));
      if(rn<(fsum1+fsum2)) Life[row][col].val=0;
    }else if(Life[row][col].val==1 || Life[row][col].val==2){
      double fsum3 = p3*CountMoore8(Life,3,row,col);
      if(fsum3 ==0 ) return;
      double rn = (genrand_real1() * (fsum3+NONE));
      if(rn<(fsum3)) Life[row][col].val=0;
    }
  }
  
}

void Death(TYPE2 **Life){
  int i,j;
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){
    if(Life[i][j].val!=0) 
      if(genrand_real2() < death) Life[i][j].val=0;
  }
}

void Update(void)
{ 
  Display(Life);
  Plot(1,Life);
  //Asynchronous();
  Synchronous(1,Life);
  Death(Life);
  //while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement
  
}
