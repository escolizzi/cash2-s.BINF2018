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
double a=1.;    //growth rate of strain 1
double b=1.;   //growth rate of strain 2

void Initial(void)
{
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 200; /* # of row (default=100)*/
  ncol = 200; /* # of column (default=100)*/
  nplane = 1; /* # of planes (default=0)*/
  scale = 4; /* size of the window (default=2)*/
  boundary = FIXED; /* the type of boundary: FIXED, WRAP, ECHO (default=WRAP). Note that
		      Margolus diffusion is not supported for ECHO. */
  ulseedG = 47329; /* random seed (default=56)*/

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
  //InitialSet(Life,1,0,1,0.3);
  //ReadSavedData("glidergun.sav",1,Life);
  //printf("\n\nLife is ready. Click or scroll your mouse to update the CA.\n\n");
  
  InitialSet(Life,2,0);
  int i,j;
  for(i=1;i<=nrow;i++) for(j=1;j<=ncol;j++){
    Life[i][j].val=0;
    if(i<nrow/2+10 && i>nrow/2-10 && j>ncol/2-10 && j<ncol/2+10){
      if(genrand_real1()<0.5) Life[i][j].val =1;
      else Life[i][j].val =2;
    }
    //if(j<10) Life[i][j].val =1;
    //if(i<=10 && j<10) Life[i][j].val =1;
    //else if(i>10 && j<10) Life[i][j].val =2;
    //else Life[i][j].val=0;
  }
  Boundaries2(Life);
}

void NextState(int row,int col)
{
  int sum1,sum2;
  double fsum1,fsum2;
  double rn;
  
  if(Life[row][col].val==0){
    //printf("Enter life[row][col]=0\n");
    
    sum1 = CountMoore8(Life,1,row,col);
    sum2 = CountMoore8(Life,2,row,col);
    
    //printf("sum1=%d,sum2=%d\n",sum1,sum2);
    
    if(sum1+sum2 > 0){
      fsum1 = a * (double)sum1;
      fsum2 = b * (double)sum2;
      //ftot=(double)tot;
      
      rn=8.*genrand_real2();
      
      //printf("fit1=%f,fit2=%f, rn=%f\n",fsum1,fsum2, rn);
      
      if(rn<fsum1) Life[row][col].val=1;
      else if( rn< (fsum1+fsum2) ) {
        Life[row][col].val=2;
        
        //if(sum2==0 || row<nrow/4) {
        //  printf("We got this Life[%d][%d].val=%d\n\n",row,col,Life[row][col].val);
        //  printf("This is error\n");
        //  exit(1);
        // }
      }
      //printf("We got this Life[%d][%d].val=%d\n\n",row,col,Life[row][col].val);
      
    }
  }else{
    //tiny movement into empty space to get things less stuck
    ;/*
    int neirow, neicol,randir;
    randir = 1+(int)(8.*genrand_real2());
    GetNeighborC(Life,row,col,randir ,&neirow,&neicol);
    
    
    if(Life[neirow][neicol].val==0 && genrand_real2()<0.00){
      Life[neirow][neicol].val=Life[row][col].val;
      Life[row][col].val=0;
    }
    */
  }
  
}

void Update(void)
{ 
  Display(Life);
  Asynchronous();
  //Synchronous(1,Life);
  int i;
  
  
  for(i=1;i<=ncol;i++){
    if(Life[nrow][i].val==1){
      printf("Overtaking ended: ncol=%d, Time=%d, speed=%f\n",i-10,Time, (i-10)/(double)Time);
      while( Mouse()==0) {};
      exit(0);
    }
    
  }
  /*
  for(i=1;i<=nrow;i++){
    if(Life[i][ncol].val!=0){
      printf("Simulation ended: ncol=%d, Time=%d, speed=%f\n",ncol-10,Time, (ncol-10)/(double)Time);
      exit(0);
    }
    
  }
  */
  //while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement
  
}


