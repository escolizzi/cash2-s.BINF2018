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
  nrow = 100; /* # of row (default=100)*/
  ncol = 100; /* # of column (default=100)*/
  nplane = 1; /* # of planes (default=0)*/
  scale = 8; /* size of the window (default=2)*/
  boundary = WRAP; /* the type of boundary: FIXED, WRAP, ECHO (default=WRAP). Note that
		      Margolus diffusion is not supported for ECHO. */
  ulseedG = 156; /* random seed (default=56)*/

  /* useally, one does not have to change the followings */
  /* the value of boundary (default=(TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.})*/
  boundaryvalue2 = (TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.}; 
}

void InitialPlane(void)
{
  int i,j;
  MakePlane(&Life);

  /* InitialSet(1,2,3,4,5)
    1: name of plane
    2: number of the state other than the background state
    3: background state or empty state
    4: state I want to put
    5: fraction of cells that get S1 state (0 to 1)
  */
  for(i=1; i<=nrow; i++) for(j=1; j<=ncol; j++){
      Life[i][j].fval=0;
      if(genrand_real1()<.99){// && i>nrow/2-25 && i<nrow/2+25 &&  j>ncol/2-25 && j<ncol/2+25){
        Life[i][j].fval=0.+8*genrand_real1();
        Life[i][j].val=Life[i][j].fval/10;
        //if(j<ncol/2) bact[i][j].fval=0.;
       }
  }
  //InitialSet(Life,1,0,1,0.3);
  //ReadSavedData("glidergun.sav",1,Life);
  //printf("\n\nLife is ready. Click or scroll your mouse to update the CA.\n\n");
  Boundaries2(Life);
  
  for(i=0;i<=15;i++){
    
    int r=i* 255/15;
    int g=r;
    int b=r;
    ColorRGB(i,r,g,b);	/*so colour indexes in the gradient start at 10*/
  }
  
}

double sigma=100.;
double delta=0.0;//1.4;
double h=25.;
double i=1.;

void NextState(int row,int col)
{
  double numerator,denominator;
  int k;
  TYPE2 *nei;
  
  numerator = sigma*pow( delta+Life[row][col].fval ,2.);
  denominator = pow(h,2.) + pow( delta+Life[row][col].fval ,2.);
  
  for(k=1;k<=8;k++){
    nei=GetNeighborP(Life,row,col,k);
    denominator += i* pow(nei->fval,2.);
  }
  Life[row][col].fval= numerator/denominator;
  Life[row][col].val=Life[row][col].fval/10;
}

void Update(void)
{ 
  int i,j;
  double avrg=0.;
  TYPE2 *nei;
  
  Display(Life);
  
  for(i=1; i<=nrow; i++) for(j=1; j<=ncol; j++){
    avrg+=Life[i][j].fval;
    if(Life[i][j].fval>=0) printf("this: %f\n",Life[i][j].fval);
  }
  printf("%f\n", avrg / (double)((nrow-1)*(ncol-1)) );
  
  //REPLICATION - works quite stably 
  for(i=1; i<=nrow; i++) for(j=1; j<=ncol; j++){
    if(Time%100==0 && genrand_real1()<0.5 && Life[i][j].fval>70){
      int rn = 1+ (int)( 8.*genrand_real2() );
      nei=GetNeighborP(Life,i,j,rn);
      nei->fval=Life[i][j].fval/2.;
      Life[i][j].fval/=2.;Life[i][j].fval++;
    }
  }
  
  //Asynchronous();
  Synchronous(1,Life);
  //if( /*Time%1==0 &&*/ Time>50){
  //  for(i=1; i<=nrow; i++) for(j=1; j<=ncol; j++){
  //    
  //    //Life[i][j].fval=0;
  //    if(genrand_real1()<1.){// && i>nrow/2-25 && i<nrow/2+25 &&  j>ncol/2-25 && j<ncol/2+25){
  //      Life[i][j].fval+=10.*genrand_real1();
  //      //Life[i][j].val=Life[i][j].fval/10;
  //     //if(j<ncol/2) bact[i][j].fval=0.;
  ////     }
  ////  }
  //}
  //sleep(0.5);
  while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement
  
}
