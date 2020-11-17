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
      Life[i][j].val=0;
      if(genrand_real1()<1){// && i>nrow/2-25 && i<nrow/2+25 &&  j>ncol/2-25 && j<ncol/2+25){
        Life[i][j].val=80+(genrand_real1()-0.5)*80/2;
        //Life[i][j].val=Life[i][j].fval/10;
        //if(j<ncol/2) bact[i][j].fval=0.;
       }
  }
  //InitialSet(Life,1,0,1,0.3);
  //ReadSavedData("glidergun.sav",1,Life);
  //printf("\n\nLife is ready. Click or scroll your mouse to update the CA.\n\n");
  Boundaries2(Life);
  
  for(i=0;i<=100;i++){
    
    int r=i* 255/100;
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
  
  numerator = sigma*pow( delta+ (double)(Life[row][col].val) ,2.);
  denominator = pow(h,2.) + pow( delta+(double)(Life[row][col].val) ,2.);
  
  for(k=1;k<=8;k++){
    nei=GetNeighborP(Life,row,col,k);
    denominator += i* pow((double)(nei->val),2.);
  }
  
  double fval= numerator/denominator;
  int val=(int)(numerator/denominator);
  double diff= fval - (double)val;
  
  //if(Life[row][col].val==0 || Life[row][col].val==1 ){
  //  printf("This is %d, num/den= %f, val %d, diff = %f\n",Life[row][col].val, numerator/denominator, val,diff);
  //}
  
  Life[row][col].val = val;
  double rn=genrand_real2();
  if(rn<diff) Life[row][col].val += 1;
  
  //printf(" because of rn=%f, I now am %d\n",rn,Life[row][col].val);
  //if(val>0 && diff<genrand_real2()) 
  //  Life[row][col].val += 1;
  
  //Life[row][col].val=Life[row][col].fval/10;
}

void Update(void)
{ 
  int i,j;
  double avrg=0.;
  TYPE2 *nei;
  
  Display(Life);
  
  if(Time%10==0) printf("Time: %d\n",Time);
  
  for(i=1; i<=nrow; i++) for(j=1; j<=ncol; j++){
    avrg+=Life[i][j].val;
    if(Life[i][j].val>=0) printf("this: %d\n",Life[i][j].val);
  }
  printf("%f\n", avrg / (double)((nrow-1)*(ncol-1)) );
  
  //REPLICATION - works quite stably 
  /*
  for(i=1; i<=nrow; i++) for(j=1; j<=ncol; j++){
    if(Time%50==0 && genrand_real1()<1.5 && Life[i][j].val>60){
      int rn = 1+ (int)( 8.*genrand_real2() );
      nei=GetNeighborP(Life,i,j,rn);
      nei->fval=Life[i][j].val/2;
      Life[i][j].val/=2;Life[i][j].val++;
    }
  }
  */
  //Asynchronous();
  Synchronous(1,Life);
  /*
  if( Time%50==0 && Time>0){
    for(i=1; i<=nrow; i++) for(j=1; j<=ncol; j++){
      
      //Life[i][j].fval=0;
      if(genrand_real1()<1.){// && i>nrow/2-25 && i<nrow/2+25 &&  j>ncol/2-25 && j<ncol/2+25){
        Life[i][j].val+=30.*genrand_real1();
        //Life[i][j].val=Life[i][j].fval/10;
       //if(j<ncol/2) bact[i][j].fval=0.;
      }
    }
  }
  */
  //sleep(0.5);
  while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement
  
}
