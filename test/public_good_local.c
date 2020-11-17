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
// double a=1.;    //growth rate of strain 1
// double b=0.25;   //growth rate of strain 2

//    THIS IS A NICE PARAMETER COMBINATION  10,2,0.4,0.1,0.,10  //
double B = 10.;
double C = 2.0;
double d = 0.4;
double beta=0.1;//0.005;
double kdiff=0.0;
double publicgood_pp=10.;

//double NONE=30;

int CountLargerMoore(TYPE2** Life,int val ,int row,int col);
void DiffBySwap(TYPE2 **bact,int row, int col);

void Initial(void)
{
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 300; /* # of row (default=100)*/
  ncol = 302; /* # of column (default=100)*/
  nplane = 1; /* # of planes (default=0)*/
  scale = 2; /* size of the window (default=2)*/
  boundary = WRAP; /* the type of boundary: FIXED, WRAP, ECHO (default=WRAP). Note that
		      Margolus diffusion is not supported for ECHO. */
  ulseedG = 363; /* random seed (default=56)*/
  display=1;
  /* useally, one does not have to change the followings */
  /* the value of boundary (default=(TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.})*/
  boundaryvalue2 = (TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.}; 
}

void InitialPlane(void)
{
  MakePlane(&Life);
  if(kdiff+d >=1.){fprintf(stderr,"InitialPlane(): Error. kdiff+d >=1.\n");exit(1);}
  /* InitialSet(1,2,3,4,5)
    1: name of plane
    2: number of the state other than the background state
    3: background state or empty state
    4: state I want to put
    5: fraction of cells that get S1 state (0 to 1)
  */
//    InitialSet(Life,2,0,1,0.1,2,0.05);
  //ReadSavedData("glidergun.sav",1,Life);
  //printf("\n\nLife is ready. Click or scroll your mouse to update the CA.\n\n");
  
  InitialSet(Life,2,0);
  int i,j;
  for(i=1;i<=nrow;i++) for(j=1;j<=ncol;j++){
    Life[i][j].val=0;
    if(i<20 || i>80) continue; 
    if(j>=2 && j<20) Life[i][j].val =1;
    if(j==1) Life[i][j].val =2;
    //if(i==nrow/2){
      //if( j<ncol/2+2 && j>ncol/2-2){ 
    //   Life[i][j].val=1;
      //}
    //}if(i==nrow/2+1){
    //  if( j==ncol/2){ 
    //    Life[i][j].val=2;
    //  }
    //}
    
    
    //if(i>nrow/2-25 && i<nrow/2+25 && j>ncol/2-25 && j<ncol/2+25){ 
    //  Life[i][j].val =1;
    //  if(j<ncol/2 ) Life[i][j].val =2;
    //
    //}
  }
  Boundaries2(Life);
}

int Checkboundary(int i, char r_or_c){
  if( i>=1 && i<=nrow ) return i;
  else{
    if(boundary==FIXED) return -1;
    if(boundary==WRAP) return (r_or_c=='r')? ((i+nrow-1)%nrow)+1 : ((i+ncol-1)%ncol)+1;
  }
}

int CountLargerMoore(TYPE2** Life,int val ,int row,int col)
{
  //let's do only the 5x5 moore neigh.
  int i,j,counter=0;
  for(i=row-2; i<=row+2;i++)for(j=col-2; j<=col+2;j++){
    //printf("old i %d, j %d\n",i,j);
    int ii=Checkboundary(i,'r');
    int jj=Checkboundary(j,'c');
    //printf("new i %d, j %d\n",ii,jj);
    if( ii==row && jj==col) continue;
    if( (ii<0 || jj<0) && boundary==FIXED ) continue;
    if(Life[ii][jj].val==val) counter++;
  }
  return counter;
}

void DiffBySwap(TYPE2 **bact,int row, int col)
{
  int rnei;
  TYPE2 *nei; 
  TYPE2 tmp;
  
  rnei = 1 + (int)( 8*genrand_real2() );
  nei=GetNeighborP(bact,row,col,rnei);
  
  tmp=bact[row][col];
  bact[row][col]=*nei;
  *nei=tmp;
  
}

void NextState(int row,int col)
{
  if(col<=1) return;
  
  int sum1,sum2;
  double fsum1,fsum2;
  double rn;
  
  if(Life[row][col].val==0){
    //printf("Enter life[row][col]=0\n");
    //printf("row %d, col %d\n",row,col);
    //sum1 = CountLargerMoore(Life,1,row,col); //this is also the total amount of public good
    //sum2 = CountLargerMoore(Life,2,row,col);
    
    sum1 = CountMoore8(Life,1,row,col); //this is also the total amount of public good
    sum2 = CountMoore8(Life,2,row,col);
    
    
    if(sum1==0 || sum1+sum2==0) return;
    //double percapita_pg = publicgood_pp*sum1/9.;    // / (double)(sum1+sum2); <- No, public good diffuses also to empty space
    double totfitness_coop = (B*publicgood_pp*sum1/9. - C*publicgood_pp)*sum1;
    totfitness_coop = (totfitness_coop>0.)?totfitness_coop:0.;
    
    double totfitness_cheat = (B*publicgood_pp*sum1/9.)*sum2;
    double totalfitness = totfitness_coop + totfitness_cheat;
    
    rn=genrand_real2() * (totalfitness)/( 1. - exp(-beta*totalfitness) );
    if(rn<totfitness_coop) 
      Life[row][col].val=1;
    else if (rn<totfitness_coop+totfitness_cheat) Life[row][col].val=2;
    
  }else{
    rn=genrand_real2();
    if(rn < d ) Life[row][col].val=0;
    else if(rn < d+kdiff) DiffBySwap(Life,row,col);
    
  }
  
}

void Update(void)
{ 
  Display(Life);
  //Asynchronous();
  Synchronous(1,Life);
  int i,j;
  int totpop=0, totcoop=0,totcheat=0;
  
//   if(Time%1000==0 && Time !=0){
//     C+=1.;
//     fprintf(stderr,"Time=%d, C=%f\n",Time,C);
//   }
  
  static int halfwayyet=0;
  for(i=1;i<=nrow;i++){
    if(Life[i][ncol/2].val!=0 && !halfwayyet){
      printf("Half pass, time=%d\n",Time);
      halfwayyet=1;
    }
    if(Life[i][ncol].val!=0){
      printf("Tot pass, time=%d\n",Time);
      exit(0);
    }
  }
  
  //while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement
  
  //fprintf(stderr,"Time=%d\n",Time);
  
  if(Time%100==0){
    //fprintf(stderr,"Time=%d\n",Time);
    for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){
      if(Life[i][j].val>0){
        totpop++;
        if(Life[nrow][i].val==1) totcoop++;
        else totcheat++;
      }
    }

    if(totpop==0) {
      fprintf(stderr,"totpop=%d, Population extinct\n",totpop);
      exit(1);
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


