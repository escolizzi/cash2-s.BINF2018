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

#define MAXNUMBER 19683

static TYPE2** Life;
// double a=1.;    //growth rate of strain 1
// double b=0.25;   //growth rate of strain 2

double B = 10;
double C = 2;
double d = 0.5;
double NONE=10;

int CountLargerMoore(TYPE2** Life,int val ,int row,int col);

void Initial(void)
{
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 400; /* # of row (default=100)*/
  ncol = 400; /* # of column (default=100)*/
  nplane = 1; /* # of planes (default=0)*/
  scale = 1; /* size of the window (default=2)*/
  boundary = WRAP; /* the type of boundary: FIXED, WRAP, ECHO (default=WRAP). Note that
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
   InitialSet(Life,2,0,1,0.4,2,0.1);
  //ReadSavedData("glidergun.sav",1,Life);
  //printf("\n\nLife is ready. Click or scroll your mouse to update the CA.\n\n");
  
//   InitialSet(Life,2,0);
//   int i,j;
//   for(i=1;i<=nrow;i++) for(j=1;j<=ncol;j++){
//     //if(j<10) Life[i][j].val =1;
//     if(i>nrow/2-50 && i<nrow/2+50 && j>ncol/2-50 && j<ncol/2+50) Life[i][j].val =1;
//     else if(i==nrow/2-20 ) Life[i][j].val =2;
//     else Life[i][j].val=0;
//   }
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

int Transition(int state)
{
  //Life[row][col].val is a number between 0 and 3^9 representing a neighbourhood
  int config[9]={0};   // 9 is neighbourhood size, includng self
  int comp[3]={0};
  int number;
  int pos=0;
  
  if(state==0) return 0; //empty neigh stays empty
  number=state;
  
  //get neighbourhood from state
  while(number>0){
    int which = number%3;
    config[pos] = which;
    comp[which] ++ ;
    number /= 3;
    pos+=1;
  }
  //now config has full information over the neighbourhood
  // note by the way that this is overhead, because we only care about 
  // number of coop, cheat and empty, i.e. comp
  
  
}

void NextState(int row,int col)
{
  int sum1,sum2;
  double fsum1,fsum2;
  double rn;
  
  fprintf(stderr,"UNDER CONSTRUCTION\n");
  exit(1);
  
  Life[row][col].val=Transition(Life[row][col].val);
  
//   if(Life[row][col].val==0){
//     //printf("Enter life[row][col]=0\n");
//     //printf("row %d, col %d\n",row,col);
//     sum1 = CountLargerMoore(Life,1,row,col); //this is also the total amount of public good
//     sum2 = CountLargerMoore(Life,2,row,col);
//     
//     if(sum1==0 || sum1+sum2==0) return;
//     double percapita_pg = sum1/8.;    // / (double)(sum1+sum2); <- No, public good diffuses also to empty space
//     double totfitness_coop = (B*percapita_pg - C)*sum1;
//     totfitness_coop = (totfitness_coop>0.)?totfitness_coop:0.;
//     double totfitness_cheat = B*percapita_pg*sum2;
//     double totalfitness = totfitness_coop + totfitness_cheat;
//     
//     rn=genrand_real2()*(totalfitness+NONE);
//     if(rn<totfitness_coop) Life[row][col].val=1;
//     else if (rn<totfitness_coop+totfitness_cheat) Life[row][col].val=2;
//     
//     //printf("row %d, col %d, sum1=%d,sum2=%d\n",row,col,sum1,sum2);
//     
// //     if(sum1+sum2 > 0){
// //       fsum1 = a * (double)sum1;
// //       fsum2 = b * (double)sum2;
// //       //ftot=(double)tot;
// //       
// //       rn=8.*genrand_real2();
// //       
// //       //printf("fit1=%f,fit2=%f, rn=%f\n",fsum1,fsum2, rn);
// //       
// //       if(rn<fsum1) Life[row][col].val=1;
// //       else if( rn< (fsum1+fsum2) ) {
// //         Life[row][col].val=2;
// //         
// //         //if(sum2==0 || row<nrow/4) {
// //         //  printf("We got this Life[%d][%d].val=%d\n\n",row,col,Life[row][col].val);
// //         //  printf("This is error\n");
// //         //  exit(1);
// //         // }
// //       }
// //       //printf("We got this Life[%d][%d].val=%d\n\n",row,col,Life[row][col].val);
// //       
// //     }
//   }else{
//     //tiny movement into empty space to get things less stuck
//     if(genrand_real2()<d) Life[row][col].val=0;; /*
//     int neirow, neicol,randir;
//     randir = 1+(int)(8.*genrand_real2());
//     GetNeighborC(Life,row,col,randir ,&neirow,&neicol);
//     
//     
//     if(Life[neirow][neicol].val==0 && genrand_real2()<0.00){
//       Life[neirow][neicol].val=Life[row][col].val;
//       Life[row][col].val=0;
//     }
//     */
//   }
//   
}

void Update(void)
{ 
  Display(Life);
  Asynchronous();
  //Synchronous(1,Life);
  int i,j;
  int totpop=0, totcoop=0,totcheat=0;
  
  if(Time%300==0 && Time !=0){
    C+=1.;
    fprintf(stderr,"Time=%d, C=%f\n",Time,C);
  }
  if(Time%100==0){
    fprintf(stderr,"Time=%d\n",Time);
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


