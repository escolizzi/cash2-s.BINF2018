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


static TYPE2** world;

double kdiff = 0.1;
double beta=2.;

void Initial(void)
{
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 100; /* # of row (default=100)*/
  ncol = 100; /* # of column (default=100)*/
  nplane = 1; /* # of planes (default=0)*/
  scale = 8; /* size of the window (default=2)*/
  boundary = WRAP; /* the type of boundary: FIXED, WRAP, ECHO (default=WRAP). Note that
		      Margolus diffusion is not supported for ECHO. */
  ulseedG = 56; /* random seed (default=56)*/

  /* useally, one does not have to change the followings */
  /* the value of boundary (default=(TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.})*/
  boundaryvalue2 = (TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.}; 
}

void InitialPlane(void)
{
  MakePlane(&world);

  /* InitialSet(1,2,3,4,5)
    1: name of plane
    2: number of the state other than the background state
    3: background state or empty state
    4: state I want to put
    5: fraction of cells that get S1 state (0 to 1)
  */
  int i,j;
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){
    world[i][j].val=0;
    if(pow(i-nrow/2,2)+pow(j-ncol/2,2) < 100)
      world[i][j].val=1;
  }
  
  //InitialSet(world,1,0,1,0.1);
  //ReadSavedData("glidergun.sav",1,world);
  //printf("\n\nworld is ready. Click or scroll your mouse to update the CA.\n\n");
  //Boundaries2(world);
}

void NextState(int row,int col)
{
  int s_sum,n_sum,nr,nc,rn;
  
  //TYPE2 *nei;
  
  rn=1+(int)(8.*genrand_real2());
  GetNeighborC(world,row,col,rn,&nr,&nc);
  
  s_sum = CountMoore8(world,1,row,col);
  n_sum = CountMoore8(world,1,nr,nc);
  
  //if(world[row][col].val==world[nr][nc].val){
  //  ;
  //}
  
  if( world[row][col].val ==0 && world[nr][nc].val==1 && s_sum==1 && n_sum==0 ){
    if(genrand_real2()< kdiff/2.){
      world[row][col].val =1;
      world[nr][nc].val =0;
    }
  }else if(world[row][col].val ==0 && world[nr][nc].val==1 && n_sum>0){
    if(s_sum > n_sum){
      if(genrand_real2()< kdiff/2.){
        world[row][col].val =1;
        world[nr][nc].val =0;
      }
    }else{
      if(genrand_real2()< kdiff/2. * exp(-beta*( n_sum - (s_sum-1) ))){
        world[row][col].val =1;
        world[nr][nc].val =0;
      }
    }
  }
  
  
  
  /*
   *   empty [0], empty [0] -> do nothing
   *   empty [0], loose lipid [0] -> free diffusion
   *   empty [n], lipid [m] -> if n>=m P_swap=1, else P_swap = exp(-beta*DG), DG = m-n
   *   lipid [n], lipid [m] -> for now do nothing, actually lipids with few neighs should move towards more neighs
   *   
  */
  
  
  //if(sum==3 || (sum==2 && world[row][col].val==1))
  //  world[row][col].val = 1;
  //else
  //  world[row][col].val = 0;
  
}

void Update(void)
{ 
  Display(world);
  Asynchronous();
  //Synchronous(1,world);

  //while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement
  
}
