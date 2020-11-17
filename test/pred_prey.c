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
static TYPE2** food;
double foodinflux;
double NONE = 0.5;

void Initial(void)
{
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 200; /* # of row (default=100)*/
  ncol = 200; /* # of column (default=100)*/
  nplane = 2; /* # of planes (default=0)*/
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
  MakePlane(&world,&food);
  printf("Hello0\n");

  /* InitialSet(1,2,3,4,5)
    1: name of plane
    2: number of the state other than the background state
    3: background state or empty state
    4: state I want to put
    5: fraction of cells that get S1 state (0 to 1)
  */
  InitialSet(world,1,0,1,0.3);
  printf("Hello1\n");
  int i,j;
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){
    if(genrand_real1()<0.3) food[i][j].val=1;
  }
  
  //InitialSet(food,1,0,1,0.3);
  printf("Hello2\n");
  //ReadSavedData("glidergun.sav",1,world);
  //printf("\n\nworld is ready. Click or scroll your mouse to update the CA.\n\n");
  //Boundaries2(world);
}

void NextState(int row,int col)
{
  int i,who,sum=0;
  TYPE2 *nei;
  TYPE2 *arrnei[8];
  if(world[row][col].val==0){
    //replication can happen here
    for(i=1;i<=8;i++){
      //grab the guy, irrespective of what it is, if it has eaten, it has a chance of reproducing
      nei = GetNeighborP(world, row,col,i);
      if(nei->val!=0 && nei->val2==1){
        arrnei[sum]=nei;
        sum++;
      }
    }
    if( sum>0 && genrand_real2()< sum/(8. + NONE) ){
      who = (int)( sum * genrand_real2() ) ;
      world[row][col]=*nei;
    }
    
  }else if(world[row][col].val==PREY){
    //prey can eat from food plane
    if( world[row][col].val2==1 ) ... decay of food inside prey
    if( food[row][col].val>0 && world[row][col].val2==0 ){
      ... eating of prey
    }
    
    
    //predator can eat if in neighbour
    ...
    
  }
  
  
  
}

void Update(void)
{ 
  Display(world,food);
  Asynchronous();
  //Synchronous(1,world);

  //while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement
  
}
