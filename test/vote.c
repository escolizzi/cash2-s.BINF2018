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

static TYPE2** Vote;
void Initial(void)
{
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 100; /* # of row (default=100)*/
  ncol = 100; /* # of column (default=100)*/
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
  MakePlane(&Vote);

  /* InitialSet(1,2,3,4,5)
    1: name of plane
    2: number o2f the state other than the background state
    3: background state or empty state
    4: state I want to put
    5: fraction of cells that get S1 state (0 to 1)
  */
  InitialSet(Vote,1,0,1,0.5);
  printf("\n\nVote is ready. Click or scroll your mouse to update the CA.\n\n");
  Boundaries2(Vote);
}

void NextState(int row,int col)
{
  int sum;

  sum = CountMoore8(Vote,1,row,col);

  if(genrand_real1()<0.75){
    Vote[row][col].val = 1-Vote[row][col].val ;
  }
  else {
  if(sum >= 0 && sum <= 3){
    Vote[row][col].val = 0;
  }
  else if(sum >= 5){
    Vote[row][col].val = 1;
  }
  }
  /*
    genrand_real1() returns a floating number between 0 and 1
    (both 0 and 1 inclusive) with a uniform probability
    distribution. For example,

    if(genrand_real1()<0.1){
       case1();
    }else{
       case2();
    }

    case1() will be executed with a chance of 10% whicle
    case2() with 90%.
  */
}

void Update(void)
{
  Display(Vote);
  Synchronous(1,Vote);
  //while( Mouse()==0) {};
  if(Time%50 ==0) printf("Time: %d\n",Time);
  //Asynchronous();
}
