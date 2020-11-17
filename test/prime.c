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

static TYPE2** Prime;

void Initial(void)
{
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 400; /* # of row (default=100)*/
  ncol = 400; /* # of column (default=100)*/
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
  MakePlane(&Prime);

  /* InitialSet(1,2,3,4,5)
    1: name of plane
    2: number o2f the state other than the background state
    3: background state or empty state
    4: state I want to put
    5: fraction of cells that get S1 state (0 to 1)
  */
  //InitialSet(Prime,1,0,1,0.1);
  ReadSavedData("big_elephant.sav",1,Prime);
  printf("\n\nPrime is ready. Click or scroll your mouse to update the CA.\n\n");
  Boundaries2(Prime);
}

void NextState(int row,int col)
{
  Prime[row][col].val = SumNeumann4(Prime,row,col)%2;
}

void Update(void)
{
  Display(Prime);

  Synchronous(1,Prime);
  while( Mouse()==0) {};
  //Asynchronous();
}
