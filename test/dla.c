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

static TYPE2** Free;
static TYPE2** Fixed;

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
  MakePlane(&Free,&Fixed);

  /* InitialSet(1,2,3,4,5)
    1: name of plane
    2: number of the state other than the background state
    3: background state or empty state
    4: state I want to put
    5: fraction of cells that get S1 state (0 to 1)
  */
  InitialSet(Free,1,0,1,0.1);
  InitialSet(Fixed,0,0);

  Fixed[100][100].val = 1;

  Boundaries2(Free);
  Boundaries2(Fixed);
}

void NextState(int row,int col)
{
  if(Fixed[row][col].val == 0){
    if(Free[row][col].val == 1 && CountMoore8(Fixed,1,row,col) > 0){
      Free[row][col].val = 0;
      Fixed[row][col].val = 1;
    }
  }else{
    ;
  }
}

void Update(void)
{
  Display(Free,Fixed);
  Synchronous(2,Free,Fixed);
  MDiffusion(Free);
}
