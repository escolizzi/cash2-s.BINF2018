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
n=3;

int WrapNei(pos,which)
{ 
  if(which==0){ 
    if(pos>0 || pos<=nrow) return pos;
    if(pos<=0) pos +=nrow;
    if(pos>mrow) pos -=
  } 
}

int CountNeiN(Life,n,row,col)
{
  int i,j,k,sum=0;
  for(k=1;k<=n;k++){
    //row on top: i=row-k, j from col-k to col+k incl.
    thisrow=WrapNei(row-k,0);     //WrapNei return the position corrected for boundary conditions (first argument is pos number, second is 0 or 1 for row or col)
    for(j=col-k;j<=col+k;j++){
      thiscol=WrapNei(j,1);
      if(Life[thisrow][thiscol].val == 1) sum++;
    }
    //row below
    thisrow=WrapNei(row+k,0);
    for(j=col-k;j<=col+k;j++){
      thiscol=WrapNei(j,1);
      if(Life[thisrow][thiscol].val == 1) sum++;
    }
    //left column -- carful not to recount
    
    //right column
  }
}

void Initial(void)
{
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 200; /* # of row (default=100)*/
  ncol = 200; /* # of column (default=100)*/
  nplane = 1; /* # of planes (default=0)*/
  scale = 2; /* size of the window (default=2)*/
  boundary = WRAP; /* the type of boundary: FIXED, WRAP, ECHO (default=WRAP). Note that
		      Margolus diffusion is not supported for ECHO. */
  ulseedG = 3452354; /* random seed (default=56)*/

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
  ReadSavedData("glidergun.sav",1,Life);
  printf("\n\nLife is ready. Click or scroll your mouse to update the CA.\n\n");
  Boundaries2(Life);
}

double thresh=0.999;
void NextState(int row,int col)
{
  int sum;

  //sum = CountMoore8(Life,1,row,col);
  sum = CountNeiN(Life,n,row,col);
  
  double rn = genrand_real1();
  if(Life[row][col].val == 0){ 
    if( sum==3 && rn< thresh)
      Life[row][col].val = 1;
  }else if (Life[row][col].val == 1){
    if( (sum<2 || sum>3) && rn< thresh ) 
      Life[row][col].val = 0;
  }
    
}

void Update(void)
{ 
  Display(Life);
  //Asynchronous();
  Synchronous(1,Life);

  //while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement
  
}
