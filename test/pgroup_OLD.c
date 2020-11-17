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
double kdeath = 0.1;
double kmut = 0.01;
double init_fval = 0.75;
double kcat=1.; // same for everybody
void Initial(void)
{
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 200; /* # of row (default=100)*/
  ncol = 200; /* # of column (default=100)*/
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
  MakePlane(&world);

  /* InitialSet(1,2,3,4,5)
    1: name of plane
    2: number of the state other than the background state
    3: background state or empty state
    4: state I want to put
    5: fraction of cells that get S1 state (0 to 1)
  */
  //InitialSet(world,1,0,1,0.3);
  //ReadSavedData("glidergun.sav",1,world);
  //printf("\n\nworld is ready. Click or scroll your mouse to update the CA.\n\n");
  
  int i,j;
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){
    world[i][j].val=0;
    world[i][j].fval=-1.;
    //if(i>nrow/2 -10 && i<nrow/2 + 10 && j>ncol/2 -10 && j<ncol/2 + 10 && genrand_real1()<0.25){
      world[i][j].val=1;
      world[i][j].fval=init_fval;
    //}
  }
  
  Boundaries2(world);
}


int Replication(TYPE2 **world, int r, int c)
{
  int rr,rc,rn,outcome=0;
  TYPE2 *icel,*nei,*parent;
  icel = &world[r][c];
  //competition for replication on focal site
  RandomMooreC8(world,r,c,&rr,&rc);  //get neigh
  nei=&world[rr][rc];
  if(nei->val!=0 && nei->crow!=-1){
    if(genrand_real2() < kcat/2.){
        //replication happens
        parent = &world[nei->crow][nei->ccol];
        
        icel->val=1;
        icel->fval = parent->fval;
        icel->crow = -1;
        icel->ccol = -1;
        
        //break complex;
        world[nei->crow][nei->ccol].crow = -1;
        world[nei->crow][nei->ccol].ccol = -1;
        nei->crow = -1;
        nei->ccol = -1;
        
        outcome=1;
    }
  }
  //else outcome=0;
  
  return outcome;
}

void Mutations(TYPE2 *icel)
{
  if(genrand_real2()<kmut)
    icel->fval += ( genrand_real1() - 0.5 )/ 20.;
  
  if(icel->fval < 0 ) icel->fval = -icel->fval;
  if(icel->fval > 1.) icel->fval = 2. - icel->fval;
}

void NextState(int row,int col)
{
  int outcome=-1,rr,rc;
  TYPE2 *nei;
  if(world[row][col].val==0){           //this will have to change when vescicles are implemented 
    outcome = Replication(world,row,col);
    if(outcome==1) Mutations(&world[row][col]);
  }
  else{
    //pick random neigh
    if(world[row][col].crow==-1){
      RandomMooreC8(world,row,col,&rr,&rc);  //get neigh
      nei = &world[rr][rc]; 
      if(nei->val!=0){
        //complex formation, guy at i,j is the replicator
        if(genrand_real2() < world[row][col].fval){
          world[row][col].crow=rr;
          world[row][col].ccol=rc;
          nei->crow = row;
          nei->ccol = col;
        }
      }
    }
  }
}

void Death(TYPE2 **world){
  int i,j;
  for(i=1; i<=nrow;i++)for(j=1; j<=ncol;j++){
    if(genrand_real2()<kdeath){
      world[i][j].val=0;
      world[i][j].fval=0.;
      world[i][j].crow=-1;
    }
  }
    
}

void PrintStuff(TYPE2**world)
{
  double avfval=0.,minfval=2.,maxfval=-1;
  int popsize=0;
  int i,j;
  for(i=1; i<=nrow;i++)for(j=1; j<=ncol;j++){
    if(world[i][j].val != 0){
      if(world[i][j].fval<minfval) minfval=world[i][j].fval;
      if(world[i][j].fval>maxfval) maxfval=world[i][j].fval;
      avfval += world[i][j].fval;
      popsize++;
    }
  }
  printf("Av fval = %f, min fval = %f, max fval = %f, popsize = %d\n", avfval/(float)popsize,minfval,maxfval, popsize);
  if(popsize==0){
    fprintf(stderr,"Population extinct, ending simulation\n");
    exit(1);
  }
}

void Update(void)
{ 
  Display(world);
  Asynchronous();
  //Synchronous(1,world);
  Death(world);
  //while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement
  if(Time % 100 == 0 ) 
    PrintStuff(world);    // also ends sim if popisze == 0
}


