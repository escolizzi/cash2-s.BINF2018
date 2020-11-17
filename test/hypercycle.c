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

#define MYBUFFER 50

static TYPE2** world;
int nspecies=6;
int hyperc1=6;
int hyperc2=0;
double uncat_birth=0.06;
double cat_birth=0.1; //7*0.1 = 0.7 who cares
double cat_birth1=0.08; //7*0.08 = 0.56
double cat_birth2=0.12; //7*0.12 = 0.84 who cares
double death = 0.03;
//double death1 = 0.022;
//double death2 = 0.015;

//static TYPE2** food;
//double foodinflux;
//double NONE = 0.5;

void Initial(void)
{
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 2048; /* # of row (default=100)*/
  ncol = 2048; /* # of column (default=100)*/
  nplane = 1; /* # of planes (default=0)*/
  scale = 4; /* size of the window (default=2)*/
  boundary = FIXED; /* the type of boundary: FIXED, WRAP, ECHO (default=WRAP). Note that
		      Margolus diffusion is not supported for ECHO. */
  ulseedG = 324798; /* random seed (default=56)*/
  display = 1;
  /* usually, one does not have to change the followings */
  /* the value of boundary (default=(TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.})*/
  boundaryvalue2 = (TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.}; 
}

void InitialPlane(void)
{
  MakePlane(&world);
  //printf("Hello0\n");

  /* InitialSet(1,2,3,4,5)
    1: name of plane
    2: number of the state other than the background state
    3: background state or empty state
    4: state I want to put
    5: fraction of cells that get S1 state (0 to 1)
  */
  
  //InitialSet(world,nspecies,0,1,0.1,2,0.1,3,0.1,4,0.1,5,0.1);
  //printf("Hello1\n");
  int i,j;
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){
    //homogeneous init.
    //world[i][j].val=(int)((1+nspecies)*genrand_real1());
    
    //split random init. for two hyperc
    //if(j<ncol/2){
      world[i][j].val = (int)((1+hyperc1)*genrand_real1());
      world[i][j].fval = death;
    //}
    /*
    else{
      world[i][j].val = hyperc1 +1 + (int)((hyperc2)*genrand_real1());
      world[i][j].fval = death;
    }*/
  }
  
  //InitialSet(food,1,0,1,0.3);
  //printf("Hello2\n");
  //ReadSavedData("glidergun.sav",1,world);
  //printf("\n\nworld is ready. Click or scroll your mouse to update the CA.\n\n");
  Boundaries2(world);
}

void NextState(int row,int col)
{
  int rn1,ival,neirow,neicol;
  int sum;
  double prepl;
  int cat_val;
  
  if( world[row][col].val == 0 ){
    rn1 =(int)(1. + 8.*genrand_real2() );   //random direction: [1,8]
    ival = GetNeighbor( world,row,col, rn1 ); //val of molecule that can replicate
    //printf("Got zero, at row,col %d %d, random dir = %d, ival = %d\n",row,col,rn1,ival);
    if( 0 != ival ){
      //printf("ival %d\n",ival);
      GetNeighborC(world, row,col,rn1,&neirow,&neicol); // where is it?
      
      // Uncomment this for single hypercycle
      cat_val = 1+ (ival)%hyperc1 ;
      sum = CountMoore8(world, cat_val , neirow,neicol); //look at its neigh, and see how many mol ival+1 are there 
      prepl = uncat_birth + cat_birth1*sum;
      /*
      cat_val = hyperc1 + 1+ (ival)%hyperc2; //goes from hyperc1 to hyperc2 both incl
      sum = CountMoore8(world, cat_val , neirow,neicol); //look at its neigh, and see how many mol ival+1 are there 
      prepl += cat_birth2*sum;
      */
      //prepl = uncat_birth + cat_birth2*(double)sum;
      
      //this is for two hypercycles
      /*
      if(ival<=hyperc1) {
        cat_val = 1+ (ival)%hyperc1; //goes from 1 to hyperc1 both incl
        sum = CountMoore8(world, cat_val , neirow,neicol); //look at its neigh, and see how many mol ival+1 are there 
        prepl = uncat_birth + cat_birth1*(double)sum;
        if(cat_val>hyperc1) printf("Error\n");
      }
      else{
        cat_val = hyperc1 + 1+ (ival)%hyperc2; //goes from hyperc1 to hyperc2 both incl
        sum = CountMoore8(world, cat_val , neirow,neicol); //look at its neigh, and see how many mol ival+1 are there 
        prepl = uncat_birth + cat_birth2*(double)sum;
        if(cat_val<=hyperc1) printf("Error\n");
      }
      */
      //printf("Because %d != 0, \n",ival);
      
      if(genrand_real1() < prepl) {
        world[row][col].val=ival;
        //printf("ival %d\n",ival);
        //world[row][col].val2=ival;
      }
    }
  }else{
    //hypercycles compete for death rate
    /*
    if( world[row][col].val <= hyperc1 ){
      if( genrand_real2() < death1 )
        world[row][col].val = 0;
    } 
    else{
      if( genrand_real2() < death2 )
        world[row][col].val = 0;
    }
    */
    
    //for single death rate
    if(genrand_real1() < world[row][col].fval) 
      world[row][col].val = 0;
  }
  
}

void Update(void)
{ 
  char moviedir[MYBUFFER] = "hyperc_competition";
  
  //Display(world);
  
  Plot(1,world);
  
  if(Time % 50 ==0) DrawSlide(world, moviedir);
  
  Asynchronous();
  //Synchronous(1,world);
  
  if( Time % 10 ==0 ) MDiffusion(world);
  
  //while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement
  
}




