#include <stdio.h> 
// #include <stdlib.h>
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

double init_p_formcomplex = 0.9;
double p_decay = 0.02; //everything dies at 0.03 

void Decay(TYPE2 **world);
int CountVal(int val);
double Avrgfval(int val);

void Initial(void)
{
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 200; /* # of row (default=100)*/
  ncol = 200; /* # of column (default=100)*/
  nplane = 1; /* # of planes (default=0)*/
  scale = 4; /* size of the window (default=2)*/
  boundary = WRAP; /* the type of boundary: FIXED, WRAP, ECHO (default=WRAP). Note that
		      Margolus diffusion is not supported for ECHO. */
  ulseedG = time(NULL); /* random seed (default=56)*/
  fprintf(stderr,"Seeding with time: %ld\n", ulseedG);
  
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
  // InitialSet(world,3,0,1,0.2,2,.2,3,.2);
  //ReadSavedData("glidergun.sav",1,world);
  //printf("\n\nworld is ready. Click or scroll your mouse to update the CA.\n\n");
  int row,col;
  for(row=1;row<=nrow;row++)for(col=1;col<=ncol;col++){
    world[row][col].val=0;
    world[row][col].crow=-1;
    
    double drn = genrand_real2();
    if(drn < 0.2) {
      world[row][col].val = 1;
      world[row][col].fval = init_p_formcomplex; //initial
      world[row][col].fval2 = init_p_formcomplex;
    }
    else if(drn < 0.4) {
      world[row][col].val = 2;
      world[row][col].fval = init_p_formcomplex;
    }
    else if(drn < 0.6) {
      world[row][col].val = 3;
      world[row][col].fval = init_p_formcomplex;
    }
    
  }
  
  Boundaries2(world);
}

void NextState(int row,int col)
{
  int rn, grow,gcol, catrow,catcol;
  // double drn;
  TYPE2 *nei; 
  TYPE2* neinei;
  
  if(world[row][col].val==0){
    rn = 1. + 8.*genrand_real2(); //pick random nei
    nei = GetNeighborP(world,row,col,rn);
    
    int neirow,neicol;      // just for checking
    GetNeighborC(world,row,col,rn, &neirow,&neicol );
    // printf("Neirow and neicol = %d,%d\n",neirow,neicol );
    if(nei->crow != -1){
      neinei = &world[nei->crow][nei->ccol]; //nei of nei
      //find who is who
      // 1 is genome
      // 2 is replicase
      // 3 is endonuclease
      if(nei->val == 1 ){
        //nei is genome
        grow = neinei->crow;
        gcol = neinei->ccol;
        catrow = nei->crow;
        catcol = nei->ccol;
      }else if(neinei->val == 1){
        //neinei is genome
        grow = nei->crow;
        gcol = nei->ccol;
        catrow = neinei->crow;
        catcol = neinei->ccol;
      }else{
        fprintf(stderr,"NextState(): Error. got complex with no genome?\n");
        fprintf(stderr,"world[row][col].val = %d, nei->val = %d, neinei->val = %d\n", world[row][col].val, nei->val,neinei->val);
        fprintf(stderr,"world[row][col].crow,ccol = %d,%d; nei->crow,ccol = %d,%d, neinei->crow,ccol = %d,%d\n", world[row][col].crow, world[row][col].ccol, nei->crow,nei->ccol, neinei->crow, neinei->ccol);
        // while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement
        exit(0);
      }
      
      // perform reaction and set complex free
      if(world[catrow][catcol].val == 2){
        //got replicase
        world[row][col].val = world[grow][gcol].val; //replicase replicates genome
        //undo complex
        world[row][col].crow =-1;
        world[grow][gcol].crow =-1;
        world[catrow][catcol].crow =-1;
        
        //Mutations in the new genome
        world[row][col].fval = world[grow][gcol].fval + (genrand_real1() - 0.5) / 10.;
        if(world[row][col].fval>1) world[row][col].fval = 2. - world[row][col].fval;
        if(world[row][col].fval<0) world[row][col].fval = -world[row][col].fval;
        
        world[row][col].fval2 = world[grow][gcol].fval2 + (genrand_real1() - 0.5) / 10.;
        if(world[row][col].fval2>1) world[row][col].fval2 = 2. - world[row][col].fval2;
        if(world[row][col].fval2<0) world[row][col].fval2 = -world[row][col].fval2;
        
      }else if(world[catrow][catcol].val == 3){
        //got endonuclease
        // eventually randomise who gets into which position
        world[row][col].val = 2; // this is the new replicase
        // world[row][col].fval = fval it's the same fval
        
        world[grow][gcol].val = 3;
        world[grow][gcol].fval = world[row][col].fval2;
        
        //undo complex
        world[row][col].crow=-1;
        world[grow][gcol].crow=-1;
        world[catrow][catcol].crow=-1;
        
      }else{
        printf("NextState(): Error. got catalyst neither repl nor transcr?\n");
      }
    }
    //if nei not in complex -> do nothing
  }else{
    //if world[row][col] is not empty, if it is NOT a genome and nei is a genome complex can happen, and neither of them are already in complex
    if( (world[row][col].val == 2 || world[row][col].val == 3) && world[row][col].crow==-1 ){
      rn = 1. + 8.*genrand_real2(); //pick random nei
      nei = GetNeighborP(world,row,col,rn);
      //check if nei is genome, if yes make complex
      if(nei->val==1 && nei->crow==-1  && genrand_real2() < world[row][col].fval ){
        //find coordinates of nei
        GetNeighborC(world,row,col,rn,&grow,&gcol);
        
        world[row][col].crow = grow; // set me as complex
        world[row][col].ccol = gcol;
        
        world[grow][gcol].crow = row;
        world[grow][gcol].ccol = col;
        
      }
    }
  }
}

void Decay(TYPE2 **world)
{
  int row,col;
  for(row=1;row<=nrow;row++)for(col=1;col<=ncol;col++){
    if(world[row][col].val != 0){
      if( genrand_real2() < p_decay ){
        world[row][col].val = 0;
        // if in complex, undo complex
        if(world[row][col].crow != -1){
          //first undo complexee
          world[ world[row][col].crow ][ world[row][col].ccol ].crow=-1;
          //and only after undo this complex
          world[row][col].crow = -1;
        }
      }
    }
  }
}

int CountVal(int val){
  int row,col,counter=0;
  for(row=1;row<=nrow;row++)for(col=1;col<=ncol;col++){
    if(world[row][col].val == val) counter++;
  }
  return counter;
}
double Avrgfval(int val)
{
  int row,col,counter=0;
  double fval;
  for(row=1;row<=nrow;row++)for(col=1;col<=ncol;col++){
    if(world[row][col].val == val){
      fval += world[row][col].fval;
      counter++;
    }
  }
  if(counter) fval /= (double)counter;
  else fval = 0.; 
  return fval;
}

void Update(void)
{
  Display(world);
  Asynchronous();
  //Synchronous(1,world);
  
  // decay can also happen
  Decay(world);
  //Diffusion!
  
  if(Time%100 ==0) {
    int howmany_genomes = CountVal(1);
    int howmany_repl = CountVal(2);
    int howmany_end = CountVal(3);
    printf("# genomes = %d, # repl = %d, # endonucl = %d, \n", howmany_genomes, howmany_repl, howmany_end);
    printf("Avrg repl: %f, Avrg end: %f\n", Avrgfval(2), Avrgfval(3));
  }
  
  //while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement
  
}
