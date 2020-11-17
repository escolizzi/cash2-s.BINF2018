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

#define SELF 0
#define NEI 1

static TYPE2** world;
//const int Nl=2;
//const int Jtable[3][3] = {{0,8, 8 },{8, 12, 20},{8, 20 , 12}};

const int Nl=5;
//Notice that interaction with environment is twice what's written, so if you see 8, it actually means 16, why this is the case, ask Renske
const int Jtable[6][6] = {{0,  12,  8,  8,  12 , 12},
                          {12,  8, 10, 20, 20, 20},
                          {8,  10, 12, 20, 20, 20},
                          {8,  20, 20, 12, 12, 12},
                          {12, 20, 20, 12, 10, 10},
                          {12, 20, 20, 12, 10, 10},
                                                }; 

double kdiff = 1.;
double beta=.3;

void Initial(void)
{
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 100; /* # of row (default=100)*/
  ncol = 100; /* # of column (default=100)*/
  nplane = 1; /* # of planes (default=0)*/
  scale = 4; /* size of the window (default=2)*/
  boundary = FIXED; /* the type of boundary: FIXED, WRAP, ECHO (default=WRAP). Note that
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
  
  /*
  int i,j;
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){
    world[i][j].val=0;
    if(pow(i-nrow/2,2)+pow(j-ncol/2,2) < 200){
      world[i][j].val=1;
      if(genrand_real2()<.5) world[i][j].val=2;

      
    }
    
  }
  */
  InitialSet(world,Nl,0,1,0.065,2,0.2,3,0.13,4,0.13,5,0.13);
  
  //ReadSavedData("glidergun.sav",1,world);
  //printf("\n\nworld is ready. Click or scroll your mouse to update the CA.\n\n");
  //Boundaries2(world);
  
  //fprintf(stderr,"DEBUG ME, I am not working very well!!!\n");
}

int SumNei(int *bla)
{
  int i,sum=0;
  for(i=0;i<Nl;i++) sum += bla[i];
  
  return sum;
}


//returns nrg of pairwise interaction of lipids
int NrgOfPair(int val1, int val2)
{
  return (Jtable[val1][val2]);
}

int GetNrgConfig(TYPE2 **world, int row, int col, int val )
{
  int i,nrg=0;
  TYPE2 *nei;
  for(i=1;i<=8;i++) {
    nei = GetNeighborP(world,row,col,i);
    nrg+= NrgOfPair(nei->val, val);
  }
  return nrg;
}

int IsBeyondBoundary(int nr,int nc)
{
  if(nr>nrow || nc>ncol || nr<1 ||nc<1) 
    return 1;
  else
    return 0;
}

void NextState(int row,int col)
{
  int nr,nc,rn;
  //int tsum[2][Nl];
  //int i,j;
  //TYPE2 *nei;
  //printf("NextState(): Hello0, row, col: %d %d\n",row,col);
  rn=1+(int)(8.*genrand_real2());
  GetNeighborC(world,row,col,rn,&nr,&nc); //pick a random neighbour
  if(IsBeyondBoundary(nr,nc) && boundary == FIXED){
    //always give nrg of exchange with beyond boundary as high (adsorbing boundary condition)
    //so exchanged is always accepted:
    world[row][col].val =0;
    
    //int nrg_now,nrg_fut;
    //nrg_now = GetNrgConfig(world,row,col, world[row][col].val); // + 0. from beyond boundary
    //nrg_fut = GetNrgConfig(world,row,col, 0); //Fixed Boundaries: the future is just the present minus the molecule at row,col
    //if( (nrg_fut < nrg_now &&  genrand_real2()< kdiff/2.) ||  (genrand_real2()< kdiff/2. * exp(-beta*( nrg_fut - nrg_now )) ) ){
    //  int temp=world[row][col].val;
    //  world[row][col].val = world[nr][nc].val;
    //  world[nr][nc].val =temp;
    //}
    return;
  }else{
  //printf("NextState(): Hello1, neirow neicol: %d %d\n",nr,nc);
  //for(j=0;j<Nl;j++){
  //  tsum[SELF][j]=CountMoore8(world,j+1,row,col); //count how many neigh ar at row,col
  //  tsum[NEI][j]=CountMoore8(world,j+1,nr,nc);    //count how many neigh ar at nei row, nei col
  //}
  
  //if(world[row][col].val==world[nr][nc].val){
  //  ;
  //}
  
  
  //As long as there is at least one particle between self and nei, movements can happen
  //Movements are accepted if free nrg_future<nrg_now, or if prob, drawn from boltzmann distr.
  
  if( world[row][col].val + world[nr][nc].val >0 ){
    //printf("NextState(): Hello1.1\n");
    int nrg_now=0,nrg_fut=0;
    nrg_now = GetNrgConfig(world,row,col, world[row][col].val) + GetNrgConfig(world,nr,nc, world[nr][nc].val);
    //printf("NextState(): Hello1.2\n");
    nrg_fut = GetNrgConfig(world,row,col, world[nr][nc].val) + NrgOfPair(world[row][col].val,world[nr][nc].val) -  NrgOfPair(world[nr][nc].val,world[nr][nc].val) 
            + GetNrgConfig(world,nr,nc, world[row][col].val) - NrgOfPair(world[row][col].val,world[row][col].val) + NrgOfPair(world[row][col].val,world[nr][nc].val);
    //printf("NextState(): Hello1.3\n");
    if( (nrg_fut < nrg_now &&  genrand_real2()< kdiff/2.) ||  (genrand_real2()< kdiff/2. * exp(-beta*( nrg_fut - nrg_now )) ) ){
      int temp=world[row][col].val;
      world[row][col].val = world[nr][nc].val;
      world[nr][nc].val =temp;
    }
    
  }
  //printf("NextState(): Hello2\n");
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
  //printf("Update(): Hello0\n");
  Display(world);
  //printf("Update(): Hello1\n");
  Asynchronous();
  //Synchronous(1,world);
  // printf("Update(): Hello2\n");
  //while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement
  
}
