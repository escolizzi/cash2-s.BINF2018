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

void AllocateJtable(int ***dJtable, int Nl);
void InitialiseJtable(int **dJtable, int Nl);
void InitialiseWorld(TYPE2 **world,int Nl);

static TYPE2** world;
//const int Nl=2;
//const int Jtable[3][3] = {{0,8, 8 },{8, 12, 20},{8, 20 , 12}};
//static int **Jtable;    //allocate this dynamically in the future

const int Nl=30;

static int **dJtable;

//Notice that interaction with environment is twice what's written, so if you see 8, it actually means 16, why this is the case, ask Renske
const int Jtable[8][8] = {{0,  12,  8,  8,  12 , 12, 12, 8},
                          {12,  8, 10, 20, 20, 20, 20, 20},
                          {8,  10, 12, 20, 20, 20, 20, 20},
                          {8,  20, 20, 12, 12, 12, 12, 20},
                          {12, 20, 20, 12, 10, 10, 10, 10},
                          {12, 20, 20, 12, 10, 10, 10, 10},
                          {12, 20, 20, 12, 10, 10, 10, 10},
                          {8,  20, 20, 20, 10, 10, 10, 12},
                                                         }; 


/*const int Nl=3;
//Notice that interaction with environment is twice what's written, so if you see 8, it actually means 16, why this is the case, ask Renske
const int Jtable[4][4] = {{0,   8,  12,  12},
                          {8,   16, 12, 12},
                          {12,  12, 10, 10},
                          {12,  12, 10, 10},
                                          }; 

*/
double kdiff = 1.;
double beta=.2;
double kloss = .0001;//00001; //the higher this term, the smaller the steady state value of the forming vescicles
                      //intuitively, if kloss=1 the system is constantly being mixed, 
                      // and if it is 0 it slowly converges to steady state
                      // WITH THIS and the above matrix I have seen vescicles splitting!!!

void Initial(void)
{
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 150; /* # of row (default=100)*/
  ncol = 150; /* # of column (default=100)*/
  nplane = 1; /* # of planes (default=0)*/
  scale = 4; /* size of the window (default=2)*/
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
  
  AllocateJtable(&dJtable, Nl);
  InitialiseJtable(dJtable, Nl);
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
  InitialiseWorld(world,Nl);
  //InitialSet(world,Nl,0,1,0.125,2,0.0525,3,0.0525,4,0.125,5,0.125,6,0.125,7,0.0525);
  //InitialSet(world,Nl,0,1,0.5,2,0.15,3,0.15);
  //ReadSavedData("glidergun.sav",1,world);
  //printf("\n\nworld is ready. Click or scroll your mouse to update the CA.\n\n");
  Boundaries2(world);
  
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
  return (dJtable[val1][val2]);
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

void NextState(int row,int col)
{
  int nr,nc,rn;
  //int tsum[2][Nl];
  //int i,j;
  //TYPE2 *nei;
  
  rn=1+(int)(8.*genrand_real2());
  GetNeighborC(world,row,col,rn,&nr,&nc); //pick a random neighbour
  
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
    int nrg_now=0,nrg_fut=0;
    nrg_now = GetNrgConfig(world,row,col, world[row][col].val) + GetNrgConfig(world,nr,nc, world[nr][nc].val);
    nrg_fut = GetNrgConfig(world,row,col, world[nr][nc].val) + NrgOfPair(world[row][col].val,world[nr][nc].val) -  NrgOfPair(world[nr][nc].val,world[nr][nc].val) 
            + GetNrgConfig(world,nr,nc, world[row][col].val) - NrgOfPair(world[row][col].val,world[row][col].val) + NrgOfPair(world[row][col].val,world[nr][nc].val);
    if( (nrg_fut < nrg_now &&  genrand_real2()< kdiff/2.) ||  (genrand_real2()< kdiff/2. * exp(-beta*( nrg_fut - nrg_now )) ) ){
      int temp=world[row][col].val;
      world[row][col].val = world[nr][nc].val;
      world[nr][nc].val =temp;
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

// FOR FUTURE: MAKE large scale swapping accepted based on energy
// for now random loss and add (into empty space) with mass conservation
// Actually - this is really mixing, there is nothing in the vapor after we are out of loop
void LossAndAdd(TYPE2 **world)
{
  int i,j;
  int l_vapor[nrow*ncol];
  int l_empty[nrow*ncol];
  int tot_vapor=0, tot_empty=0;
  
  //initialise arrays
  memset(l_vapor, 0, sizeof(int) * nrow * ncol);
  memset(l_empty, 0, sizeof(int) * nrow * ncol);
  //int flag=0;
  
  //we do lipids evaporation (into l_vapor), and we can also keep track of who is there
  for(i=1;i<=nrow;i++){
    for(j=1;j<=ncol;j++){
     if(world[i][j].val != 0){
      if(genrand_real2() < kloss){
        l_vapor[ tot_vapor ] = world[i][j].val;
        tot_vapor++;
        world[i][j].val = 0;
        //if this has become zero -> add this spot to the possible ones for lipid condensation onto surface
        l_empty[tot_empty]= nrow*i + j;
        tot_empty++;
        
        // FOR THE SAKE OF TESTING
        //fprintf(stderr,"LossAndAdd(): Warning placed a break, remove it. row,col=%d %d, code=%d\n", i,j, nrow*i + j);
        //flag=1;
        //break;
        
      }
    }else{
      //else this spot is zero, so it should be added for condensation anyway
      //fprintf(stderr,"LossAndAdd(): row,col=%d %d, code=%d... is this empty? val = %d\n", i,j, nrow*i + j,world[i][j].val);
      l_empty[tot_empty]= nrow*i + j;
      tot_empty++;
    }
      
  }
  //if(flag==1){fprintf(stderr,"Another break\n"); break;}
  }
  
  //fprintf(stderr,"Tot empty space= %d, tot lipid evaporated = %d\n", tot_empty, tot_vapor);
  
  //condensation
  //this can be randomised on the fly:
  //also notice that tot_vapor <= tot_empty
  for(i=tot_vapor-1;i>=0;i--){
    //take random lipid
    int ran_lip = (int)( (i+1)*genrand_real2() );
    int lipid = l_vapor[ran_lip];
    l_vapor[ran_lip] = l_vapor[i];
    
    //take a random position from l_empty list
    int ran_pos = (int)( tot_empty*genrand_real2() );
    int bla = l_empty[ran_pos];
    int row= l_empty[ran_pos] / nrow;
    int col= l_empty[ran_pos] % nrow;
    if(col==0){
      col=ncol;
      row--;
    }
    
    //swap that random pos with the last tot_empty list
    l_empty[ran_pos] = l_empty[tot_empty-1];
    tot_empty--;
    if(world[row][col].val != 0 ) 
      fprintf(stderr,"LossAndAdd(): Error. placing lipid on non empty spot, code=%d, decoded as row col = %d %d\n",bla,row,col);
    world[row][col].val= lipid;
    
  }

  //fprintf(stderr,"LossAndAdd(): this is it, i=%d\n",i);
  //exit(1);
}

int CountPop(TYPE2 **world)
{
  int i,j, tot=0;
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){
    if(world[i][j].val != 0) tot++;
  }
  return tot;
}

void Update(void)
{ 
  Display(world);
  Asynchronous();
  //Synchronous(1,world);
  LossAndAdd(world);
  //while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement
//   if(!(Time%100)){
//     printf("Tot pop: %d\n", CountPop(world));
//   }
}

void InitialiseWorld(TYPE2 **world,int Nl)
{
  int i,j;
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){
    if(genrand_real2() < 0.4) world[i][j].val=0;
    else world[i][j].val = 1 + (int)( Nl * genrand_real2() ) ;
  }
}

void InitialiseJtable(int **dJtable, int Nl)
{
  int i,j,Jval;
  dJtable[0][0]=0;
  //initialise the top half with some values
  for(i=0;i<=Nl;i++)for(j=i;j<=Nl;j++){
    if(i==0 && j==0) continue;  //already initialised to zero
    Jval=5+(int)(20*genrand_real2());
    if(i==0) Jval/=2; //half the number
    dJtable[i][j]=Jval;
  }
  
  //copy those values to the other half
  for(i=1;i<=Nl;i++)for(j=0;j<i;j++){
    dJtable[i][j] = dJtable[j][i];
  }
  
  //print table
  printf("Jtable thus created\n");
  for(i=0;i<=Nl;i++){
    for(j=0;j<=Nl;j++){
      printf("%d\t", dJtable[i][j]);
    }
    printf("\n\n");
  }
}

//dynamically allocate dJtable
void AllocateJtable(int ***dJtable, int Nl)
{
  int i;
  
  *dJtable = (int **)malloc( (Nl+1) * sizeof(int *) );
  for(i=0;i<Nl+1;i++){
    (*dJtable)[i]=(int *)malloc( (Nl+1) * sizeof(int) );
  }
  (*dJtable)[0][0]=1;
  //printf("%d\n",(*dJtable)[0][0]);
  //exit(1);
}





