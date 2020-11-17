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
int MAXSIZE=256; //keep this same size as seq is cash.h
double beta=0.05; //scaling factor for probability of replication
int MAX_ANTIB = 16; // max number of antibiotics

static int init=0;
static double*** antib;
void UpdateAntibPlane(double ****pantib, TYPE2 **world); //antib generated, degraded, diffusion
void MakeAntibioticPlane(double ****a,int nrow,int ncol,int MAX_ANTIB);

//antib generated, degraded, diffusion
void UpdateAntibPlane(double ****pantib, TYPE2 **world)
{
  //antibiotics are generated


}

int Genome2genenumber(char *seq, char gene)
{
  int i,genecount=0;
  int seqlen=strlen(seq);
  for(i=0;i<seqlen;i++) if(seq[i]==gene) genecount++;
  return genecount;
}

// void Replicate( int row, int col,int dir)
// {
//   world[row][col] = *GetNeighborP(world,row,col,dir);
// }

int Mutate(int row, int col)
{
  return 0;
}

void Initial(void)
{
  fprintf(stderr,"Use the file in projects/projectStrepto, not this one");
  exit(1);
  
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 200; /* # of row (default=100)*/
  ncol = 200; /* # of column (default=100)*/
  nplane = 1; /* # of planes (default=0)*/
  scale = 4; /* size of the window (default=2)*/
  boundary = WRAP; /* the type of boundary: FIXED, WRAP, ECHO (default=WRAP). Note that
		      Margolus diffusion is not supported for ECHO. */
          
  ulseedG = time(NULL); /* random seed ... if you don't know the time */
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
  int i,j,k;
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){
    world[i][j].val=0;
    for(k=0;k<MAXSIZE;k++) world[i][j].seq[k]='\0';
    if(genrand_real2()<0.001){
      world[i][j].val=1;
      world[i][j].val=10.*genrand_real2(); // just for coloring
      world[i][j].seq[0]='g'; //gene for growth
      world[i][j].val2 = Genome2genenumber(world[i][j].seq,'g');
    }
  }

  //InitialSet(world,1,0,1,0.001);
  //ReadSavedData("glidergun.sav",1,world);
  printf("\n\nworld is ready. Click or scroll your mouse to update the CA.\n\n");
  Boundaries2(world);
}

void NextState(int row,int col)
{
  //int sum1,sum0;
  int k;
  TYPE2 *nei;
  if(world[row][col].val==0){
    int howmany_emptynei = CountMoore8(world,0,row,col); //counts how many neigh are empty
    if(howmany_emptynei<7 || howmany_emptynei==8) return; // if number is less than this -> not enough resources, returns

    int howmanynei = 8 - howmany_emptynei; //so this is number of competing cells (1 or 2)
    // printf("howmanynei: %d\n", howmanynei);
    int dirarray[howmanynei];
    int grarray[howmanynei];

    int counter=0;
    int totgrowthgenes=0;

    //we find who these neighbors are:
    for(k=1;k<9;k++){
      nei = GetNeighborP(world,row,col,k);
      if(nei->val!=0){
        //save direction
        dirarray[counter]=k;
        //save number of growth genes
        grarray[counter]=nei->val2;
        totgrowthgenes+=nei->val2;

        counter++;
      }
      if(counter==howmanynei) break;
    }

    double totalgrowth=1. - exp(-beta*totgrowthgenes);  //total prob of replication
    double rn = genrand_real2();
    counter=0;
    double cumprob=0.0;
    if(rn < totalgrowth){
      while (rn>cumprob && counter<howmanynei){
        cumprob+=totalgrowth*grarray[counter]/(double)totgrowthgenes;
        counter++;
      }
      //now we replicate the individual at direction dirarray[counter] into the focal pixel
      //Replicate( row,col,dirarray[counter]);
      int direction=dirarray[counter-1];
      // printf("counter: %d, direction: %d\n",counter, direction);
      world[row][col] = *GetNeighborP(world,row,col,direction);
      Mutate(row,col);

    }
    // else no one replicates.


    //int neival=GetNeighbor(world,row,col, 1+(int)(8*genrand_real2()) );
    //if(neival>0 && CountMoore8(world,0,row,col)>6 && genrand_real2() < 0.5 ) world[row][col].val=neival;
  }

  // sum = CountMoore8(world,1,row,col);
  //
  //
  // if(sum==3 || (sum==2 && world[row][col].val==1))
  //   world[row][col].val = 1;
  // else
  //   world[row][col].val = 0;

}

void Update(void)
{
  //done only the first time step
  if(!init){
    MakeAntibioticPlane(&antib,nrow,ncol,MAX_ANTIB); //also initilizes it to zero
    init=1;
  }

  UpdateAntibPlane(&antib,world); //antib generated, degraded, diffusion

  Display(world);
  Asynchronous();
  //Synchronous(1,world);

  //while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement

}

void MakeAntibioticPlane(double ****a,int nrow,int ncol,int MAX_ANTIB)
{
  //double ***b;
  int i,j,k;
  (*a)=(double ***)malloc((nrow+2)*sizeof(double **));
  if((*a)==NULL){printf("Error at allocating antibiotics plane\n");exit(1);}

  for(i=0; i<=nrow+1;i++){
    (*a)[i]=(double **)malloc((ncol+2)*sizeof(double *));
    if((*a)[i]==NULL){printf("Error at allocating antibiotics plane\n");exit(1);}
  }

  for(i=0; i<=nrow+1;i++)for(j=0; j<=ncol+1;j++){
    (*a)[i][j]=(double *)malloc(MAX_ANTIB*sizeof(double));
    if((*a)[i][j]==NULL){printf("Error at allocating antibiotics plane\n");exit(1);}
  }
  // printf("Hello\n" );
  //(*a)[nrow][ncol][MAX_ANTIB -1] = 1.;
  //printf("%f\n", (*a)[nrow][ncol][MAX_ANTIB-1]);

  //a=&b;
  for(i=0; i<=nrow+1;i++)for(j=0; j<=ncol+1;j++)for(k=0;k<MAX_ANTIB;k++) (*a)[i][j][k]=0.;

  return;
}
