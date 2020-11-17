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

double NONE = 0.5;
double pdeath = 0.2;
double init_pmut = 0.05;
int init_gs=50;
double fdup = 0.49;
double fdel = 0.49;
double flof = 0.01;

double p_BGlof=0.001;

double tar=50;

static TYPE2** Life;
static TYPE2** Genm;
static TYPE2** Mutr;

void Mutations(TYPE2* ic);
//double Fitness(TYPE2* ic);

void Fitness(TYPE2* ic){
  int a=0,i=0;
  double fitness;
  int gs = ic->val2;
  //printf("gs = %d\n",gs);
  
  for(int k=0;k<gs;k++){
    //printf("%c ",ic->seq[k]);
    if(ic->seq[k]=='1') a++;
    else if(ic->seq[k]=='0') i++;
  }
  //printf("\n");
  
  //printf("genome = %s\n",ic->seq);
  //printf("a = %d, i =% d\n",a,i);
  
  if(a>=tar) fitness = 1.;
  else fitness = a/tar;
  
  fitness *= exp(-0.01*i);
  fitness *= exp(-pow(0.001*gs,4.));
  
  ic->fval = fitness;
}

void Initial(void)
{
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 200; /* # of row (default=100)*/
  ncol = 200; /* # of column (default=100)*/
  nplane = 3; /* # of planes (default=0)*/
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
  MakePlane(&Life,&Genm,&Mutr);

  /* InitialSet(1,2,3,4,5)
    1: name of plane
    2: number of the state other than the background state
    3: background state or empty state
    4: state I want to put
    5: fraction of cells that get S1 state (0 to 1)
  */
  //InitialSet(Life,1,0,1,0.3);
  
  TYPE2* ic;
  
  for(int i=1;i<=nrow;i++)for(int j=1;j<=ncol;j++){
    ic = &Life[i][j];
    ic->val=1;
    ic->val2=init_gs; // initial genome size
    for(int k=0; k<init_gs;k++){
      ic->seq[k]='1';
    }
    for(int k=init_gs; k<256;k++){
      ic->seq[k]='\0';
    }
    
    ic->fval2=init_pmut;
    Fitness(ic);
    //printf("fitness = %f\n",ic->fval);
  }
  
  InitialSet(Genm,0,0);
  InitialSet(Mutr,0,0);
  Boundaries2(Life);
  Boundaries2(Genm);
  Boundaries2(Mutr);
  
  // Initialise colorRGB table (for colouring 2nd and 3rd plane)
  int r=0,g=0,b=255;
  double nr=102.;    //nr of ColorRGB's you want to create
  double range=1275.;  // range of coloursteps: 255*5= 1275
  double x = range/nr;  	 
  int y=0,c;

  for(c=0;c<(int)range;c++){
          if(c<255){			//starts blue
                  r = r + 1;		//goes magenta
          }else if(c<255*2){	
                  b = b - 1;		//goes red
          }else if(c<255*3){
                  g = g + 1;		//goes yellow
          }else if(c<255*4){
                  r = r -1;    		//goes green
          }else if(c<255*5){
                  b = b + 1;		//goes cyan
          }

          if(c == (int)(y*x+0.5)){
                  ColorRGB(10+y,r,g,b);	//so colour indexes in the gradient start at 10
                  y++;
          }
  }
  
}

void NextState(int row,int col)
{
  
  if(Life[row][col].val == 0){
    double sumfitness = NONE;
    double counterfitness=0.;
    int i;
    for(i=1;i<=8;i++) sumfitness += (GetNeighborP(Life,row,col, i)->val)? GetNeighborP(Life,row,col, i)->fval: 0.;
    //printf("sumfitness = %f\n",sumfitness);
    //exit(1);
    double rn = sumfitness * genrand_real2();
    for(i=1;i<=8;i++){
      counterfitness += GetNeighborP(Life,row,col, i)->fval;
      if(counterfitness> rn){
        Life[row][col] = *GetNeighborP(Life,row,col, i);
        // mutations
        Mutations(&Life[row][col]); // also sets fitness
        
        break;
      }
    }
    
  }else{
    if(genrand_real2() < pdeath){ 
      Life[row][col].val = 0;
      Life[row][col].fval = 0.;
    }
  }
  
  if(Life[row][col].val == 0) Mutr[row][col].val = 0;
  else Mutr[row][col].val = 10+(int)(100.*Life[row][col].fval2);
  
}

void Update(void)
{ 
  double avpmut=0,varpmut=0,ava=0,avi=0,avfitness=0;
  int popsize=0;
  if(Time%10){
    for(int i=1;i<=nrow;i++)for(int j=1;j<=ncol;j++){
      if(Life[i][j].val){
        popsize++;
        avpmut+=Life[i][j].fval2;
        varpmut+=Life[i][j].fval2*Life[i][j].fval2;
        avfitness+=Life[i][j].fval;
      }
    }
    
    printf("Av fitness = %f, Av mutrate = %f, std mutrate = %f\n", avfitness/(double)popsize, avpmut/(double)popsize, sqrt(  varpmut/(double)popsize - pow(avpmut/(double)popsize , 2.) )  );
  }
  Display(Life,Genm,Mutr);
  Asynchronous();
  //Synchronous(1,Life);

  //while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement
  
  
  
}

void Mutations(TYPE2* ic)
{
  
  //ic->val2 : genome size
  //ic->val2 : mutation rates
  int howmany_BG_mut = bnldev(ic->fval2 + p_BGlof,ic->val2);
  //printf("Nr. mutations: %d, genome size: %d\n",howmany_BG_mut,ic->val2);
  // first finds how many mutations happen
  //int howmany_BG_mut = BINOMIAL(gs, totpmut);
  //howmany_mut+=howmany_BG_mut;
  if(howmany_BG_mut>0){
    //printf("Genome bef: %s, mut: ",ic->seq);
    int gs = ic->val2;
    for(int i=0;i<howmany_BG_mut;i++){
      //choose site
      int pos=(int)(gs*genrand_real2());
      //choose type of mutation
      double rn = (ic->fval2+p_BGlof)*RANDOM();
      if(rn<fdup * ic->fval2){
        //printf("D");
        //duplication happens - MAX genes 256
        if(gs>=255) continue;
        for(int j=gs; j>pos; j-- ) ic->seq[j]=ic->seq[j-1];
        gs++;
      }else if(rn < (fdup+fdel)*ic->fval2 ){
        //printf("d");
        //delete - MIN genes 0
        if(gs==0) continue;
        for(int j=pos; j<gs; j++ ) ic->seq[j]=ic->seq[j+1];
        gs--;
      }
      else{
        //printf("i");
        ic->seq[pos]='0';
      }
    }
    
    ic->val2 = gs;
    //printf("\nGenome aft: %s\n",ic->seq);
  }
  
  if(genrand_real2()<0.02){
    ic->fval2 += 2.*( genrand_real1() - 0.5 )/100.;
    if(ic->fval2 > 1.) ic->fval2 =1.;
    if(ic->fval2 <=0.) ic->fval2 *= -1.;
  }
  if(howmany_BG_mut) Fitness(ic); // set fitness
  
}




