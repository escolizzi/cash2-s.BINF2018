/*PUBLIC GOOD MODEL
 Version 5.0:
 include command line backup file name
*/

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

#define BUFSIZE 128

static TYPE2** bact;
static TYPE2** good;
static TYPE2 empty={0,0,0,0,0,0.,0.,0.,0.,0.};

int MyArguments(int argc_g, char *argv_g[]);
void SetAllBack(int howmanycol,TYPE2 **bact);
int howmanycol=100;
double UpdateBactFitness(TYPE2 *nei,int neirow, int neicol);
double UpdateBactFitness_Nowak(TYPE2 **bact,TYPE2 *nei, int neirow, int neicol);
void PrintTheWelcomeStatement(void);

void ToBackup(TYPE2 **bact);
int GetColorIndexFrom(int val,double fval);
double ToMovie(TYPE2 **bact);
void ToSampleData(TYPE2 **bact);

void InitialiseFromBackup(TYPE2 **bact, TYPE2 empty);

double w=0.1; //intensity of selection for altruism fitness of individuals in Nowak's model is 
			  // f_i = 1 - w + w(  b*sum_j(fval_j) - sum_j(fval_i) )
double k_max=10.;// max possible production rate
double init_fval=5.;	//initial production rate of public good
double cst=4.2;	// production of good from resources: conversion rate = cost
double bnf=10.;	// intake of good from environment: conversion rate = benefit
double NONE=20.; //propbability that nothing happens ~ sort of

static double maxprob;  //bnf*k_max*7/8 ~ to let some probability of nothing happening, we just do b*k_max
double death=0.2;	//death rate
double kmut=0.0;	//mutation rate
double kdelta=.1;	//parameter to determine the mutation step: kdelta*(genrand_real1()-0.5 )
double kdiff=0.02;	//diffusion rate

int time_save_sample=5000;
int time_save_movie=100;
int time_backup=1000;
int sample_size=10000;
char par_movie_directory_name[BUFSIZE]="moviePB";
char par_data_file_name[BUFSIZE]="dataPB.txt";
char par_backup_directory_name[BUFSIZE]="backupPB";
char par_input_file_name1[BUFSIZE]="";

void Initial(void)
{
  int status;
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 729; /* # of row (default=100)*/
  ncol = 729; /* # of column (default=100)*/
  nplane = 2; /* # of planes (default=0)*/
  scale = 1; /* size of the window (default=2)*/
  boundary = WRAP; /* the type of boundary: FIXED, WRAP, ECHO (default=WRAP). Note that
		      Margolus diffusion is not supported for ECHO. */
  ulseedG = 2407824; /* random seed (default=56)*/
  
  /* useally, one does not have to change the followings */
  /* the value of boundary (default=(TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.})*/
  boundaryvalue2 = (TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.}; 
  display=0; //1 for display, 0 no display
  
  /*
    nrow = 80; 
    ncol = 400;
    boundary = FIXED;
  */
  
  status=MyArguments(argc_g, argv_g);
  if(status==1){ 
    PrintTheWelcomeStatement();
    exit(1);
  }
  
  OpenPNG(par_movie_directory_name, nrow, ncol);
  
}


void InitialPlane(void)
{
  MakePlane(&bact,&good);
  
  /* InitialSet(1,2,3,4,5)
    1: name of plane
    2: number of the state other than the background state
    3: background state or empty state
    4: state I want to put
    5: fraction of cells that get S1 state (0 to 1)
  */
  //InitialSet(Free,1,0,1,0.1);
  //InitialSet(Fixed,0,0);
  
  InitialSetS(bact,0,boundaryvalue2);
  //InitialSetS(good,0,boundaryvalue2);
  
  //InitGradationColor(1,0,1);	//
  
  //Fixed[100][100].val = 1;
  int i,j;
  
  /*
  for(i=nrow/2; i<nrow/2+40; i++) for(j=ncol/2; j<ncol/2 +40; j++){
    bact[i][j].val=1;
    bact[i][j].fval=init_fval;
    bact[i][j].fval2=0.;
  }
  */
  /*for(i=nrow/2; i<nrow/2+20; i++) for(j=ncol/2; j<ncol/2 +20; j++){
    bact[i][j].val=1;
    bact[i][j].fval=0.1;
    bact[i][j].fval2=0.;
  }
  */
  /*
  for(i=1; i<nrow; i++) for(j=21; j<40; j++){
    bact[i][j].val=1;
    bact[i][j].fval=init_fval;
    bact[i][j].fval2=0.;
  }
  */
  
  if(strlen(par_input_file_name1)!=0) InitialiseFromBackup(bact,empty);
  else{
    for(i=1; i<=nrow; i++) for(j=1; j<=ncol; j++){
      bact[i][j]=empty;
      if(genrand_real1()<0.15){// && i>nrow/2-25 && i<nrow/2+25 &&  j>ncol/2-25 && j<ncol/2+25){
        bact[i][j].val=1;
        bact[i][j].fval=init_fval;
        if(j<ncol/2) bact[i][j].fval=0.;
       }
    }
  }
  
  Boundaries2(bact);
  Boundaries2(good);
  
  int r=0,g=0,b=255;
  double nr=102.;    //nr of ColorRGB's you want to create
  double range=1530.;  // range of coloursteps: 255*6= 1530
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
    }else if(c<255*6){
      g = g - 1;		//back to blue
    }
    
    if(c == (int)(y*x+0.5)){
      ColorRGB(10+y,r,g,b);	/*so colour indexes in the gradient start at 10*/
      y++;
    }
  }
  
  maxprob = bnf*k_max;
  //maxprob = bnf*k_max*7./8.;
  //printf("\n\n\tDID YOU SET maxprob to b*max(k)*7/8? \n\n\tCHECK IF %f == %f \n\n", (bnf*k_max) , maxprob);
   
}

void DiffBySwap(TYPE2 **bact,int row, int col)
{
  int rnei;
  TYPE2 *nei; 
  TYPE2 tmp;
   
  rnei = 1 + (int)( 8*genrand_real2() );
  nei=GetNeighborP(bact,row,col,rnei);
  
  tmp=bact[row][col];
  bact[row][col]=*nei;
  *nei=tmp;
  
}

/*Same as UpdateBactFitness(), except you do not get any benefit
  from your own production, and still you pay the costs.
  Notice that this method is VERY similar to Nowak's.
*/
double UpdateBactFitness_no_own_gain(TYPE2 *nei, int neirow, int neicol)
{
  
  int k;
  double krepl;
  TYPE2 *neinei;
  
  krepl = -cst*nei->fval;
  //loop through all nei  ** NOT ** including self (k=0) and add everybody's contribution to you
  for(k=1;k<=8;k++){
    neinei=GetNeighborP(bact,neirow,neicol,k);
    if( neinei->val==1 ) krepl += (1./9.)*bnf*neinei->fval;
  }
  //printf("krepl=%f\n", krepl);
  krepl=( krepl>0. )?krepl:0;
  return krepl;
  
}

// fval is public good produced per time step, and it's evolvable, fval2 is fitness
double UpdateBactFitness(TYPE2 *nei, int neirow, int neicol)
{
  int k;
  double krepl;
  TYPE2 *neinei;
  
  krepl = -cst*nei->fval;
  //printf("krepl=%f\n", krepl);
  //loop through all nei including self and add everybody's contribution to you
  for(k=0;k<=8;k++){
    neinei=GetNeighborP(bact,neirow,neicol,k);
    if( neinei->val==1 ) krepl += (1./9.)*bnf*neinei->fval;
  }
  //printf("krepl=%f\n", krepl);
  krepl=( krepl>0. )?krepl:0;
  return krepl;
  
}

double UpdateBactFitness_Nowak(TYPE2** bact,TYPE2 *nei, int neirow, int neicol)
{
  int k;
  double krepl=0.;
  TYPE2 *neinei;
  
  //int neineirow,neineicol;
  //krepl = -cst*nei->fval;
  //printf("krepl=%f\n", krepl);
  
  //loop through all nei including self and add everybody's contribution to you
  for(k=1;k<=8;k++){
    neinei=GetNeighborP(bact,neirow,neicol,k);
    if( neinei->val==1 ){
      //GetNeighborC(bact,neirow,neicol,k,&neineirow,&neineicol);
      //printf("\n\t\tneineirow %d neineicol %d; ", neineirow, neineicol);
      
      krepl += (bnf*neinei->fval - cst*nei->fval);
      
    }
  }
  
  //printf("krepl=%f\n", krepl);
  krepl=( krepl>0. )?krepl:0;
  krepl = 1. - w + w*krepl;
  return krepl;
}

/********   ---   AGENT BASED   ---   ************/
// Agent based replication works like this:
// an individual gets chosen for replication with probability:
// fval2_i/sum_nei(fval2) * ( 1 - exp( -b*sum_nei(fval2) ) )
// essentially like SIR model

int AgentBasedReplication(TYPE2 **bact, int row,int col, TYPE2 *child)
{
  int k,neirow,neicol,l=0;
  double rsum,rn;
  double totsum=NONE;
  TYPE2 *neis[8];
  int success=0;
  
  for(k=1;k<=8;k++){
    //Notice that if this is empty, next time we are overwriting this address.
    //If it is all empty, then l=0, which is excluded later
    neis[l]=GetNeighborP(bact,row,col,k);
    if( neis[l]->val==1 ){
      GetNeighborC(bact,row,col,k,&neirow,&neicol);
      //printf("\n\tneirow %d neicol %d; ", neirow, neicol);
      
      //goes through the neighbourhood of every individual and assign them fitness
      //neis[l]->fval2=UpdateBactFitness_Nowak(bact,neis[l],neirow,neicol);
      //neis[l]->fval2=UpdateBactFitness(neis[l],neirow,neicol);  
      neis[l]->fval2=UpdateBactFitness(neis[l],neirow,neicol);  
      totsum+=neis[l]->fval2;
      l++;
    }
  }
  
  //l now says how many neighbours are competing
  if(l!=0 && totsum>0.){
    //replication happens: who's the lucky one?
    //rn= totsum*genrand_real2();
    //rn= 8.*maxprob*genrand_real2(); 					//This is problematic!
    //rn= NONE+8.*maxprob*genrand_real2(); 					//THIS too
    rn = genrand_real2() * totsum/( 1. - exp(-totsum) );			//Totally CLAIM !!!  -- totsum already includes NONE
    rsum=0.;
    for(k=0;k<l;k++){
      rsum+=neis[k]->fval2;
      //this is the lucky one
      if(rn<rsum){
        *child=*neis[k];
        success=1;
        break;
      }
    }
    /*
    if(rsum>=l*maxprob) {
      printf("rsum too large =%f\n",rsum);
      printf("%d competitors, fvals2: ",l);
      for(k=0;k<l;k++) printf("%f ", neis[k]->fval2);
      printf("\n");
    }
    */
  }
  //printf("Success: %d, Claim: %d %f\n", success,child->val,child->fval);
  return success;
}

		
		
/*
		RANDOM NEIGHBOUR   --- UNTESTED !!!
*/		
int RandomNeiReplication(TYPE2 **bact, int row,int col, TYPE2 *child)
{
  int k,neirow,neicol,success=0;
  double rn;
  TYPE2 *nei;
  
  k=1 + (int)( 8.*genrand_real2());		// get random direction from focal point
  GetNeighborC(bact,row,col,k, &neirow,&neicol);	// get coordinates
  nei=GetNeighborP(bact,row,col,k);					// and pointer to neighbour
  if(nei->val==1){
    nei->fval2=UpdateBactFitness(nei,row,col); 		// calculate fintess
    rn=maxprob*genrand_real2(); 
    if(rn<nei->fval2){
      *child=*nei;									//replication happens
      success=1;
    }
  }
  return success;
}


/***********************/
/* ------ CLAIM ------ */
/***********************/
int ClaimReplication(TYPE2 **bact, int row,int col, TYPE2 *child)
{
  int k,neirow,neicol,l=0;
  double rsum,rn;
  double totsum=NONE;
  TYPE2 *neis[8];
  int success=0;
  
  //printf("\nrow %d col %d; ", row, col);
  //collect everybody's fitness (except focal grid-cell, i.e. k=0) and keep pointer to who they are
  for(k=1;k<=8;k++){
    //Notice that if this is empty, next time we are overwriting this address.
    //If it is all empty, then l=0, which is excluded later
    neis[l]=GetNeighborP(bact,row,col,k);
    if( neis[l]->val==1 ){
      GetNeighborC(bact,row,col,k,&neirow,&neicol);
      //printf("\n\tneirow %d neicol %d; ", neirow, neicol);
      
      //goes through the neighbourhood of every individual and assign them fitness
      //neis[l]->fval2=UpdateBactFitness_Nowak(bact,neis[l],neirow,neicol);
      neis[l]->fval2=UpdateBactFitness(neis[l],neirow,neicol);  
      //neis[l]->fval2=UpdateBactFitness_no_own_gain(neis[l],neirow,neicol);  
      
      totsum+=neis[l]->fval2;
      l++;
    }
  }
  //l now says how many neighbours are competing
  if(l!=0 && totsum>0.){
    //replication happens: who's the lucky one?
    //rn= totsum*genrand_real2();
    //rn= 8.*maxprob*genrand_real2(); 					//This is problematic!
    //rn= NONE+8.*maxprob*genrand_real2(); 					//THIS too
    rn= totsum * genrand_real2();			//Totally CLAIM !!!  -- totsum already includes NONE
    rsum=0.;
    for(k=0;k<l;k++){
      rsum+=neis[k]->fval2;
      //this is the lucky one
      if(rn<rsum){
        *child=*neis[k];
        success=1;
        break;
      }
    }
    /*
    if(rsum>=l*maxprob) {
      printf("rsum too large =%f\n",rsum);
      printf("%d competitors, fvals2: ",l);
      for(k=0;k<l;k++) printf("%f ", neis[k]->fval2);
      printf("\n");
    }
    */
  }
  //printf("Success: %d, Claim: %d %f\n", success,child->val,child->fval);
  return success;
}


double Mutations(double fval)
{
  if(genrand_real1() < kmut) fval += kdelta*(genrand_real1() - 0.5);
  if( fval<0. ) fval *= -1.;
  if( fval>k_max ) fval = 2.*k_max - fval;
  
  return fval;
}    



/*********************************/
/****** NEXT STATE FUNCTION ******/
/*********************************/

void NextState(int row,int col)
{
  
  //static int count=0;
  //printf("%d %d\n",count,row);
  //count++;
  
  //if empty spot replication happens
  if(bact[row][col].val==0){
    if( AgentBasedReplication(bact,row,col,&bact[row][col]) == 1 ) bact[row][col].fval = Mutations(bact[row][col].fval);  
    //if( ClaimReplication(bact,row,col,&bact[row][col]) == 1 ) bact[row][col].fval = Mutations(bact[row][col].fval);  
  }
  //diffusion
  else if(bact[row][col].val==1 && genrand_real1()<kdiff) DiffBySwap(bact,row,col);
  
  //death is done outside next state function
}



/****************************/
/****** DEATH FUNCTION ******/
/****************************/

void Death(void)
{
  int i,j;
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){
    if(bact[i][j].val==1){
      if(genrand_real1()<death) bact[i][j]=empty;
    }
  }
  
  /*
    for(i=1;i<=nrow;i++){
      if(bact[i][ncol].val==1) {
        SetAllBack(howmanycol,bact);
        break;
      }
    }
  */
}

/*****************************/
/****** UPDATE FUNCTION ******/
/*****************************/
void Update(void)
{
  double avfval;
  
  if(Time % time_save_movie==0){
    avfval=ToMovie(bact);
    if(avfval==-1.) {
      fprintf(stderr,"Simulation finished due to extinction\n");
      exit(1);
    }else if(avfval==-2.){
      fprintf(stderr,"Update(): Warning, memory for png not allocated\n");
    }
  }
  if(Time % time_save_sample==0){
    //printf("%d %f\n",Time,avfval); //does it really matter if in sync or not?
    ToSampleData(bact);
  }
  
  
  if(Time % time_backup==0 && Time >0){
    ToBackup(bact);
  }
  
  
  //DiffBySwap(bact);
  //Synchronous(1,bact); //Reproduction(bact);
  Asynchronous();
  Death();
  //while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement
  
}






/******************************************/
/*******                            *******/
/*   ---   \\\   OUTPUTTING   |||   ---   */
/*******                            *******/
/******************************************/
void ToBackup(TYPE2 **bact)
{
  int i,j;
  //double max_fval_we_color = (k_max*100. - 10.)/100.; /*Remember that colour indexes in the gradient start at 10*/
  
  FILE *fp;
  char fname[BUFSIZE];
  char command[512];
  
  sprintf(command,"%s %s","mkdir -p",par_backup_directory_name);
  if(system(command)==-1){
    fprintf(stderr,"ToMovie: Failed to mkdir %s. Save here\n",par_backup_directory_name);
    sprintf(fname,"backup_t%d.txt",Time);
  }
  else sprintf(fname,"%s/backup_t%d.txt",par_backup_directory_name,Time);
  
  fp=fopen(fname,"w");
  fprintf(fp,"Time %d nrow %d ncol %d k_max %f\n", Time,nrow,ncol,k_max);
  
  for(i=1; i<=nrow; i++) for(j=1; j<=ncol; j++){
      if(bact[i][j].val==0) fprintf(fp,"-1.\n");
      else fprintf(fp,"%f\n",bact[i][j].fval);
  }
  fclose(fp);
}


int GetColorIndexFrom(int val,double fval)
{
  int color;
  double max_fval_we_color;
  
  if(val==0) return 0;
  max_fval_we_color = (k_max*100. - 10.)/100.; /*Remember that colour indexes in the gradient start at 10*/
  if(fval> max_fval_we_color ) color=100;
  else color= (int)( 100.* fval/k_max )+10 ;
  return color;
}

double ToMovie(TYPE2 **bact)
{
  int popcount=0;
  double avfval=0.;
  //double max_fval_we_color = (k_max*100. - 10.)/100.; /*Remember that colour indexes in the gradient start at 10*/
  
  int i,j,color_index;
  int **tomovie;
  
  if(  NULL == ( tomovie=(int **)malloc( (nrow+1) *sizeof(int*)) )  ){
   fprintf(stderr,"SaveMovie(): Error. Memory allocation unsuccessful.\n");
   return -2.;
  }
  for(i=0; i<=nrow;i++){
    if(NULL==(tomovie[i]=(int *)malloc( (ncol+1) *sizeof(int)))) {
      fprintf(stderr,"SaveMovie(): Error. Memory allocation unsuccessful.\n");
      return -2.;
    }
  }
  //here you'll save data
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){ 
    if(bact[i][j].val==0) color_index=0;
    else{
      popcount++;
      avfval+=bact[i][j].fval;
      //now we colour krec
      color_index=GetColorIndexFrom(bact[i][j].val,bact[i][j].fval); //this is function that makes colours from some features...
    }
    //color_index=GetColorIndexFrom(world[i+1][j+1].val,world[i+1][j+1].kcat); //this is function that makes colours from some features...
    tomovie[i][j]=color_index;
  }
  
  PlanePNG(tomovie,0);
  
  for(i=0; i<nrow+1;i++) free(tomovie[i]);
  free(tomovie);
  
  //printf("%d\n",popcount);
  
  return (popcount>0)?(avfval/(double)popcount):-1. ;	// Return average fval, or -1 if extinction
}

void ToSampleData(TYPE2 **bact)
{
  int i,j,k=0;
  double *presample;
  FILE *fp;
    
  if(  NULL == ( presample=(double *)malloc( (nrow*ncol) *sizeof(double)) )  ){
   fprintf(stderr,"ToSampleData(): Error. Memory allocation unsuccessful.\n");
   return;
  }
  
  for(i=1; i<=nrow; i++) for(j=1; j<=ncol; j++){
    if(bact[i][j].val>0){
      presample[k]=bact[i][j].fval;
      k++;
    }
  }
  
  // Now k is the number of individuals in the field
  // if there are fewer individuals than the sample size, we dump all of them
  fp = fopen(par_data_file_name,"a");
  if(k<=sample_size){
    for(i=0;i<k;i++) fprintf(fp,"%d %f\n",Time,presample[i]);
  }else{
    //else we have to sample the file:
    for(i=0;i<sample_size;i++){
      //generate random integer number in [0,k[
      j=(int)(k*genrand_real2());
      //pick and print that guy to file
      fprintf(fp,"%d %f\n",Time,presample[j]);
      //copy last element of presample to that spot
      presample[j]=presample[k];
      //decrease counter of k
      k--;
    }
  }
  
  free(presample);
  fclose(fp);
  
}

void SetAllBack(int howmanycol,TYPE2 **bact)
{
  int row,col;
  for(row=1; row<=nrow; row++) for(col=1; col<=ncol; col++){
    if(col<=howmanycol) bact[row][col]=bact[row][ncol-howmanycol+col];
    else bact[row][col]=empty;
  }  
  for(col=1; col<=ncol; col++){
    bact[0][col]=boundaryvalue2;
    bact[nrow+1][col]=boundaryvalue2;
  }
}



/**********************************************/
/*******                                *******/
/*   ---   \\\   INITIALISATION   |||   ---   */
/*******                                *******/
/**********************************************/

/*Individuals are copied rom backup to current field starting from upper left corner
  if file dimension is larger than current field, we cut out the excedence.
  if file dimension is smaller, the rest will be left empty.
*/
void InitialiseFromBackup(TYPE2 **bact, TYPE2 empty)
{
  int fr,fc,nr,nc,i,j;
  double init_fval;
  char instr[BUFSIZE];
  FILE *fin;
  
  fin=fopen(par_input_file_name1,"r");
  if(fin==NULL){
    fprintf(stderr, "InitialiseFromBackup(): Error. Can't open %s\n",par_input_file_name1);
    exit(1);
  }
  
  // This is the header of backup file 
  //fprintf(fp,"Time %d nrow %d ncol %d k_max %f\n", Time,nrow,ncol,k_max);
  
  //read heading of backup file
  fgets(instr,BUFSIZE,fin);
  strtok(instr," ");	//strip #Time
  strtok(NULL," ");
  strtok(NULL," ");
  fr = atoi(strtok(NULL," ")); //read # of row
  strtok(NULL," "); // strip 
  fc = atoi(strtok(NULL," \n")); //read # of col
    
  fprintf(stderr,"Initialise from backup file. #row: %d, #col: %d\n",fr,fc);
  fprintf(stderr,"Into a field of size. #row: %d, #col: %d\n",nrow,ncol);
  //nr is the largest between the grid dimension and the file dimension
  // we iterate over nr and nc making sure we do not iterate outside the file.
  // however we cannot excede the current field size either.
  nr = (nrow>fr) ? nrow : fr;
  nc = (ncol>fc) ? ncol : fc;
  for(i=1;i<=nr;i++){
    for(j=1;j<=nc;j++){
      instr[0]='\0';
      //we must read all the file lines, even if we do nothing with them, otherwise we lose track of where they should go in the actual field
      if(i<=fr && j<=fc){
        //read a line 
        if(fgets(instr,BUFSIZE,fin)==NULL){
	      fprintf(stderr, "InitialiseFromBackup(): Error at reading file\n");
	      exit(1);
	    }
	  }
	  //if we are within field size we may do something with them or not
	  //in any case we must initialise the place to empty
      if(i<=nrow && j<=ncol){
	    bact[i][j]=empty; //initialise spot
  		//file line is eihter -1. -for empty- or fval
  		if(instr[0]!='\0'){
	      init_fval = atof(strtok(instr,"\n"));
	      //printf("Hello\n");
	      if( init_fval>=0. ){
	        bact[i][j].val = 1;
	        bact[i][j].fval = init_fval;
	      }
	    }
	  }
	} 
  }	
}

// Read arguments, if any, displays defaults if mistakes
int MyArguments(int argc_g, char *argv_g[])
{
  int i,iarg;
  double darg;
  char *sarg;
  long int liarg;
  
  //printf("argc: %d\n\n\n", argc_g);
  
  for (i = 1; i < argc_g; i++) {
    if(strcmp(argv_g[i],"-rowcol") == 0){
      i++; if(i == argc_g) return 1;
      iarg = atoi(argv_g[i]);
      if(iarg<=1) { 
        fprintf(stderr,"MyArguments(): Error, nrow sould be larger than 1\n");
        exit(1);
      }else nrow=iarg;
      
      i++; if(i == argc_g) return 1;
      iarg = atoi(argv_g[i]);
      if(iarg<=1){ 
        fprintf(stderr,"MyArguments(): Error, ncol should be larger than 1\n");
        exit(1);
      }else ncol=iarg;
    }
    else if(strcmp(argv_g[i],"-bnf") == 0){
      i++; if (i == argc_g) return 1;
      darg = atof(argv_g[i]);
      if(darg<=0.){
        fprintf(stderr,"MyArguments(): Error, benefit should be larger than 0\n");
        exit(1);
      }else bnf=darg;
    }
    else if(strcmp(argv_g[i],"-cst") == 0){
      i++; if (i == argc_g) return 1;
      darg = atof(argv_g[i]);
      if(darg<=0.){
        fprintf(stderr,"MyArguments(): Error, cost should be larger than 0\n");
        exit(1);
      }else cst=darg;
    }
    else if(strcmp(argv_g[i],"-init_k") == 0){
      i++; if (i == argc_g) return 1;
      darg = atof(argv_g[i]);
      if(darg<=0. || darg>=10.){
        fprintf(stderr,"MyArguments(): Error, initial production should be larger than 0 and smaller than 10\n");
        exit(1);
      }else init_fval=darg;
    }
    else if(strcmp(argv_g[i],"-NONE") == 0){
      i++; if (i == argc_g) return 1;
      darg = atof(argv_g[i]);
      if(darg<0.){
        fprintf(stderr,"MyArguments(): Error, probability that nothing happens should be larger than or equal to 0\n");
        exit(1);
      }else NONE=darg;
    }
    else if(strcmp(argv_g[i],"-death") == 0){
      i++; if (i == argc_g) return 1;
      darg = atof(argv_g[i]);
      if(darg<0.){
        fprintf(stderr,"MyArguments(): Error, death rate should be larger than or equal to 0\n");
        exit(1);
      }else death=darg;
    }
    else if(strcmp(argv_g[i],"-kdiff") == 0){
      i++; if (i == argc_g) return 1;
      darg = atof(argv_g[i]);
      if(darg<0.){
        fprintf(stderr,"MyArguments(): Error, rate of diffusion should be larger than or equal to 0\n");
        exit(1);
      }else kdiff=darg;
    }
    else if(strcmp(argv_g[i],"-seed") == 0){
      i++; if (i == argc_g) return 1;
      liarg = atol(argv_g[i]);
      if(liarg<=0){
        fprintf(stderr,"MyArguments(): Error, seed for random number generator should be larger than 0\n");
        exit(1);
      }else ulseedG=liarg;
    }
    else if(strcmp(argv_g[i],"-selection") == 0){
      i++; if (i == argc_g) return 1;
      darg = atof(argv_g[i]);
      if(darg<0. || darg>1.){
        fprintf(stderr,"MyArguments(): Error, chose intensity of selection 0 <= w <= 1\n");
        exit(1);
      }else w=darg;
    }
    
    else if(strcmp(argv_g[i],"-time_save_sample") == 0){
      i++; if (i == argc_g) return 1;
      iarg = atoi(argv_g[i]);
      if(iarg<0){
        fprintf(stderr,"MyArguments(): Error, chose sample saving period > 0\n");
        exit(1);
      }else time_save_sample=iarg;
    }
    else if(strcmp(argv_g[i],"-time_save_movie") == 0){
      i++; if (i == argc_g) return 1;
      iarg = atoi(argv_g[i]);
      if(iarg<0){
        fprintf(stderr,"MyArguments(): Error, chose movie saving period > 0\n");
        exit(1);
      }else time_save_movie=iarg;
    }
    else if(strcmp(argv_g[i],"-time_backup") == 0){
      i++; if (i == argc_g) return 1;
      iarg = atoi(argv_g[i]);
      if(iarg<0){
        fprintf(stderr,"MyArguments(): Error, chose backup saving period > 0\n");
        exit(1);
      }else time_backup=iarg;
    }
    
    else if(strcmp(argv_g[i],"-backup") == 0){
      i++; if(i==argc_g) return 1;
      sarg = argv_g[i];
      if ((strlen(sarg) > 0) && (strlen(sarg) <= BUFSIZE-1)) 
	strcpy(par_input_file_name1, sarg);
      else {
	fprintf(stderr,"Args: file name should be 1 to 128 in length\n");
	exit(1);
      }
      continue;
    }
    else{
      return 1;
    }
  }
  return 0;
}

void PrintTheWelcomeStatement(void)
{
  fprintf(stderr,"\n\n\t ---||| PUBLIC GOOD |||---\n\n ---|||  Spatial Evolutionary Dynamics of altruistic traits |||---\n\n");
  fprintf(stderr,"List of options:\n");
  fprintf(stderr,"-rowcol : default=%d %d\n",nrow,ncol);
  fprintf(stderr,"-bnf (benefit) : default=%f\n",bnf);
  fprintf(stderr,"-cst (cost) : default=%f\n",cst);
  fprintf(stderr,"-init_k (initial_production) : default=%f\n",init_fval);
  fprintf(stderr,"-NONE (probability that nothing happens) : default=%f\n",NONE);
  fprintf(stderr,"-death (death rate) : default=%f\n",death);
  fprintf(stderr,"-kdiff (diffusion rate) : default=%f\n",kdiff);
  fprintf(stderr,"-seed (for random number gerenator) : default=%ld\n",ulseedG);
  fprintf(stderr,"-selection (intensity of selection) : default=%f\n",w);
  fprintf(stderr,"-backup (start simulations from bakup file) : Not default\n");
  fprintf(stderr,"-time_save_sample : default=%d\n",time_save_sample);
  fprintf(stderr,"-time_save_movie : default=%d\n",time_save_movie);
  fprintf(stderr,"-time_backup : default=%d\n",time_backup);
  fprintf(stderr,"\n");
}






