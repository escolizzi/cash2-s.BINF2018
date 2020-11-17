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

#define MAXR 10     //This must be always less than or equal to the size of replicator array in cash2.h

void CleanGridPointInfo(TYPE2 *ic);
double SplittingRule(TYPE2 *ic);
int UpdateGridPointInfo(TYPE2 *ic);
double Mutate_krec( double );
double Mutate_tar_n( double );
int ExternalReplication(TYPE2 **world, int row, int col);
int InternalReplication(TYPE2 **world, int row, int col);
int ComplexFormationDissociation(TYPE2 **world, int row, int col);
void PrintCellInfo(TYPE2 *ic,int i,int j);
void ScrambleIntArray(int *arr,int size);
int IsInArray(int *arr, int size ,int id);
void InitialiseUpdorderarray(struct updorder **updorderDiff,int nrow,int ncol);
void ScrambleUpdorderarray(struct updorder *updorderDiff);


static TYPE2** world;
double kmut_krec = 0.; //0.01;
double kmut_tar_n = 0.; //0.05;
double kdeath = 0.05;  //death rate
double kdiss = 0.2;   //complex dissociation rate
double kdiff = 0.3;   // diffusion rate for single molecules (for vescicle it is kdiff/ic->val)
double k_compl_factor_diff = 0.85;  //factor to multiply kdiff for complexes

static struct updorder *updorderDiff;

void Initial(void)
{
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 10; /* # of row (default=100)*/
  ncol = 10; /* # of column (default=100)*/
  nplane = 1; /* # of planes (default=0)*/
  scale = 2; /* size of the window (default=2)*/
  display = 0;  //display=0 -> no display
  boundary = WRAP; /* the type of boundary: FIXED, WRAP, ECHO (default=WRAP). Note that
		      Margolus diffusion is not supported for ECHO. */
  ulseedG = 56; /* random seed (default=56)*/

  /* useally, one does not have to change the followings */
  /* the value of boundary (default=(TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.})*/
  boundaryvalue2 = (TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.}; 
}

void InitialPlane(void)
{
  int i,j;
  TYPE2 *ic;
  MakePlane(&world);
  
  InitialiseUpdorderarray(&updorderDiff,nrow,ncol); //we pass the address of the pointer to modify it inside the function
  //printf("CHECK INIT AND USAGE OF updorderDiff... %d\n",updorderDiff[0].col);
  //exit(1);
  
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
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){
    CleanGridPointInfo(&world[i][j]);
  }
  //for the sake of initial testing:
  //initilise a vescicle by hand, which is ready to split:
  ic = &world[nrow/2][ncol/2];
  
  for(i=0;i<5;i++){
    ic->replicator[i].id=i;
    ic->replicator[i].krec=0.15*i;
    ic->replicator[i].tar_n=5.;
  }
  //make them in complex by hand:
  ic->replicator[0].c_id=1; ic->replicator[0].replicase = 1;
  ic->replicator[1].c_id=0; ic->replicator[1].replicase = 0;
  
  ic->replicator[2].c_id=3; ic->replicator[2].replicase = 1;
  ic->replicator[3].c_id=2; ic->replicator[3].replicase = 0;
  
  int bla=UpdateGridPointInfo(ic);
  printf("CHECK: vescicle val=%d and fval=%f\n",ic->val,ic->fval);
  //exit(0);
  Boundaries2(world);
}

// splitting rule, what it the fval at which a vescicle splits?
//the minimum? the average? a constant?
// ... for now, the average
double SplittingRule(TYPE2 *ic)
{
  int i;
  double splitsize=0.;
  for(i=0; i < ic->val; i++) splitsize += (ic->replicator[i].tar_n > 0.)?ic->replicator[i].tar_n:0.;
  return splitsize/(double)(ic->val);
}

//just set a cell of the CA to empty
void CleanGridPointInfo(TYPE2 *ic)
{
  int i;
  ic->val=0;
  ic->fval=0;
  for(i=0;i<MAXR;i++){
    ic->replicator[i].id=-1;
    ic->replicator[i].c_id=-1;
    //ic->replicator[i].replicase=-1;
  }
}

//Updates val and fval (i.e. vescicle size and average target splitting size)
//Updates id inside the replicators array
//be careful with complexes
int UpdateGridPointInfo(TYPE2 *ic)
{
  int i,j,howmany=0;
  //first count how many replicators are in there
  for(i=0;i<MAXR;i++){
    if( ic->replicator[i].id != -1 ) howmany++;
  }
  ic->val=howmany;  //set vescicle size val
  
  //printf("UpdateGridPointInfo()\n");
  
  //After replication there will be holes in the array. 
  //We remove them by restacking the array, being careful about complexes.
  j=MAXR-1; // define it here, so that we can keep going with the loop from the right
  for(i=0;i<ic->val;i++){
    if( ic->replicator[i].id == -1 ){
      //this is a hole, we fill it by taking the last guy of the array an put it here
      for(;j>0;j--){
        if(ic->replicator[j].id != -1){
          //printf("UpdateGridPointInfo(): check id, krec and tar_n of guy at pos j=%d: %d, %f, %f\n",j,ic->replicator[j].id,ic->replicator[j].krec,ic->replicator[j].tar_n);
          //copy this to position i, reassign id, careful with complex
          ic->replicator[i] = ic->replicator[j];
          ic->replicator[i].id=i;
          ic->replicator[j].id=-1; ic->replicator[j].c_id=-1;
          if(ic->replicator[i].c_id != -1){
            int complex_id;
            complex_id = ic->replicator[i].c_id;
            ic->replicator[complex_id].c_id = i;
          }
          break;  //break out of the inner loop to let the outer one can proceed
        }
      }
    }
  }
  
  ic->fval=SplittingRule(ic); //This re-determines the size at which a vescicle will break
  return 0;
}

double Mutate_krec( double krec )
{
  if(genrand_real2()<kmut_krec)
    krec += (genrand_real1() - 0.5)/10.;
  return krec;
}

double Mutate_tar_n(double tar_n)
{
  if(genrand_real2()<kmut_tar_n)
    tar_n += (genrand_real1() - 0.5)/5.;
  return tar_n;
}

//Replication works that you look only in your own grid point, 
// see who is in complex, and have some probability of replication depending on that
// if replication happens in a bin that exceeds some function of the average target number for splitting, the whole thing splits.
int ExternalReplication(TYPE2 **world, int row, int col)
{
  int outcome=0,global_replincomplex=0,neicounter,localcounter,which;
  int nr,nc;
  TYPE2 *ic,*nei;
  
  //first check if there are vescicles that are ready to split 
  // this includes "lone" pairs of replicators
  // if there are, check if anyone is in complex inside and can replicate
  // for now we just count how many???
  for(neicounter=1;neicounter<=8;neicounter++){
    nei=GetNeighborP(world,row,col,neicounter);
    if(nei->val >= (int)(0.5+nei->fval) ){ 
      int localcounter;
      for(localcounter=0; nei->replicator[localcounter].id != -1; localcounter++){
        if(nei->replicator[localcounter].c_id != -1 && nei->replicator[localcounter].replicase) 
          global_replincomplex ++; 
      }
    }
  }
  
  //second pass: make one of them replicate (if there is at least on of them)
  which = 1+ (int)(global_replincomplex * genrand_real2());
  if(global_replincomplex>0){
    for(neicounter=1;neicounter<=8;neicounter++){
      GetNeighborC(world,row,col,neicounter, &nr, &nc);
      nei = &world[nr][nc];
      if(nei->val >= (int)(0.5+nei->fval) ){
        for(localcounter=0; nei->replicator[localcounter].id != -1; localcounter++){
          if(nei->replicator[localcounter].c_id != -1 && nei->replicator[localcounter].replicase) 
            which --;
          if(which == 0){
            outcome=1;
            break;
          }
        }
      }
      if(outcome==1) break;
    }
  }
  
  if(outcome==1){
    //printf("External replication happening into position %d, %d\n", row,col);
    //this one replicates whoever it is in complex with, and carries along everybody else:
    
    //1) internal replication - we do not update all the statistics because we will do it after we break cell
    int template_pos = nei->replicator[localcounter].c_id;
    int child_pos = nei->val; //we put child molecule at last position of array
    nei->replicator[child_pos] = nei->replicator[template_pos];
    //2) Mutations
    nei->replicator[child_pos].krec = Mutate_krec( nei->replicator[child_pos].krec );
    nei->replicator[child_pos].tar_n = Mutate_tar_n( nei->replicator[child_pos].tar_n );
    //3) Undo complex
    nei->replicator[child_pos].c_id=-1;   //child replicator not in complex
    nei->replicator[template_pos].c_id=-1;//parent replicator not in complex
    nei->replicator[localcounter].c_id =-1; //replicase not in complex
    //printf("Replicated id %d into position %d\n",localcounter,nei->val);
    //break cell: I have an array long val+1 -> I want to split it binomially, 
    //do I just do the bernoulli trials?... for now yes
    //also, do I keep complexes together? I should...
    ic=&world[row][col];
    CleanGridPointInfo(ic);   //this should be already clean
    for(localcounter=0;localcounter < nei->val+1; localcounter++){
      //it can be that id ==-1 even though we have not reached the end of the replicator array 
      // because the molecule that was there was in complex with one that moved
      // also if you're at the template end of complex you move only if the replicase end moves
      if(nei->replicator[localcounter].id != -1 && genrand_real2() < 0.5){
        //move it if not in complex
        if( nei->replicator[localcounter].c_id == -1){
          //printf("copying cell not in complex: localcounter %d\n",localcounter);
          ic->replicator[localcounter] = nei->replicator[localcounter];
          //printf("Check that something has been copied, copied krec: %f\n",ic->replicator[localcounter].krec);
          nei->replicator[localcounter].id=-1;  // and this guy has also moved
        }
        //we check if in complex: in order not to give complexes twice the chances of moving
        //we move complex only by the replicase
        if( nei->replicator[localcounter].c_id != -1 && nei->replicator[localcounter].replicase){
          //we copy also the other member of the complex
          //printf("copying cell in complex: localcounter %d\n",localcounter);
          ic->replicator[localcounter] = nei->replicator[localcounter];
          
          int complex_id = nei->replicator[localcounter].c_id;
          ic->replicator[complex_id] = nei->replicator[complex_id];
          nei->replicator[complex_id].id=-1;    //and we set id to -1 to signal that this has moved
          nei->replicator[localcounter].id=-1;  // and this guy has also moved
        }
      }
    }
    // update all statistics
    UpdateGridPointInfo(&world[row][col]);
    UpdateGridPointInfo(&world[nr][nc]);
    
    
    PrintCellInfo(&world[row][col],row,col);
    PrintCellInfo(&world[nr][nc],nr,nc);
    //exit(1);
  }
  
  
  return outcome;
}

//This function makes replication of a lone replicator into a spot occupied by only a lone replicator
// we take all lone replicators on the template side of complex (if any)
// and make them replicate
int LoneReplExternalReplication(world, row,col)
{
  int outcome;
  
  //first pass check who is there
  
  
  //second pass replicate one of them randomly
  printf("YOU ARE HERE\n");
  exit(1);
  
  return outcome;
}

int InternalReplication(TYPE2 **world,int row, int col)
{
  int localcounter,howmanyincompl=0,howmuchspace,i,outcome; 
  int l_id[]={[0 ... (MAXR/2 - 1) ] = -1}; // This works - tested!!!
  TYPE2 *ic;
  ic = &world[row][col];
  
  howmuchspace = (int)(0.5+ic->fval) - ic->val;
  
  if(howmuchspace <= 0) return 0;   //if the vescicle is full, just return
  
  printf("Internal replication might happen\n");
  //else there is space, and people in complex will replicate
  //1) see who can replicate - because it is in complex
  // save the id of these people in l_id, 
  // and later randomize where they come from
  for(localcounter=0;localcounter < ic->val; localcounter++){
    if(ic->replicator[localcounter].c_id !=-1 && ic->replicator[localcounter].replicase==0){
      l_id[howmanyincompl]=ic->replicator[localcounter].id;
      howmanyincompl++;
    }
  }
  
  //second pass, copy all of them until either the vescicle is full, or until we run out of possibilities
  //if there is more space than people in complex then everybody can replicate
  if(howmanyincompl <= howmuchspace){
    for(i=0;i<howmanyincompl;i++){
      //everybody replicates
      int template_id = l_id[i];
      int child_pos = ic->val;
      ic->replicator[child_pos] = ic->replicator[template_id];
      //mutations
      ic->replicator[child_pos].krec = Mutate_krec( ic->replicator[child_pos].krec );
      ic->replicator[child_pos].tar_n = Mutate_tar_n( ic->replicator[child_pos].tar_n );
      //undo complex
      ic->replicator[child_pos].c_id=-1;
      int replicase_pos = ic->replicator[template_id].c_id;
      ic->replicator[replicase_pos].c_id=-1;
      ic->replicator[template_id].c_id=-1;
      //icrement count of people in there 
      ic->val++;
    }
  }
  else{
    //else we have to randomly choose "howmuchspace"-people out of the l_id array
    //scramble the array l_id
    for(i=0; i<howmuchspace; i++){
      int ranpos= i + (int)((float)(howmanyincompl-i)*genrand_real2());
      int tmp = l_id[i];
      l_id[i] = l_id[ranpos];
      l_id[ranpos] = tmp;
    }
    for(i=0; i<howmuchspace; i++){
      //everybody replicates
      int template_id = l_id[i];
      int child_pos = ic->val;
      ic->replicator[child_pos] = ic->replicator[template_id];
      //mutations
      ic->replicator[child_pos].krec = Mutate_krec( ic->replicator[child_pos].krec );
      ic->replicator[child_pos].tar_n = Mutate_tar_n( ic->replicator[child_pos].tar_n );
      //undo complex
      ic->replicator[child_pos].c_id=-1;
      int replicase_pos = ic->replicator[template_id].c_id;
      ic->replicator[replicase_pos].c_id=-1;
      ic->replicator[template_id].c_id=-1;
      //icrement count of people in there 
      ic->val++;
    }
  }
  
  outcome = UpdateGridPointInfo(ic);
  return 0;
  
}

int ComplexFormationDissociation(TYPE2 **world, int row, int col)
{
  int i,j;
  TYPE2 *ic;
  int l_id[] = {[0 ... MAXR] = -1};   //crazy initialisation that works
  
  ic=&world[row][col];
  if(ic->val <= 1) 
    return 0; //if there's not enough replicators to make complex, there's no point in any of this
  
  //First pass: find who can be in complex, save its id in l_id and scramble l_id
  j=0;
  for(i=0; i < ic->val; i++){
    if(ic->replicator[i].c_id == -1){
      l_id[j]=i;  //save id in l_id
      j++;
    }
    else{
      //DISSOCIATION ON THE FLY
      // if someone is in complex, in order not to upset what we are doing above, 
      // we check also the complexee
      int complex_id = ic->replicator[i].c_id;
      if( ic->replicator[complex_id].c_id == -1 ){
        //then we have already undone this complex from the other end, 
        // thus we just undo the complex at position i
        // notice that in this case it must be that i > complex_id
        // because otherwise the case below has already properly taken care of the complex
        ic->replicator[i].c_id =-1;
      }else{
        //we have not undone the complex yet, so we attempt this
        // (kdiss/2 because there are two ways of undoing complex)
        if(genrand_real2() < kdiss/2.){
          if(i>complex_id) 
            ic->replicator[complex_id].c_id=-1;  // we should undo the complexee if i > complex_id, 
                                                 // because otherwise we'll have a mismatching complex
          ic->replicator[i].c_id =-1; 
        }
      }
    }
  }
  
  ScrambleIntArray(&l_id[0],j);  //scramble
  printf("Scrambled array of who can form complex\n");
  for(i=0;i<j;i++) printf("%d ",l_id[i]);
  printf("\n");
  
  if(j<=1) 
    return 0; // if there's not enough free to make complex, return
  
  //second pass, actually make complexes - l_id is scrambled and j is its size
  //we take the id's of people that attempt complex formation from l_id 
  for(i=0; i < j; i++){
    //everyone that is not in complex attempts to be the replicase with its own krec
    int me_try = l_id[i];
    int complex_id = ( (ic->val -1) * genrand_real2() ); //with whom me_try tries to make complex (complex_id is index over ic->replicator[])
    if(complex_id >= me_try ) complex_id ++;             // because we exclude me_try itself
    //if  it is not in complex,              and RAND() < me_try.krec                         and complex_id is in l_id
    if( ic->replicator[complex_id].c_id ==-1 && genrand_real2() < ic->replicator[me_try].krec && IsInArray(&l_id[0],j,complex_id)){
      //complex formation happens
      ic->replicator[me_try].c_id = complex_id;
      ic->replicator[me_try].replicase = 1;
      ic->replicator[complex_id].c_id = me_try;
      ic->replicator[complex_id].replicase = 0;
    }
  }
  printf("Who is in complex? printing vescicle info\n");
  PrintCellInfo(ic,row,col);
  return 0;
}

void NextState(int row,int col)
{
  int outcome=0;
  
  //let's have that val indicates number of individuals per pixel (i.e. length of repl[])
  if(world[row][col].val==0){
    //RandomMooreC8(world,row,col,&rr,&rc);  //get neigh
    //if(world[rr][rc].val>1) 
    outcome=ExternalReplication(world, row,col);  //Mutations happen inside Replication functions
  }
  else{
    if(world[row][col].val==1 && world[row][col].fval<2.5 ){
      outcome=LoneReplExternalReplication(world, row,col);    // this and following replication do not interfere with each other
    }else{
      outcome=InternalReplication(world, row,col);
    }
    outcome=ComplexFormationDissociation(world, row,col);
    
  }
}

//diffusion is very important in this model, because it is the only way for a lone replicator
//to meet another one and make complex
void MyDiffusion(TYPE2 **world)
{
  int l,i,j,k,placer_for_ic,outcome;
  TYPE2 *ic,*nei;
  //fprintf(stderr,"Time = %d, ORDER OF DIFFUSION SHOULD BE RANDOMISED\n",Time);
  
  ScrambleUpdorderarray(updorderDiff);
  
  //printf("AFTER FUNCTION CALL\n");
  //int l;
  //for(l=0;l<nrow*ncol;l++)
  //  //printf("%d %d %c",updorderDiff[l].row, updorderDiff[l].col, (l%nrow)?'\n':'\t');
  //  printf("%d %d\t",updorderDiff[l].row, updorderDiff[l].col);
  //exit(1);
  
  //for(i=1;i<=nrow; i++)for(j=1;j<=nrow; j++){       // <==== THIS MUST BE RANDOMIZED !!!!!!!!!!!!!
  for(l=0;l<nrow*ncol;l++){
    i=updorderDiff[l].row;
    j=updorderDiff[l].col;
    
    ic=&world[i][j];
    //empty spot, everybody can diffuse with prob. proportional to their size
    int random_pos = 1+(int)(8.*genrand_real2());
    int rr,rc;
    GetNeighborC(world,i,j,random_pos, &rr,&rc);
    nei=&world[rr][rc];
    
    // -------------------------
    //first case: if ic is empty
    // -------------------------
    if(ic->val == 0){
      //CASE 1.0 --- if nothing in neigh, there is nothing to do
      if(nei->val==0) continue;   
      //CASE 1.1 --- if somehitng in neigh, then it depends whether there are one or two independent replicators - in compex or not, or a vescicle
      //if there is something, nei->val > 0:
      //if (int)(0.5+fval) < 2 this is not a vescicle, we move replicators independently (unless they are in complex)
      else if(nei->fval < 2.5 ){
        //in this case no vescicle is formed (yet), so replicators move independently
        //also, there should not be more than two replicators here (this is just a check)
        if(nei->val >2) printf("MyDiffusion(): Warning. ic->fval = %f but ic->val = %d. Should there be more replicators than fval says?\n", ic->fval, ic->val);
        for(k=0; k<nei->val;k++){
          placer_for_ic = 0;
          if(nei->replicator[k].id==-1) continue; //because it could be someone in complex
          //if not in complex
          if( nei->replicator[k].c_id == -1 && genrand_real2() < kdiff ){
            //simply copy nei->replicator[k] into ic->replicator[placer_for_ic]
            ic->replicator[placer_for_ic] = nei->replicator[k];
            nei->replicator[k].id = -1;
            placer_for_ic++;
            //if in complex
          }else if(nei->replicator[k].c_id != -1 && nei->replicator[k].replicase ==1 && genrand_real2() < k_compl_factor_diff*kdiff){
            //complex formation diffusion
            ic->replicator[placer_for_ic] = nei->replicator[k];
            int complex_id = nei->replicator[k].c_id;
            ic->replicator[placer_for_ic+1] = nei->replicator[complex_id];
            //flags for complex - replicase vs. template flags are taken care of because we move comple from replicase perspective
            ic->replicator[placer_for_ic].c_id = placer_for_ic+1;
            ic->replicator[placer_for_ic +1 ].c_id = placer_for_ic;
            //delete information from nei
            nei->replicator[ nei->replicator[k].c_id ].id =-1;
            nei->replicator[ nei->replicator[k].c_id ].c_id =-1;
            nei->replicator[k].c_id = -1;
            nei->replicator[k].id = -1;
            //increase placer_for_ic counter
            placer_for_ic += 2;
          }
        }
      }
      //else there is a vescicle, and should be moved all at once
      else{
        if( genrand_real2() < kdiff/(float)(nei->val) ){
          *ic = *nei;
          CleanGridPointInfo(nei);//delete from previus position
        }
      }
    // -------------------------
    //second case: if ic->val >0 and < 2 ...this is same as ic->fval==1 and it is not a vescicle
    // -------------------------
    }else if(ic->val ==1 && ic->fval < 2.5 && nei->fval < 2.5){
      //if nei is not a vescicle itself, it has at least one member and if it has two they are not in complex
      if( nei->val ==2 && nei->replicator[0].c_id==-1 && genrand_real2() < kdiff){
        int which = (int)(2.*genrand_real2());
        ic->replicator[ ic->val ] = nei->replicator[which];
        nei->replicator[which].id=-1;
      }else if(nei->val == 1 && genrand_real2() < kdiff){
        ic->replicator[ ic->val ] = nei->replicator[0];
        nei->replicator[0].id=-1;
      }
    // ------------------------
    // third case: if ic is occupied by vescicle -- should be maybe moved to top
    // ------------------------
    }else if(ic->fval >= 2.5){
      continue;
    }
    
    outcome = UpdateGridPointInfo(ic);
    outcome = UpdateGridPointInfo(nei);
    
  }
}

//DEATH OF SOME
//you'll have to re-stack the replicators (and re-assign complexes) everytime someone dies
void Death(TYPE2 **world)
{
  int i,j,k,outcome;
  TYPE2 *ic;
  
  printf("\nDeath takes us all, but what is death? Death is this routine:\n");
  
  for(i=1;i<=nrow; i++)for(j=1;j<=nrow; j++){
    ic=&world[i][j];
    for(k=0; k < ic->val; k++){
      if(genrand_real2() < kdeath){
        printf("A fallen replicator: pos %d %d, id %d, krec = %f\n", i,j,ic->replicator[k].id, ic->replicator[k].krec);
        //eliminate
        ic->replicator[k].id=-1;
        //if in complex, break it
        if(ic->replicator[k].c_id != -1){
          int complex_id = ic->replicator[k].c_id;
          ic->replicator[complex_id].c_id = -1; //break complex from the other member
          ic->replicator[k].c_id = -1;          //break complex here
        }
      }
    }
    
    outcome = UpdateGridPointInfo(ic);
  }
}

void PrintCellInfo(TYPE2 *ic,int row,int col)
{
  int k;
  printf("I'm at pos %d %d\n",row,col);
  printf("Vescicle size = %d, target split size = %f\nThis is inside me:\n",ic->val,ic->fval);
  for(k=0; k < ic->val; k++){
    printf("id: %d, krec = %f, tar_n = %f, c_id = %d, replicase = %d\n",ic->replicator[k].id,ic->replicator[k].krec, ic->replicator[k].tar_n, ic->replicator[k].c_id,ic->replicator[k].replicase);
  }
  //printf("just a check that there is nothing in replicator[] at position val: %d\n",ic->replicator[ic->val].id);
    
}

void PrintStuff(TYPE2 **world)
{
  int i,j,popsize=0;
  TYPE2 *ic;
  
  printf("\n\n\nTime: %d\n",Time);
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){
    ic = &world[i][j];
    if(ic->val > 0) {
      PrintCellInfo(ic,i,j);
      popsize+=ic->val;
    }
  }
  
  if(popsize==0){
    fprintf(stderr, "Global Extinction. Simulation terminated.\n");
    exit(0);
  }
}

void Update(void)
{ 
  PrintStuff(world);
  
  //Display(world);
  Asynchronous();
  //Synchronous(1,world);
  MyDiffusion(world);
  Death(world);
  //while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement
  
  //if(Time>5) exit(0);
  
}


// USEFUL ROUTINES
void ScrambleIntArray(int *arr,int size)
{
  int i,tmp;
  for(i=size-1; i>0; i--){
    int ran_pos = (int)((i+1)*genrand_real2() );  //random pos between 0 and i both included (i = size-1)
    tmp=arr[i];
    arr[i]=arr[ran_pos];
    arr[ran_pos]=tmp;
  }
}

int IsInArray(int *arr, int size ,int id)
{
  int i;
  for(i=0;i<size;i++){
    if(arr[i]==id) 
      return 1;
  }
  return 0;
}

//Allocate memory for the array, and initialise its values 
void InitialiseUpdorderarray(struct updorder **updorderDiff,int nrow, int ncol)
{
  int i,j;
  //struct updorder *a;
  *updorderDiff = (struct updorder *)malloc(nrow*ncol*sizeof(struct updorder));
  
  //THERE ARE CORRECTIONS TO MAKE BECAUSE i and j should go from 1 to nrow and from 1 to ncol
  for(i=0;i<nrow;i++)for(j=0;j<ncol;j++){
    (*updorderDiff)[nrow*i+j].row=i+1;
    (*updorderDiff)[nrow*i+j].col=j+1;
  }
  
  //for(i=0;i<nrow*ncol;i++)
  //  printf("%d %d\n",(*updorderDiff)[i].row, (*updorderDiff)[i].col);
  //printf("Are there corrections to make about this?\n");
  //exit(1);
}
//scramble array
void ScrambleUpdorderarray(struct updorder *updorderDiff)
{
  int l,rpos;
  l=nrow*ncol;
  struct updorder tmp;
  
  //printf("\nBEFORE SCRAMBLE\n");
  //for(l=0;l<nrow*ncol;l++)
  //  printf("%d %d\t",updorderDiff[l].row, updorderDiff[l].col);
  
  for(l=nrow*ncol-1; l>0; l--){
    rpos=(int)( l * genrand_real2()  );
    tmp = updorderDiff[rpos];
    updorderDiff[rpos] = updorderDiff[l];
    updorderDiff[l] = tmp;
  }
  
  //printf("\nAFTER SCRAMBLE\n");
  //for(l=0;l<nrow*ncol;l++)
  //  printf("%d %d\t",updorderDiff[l].row, updorderDiff[l].col);
  
}














