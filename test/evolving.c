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

static TYPE2** Comp;		// Normal Plane
static TYPE2** ColorBirth;	// Plane that displays evolved birthrates
static TYPE2** ColorDeath;	// Plane that displays evolved deathrates
static TYPE2** ColorMap;	// Plane that displays colorbar

static TYPE2 empty={0,0,0,0,0,0.,0.,0.,0.,0.}; // Predefine empty cell

double init_birth_rate=0.35;	// Starting birthrate
double init_death_rate=0.2;	// Startint deathrate
double deathrate = 0.1;		// Minimal deathrate 
double mut_rate= 0.005;		// Chance of mutations
double mut_step= 0.2;		// Size of mutations
double max_fval=1.;		// Upper bound (if applicable) for your evolvable parameter

void Initial(void)
{
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = 200; /* # of row (default=400)*/
  ncol = 200; /* # of column (default=400)*/
  nplane = 4; /* # of planes (default=0)*/
  margin = 10;
  scale = 1; /* size of the window (default=2)*/
  boundary = WRAP; /* the type of boundary: FIXED, WRAP, ECHO (default=WRAP). Note that
		      Margolus diffusion is not supported for ECHO. */
  ulseedG = 56; /* random seed (default=56)*/

  /* useally, one does not have to change the followings */
  /* the value of boundary (default=(TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.})*/
  boundaryvalue2 = (TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.}; 
}

void InitialPlane(void)
{
  MakePlane(&Comp,&ColorBirth,&ColorDeath,&ColorMap); // Init all planes

  /* InitialSet(1,2,3,4,5)
    1: name of plane
    2: number of the state other than the background state
    3: background state or empty state
    4: state I want to put
    5: fraction of cells that get S1 state (0 to 1)
  */
  
  //InitialSet(Comp,1,0,1,0.5);
  //InitialSet(Color,1,0);
  Boundaries2(Comp);			// Init boundaries for planes 
  Boundaries2(ColorBirth);		// Init boundaries for planes 
  Boundaries2(ColorDeath); 		// Init boundaries for planes 
  Boundaries2(ColorMap);		// Init boundaries for planes 
  

  // PART I: Store everyones birth/deathrates in fval/fval2
  // PART II: Make colorbar
  int i,j;
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++)
    {
// PART I
    Comp[i][j]=empty;
    if(genrand_real1()<0.2){
      Comp[i][j].val=1;
      Comp[i][j].fval=init_birth_rate;
      Comp[i][j].fval2=init_death_rate;
    }
// PART II
    if(j>10 && j < 30) ColorMap[i][j].val= 110 - (float)i/(float)nrow*100;	//maps fval on Comp plane to val on Color plane
  }


  
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
  TYPE2 nei=empty;
  nei = RandomMooreS8(Comp,row,col);	//this can be a pointer, see RandomMooreP8(Comp,row,col);
  
  if(Comp[row][col].val==0){
    if(nei.val==1 &&  genrand_real1()<nei.fval ) {
      Comp[row][col]=nei;				//birth

	// Mutate birth
      if(genrand_real1()<mut_rate) Comp[row][col].fval += mut_step*(genrand_real2()-0.5);	//mutations
      if(Comp[row][col].fval<0.) Comp[row][col].fval = (-1)*Comp[row][col].fval;		//bounds mutations to decent values
      if(Comp[row][col].fval>1.) Comp[row][col].fval = 2. - Comp[row][col].fval;
	// Mutate death
      if(genrand_real1()<mut_rate) Comp[row][col].fval2 += mut_step*(genrand_real2()-0.5);	//mutations
      if(Comp[row][col].fval2<0.) Comp[row][col].fval2 = (-1)*Comp[row][col].fval2;		//bounds mutations to decent values
      if(Comp[row][col].fval2>1.) Comp[row][col].fval2 = 2. - Comp[row][col].fval2;

    }

  }
  else 
    {
     if(genrand_real1()<(Comp[row][col].fval2+deathrate)) Comp[row][col]=empty; //death
    }
  

}
// This function maps fvals stored in individuals to a colorvalue in the colorRGB table
int GetColorIndexFrom(int val,double fval)
{
  int color;
  double max_fval_we_color;
  
  if(val==0) return 0;
  max_fval_we_color = (max_fval*100. - 10.)/100.; //Remember that colour indexes in the gradient start at 10
  if(fval> max_fval_we_color ) color=100;
  else color= (int)( 100.* fval/max_fval )+10 ;
  return color;
}

void Update(void)
{
  Display(Comp,ColorBirth,ColorDeath,ColorMap);  
  Synchronous(1,Comp);
  


// ************************************************************************************//
// GRAPHING Below is some code to plot fvals...
	// 1) ...over time for 120 random individuals every 10 timesteps
	// 2) ...as a color in the second plane
// ************************************************************************************//
  int i, j;
if(Time%10==0) {

   int r;
  for(r=1;r<=120;r++)
	{
          i = genrand_int(1,nrow);	
	  j = genrand_int(1,ncol);	
	  if(Comp[i][j].val > 0)
	  {
            PlotXY(Time,Comp[i][j].fval); 
            PlotXY(Time,Comp[i][j].fval2); 
	  }
	}
  // Setup display of graph (once)

}


if(Time == 0) 
{
	GracePrintf("g0 on");
	GracePrintf("with g0");
	GracePrintf("focus g0");
	GracePrintf("s1 on");
	GracePrintf("s1 color 5");
	GracePrintf("s1 line linestyle 0");
	GracePrintf("s1 symbol 1");
	GracePrintf("s1 symbol size 0.2");
}

for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){
    ColorBirth[i][j].val=GetColorIndexFrom(Comp[i][j].val, Comp[i][j].fval );	//maps fval on Comp plane to val on Color plane
    ColorDeath[i][j].val=GetColorIndexFrom(Comp[i][j].val, Comp[i][j].fval2 );	//maps fval on Comp plane to val on Color plane
}

// ************************************************************************************//
// END PLOTTING PART
// ************************************************************************************//

}
