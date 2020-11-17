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

static TYPE2** Cheater;		// Normal Plane
static TYPE2** ColorPars;	// Plane that displays evolved birthrates
static TYPE2** ColorMap;	// Plane that displays colorbar

static TYPE2 empty={0,0,0,0,0,0.,0.,0.,0.,0.}; // Predefine empty cell

double init_birth_rate=0.35;	// Starting birthrate
double init_death_rate=0.2;	// Startint deathrate
double deathrate = 0.01;		// Minimal deathrate 
double max_fitness = 0.0;
double mut_rate= 0.0;		// Chance of mutations
double mut_step= 0.1;		// Size of mutations
double max_fval=3.;		// Upper bound (if applicable) for plotting the fval on the grid

void Initial(void)
{
        MaxTime = 2147483647; /* default=2147483647 */
        nrow = 200; /* # of row (default=400)*/
        ncol = 200; /* # of column (default=400)*/
        nplane = 3; /* # of planes (default=0)*/
        margin = 10;
        scale = 3; /* size of the window (default=2)*/
        boundary = WRAP; /* the type of boundary: FIXED, WRAP, ECHO (default=WRAP). Note that
                            Margolus diffusion is not supported for ECHO. */
        ulseedG = 58; /* random seed (default=56)*/

        /* useally, one does not have to change the followings */
        /* the value of boundary (default=(TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.})*/
        boundaryvalue2 = (TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.}; 
}

void InitialPlane(void)
{
        MakePlane(&Cheater,&ColorPars,&ColorMap); // Init all planes

        /* InitialSet(1,2,3,4,5)
1: name of plane
2: number of the state other than the background state
3: background state or empty state
4: state I want to put
5: fraction of cells that get S1 state (0 to 1)
         */

	ColorRGB(3,0,0,255);	// I make the cheater blue since I'm RG-colourblind :D
        // PART I: Store everyones birth/deathrates in fval/fval2
        // PART II: Make colorbar
        int i,j;
        for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++)
        {
                // PART I
                Cheater[i][j]=empty;
		double kans = genrand_real1();
                if(kans<0.3333){
                        Cheater[i][j].val=1;
                        Cheater[i][j].fval=0.4;		// Helping rate 	(nonevolvable)
			Cheater[i][j].fval2=0.4;		// Receiving rate	(evolvable)
                }
                else if(kans<0.6666){
                        Cheater[i][j].val=2;
                        Cheater[i][j].fval=0.4;		// Helping rate		(nonevolvable)
			Cheater[i][j].fval2=0.4;		// Receiving rate	(evolvable)
                }
                else{
                        Cheater[i][j].val=3;
			Cheater[i][j].fval2=0.6;	// Receiving rate cheater 1	(evolvable)
                }
                // PART II
                if(j>10 && j < 30) ColorMap[i][j].val= 110 - (float)i/(float)nrow*100;	//maps fval on Cheater plane to val on Color plane
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
        if (Cheater[row][col].val==0)
        {
		//printf("Start Cheater\n");
         	TYPE2 sum = SumMooreS8(Cheater, row, col);
	        TYPE2 nei=empty;
		double sum_fit = sum.fval3;
		double kans = genrand_real1();
		double pinv = 0.0;
		for(int i =1; i <9; i++)	
                {
			nei = GetNeighborS(Cheater,row,col,i);
			pinv += nei.fval3 / (sum_fit + 2.0);
			if(kans < pinv)
			{
		                Cheater[row][col].val = nei.val;
			        Cheater[row][col].fval = nei.fval;
				Cheater[row][col].fval2 = nei.fval2;
				if(Cheater[row][col].val < 3)
				{
                                	if(genrand_real1()<mut_rate) Cheater[row][col].fval2 += mut_step*(genrand_real2()-0.5);	//mutations
	                                if(Cheater[row][col].fval2<0.) Cheater[row][col].fval2 = (-1)*Cheater[row][col].fval2;		//bounds mutations to decent values
	                                
				}
				else if(Cheater[row][col].val >= 3)
				{
                                	if(genrand_real1()<mut_rate) Cheater[row][col].fval2 += mut_step*(genrand_real2()-0.5);	//mutations
	                                if(Cheater[row][col].fval2<0.) Cheater[row][col].fval2 = (-1)*Cheater[row][col].fval2;		//bounds mutations to decent values
	                                
				}
				break;
			}
                }
        }
      	else{
	        if (genrand_real1()<0.20)
                {
                        Cheater[row][col].val=0;
			Cheater[row][col].fval=0.;
			Cheater[row][col].fval2=0.;	
			Cheater[row][col].fval3=0.;
                }              
        }

}

int GetColorIndexFrom(int val,double fval)
{
        int color;
        double max_fval_we_color;

        if(val==0) return 0;
        max_fval_we_color = (max_fval*100.)/100.; //Remember that colour indexes in the gradient start at 10
        if(fval> max_fval_we_color ) color=100;
        else color= (int)( 100.* fval/max_fval )+10 ;
        return color;
}

void StoreFitness(TYPE2** p)
{
	max_fitness = 0.0;
	for(int row = 1; row <= nrow; row++)
		for(int col = 1; col<=ncol; col++)
		{
			p[row][col].fval3 = 0.0; 		// Reset fitness to 0.0
			if(p[row][col].val >0)
			{			
				for(int i =1; i <9; i++)	// First get all neighbors to check who's there
				{
				        TYPE2 nei = GetNeighborS(Cheater,row,col,i);
					if(p[row][col].val == 1 && nei.val == 2)
						p[row][col].fval3 += p[row][col].fval2*nei.fval;	// X gets fitness from Y
					else if(p[row][col].val == 2 && nei.val == 1)
						p[row][col].fval3 += p[row][col].fval2*nei.fval;	// Y gets fitness from X
					else if (p[row][col].val == 3 && nei.val == 2)
						p[row][col].fval3 += p[row][col].fval2*nei.fval;	// Cheater gets fitness from Y
					else if (p[row][col].val == 4 && nei.val == 1)
						p[row][col].fval3 += p[row][col].fval2*nei.fval;	// Cheater gets fitness from Y
				}
			}	
		}
}

void Update(void)
{

	//if(Time%3==0)MDiffusion(Cheater);


        // ************************************************************************************//
        // GRAPHING Below is some code to plot fvals...
        // 1) ...over time for 120 random individuals every 10 timesteps
        // 2) ...as a color in the second plane
        // ************************************************************************************//
        int i, j;
	//Plot(1,Cheater);
        if(Time%100==0) {
		int sum1 = 0;
		int sum2 = 0;
		int sum3 = 0;

		double sumrec1 = 0.0;
		double sumrec2 = 0.0;
		double sumrec3 = 0.0;

                for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++)
		{	
	                //if(Cheater[i][j].val > 0 && Cheater[i][j].val < 3) { PlotXY(Time,Cheater[i][j].fval); }
	                //PlotXY(Time,Cheater[i][j].fval2);                         
			if(Cheater[i][j].val == 1) sumrec1+=Cheater[i][j].fval2, sum1++;
			if(Cheater[i][j].val == 2) sumrec2+=Cheater[i][j].fval2, sum2++;
			if(Cheater[i][j].val == 3) sumrec3+=Cheater[i][j].fval2, sum3++;
                }
		if(Time==0)printf("Alpha 1\t\tAlpha 2\t\tBeta\t\t(you can use these data to make a pretty plot if you want)\n");
		
		printf("%f\t%f\t%f\n",sumrec1/sum1,sumrec2/sum2,sumrec3/sum3);
                // Setup display of graph (once)
        }
	Plot(1,Cheater);

        for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++)
	{
                ColorPars[i][j].val=GetColorIndexFrom(Cheater[i][j].val, Cheater[i][j].fval );	//maps fval on Cheater plane to val on Color plane
        }

        // ************************************************************************************//
        // END PLOTTING PART
        // ************************************************************************************//

       	if(Time%3==0)Display(Cheater,ColorPars,ColorMap);  				// Only display every 3 timesteps
	//if(Time%100 == 0) DrawSlide2("images/",4,Cheater,ColorPars,ColorDeath,ColorMap);	// Saves a nice movie :)
	StoreFitness(Cheater);	// Stores fitness for everyone before updating
        Synchronous(1,Cheater);	// Update on the basis of fitness
}
