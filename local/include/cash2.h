/******************************************************
 ************Here you define your own type*************
 ******************************************************/

// THIS IS A REPLICATOR
typedef struct __repl {
  int id;       //id
  double krec;  //association rate
  double tar_n; //target number of replicators like self that would destabilise group
  int c_id;     //if in complex, tells id of whom it is in complex with, if not == -1
  int replicase;//if in complex: 1 if it is replicase, 0 if template
} REPLICATOR;

//THIS IS A GRID POINT, may contain many replicators (max=256)
typedef struct __type2{
  int val;
  int val2;
  int val3;
  int val4;
  int val5;
  double fval;
  double fval2;
  double fval3;
  double fval4;
  double fval5;
  char seq[1024];
  char str[1024];
  int crow;
  int ccol;
  int valarray[1024];
//  REPLICATOR replicator[256];
} TYPE2;

/*****************************************************
 *************Do not change from here*****************
 *****************************************************/

/* This structure is used in asynchronous updating to determine
   the order of the cells to update. Look at the functions in
   chash2.c which use "struct updorder" for details and the usage
   of this structure. Look at UpdOrdReset() in cash2.c  */
struct updorder {
  int row;
  int col;
};

/* This structure is used in neighborhoods retrieving. Look
   at NeighSet() in cash2.c */
struct point {
  int row;
  int col;
};

/* This structure is used in Margolus neighborhood
   retrieving. Read MarGolusNeigh(). If you say "MARGOLUS
   something_you_defined[*][i][j]" where [*] denote even-or-odd
   phase, then "----[*][i][j].m[CW].row" and
   "----[*][i][j].m[CW].col" give you the coordinate of the CW
   neighbor (see below) of Margolus neighborhood of the [i][j]
   cell in the phase [*] (even/odd, 0 or 1).

   -----------
   |HERE| CW |
   -----------
   | CCW| OPP|
   -----------

   Note that the exact place of "HERE" does not matter in that
   "CW", "OPP" and "CCW" notations are invariant for the position
   of "HERE".

   What is even and what is odd? Let "SYD" denote
   "something_you_defined". Then, SYD[0] and SYD[1] gives you tow
   alternative Margolus partitioning of the plane. Then, SYD[0]'s
   partition will be such that [0][0] cell of your plane will be
   the upper-left corner (UL) of the 2x2 square. SYD[1]'s
   partition will be such that [0][0] cell  of your plane will be
   the lower-right corner (LR) of the 2x2 square. See the
   followings too.

   UL,UR,LL,LR are like,
   -----------
   | UL | UR |
   -----------
   | LL | LR |
   -----------


   in even phase (SYD[0]), the plane is divided as follows:
   -----------------
   | 1,1 | 1,2 | etc.
   -------------
   | 2,1 | 2,2 | etc.
   -------------
   | etc. etc.   etc.

   in odd phase (SYD[1]) and boundary=WRAP, the plane is divided as follows:
   ------------------
   |nr,nc| nr,1| etc.
   -------------
   | 1,nc|  1,1| etc.
   -------------
   | etc. etc.   etc.

   where "nr" is the number of row, and "nc" is the number of
   column. If boundary==FIXED, then read the avobe table such
   that "nr=nc=0".
*/
#define CW (0)
#define OPP (1)
#define CCW (2)
typedef struct __margolus{
  struct point m[3];/* 0:clockwise 1:opposite 2:counter-clockwise */
} MARGOLUS;


/*********************************************basic*/
TYPE2 **NewP2(void);
TYPE2 **New2(void);
double **NewDB(void);
int PlaneFree2(TYPE2**);
void UpdOrdReset(struct updorder**);
TYPE2 **Copy2(TYPE2 **,TYPE2 **);
TYPE2 **Fill2(TYPE2 **,TYPE2);
void DiffusionDB(double **,double,struct point**[]);
int InitialPlaneSet(TYPE2**,int,TYPE2,...);

/*********************************************shift*/
TYPE2 **Boundaries2(TYPE2**);
double **BoundariesDB(double **,double);
/*****************************************neighbors*/
struct point **NeighSet(int);
void NeighFree(struct point**);
TYPE2 ***NeighSetP(TYPE2**,int);
void NeighFreeP(TYPE2***);

/******************************************margolus*/
MARGOLUS ***MargolusNeigh(void);
void MargolusFree(MARGOLUS***);
void MargolusDiffusion(TYPE2**,MARGOLUS***,int);

/*********************************************noise*/
void PerfectMix(TYPE2**); /* Shake() in CASH */
