/*
  CASH2-student

  CASH is "Cellular Automata Simulated Hardware" created by R de
  Boer (http://theory.bio.uu.nl/rdb/). CASH2 was created as an
  extension of CASH by N Takeuchi. CASH2-student is an easy to
  use version of CASH2 prepared for the use in the course
  Bioinformatic Processes given by P Hogeweg at Utrecht
  University. Authors hold copy right. No warranty.
*/

/* Version
   0.1 change for-loop to duff's device in Asynchronous().
   0.0 the first.
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
#include <stdarg.h>

#include "cash2003.h"
#include "cash2.h"
#include "mersenne.h"

/* system parameters */
#define BUFSIZE 250

/* PlaneNei structure will keep a track of planes and their
   neigobrs: we do not know how many planes user makes and which
   plane user is referring when user want a neighbor. */
struct PlaneNei {
  TYPE2 **plane;
  TYPE2 ***neighbor[9];
};

int MaxTime=INT_MAX;
int Time=0;
int nplane=1;
unsigned long ulseedG=56;

int argc_g;
char **argv_g;

int display=1;
static int movie_g=0;

int margin=0;
static MARGOLUS ***margolusG;
static struct updorder *updorderG;
static struct point **neighborG[9];
static struct PlaneNei *planeneiG;

extern int nrow,ncol,scale,boundary;
extern TYPE2 boundaryvalue2;

static TYPE **viewG;
static int poseG=0;
static int stopG=0;

/* function prototype of cash2-student.c */
void InitXmgrace(void);
void InitialMain(void);
int Display(TYPE2**,...);
int InitGradationColor(int,int,int);
int UpdateGradationColor(TYPE2 **,double,double,int);
void MDiffusion(TYPE2**);
void DiffusionFVAL(TYPE2**,double,int);
void Plot(int,...);
void PlotArray(double []);
void PlotXY(double x,double y);
void Synchronous(int, ...);
void Asynchronous(void);

int GetNeighbor(TYPE2**,int,int,int);
int RandomMoore8(TYPE2**,int,int);
int RandomMoore9(TYPE2**,int,int);
int RandomNeumann4(TYPE2**,int,int);
int RandomNeumann5(TYPE2**,int,int);

int CountMoore8(TYPE2 **,int,int,int);
int CountMoore9(TYPE2 **,int,int,int);
int CountNeumann4(TYPE2 **,int,int,int);
int CountNeumann5(TYPE2 **,int,int,int);

TYPE2 GetNeighborS(TYPE2**,int,int,int);
TYPE2 RandomMooreS8(TYPE2**,int,int);
TYPE2 RandomMooreS9(TYPE2**,int,int);
TYPE2 RandomNeumannS4(TYPE2**,int,int);
TYPE2 RandomNeumannS5(TYPE2**,int,int);

TYPE2 CountMooreS8(TYPE2**,int,int,int);
TYPE2 CountMooreS9(TYPE2**,int,int,int);
TYPE2 CountNeumannS4(TYPE2**,int,int,int);
TYPE2 CountNeumannS5(TYPE2**,int,int,int);

int SumMoore8(TYPE2**,int,int);
int SumMoore9(TYPE2**,int,int);
int SumNeumann4(TYPE2**,int,int);
int SumNeumann5(TYPE2**,int,int);

TYPE2 SumMooreS8(TYPE2**,int,int);
TYPE2 SumMooreS9(TYPE2**,int,int);
TYPE2 SumNeumannS4(TYPE2**,int,int);
TYPE2 SumNeumannS5(TYPE2**,int,int);

void GetNeighborC(TYPE2**,int,int,int,int*,int*);
void RandomMooreC8(TYPE2**,int,int,int*,int*);
void RandomMooreC9(TYPE2**,int,int,int*,int*);
void RandomNeumannC4(TYPE2**,int,int,int*,int*);
void RandomNeumannC5(TYPE2**,int,int,int*,int*);

TYPE2* GetNeighborP(TYPE2**,int,int,int);
TYPE2* RandomMooreP8(TYPE2**,int,int);
TYPE2* RandomMooreP9(TYPE2**,int,int);
TYPE2* RandomNeumannP4(TYPE2**,int,int);
TYPE2* RandomNeumannP5(TYPE2**,int,int);

int gotMouse(void);
void MakePlane(TYPE2***,...);
void SavePlane(char *,TYPE2 **,...);
void ReadSavedData(char *,int,TYPE2**,...);
void SpaceTimePlot(TYPE2**,TYPE2**);

int InitialSet(TYPE2**,int,int,...);
int InitialSetS(TYPE2**,int,TYPE2,...);


/* function prototype of model.c (made by students) */
void Initial(void);
void InitialPlane(void);
void NextState(int,int);
void Update(void);

int main (int argc, char *argv[])
{

  /* Setting default values. */
  nrow = 100;
  ncol = 100;
  scale = 2;
  boundary = WRAP;
  boundaryvalue2 = (TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.};

  /* make arguments available */
  argc_g = argc;
  argv_g = argv;

  /* initilize */
  InitialMain();

  for(Time=0;Time<=MaxTime;Time++){
    Update();
    
    if(display==1)
      if(Mouse()==1||poseG==1)gotMouse();

    if(stopG==1)
      break;
  }

  if(display==1)
    CloseDisplay();
  if(GraceIsOpen()&&display==0)
    GraceClose();
  if(movie_g==1)
    ClosePNG();
  return 0;
}

void InitXmgrace(void)
{
  int i;
  if(display==1){
    if (GraceOpen(16384) == -1) {
      fprintf(stderr,"I cannot run Xmgrace.\n");
      exit(-1);
    }
  }else{
    if (GraceOpenVA("gracebat",16384,"-noask","-nosafe",NULL) == -1) {
      fprintf(stderr,"I cannot run Grace. \n");
      exit(-1);
    }
  }
  /* g0 is a graph for population. */
  
  GracePrintf("g0 on");
  GracePrintf("with g0");

  GracePrintf("focus g0");

  for(i=0;i<15;i++){
    GracePrintf("s%d on",i);
    GracePrintf("s%d color %d",i,i+1);
  }
}

void InitialMain(void)
{
  /* set the value of nrow,ncol,boundary,scale,nplane,ulseed */
  Initial();   /* defined by student */

  init_genrand(ulseedG);
  
  ColorRGB(0,0,0,0);
  ColorRGB(1,255,255,255);
  ColorRGB(2,255,0,0);
  ColorRGB(3,0,255,0);
  ColorRGB(4,0,0,255);
  ColorRGB(5,255,255,0);
  ColorRGB(6,188, 143, 143);
  ColorRGB(7,220, 220, 220);
  ColorRGB(8,148, 0, 211);
  ColorRGB(9,0, 255, 255);
  ColorRGB(10,255, 0, 255);
  ColorRGB(11,255, 165, 0);
  ColorRGB(12,114, 33, 188);
  ColorRGB(13,103, 7, 72);
  ColorRGB(14,64, 224, 208);
  ColorRGB(15,0, 139, 0);

  /*  ColorRGB(1,254,254,254);
      ColorRGB(2,254,0,0);
      ColorRGB(3,254,254,0);
      ColorRGB(4,0,0,254);
      ColorRGB(5,255,0,255);
      ColorRGB(6,0,255,255);
      ColorRGB(7,0,255,0);
      ColorRGB(8,85,85,85);
      ColorRGB(9,198,113,113);
      ColorRGB(10,113,198,113);
      ColorRGB(11,142,142,56);
      ColorRGB(12,113,113,198);
      ColorRGB(13,142,56,142);
      ColorRGB(14,56,142,142);
      ColorRGB(15,170,170,170); */
  ColorRGB(16,76,113,158);
  ColorRGB(17,73,91,91);
  ColorRGB(18,35,43,43);
  ColorRGB(19,164,173,173);
  ColorRGB(20,0,55,91);
  ColorRGB(21,115,145,165);
  ColorRGB(22,255,191,0);
  ColorRGB(23,0,24,40);
  ColorRGB(24,153,145,58);
  ColorRGB(25,193,193,193);
  ColorRGB(26,39,39,39);
  ColorRGB(27,232,232,232);
  ColorRGB(28,112,112,112);
  ColorRGB(29,190,190,190);
  ColorRGB(30,185,142,142);
  ColorRGB(31,114,156,156);
  ColorRGB(32,86,128,171);
  ColorRGB(33,226,208,208);
  ColorRGB(34,102,78,78);
  ColorRGB(35,213,213,213);
  ColorRGB(36,153,153,153);
  ColorRGB(37,0,0,63);
  ColorRGB(38,214,214,214);
  ColorRGB(39,173,173,173);
  ColorRGB(40,204,204,204);
  ColorRGB(41,115,115,115);
  ColorRGB(42,76,76,76);
  ColorRGB(43,128,128,128);
  ColorRGB(44,239,239,239);
  ColorRGB(45,141,227,141);
  ColorRGB(46,155,155,155);
  ColorRGB(47,134,227,227);
  ColorRGB(48,156,156,156);
  ColorRGB(49,191,191,191);
  ColorRGB(50,131,131,131);
  ColorRGB(51,53,53,53);
  ColorRGB(52,176,0,0);
  ColorRGB(53,33,33,33);
  ColorRGB(54,8,91,8);
  ColorRGB(55,187,187,187);
  ColorRGB(56,93,93,93);
  ColorRGB(57,107,107,107);
  ColorRGB(58,7,7,7);
  ColorRGB(59,22,22,22);
  ColorRGB(60,45,45,45);
  ColorRGB(61,73,73,73);
  ColorRGB(62,146,146,146);
  ColorRGB(63,125,125,125);
  ColorRGB(64,0,0,100);
  ColorRGB(65,0,91,0);
  ColorRGB(66,105,32,172);
  ColorRGB(67,178,192,220);
  ColorRGB(68,198,213,226);
  ColorRGB(69,139,153,181);
  ColorRGB(70,55,5,29);
  ColorRGB(71,142,130,154);
  ColorRGB(72,255,0,102);
  ColorRGB(73,255,56,131);
  ColorRGB(74,255,161,181);
  ColorRGB(75,255,153,0);
  ColorRGB(76,168,97,12);
  ColorRGB(77,135,61,56);
  ColorRGB(78,109,211,57);
  ColorRGB(79,171,77,0);
  ColorRGB(80,5,122,11);
  ColorRGB(81,204,255,102);
  ColorRGB(82,255,102,204);
  ColorRGB(83,25,39,133);
  ColorRGB(84,140,255,216);
  ColorRGB(85,173,217,121);
  ColorRGB(86,99,61,255);
  ColorRGB(87,106,169,216);
  ColorRGB(88,198,101,255);
  ColorRGB(89,204,102,102);
  ColorRGB(90,249,198,100);
  ColorRGB(91,46,140,104);
  ColorRGB(92,53,70,30);
  ColorRGB(93,166,36,161);
  ColorRGB(94,57,70,175);
  ColorRGB(95,162,0,255);
  ColorRGB(96,125,168,175);
  ColorRGB(97,204,153,255);
  ColorRGB(98,255,204,153);
  ColorRGB(99,177,177,21);
  ColorRGB(100,185,109,173);
  ColorRGB(101,194,206,255);
  ColorRGB(102,155,46,255);
  ColorRGB(103,236,9,73);
  ColorRGB(104,255,243,82);
  ColorRGB(105,255,166,35);
  ColorRGB(106,255,240,98);
  ColorRGB(107,148,144,84);
  ColorRGB(108,173,236,192);
  ColorRGB(109,40,167,167);
  ColorRGB(110,102,200,143);
  ColorRGB(111,255,85,71);
  ColorRGB(112,67,255,63);
  ColorRGB(113,100,108,5);
  ColorRGB(114,99,165,179);
  ColorRGB(115,189,50,90);
  ColorRGB(116,25,26,61);
  ColorRGB(117,117,107,208);
  ColorRGB(118,211,255,216);
  ColorRGB(119,82,222,197);
  ColorRGB(120,204,204,51);
  ColorRGB(121,255,51,51);
  ColorRGB(122,21,61,193);
  ColorRGB(123,218,48,185);
  ColorRGB(124,27,87,111);
  ColorRGB(125,40,140,180);
  ColorRGB(126,93,88,203);
  ColorRGB(127,211,170,217);
  ColorRGB(128,114,229,88);
  ColorRGB(129,147,181,255);
  ColorRGB(130,44,89,162);
  ColorRGB(131,88,37,47);
  ColorRGB(132,39,67,115);
  ColorRGB(133,218,153,64);
  ColorRGB(134,155,94,170);
  ColorRGB(135,43,165,64);
  ColorRGB(136,44,69,216);
  ColorRGB(137,181,255,97);
  ColorRGB(138,27,51,167);
  ColorRGB(139,146,137,226);
  ColorRGB(140,86,129,185);
  ColorRGB(141,153,90,255);
  ColorRGB(142,94,146,183);
  ColorRGB(143,190,81,44);
  ColorRGB(144,149,86,66);
  ColorRGB(145,111,55,188);
  ColorRGB(146,207,106,193);
  ColorRGB(147,39,8,154);
  ColorRGB(148,76,58,255);
  ColorRGB(149,91,89,165);
  ColorRGB(150,225,225,225);
  ColorRGB(151,159,40,108);
  ColorRGB(152,235,255,215);
  ColorRGB(153,43,26,106);
  ColorRGB(154,110,89,179);
  ColorRGB(155,255,221,8);
  ColorRGB(156,235,186,78);
  ColorRGB(157,57,47,149);
  ColorRGB(158,245,74,152);
  ColorRGB(159,217,15,101);
  ColorRGB(160,159,155,27);
  ColorRGB(161,94,67,39);
  ColorRGB(162,157,127,72);
  ColorRGB(163,197,54,120);
  ColorRGB(164,193,69,69);
  ColorRGB(165,163,52,55);
  ColorRGB(166,175,38,97);
  ColorRGB(167,166,209,255);
  ColorRGB(168,48,31,121);
  ColorRGB(169,134,158,74);
  ColorRGB(170,60,52,138);
  ColorRGB(171,186,93,0);
  ColorRGB(172,255,0,13);
  ColorRGB(173,154,221,122);
  ColorRGB(174,43,83,145);
  ColorRGB(175,41,80,117);
  ColorRGB(176,57,54,149);
  ColorRGB(177,137,165,179);
  ColorRGB(178,162,181,205);
  ColorRGB(179,175,175,175);
  ColorRGB(180,159,159,159);
  ColorRGB(181,145,0,0);
  ColorRGB(182,91,91,91);
  ColorRGB(183,221,126,126);
  ColorRGB(184,185,185,185);
  ColorRGB(185,247,187,62);
  ColorRGB(186,231,186,129);
  ColorRGB(187,90,31,144);
  ColorRGB(188,230,88,216);
  ColorRGB(189,0,0,79);
  ColorRGB(190,210,180,140);
  ColorRGB(191,139,91,122);
  ColorRGB(192,179,179,179);
  ColorRGB(193,139,119,101);
  ColorRGB(194,238,130,238);
  ColorRGB(195,217,217,217);
  ColorRGB(196,114,119,133);
  ColorRGB(197,77,77,77);
  ColorRGB(198,51,255,255);
  ColorRGB(199,0,221,0);
  ColorRGB(200,0,0,221);
  ColorRGB(201,153,0,255);
  ColorRGB(202,255,255,51);
  ColorRGB(203,255,0,153);
  ColorRGB(204,248,105,166);
  ColorRGB(205,20,41,248);
  ColorRGB(206,177,41,248);
  ColorRGB(207,145,145,145);
  ColorRGB(208,47,47,47);
  ColorRGB(209,37,37,37);
  ColorRGB(210,192,192,192);
  ColorRGB(211,224,224,224);
  ColorRGB(212,139,106,106);
  ColorRGB(213,172,172,172);
  ColorRGB(214,153,0,0);
  ColorRGB(215,0,102,102);
  ColorRGB(216,0,204,255);
  ColorRGB(217,102,0,102);
  ColorRGB(218,102,102,0);
  ColorRGB(219,102,102,204);
  ColorRGB(220,102,204,102);
  ColorRGB(221,204,255,0);
  ColorRGB(222,255,0,204);
  ColorRGB(223,51,51,51);
  ColorRGB(224,51,153,153);
  ColorRGB(225,153,51,153);
  ColorRGB(226,153,153,51);
  ColorRGB(227,153,255,204);
  ColorRGB(228,0,51,204);
  ColorRGB(229,0,153,204);
  ColorRGB(230,0,204,51);
  ColorRGB(231,51,0,204);
  ColorRGB(232,51,51,153);
  ColorRGB(233,51,51,255);
  ColorRGB(234,51,153,51);
  ColorRGB(235,51,153,255);
  ColorRGB(236,51,204,0);
  ColorRGB(237,51,204,204);
  ColorRGB(238,51,255,51);
  ColorRGB(239,102,255,153);
  ColorRGB(240,153,51,51);
  ColorRGB(241,153,102,255);
  ColorRGB(242,153,204,0);
  ColorRGB(243,153,255,51);
  ColorRGB(244,204,0,51);
  ColorRGB(245,204,0,153);
  ColorRGB(246,204,51,0);
  ColorRGB(247,204,51,204);
  ColorRGB(248,255,51,153);
  ColorRGB(249,255,51,255);
  ColorRGB(250,255,153,102);
  ColorRGB(251,0,0,51);
  ColorRGB(252,0,0,102);
  ColorRGB(253,0,0,204);
  ColorRGB(254,0,51,0);
  //  ColorRGB(255,0,51,51);


  /* Create planes by New2() and create neighborhood planes by
     NeighSet(), NeighSetP() and MargolusNeigh(). Also initilize
     them */
  InitialPlane();   /* this must be defined by student */

  /* Create a plane for display */
  viewG = New();

  /* OpenDisplay() must be after the definition of the color
     table, thus after InitialPlane() */
  if(display==1)
    OpenDisplay("CASH2-student",nrow+2*margin,nplane*ncol+(1+nplane)*margin);
}

void Synchronous(int npl,...)
{
  int i,j,k,nr,nc;
  va_list ap;
  static TYPE2 ***prev=NULL; /* they will get planes of user */
  static TYPE2 ***next=NULL; /* they will be buffer to make next
				state, which will be copied into
				the user's planes in the end */
  static TYPE2 *copy=NULL;


  if(npl == 0){
    fprintf(stderr,"Synchronous() npl==0, that's wrong. See manual.\n");
    exit(-1);
  }

  /* If the first call */
  if(next==NULL){
    prev = (TYPE2***) malloc(sizeof(TYPE2**)*npl);

    next = (TYPE2***) malloc(sizeof(TYPE2**)*npl);
    for(i=0;i<npl;i++)
      *(next+i) = New2();

    copy = (TYPE2*) malloc(sizeof(TYPE2)*npl);
  }

  /* get arguments */
  va_start(ap,npl);
  for(i=0;i<npl;i++){
    *(prev+i) = va_arg(ap,TYPE2**);
  }

  nr = nrow;
  nc = ncol;

  for(i=1;i<=nr;i++)
    for(j=1;j<=nc;j++){
      /* take a copy of previous states at position [i,j]*/
      for(k=0;k<npl;k++)
	*(copy+k) = (*(prev+k))[i][j];

      /* This will modify the state of the [i][j] cell in
	 "*((prev+k))[i][j]". */
      NextState(i,j);

      for(k=0;k<npl;k++){
	/* Copy the new state of the [i][j] cells into the "next" */
	(*(next+k))[i][j] = (*(prev+k))[i][j];

	/* Revert the state of "prev" because this is synchronous
	   update. */
	(*(prev+k))[i][j] = *(copy+k);
      }
    }
  
  /* Then, update all the cells in user's plane, this makes
     CASH-student slower than CASH(or 2). */
  for(i=0;i<npl;i++)    
    Copy2(*(prev+i),*(next+i));
}

void Asynchronous(void)
{
  int i;
  int length;
  int irow,icol;
  int testVal=0;

  UpdOrdReset(&updorderG);
  
  length = nrow*ncol;

  i=(length+7)/8;
  switch(length%8)
    {
    case 0: do { irow=(updorderG+testVal)->row;
                 icol=(updorderG+testVal)->col;
                 NextState(irow,icol);
                 testVal++;
    case 7:      irow=(updorderG+testVal)->row;
                 icol=(updorderG+testVal)->col;
                 NextState(irow,icol);
                 testVal++;
    case 6:      irow=(updorderG+testVal)->row;
                 icol=(updorderG+testVal)->col;
                 NextState(irow,icol);
                 testVal++;
    case 5:      irow=(updorderG+testVal)->row;
                 icol=(updorderG+testVal)->col;
                 NextState(irow,icol);
                 testVal++;
    case 4:      irow=(updorderG+testVal)->row;
                 icol=(updorderG+testVal)->col;
                 NextState(irow,icol);
                 testVal++;
    case 3:      irow=(updorderG+testVal)->row;
                 icol=(updorderG+testVal)->col;
                 NextState(irow,icol);
                 testVal++;
    case 2:      irow=(updorderG+testVal)->row;
                 icol=(updorderG+testVal)->col;
                 NextState(irow,icol);
                 testVal++;
    case 1:      irow=(updorderG+testVal)->row;
                 icol=(updorderG+testVal)->col;
                 NextState(irow,icol);
                 testVal++;
               }while(--i>0);
    }
}

/* Usage: Display(Prey,Redator) */
int Display(TYPE2** Arg1,...)
{
  int i;
  int j,k,nr,nc;
  va_list ap;
  TYPE2** a;

  if(display==0)
    return 1;

  nr = nrow;
  nc = ncol;

  if(nplane > 1){
    va_start(ap,Arg1);
  }
  for(i=0;i<nplane;i++){
    if(i==0){
      a = Arg1;
    }else{
      a = va_arg(ap,TYPE2**);
    }
    for(j=1;j<=nr;j++)
      for(k=1;k<=nc;k++){
	viewG[j][k] = a[j][k].val;
      }

    PlaneDisplay(viewG,margin,margin+i*(ncol+margin),0);
  }

  return 0;
}

int DrawSlide(TYPE2 **Arg1,char *dirname)
{
  int j,k,nr,nc;
  static int init_flag=0;

  if(init_flag==0){
    OpenPNG(dirname,nrow,ncol);
    init_flag = 1;
    movie_g = 1;
  }
  
  nr = nrow;
  nc = ncol;
  
  for(j=1;j<=nr;j++)
    for(k=1;k<=nc;k++){
      viewG[j][k] = Arg1[j][k].val;
    }
  PlanePNG(viewG,0);
  
  return 0;

}

/* InitGradationColor() & UpdateGradationColor() are for coloring
   of floating number states of the cell */
static int minGradColor;
static int maxGradColor;
int InitGradationColor(int ncolor,int min,int max)
{
  int i;

  if(min<0 || max > 255 || min >= max){
    fprintf(stderr,"SetGradationColor: must be 0<min<max<255\n");
    exit(1);
  }

  /* # of color is 1 */
  if(ncolor<2){
    for(i=min;i<=max;i++){
      ColorRGB(i,(int)(255.*(double)(i-min)/(double)(max-min)),0,0);
    }
    minGradColor = min;
    maxGradColor = max;
  }
  /* # of color is 2 */
  else{
    if((max-min)%2!=0){
      fprintf(stderr,"SetGradationColor: if ncolor=2, max-min must be even\n");
      exit(1);
    }

    for(i=min;i<min+(max-min)/2;i++){
      ColorRGB(i,0,(int)(255.*(double)(min+(max-min)/2-i)/(double)((max-min)/2)),0);
    }
    ColorRGB(i,0,0,0);
    for(i=min+(max-min)/2+1;i<=max;i++){
      ColorRGB(i,(int)(255.*(double)(i-min-(max-min)/2)/(double)((max-min)/2)),0,0);
    }
    minGradColor = min;
    maxGradColor = max;
  }
  return 0;
}

int UpdateGradationColor(TYPE2 **plane,double min,double max,int fnum)
{
  int i,j;
  int nr=nrow,nc=ncol;
  double clen = (double)maxGradColor-minGradColor; /* color length */
  double cmin = (double)minGradColor;
  int cmini = minGradColor;
  int cmaxi = maxGradColor;
  double vlen = max - min; /* value length */
  double absMax,absMin,abs;

  if(max < min){
    fprintf(stderr,"UpdateGradationColor: must be max > min\n");
    exit(1);
  }

  absMax = (max > -1.*max) ? max : -1.*max;
  absMin = (min > -1.*min) ? min : -1.*min;
  abs = (absMax>absMin) ? absMax : absMin;
  abs = (abs>DBL_MIN) ? abs : DBL_MIN;

  if(vlen/abs < DBL_EPSILON){
    fprintf(stderr,"UpdateGradationColor: must be max != min\n");
    exit(1);
  }


  switch(fnum){
  case 1:
    for(i=1;i<=nr;i++)
      for(j=1;j<=nc;j++){
	/* +0.5 makes x.4 to x, x.6 to x+1*/
	plane[i][j].val=clen*(plane[i][j].fval-min)/vlen + cmin +0.5;
	if(plane[i][j].val<cmini)
	  plane[i][j].val=cmini;
	else if(plane[i][j].val>cmaxi)
	  plane[i][j].val=cmaxi;
      }
    break;
  case 2:
    for(i=1;i<=nr;i++)
      for(j=1;j<=nc;j++){
	plane[i][j].val=clen*(plane[i][j].fval2-min)/vlen + cmin +0.5;
	if(plane[i][j].val<cmini)
	  plane[i][j].val=cmini;
	else if(plane[i][j].val>cmaxi)
	  plane[i][j].val=cmaxi;
      }
    break;
  case 3:
    for(i=1;i<=nr;i++)
      for(j=1;j<=nc;j++){
	plane[i][j].val=clen*(plane[i][j].fval3-min)/vlen + cmin +0.5;
	if(plane[i][j].val<cmini)
	  plane[i][j].val=cmini;
	else if(plane[i][j].val>cmaxi)
	  plane[i][j].val=cmaxi;
      }
    break;
  case 4:
    for(i=1;i<=nr;i++)
      for(j=1;j<=nc;j++){
	plane[i][j].val=clen*(plane[i][j].fval4-min)/vlen + cmin +0.5;
	if(plane[i][j].val<cmini)
	  plane[i][j].val=cmini;
	else if(plane[i][j].val>cmaxi)
	  plane[i][j].val=cmaxi;
      }
    break;
  case 5:
    for(i=1;i<=nr;i++)
      for(j=1;j<=nc;j++){
	plane[i][j].val=clen*(plane[i][j].fval5-min)/vlen + cmin +0.5;
	if(plane[i][j].val<cmini)
	  plane[i][j].val=cmini;
	else if(plane[i][j].val>cmaxi)
	  plane[i][j].val=cmaxi;
      }
    break;
  default:
    break;
  }

    return 0;
}

void Plot(int npl,...)
{
  static int count;
  int i,j,k;
  int nc,nr;
  int pop[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  va_list ap;
  TYPE2 **world;
  static int initial_flag = 0;

  if(initial_flag == 0){
    InitXmgrace();
    initial_flag =1 ;
  }
  

  nr = nrow;
  nc = ncol;
  
  va_start(ap,npl);
  for(k=0;k<npl;k++){
    world = va_arg(ap,TYPE2**);
    
    for(i=1; i<=nr; i++)
      for(j=1; j<=nc; j++){
	if(world[i][j].val>0 && world[i][j].val<=15)
	  pop[world[i][j].val-1]++;
      }
  }
  
  for(i=0;i<15;i++)
    GracePrintf("g0.s%d point %d, %d",i,Time,pop[i]);


  if(count%10==0 && display==1)
    GracePrintf("redraw");
  
  if(count%100==0 && display==1){
    GracePrintf("focus g0");
    GracePrintf("autoscale");
    GracePrintf("redraw");
  }

  count++;
}

void PlotArray(double data[15])
{
  int i;
  static int count;
  static int initial_flag = 0;

  if(initial_flag == 0){
    InitXmgrace();
    initial_flag =1 ;
  }


  for(i=0;i<15;i++)
    GracePrintf("g0.s%d point %d, %f",i,Time,data[i]);

  if(count%10==0 && display==1)
    GracePrintf("redraw");
  
  if(count%100==0 && display==1){
    GracePrintf("focus g0");
    GracePrintf("autoscale");
    GracePrintf("redraw");
  }

  count++;
}

void SavePlot(char *fname)
{
  if(GraceIsOpen()){
    GraceFlush();
    GracePrintf("saveall \"%s\"",fname);
    GraceClose();
  }
}

/* phase_list will keep track of the phase of the planes: we
   don't know how many planes user have made and which plane user
   wants to do diffusion. */
struct phase_list {
  int phase;
  TYPE2 **plane;
};

void  MDiffusion(TYPE2** a)
{
  int i;
  int planeID;
  static int init_flag=0;
  static struct phase_list *list;
 
  /* initilization */
  if(init_flag==0){
    list = (struct phase_list *)malloc(sizeof(struct phase_list)*nplane);
    for(i=0;i<nplane;i++){
      (list+i)->phase = 0;
      (list+i)->plane = NULL;
    }
    init_flag = 1;
  }

  /* search "a_name" in the list */
  planeID = -1;
  i=0;
  while(i<nplane){
    if((list+i)->plane == a){
      planeID = i;
      break;
    }
    else if((list+i)->plane==NULL){
      /* insert this plane */
      (list+i)->plane = a;
      planeID = i;
      break;
    }
    i++;
  }

  /* weak error cheack */
  if(planeID == -1){
    fprintf(stderr,"MDiffusion: a bug. please report this\n");
    exit(1);
  }

  /* do Margolus diffusion */
  MargolusDiffusion(a,margolusG,(list+planeID)->phase);

  /* incliment the phase */
  (list+planeID)->phase = ((list+planeID)->phase+1)%2;
}

/* diffusion of continuous values. DiffusionFVAL() function is
   specilized version of DiffusionDBP() for the use in
   cash2-student version */
void DiffusionFVAL(TYPE2** plane,double diff_rate,int nval)
{
  int iplane;
  static double **influx=NULL;
  int i, j, nc, nr;

  /* if first call */
  if(influx==NULL){
    influx = NewDB();
  }

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  nr = nrow;
  nc = ncol;

  switch(nval){
  case 1: /* if fval */
    for (i=1;i<=nr;i++)
      for (j=1;j<=nc;j++){
	influx[i][j] = (planeneiG[iplane].neighbor[1][i][j])->fval+(planeneiG[iplane].neighbor[2][i][j])->fval+(planeneiG[iplane].neighbor[3][i][j])->fval+(planeneiG[iplane].neighbor[4][i][j])->fval;
      }
    for (i=1; i <= nr; i++)
      /*$dir force_vector*/
      for (j=1; j <= nc; j++){
	plane[i][j].fval = plane[i][j].fval + diff_rate*influx[i][j] - 4*diff_rate*plane[i][j].fval;
      }
    break;
  case 2: /* if fval2 */
    for (i=1;i<=nr;i++)
      for (j=1;j<=nc;j++){
	influx[i][j] = (planeneiG[iplane].neighbor[1][i][j])->fval2+(planeneiG[iplane].neighbor[2][i][j])->fval2+(planeneiG[iplane].neighbor[3][i][j])->fval2+(planeneiG[iplane].neighbor[4][i][j])->fval2;
      }
    for (i=1; i <= nr; i++)
      /*$dir force_vector*/
      for (j=1; j <= nc; j++){
	plane[i][j].fval2 = plane[i][j].fval2 + diff_rate*influx[i][j] - 4*diff_rate*plane[i][j].fval2;
      }
    break;
  case 3: /* if fval3 */
    for (i=1;i<=nr;i++)
      for (j=1;j<=nc;j++){
	influx[i][j] = (planeneiG[iplane].neighbor[1][i][j])->fval3+(planeneiG[iplane].neighbor[2][i][j])->fval3+(planeneiG[iplane].neighbor[3][i][j])->fval3+(planeneiG[iplane].neighbor[4][i][j])->fval3;
      }
    for (i=1; i <= nr; i++)
      /*$dir force_vector*/
      for (j=1; j <= nc; j++){
	plane[i][j].fval3 = plane[i][j].fval3 + diff_rate*influx[i][j] - 4*diff_rate*plane[i][j].fval3;
      }
    break;
  case 4: /* if fval4 */
    for (i=1;i<=nr;i++)
      for (j=1;j<=nc;j++){
	influx[i][j] = (planeneiG[iplane].neighbor[1][i][j])->fval4+(planeneiG[iplane].neighbor[2][i][j])->fval4+(planeneiG[iplane].neighbor[3][i][j])->fval4+(planeneiG[iplane].neighbor[4][i][j])->fval4;
      }
    for (i=1; i <= nr; i++)
      /*$dir force_vector*/
      for (j=1; j <= nc; j++){
	plane[i][j].fval4 = plane[i][j].fval4 + diff_rate*influx[i][j] - 4*diff_rate*plane[i][j].fval4;
      }
    break;
  case 5: /* if fval5 */
    for (i=1;i<=nr;i++)
      for (j=1;j<=nc;j++){
	influx[i][j] = (planeneiG[iplane].neighbor[1][i][j])->fval5+(planeneiG[iplane].neighbor[2][i][j])->fval5+(planeneiG[iplane].neighbor[3][i][j])->fval5+(planeneiG[iplane].neighbor[4][i][j])->fval5;
      }
    for (i=1; i <= nr; i++)
      /*$dir force_vector*/
      for (j=1; j <= nc; j++){
	plane[i][j].fval5 = plane[i][j].fval5 + diff_rate*influx[i][j] - 4*diff_rate*plane[i][j].fval5;
      }
    break;
  default:
    break;
  }
}

int GetNeighbor(TYPE2** a,int row,int col,int nei)
{
  int iplane;
  if(nei<0 || nei>8){
    fprintf(stderr,"GetNeighbor(): wrong value in the forth argument\n");
    exit(1);
  }

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  return((planeneiG[iplane].neighbor[nei][row][col])->val);
}

int RandomMoore8(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(1,8);
  return((planeneiG[iplane].neighbor[nei][row][col])->val);
}

int RandomMoore9(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(0,8);
  return((planeneiG[iplane].neighbor[nei][row][col])->val);
}

int RandomNeumann4(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(1,4);
  return((planeneiG[iplane].neighbor[nei][row][col])->val);
}

int RandomNeumann5(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(0,4);
  return((planeneiG[iplane].neighbor[nei][row][col])->val);
}

int CountMoore8(TYPE2 **plane,int aval,int row,int col)
{
  int iplane;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  return(((planeneiG[iplane].neighbor[1][row][col])->val==aval)+((planeneiG[iplane].neighbor[2][row][col])->val==aval)+((planeneiG[iplane].neighbor[3][row][col])->val==aval)+((planeneiG[iplane].neighbor[4][row][col])->val==aval)+((planeneiG[iplane].neighbor[5][row][col])->val==aval)+((planeneiG[iplane].neighbor[6][row][col])->val==aval)+((planeneiG[iplane].neighbor[7][row][col])->val==aval)+((planeneiG[iplane].neighbor[8][row][col])->val==aval));
}

int CountMoore9(TYPE2 **plane,int aval,int row,int col)
{
  int iplane;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  return(((planeneiG[iplane].neighbor[0][row][col])->val==aval)+((planeneiG[iplane].neighbor[1][row][col])->val==aval)+((planeneiG[iplane].neighbor[2][row][col])->val==aval)+((planeneiG[iplane].neighbor[3][row][col])->val==aval)+((planeneiG[iplane].neighbor[4][row][col])->val==aval)+((planeneiG[iplane].neighbor[5][row][col])->val==aval)+((planeneiG[iplane].neighbor[6][row][col])->val==aval)+((planeneiG[iplane].neighbor[7][row][col])->val==aval)+((planeneiG[iplane].neighbor[8][row][col])->val==aval));

}

int CountNeumann4(TYPE2 **plane,int aval,int row,int col)
{
  int iplane;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  return(((planeneiG[iplane].neighbor[1][row][col])->val==aval)+((planeneiG[iplane].neighbor[2][row][col])->val==aval)+((planeneiG[iplane].neighbor[3][row][col])->val==aval)+((planeneiG[iplane].neighbor[4][row][col])->val==aval));
}

int CountNeumann5(TYPE2 **plane,int aval,int row,int col)
{
  int iplane;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  return(((planeneiG[iplane].neighbor[0][row][col])->val==aval)+((planeneiG[iplane].neighbor[1][row][col])->val==aval)+((planeneiG[iplane].neighbor[2][row][col])->val==aval)+((planeneiG[iplane].neighbor[3][row][col])->val==aval)+((planeneiG[iplane].neighbor[4][row][col])->val==aval));
}

int SumMoore8(TYPE2**plane,int row,int col)
{
  int iplane;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  return(((planeneiG[iplane].neighbor[1][row][col])->val)+((planeneiG[iplane].neighbor[2][row][col])->val)+((planeneiG[iplane].neighbor[3][row][col])->val)+((planeneiG[iplane].neighbor[4][row][col])->val)+((planeneiG[iplane].neighbor[5][row][col])->val)+((planeneiG[iplane].neighbor[6][row][col])->val)+((planeneiG[iplane].neighbor[7][row][col])->val)+((planeneiG[iplane].neighbor[8][row][col])->val));
}

int SumMoore9(TYPE2**plane,int row,int col){
  int iplane;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  return(((planeneiG[iplane].neighbor[0][row][col])->val)+((planeneiG[iplane].neighbor[1][row][col])->val)+((planeneiG[iplane].neighbor[2][row][col])->val)+((planeneiG[iplane].neighbor[3][row][col])->val)+((planeneiG[iplane].neighbor[4][row][col])->val)+((planeneiG[iplane].neighbor[5][row][col])->val)+((planeneiG[iplane].neighbor[6][row][col])->val)+((planeneiG[iplane].neighbor[7][row][col])->val)+((planeneiG[iplane].neighbor[8][row][col])->val));
}

int SumNeumann4(TYPE2**plane,int row,int col)
{
  int iplane;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  return(((planeneiG[iplane].neighbor[1][row][col])->val)+((planeneiG[iplane].neighbor[2][row][col])->val)+((planeneiG[iplane].neighbor[3][row][col])->val)+((planeneiG[iplane].neighbor[4][row][col])->val));
}
int SumNeumann5(TYPE2**plane,int row,int col)
{
  int iplane;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  return(((planeneiG[iplane].neighbor[0][row][col])->val)+((planeneiG[iplane].neighbor[1][row][col])->val)+((planeneiG[iplane].neighbor[2][row][col])->val)+((planeneiG[iplane].neighbor[3][row][col])->val)+((planeneiG[iplane].neighbor[4][row][col])->val));
}

TYPE2 SumMooreS8(TYPE2** plane,int row,int col)
{
  int iplane;
  TYPE2 result;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  result.val = (((planeneiG[iplane].neighbor[1][row][col])->val)+((planeneiG[iplane].neighbor[2][row][col])->val)+((planeneiG[iplane].neighbor[3][row][col])->val)+((planeneiG[iplane].neighbor[4][row][col])->val)+((planeneiG[iplane].neighbor[5][row][col])->val)+((planeneiG[iplane].neighbor[6][row][col])->val)+((planeneiG[iplane].neighbor[7][row][col])->val)+((planeneiG[iplane].neighbor[8][row][col])->val));

  result.val2 = (((planeneiG[iplane].neighbor[1][row][col])->val2)+((planeneiG[iplane].neighbor[2][row][col])->val2)+((planeneiG[iplane].neighbor[3][row][col])->val2)+((planeneiG[iplane].neighbor[4][row][col])->val2)+((planeneiG[iplane].neighbor[5][row][col])->val2)+((planeneiG[iplane].neighbor[6][row][col])->val2)+((planeneiG[iplane].neighbor[7][row][col])->val2)+((planeneiG[iplane].neighbor[8][row][col])->val2));

  result.val3 = (((planeneiG[iplane].neighbor[1][row][col])->val3)+((planeneiG[iplane].neighbor[2][row][col])->val3)+((planeneiG[iplane].neighbor[3][row][col])->val3)+((planeneiG[iplane].neighbor[4][row][col])->val3)+((planeneiG[iplane].neighbor[5][row][col])->val3)+((planeneiG[iplane].neighbor[6][row][col])->val3)+((planeneiG[iplane].neighbor[7][row][col])->val3)+((planeneiG[iplane].neighbor[8][row][col])->val3));

  result.val4 = (((planeneiG[iplane].neighbor[1][row][col])->val4)+((planeneiG[iplane].neighbor[2][row][col])->val4)+((planeneiG[iplane].neighbor[3][row][col])->val4)+((planeneiG[iplane].neighbor[4][row][col])->val4)+((planeneiG[iplane].neighbor[5][row][col])->val4)+((planeneiG[iplane].neighbor[6][row][col])->val4)+((planeneiG[iplane].neighbor[7][row][col])->val4)+((planeneiG[iplane].neighbor[8][row][col])->val4));

  result.val5 = (((planeneiG[iplane].neighbor[1][row][col])->val5)+((planeneiG[iplane].neighbor[2][row][col])->val5)+((planeneiG[iplane].neighbor[3][row][col])->val5)+((planeneiG[iplane].neighbor[4][row][col])->val5)+((planeneiG[iplane].neighbor[5][row][col])->val5)+((planeneiG[iplane].neighbor[6][row][col])->val5)+((planeneiG[iplane].neighbor[7][row][col])->val5)+((planeneiG[iplane].neighbor[8][row][col])->val5));

  result.fval = (((planeneiG[iplane].neighbor[1][row][col])->fval)+((planeneiG[iplane].neighbor[2][row][col])->fval)+((planeneiG[iplane].neighbor[3][row][col])->fval)+((planeneiG[iplane].neighbor[4][row][col])->fval)+((planeneiG[iplane].neighbor[5][row][col])->fval)+((planeneiG[iplane].neighbor[6][row][col])->fval)+((planeneiG[iplane].neighbor[7][row][col])->fval)+((planeneiG[iplane].neighbor[8][row][col])->fval));

  result.fval2 = (((planeneiG[iplane].neighbor[1][row][col])->fval2)+((planeneiG[iplane].neighbor[2][row][col])->fval2)+((planeneiG[iplane].neighbor[3][row][col])->fval2)+((planeneiG[iplane].neighbor[4][row][col])->fval2)+((planeneiG[iplane].neighbor[5][row][col])->fval2)+((planeneiG[iplane].neighbor[6][row][col])->fval2)+((planeneiG[iplane].neighbor[7][row][col])->fval2)+((planeneiG[iplane].neighbor[8][row][col])->fval2));

  result.fval3 = (((planeneiG[iplane].neighbor[1][row][col])->fval3)+((planeneiG[iplane].neighbor[2][row][col])->fval3)+((planeneiG[iplane].neighbor[3][row][col])->fval3)+((planeneiG[iplane].neighbor[4][row][col])->fval3)+((planeneiG[iplane].neighbor[5][row][col])->fval3)+((planeneiG[iplane].neighbor[6][row][col])->fval3)+((planeneiG[iplane].neighbor[7][row][col])->fval3)+((planeneiG[iplane].neighbor[8][row][col])->fval3));

  result.fval4 = (((planeneiG[iplane].neighbor[1][row][col])->fval4)+((planeneiG[iplane].neighbor[2][row][col])->fval4)+((planeneiG[iplane].neighbor[3][row][col])->fval4)+((planeneiG[iplane].neighbor[4][row][col])->fval4)+((planeneiG[iplane].neighbor[5][row][col])->fval4)+((planeneiG[iplane].neighbor[6][row][col])->fval4)+((planeneiG[iplane].neighbor[7][row][col])->fval4)+((planeneiG[iplane].neighbor[8][row][col])->fval4));

  result.fval5 = (((planeneiG[iplane].neighbor[1][row][col])->fval5)+((planeneiG[iplane].neighbor[2][row][col])->fval5)+((planeneiG[iplane].neighbor[3][row][col])->fval5)+((planeneiG[iplane].neighbor[4][row][col])->fval5)+((planeneiG[iplane].neighbor[5][row][col])->fval5)+((planeneiG[iplane].neighbor[6][row][col])->fval5)+((planeneiG[iplane].neighbor[7][row][col])->fval5)+((planeneiG[iplane].neighbor[8][row][col])->fval5));

  return result;
}
TYPE2 SumMooreS9(TYPE2**plane,int row,int col)
{
  int iplane;
  TYPE2 result;
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }
  result.val = (((planeneiG[iplane].neighbor[0][row][col])->val)+((planeneiG[iplane].neighbor[1][row][col])->val)+((planeneiG[iplane].neighbor[2][row][col])->val)+((planeneiG[iplane].neighbor[3][row][col])->val)+((planeneiG[iplane].neighbor[4][row][col])->val)+((planeneiG[iplane].neighbor[5][row][col])->val)+((planeneiG[iplane].neighbor[6][row][col])->val)+((planeneiG[iplane].neighbor[7][row][col])->val)+((planeneiG[iplane].neighbor[8][row][col])->val));

  result.val2 = (((planeneiG[iplane].neighbor[0][row][col])->val2)+((planeneiG[iplane].neighbor[1][row][col])->val2)+((planeneiG[iplane].neighbor[2][row][col])->val2)+((planeneiG[iplane].neighbor[3][row][col])->val2)+((planeneiG[iplane].neighbor[4][row][col])->val2)+((planeneiG[iplane].neighbor[5][row][col])->val2)+((planeneiG[iplane].neighbor[6][row][col])->val2)+((planeneiG[iplane].neighbor[7][row][col])->val2)+((planeneiG[iplane].neighbor[8][row][col])->val2));

  result.val3 = (((planeneiG[iplane].neighbor[0][row][col])->val3)+((planeneiG[iplane].neighbor[1][row][col])->val3)+((planeneiG[iplane].neighbor[2][row][col])->val3)+((planeneiG[iplane].neighbor[3][row][col])->val3)+((planeneiG[iplane].neighbor[4][row][col])->val3)+((planeneiG[iplane].neighbor[5][row][col])->val3)+((planeneiG[iplane].neighbor[6][row][col])->val3)+((planeneiG[iplane].neighbor[7][row][col])->val3)+((planeneiG[iplane].neighbor[8][row][col])->val3));

  result.val4 = (((planeneiG[iplane].neighbor[0][row][col])->val4)+((planeneiG[iplane].neighbor[1][row][col])->val4)+((planeneiG[iplane].neighbor[2][row][col])->val4)+((planeneiG[iplane].neighbor[3][row][col])->val4)+((planeneiG[iplane].neighbor[4][row][col])->val4)+((planeneiG[iplane].neighbor[5][row][col])->val4)+((planeneiG[iplane].neighbor[6][row][col])->val4)+((planeneiG[iplane].neighbor[7][row][col])->val4)+((planeneiG[iplane].neighbor[8][row][col])->val4));

  result.val5 = (((planeneiG[iplane].neighbor[0][row][col])->val5)+((planeneiG[iplane].neighbor[1][row][col])->val5)+((planeneiG[iplane].neighbor[2][row][col])->val5)+((planeneiG[iplane].neighbor[3][row][col])->val5)+((planeneiG[iplane].neighbor[4][row][col])->val5)+((planeneiG[iplane].neighbor[5][row][col])->val5)+((planeneiG[iplane].neighbor[6][row][col])->val5)+((planeneiG[iplane].neighbor[7][row][col])->val5)+((planeneiG[iplane].neighbor[8][row][col])->val5));

  result.fval = (((planeneiG[iplane].neighbor[0][row][col])->fval)+((planeneiG[iplane].neighbor[1][row][col])->fval)+((planeneiG[iplane].neighbor[2][row][col])->fval)+((planeneiG[iplane].neighbor[3][row][col])->fval)+((planeneiG[iplane].neighbor[4][row][col])->fval)+((planeneiG[iplane].neighbor[5][row][col])->fval)+((planeneiG[iplane].neighbor[6][row][col])->fval)+((planeneiG[iplane].neighbor[7][row][col])->fval)+((planeneiG[iplane].neighbor[8][row][col])->fval));

  result.fval2 = (((planeneiG[iplane].neighbor[0][row][col])->fval2)+((planeneiG[iplane].neighbor[1][row][col])->fval2)+((planeneiG[iplane].neighbor[2][row][col])->fval2)+((planeneiG[iplane].neighbor[3][row][col])->fval2)+((planeneiG[iplane].neighbor[4][row][col])->fval2)+((planeneiG[iplane].neighbor[5][row][col])->fval2)+((planeneiG[iplane].neighbor[6][row][col])->fval2)+((planeneiG[iplane].neighbor[7][row][col])->fval2)+((planeneiG[iplane].neighbor[8][row][col])->fval2));

  result.fval3 = (((planeneiG[iplane].neighbor[0][row][col])->fval3)+((planeneiG[iplane].neighbor[1][row][col])->fval3)+((planeneiG[iplane].neighbor[2][row][col])->fval3)+((planeneiG[iplane].neighbor[3][row][col])->fval3)+((planeneiG[iplane].neighbor[4][row][col])->fval3)+((planeneiG[iplane].neighbor[5][row][col])->fval3)+((planeneiG[iplane].neighbor[6][row][col])->fval3)+((planeneiG[iplane].neighbor[7][row][col])->fval3)+((planeneiG[iplane].neighbor[8][row][col])->fval3));

  result.fval4 = (((planeneiG[iplane].neighbor[0][row][col])->fval4)+((planeneiG[iplane].neighbor[1][row][col])->fval4)+((planeneiG[iplane].neighbor[2][row][col])->fval4)+((planeneiG[iplane].neighbor[3][row][col])->fval4)+((planeneiG[iplane].neighbor[4][row][col])->fval4)+((planeneiG[iplane].neighbor[5][row][col])->fval4)+((planeneiG[iplane].neighbor[6][row][col])->fval4)+((planeneiG[iplane].neighbor[7][row][col])->fval4)+((planeneiG[iplane].neighbor[8][row][col])->fval4));

  result.fval5 = (((planeneiG[iplane].neighbor[0][row][col])->fval5)+((planeneiG[iplane].neighbor[1][row][col])->fval5)+((planeneiG[iplane].neighbor[2][row][col])->fval5)+((planeneiG[iplane].neighbor[3][row][col])->fval5)+((planeneiG[iplane].neighbor[4][row][col])->fval5)+((planeneiG[iplane].neighbor[5][row][col])->fval5)+((planeneiG[iplane].neighbor[6][row][col])->fval5)+((planeneiG[iplane].neighbor[7][row][col])->fval5)+((planeneiG[iplane].neighbor[8][row][col])->fval5));

  return result;
}
TYPE2 SumNeumannS4(TYPE2**plane,int row,int col)
{
  int iplane;
  TYPE2 result;
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }
  result.val = (((planeneiG[iplane].neighbor[1][row][col])->val)+((planeneiG[iplane].neighbor[2][row][col])->val)+((planeneiG[iplane].neighbor[3][row][col])->val)+((planeneiG[iplane].neighbor[4][row][col])->val));

  result.val2 = (((planeneiG[iplane].neighbor[1][row][col])->val2)+((planeneiG[iplane].neighbor[2][row][col])->val2)+((planeneiG[iplane].neighbor[3][row][col])->val2)+((planeneiG[iplane].neighbor[4][row][col])->val2));

  result.val3 = (((planeneiG[iplane].neighbor[1][row][col])->val3)+((planeneiG[iplane].neighbor[2][row][col])->val3)+((planeneiG[iplane].neighbor[3][row][col])->val3)+((planeneiG[iplane].neighbor[4][row][col])->val3));

  result.val4 = (((planeneiG[iplane].neighbor[1][row][col])->val4)+((planeneiG[iplane].neighbor[2][row][col])->val4)+((planeneiG[iplane].neighbor[3][row][col])->val4)+((planeneiG[iplane].neighbor[4][row][col])->val4));

  result.val5 = (((planeneiG[iplane].neighbor[1][row][col])->val5)+((planeneiG[iplane].neighbor[2][row][col])->val5)+((planeneiG[iplane].neighbor[3][row][col])->val5)+((planeneiG[iplane].neighbor[4][row][col])->val5));

  result.fval = (((planeneiG[iplane].neighbor[1][row][col])->fval)+((planeneiG[iplane].neighbor[2][row][col])->fval)+((planeneiG[iplane].neighbor[3][row][col])->fval)+((planeneiG[iplane].neighbor[4][row][col])->fval));

  result.fval2 = (((planeneiG[iplane].neighbor[1][row][col])->fval2)+((planeneiG[iplane].neighbor[2][row][col])->fval2)+((planeneiG[iplane].neighbor[3][row][col])->fval2)+((planeneiG[iplane].neighbor[4][row][col])->fval2));

  result.fval3 = (((planeneiG[iplane].neighbor[1][row][col])->fval3)+((planeneiG[iplane].neighbor[2][row][col])->fval3)+((planeneiG[iplane].neighbor[3][row][col])->fval3)+((planeneiG[iplane].neighbor[4][row][col])->fval3));

  result.fval4 = (((planeneiG[iplane].neighbor[1][row][col])->fval4)+((planeneiG[iplane].neighbor[2][row][col])->fval4)+((planeneiG[iplane].neighbor[3][row][col])->fval4)+((planeneiG[iplane].neighbor[4][row][col])->fval4));

  result.fval5 = (((planeneiG[iplane].neighbor[1][row][col])->fval5)+((planeneiG[iplane].neighbor[2][row][col])->fval5)+((planeneiG[iplane].neighbor[3][row][col])->fval5)+((planeneiG[iplane].neighbor[4][row][col])->fval5));

  return result;
}
TYPE2 SumNeumannS5(TYPE2**plane,int row,int col)
{
  int iplane;
  TYPE2 result;
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }
  result.val = (((planeneiG[iplane].neighbor[0][row][col])->val)+((planeneiG[iplane].neighbor[1][row][col])->val)+((planeneiG[iplane].neighbor[2][row][col])->val)+((planeneiG[iplane].neighbor[3][row][col])->val)+((planeneiG[iplane].neighbor[4][row][col])->val));

  result.val2 = (((planeneiG[iplane].neighbor[0][row][col])->val2)+((planeneiG[iplane].neighbor[1][row][col])->val2)+((planeneiG[iplane].neighbor[2][row][col])->val2)+((planeneiG[iplane].neighbor[3][row][col])->val2)+((planeneiG[iplane].neighbor[4][row][col])->val2));

  result.val3 = (((planeneiG[iplane].neighbor[0][row][col])->val3)+((planeneiG[iplane].neighbor[1][row][col])->val3)+((planeneiG[iplane].neighbor[2][row][col])->val3)+((planeneiG[iplane].neighbor[3][row][col])->val3)+((planeneiG[iplane].neighbor[4][row][col])->val3));

  result.val4 = (((planeneiG[iplane].neighbor[0][row][col])->val4)+((planeneiG[iplane].neighbor[1][row][col])->val4)+((planeneiG[iplane].neighbor[2][row][col])->val4)+((planeneiG[iplane].neighbor[3][row][col])->val4)+((planeneiG[iplane].neighbor[4][row][col])->val4));

  result.val5 = (((planeneiG[iplane].neighbor[0][row][col])->val5)+((planeneiG[iplane].neighbor[1][row][col])->val5)+((planeneiG[iplane].neighbor[2][row][col])->val5)+((planeneiG[iplane].neighbor[3][row][col])->val5)+((planeneiG[iplane].neighbor[4][row][col])->val5));

  result.fval = (((planeneiG[iplane].neighbor[0][row][col])->fval)+((planeneiG[iplane].neighbor[1][row][col])->fval)+((planeneiG[iplane].neighbor[2][row][col])->fval)+((planeneiG[iplane].neighbor[3][row][col])->fval)+((planeneiG[iplane].neighbor[4][row][col])->fval));

  result.fval2 = (((planeneiG[iplane].neighbor[0][row][col])->fval2)+((planeneiG[iplane].neighbor[1][row][col])->fval2)+((planeneiG[iplane].neighbor[2][row][col])->fval2)+((planeneiG[iplane].neighbor[3][row][col])->fval2)+((planeneiG[iplane].neighbor[4][row][col])->fval2));

  result.fval3 = (((planeneiG[iplane].neighbor[0][row][col])->fval3)+((planeneiG[iplane].neighbor[1][row][col])->fval3)+((planeneiG[iplane].neighbor[2][row][col])->fval3)+((planeneiG[iplane].neighbor[3][row][col])->fval3)+((planeneiG[iplane].neighbor[4][row][col])->fval3));

  result.fval4 = (((planeneiG[iplane].neighbor[0][row][col])->fval4)+((planeneiG[iplane].neighbor[1][row][col])->fval4)+((planeneiG[iplane].neighbor[2][row][col])->fval4)+((planeneiG[iplane].neighbor[3][row][col])->fval4)+((planeneiG[iplane].neighbor[4][row][col])->fval4));

  result.fval5 = (((planeneiG[iplane].neighbor[0][row][col])->fval5)+((planeneiG[iplane].neighbor[1][row][col])->fval5)+((planeneiG[iplane].neighbor[2][row][col])->fval5)+((planeneiG[iplane].neighbor[3][row][col])->fval5)+((planeneiG[iplane].neighbor[4][row][col])->fval5));

  return result;
}

TYPE2 GetNeighborS(TYPE2** a,int row,int col,int nei)
{
  int iplane;
  if(nei<0 || nei>8){
    fprintf(stderr,"GetNeighbor(): wrong value in the forth argument\n");
    exit(1);
  }

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  return(*(planeneiG[iplane].neighbor[nei][row][col]));
}

TYPE2 RandomMooreS8(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(1,8);
  return(*(planeneiG[iplane].neighbor[nei][row][col]));
}

TYPE2 RandomMooreS9(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(0,8);
  return(*(planeneiG[iplane].neighbor[nei][row][col]));
}

TYPE2 RandomNeumannS4(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(1,4);
  return(*(planeneiG[iplane].neighbor[nei][row][col]));
}

TYPE2 RandomNeumannS5(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(0,4);
  return(*(planeneiG[iplane].neighbor[nei][row][col]));
}

TYPE2 CountMooreS8(TYPE2 **plane,int aval,int row,int col)
{
  int iplane;
  TYPE2 result;
  
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }

  result.val = (((planeneiG[iplane].neighbor[1][row][col])->val==aval)+((planeneiG[iplane].neighbor[2][row][col])->val==aval)+((planeneiG[iplane].neighbor[3][row][col])->val==aval)+((planeneiG[iplane].neighbor[4][row][col])->val==aval)+((planeneiG[iplane].neighbor[5][row][col])->val==aval)+((planeneiG[iplane].neighbor[6][row][col])->val==aval)+((planeneiG[iplane].neighbor[7][row][col])->val==aval)+((planeneiG[iplane].neighbor[8][row][col])->val==aval));

  result.val2 = (((planeneiG[iplane].neighbor[1][row][col])->val2==aval)+((planeneiG[iplane].neighbor[2][row][col])->val2==aval)+((planeneiG[iplane].neighbor[3][row][col])->val2==aval)+((planeneiG[iplane].neighbor[4][row][col])->val2==aval)+((planeneiG[iplane].neighbor[5][row][col])->val2==aval)+((planeneiG[iplane].neighbor[6][row][col])->val2==aval)+((planeneiG[iplane].neighbor[7][row][col])->val2==aval)+((planeneiG[iplane].neighbor[8][row][col])->val2==aval));

  result.val3 = (((planeneiG[iplane].neighbor[1][row][col])->val3==aval)+((planeneiG[iplane].neighbor[2][row][col])->val3==aval)+((planeneiG[iplane].neighbor[3][row][col])->val3==aval)+((planeneiG[iplane].neighbor[4][row][col])->val3==aval)+((planeneiG[iplane].neighbor[5][row][col])->val3==aval)+((planeneiG[iplane].neighbor[6][row][col])->val3==aval)+((planeneiG[iplane].neighbor[7][row][col])->val3==aval)+((planeneiG[iplane].neighbor[8][row][col])->val3==aval));

  result.val4 = (((planeneiG[iplane].neighbor[1][row][col])->val4==aval)+((planeneiG[iplane].neighbor[2][row][col])->val4==aval)+((planeneiG[iplane].neighbor[3][row][col])->val4==aval)+((planeneiG[iplane].neighbor[4][row][col])->val4==aval)+((planeneiG[iplane].neighbor[5][row][col])->val4==aval)+((planeneiG[iplane].neighbor[6][row][col])->val4==aval)+((planeneiG[iplane].neighbor[7][row][col])->val4==aval)+((planeneiG[iplane].neighbor[8][row][col])->val4==aval));

  result.val5 = (((planeneiG[iplane].neighbor[1][row][col])->val5==aval)+((planeneiG[iplane].neighbor[2][row][col])->val5==aval)+((planeneiG[iplane].neighbor[3][row][col])->val5==aval)+((planeneiG[iplane].neighbor[4][row][col])->val5==aval)+((planeneiG[iplane].neighbor[5][row][col])->val5==aval)+((planeneiG[iplane].neighbor[6][row][col])->val5==aval)+((planeneiG[iplane].neighbor[7][row][col])->val5==aval)+((planeneiG[iplane].neighbor[8][row][col])->val5==aval));

  return result;
}

TYPE2 CountMooreS9(TYPE2 **plane,int aval,int row,int col)
{
  int iplane;
  TYPE2 result;
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }
  result.val = (((planeneiG[iplane].neighbor[0][row][col])->val==aval)+((planeneiG[iplane].neighbor[1][row][col])->val==aval)+((planeneiG[iplane].neighbor[2][row][col])->val==aval)+((planeneiG[iplane].neighbor[3][row][col])->val==aval)+((planeneiG[iplane].neighbor[4][row][col])->val==aval)+((planeneiG[iplane].neighbor[5][row][col])->val==aval)+((planeneiG[iplane].neighbor[6][row][col])->val==aval)+((planeneiG[iplane].neighbor[7][row][col])->val==aval)+((planeneiG[iplane].neighbor[8][row][col])->val==aval));

  result.val2 = (((planeneiG[iplane].neighbor[0][row][col])->val2==aval)+((planeneiG[iplane].neighbor[1][row][col])->val2==aval)+((planeneiG[iplane].neighbor[2][row][col])->val2==aval)+((planeneiG[iplane].neighbor[3][row][col])->val2==aval)+((planeneiG[iplane].neighbor[4][row][col])->val2==aval)+((planeneiG[iplane].neighbor[5][row][col])->val2==aval)+((planeneiG[iplane].neighbor[6][row][col])->val2==aval)+((planeneiG[iplane].neighbor[7][row][col])->val2==aval)+((planeneiG[iplane].neighbor[8][row][col])->val2==aval));

  result.val3 = (((planeneiG[iplane].neighbor[0][row][col])->val3==aval)+((planeneiG[iplane].neighbor[1][row][col])->val3==aval)+((planeneiG[iplane].neighbor[2][row][col])->val3==aval)+((planeneiG[iplane].neighbor[3][row][col])->val3==aval)+((planeneiG[iplane].neighbor[4][row][col])->val3==aval)+((planeneiG[iplane].neighbor[5][row][col])->val3==aval)+((planeneiG[iplane].neighbor[6][row][col])->val3==aval)+((planeneiG[iplane].neighbor[7][row][col])->val3==aval)+((planeneiG[iplane].neighbor[8][row][col])->val3==aval));

  result.val4 = (((planeneiG[iplane].neighbor[0][row][col])->val4==aval)+((planeneiG[iplane].neighbor[1][row][col])->val4==aval)+((planeneiG[iplane].neighbor[2][row][col])->val4==aval)+((planeneiG[iplane].neighbor[3][row][col])->val4==aval)+((planeneiG[iplane].neighbor[4][row][col])->val4==aval)+((planeneiG[iplane].neighbor[5][row][col])->val4==aval)+((planeneiG[iplane].neighbor[6][row][col])->val4==aval)+((planeneiG[iplane].neighbor[7][row][col])->val4==aval)+((planeneiG[iplane].neighbor[8][row][col])->val4==aval));

  result.val5 = (((planeneiG[iplane].neighbor[0][row][col])->val5==aval)+((planeneiG[iplane].neighbor[1][row][col])->val5==aval)+((planeneiG[iplane].neighbor[2][row][col])->val5==aval)+((planeneiG[iplane].neighbor[3][row][col])->val5==aval)+((planeneiG[iplane].neighbor[4][row][col])->val5==aval)+((planeneiG[iplane].neighbor[5][row][col])->val5==aval)+((planeneiG[iplane].neighbor[6][row][col])->val5==aval)+((planeneiG[iplane].neighbor[7][row][col])->val5==aval)+((planeneiG[iplane].neighbor[8][row][col])->val5==aval));

  return result;
}

TYPE2 CountNeumannS4(TYPE2** plane,int aval,int row,int col)
{
  int iplane;
  TYPE2 result;
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }
  result.val = (((planeneiG[iplane].neighbor[1][row][col])->val==aval)+((planeneiG[iplane].neighbor[2][row][col])->val==aval)+((planeneiG[iplane].neighbor[3][row][col])->val==aval)+((planeneiG[iplane].neighbor[4][row][col])->val==aval));

  result.val2 = (((planeneiG[iplane].neighbor[1][row][col])->val2==aval)+((planeneiG[iplane].neighbor[2][row][col])->val2==aval)+((planeneiG[iplane].neighbor[3][row][col])->val2==aval)+((planeneiG[iplane].neighbor[4][row][col])->val2==aval));

  result.val3 = (((planeneiG[iplane].neighbor[1][row][col])->val3==aval)+((planeneiG[iplane].neighbor[2][row][col])->val3==aval)+((planeneiG[iplane].neighbor[3][row][col])->val3==aval)+((planeneiG[iplane].neighbor[4][row][col])->val3==aval));

  result.val4 = (((planeneiG[iplane].neighbor[1][row][col])->val4==aval)+((planeneiG[iplane].neighbor[2][row][col])->val4==aval)+((planeneiG[iplane].neighbor[3][row][col])->val4==aval)+((planeneiG[iplane].neighbor[4][row][col])->val4==aval));

  result.val5 = (((planeneiG[iplane].neighbor[1][row][col])->val5==aval)+((planeneiG[iplane].neighbor[2][row][col])->val5==aval)+((planeneiG[iplane].neighbor[3][row][col])->val5==aval)+((planeneiG[iplane].neighbor[4][row][col])->val5==aval));

  return result;
}

TYPE2 CountNeumannS5(TYPE2** plane,int aval,int row,int col)
{
  int iplane;
  TYPE2 result;
  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == plane)
      break;
  }
  result.val = (((planeneiG[iplane].neighbor[0][row][col])->val==aval)+((planeneiG[iplane].neighbor[1][row][col])->val==aval)+((planeneiG[iplane].neighbor[2][row][col])->val==aval)+((planeneiG[iplane].neighbor[3][row][col])->val==aval)+((planeneiG[iplane].neighbor[4][row][col])->val==aval));

  result.val2 = (((planeneiG[iplane].neighbor[0][row][col])->val2==aval)+((planeneiG[iplane].neighbor[1][row][col])->val2==aval)+((planeneiG[iplane].neighbor[2][row][col])->val2==aval)+((planeneiG[iplane].neighbor[3][row][col])->val2==aval)+((planeneiG[iplane].neighbor[4][row][col])->val2==aval));

  result.val3 = (((planeneiG[iplane].neighbor[0][row][col])->val3==aval)+((planeneiG[iplane].neighbor[1][row][col])->val3==aval)+((planeneiG[iplane].neighbor[2][row][col])->val3==aval)+((planeneiG[iplane].neighbor[3][row][col])->val3==aval)+((planeneiG[iplane].neighbor[4][row][col])->val3==aval));

  result.val4 = (((planeneiG[iplane].neighbor[0][row][col])->val4==aval)+((planeneiG[iplane].neighbor[1][row][col])->val4==aval)+((planeneiG[iplane].neighbor[2][row][col])->val4==aval)+((planeneiG[iplane].neighbor[3][row][col])->val4==aval)+((planeneiG[iplane].neighbor[4][row][col])->val4==aval));

  result.val5 = (((planeneiG[iplane].neighbor[0][row][col])->val5==aval)+((planeneiG[iplane].neighbor[1][row][col])->val5==aval)+((planeneiG[iplane].neighbor[2][row][col])->val5==aval)+((planeneiG[iplane].neighbor[3][row][col])->val5==aval)+((planeneiG[iplane].neighbor[4][row][col])->val5==aval));

  return result;
}


void GetNeighborC(TYPE2** a,int row,int col,int nei,int*neirow,int*neicol)
{
  if(nei<0 || nei>8){
    fprintf(stderr,"GetNeighbor(): wrong value in the forth argument\n");
    exit(1);
  }

  if(neirow!=NULL)
    *(neirow)=neighborG[nei][row][col].row;
  if(neicol!=NULL)
    *(neicol)=neighborG[nei][row][col].col;
}

void RandomMooreC8(TYPE2** a,int row,int col,int*neirow,int*neicol)
{
  int nei;
  nei = genrand_int(1,8);

  if(neirow!=NULL)
    *(neirow)=neighborG[nei][row][col].row;
  if(neicol!=NULL)
    *(neicol)=neighborG[nei][row][col].col;
}

void RandomMooreC9(TYPE2** a,int row,int col,int*neirow,int*neicol)
{
  int nei;
  nei = genrand_int(0,8);

  if(neirow!=NULL)
    *(neirow)=neighborG[nei][row][col].row;
  if(neicol!=NULL)
    *(neicol)=neighborG[nei][row][col].col;
}

void RandomNeumannC4(TYPE2** a,int row,int col,int*neirow,int*neicol)
{
  int nei;
  nei = genrand_int(1,4);

  if(neirow!=NULL)
    *(neirow)=neighborG[nei][row][col].row;
  if(neicol!=NULL)
    *(neicol)=neighborG[nei][row][col].col;
}

void RandomNeumannC5(TYPE2** a,int row,int col,int*neirow,int*neicol)
{
  int nei;
  nei = genrand_int(0,4);

  if(neirow!=NULL)
    *(neirow)=neighborG[nei][row][col].row;
  if(neicol!=NULL)
    *(neicol)=neighborG[nei][row][col].col;
}

TYPE2* GetNeighborP(TYPE2** a,int row,int col,int nei)
{
  int iplane;
  if(nei<0 || nei>8){
    fprintf(stderr,"GetNeighbor(): wrong value in the forth argument\n");
    exit(1);
  }

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  return(planeneiG[iplane].neighbor[nei][row][col]);
}

TYPE2* RandomMooreP8(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(1,8);
  return(planeneiG[iplane].neighbor[nei][row][col]);
}

TYPE2* RandomMooreP9(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(0,8);
  return(planeneiG[iplane].neighbor[nei][row][col]);
}

TYPE2* RandomNeumannP4(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(1,4);
  return(planeneiG[iplane].neighbor[nei][row][col]);
}

TYPE2* RandomNeumannP5(TYPE2** a,int row,int col)
{
  int iplane;
  int nei;

  /* we search through planeneiG[] to find which plane user is
     refering to */
  for(iplane=0;iplane<nplane;iplane++){
    if(planeneiG[iplane].plane == a)
      break;
  }

  nei = genrand_int(0,4);
  return(planeneiG[iplane].neighbor[nei][row][col]);
}

int gotMouse(void)
{
  char answer[BUFSIZE];

  poseG = 0;
  while(1){
    printf("Enter: Contine (c), OneStep (o), Quit (q): ?");
    fgets(answer,BUFSIZE,stdin);
    switch(answer[0]){
    case 'c':
      return 0;
    case 'o':
      poseG = 1;
      return 0;
    case 'q':
      stopG = 1;
      return 0;
    default:
      printf("Unknown command: %s\n",answer);
      break;
    }
  }


}


/* Make planes by calling New2() and make also neighborhood
   planes by calling NeighSet(), NeighSetP() and
   MargolusNeigh(). */
void MakePlane(TYPE2*** Arg1,...)
{
  int i,j;
  va_list ap;
  TYPE2*** planeP;

  planeneiG = (struct PlaneNei*) malloc(sizeof(struct PlaneNei)*nplane);

  /* if there is more plane, then there should be more arguments */
  if(nplane>1)
    va_start(ap,Arg1);

  for(i=0;i<nplane;i++){
    if(i==0){
      planeP = Arg1;
    }
    else{
      planeP = va_arg(ap,TYPE2***);
    }

    *planeP = New2();

    /* here, we keep track of planes created. */
    planeneiG[i].plane = *planeP;

    /* we also make address based neighborhood planes */
    for(j=0;j<9;j++)
      planeneiG[i].neighbor[j] = NeighSetP(*planeP,j);
  }

  /* we also make coordinate based neighborhood planes */
  for(i=0;i<9;i++)
    neighborG[i]=NeighSet(i);
  margolusG = MargolusNeigh();

}

void SavePlane(char *fname,TYPE2 **Arg1,...)
{
  int i,j,k;
  FILE *fout;
  va_list ap;
  TYPE2** plane;


  fout = fopen(fname,"w");
  if(fout == NULL){
    fprintf(stderr,"SavePlane: can't open %s\n",fname);
    exit(1);
  }

  fprintf(fout,"nrow=%d\nncol=%d\nplane=%d\n",nrow,ncol,nplane);

  if(nplane > 1){
    va_start(ap,Arg1);
  }
  for(k=0;k<nplane;k++){
    if(k==0){
      plane = Arg1;
    }else{
      plane=va_arg(ap,TYPE2**);
    }
    for(i=1;i<=nrow;i++)
      for(j=1;j<=ncol;j++){
	fprintf(fout,"%d %d %d %d %d %f %f %f %f %f\n",plane[i][j].val,plane[i][j].val2,plane[i][j].val3,plane[i][j].val4,plane[i][j].val5,plane[i][j].fval,plane[i][j].fval2,plane[i][j].fval3,plane[i][j].fval4,plane[i][j].fval5);
      }
  }

  fclose(fout);
}


void ReadSavedData(char *fname,int npl,TYPE2** Arg2,...)
{
  int i,j,k;
  FILE *fin;
  char line[1000];
  int nr,nc;
  TYPE2** plane;
  va_list ap;
  char *val;

  fin = fopen(fname,"r");
  if(fin == NULL){
    fprintf(stderr,"ReadSavedData: can't open %s\n",fname);
    exit(1);
  }

  /* nrow */
  fgets(line,1000,fin);
  nr=atoi(strpbrk(line,"=")+1);
  /* ncol */
  fgets(line,1000,fin);
  nc=atoi(strpbrk(line,"=")+1);
  /* nplane */
  fgets(line,1000,fin);
  if(npl!=atoi(strpbrk(line,"=")+1)){
    fprintf(stderr,"ReadSavedData: the number of planes is different?\n");
    exit(1);
  }



  if(nrow<nr||ncol<nc){
    fprintf(stderr,"ReadSavedData: nrow or ncol are too small to read the file\n");
    exit(1);
  }

  if(npl > 1){
    va_start(ap,Arg2);
  }
  for(k=0;k<npl;k++){
    if(k==0){
      plane = Arg2;
    }else{
      plane = va_arg(ap,TYPE2**);
    }

    for(i=1;i<=nr;i++)
      for(j=1;j<=nc;j++){
	fgets(line,1000,fin);

	val = strtok(line," ");
	plane[i][j].val = atoi(val);

	val = strtok(NULL," ");
	plane[i][j].val2 = atoi(val);

	val = strtok(NULL," ");
	plane[i][j].val3 = atoi(val);

	val = strtok(NULL," ");
	plane[i][j].val4 = atoi(val);

	val = strtok(NULL," ");
	plane[i][j].val5 = atoi(val);

	val = strtok(NULL," ");
	plane[i][j].fval = atof(val);

	val = strtok(NULL," ");
	plane[i][j].fval2 = atof(val);

	val = strtok(NULL," ");
	plane[i][j].fval3 = atof(val);

	val = strtok(NULL," ");
	plane[i][j].fval4 = atof(val);

	val = strtok(NULL," ");
	plane[i][j].fval5 = atof(val);
      }
  }
}

/*
void ReadCrap(char* fname,TYPE2** plane)
{
  int i,j;
  FILE *fin;
  char val[2]="";
  fin = fopen(fname,"r");

  for(i=1;i<=130;i++){
    for(j=1;j<=130;j++){
      val[0]=(char)fgetc(fin);
      plane[i][j].val = atoi(val);
    }
    fgetc(fin);
  }
}
*/

struct space_time_rec {
  TYPE2 **plane;  
  int current; /* the row of the destination, where we should draw somehting */
};
void SpaceTimePlot(TYPE2 **dest, TYPE2 **source)
{
  int nr=nrow,nc=ncol;
  int i;
  static struct space_time_rec *list;
  static int init_flag = 0;
  int planeID;

  if(init_flag == 0){
    list = (struct space_time_rec*) malloc(sizeof(struct space_time_rec)*nplane);
    for(i=0;i<nplane;i++){
      (list+i)->plane=NULL;
      (list+i)->current = 1;
    }
    init_flag = 1;
  }

  /* search through the list */
  planeID = -1;
  for(i=0;i<nplane;i++){
    if((list+i)->plane == dest){
      planeID = i;
      break;
    }
    else if((list+i)->plane==NULL){
      /* insert this plane */
      (list+i)->plane = dest;
      planeID = i;
      break;
    }
  }

  /* weak error cheack */
  if(planeID == -1){
    fprintf(stderr,"SpaceTimePlot: a bug. please report this\n");
    exit(1);
  }

  memcpy(dest[(list+planeID)->current],source[nr/2],sizeof(TYPE2)*(nc+2));
  (list+planeID)->current = ((list+planeID)->current % nr)+1;
  
}


/* A easier variant of InitailSetS() */
int InitialSet(TYPE2** plane,int ntype,int empty_int,...)
{
  va_list ap;

  TYPE2 *type_p; /* the types of cell */
  double *port_p; /* portion */
  int *numb_p; /* acutual number calculated from the portion */

  double port_total=0.;
  int nr,nc,pos,maxpos;
  int i,ii,jj;
  TYPE2 empty={0,0,0,0,0,0.0,0.0,0.0,0.0,0.0};

  /* If there is no other type thant "empty" */
  if(ntype<=0){
    empty.val = empty_int;
    Fill2(plane,empty);
    return 0;
  }

  nr = nrow;
  nc = ncol;

  /* get argument */
  type_p = (TYPE2*)malloc(sizeof(TYPE2)*ntype);
  port_p = (double*)malloc(sizeof(double)*ntype);
  numb_p = (int*)malloc(sizeof(int)*ntype);

  va_start(ap,empty_int);
  empty.val = empty_int;
  for(i=0;i<ntype;i++){
    *(type_p+i) = empty;
    (type_p+i)->val = va_arg(ap,int);
    *(port_p+i) =  va_arg(ap,double);

    /* if portion<0.0, make it 0.0 */
    if(*(port_p+i)<0.0)
      *(port_p+i) = 0.0;
    port_total += *(port_p+i);
  }
  
  /* if P1+P2 > 1, then P1=P1/(P1+P2) */
  if(port_total>1.0){
    for(i=0;i<ntype;i++){
      *(port_p+i) = (*(port_p+i))/port_total;
    }
  }

  /* 
     Set the number of cells which will be S1, S2, ... 
     The fractional part will be truncated so that the sum of the
     number will not be nrow*ncol.
  */
  for(i=0;i<ntype;i++){
    *(numb_p+i) = (*(port_p+i))*nr*nc;
  }


  /* From here, we will set the state of the cell to S1, S2 ...*/

  /* "pos" is the number of the current cell [i][j]. */
  pos = nc; /* if pos==ncol, then [i][j]=[1][1].
	       if pos==ncol+1, then [i][j]=[1][2]
	       if pos==ncol+ncol+2, then [i][j]=[2][3]
	       if pos==x, then [i][j]=[x/ncol][x%ncol+1]
	       if pos==ncol+ncol*nrow-1, then [i][j]=[nrow][ncol]
	    */
  maxpos = nc+nc*nr-1;

  /* fill the plane with the given state, S0 S1 S2 ..., we will
     shafful the plane later */
  for(i=0;i<ntype;i++){
    while(*(numb_p+i)>0){
      ii = pos/nc;
      jj = pos%nc + 1;

      plane[ii][jj] = *(type_p+i);

      pos++;
      (*(numb_p+i))--;
    }
  }
  if(pos>maxpos+1){
    fprintf(stderr,"InitialSet(): this is a bug. Report please");
    exit(1);
  }

  while(pos<=maxpos){
    ii = pos/nc;
    jj = pos%nc + 1;
    
    plane[ii][jj] = empty;
    
    pos++;
  }

  PerfectMix(plane);

  free(type_p);
  free(port_p);
  free(numb_p);

  return 0;
}

/* A easier variant of InitailSetS() */
int InitialSpot(TYPE2** plane,int type, int size, int offset)
{
	for (int row=1; row <= nrow; row++)
	{
		for (int col=1; col <= ncol; col++)
		{
			// Draw a circle
			int rowdraw, coldraw;
			rowdraw = row+offset;
			coldraw = col;
			if( (pow((col-ncol/2),2) + pow((row-nrow/2),2) ) < size)
			{
				plane[rowdraw % nrow+1][coldraw % ncol+1].val = type;
			}
		}
	}

  return 0;
}

/* This is exaclty the same as InitialPlaneSet() in cash2.c */
int InitialSetS(TYPE2** a,int ntype,TYPE2 empty,...)
{
  va_list ap;
  TYPE2 *type_p; /* the types of cell */
  double *port_p; /* portion */
  int *numb_p; /* acutual number calculated from the portion */
  double port_total=0.;
  int nr,nc,pos,maxpos;
  int i,ii,jj;

  /* If there is no other type thant "empty" */
  if(ntype<=0){
    Fill2(a,empty);
    return 0;
  }

  nr = nrow;
  nc = ncol;

  /* get argument */
  type_p = (TYPE2*)malloc(sizeof(TYPE2)*ntype);
  port_p = (double*)malloc(sizeof(double)*ntype);
  numb_p = (int*)malloc(sizeof(int)*ntype);

  va_start(ap,empty);

  for(i=0;i<ntype;i++){
    *(type_p+i) = va_arg(ap,TYPE2);
    *(port_p+i) =  va_arg(ap,double);

    /* if portion<0.0, make it 0.0 */
    if(*(port_p+i)<0.0)
      *(port_p+i) = 0.0;
    port_total += *(port_p+i);
  }
  
  /* if P1+P2 > 1, then P1=P1/(P1+P2) */
  if(port_total>1.0){
    for(i=0;i<ntype;i++){
      *(port_p+i) = (*(port_p+i))/port_total;
    }
  }

  /* 
     Set the number of cells which will be S1, S2, ... 
     The fractional part will be truncated so that the sum of the
     number will not be nrow*ncol.
  */
  for(i=0;i<ntype;i++){
    *(numb_p+i) = (*(port_p+i))*nr*nc;
  }


  /* From here, we will set the state of the cell to S1, S2 ...*/

  /* "pos" is the number of the current cell [i][j]. */
  pos = nc; /* if pos==ncol, then [i][j]=[1][1].
	       if pos==ncol+1, then [i][j]=[1][2]
	       if pos==ncol+ncol+2, then [i][j]=[2][3]
	       if pos==x, then [i][j]=[x/ncol][x%ncol+1]
	       if pos==ncol+ncol*nrow-1, then [i][j]=[nrow][ncol]
	    */
  maxpos = nc+nc*nr-1;

  /* fill the plane with the given state, S0 S1 S2 ..., we will
     shafful the plane later */
  for(i=0;i<ntype;i++){
    while(*(numb_p+i)>0){
      ii = pos/nc;
      jj = pos%nc + 1;

      a[ii][jj] = *(type_p+i);

      pos++;
      (*(numb_p+i))--;
    }
  }
  if(pos>maxpos+1){
    fprintf(stderr,"InitialPlaneSet(): this is a bug. Report please");
    exit(1);
  }

  while(pos<=maxpos){
    ii = pos/nc;
    jj = pos%nc + 1;
    
    a[ii][jj] = empty;
    
    pos++;
  }

  PerfectMix(a);

  free(type_p);
  free(port_p);
  free(numb_p);

  return 0;
}

/***********************************
 Margolus diffusion with obstacles by Anton
**********************************/

/* Some auxilary functions for rotating the particles
 * 
 * These functions should be inlined.... as if they were macros.
 */
// inline void rotate_square( TYPE2 *a, TYPE2 *b, TYPE2 *c, TYPE2 *d ) {
//     TYPE2 aux;
// 
//     if( genrand_real1() < 0.5 ) {
//         aux = *a;
//         *a = *b;
//         *b = *c;
//         *c = *d;
//         *d = aux;
//     } else {
//         aux = *d;
//         *d = *c;
//         *c = *b;
//         *b = *a;
//         *a = aux;
//     }
// }
// 
// inline void rotate_triangle( TYPE2 *a, TYPE2 *b, TYPE2 *c ) {
//     TYPE2 aux;
// 
//     if( genrand_real1() < 0.5 ) {
//         aux = *a;
//         *a = *b;
//         *b = *c;
//         *c = aux;
//     } else {
//         aux = *c;
//         *c = *b;
//         *b = *a;
//         *a = aux;
//     }
// }
// 
// inline void rotate_pair( TYPE2 *a, TYPE2 *b ) {
//     TYPE2 aux;
// 
//     aux = *a;
//     *a = *b;
//     *b = aux;
// }
// 
// 
// /* Modified Margolis diffusion. Now we can diffuse with obstacles
//  * in our grid as well.
//  * 
//  * The function expects two planes, a diffusion plane and an
//  * obstacle plane. For the obstacle plane one needs to supply an
//  * example obstacle. (There may be more states in the plane of
//  * which only one is an obstacle.)
//  *
//  * Assuming both planes have equal dimensions.
//  */
// void MargolusWithObstacle( TYPE2** d_plane, TYPE2** o_plane, int obstacle, 
//                       MARGOLUS*** margolus, int phase ) {
//   int my_nrow, my_ncol;
//   int i, j, ii;
//   unsigned int mask;
//   int sq_row[ 4 ], sq_col[ 4 ];
//   
//   my_nrow = nrow - 1;
//   my_ncol = ncol - 1;
// 
//   if( boundary == WRAP ) {
//     // use margolis neighbourhood and rotate according to obstacles
//     for( i = 1; i <= my_nrow; i += 2 ) {
//       for( j = 1; j <= my_ncol; j += 2 ) {
// 
// 	// 1. determine nbh
// 
// 	/* (i,j)=HERE */
// 	sq_row[ 0 ] = i;
// 	sq_col[ 0 ] = j;
// 	/* (cw_r,cw_c)=CW_of_(i,j) */
// 	sq_row[ 1 ] = margolus[phase][i][j].m[CW].row;
// 	sq_col[ 1 ] = margolus[phase][i][j].m[CW].col;
// 	/* (ccw_r,ccw_c)=CCW_f_(i,j) */
// 	sq_row[ 2 ] = margolus[phase][i][j].m[CCW].row;
// 	sq_col[ 2 ] = margolus[phase][i][j].m[CCW].col;
// 	/* (opp_r,opp_c)=OPP_of_(i,j) */
// 	sq_row[ 3 ] = margolus[phase][i][j].m[OPP].row;
// 	sq_col[ 3 ] = margolus[phase][i][j].m[OPP].col;
// 
// 	// 2. look for obstacles
// 	mask = 0;
// 	for( ii = 0; ii < 4; ++ii ) {
// 	  // FIX whole structure comparison?
// 	  if( o_plane[ sq_row[ ii ] ][ sq_col[ ii ] ].val == obstacle ) {
// 	    mask = mask | ( 1 << ii );
// 	  }
// 	}
//         
// 	// 3. rotate depending on nr and position of obstacles
// 	/* The neighbourhood looks like this:
// 	 * 
// 	 * ---------
// 	 * | 0 | 1 |
// 	 * ---------
// 	 * | 2 | 3 |
// 	 * ---------
// 	 */
// 	switch( mask ) {
// 	case 0:
// 	  // rotate 0, 1, 2, 3
// 	  rotate_square( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
// 			 &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ],
// 			 &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ],
// 			 &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ] );
// 	  break;
// 	case 1:
// 	  // rotate 0, 1, 2
// 	  rotate_triangle( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
// 			   &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ],
// 			   &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ] );
// 	  break;
// 	case 2:
// 	  // rotate 0, 1, 3
// 	  rotate_triangle( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
// 			   &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ],
// 			   &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ] );
// 	  break;
// 	case 3:
// 	  // rotate 0, 1
// 	  rotate_pair( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
// 		       &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ] );
// 	  break;
// 	case 4:
// 	  // rotate 0, 2, 3
// 	  rotate_triangle( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
// 			   &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ],
// 			   &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ] );
// 	  break;
// 	case 5:
// 	  // rotate 0, 2
// 	  rotate_pair( &d_plane[ sq_row[ 0 ] ][ sq_col[ 0 ] ],
// 		       &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ] );
// 	  break;
// 	case 8:
// 	  // rotate 1, 2, 3
// 	  rotate_triangle( &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ],
// 			   &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ],
// 			   &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ] );
// 	  break;
// 	case 10:
// 	  // rotate 1, 3
// 	  rotate_pair( &d_plane[ sq_row[ 1 ] ][ sq_col[ 1 ] ],
// 		       &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ] );
// 	  break;
// 	case 12:
// 	  // rotate 2, 3
// 	  rotate_pair( &d_plane[ sq_row[ 2 ] ][ sq_col[ 2 ] ],
// 		       &d_plane[ sq_row[ 3 ] ][ sq_col[ 3 ] ] );
// 	  break;
// 	}
//       }
//     }
//   } else if( boundary == FIXED ) {
//     fprintf(stderr,"MargolusWithObstacle() Use WRAP boundary\n");
//     exit(-1);
//   }
// }

// /*
// void  ObstacleMargolus( TYPE2 **a, TYPE2 **o_plane, int obstacle )
// {
//   int i;
//   int planeID;
//   static int init_flag=0;
//   static struct phase_list *list;
//  
//   /* initilization */
//   if(init_flag==0){
//     list = (struct phase_list *)malloc(sizeof(struct phase_list)*nplane);
//     for(i=0;i<nplane;i++){
//       (list+i)->phase = 0;
//       (list+i)->plane = NULL;
//     }
//     init_flag = 1;
//   }
// 
//   /* search "a_name" in the list */
//   planeID = -1;
//   i=0;
//   while(i<nplane){
//     if((list+i)->plane == a){
//       planeID = i;
//       break;
//     }
//     else if((list+i)->plane==NULL){
//       /* insert this plane */
//       (list+i)->plane = a;
//       planeID = i;
//       break;
//     }
//     i++;
//   }
// 
//   /* weak error cheack */
//   if(planeID == -1){
//     fprintf(stderr,"ObstacleMargolis: a bug. please report this\n");
//     exit(1);
//   }
// 
//   /* do Margolus diffusion */
//   MargolusWithObstacle( a, o_plane, obstacle, margolusG, (list+planeID)->phase );
// 
//   /* incliment the phase */
//   (list+planeID)->phase = ((list+planeID)->phase+1)%2;
// }*/

int countGlobal(TYPE2 **world, int num){

  int pop = 0;
  int i,j;
  for(i=1;i<=nrow;i++){                     
    for(j=1;j<=ncol;j++){                     
     if(world[i][j].val == num)
          pop++;
     }
  }
  return pop;
}

void PlotXY(double x, double y)
{
  static int count;
  static int initial_flag = 0;

  if(initial_flag == 0){
    InitXmgrace();
    initial_flag =1 ;
  }
   
  GracePrintf("g0.s%d point %f, %f",1,x,y);

  if(count%10==0)
    GracePrintf("redraw");
  
  if(count%100==0){
    GracePrintf("focus g0");
    GracePrintf("autoscale");
    GracePrintf("redraw");
  }

  count++;
}
