/*
caldep + fragment + sphere
*/
//#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include "struct.h"
#include "func.h"
#include "mrc.h"
//#include "scoring.h"

#define PDB_STRLEN 55

void malloc_error(char *a){
 fprintf(stderr,"malloc error in %s\n",a);
 exit(0);
}
double gettimeofday_sec()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

int readlist(char *fname,char **list){
 int num=0;
 FILE *fp;
 int len;
 if((fp=fopen(fname,"r"))==NULL)
  return FALSE;
 while(fgets(list[num],LIN,fp)!=NULL){
  len=strlen(list[num]);
  list[num][len-1]='\0';//ignore terminal \n
  num++;
 }
 fclose(fp);
 return TRUE;
}

int line_num(char *fname){
 int num=0;
 FILE *fp;
 char line[LIN];
 if((fp=fopen(fname,"r"))==NULL)
  return FALSE;
 while(fgets(line,LIN,fp)!=NULL){
  num++;
 }
 fclose(fp);
 return num;
}


CMD cmd;

int main(int argc, char **argv)
{
 double t1=gettimeofday_sec();
 double t4;
 POINTS pt;
 MRC mrc;
 GRAPH g;
 TREE mst;
 //Get Options
 if(chkcmdline(argc,argv,&cmd)==FALSE)
  return(0);

 //Set threads
 if(cmd.Nthr < omp_get_num_procs()){
  omp_set_num_threads(cmd.Nthr);
 }else{
  omp_set_num_threads(omp_get_num_procs());
 }
 
 if(readmrc(&mrc,cmd.filename))
  return(0);

 //Convert to cubic
 // if(ToCubic(&mrc))
 // return(0);

 //return 0;
 //upsampling
 if(upsampling(&mrc,cmd.map_t))
  return(0);

 printf("#Nact= %d\n",mrc.Nact);

 //out_situs(&mrc);

 if(cmd.Mode==3){//fast LDPs only
  if(fastLDP(&mrc,&pt))
   return(0);
  //if(MergePoints(&mrc,&pt))
  // return(0);
  ShowLDP(&mrc,&pt);
  t4=gettimeofday_sec();
  printf("#FINISHED TOTAL TIME= %f\n",t4-t1);
  return 0;
 }

 //Mean Shifting
 if(meanshift(&mrc,&pt))
  return(0);

 if(MergePoints(&mrc,&pt))
  return(0);

 //Set Up Graph and MST
 if(SetUpGraph(&pt,&g,&mrc,&mst))
  return(0);

 if(cmd.Mode==1){//Graph
  ShowModel(&mrc,&pt);
  ShowGraph(&g);
  return 0;
 }
 if(cmd.Mode==2){//MST
  ShowModel(&mrc,&pt);
  //ShowOri(&mrc,&pt);
  ShowTree(&g,&mst);
  return 0;
 }

 //Optimize Graph by Tabu seatch
 TREE results[100];//Max 100 simu
 if(Tabu(&g,&mst,results))
  return 0;

 //if(PairExhaust(&g,&mst,results))
 // return 0;
 
 //Volume and density data
 ShowPath2(&mrc,&pt,&g,results,cmd.Nsim);

 t4=gettimeofday_sec();
 printf("#FINISHED TOTAL TIME= %f\n",t4-t1);
 return 0;

}

