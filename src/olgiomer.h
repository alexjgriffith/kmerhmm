#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int bp_options(int *n,char ** nmers);
int baseToValue(char * inchar,int * outint);
int baseScore(int *score,char * string,int start,int end);
int reverseBase(int * comp, int * modul);
int findCompliment(int * inValue, int * compliment,int * n);
int kmerFreq(double * weight, double * score, char * fastaFilename,int *n);
int kmerPrint(int * n, char ** nmers,  double * score);
void countNMers(char ** fastaFilename, int * n ,double * weight,double * retmer , char ** charmer );
int compMer(double * score, int * counts, double* retmer,int * n2, int * length);

R_CMethodDef cMethods[]={
{"countNMers",(DL_FUNC) &countNMers,5,{STRSXP,INTSXP,REALSXP,REALSXP,STRSXP}},
{NULL,NULL,0}
};

void R_init_myLib(DllInfo *info);
