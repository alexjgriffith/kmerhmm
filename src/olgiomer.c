#include "olgiomer.h"

int bp_options(int *n,char ** nmers)
{
  const char ntides[5]= "ACGT";
  int i,j,k,c,cp,acum,outer,t;
  c=1;
  outer=pow(4,(*n-1));
  acum=0;
  cp=0;
  for(t=0;t<*n;t++)
    {
      for(j=0;j<outer;j++)
	{
	  for(i=0;i<4;i++)
	    {
	      for(k=0;k<c;k++){	    
		nmers[acum][t]=ntides[i];
		acum=acum+1;	
	      }
	    }
	}
      outer=pow(4,(*n-2-t));
      c=pow(4,(t+1));      
      acum=0;
    }  	   
  return 0;
}

int baseToValue(char * inchar,int * outint)
{
  if(* inchar=='A' || * inchar=='a')
    *outint= 0;
  else if(* inchar=='C' || * inchar=='c')
    *outint= 1;
  else if(* inchar=='G' || * inchar=='g')
    *outint = 2;
  else if(* inchar=='T' || * inchar=='t')
    *outint =3;
  return 0;
}

int baseScore(int *score,char * string,int start,int end)
{
  int i;
  int k=0;
  int temp=0;
  int value=0;
  *score=0;
  for(i=start;i<end;i++)
    {      
      baseToValue(&string[i],&value);
      temp=*score+ (pow(4,k) * value) ;
      *score=temp;
      k++;
    }
  return 0;
}

int reverseBase(int * comp, int * modul)
{
  if(*modul==3)
    *comp=0;
  if(*modul==2)
    *comp=1;
  if(*modul==1)
    *comp=2;
  if(*modul==0)
    *comp=3;
  return 0;
}

int findCompliment(int * inValue, int * compliment,int * n)
{
  int i;
  int acum;
  int temp;
  int comp;
  int modul;
  *compliment=0;
  acum=* inValue;
  for(i=0;i<*n;i++)
    {
      modul=acum%4;
      reverseBase(&comp,&modul);
      temp=*compliment+ (pow(4,*n-i-1) * comp );
      *compliment=temp;    
      temp=floor(acum/4);
      acum=temp;
    }
  return 0;
}

int stitch_buff(char * buffer_2, int a, int b, char * buffer , int a2 , int b2,char *  stitch){
  int i,j,k;
  for(i=a;i<b;i++)
    {
      stitch[k]=buffer_2[i];
      k++; 
    }
  
  for(i=a2;i<b2;i++)
    {
      stitch[k]=buffer[i];
      k++; 
    }
    
  return 0;
}


int kmerFreq(double * weight, double * score, char * fastaFilename,int *n){
  int i,j,k,r;
  int been_in;
  char buffer[2048];
  char buffer_2[2048];
  char stitch[200];
  int value;
  double temp;
  int cutoff=100;
  FILE * fasta = fopen(fastaFilename,"r");   
  //reads two lines, first >chr:start-end second {AGCTagct}
  // Ns will break it
  // the string line must be a single line
  ///bug here somewhere
  if(fasta==NULL){
    Rprintf("failed to open %s \n",fastaFilename);
    exit(0);
  }
  k=0;
  been_in=0;
  while(fgets(buffer,2048,fasta))
    {
 
      if(buffer[0]!='>'){
	if(been_in==1)
	  {
	    for(i=0;i<*n;i++)
	      {
		stitch_buff(buffer_2,r-*n+i ,r,buffer,0,i,stitch);
		value=0;
		baseScore(&value, stitch,0,*n);
		temp=score[value]+weight[k];
		score[value]=temp;
	      }
	  }
	for(i=0;i<(strlen(buffer)-*n);i++)
	  {
	    value=0;
	    baseScore(&value, buffer,i,i+*n);
	    temp=score[value]+weight[k];
	    score[value]=temp;
	    
	  }
	r=i;
	strcpy(buffer_2,buffer);
	been_in=1;
      }
      if(buffer[0]=='>')
	{
	  if(been_in==1){
	    k++;
	    been_in=0;
	  }
	}
    }   
  fclose(fasta);
  return 0;
}

int kmerPrint(int * n, char ** nmers,  double * score)
{
  int i;
  for(i=0;i<pow(4,*n);i++)
    if(score[i]>0)
      printf("%s:%lf\n",nmers[i],score[i]);
  return 0;
}


void countNMers(char ** fastaFilename, int * n ,double * weight,double * retmer , char ** charmer ){
  int i;
  int length;
  double * score = (double *) calloc(pow(4,*n)+1,sizeof(double));
  char ** nmers=  calloc(pow(4,*n)+1,sizeof(char*)); 
  int * counts = (int *) calloc(pow(4,*n)+1,sizeof(int));
  for(i=0 ;i<(pow(4,*n)+1);i++){
    nmers[i]=calloc(1,sizeof(char)*(*n+1));
    score[i]=0;
  }
  bp_options(n,nmers);
  kmerFreq(weight,score,*fastaFilename, n);
  compMer(score,counts,retmer,n,&length);

  for(i=0;i<length;i++){
    strcpy(charmer[i],nmers[counts[i]]);

  }
  
}


int compMer(double * score, int * counts, double* retmer,int * n2, int * length){
  int save;
  int k,j,i;
  double * outmer;
  int * select;
  outmer = (double *) calloc(pow(4,*n2)+1,sizeof(double));
  select = (int *) calloc(pow(4,*n2)+1,sizeof(int));
  for(i=0;i<pow(4,*n2)+1;i++){
    outmer[i]=0;
    select[i]=0;
  }
  for(i=0;i<pow(4,*n2);i++){
    findCompliment(&i,&save,&*n2);
    if(i<save){
      outmer[i]=score[i]+score[save];
      select[i]=1;
    }
    if(i==save)
      {
	outmer[i]=score[i];
	select[i]=1;	
      }
  }
  k=0;
  for(i=0;i<pow(4,*n2)+1;i++){
    if(select[i]==1){
      counts[k]=i;
      retmer[k]=outmer[i];      
      k++;      
    }
  }
  *length=k;
  return 0;
}


/*int main(int argc, char ** argv)
{
  int i;
  char * fastaFilename = argv[1];
  int n=4;
  double * weight =(double * ) calloc(200000,sizeof(double));
  for(i=0;i<200000;i++)
    weight[i]=1;
  double * score = (double *) calloc(pow(4,n)+1,sizeof(double));
  char ** nmers=  calloc(pow(4,n)+1,sizeof(char*)); 


  char ** charmers=  calloc(pow(4,n)+1,sizeof(char*)); 
  for(i=0 ;i<(pow(4,n)+1);i++){
    charmers[i]=calloc(1,sizeof(char)*(n+1));
    nmers[i]=calloc(1,sizeof(char)*(n+1));
    score[i]=0;
  }
  bp_options(&n,nmers);
  kmerFreq(weight,score,fastaFilename, &n);
  int length;
  int * counts = (int *) calloc(pow(4,n)+1,sizeof(int));
  double * retmer = (double *) calloc(pow(4,n)+1,sizeof(double));
  compMer(score,counts,retmer,&n,&length);
  //kmerPrint(&n,nmers,score);    
  for(i=0;i<length;i++)
    strcpy(charmers[i],nmers[counts[i]]);
      //printf("%s:%lf\n",nmers[counts[i]],retmer[i]);
  return 0;
  }*/

void R_init_myLib(DllInfo *info)
{
  R_registerRoutines(info,cMethods,NULL,NULL,NULL);
}
