#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SQRT3 1.732050808
#define SQRT6 2.449489743
#define LEN1 60
#define LEN 250

int main(int argc,char *argv[]){

FILE *FPin;
char **atm,*filein, *tmpstr, *comment;
int i,j,ind,nsp,natm,adv, l1,l2,l3, n1,n2,n3;
double x,y,z, xl,yl,zl, v1[3],v2[3],v3[3], or[3], p[3];
char ch;
int fold;

if (argc<2)
{ 
  printf("Usage: gen2xyz filein.gen -s n1 n2 n3 -c 'comment line'\n");
  printf(" -s  supercell option, repeat n1 n2 n3 \n");
  printf("     0 0 0 folds the supercell in a box \n");
  printf(" -c  'comment line within apices' \n");
  return;
}

n1 = 1; n2 = 1; n3 = 1;
comment=malloc(LEN1*sizeof(char));
strcpy(comment," ");

for(i=1; i<argc; i++)
{
   if (strcmp(argv[i],"-s")==0)
   {
     n1 = atoi(argv[i+1]);
     n2 = atoi(argv[i+2]);
     n3 = atoi(argv[i+3]);
   }
   if (strcmp(argv[i],"-c")==0)
   {
     strcpy(comment,argv[i+1]);
   }

fold=0;
if (n1==0 && n2==0 && n3==0) fold=1;
if (fold){n1 = 1; n2 = 1; n3 = 1;}

nsp=20;
filein=malloc(LEN1*sizeof(char));
tmpstr=malloc(LEN*sizeof(char));

strcpy(filein,argv[1]);

atm=malloc(nsp*sizeof(char *));
for (i=0;i<nsp;i++)
  {
    atm[i]=malloc(3*sizeof(char));
    //atm[i]=NULL;
  }
 
FPin=fopen(filein,"r");
if (!FPin)
  { 
    printf("File error or doesn't exist\n");
    free(filein);
    return;
  }

fgets(tmpstr,LEN,FPin);
sscanf(tmpstr," %d  %c ",&natm,&ch);
 
if (ch=='S')
{
  fgets(tmpstr,LEN,FPin);
  for(i=1;i<=natm;i++)
  {
    fgets(tmpstr,LEN,FPin);
  }
  fgets(tmpstr,LEN,FPin);
  sscanf(tmpstr," %lg %lg %lg",&(or[0]),&(or[1]),&(or[2]));
  fgets(tmpstr,LEN,FPin);
  sscanf(tmpstr," %lg %lg %lg",&(v1[0]),&(v1[1]),&(v1[2]));
  fgets(tmpstr,LEN,FPin);
  sscanf(tmpstr," %lg %lg %lg",&(v2[0]),&(v2[1]),&(v2[2]));
  fgets(tmpstr,LEN,FPin);
  sscanf(tmpstr," %lg %lg %lg",&(v3[0]),&(v3[1]),&(v3[2]));
  //printf(" %f %f %f \n",or[0],or[1],or[2]);
  //printf(" %f %f %f \n",v1[0],v1[1],v1[2]);
  //printf(" %f %f %f \n",v2[0],v2[1],v2[2]);
  //printf(" %f %f %f \n",v3[0],v3[1],v3[2]);
  rewind(FPin);
  fgets(tmpstr,LEN,FPin); //first line
}

fgets(tmpstr,LEN,FPin); //atoms line
adv = strcspn(tmpstr,"\n");
tmpstr[adv]=' ';

adv=0;
i=0;
while(sscanf(tmpstr+adv," %2c",atm[i]) > 0)
{
   adv=(int) ( strstr(tmpstr,atm[i])- tmpstr) +2;
   i++;
}

printf("%d\n",natm*n1*n2*n3);

nsp=i;
printf("%s \n",comment);

for(l3=0; l3<n3; l3++)
{
  for(l2=0; l2<n2; l2++)
  {
    for(l1=0; l1<n1; l1++)
    { 
      rewind(FPin);
      fgets(tmpstr,LEN,FPin); //first line
      fgets(tmpstr,LEN,FPin); //atoms line
      for(i=1;i<=natm;i++)
      {
        fgets(tmpstr,LEN,FPin);
        sscanf(tmpstr," %d %d %lg %lg %lg",&j,&ind,&x,&y,&z);
        if(ind>nsp){printf("Error in file: at row %d, type %d > nsp %d\n",j,ind,nsp);return;}
       
        xl = x +  v1[0]*l1 + v2[0]*l2 + v3[0]*l3;
        yl = y +  v1[1]*l1 + v2[1]*l2 + v3[1]*l3;
        zl = z +  v1[2]*l1 + v2[2]*l2 + v3[2]*l3;

        if (fold) 
        {  
          p[0] = v1[0] + v2[0] + v3[0];
          p[1] = v1[1] + v2[1] + v3[1];
          p[2] = v1[2] + v2[2] + v3[2];
          if (xl > or[0] +  p[0]){ xl -= p[0]; } 
          else if (xl < or[0]){   xl += p[0]; }
          if (yl > or[1] +  p[1]){  yl -= p[1]; }
          else if (yl < or[1]){  yl += p[1]; }
          if (zl > or[2] +  p[2]){ zl -= p[2]; }
          else if (zl < or[2]){ zl += p[2]; }
        }
        printf("%s %f %f %f \n",atm[ind-1],xl,yl,zl); 
      }
    }
  }
}

fclose(FPin);

free(filein);
free(tmpstr);
free(comment);
for(i=0;i<nsp;i++){ free(atm[i]); }
free(atm);

}
	 	
    /* printf("HETATM%5d%3s%12d     %6.3f  %6.3f  %6.3f\n",i,atm[ind],i,x,y,z);*/	 	
