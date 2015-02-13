#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LEN 250
#define MAXPNTS 3000
#define ANG 1e-10


int main(int argc,char *argv[])
{

  FILE *FPin;
  char *filein, *tmpstr, *name1;
  int k,p1,p2,p3;
  int *pnts;
  double *x,*y,*z, xx,yy,zz,conv;
  double r1[3],r2[3],d1,d2,r[3],d;
  int *line1,*line2;

  int echo;

 
echo = 0; 

if (argc<3)
{ 
  printf("Usage: drawlines unrolled.geo units(meters)\n");
  return;
}

filein=malloc(60*sizeof(char));
tmpstr=malloc(LEN*sizeof(char));
name1=malloc(LEN*sizeof(char)); 

 pnts=malloc(MAXPNTS*sizeof(int));
 line1=malloc(MAXPNTS*sizeof(int));
 line2=malloc(MAXPNTS*sizeof(int));
 x=malloc(MAXPNTS*sizeof(double));
 y=malloc(MAXPNTS*sizeof(double));
 z=malloc(MAXPNTS*sizeof(double));


strcpy(filein,argv[1]);

 conv = atof(argv[2])/ANG;

FPin=fopen(filein,"r");
if (!FPin)
  { 
    printf("File error or doesn't exist\n");
    free(filein);
    return;
  }

	while (feof(FPin)==0)
        {

	        name1[0]='\0';
		fgets(tmpstr,LEN,FPin);

		//if(echo){printf("\n%s",tmpstr);}

                sscanf(tmpstr," %5c %d",name1,&k);
		
		if(strcmp(name1,"Point") == 0) 
		{
		  
		  if(echo){printf("\n%s%d) = ",name1,k);}
		  
		  sscanf(tmpstr," Point( %d ) = { %lg , %lg , %lg",&k,&xx,&yy,&zz);		  
		  
		  x[k] = xx*conv; y[k] = yy*conv; z[k] = zz*conv;

		  if(echo){printf("%f %f %f \n",x[k],y[k],z[k]);}

		}
		
		if(strcmp(name1,"Line(") == 0) 
		{
		
		  if(echo){printf("\n%s%d) = ",name1,k);}		  
		  
		  sscanf(tmpstr," Line( %d ) = { %d , %d ",&k,&p1,&p2);

		  printf("draw l%d {%f %f %f}{%f %f %f}\n",k,x[p1],y[p1],z[p1],x[p2],y[p2],z[p2]);

	
		}

		if(strcmp(name1,"Circl") == 0) 
		{
		
		  if(echo){printf("\n%s%d) = ",name1,k);}		  
		  
		  sscanf(tmpstr," Circle( %d ) = { %d , %d , %d ",&k,&p1,&p2,&p3);

      r1[0] = x[p1] - x[p2];
      r1[1] = y[p1] - y[p2];
      r1[2] = z[p1] - z[p2];
      r2[0] = x[p3] - x[p2];
      r2[1] = y[p3] - y[p2];
      r2[2] = z[p3] - z[p2];


      r[0] = r1[0]+r2[0];
      r[1] = r1[1]+r2[1];
      r[2] = r1[2]+r2[2];

      d1 = sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]);
      d2 = sqrt(r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2]);
      d = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

      r[0] = x[p2] + (d1+d2)/2.0 * r[0]/d;
      r[1] = y[p2] + (d1+d2)/2.0 * r[1]/d;
      r[2] = z[p2] + (d1+d2)/2.0 * r[2]/d;

		  printf("draw l%d curve {%f %f %f}{%f %f %f}{%f %f %f}\n",k,x[p1],y[p1],z[p1], 
			 r[0],r[1],r[2], x[p3],y[p3],z[p3]);

	
		}


	}

   

}
