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
  int *line1,*line2;
  double v1,v2,v3,diameter,angle;	

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

	          v1 = (y[p1]*z[p3] - y[p3]*z[p1]);
	          v2 = (z[p1]*x[p3] - z[p3]*x[p1]);
	          v3 = (x[p1]*y[p3] - x[p3]*y[p1]);
		  v1 = v1/sqrt(v1*v1+v2*v2+v3*v3);
		  v2 = v2/sqrt(v1*v1+v2*v2+v3*v3);
		  v3 = v3/sqrt(v1*v1+v2*v2+v3*v3);
		  diameter = 2*sqrt((x[p1]-x[p2])*(x[p1]-x[p2])+(y[p1]-y[p2])*(y[p1]-y[p2])+(z[p1]-z[p2])*(z[p1]-z[p2]));	
                  angle = acos(x[p1]*x[p3]+y[p1]*y[p3]+z[p1]*z[p3]);

		  printf("draw l%d arc {%f %f %f}{%f %f %f}{%f %f %f}{0 %f 0} SCALE %f\n",k,x[p2],y[p2],z[p2], 
			 v1,v2,v3, x[p1],y[p1],z[p1],angle*180/3.141593,diameter);

	
		}


	}

   

}
