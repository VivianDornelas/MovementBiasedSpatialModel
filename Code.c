/* 29/07/22                                                                 	*/
/* Critical habitat (Home range + diffusion)         			        */
/* Equation: d_t u = ru (1-u/K) + tau^(-1)d_x [(x-lambda)u]  + D d_xx u 	*/
/*                          By Vivian Dornelas ;) 	    			*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//***********************  funtions  *********************************
int alocar();
int libera();
void evolution();
int seed();
double ranf();
int int_ranf(int a,int b);
void medidas();
//********************************************************************


//********************  Global variables  ****************************
int i, j, x, t, tmax, dxs, dts, disc, iseed, flag, cont;
double dt, Ll, Lr, L, dx, dt, dx2, x0, g, K, r, Diff, u0;
double *ubefore, *ucurrent, d1u, d2u; 
double rand1;
double mu1, mu, mx, sigma2x, umax, xumax, dist; 
//********************************************************************

FILE *out1;

int main(int argc, char *argv[]){
	if(argc != 10){
		printf("Error de sintaxis: %s >\n",argv[0]);exit(1);}
	if(sscanf(argv[1],"%lf",&Ll) != 1){
		printf("Erro no parametro A\n");exit(1);}
	if(sscanf(argv[2],"%lf",&Lr) != 1){
		printf("Erro no parametro A\n");exit(1);}
	if(sscanf(argv[3],"%d",&tmax) != 1){
		printf("Erro no parametro A\n");exit(1);}
	if(sscanf(argv[4],"%lf",&g) != 1){
		printf("Erro no parametro A\n");exit(1);}        
	if(sscanf(argv[5],"%d",&dxs) != 1){
		printf("Erro no parametro A\n");exit(1);}
	if(sscanf(argv[6],"%d",&dts) != 1){
		printf("Erro no parametro A\n");exit(1);}
	if(sscanf(argv[7],"%lf",&K) != 1){
		printf("Erro no parametro A\n");exit(1);}
	if(sscanf(argv[8],"%lf",&r) != 1){
		printf("Erro no parametro A\n");exit(1);}
	if(sscanf(argv[9],"%lf",&Diff) != 1){
		printf("Erro no parametro A\n");exit(1);}

	//  ***** inicializa variaveis *****
	char dest1[256];
	iseed=seed();
	dx=1.0/pow(10,dxs);
	L=Lr+Ll;
	x=L/dx;
	dx2=dx*dx;
	x0=Ll;
	disc=pow(10,dts);
	dt=1.0/disc;
	u0=0.0;
	if ((2.0*Diff*dt/(dx2)) >= 1.0){
		printf ("Errou feio, errou rude!\n\n");
		exit(0);
	}
	if ((g*dt/dx)>= 1.0){
		printf ("Errou feio, errou rude!\n\n");
		exit(0);
	}

	alocar();

	//  ************* c.i. *************
	ubefore[0]=0;		
	ubefore[x]=0;	
	ucurrent[0]=0;		
	ucurrent[x]=0;
	for (i=1;i<x;i++){
		rand1=ranf();
		rand1/=pow(10,4);     
		ubefore[i]=u0+rand1;
	}
	sprintf(dest1, "meas-Ll%.2lf_Lr%.2lf_r%.2lf_K%.1lf_g%.2lf_D%.4lf_t%d_dxe-%d_dte-%d.txt", Ll, Lr, r, K, g, Diff, tmax, dxs, dts);
	out1 = fopen( dest1, "w");
	fprintf(out1,"#t, g, umax, xumax, mx, mu, sigma2x\n"); 

	mu1=0.0;
	mu=0.0;

	for (t=0; t < tmax+1; t++){
	        for (j=0;j<disc;j++)
			    evolution();	
		if ((t%100)==0){
			medidas();
			fprintf(out1,"%d\t%.4lf\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n", t, g, umax, xumax, mx, mu, sigma2x);
			if ( (mu == mu1) || ( mu < pow(10,-100) ) )
				break;
			mu1=mu;
		}
	}
	fclose(out1);

	medidas();
	printf("%.3lf\t%.3lf\t%.17g\t%.17g\t%d\t%.17g\t%.17g\t%.17g\n", Ll, Lr, umax, xumax, tmax, mx, mu, sigma2x);

	sprintf(dest1, "pop-Ll%.2lf_Lr%.2lf_r%.2lf_K%.1lf_g%.2lf_D%.4lf_t%d_dxe-%d_dte-%d.txt", Ll, Lr, r, K, g, Diff, tmax, dxs, dts);
	out1 = fopen( dest1, "w");
	fprintf(out1,"#x\tu\n"); 
	for (i=0;i<x+1;i++)
        	fprintf(out1,"%.3lf\t%.17g\n", (i*dx-x0), ubefore[i]);
	fclose(out1);
    
	libera();
	return 0;
}    

//********************************************************************

void evolution(){
	ucurrent[0]=0;
	ucurrent[x]=0;
	for (i=1;i<x;i++){
		d1u= ( ubefore[i+1] - ubefore[i-1])/(2.0*dx);
		d2u= ( ubefore[i+1] + ubefore[i-1] - 2*ubefore[i])/dx2;
		ucurrent[i]= ubefore[i] + dt * ( r*ubefore[i] *( 1.0- (ubefore[i]/K)) + g* (ubefore[i] + ( (i*dx)-x0  )*d1u) + Diff*d2u );
		if (ucurrent[i]<0)
			ucurrent[i]=0;
	}
	for (i=1;i<x;i++)
    		ubefore[i]=ucurrent[i];	
}

void medidas(){
	umax = 0.0;
	mu = 0.0;
	mx = 0.0;
	sigma2x = 0.0;
	for (i=0;i<x+1;i++){
		if (ucurrent[i]>umax){
			umax = ucurrent[i];
			xumax= i*dx-x0;
		}
		mu += ucurrent[i];
		mx += (i*dx-x0) * ucurrent[i];
	}
	mx/=mu;
	for (i=0;i<x+1;i++){
		sigma2x += ( pow((i*dx-x0 - mx), 2.0)*ucurrent[i] );
	}
	sigma2x /= mu;
}


int alocar(){
	if(((ubefore = (double*)calloc(x+1,sizeof(double))) == NULL)){
		printf("Sem memoria\n"); exit(1);}
	if(((ucurrent = (double*)calloc(x+1,sizeof(double))) == NULL)){
		printf("Sem memoria\n"); exit(1);}
	return 0;
}

int libera(){
	free(ubefore);
	free(current);
	return 0;
}

// Initial seed from the system time and forced to be odd 
int seed(){
	int t1,t2,t3,t4,t5,t6;
	time_t tp;
	struct tm *t;
	time ( &tp);  
	t  = localtime(&tp);
	t1 = t->tm_sec;
	t2 = t->tm_min;
	t3 = t->tm_hour;
	t4 = t->tm_mday;
	t5 = t->tm_mon;
	t6 = t->tm_year;
	iseed = t6+70*(t5+12*(t4+31*(t3+23*(t2+59*t1))));
	if ((iseed%2) == 0)
		iseed--;
	return iseed;
}

// Random number generator between 0 and 1
double ranf(){
	const int ia=16807,ic=2147483647,iq=127773,ir=2836;
	int il,ih,it;
	double rc;
	ih = iseed/iq;
	il = iseed%iq;
	it = ia*il-ir*ih;
	if (it > 0)
		iseed = it;
	else
		iseed = ic+it;
	rc = ic;
	return iseed/rc;
}

// Generates integer random numbers between a and b, including them 
int int_ranf(int a,int b){
	return a+(b-a)*ranf()+0.5; 
}
