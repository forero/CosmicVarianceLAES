#include<mass_min_max.h>
/*This routine assign LAEs to halos in one of two differents models.

  The first one fixes the mininal mass (double *mass_min) for an halo to host a LAEs is fixed and maximal mass is obtained in ascendant order till the observational  number  of LAEs is reproduced (int N_laes). The mass of the last halo is the maximal mass(double *mass_max)  of the purposed model. The halos that within that interval (m_min, m_max) that host LAEs are chosen in an aleatory way. The fraction between the number of LAEs and the number of Halos is also measured (double *duty_cycle).
   

The second model fixes both the minimal mass (double *mass_min) and the maximal mass (double *mass_max).  The fraction between the observational number of LAEs (int N_laes) and the number of Halos between mass_min and mass_max  is  measured (double *duty_cycle). LAEs are assigned randomly  to halos with mass within that interval.

To choose the first model please set int model to 0 and 1 for the second one.
   

*/



/*Numero de objetos sin leer en el subcatalogo:  N*/
int LAE_catalog_generator(int model,int N_laes,double mass[],double x[],double y[],double z[],double *mass_min,double halo_m_dev,int *N,int *T,double *mass_max,double *duty_cycle,double x_cat[],double y_cat[],double z_cat[],int *catalog)
{
   //1-DECLARA
   const gsl_rng_type *tipo;
  gsl_rng *generador;
  
  //2-ALOCA E INICIALIZA
  tipo=gsl_rng_mt19937;
  generador=gsl_rng_alloc(tipo);
  gsl_rng_set(generador,time(NULL));   
  int rand;
  int k,l;
  double exp_mass_min;
  double mass1;
  l=1;
  while( mass[*N-l]<=*mass_min)
    {
      l++;
    }
  
  *mass_min=mass[*N-l];
  exp_mass_min=pow(10,*mass_min);
  *N=*N-l;

  *T=*N;		
  //printf("N=%d \n",*N);
  if(*N-N_laes>0)
    { 
      // printf("N-N_LAES>0, m_min=%lf m_max=%lf \n",*mass_min,*mass_max);
      *catalog=1;
      if(model==0)
	{
	  printf("model=0 \n");
	  if( ( pow(10,mass[*N-N_laes]) - exp_mass_min  ) <= halo_m_dev )
	    {	  
	      *mass_max=log10( exp_mass_min + halo_m_dev );	  
	    }
	  else
	    {
	      *mass_max=mass[*N-N_laes]; 
	    }
	  l=1;			  
	  
	  while( mass[*T-l]<=(*mass_max) )
	    {
	      l++;
	    }
	  
	  for(k=0;k<N_laes;k++)
	    {
	      rand=gsl_rng_uniform_int(generador, l );
	      x_cat[k]=x[*N-rand];
	      y_cat[k]=y[*N-rand];
	      z_cat[k]=z[*N-rand];
	    }
	  
	  *duty_cycle=((double)N_laes)/((double)l);
	}
      if(model==1)
	{
	  //printf("model=1 \n");
	  l=1;
	  *duty_cycle=2;	
	  //k=0;
	  while( mass[*T-l]<=(*mass_max) )
	    {
	      l++;
	    }
	  
	  *duty_cycle=((double)N_laes)/((double)l);	  
	  if(*duty_cycle>1.0) *catalog=0;
	  //*mass_max+=k*0.1;
	  //k++;
	}
      //printf("l=%d \n",l);
      
      for(k=0;k<N_laes;k++)
	    {
	      rand=gsl_rng_uniform_int(generador, l );
	      //printf("x[N-rand]=%lf y[N-rand]=%lf z[N-rand]=%lf \n",x[*N-rand],y[*N-rand],z[*N-rand]);
	      //printf("rand=%d \n",rand);
	      x_cat[k]=x[*N-rand];
	      //printf("x[N-l]=%lf  \n",x[*N-rand]);
	      y_cat[k]=y[*N-rand];
	      //printf("y[N-l]=%lf  \n",y[*N-rand]);
	      z_cat[k]=z[*N-rand];
	      //printf("x[N-l]=%lf y[N-l]=%lf z[N-l]=%lf \n",x[*N-rand],y[*N-rand],z[*N-rand]);
	    }
      // printf("x[N-l]=%lf y[N-l]=%lf z[N-l]=%lf \n",x[*N-l],y[*N-l],z[*N-l]);
      *duty_cycle=((double)N_laes)/((double)l);	  
	  //printf("duty_cycle=%lf \n",*duty_cycle);
    
    }else{*catalog=0;}
  //  printf("saliendo \n");
  return 0;
}
