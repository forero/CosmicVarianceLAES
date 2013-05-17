#include<model_correlation.h>

int DD_histogram(int N_cat,double X_cat[],double Y_cat[],double histogram_cat[])
{
  
  double theta;
  double d1;
  gsl_histogram *hist_theta_cat;
  
  int i,j,w;
  double scale_cat=2.0/(1.0*N_cat*(N_cat-1));
  
  
  hist_theta_cat=gsl_histogram_alloc(THETA_PARTITIONS);
  gsl_histogram_set_ranges_uniform(hist_theta_cat,THETA_MIN,THETA_MAX);
  
  for(i=0;i<N_cat;i++)
    {
      
      for(j=i+1;j<N_cat;j++)
	{
	  d1=sqrt(  (X_cat[i]-X_cat[j])*(X_cat[i]-X_cat[j]) + (Y_cat[i]-Y_cat[j])*(Y_cat[i]-Y_cat[j]) );	      
	  theta=206265*d1/(D_z);
          //theta=log10(theta); 
	  gsl_histogram_increment(hist_theta_cat,theta);
	}
    }  
  
  gsl_histogram_scale (hist_theta_cat,scale_cat);
  
  for(w=0;w<THETA_PARTITIONS;w++)
    {
      histogram_cat[w]=gsl_histogram_get(hist_theta_cat, w);
    }
  
  return 0;
  
}



int RR_histogram(double histogram_rand[],int random_points,double X[],double Y[],double std[])
{
  
  
  double *hist_all;
  hist_all=malloc( N_CATS_RANDOM*THETA_PARTITIONS*sizeof(double) );
  double theta;
  double d1;
  gsl_histogram *hist_theta_rand;
  
  int i,j,k,w,m;
  double scale_rand=2.0/(1.0*random_points*(random_points-1));
  
  
  hist_theta_rand=gsl_histogram_alloc(THETA_PARTITIONS);
  gsl_histogram_set_ranges_uniform(hist_theta_rand,THETA_MIN,THETA_MAX);
  
  //1-DECLARA
  const gsl_rng_type *tipo;
  gsl_rng *generador;
  
  //2-ALOCA E INICIALIZA
  tipo=gsl_rng_mt19937;
  generador=gsl_rng_alloc(tipo);
  gsl_rng_set(generador,time(NULL));
  
  
  for(i=0;i<THETA_PARTITIONS;i++)
    {
      histogram_rand[i]=0;
      std[i]=0;
    }
  
  
 
  for(k=0;k<N_CATS_RANDOM;k++)
    {
      
      for(i=0;i<random_points;i++)
	{
	  
	  for(j=i+1;j<random_points;j++)
	    {	
	      d1=sqrt(  ( X[i+random_points*k]-X[j+random_points*k] )*( X[i+random_points*k]-X[j+random_points*k] ) + ( Y[i+random_points*k]-Y[j+random_points*k] )*( Y[i+random_points*k]-Y[j+random_points*k] ));	      
	      theta=206265*d1/(D_z);
	      
	      //theta=log10(theta); 
	      gsl_histogram_increment(hist_theta_rand,theta);
	    } 
	}
      gsl_histogram_scale(hist_theta_rand,scale_rand);
      for(w=0;w<THETA_PARTITIONS;w++)
	{
	  histogram_rand[w]+=gsl_histogram_get(hist_theta_rand, w);
	  hist_all[k+N_CATS_RANDOM*w]=gsl_histogram_get(hist_theta_rand, w);
	}
      
      gsl_histogram_reset (hist_theta_rand);
    }
  
  
  for(k=0;k<N_CATS_RANDOM;k++)
    {
      for(w=0;w<THETA_PARTITIONS;w++)
	{
	  std[w]+= ( histogram_rand[w] - hist_all[k+N_CATS_RANDOM*w] )*( histogram_rand[w] - hist_all[k+N_CATS_RANDOM*w] ) ;
	}
      
    }
  
  for(w=0;w<THETA_PARTITIONS;w++)
    {
      histogram_rand[w]/=N_CATS_RANDOM;      
      std[w]/=N_CATS_RANDOM;
      std[w]=sqrt( std[w] );
    }
  
  return 0;      
}

int DR_histogram(int N_cat,int N_cats_rand,double X_cat[],double Y_cat[],double X[],double Y[],double histogram_dr[])
{
  
  
  
  double theta;
  double d1;
  gsl_histogram *hist_theta_dr;
  
  int i,j,k;
  double scale_dr=N_cat*N_cat;
  
  
  hist_theta_dr=gsl_histogram_alloc(THETA_PARTITIONS);
  gsl_histogram_set_ranges_uniform(hist_theta_dr,THETA_MIN,THETA_MAX);
  
  //1-DECLARA
  const gsl_rng_type *tipo;
  gsl_rng *generador;
  
  //2-ALOCA E INICIALIZA
  tipo=gsl_rng_mt19937;
  generador=gsl_rng_alloc(tipo);
  gsl_rng_set(generador,time(NULL));

  
  for(i=0;i<THETA_PARTITIONS;i++)
    {
      histogram_dr[i]=0;
    }
  

  
      
      
      for(i=0;i<N_cat;i++)
	{
	  for(k=0;k<N_cat;k++)
	    {
	      d1=sqrt(  (X[i]-X_cat[k])*(X[i]-X_cat[k]) + (Y[i]-Y_cat[k])*(Y[i]-Y_cat[k]) );
	      theta=206265*d1/(D_z);
	     
	      //theta=log10(theta); 
	      gsl_histogram_increment(hist_theta_dr,theta);
	    }
	}
      
      for(i=0;i<THETA_PARTITIONS;i++)
	{
	  //histogram_dr[i]+=gsl_histogram_get(hist_theta_dr, i);
	  histogram_dr[i]=gsl_histogram_get(hist_theta_dr, i);
	}
      gsl_histogram_reset(hist_theta_dr);
      //   }
  
  for(i=0;i<THETA_PARTITIONS;i++)
	{
	  // histogram_dr[i]/=((double)(N_cats_rand*scale_dr));	 
	  histogram_dr[i]/=((double)(scale_dr));	 
	}
  
  return 0;
  
}

int correlation_landy(double histogram_cat[],double histogram_rand[], double histogram_dr[],double correlation[])
{
  int i;
  for(i=0;i<THETA_PARTITIONS;i++)
    {
      correlation[i]=(histogram_cat[i]-2.0*histogram_dr[i]+histogram_rand[i])/histogram_rand[i];
    }
  
   return 0;
}

int correlation_standard(double histogram_cat[],double histogram_rand[],double correlation[])
{
  int i;
  for(i=0;i<THETA_PARTITIONS;i++)
    {
      correlation[i]=histogram_cat[i]/histogram_rand[i]-1.0;
    }
  
   return 0;
}

int correlation_peebles(double histogram_cat[],double histogram_dr[],double correlation[])
{
 int i;
  for(i=0;i<THETA_PARTITIONS;i++)
    {
      correlation[i]=histogram_cat[i]/histogram_dr[i]-1.0;
    }
  
   return 0;
}
