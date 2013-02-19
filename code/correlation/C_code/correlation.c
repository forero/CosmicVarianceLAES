#include<mass_min_max.h>

int DD_histogram(int N_cat,double X_cat[],double Y_cat[],double Z_cat[],double histogram_cat[])
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
      
      for(j=0;j<i;j++)
	{
	  d1=sqrt(  (X_cat[i]-X_cat[j])*(X_cat[i]-X_cat[j]) + (Y_cat[i]-Y_cat[j])*(Y_cat[i]-Y_cat[j]) );	      
	  theta=206265*d1/(D_z);
          theta=log10(theta); 
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



int RR_histogram(double histogram_rand[],int N_POINTS_RANDOM,double std[])
{
  
  //  double X[N_POINTS_RANDOM*N_CATS_RANDOM],Y[N_POINTS_RANDOM*N_CATS_RANDOM],Z[N_POINTS_RANDOM*N_CATS_RANDOM];
  double *X,*Y,*Z;
  X=malloc( N_POINTS_RANDOM*N_CATS_RANDOM*sizeof(double) );
  Y=malloc( N_POINTS_RANDOM*N_CATS_RANDOM*sizeof(double) );
  Z=malloc( N_POINTS_RANDOM*N_CATS_RANDOM*sizeof(double) );
  double *hist_all;
  hist_all=malloc( N_CATS_RANDOM*THETA_PARTITIONS*sizeof(double) );
  double theta;
  double d1;
  gsl_histogram *hist_theta_rand;
  
  int i,j,k,w,m;
  double scale_rand=2.0/(1.0*N_POINTS_RANDOM*(N_POINTS_RANDOM-1));
  
  
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
  
  for(m=0;m<N_POINTS_RANDOM*N_CATS_RANDOM;m++)
    {
      
      X[m]= 0.0 + SUB_BOX_SIZE_X*gsl_rng_uniform (generador);
      Y[m]= 0.0 + SUB_BOX_SIZE_Y*gsl_rng_uniform (generador);
      Z[m]= 0.0 + SUB_BOX_SIZE_Z*gsl_rng_uniform (generador);
      
    }
 
  for(k=0;k<N_CATS_RANDOM;k++)
    {
      
      for(i=0;i<N_POINTS_RANDOM;i++)
	{
	  
	  for(j=0;j<i;j++)
	    {	
	      d1=sqrt(  ( X[i+N_POINTS_RANDOM*k]-X[j+N_POINTS_RANDOM*k] )*( X[i+N_POINTS_RANDOM*k]-X[j+N_POINTS_RANDOM*k] ) + ( Y[i+N_POINTS_RANDOM*k]-Y[j+N_POINTS_RANDOM*k] )*( Y[i+N_POINTS_RANDOM*k]-Y[j+N_POINTS_RANDOM*k] ));	      
	      theta=206265*d1/(D_z);
	      
	      theta=log10(theta); 
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
  
  for(i=0;i<THETA_PARTITIONS;i++)
    {
      histogram_rand[i]/=N_CATS_RANDOM;      
    }
  
  for(k=0;k<N_CATS_RANDOM;k++)
    {
      for(w=0;w<THETA_PARTITIONS;w++)
	{
	  std[w]+=(histogram_rand[w]-hist_all[k+N_CATS_RANDOM*w])*(histogram_rand[w]-hist_all[k+N_CATS_RANDOM*w]);
	}
      
    }
  
  for(w=0;w<THETA_PARTITIONS;w++)
    {
      std[w]/=N_CATS_RANDOM;
    }
  
  return 0;      
}

int DR_histogram(int N_cat,int p,int q,int r,int N_cats_rand,double X_cat[],double Y_cat[],double Z_cat[],double histogram_dr[])
{
  double X[N_cat],Y[N_cat],Z[N_cat];
  
  
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
  

  //for(j=0;j<N_cats_rand;j++)
  //{
      for(i=0;i<N_cat;i++)
	{ 
	  X[i]= p*SUB_BOX_SIZE_X +  SUB_BOX_SIZE_X*(p+1)*gsl_rng_uniform (generador);
	  Y[i]= q*SUB_BOX_SIZE_Y + SUB_BOX_SIZE_Y*(q+1)*gsl_rng_uniform (generador);
	  Z[i]= r*SUB_BOX_SIZE_Z + SUB_BOX_SIZE_Z*(r+1)*gsl_rng_uniform (generador);
	} 
      
      
      for(i=0;i<N_cat;i++)
	{
	  for(k=0;k<N_cat;k++)
	    {
	      d1=sqrt(  (X[i]-X_cat[k])*(X[i]-X_cat[k]) + (Y[i]-Y_cat[k])*(Y[i]-Y_cat[k]) );
	      theta=206265*d1/(D_z);
	     
	      theta=log10(theta); 
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
