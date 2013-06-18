#include<model_correlation.h>

int main(int argc, char *argv[])
{

 
  
  
  int model, survey,i;
  
  double correlation[THETA_PARTITIONS];
  double correlation_s[THETA_PARTITIONS];
 double correlation_p[THETA_PARTITIONS];
  char filename[100];
  double theta;
  double xc,yc,xr,yr;
  double *x_cat,*y_cat,*x_rand,*y_rand;
  double histogram_cat[THETA_PARTITIONS],histogram_rand[THETA_PARTITIONS],histogram_dr[THETA_PARTITIONS];
  double std_rand[THETA_PARTITIONS];

  FILE *file;
  int ch;
  int n_laes;
  
  for(model=0;model<N_MODELS;model++)
    {
      for(survey=0;survey<MOCK_SURVEYS;survey++)
	{
	  
	  sprintf(filename,"/home/jemejia/CosmicVarianceLAES/data/laes/mock_cat/model_%d_mock_%d.txt",model,survey);
	  file=fopen(filename,"r");
	  
	  n_laes=0;
	  while(!feof(file))
	    {
	      fscanf(file,"%lf %lf %lf %lf",&xc,&yc,&xr,&yr);
	      n_laes++;
	    }
	  
	  fclose(file);
	  
  
	  x_cat=malloc(n_laes*sizeof(double));
	  y_cat=malloc(n_laes*sizeof(double));
	  x_rand=malloc(n_laes*sizeof(double));
	  y_rand=malloc(n_laes*sizeof(double));
	  
	  file=fopen(filename,"r");
	  for(i=0;i<n_laes;i++)
	    {
	      fscanf(file,"%lf %lf %lf %lf",&xc,&yc,&xr,&yr);
	      x_cat[i]=xc;
	      y_cat[i]=yc;
	      x_rand[i]=xr;
	      y_rand[i]=yr;	       
	      
	    }
	  printf("nlaes=%d\n",n_laes);
	  fclose(file);
	  printf("Computing RR \n");
	  RR_histogram(histogram_rand,n_laes,x_rand,y_rand,std_rand);
	  printf(" RR Computed model %d survey %d \n",model,survey);
	  //	  for(i=0;i<THETA_PARTITIONS;i++)
	  //{
	  //  printf("hist[%d]=%lf std[%d]=%lf  \n",i,histogram_rand[i],i,std_rand[i]);
	  //  }
	  printf("Computing DD \n");
	  DD_histogram(n_laes,x_cat,y_cat,histogram_cat);
	  printf(" DD Computed \n");
	  printf("Computing DR \n");
	  DR_histogram(n_laes,(int)N_CATS_DR, x_cat, y_cat, x_rand,y_rand,histogram_dr);
	  printf(" DR Computed  \n");
	  
	  printf("computing Correlation funtion  \n");
	  correlation_landy(histogram_cat,histogram_rand,histogram_dr,correlation);
	  //correlation_peebles(histogram_cat,histogram_dr,correlation_p);
	  //correlation_standard(histogram_cat,histogram_rand,correlation_s);
	  printf("Correlation funtion computed \n");
	  sprintf(filename,"/home/jemejia/CosmicVarianceLAES/data/laes/correlation/model_%d_mock_%d.txt",model,survey);
	  file=fopen(filename,"w");
	  for(i=0;i<THETA_PARTITIONS;i++)
	    {
	      theta=THETA_MIN + (THETA_MAX - THETA_MIN)*i/(1.0*THETA_PARTITIONS);
	      fprintf(file,"%lf %lf\n",theta,correlation[i]);
	      printf("%lf %lf\n",theta,correlation[i]);
	      
	      
	    }
	  fclose(file);
	  free(x_cat);
	  free(y_cat);
	  free(x_rand);
	  free(y_rand);
	}
    }		 
  
  
  return 0;
}

