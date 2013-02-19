#include<mass_min_max.h>
/*
Ingresar a=0 si z=5.7 ó a=1 para otro redshift
observational_filename es la ruta del archivo con la funcion de correlacion para cada uno de los redshift
*/

int correlation_statistic(double correlation_data[],char observational_file[],int a,int catalog[])
{
  double correlation1;
  double *correlation_mean;
  double *correlation_std;
  double *correlation_data_obs;
  double *correlation_std_obs; 
  double *chisquare_all_models;
  double *chisquare_mean_models;

  int x_div,y_div,z_div;
  x_div=floor(BOXSIZE/SUB_BOX_SIZE_X);
  y_div=floor(BOXSIZE/SUB_BOX_SIZE_Y);
  z_div=floor(BOXSIZE/SUB_BOX_SIZE_Z);

  
  double theta[THETA_PARTITIONS];
  int s,m,i,K,g;
  char filename[300];
  FILE *file3;
  double trans;
  K=N_STEPS;
  
  
	 
  ordering order_mws[K*(x_div*y_div*z_div)];
 
  for(m=0;m<K;m++)
    { 
      
      for(s=0;s<(x_div*y_div*z_div);s++)
	{
	  g=m+K*s;
	  
	  order_mws[g].p=m;
	  
	  order_mws[g].q=s;
	  
	  
	}
      
    }
  
  
  
  correlation_mean=malloc(K*THETA_PARTITIONS*sizeof(double));
  correlation_std=malloc(K*THETA_PARTITIONS*sizeof(double));
  
  correlation_data_obs=malloc(THETA_PARTITIONS*sizeof(double));
  correlation_std_obs=malloc(THETA_PARTITIONS*sizeof(double));
  
  chisquare_all_models=malloc((x_div*y_div*z_div)*K*sizeof(double));
  chisquare_mean_models=malloc(K*sizeof(double));
  
  file3=fopen(observational_file,"r");
  printf("hola_2 \n");
  
  if(a==1)
    {
      for(i=0;i<THETA_PARTITIONS;i++)
	{
	  fscanf(file3,"%lf %lf %lf",&theta[i],&correlation_data_obs[i],&correlation_std_obs[i]); 
	  printf("%lf %lf %lf \n",theta[i],correlation_data_obs[i],correlation_std_obs[i]);
	}
    } 
  else
    {
      for(i=0;i<THETA_PARTITIONS;i++)
	{
	  fscanf(file3,"%lf %lf %lf %lf",&theta[i],&trans,&correlation_data_obs[i],&correlation_std_obs[i]); 
	  printf("%lf %lf %lf \n",theta[i],correlation_data_obs[i],correlation_std_obs[i]);
	}   
    }
  fclose(file3);
  
  

  for(m=0;m<K;m++)
    {
      for(i=0;i<THETA_PARTITIONS;i++)
	{	  
	  correlation_mean[m+K*i]=0;
	  correlation_std[m+K*i]=0;
	  
	}
    }

  
  
  for(m=0;m<K;m++)
    {
      
      chisquare_mean_models[m]=0;
      for(s=0;s<(x_div*y_div*z_div);s++)
	{
	  chisquare_all_models[m+K*s]=0;
	}
    }
  
  
  printf("hola 3 \n");
  int w[K];
  for(m=0;m<K;m++)
    {
      w[m]=0;
      for(s=0;s<(x_div*y_div*z_div);s++)
	{
	  if(catalog[m+N_STEPS*s])
	    {
	      w[m]++;
	      for(i=0;i<THETA_PARTITIONS;i++)
		{
		  correlation1=correlation_data[m+K*(s+(x_div*y_div*z_div)*i)];	      	      
		  correlation_mean[m+K*i]+=correlation1;
		  printf("c[%d_%d_%d]=%lf \n",m,s,i,correlation1);
		}
	    }
	  else
	    {
	      
	      for(i=0;i<THETA_PARTITIONS;i++)
		{
		  correlation_data[m+K*(s+(x_div*y_div*z_div)*i)]=-1000.0;   
		}
	      
	    }
	}
    }
  
  
  
 
  for(m=0;m<K;m++)
    {
      for(i=0;i<THETA_PARTITIONS;i++)
	{	  
	  //  correlation_mean[m+K*i]/=((x_div*y_div*z_div));
	  correlation_mean[m+K*i]=correlation_mean[m+K*i]/(w[m]);
	  printf("c_mean[%d_%d]=%lf \n",m,i,correlation_mean[m+K*i]);
	  
	}
    }
  
  
  for(m=0;m<K;m++)
    {
      for(s=0;s<(x_div*y_div*z_div);s++)
	{
	  if(catalog[m+K*s])
	    {
	      for(i=0;i<THETA_PARTITIONS;i++)
		{   
		  correlation_std[m+K*i]+=( correlation_data[m+K*(s+(x_div*y_div*z_div)*i)]-correlation_mean[m+K*i] )*( correlation_data[m+K*(s+(x_div*y_div*z_div)*i)]-correlation_mean[m+K*i] );
		  
		}
	    }else{continue;}
	}
      
    }
  printf("hola 6 \n");
  FILE *file2[K];
  
  for(m=0;m<K;m++)
    {
      
      sprintf(filename,"./correlation_files/correlation_%d.data",m);
      file2[m]=fopen(filename,"w");
      for(i=0;i<THETA_PARTITIONS;i++)
	{
	  //correlation_std[m+K*i]=sqrt(correlation_std[m+K*i]/(x_div*y_div*z_div));
	  correlation_std[m+K*i]=sqrt( correlation_std[m+K*i]/w[m] );
	  fprintf(file2[m],"%lf %lf %lf \n",theta[i],correlation_mean[m+K*i],correlation_std[m+K*i]);
	  
	}
      
     
    }

  printf("hola 7 \n");
  

  
  
  
  for(m=0;m<K;m++)
    {  
      
      for(i=0;i<THETA_PARTITIONS;i++)
	{	  
	  
	  chisquare_mean_models[m]+=(correlation_data_obs[i]-correlation_mean[m+K*i])*(correlation_data_obs[i]-correlation_mean[m+K*i])/(THETA_PARTITIONS*(correlation_std_obs[i]*correlation_std_obs[i]+correlation_std[m+K*i]*correlation_std[m+K*i]));
	  
	}
      
      
    }
  
  printf("hola 9 \n");
 for(m=0;m<K;m++)
    {
      for(s=0;s<(x_div*y_div*z_div);s++)
	{
	  
	  for(i=0;i<THETA_PARTITIONS;i++)
	    {
	      chisquare_all_models[m+K*s]+=(correlation_data_obs[i]- correlation_data[m+K*(s+(x_div*y_div*z_div)*i)])*(correlation_data_obs[i]- correlation_data[m+K*(s+(x_div*y_div*z_div)*i)])/(THETA_PARTITIONS*correlation_std_obs[i]*correlation_std_obs[i]);
	    }
	}
    }
 

 printf("hola 10 \n");
 FILE *file5,*file6;
 file5=fopen("./correlation_files/chisquare_models.data","w");
 file6=fopen("./correlation_files/chisquare_all.data","w");
 
 
 size_t permutation_mean[K];
 size_t permutation_all[K*(x_div*y_div*z_div)];
 
 printf("hola 11 \n");
 gsl_sort_smallest_index(permutation_mean,K,chisquare_mean_models,1,K);
 gsl_sort_smallest_index(permutation_all,K*(x_div*y_div*z_div),chisquare_all_models,1,K*(x_div*y_div*z_div));
 
 printf("hola 12 \n");
 i=0;
 
 for(m=0;m<K;m++)
   {
     fprintf(file5,"%d %lf %lf %d \n",(int)permutation_mean[m],9.95+0.05*permutation_mean[m],chisquare_mean_models[permutation_mean[m]],w[permutation_mean[m]]);
     
     for(s=0;s<(x_div*y_div*z_div);s++)
       {
	 fprintf(file6,"%d %d  %lf %lf \n",order_mws[permutation_all[i]].p,order_mws[permutation_all[i]].q,9.95+0.05*order_mws[permutation_all[i]].p,chisquare_all_models[permutation_all[i]]);
	 i++;
       }
   }
 printf("hola 13 \n");
 
 return 0;
    
 
}


      
int correlation_statistic_m_min_max(double correlation_data[],char observational_file[],int a,int catalog[])
{
  double correlation1;
  double *correlation_mean;
  double *correlation_std;
  double *correlation_data_obs;
  double *correlation_std_obs; 
  double *chisquare_all_models;
  double *chisquare_mean_models;
  double *chisquare_mean_models_obs;
  
  double theta[THETA_PARTITIONS];
  int s,m,n,i,j,k,K,g;
  char filename[300];
  FILE *file3;
  double trans;
  K=N_STEPS;
  
  int x_div,y_div,z_div;
  x_div=floor(BOXSIZE/SUB_BOX_SIZE_X);
  y_div=floor(BOXSIZE/SUB_BOX_SIZE_Y);
  z_div=floor(BOXSIZE/SUB_BOX_SIZE_Z);

	 
  ordering order_mws[K*(x_div*y_div*z_div)*(K+1)];
  ordering order_mn[K*(K+1)];
 
  for(m=0;m<K;m++)
    { 
      for(n=1;n<K+1-m;n++)
	{
	  for(s=0;s<(x_div*y_div*z_div);s++)
	    {
	      g=m+K*(s+(x_div*y_div*z_div)*n);
	  
	      order_mws[g].p=m;
	      
	      order_mws[g].q=s;
	      
	      order_mws[g].r=n;
	    }
	  
	}
    }
  

for(m=0;m<K;m++)
    { 
      for(n=1;n<K+1-m;n++)
	{	
	  g=m+K*n;
	  order_mn[g].p=m;	      
	  order_mn[g].q=n;
	}
    }

  
  correlation_mean=malloc((K*(K+1)*THETA_PARTITIONS)*sizeof(double));
  correlation_std=malloc((K*(K+1)*THETA_PARTITIONS)*sizeof(double));
  
  correlation_data_obs=malloc(THETA_PARTITIONS*sizeof(double));
  correlation_std_obs=malloc(THETA_PARTITIONS*sizeof(double));
  
  chisquare_all_models=malloc((x_div*y_div*z_div)*(K*(K+1))*sizeof(double));
  chisquare_mean_models=malloc((K*(K+1))*sizeof(double));
 chisquare_mean_models_obs=malloc((K*(K+1))*sizeof(double));
  
  file3=fopen(observational_file,"r");
  printf("hola_2 \n");
  
  if(a==1)
    {
      for(i=0;i<THETA_PARTITIONS;i++)
	{
	  fscanf(file3,"%lf %lf %lf",&theta[i],&correlation_data_obs[i],&correlation_std_obs[i]); 
	  printf("%lf %lf %lf \n",theta[i],correlation_data_obs[i],correlation_std_obs[i]);
	}
    }
  else
    {
      for(i=0;i<THETA_PARTITIONS;i++)
	{
	  fscanf(file3,"%lf %lf %lf %lf",&theta[i],&trans,&correlation_data_obs[i],&correlation_std_obs[i]); 
	  printf("%lf %lf %lf \n",theta[i],correlation_data_obs[i],correlation_std_obs[i]);
	}   
    }
  fclose(file3);



  for(m=0;m<K;m++)
    {
       for(n=1;n<K+1-m;n++)
	{
	  for(i=0;i<THETA_PARTITIONS;i++)
	    {	  
	      correlation_mean[m+K*(i+THETA_PARTITIONS*n)]=0;
	      correlation_std[m+K*(i+THETA_PARTITIONS*n)]=0;
	  
	    }
	}
    }


  for(m=0;m<K;m++)
    {
      for(n=0;n<K+1;n++) 
	{ 
	  //aqui se corre entre 0 y K+1 pues debemos llenar todas las entradas de 0 incluso las que no se van a usar 
	  chisquare_mean_models[m+K*n]=0;
	  chisquare_mean_models_obs[m+K*n]=0;
	  for(s=0;s<(x_div*y_div*z_div);s++)
	    {
	      chisquare_all_models[m+K*(s+x_div*y_div*z_div*n)]=0;
	    }
	}
    }
  
  printf("hola 3 \n");
  int w[K*(K+1)];
  for(m=0;m<K;m++)
    {
      for(n=1;n<K+1-m;n++)
	{
	  w[m+K*n]=0;
	  for(s=0;s<(x_div*y_div*z_div);s++)
	    {
	      if(catalog[m+N_STEPS*(s+x_div*y_div*z_div*n)])
		{
		  w[m+K*n]++;
		  for(i=0;i<THETA_PARTITIONS;i++)
		    {
		      correlation1=correlation_data[m+K*(s+(x_div*y_div*z_div)*(i+THETA_PARTITIONS*n))];	      	      
		      correlation_mean[m+K*(i+THETA_PARTITIONS*n)]+=correlation1;
		      //printf("c[%d_%d_%d_%d]=%lf \n",m,n,s,i,correlation1);
		    }
		}
	      else
		{
		  
		  for(i=0;i<THETA_PARTITIONS;i++)
		    {
		      correlation_data[m+K*(s+(x_div*y_div*z_div)*(i+THETA_PARTITIONS*n))]=-1000.0;   
		    }
		  
		}
	    }
	}
    }


 
  for(m=0;m<K;m++)
    {
      for(n=1;n<K+1-m;n++)
	{
	  for(i=0;i<THETA_PARTITIONS;i++)
	    {	  
	      //  correlation_mean[m+K*i]/=((x_div*y_div*z_div));
	      correlation_mean[m+K*(i+THETA_PARTITIONS*n)]=correlation_mean[m+K*(i+THETA_PARTITIONS*n)]/(w[m+K*n]);
	      //printf("c_mean[%d_%d_%d]=%lf \n",m,n,i,correlation_mean[m+K*(i+THETA_PARTITIONS*n)]);
	
	    }
	}
    }
      
      for(m=0;m<K;m++)
	{
	   for(n=1;n<K+1-m;n++)
	     {
	       for(s=0;s<(x_div*y_div*z_div);s++)
		 {
		   if(catalog[m+K*(s+x_div*y_div*z_div*n)])
		     {
		       for(i=0;i<THETA_PARTITIONS;i++)
			 {   
			   correlation_std[m+K*(i+THETA_PARTITIONS*n)]+=( correlation_data[m+K*(s+(x_div*y_div*z_div)*(i+THETA_PARTITIONS*n))]-correlation_mean[m+K*(i+THETA_PARTITIONS*n)] )*( correlation_data[m+K*(s+(x_div*y_div*z_div)*(i+THETA_PARTITIONS*n))]-correlation_mean[m+K*(i+THETA_PARTITIONS*n)] );
			   
			 }
		     }else{continue;}
		 }
	       
	     }
	}
  printf("hola 6 \n");
  FILE *file2[K*(K+1)];

    for(m=0;m<K;m++)
      {
	for(n=1;n<K+1-m;n++)
	  {
	    sprintf(filename,"./correlation_files/correlation_%d_%d.data",m,n);
	    file2[m+K*n]=fopen(filename,"w");
	    for(i=0;i<THETA_PARTITIONS;i++)
	      {
		//correlation_std[m+K*i]=sqrt(correlation_std[m+K*i]/(x_div*y_div*z_div));
		correlation_std[m+K*(i+THETA_PARTITIONS*n)]=sqrt( correlation_std[m+K*(i+THETA_PARTITIONS*n)]/w[m+K*n] );
		fprintf(file2[m+K*n],"%lf %lf %lf \n",theta[i],correlation_mean[m+K*(i+THETA_PARTITIONS*n)],correlation_std[m+K*(i+THETA_PARTITIONS*n)]);
		
	      }
	    fclose(file2[m+K*n]);
	  } 
	
      }
    
    printf("hola 7 \n");
    
   
    for(m=0;m<K;m++)
      {  
	for(n=1;n<K+1-m;n++)
	  {
	    for(i=0;i<THETA_PARTITIONS;i++)
	      {	  
		
		chisquare_mean_models[m+K*n]+=(correlation_data_obs[i]-correlation_mean[m+K*(i+THETA_PARTITIONS*n)])*(correlation_data_obs[i]-correlation_mean[m+K*(i+THETA_PARTITIONS*n)])/(THETA_PARTITIONS*(correlation_std_obs[i]*correlation_std_obs[i]+correlation_std[m+K*(i+THETA_PARTITIONS*n)]*correlation_std[m+K*(i+THETA_PARTITIONS*n)]));
	       

		
		chisquare_mean_models_obs[m+K*n]+=(correlation_data_obs[i]-correlation_mean[m+K*(i+THETA_PARTITIONS*n)])*(correlation_data_obs[i]-correlation_mean[m+K*(i+THETA_PARTITIONS*n)])/(THETA_PARTITIONS*(correlation_std_obs[i]*correlation_std_obs[i]));
	      }
	    
	    
	  }
      }




  

 printf("hola 9 \n");
 for(m=0;m<K;m++)
    {
      for(n=1;n<K+1-m;n++)
	{
	  for(s=0;s<(x_div*y_div*z_div);s++)
	    {
	      if(catalog[m+N_STEPS*(s+x_div*y_div*z_div*n)])
		{
		  for(i=0;i<THETA_PARTITIONS;i++)
		    {
		      
		      chisquare_all_models[m+K*(s+x_div*y_div*z_div*n)]+=(correlation_data_obs[i]- correlation_data[m+K*(s+(x_div*y_div*z_div)*(i+THETA_PARTITIONS*n))])*(correlation_data_obs[i]- correlation_data[m+K*(s+(x_div*y_div*z_div)*(i+THETA_PARTITIONS*n))])/(THETA_PARTITIONS*correlation_std_obs[i]*correlation_std_obs[i]);
		    }
		}
	      //printf("chisquiare[%d_%d_%d]=%lf \n",m,n,s,chisquare_all_models[m+K*(s+x_div*y_div*z_div*n)]);
	    }
	}
    }

printf("hola 10 \n");
 FILE *file5,*file6,*file7;
 file5=fopen("./correlation_files/chisquare_models.data","w");
 file7=fopen("./correlation_files/chisquare_models_obs.data","w");
 file6=fopen("./correlation_files/chisquare_all.data","w");
 
 
 size_t permutation_mean[K*(K+1)];
 size_t permutation_mean_obs[K*(K+1)];
 size_t permutation_all[K*(x_div*y_div*z_div)*(K+1)];
 
printf("hola 11 \n");
 gsl_sort_smallest_index(permutation_mean,K*(K+1),chisquare_mean_models,1,K*(K+1));
gsl_sort_smallest_index(permutation_mean_obs,K*(K+1),chisquare_mean_models_obs,1,K*(K+1));
 gsl_sort_smallest_index(   permutation_all, K*(x_div*y_div*z_div)*(K+1), chisquare_all_models, 1, K*(x_div*y_div*z_div)*(K+1)   );
 
printf("hola 12 \n");
 i=0;
 j=0;
 int d,b,c,e=0;
 FILE *corr_cats[x_div*y_div*z_div];
 for(m=0;m<K;m++)
   {
     for(n=0;n<K+1;n++)
       {
	 if(chisquare_mean_models[permutation_mean[i]]!=0)
	   {
	     fprintf(file5,"%d %d %lf %lf %lf %d \n",order_mn[permutation_mean[i]].p,order_mn[permutation_mean[i]].q,LOWEST_MASS+MASS_STEP*order_mn[permutation_mean[i]].p,LOWEST_MASS + MASS_STEP*(order_mn[permutation_mean[i]].p + order_mn[permutation_mean[i]].q  ),chisquare_mean_models[permutation_mean[i]],w[permutation_mean[i]]);

fprintf(file7,"%d %d %lf %lf %lf %d \n",order_mn[permutation_mean_obs[i]].p,order_mn[permutation_mean_obs[i]].q,LOWEST_MASS+MASS_STEP*order_mn[permutation_mean_obs[i]].p,LOWEST_MASS + MASS_STEP*(order_mn[permutation_mean_obs[i]].p + order_mn[permutation_mean_obs[i]].q  ),chisquare_mean_models_obs[permutation_mean_obs[i]],w[permutation_mean_obs[i]]);
	   }
	 for(s=0;s<(x_div*y_div*z_div);s++)
	   {
	     if(chisquare_all_models[permutation_all[j]]!=0)
	       {
		 d=order_mws[permutation_all[j]].p;
		 b=order_mws[permutation_all[j]].q;
		 c=order_mws[permutation_all[j]].r;
		 fprintf(file6,"%d %d %d  %lf %lf %lf \n",order_mws[permutation_all[j]].p,order_mws[permutation_all[j]].r,order_mws[permutation_all[j]].q,LOWEST_MASS+MASS_STEP*order_mws[permutation_all[j]].p,LOWEST_MASS+MASS_STEP*(order_mws[permutation_all[j]].p + order_mws[permutation_all[j]].r),chisquare_all_models[permutation_all[j]]);
		 if(e<x_div*y_div*z_div)
		   {
		     sprintf(filename,"./correlation_files/correlation_%d_%d_%d.data",d,c,b);
		     corr_cats[e]=fopen(filename,"w");
		     for(k=0;k<THETA_PARTITIONS;k++)
		       {
			 fprintf(corr_cats[e],"%lf %lf \n",theta[k], correlation_data[m+K*(s+(x_div*y_div*z_div)*(k+THETA_PARTITIONS*n))]);
		       }
		     fclose(corr_cats[e]);
		     e++;
		   }
	       }
	    
	     j++;
	     
	   }
	 i++;
       }
   }
 
 printf("hola 13 \n");
 fclose(file5);
 fclose(file6);
 fclose(file7);
 return 0;
 
}




int correlation_statistic_models(double correlation_data[],int catalog[])
{
  
  double *chisquare_all_models;
  
  
  
  int s,m,n,i,j,K,g;
  int o,p,t;
  int d,b,c,f,h;
  
  

  
  FILE *file6; 
  size_t *permutation_all;
  

  int x_div,y_div,z_div;
  x_div=floor(BOXSIZE/SUB_BOX_SIZE_X);
  y_div=floor(BOXSIZE/SUB_BOX_SIZE_Y);
  z_div=floor(BOXSIZE/SUB_BOX_SIZE_Z);
  
  
  K=N_STEPS;
  
  
  
  ordering *order_mws;
  order_mws=malloc(K*(x_div*y_div*z_div)*(K+1)*K*(x_div*y_div*z_div)*(K+1)*sizeof(ordering));
    
  
  printf("hola parce 1 \n");
  
  file6=fopen("./correlation_files/chisquare_between_boxes.data","w");
  
  
  printf("hola parce 2 \n");
  
  for(m=0;m<K;m++)
    { 
      for(n=1;n<K+1-m;n++)
	{
	  // printf("hola parcerito %d \n",n);
	  for(s=0;s<(x_div*y_div*z_div);s++)
	    {
	      
	      for(o=0;o<K;o++)
		{
		  for(p=1;p<K+1-m;p++)
		    {
		      
		      for(t=0;t<(x_div*y_div*z_div);t++)
			{
			  g=m + K*( s + x_div*y_div*z_div*( n + ( K + 1 )*( o + K*(t + x_div*y_div*z_div*p ) ) ) );
			  //  printf("hola parcerito %d %d %d %d \n",n,m,o,t);
			  order_mws[g].p=m;			  
			  order_mws[g].q=n;
			  order_mws[g].r=s;
			  order_mws[g].w=o;
			  order_mws[g].v=p;
			  order_mws[g].u=t;
			  
			}
		    }
		}
	    }
	}
      
    }
  
  printf("hola parce 3 \n");

  
  

  
  chisquare_all_models=malloc((x_div*y_div*z_div)*(K*(K+1))*(x_div*y_div*z_div)*(K*(K+1))*sizeof(double));
  
 
  printf("hola parce 4 \n");
 
 
  
 



  for(m=0;m<K;m++)
    {
      for(n=0;n<K+1;n++) 
	{ 
	  for(s=0;s<(x_div*y_div*z_div);s++)
	    {
	      for(o=0;o<K;o++)
		{
		  for(p=1;p<K+1-m;p++)
		    {
	  
		      for(t=0;t<(x_div*y_div*z_div);t++)
			{
			  g=m + K*( s + x_div*y_div*z_div*( n + ( K + 1 )*( o + K*(t + x_div*y_div*z_div*p ) ) ) );
			  if(t<=s)
			    {
			      chisquare_all_models[g]=10000;  
			    }
			  else
			    {
			      chisquare_all_models[g]=0;
			    }
			}
		    }
		}
	    }
	}
    }
 

  printf("hola parce 5 \n");
  
 

  for(m=0;m<K;m++)
    {
      for(n=1;n<K+1-m;n++)
	{
	  
	  for(s=0;s<(x_div*y_div*z_div);s++)
	    {
	      for(o=0;o<K;o++)
		{
		  for(p=1;p<K+1-m;p++)
		    {
	  
		      for(t=0;t<(x_div*y_div*z_div);t++)
			{
			  if(catalog[m+N_STEPS*(s+x_div*y_div*z_div*n)] && catalog[o+N_STEPS*(t+x_div*y_div*z_div*p)])
			    {
			      for(i=0;i<THETA_PARTITIONS;i++)
				{
				  g=m + K*( s + x_div*y_div*z_div*( n + ( K + 1 )*( o + K*(t + x_div*y_div*z_div*p ) ) ) );
				  f=m+K*(s+(x_div*y_div*z_div)*(i+THETA_PARTITIONS*n));
				  h=o+K*(t+(x_div*y_div*z_div)*(i+THETA_PARTITIONS*p));
			          if(correlation_data[f]<90000 && correlation_data[h]<90000)
                                    {
				     chisquare_all_models[g]+=( correlation_data[f] - correlation_data[h] ) * ( correlation_data[f] - correlation_data[h] );
                                    }
				}
			    }
			  //printf("chisquiare=%lf \n",chisquare_all_models[g]);
			}
		    }
		}
	    }
	}
    }

 
 permutation_all=malloc(K*(x_div*y_div*z_div)*(K+1)*K*(x_div*y_div*z_div)*(K+1)*sizeof(size_t));


 printf("hola parce 6 \n"); 
 
 
 // gsl_sort_smallest_index(   permutation_all,K*K*(x_div*y_div*z_div)*(K+1)*(K+1)*(x_div*y_div*z_div) , chisquare_all_models, 1, K*K*(x_div*y_div*z_div)*(K+1)*(K+1)*(x_div*y_div*z_div)   );
 
 gsl_sort_index (permutation_all,chisquare_all_models,1, K*K*(x_div*y_div*z_div)*(K+1)*(K+1)*(x_div*y_div*z_div) );

 printf("hola parce 7 \n"); 

 j=0;

 
 for(m=0;m<K;m++)
   {
     for(n=0;n<K+1;n++)
       {
	 for(s=0;s<(x_div*y_div*z_div);s++)
	   {
	     for(o=0;o<K;o++)
	       {
		 for(p=0;p<K+1;p++)
		   {
		     for(t=0;t<(x_div*y_div*z_div);t++)
		       {
			 if(chisquare_all_models[permutation_all[j]]!=0 && chisquare_all_models[permutation_all[j]]<9000)
			   {
			     d=order_mws[permutation_all[j]].p;
			     b=order_mws[permutation_all[j]].q;
			     c=order_mws[permutation_all[j]].r;
			     f=order_mws[permutation_all[j]].w;
			     g=order_mws[permutation_all[j]].v;
			     h=order_mws[permutation_all[j]].u;
			     fprintf(file6,"%5d %5d %5d  %10d %5d %5d %10.3lf \n", d, b, c, f, g, h, chisquare_all_models[permutation_all[j]]);
		
			   }
			 
			 j++;
			 
		       }
		   }
	       }
	   }
       }
   }
 printf("hola parce 8 \n"); 

 fclose(file6);

 return 0;
 
}
