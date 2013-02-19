#include<mass_min_max.h>


int distance_mass_halos_counter(char catalog_file[],double xmin,double xmax,double ymin,double ymax,double zmin,double zmax)
{
  int j=0,n,n_particles,i=0;
  double xtest,ytest,ztest,vx,vy,vz,mass,sigma_r, sigma_v, r_sph,delta, xlamba;//distance;
  mass=1e15;
  
  FILE *f1;

  f1=fopen(catalog_file,"r");

  while(log10(mass)>LOWEST_MASS-0.1)
    {
      
      fscanf(f1,"%d %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf",&n,&xtest,&ytest,&ztest,&vx,&vy,&vz,&n_particles,&mass, &sigma_r, &sigma_v, &r_sph, &delta, &xlamba);
      
      i++;
      if(xtest>=xmin &&xtest<xmax && ytest>=ymin && ytest<ymax && ztest>=zmin && ztest<zmax )
	{
	  j++; 
	}	     
      
    }
  printf("El numero de objetos en el subcatálogo es de %d \n",j);
  
  
  return j;
  
}




int distance_mass_halos_selector(char catalog_file[],int N,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,double mas[],double x[],double y[],double z[])
{

  int n,n_particles,i=0,j=0;
  double xtest,ytest,ztest,vx,vy,vz,mass=1e15,sigma_r, sigma_v, r_sph,delta, xlamba;//distance;
  
  FILE *f1;
  f1=fopen(catalog_file,"r");

 while( log10(mass)>LOWEST_MASS-0.1)
    {
      fscanf(f1,"%d %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf",&n,&xtest,&ytest,&ztest,&vx,&vy,&vz,&n_particles,&mass, &sigma_r, &sigma_v, &r_sph, &delta, &xlamba);
      
      i++;
     
      if(xtest>=xmin && xtest<xmax  && ytest>=ymin && ytest<ymax && ztest>=zmin && ztest<zmax )
	{
	  
	  mas[j]=log10(mass);
	  x[j]=xtest;
	  y[j]=ytest;
	  z[j]=ztest;
	  j++;
	  
	}	     
    }
 
  fclose(f1);
  return 0;
  
}


