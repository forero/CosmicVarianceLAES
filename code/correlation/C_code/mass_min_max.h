//cabeceras del C
#include <stdlib.h>
#include<stdio.h>
#include<math.h>
#include <time.h>
#include <gsl/gsl_roots.h>
#include<string.h>




//Cabeceras de gsl
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_rng.h>
#include<gsl/gsl_sort_double.h>



//Cabeceras propias
#include<randomsample.h>
#include<mass_distance_selector.h>
#include<theta_histogram.h>
#include<correlation.h>
#include<angular_distribution.h>
#include<density_cluster_calculation.h>
#include<mass_function.h>
#include<LAE_catalog_generator.h>
#include <correlation_statistic.h>
//VARIABLES GLOBALES
#define RMIN 50.0
#define RMAX 51.25
#define N_CATS_RANDOM 20000
#define N_CATS_DR 100
#define THETA_MIN 1.30 //arcsec
#define THETA_MAX 3.50 //arcsec
#define THETA_PARTITIONS 11
#define D_z 8275
#define H 0.709
#define PHO_M 3.7e10
#define BOXSIZE 250.0
#define GRID_SIZE 250.0
#define N_LAES_OBS 111
//#define SUB_BOX_SIZE 40.0
#define SUB_BOX_SIZE_X_OBS 46.0
#define SUB_BOX_SIZE_Y_OBS 35.0
#define SUB_BOX_SIZE_Z_OBS 41.0
#define SUB_BOX_SIZE_X 46.0
#define SUB_BOX_SIZE_Y 35.0
#define SUB_BOX_SIZE_Z 41.0
#define LOWEST_MASS 10.0
#define HIGHEST_MASS 13.0
#define MASS_STEP 0.10
#define N_STEPS 25
#define HALO_MASS_DEV 0.127725e10

typedef struct corr_func
{
  double theta[THETA_PARTITIONS];
}corr_func;

typedef struct ordering
{
  int p;
  int q;
  int r;
  int w;
  int v;
  int u;
}ordering;


