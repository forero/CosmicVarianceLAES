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

#include<correlation.h>

//VARIABLES GLOBALES
#define RMIN 50.0
#define RMAX 51.25
#define N_CATS_RANDOM 1
#define N_CATS_DR 1
#define THETA_MIN 0.0 //arcsec
#define THETA_MAX 1040.0 //arcsec
#define THETA_PARTITIONS 13
#define D_z 6558.3
#define H 0.709
#define PHO_M 3.7e10
#define BOXSIZE 250.0
#define GRID_SIZE 250.0

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
#define MOCK_SURVEYS 15
#define N_MODELS 101
#define HALO_MASS_DEV 0.127725e10





