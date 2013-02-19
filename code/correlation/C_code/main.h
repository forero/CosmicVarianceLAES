#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#define R_MIN 140.0 
#define R_MAX 140.0
#define DELTA_LIM 1.6
#define N_OBS_SAMPLE 246
#define N_STEPS 160
//#define nlaes[27][200]

#define N_POINTS_RANDOM 1000
#define THETA_MIN 0.0//arcseconds
#define THETA_MAX 1520//arcseconds
#define THETA_PARTITIONS 19
#define D_z 4747500
#define H 0.719
#define PHO_M 3.7e10
#define BOXSIZE 50000.0
#define N_GRID_X 40
#define N_GRID_Y 40
#define N_GRID_Z 40
#define GRID_SIZE 50000.0
#define CELL_SIZE 1250.0
#define N_LAES 1000
#define USE "./a.out density_catalog"

//cabeceras del gsl
#include<gsl/gsl_sort_double.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_rng.h>

//estructuras del propgrama

  typedef struct position{
    double x;
    double y;
    double z;
  }position;

typedef struct n_laes_grid
{
  int p;
  int q;
  int r;
  int n_laes;
}n_laes_grid;

typedef struct grid_str
{    
  int N_cells;
  double *cargo;
  float cell_size;
  float smoothing;
  float x_level;
  float y_level;
  float z_level;
  int N_grid_x; 
  int N_grid_y;
  int N_grid_z;
  float a_exp;
  
} grid;

#include<load_cube_file_cic.h>
#include<lowest_to_highest.h>
#include<laes_selector.h>
#include<compute_correlation.h>
#include<randomsample.h>
