int DD_histogram(int N_cat,double X_cat[],double Y_cat[],double Z_cat[],double histogram_cat[]);
int RR_histogram(double histogram_rand[],int N_POINTS_RANDOM,double std[]);
int DR_histogram(int N_cat,int p,int q,int r,int N_cats_rand,double X_cat[],double Y_cat[],double Z_cat[],double histogram_dr[]);
int correlation_landy(double histogram_cat[],double histogram_rand[], double histogram_dr[],double correlation[]);
int correlation_standard(double histogram_cat[],double histogram_rand[],double correlation[]);
