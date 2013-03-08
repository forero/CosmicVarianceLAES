import numpy as np 
#import correlation_lib as corr

D_z=6558.3

'''
DD histogram
Computes the angular distribution of points in the observational catalog

'''
def DD_histogram(X,Y,distance,th_min,th_max,theta_bins):

    n_points=len(X)
    d_arr=[]
    for i in range(n_points):
        for j in range(i,n_points):
            d=(X[i] - X[j])*(X[i] - X[j]) + (Y[i] - Y[j])*(Y[i] - Y[j])
            d=np.sqrt(d)/distance
            d_arr=np.append(d_arr,d)
            

    theta=d_arr/distance           
    DD=np.histogram(theta,bins=theta_bins, range=(th_min,th_max))
    N=1.0*n_points*(1.0*n_points-1.0)/2.0
    DD=DD/N
    return DD

def RR_histogram(X,Y,distance,th_min,th_max,theta_bins,cat_number=1000):
    
    n_points=len(X)
    Xmin=1.0*np.floor( np.amin(X) )
    Xmax=1.0*np.ceil( np.amin(X) )
    Ymin=1.0*np.floor( np.amin(Y) )
    Ymax=1.0*np.ceil( np.amin(Y) )
    

    Xr= Xmin +  ( Xmax - Xmin )*np.random.random_sample(n_points*cat_number)
    Yr=Ymin +   ( Ymax - Ymin )*np.random.random_sample(n_points*cat_number)
    d_arr=[]
    for m in range(cat_number):
        for i in range(n_points):
            for j in range(i,n_points):
                d=(Xr[i + n_points*m ] - Xr[j + n_points*m ])*(Xr[i + n_points*m ] - Xr[j + n_points*m ]) + (Yr[i + n_points*m ] - Yr[j + n_points*m ])*(Yr[i + n_points*m ] - Yr[j + n_points*m ])
                d=np.sqrt(d)/distance
                d_arr=np.append(d_arr,d)
    theta=d_arr/distance
    RR=np.histogram(d_arr,bins=theta_bins, range=(th_min,th_max))
    N=1.0*n_points*(1.0*n_points-1.0)/2.0
    RR=RR/(N*cat_number)            
    return RR

def DR_histogram(X,Y,distance,th_min,th_max,theta_bins,cat_number=1000):
    n_points=len(X)
    Xmin=1.0*np.floor( np.amin(X) )
    Xmax=1.0*np.ceil( np.amin(X) )
    Ymin=1.0*np.floor( np.amin(Y) )
    Ymax=1.0*np.ceil( np.amin(Y) )
    

    Xr= Xmin +  ( Xmax - Xmin )*np.random.random_sample(n_points*cat_number)
    Yr=Ymin +   ( Ymax - Ymin )*np.random.random_sample(n_points*cat_number)
    d_arr=[]
    for m in range(cat_number): 
            for j in range(i,n_points):
                for i in range(n_points):
                    d=( X[i] - Xr[j + n_points*m] )*( X[i] - Xr[j + n_points*m] ) + ( Y[i] - Yr[j + n_points*m] )*( Y[i] - Yr[j + n_points*m] )
                    d=np.sqrt(d)/distance
                    d_arr=np.append(d_arr,d)
    theta=d_arr/distance  
    DR=np.histogram(theta,bins=theta_bins, range=(th_min,th_max))
    N=1.0*n_points*(n_points)
    DR=DR/(N*cat_number)            
    return DR


def landy_correlation(DD,RR,DR):
    return (DD - 2.0*DR + RR)/RR
def peebles_correlation(DD,DR):    
    return DD/DR - 1.0
def standard_correlation(DD,RR):
    return DD/RR - 1.0






'''--------------------------------------------
best_model_sel

Select the models with a KS prob greater that prob_theshold

Params:


prob_treshold=> Ks test probability to consider a model feasible

survey_type => Indicate the kind of mock survey
to consider for the model.

  Description:
 
  match: The mock-field distribution exactly takes the observational distribution
         in position and number. 

  random: The number of fields is the same as the observational fields but  
          mock-field positions are taken randomly 
 
  full: All the mock-fields are taken into account for the computation

pro_path==> Path with the KS test file and ID file

return:

best models => 5 column array with the information of the models  models
where  ks_prob>=p_threshold 

   column 1: m_min 
   column 2: m_max
   column 3: f_occ
   column 4: ID_survey
   column 5: KS test probability

'''



def best_model_sel(prob_treshold,survey_type="match", pro_path="/home/jemejia/CosmicVarianceLAES/"):
   
    ID_file=pro_path + "data/mock_survey/" + "ID_" + survey_type + "_surveys.dat"
    p_values_file= pro_path + "data/mock_survey/" + "p_values_ID_" + survey_type + "_surveys.dat"

    ks_data=np.loadtxt(p_values_file) 
    

    m_min_arr = ks_data[:,0]
    m_max_arr =ks_data[:,1]
    f_occ_arr =ks_data[:,2]
    ID_survey_arr = ks_data[:,3]

    model_prob_arr = ks_data[:,4]
    
    #choosing th models with ks test probabilities greater than prob_treshold
    index=np.where(model_prob_arr>=prob_treshold)
    
    best_models=np.empty( [ np.size(index) , np.size( ks_data[0,:] ) ] )
    del(ks_data)

    best_models[:,0]=m_min_arr[index]
    best_models[:,1]=m_max_arr[index]
    best_models[:,2]=f_occ_arr[index]
    best_models[:,3]= ID_survey_arr[index] 
    best_models[:,4]=model_prob_arr[index]

    return best_models
    




'''--------------------------------------------
Best_model_correlation

Params:

best_model_array ==> It takes as input the output of the function best_model_sel
 ,i.e, this is a  5 column array with the information of the models  models
where  ks_prob>=p_threshold 

   column 1: m_min 
   column 2: m_max
   column 3: f_occ
   column 4: ID_survey
   column 5: KS test probability


theta_min ==> Minimal value of the separation of galaxies angle to compute the two angularcorrelation function
theta_max ==> maximal value of the separation of galaxies angle to compute the two angularcorrelation function
theta_bins ==> number of bins to divide the interval (theta_min, theta_max)

survey_type => Indicate the kind of mock survey
to consider for the model.

  Description:
 
  match: The mock-field distribution exactly takes the observational distribution
         in position and number. 

  random: The number of fields is the same as the observational fields but  
          mock-field positions are taken randomly 
 
  full: All the mock-fields are taken into account for the computation


distance ==> Comoving mean distance to the area of survey. 
obs_surveys ==> Number of fields where the observational survey is computed
x_width ==> Transversal observational width of the surveys in Mpc (comoving)
y_width ==> Transversal observational width of the surveys in Mpc (comoving)
z_depth ==> redshift depth  of the survey in Mpc (comoving)
box_length ==> Length of the DM simulation


pro_path==> Path with the KS test file and ID file


return: 

best_correlation ---> A two dimensional array with the mean correlation function
                      of the best models.

                      If w is the number of the model then best_correlation[w,:]
                      is the mean correlation function  of such a model evaluated
                      in the range (theta_min,theta_max) in theta_bins bins

std_correlation  --->  A two dimensional array with the standard deviation of the 
                      correlation function of the best models. This standard deviation
                      comes from the diferent values that the correlation funtion
                      take in the diferent subcatalogs that is computed

                      If w is the number of the model then best_correlation[w,:]
                      is the mean correlation function  of such a model evaluated
                      in the range (theta_min,theta_max) in theta_bins bins

angles           -->  angels where the correlation function has been computed 
'''


def best_model_correlation(best_model_array, theta_min,theta_max,theta_bins, survey_type="match", distance=6558.3, obs_surveys=12,x_width=46.0,y_width=35.0, z_depth=41.0 ,box_length=250, pro_path="/home/jemejia/CosmicVariance/"):

    
    
    dmh_path=project_path+"data/dark_matter/FOF/"
    laes_path=project_path+"data/laes/FOF/"
    n_i= int( box_length/x_width)
    n_j= int( box_length/y_width)
    n_k= int( box_length/z_depth)
    n_models = n_i * n_j * n_k
    
    ID_file=pro_path + "data/mock_survey/" + "ID_" + survey_type + "_surveys.dat"
    

   
   
    
    ID_data=np.loadtxt(ID_file,dtype='int') 
    

    survey_ID=ID_data[:,0]
    field_ID=ID_data[:,1]
    i_field = ID_data[:,2]
    j_field = ID_data[:,3]
    k_field = ID_data[:,4]

    moc_surveys=survey_ID[-1]

    ID_arr=best_models[:,3]    
    index_eq_ID=np.where(ID_arr== 1)

    cat_number= index_eq_ID[0]-1
 
    i_fields_to_measure=[]
    j_fields_to_measure=[]
    k_fields_to_measure=[]
    m_min_to_measure=[]
    m_max_to_measure=[]
    f_occ_to_measure=[]
    #choosing the subcatalogs of the best fields.
    best_correlation=np.empty([len(ID_arr),theta_bins])
    std_correlation=np.empty([len(ID_arr),theta_bins])
    for w in range( len(ID_arr) ):
        index=np.where( survey_ID == int(ID_arr[w]) )
        
        S_ID=survey_ID[index]
        ID_ini=S_ID[0]
        ID_end=int(ID_ini+obs_surveys)
        m_min=best_models[w,0]
        m_max=best_models[w,1]
        f_occ=best_models[w,2]
        print m_min, m_max, f_occ
        i_s=i_field[ID_ini:ID_end]
        j_s=j_field[ID_ini:ID_end]
        k_s=k_field[ID_ini:ID_end]
        
        
        
        corr=np.zeros( (len(i_s)*len(j_s)*len(k_s),theta_bins) )
        corr_laes=np.zeros(theta_bins)
        for i in i_s:
            for j in j_s:
                for k in k_s:
                    s=i + len(j_s)*( j + len(k_s)*k )
                    dmh_filename=dmh_path+"halos_bolshoi_"+str(i)+"-"+str(j)+"-"+str(k)+".csv"
                    halos_prop=np.loadtxt(dmh_filename,delimiter=",",skiprows=12)
                    
                    halo_mass=halos_prop[:,4]
                    print halo_mass
                    x_halos=halos_prop[:,0]
                    y_halos=halos_prop[:,1]
                    z_halos=halos_prop[:,2]

                    halo_index=np.where( (halo_mass< m_max) & (halo_mass> m_min) )
                    halo_mass_sel=halo_mass[halo_index]
                    print halo_index
                    a=np.random.shuffle(halo_index)
                    print halo_index,a
                    x_sel=x_halos[halo_index]
                    y_sel=y_halos[halo_index]
                    
                    n_halos=np.size(halo_mass_sel)
                    n_laes=int(f_occ*n_halos)
                    #print halo_index, n_laes
                    lae_index=halo_index[0:n_laes]
                    x_laes=x_halos[lae_index]
                    y_laes=y_halos[lae_index]
                    
                    
                    
                    DD_histogram(x_laes,y_laes,distance,th_min,th_max,theta_bins)
 
                    RR_histogram(x_laes,y_laes,distance,theta_min,theta_max,theta_bins,cat_number=1000)

                    DR_histogram(x_laes,y_laes,distance,theta_min,theta_max,theta_bins,cat_number=1000)
                    corr_laes[s,:]=landy_correlation(DD,RR,DR)
                    
        corr_laes=np.mean(corr_laes,axis=0)
        std_corr=np.std(corr_laes,axis=0)
        best_correlation[w,:]=corr_laes
        std_correlation[w,:]=std_corr
        dtheta=(theta_min - theta_max)/theta_bins
        angles=np.linspace(theta_min + dtheta/2.0 ,theta_max - dtheta/2.0,theta_bins -1)  
    return best_correlation,std_correlation, angles


project_path = "/home/jemejia/CosmicVarianceLAES/"
p_treshold=0.1
theta_min=40
theta_max=1320
theta_bins=13

best_models=best_model_sel(p_treshold, pro_path=project_path)
best_correlation, std_correlation=best_model_correlation(best_models, theta_min,theta_max,theta_bins,pro_path=project_path)
