#!/usr/bin/env python

#import LAE_correlation as correlation
import numpy as np
import pylab as P

field="large"
survey_type="match"
theta_max=1040
theta_min=0
theta_bins=13
pro_path="/home/jemejia/CosmicVarianceLAES/"
box_length=250.0
x_width=46.0
y_width=35.0
z_depth=41.0
obs_surveys=12
random_cat_number=1
distance=6558.3
randcats=16


dmh_path=pro_path+"data/dark_matter/FOF/"
laes_path=pro_path+"data/laes/FOF/"
n_i= int( box_length/x_width)
n_j= int( box_length/y_width)
n_k= int( box_length/z_depth)
n_models = n_i * n_j * n_k
    
ID_file=pro_path + "data/mock_survey/" + "ID_" + survey_type + "_surveys.dat"
    

models_file=pro_path + "data/mock_survey/" + "BestModels_match_mock.dat" 
model_info=np.loadtxt(models_file)

m_min_m=model_info[:,0] #_m refers to model
m_max_m=model_info[:,1]
f_occ_m=model_info[:,2]
ID_file=pro_path + "data/mock_survey/" + "ID_" + survey_type + "_surveys.dat"
p_values_file= pro_path + "data/mock_survey/" + "p_values_FOF_ID_" + survey_type + "_surveys.dat"





ID_data=np.loadtxt(ID_file,dtype='int')




ks_data=np.loadtxt(p_values_file) 

survey_ID=ID_data[:,0]
field_ID=ID_data[:,1]
i_field = ID_data[:,2]
j_field = ID_data[:,3]
k_field = ID_data[:,4]

del ID_data

m_min_arr = ks_data[:,0]
m_max_arr =ks_data[:,1]
f_occ_arr =ks_data[:,2]
ID_survey_arr = ks_data[:,3]
model_prob_arr = ks_data[:,4]


del ks_data
best_correlation=np.empty([len(m_min_m),theta_bins])
best_corr_laes_max=np.empty([len(m_min_m),theta_bins])
std_correlation=np.empty([len(m_min_m),theta_bins])
std_correlation_max=np.empty([len(m_min_m),theta_bins])
n_models=len(m_min_m)
#n_models=1 #to delete




for w in range( n_models ):

    model_index=np.where( ( m_min_arr == m_min_m[w] ) & ( m_max_arr == m_max_m[w] ) & ( f_occ_arr == f_occ_m[w] )  )  
    m_min_s=m_min_arr[model_index] #_s reffers to sorted
    m_max_s=m_max_arr[model_index]
    f_occ_s=f_occ_arr[model_index]
    print m_min_s[0], m_max_s[0], f_occ_s[0]
    ID_survey_s=ID_survey_arr[model_index]
    model_prob_s=model_prob_arr[model_index]
    print "number of mocks=",np.size(model_prob_s)
    index_prob=np.argsort(model_prob_s)
    m_min_s=m_min_s[index_prob]
    m_max_s=m_max_s[index_prob]
    f_occ_s=f_occ_s[index_prob]
    ID_survey_s=ID_survey_s[index_prob]
    model_prob_s=model_prob_s[index_prob]
    
    ID_index=np.where(survey_ID==ID_survey_s[-1])
    S_ID=survey_ID[ID_index]
    #n_mocks=np.size(S_ID)
    n_mocks=np.size(model_prob_s)
    print n_mocks
    print "n_mocks=" + str(n_mocks)
    
    #n_mocks=1 #to delete
    for j in range( n_mocks ): 
        
        ID_ini=j*obs_surveys
        ID_end=int(ID_ini+obs_surveys)
        m_min=m_min_s[-1]
        m_max=m_max_s[-1]
        f_occ=f_occ_s[-1]
        model_prob=model_prob_s[j]
        print "model:",w,"parameters:" ,m_min, m_max, f_occ,model_prob
        i_s=i_field[ID_ini:ID_end]
        j_s=j_field[ID_ini:ID_end]
        k_s=k_field[ID_ini:ID_end]

    
        x_laes=[]
        y_laes=[]
        x_random=[]
        y_random=[]
        x_drand=[]
        y_drand=[]
        n_laes=[] 
        
    
        if(field=="large"):
            i_range=7
            print "large field"
            
        else:
            i_range=np.size(i_s)
            print "full field"
            print "number of sub-catalogs=",i_range
        
        
        for i in range( i_range ):
            
            dmh_filename=dmh_path+"halos_bolshoi_"+str(i_s[i])+"-"+str(j_s[i])+"-"+str(k_s[i])+".csv"
            halos_prop=np.loadtxt(dmh_filename,delimiter=",",skiprows=12)
            
            halo_mass=halos_prop[:,4]
        
            x_halos=halos_prop[:,0]
            y_halos=halos_prop[:,1]
            z_halos=halos_prop[:,2]
            del halos_prop
            numbers=np.arange(len(halo_mass))
            halo_index=np.where( (halo_mass< m_max) & (halo_mass> m_min) )
            halo_mass_sel=halo_mass[halo_index]
            n_halos=np.size(halo_mass_sel)
            number_laes=int(f_occ*n_halos)
            n_laes=np.append(n_laes, number_laes )    
            halo_index=numbers[halo_index]
            for l in range(randcats): 
                np.random.shuffle(halo_index)
                lae_index=halo_index[0:number_laes]
                x_laes=np.append( x_laes, x_halos[lae_index] )
                y_laes=np.append( y_laes, y_halos[lae_index] )
            
            del halo_mass_sel
            del x_halos
            del y_halos
            del z_halos
            del halo_mass
            
            print "subcat number ",i,"i j k=",i_s[i],j_s[i],k_s[i]

            #random-survey catalog generation
            Xmin=x_width*i_s[i]
            Xmax=Xmin + x_width
            Ymin=y_width*j_s[i]
            Ymax=Ymin + y_width

            for l in range(randcats): 
                x_drand= np.append( x_drand , Xmin +  ( Xmax - Xmin )*np.random.random_sample(number_laes) )
                y_drand= np.append( y_drand , Ymin +   ( Ymax - Ymin )*np.random.random_sample(number_laes) )
            
            
        
        model_catalog_filename=pro_path +  "data/laes/mock_cat/model_"+str(w)+"_mock_"+str(j)+".txt"
    



        lae_cat=np.empty( [ np.size(x_laes) , 4 ] )
        lae_cat[:,0]=x_laes
        lae_cat[:,1]=y_laes
        lae_cat[:,2]=x_drand
        lae_cat[:,3]=y_drand
        np.savetxt(model_catalog_filename,lae_cat)
        
        
        
        

        max_density_index=np.argmax(n_laes)
        number_laes=n_laes[max_density_index]
        lae_pos_ini=randcats*np.sum( n_laes[ 0 : max_density_index ] )
        lae_pos_end=lae_pos_ini + randcats*number_laes
        x_laes_max = x_laes[lae_pos_ini:lae_pos_end]
        y_laes_max = y_laes[lae_pos_ini:lae_pos_end]
        #x_random_max = x_random[lae_pos_ini:lae_pos_end]
        #y_random_max = y_random[lae_pos_ini:lae_pos_end]
        x_drand_max = x_drand[lae_pos_ini:lae_pos_end]
        y_drand_max = y_drand[lae_pos_ini:lae_pos_end]
        


        lae_max=np.empty( [ np.size(x_laes_max) , 4 ] )
        lae_max[:,0]= x_laes_max
        lae_max[:,1]= y_laes_max
        lae_max[:,2]= x_drand_max
        lae_max[:,3]= y_drand_max
        

        model_catalog_filename=pro_path +  "data/laes/mock_cat/maxden_model_" + str(w) + "_mock_" + str(j) + ".txt"
        
        
        np.savetxt(model_catalog_filename,lae_max)


        n_mean=np.mean(n_laes)    
        mean_density_index=np.argmin( np.abs( n_laes-n_mean ) )
        number_laes=n_laes[mean_density_index]
        lae_pos_ini=randcats*np.sum( n_laes[ 0 : mean_density_index ] )
        lae_pos_end=lae_pos_ini + randcats*number_laes
        x_laes_mean = x_laes[lae_pos_ini:lae_pos_end]
        y_laes_mean = y_laes[lae_pos_ini:lae_pos_end]
        #x_random_max = x_random[lae_pos_ini:lae_pos_end]
        #y_random_max = y_random[lae_pos_ini:lae_pos_end]
        x_drand_mean = x_drand[lae_pos_ini:lae_pos_end]
        y_drand_mean = y_drand[lae_pos_ini:lae_pos_end]


        lae_mean=np.empty( [ np.size(x_laes_mean) , 4 ] )
        lae_mean[:,0]= x_laes_mean
        lae_mean[:,1]= y_laes_mean
        lae_mean[:,2]= x_drand_mean
        lae_mean[:,3]= y_drand_mean
        

        model_catalog_filename=pro_path +  "data/laes/mock_cat/meanden_model_" + str(w) + "_mock_" + str(j) + ".txt"
        
        
        np.savetxt(model_catalog_filename,lae_mean)
