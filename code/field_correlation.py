#!/usr/bin/env python

import LAE_correlation as correlation
import numpy as np
import pylab as P

field="full"
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
random_cat_number=16
distance=6558.3

print "computing correlation functions of the selected models"
dmh_path=pro_path+"data/dark_matter/FOF/"
laes_path=pro_path+"data/laes/FOF/"
n_i= int( box_length/x_width)
n_j= int( box_length/y_width)
n_k= int( box_length/z_depth)
n_models = n_i * n_j * n_k
    
ID_file=pro_path + "data/mock_survey/" + "ID_" + survey_type + "_surveys.dat"
    


model_info=np.loadtxt("./model_data.dat",skiprows=2)

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
std_correlation=np.empty([len(m_min_m),theta_bins])
    
for w in range( len(m_min_m) ):

    model_index=np.where( ( m_min_arr == m_min_m[w] ) & ( m_max_arr == m_max_m[w] ) & ( f_occ_arr == f_occ_m[w] )  )  
    m_min_s=m_min_arr[model_index] #_s reffers to sorted
    m_max_s=m_max_arr[model_index]
    f_occ_s=f_occ_arr[model_index]
    print m_min_s[0], m_max_s[0], f_occ_s[0]
    ID_survey_s=ID_survey_arr[model_index]
    model_prob_s=model_prob_arr[model_index]
    index_prob=np.argsort(model_prob_s)
    m_min_s=m_min_s[index_prob]
    m_max_s=m_max_s[index_prob]
    f_occ_s=f_occ_s[index_prob]
    ID_survey_s=ID_survey_s[index_prob]
    model_prob_s=model_prob_s[index_prob]
    
    ID_index=np.where(survey_ID==ID_survey_s[-1])
    S_ID=survey_ID[ID_index]
    
    ID_ini=S_ID[0]
    ID_end=int(ID_ini+obs_surveys)
    m_min=m_min_s[-1]
    m_max=m_max_s[-1]
    f_occ=f_occ_s[-1]
    model_prob=model_prob_s[-1]
    print "model:",w,"parameters:" ,m_min, m_max, f_occ,model_prob
    i_s=i_field[ID_ini:ID_end]
    j_s=j_field[ID_ini:ID_end]
    k_s=k_field[ID_ini:ID_end]
    corr=np.zeros( (len(i_s),theta_bins) )
    corr_peebles=np.zeros( (len(i_s),theta_bins) )
    corr_standard=np.zeros( (len(i_s),theta_bins) )
    corr_laes=np.zeros(theta_bins)
    
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
        halo_index=numbers[halo_index]
                    
        np.random.shuffle(halo_index)
            
            
        n_halos=np.size(halo_mass_sel)
        del halo_mass_sel
        n_laes=int(f_occ*n_halos)
            
        lae_index=halo_index[0:n_laes]
        x_laes=x_halos[lae_index]
        y_laes=y_halos[lae_index]
        
        del x_halos
        del y_halos
        del z_halos
        del halo_mass
            
            #random cat histogram generation
        P.xlabel(r'$\theta$', fontsize=16)
        P.ylabel(r"$\xi(\theta)$",fontsize=16)
        if(w==0):
            if(i==0):
                print w,i
                print "computing RR (it takes much time but it is only computed once)"
                x_random= x_width*np.random.random_sample(n_laes*random_cat_number)
                y_random=y_width*np.random.random_sample(n_laes*random_cat_number)
                RR,bins=correlation.RR_histogram(x_laes,y_laes,x_random,y_random,distance,theta_min,theta_max,theta_bins,cat_number=random_cat_number)
                del x_random
                del y_random
                print RR
             
        print "subcat number ",i,"i j k=",i_s[i],j_s[i],k_s[i]

        #random-survey histogram generation
        Xmin=x_width*i_s[i]
        Xmax=Xmin + x_width
        Ymin=y_width*j_s[i]
        Ymax=Ymin + y_width
                
                
        x_random= Xmin +  ( Xmax - Xmin )*np.random.random_sample(n_laes)
        y_random=Ymin +   ( Ymax - Ymin )*np.random.random_sample(n_laes)
            
        DR,bins=correlation.DR_histogram(x_laes,y_laes,x_random,y_random,distance,theta_min,theta_max,theta_bins,cat_number=1)
        
            #survey histogram generation
        DD,bins=correlation.DD_histogram(x_laes,y_laes,distance,theta_min,theta_max,theta_bins)
        
        corr[i,:]=correlation.landy_correlation(DD,RR,DR)
        print "CORR_landy=",corr[i,:]
            
    corr_laes=np.mean(corr,axis=0)
    std_corr=np.std(corr,axis=0)
    print "corr_landy=",corr_laes, "std_landy=",std_corr
        
    best_correlation[w,:]=corr_laes
    std_correlation[w,:]=std_corr
    dtheta=(theta_max - theta_min)/theta_bins
        
    correlation_data=np.empty(( np.size(corr_laes) , 3 ) )
    model='(Mmin,Mmax,focc)=({0},{1},{2})'.format(m_min, m_max, f_occ)
    model_name = 'model_{0}_{1}_{2}'.format(m_min, m_max, f_occ)
    filename=pro_path + "data/mock_survey/" + "correlation_best_models/" + survey_type + "_correlation_" + model_name + ".dat"
        
    angles = np.linspace( theta_min + dtheta/2.0 , theta_max - dtheta/2.0, theta_bins )
    correlation_data[:,0]=angles
    correlation_data[:,1]=best_correlation[w,:]
    correlation_data[:,2]=std_correlation[w,:]
        
    np.savetxt(filename,correlation_data)
        
    P.errorbar(correlation_data[:,0]+3.0*w, correlation_data[:,1], correlation_data[:,2],label=model,elinewidth=1.5)
    
    

file_plot=pro_path + "data/mock_survey/" + "correlation_best_models/" + survey_type + "_" + field  +"_"+ "correlation_selected_models" + ".png"
P.legend(shadow=False)
obs_correlation_file=pro_path + "data/obs/hayashino_whole_SSA22_field.txt"
obs_correlation=np.loadtxt(obs_correlation_file,skiprows=4)
#P.ylim(ymax=0.6)
P.xlim(xmax=1040)
    
P.errorbar(obs_correlation[0:theta_bins,0]-3.0, obs_correlation[0:theta_bins,1], obs_correlation[0:theta_bins,2],label="Hayashino et al 2004",elinewidth=3.0,fmt="o-")
P.legend(shadow=False)
P.title(survey_type)
P.savefig(file_plot)
P.figure()
