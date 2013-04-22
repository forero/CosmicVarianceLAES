#!/usr/bin/env python

import LAE_correlation as correlation
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
best_corr_laes_max=np.empty([len(m_min_m),theta_bins])
std_correlation=np.empty([len(m_min_m),theta_bins])
std_correlation_max=np.empty([len(m_min_m),theta_bins])
n_models=len(m_min_m)
#n_models=1 #to delete

fig1=P.figure()
fig2=P.figure()
ax_mean=fig1.add_subplot(111)
ax_max=fig2.add_subplot(111)

for w in range( n_models ):

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
    n_mocks=np.size(S_ID)
    print "n_mocks=" + str(n_mocks)
    corr=np.zeros( (n_mocks,theta_bins) )
    
    #corr_peebles=np.zeros( (n_mocks,theta_bins) )
    #corr_standard=np.zeros( (n_mocks,theta_bins) )
    corr_laes=np.zeros(theta_bins)
    
    
    corr_max_density=np.zeros( (n_mocks,theta_bins) )
    corr_laes_max=np.zeros(theta_bins)
    #n_mocks=1 #to delete
    for j in range( n_mocks ): 
        ID_ini=S_ID[j]
        ID_end=int(ID_ini+obs_surveys)
        m_min=m_min_s[-1]
        m_max=m_max_s[-1]
        f_occ=f_occ_s[-1]
        model_prob=model_prob_s[-1]
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
            halo_index=numbers[halo_index]
                    
            np.random.shuffle(halo_index)
            
            
            n_halos=np.size(halo_mass_sel)
            del halo_mass_sel
        
            number_laes=int(f_occ*n_halos)
            #number_laes=100 #to delete
            n_laes=np.append(n_laes, number_laes )    
            
            lae_index=halo_index[0:number_laes]
            x_laes=np.append( x_laes, x_halos[lae_index] )
            y_laes=np.append( y_laes, y_halos[lae_index] )
        

            del x_halos
            del y_halos
            del z_halos
            del halo_mass
            
            #random cat histogram generation

        
        
        
            #print "random catalog generator
            x_random=np.append( x_random, x_width*np.random.random_sample(number_laes*random_cat_number) )
            y_random=np.append( y_random, y_width*np.random.random_sample(number_laes*random_cat_number) )
        
        
             
            print "subcat number ",i,"i j k=",i_s[i],j_s[i],k_s[i]

            #random-survey catalog generation
            Xmin=x_width*i_s[i]
            Xmax=Xmin + x_width
            Ymin=y_width*j_s[i]
            Ymax=Ymin + y_width          
            x_drand= np.append( x_drand , Xmin +  ( Xmax - Xmin )*np.random.random_sample(number_laes) )
            y_drand= np.append( y_drand , Ymin +   ( Ymax - Ymin )*np.random.random_sample(number_laes) )
            
        
        
        
        
        
        
        
        #RR,bins=correlation.RR_histogram(x_laes,y_laes,x_random,y_random,distance,theta_min,theta_max,theta_bins,cat_number=random_cat_number)    
        print "computing DD"
        DD,bins=correlation.DD_histogram(x_laes,y_laes,distance,theta_min,theta_max,theta_bins)
        print "computed DD"
        print "computing DR"
        DR,bins=correlation.DR_histogram(x_laes,y_laes,x_drand,y_drand,distance,theta_min,theta_max,theta_bins,cat_number=1)
        print "DR computed"
        
        corr[j,:]=correlation.peebles_correlation(DD,DR)

        max_density_index=np.argmax(n_laes)
        number_laes=n_laes[max_density_index]
        lae_pos_ini=np.sum( n_laes[ 0 : max_density_index ] )
        lae_pos_end=lae_pos_ini + number_laes
        x_laes_max = x_laes[lae_pos_ini:lae_pos_end]
        y_laes_max = y_laes[lae_pos_ini:lae_pos_end]
        x_random_max = x_random[lae_pos_ini:lae_pos_end]
        y_random_max = y_random[lae_pos_ini:lae_pos_end]
        x_drand_max = x_drand[lae_pos_ini:lae_pos_end]
        y_drand_max = y_drand[lae_pos_ini:lae_pos_end]
        #RR,bins=correlation.RR_histogram(x_laes_max,y_laes_max,x_random_max,y_random_max,distance,theta_min,theta_max,theta_bins,cat_number=random_cat_number)    
        DR,bins=correlation.DR_histogram(x_laes_max,y_laes_max,x_drand_max,y_drand_max,distance,theta_min,theta_max,theta_bins,cat_number=1)
        DD,bins=correlation.DD_histogram(x_laes_max,y_laes_max,distance,theta_min,theta_max,theta_bins)
        corr_max_density[j,:]=correlation.peebles_correlation(DD,DR)


        print "CORR_landy=",corr[j,:]

    corr_laes=np.mean(corr,axis=0)
    corr_laes_max=np.mean(corr_max_density,axis=0)
    std_corr=np.std(corr,axis=0)
    std_corr_max=np.std(corr_max_density,axis=0)
    print "corr_landy=",corr_laes, "std_landy=",std_corr
    print "corr_max=",  corr_laes_max, "std_max=", std_corr_max
    best_correlation[w,:]=corr_laes
    best_corr_laes_max[w,:]=corr_laes_max
    std_correlation[w,:]=std_corr
    std_correlation_max[w,:]=std_corr_max
    dtheta=(theta_max - theta_min)/theta_bins
        
    correlation_data=np.empty(( np.size(corr_laes) , 3 ) )
    print np.size(corr_laes) 
    model='(Mmin,Mmax,focc)=({0},{1},{2})'.format(m_min, m_max, f_occ)
    model_name = 'model_{0}_{1}_{2}'.format(m_min, m_max, f_occ)
    filename=pro_path + "data/mock_survey/" + "correlation_best_models/" + survey_type + "_correlation_" + str(w) + ".dat"
        
    angles = np.linspace( theta_min + dtheta/2.0 , theta_max - dtheta/2.0, theta_bins )
    

    correlation_data[:,0]=angles
    correlation_data[:,1]=best_correlation[w,:]
    correlation_data[:,2]=std_correlation[w,:]
    
    np.savetxt(filename,correlation_data)
    
    #fig1=P.figure()
    #fig2=P.figure()
    #ax_mean=fig1.add_subplot(111)
    #ax_max=fig2.add_subplot(111)

    ax_mean.errorbar(correlation_data[:,0]+3.0*w, correlation_data[:,1], correlation_data[:,2],label=model,elinewidth=1.5)

    correlation_data[:,0]=angles
    correlation_data[:,1]=corr_laes_max
    correlation_data[:,2]=std_corr_max
    
    filename=pro_path + "data/mock_survey/" + "correlation_best_models/" + survey_type + "_maxcorrelation_" + str(w) + ".dat"
    np.savetxt(filename,correlation_data)

    
    
    ax_max.errorbar(angles + 3.0*w , corr_laes_max , std_corr_max,label=model, elinewidth=1.5)
    

file_plot=pro_path + "data/mock_survey/" + "correlation_best_models/" + survey_type + "_" + field  +"_"+ "correlation_selected_models" + ".png"
ax_mean.legend(shadow=False)
ax_max.legend(shadow=False)
obs_correlation_file=pro_path + "data/obs/hayashino_whole_SSA22_field.txt"
obs_correlation=np.loadtxt(obs_correlation_file,skiprows=4)
#P.ylim(ymax=0.6)
ax_mean.set_xlim(xmax=1040)
ax_max.set_xlim(xmax=1040)
ax_mean.set_xlabel(r'$\theta$', fontsize=16)
ax_max.set_xlabel(r'$\theta$', fontsize=16)
ax_mean.set_ylabel(r"$\xi(\theta)$",fontsize=16)    
ax_max.set_ylabel(r"$\xi(\theta)$",fontsize=16)    
ax_mean.errorbar(obs_correlation[0:theta_bins,0]-3.0, obs_correlation[0:theta_bins,1], obs_correlation[0:theta_bins,2],label="Hayashino et al 2004",elinewidth=3.0,fmt="o-")
ax_max.errorbar(obs_correlation[0:theta_bins,0]-3.0, obs_correlation[0:theta_bins,1], obs_correlation[0:theta_bins,2],label="Hayashino et al 2004",elinewidth=3.0,fmt="o-")
ax_mean.legend(shadow=False)
ax_max.legend(shadow=False)
ax_mean.set_title(survey_type)
ax_max.set_title("highest density region")

fig1.savefig("mean.png")

fig2.savefig("max.png")
P.show()
#P.figure()


