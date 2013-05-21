import numpy as np
import pylab as P

n_models=101
n_mocks=12
theta_bins=13
pro_path="/home/jemejia/CosmicVarianceLAES/"


obs_correlation_file=pro_path + "data/obs/hayashino_whole_SSA22_field.txt"
obs_correlation=np.loadtxt(obs_correlation_file,skiprows=4)
models_file=pro_path + "data/mock_survey/" + "BestModels_match_mock.dat" 
model_info=np.loadtxt(models_file)
m_min_m=model_info[:,0] #_m refers to model
m_max_m=model_info[:,1]
f_occ_m=model_info[:,2]


mean_correlation= np.zeros( [n_models,theta_bins] )
std_correlation = np.zeros( [n_models,theta_bins] )
chi_individual  = np.zeros( [n_models,n_mocks] )
for i in range(n_models):
    corr=np.zeros( (n_mocks,theta_bins) )
    
    for j in range(n_mocks):
        model_correlation_filename=pro_path +  "data/laes/correlation/maxden_model_" + str(i) + "_mock_" + str(j) + ".txt"
        
        lae_correlation  = np.loadtxt(model_correlation_filename)
        corr[j,:] = lae_correlation[:,1]
        #chi_individual[i,j]= chi_individual[i,j] + (mean_correlation[i,j] - obs_correlation[j,1])*(mean_correlation[i,j] - obs_correlation[j,1])/(std_correlation[i,j]*std_correlation[i,j] + obs_correlation[j,2]*obs_correlation[j,2])
    std_correlation[i,:]  = np.std(corr,axis=0)
    mean_correlation[i,:] = np.mean(corr,axis=0)

    chisquare=np.zeros([n_models,theta_bins])
    chi=np.zeros(n_models)
for i in range(n_models):
    for j in range(theta_bins):
        chisquare[i,j]= chisquare[i,j] + (mean_correlation[i,j] - obs_correlation[j,1])*(mean_correlation[i,j] - obs_correlation[j,1])/(std_correlation[i,j]*std_correlation[i,j] + obs_correlation[j,2]*obs_correlation[j,2])
        chi[i]= chi[i] + chisquare[i,j] 
corr_stat_file=pro_path + "data/mock_survey/" + "correlation_stat.dat"
corr_stat=np.empty( [ n_models , 4 ] )
chi=chi/theta_bins
sort_index=np.argsort(chi)
corr_stat[:,0]=m_min_m[sort_index]
corr_stat[:,1]=m_max_m[sort_index]
corr_stat[:,2]=f_occ_m[sort_index]
corr_stat[:,3]=chi[sort_index]
np.savetxt(corr_stat_file,corr_stat,fmt='%2.3lf')
fig1=P.figure()
ax_mean=fig1.add_subplot(111)
model='(Mmin,Mmax,focc)=({0},{1},{2})'.format(m_min_m[sort_index[12]], m_max_m[sort_index[12]], f_occ_m[sort_index[12]])
ax_mean.errorbar(obs_correlation[0:theta_bins,0],mean_correlation[sort_index[12],:], std_correlation[sort_index[12],:], label=model,elinewidth=1.5)

ax_mean.set_xlabel(r'$\theta$', fontsize=16)
ax_mean.set_ylabel(r"$\xi(\theta)$",fontsize=16)    
ax_mean.errorbar(obs_correlation[0:theta_bins,0]-3.0, obs_correlation[0:theta_bins,1], obs_correlation[0:theta_bins,2],label="Hayashino et al 2004",elinewidth=3.0,fmt="o-")
ax_mean.set_title("highest density region")
ax_mean.legend(shadow=False)
fig1.savefig("mean_12.png")
        
        
               
        
