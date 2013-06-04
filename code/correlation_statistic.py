import numpy as np
import pylab as P
import matplotlib.collections as collections
from scipy import stats
from scipy import optimize as op

#defining the power-law function to fit
powerlaw = lambda x, amp, index: (x/amp)**(index)
fitfunc = lambda p, x: (x/p[0])**p[1]
errfunc = lambda p, x, y,err: (y - fitfunc(p, x))/err

pinit = [1.0, -1.0] 

n_models=101
n_mocks=12
theta_bins=13
pro_path="/home/jemejia/CosmicVarianceLAES/"

obs_correlation_file=pro_path + "data/obs/hayashino_whole_SSA22_field.txt"
obs_correlation=np.loadtxt(obs_correlation_file,skiprows=4)
theta=obs_correlation[0:theta_bins,0]
correlation=obs_correlation[0:theta_bins,1]
std=obs_correlation[0:theta_bins,2]
out = op.leastsq(errfunc, pinit,args=(theta, correlation,std), full_output=1)

pfinal = out[0]
covar = out[1]

    
print pfinal
print "covariance matriz \n\n", covar, "\n\n"

ro_obs=pfinal[0]
ro_err_obs=np.sqrt(covar[0][0])
slope_obs=pfinal[1]
slope_err_obs=np.sqrt(covar[1][1])
ro_min=ro_obs - ro_err_obs
ro_max=ro_obs + ro_err_obs
slope_min=slope_obs - slope_err_obs
slope_max=slope_obs + slope_err_obs
print "ro=",ro_obs,"+-" ,ro_err_obs, "\n"
print "slope=",slope_obs,"+-" ,slope_err_obs, "\n"

theta=obs_correlation[0:theta_bins,0]
models_file=pro_path + "data/mock_survey/" + "BestModels_match_mock.dat" 
model_info=np.loadtxt(models_file)
m_min_m=model_info[:,0] #_m refers to model
m_max_m=model_info[:,1]
f_occ_m=model_info[:,2]


mean_correlation= np.zeros( [n_models,theta_bins] )
std_correlation = np.zeros( [n_models,theta_bins] )
chi_individual  = np.zeros( [n_models,n_mocks] )
slope = np.zeros( n_models )
slope_err = np.zeros( n_models )
ro = np.zeros( n_models )
ro_err = np.zeros( n_models )

for i in range(n_models):
    corr=np.zeros( (n_mocks,theta_bins) )
    
    for j in range(n_mocks):
        model_correlation_filename=pro_path +  "data/laes/correlation/maxden_model_" + str(i) + "_mock_" + str(j) + ".txt"
        
        lae_correlation  = np.loadtxt(model_correlation_filename)
        corr[j,:] = lae_correlation[:,1]

        #slope[i,j], intercept, r_value, p_value, std_err = stats.linregress( log_theta[corr_index] , np.log10(correlation[corr_index]) )
        #log_ro[i,j]=-1*intercept
        #print log_theta, corr[j,:]
        #print slope[i,j] , np.power(10,log_ro[i,j]), r_value
        #chi_individual[i,j]= chi_individual[i,j] + (mean_correlation[i,j] - obs_correlation[j,1])*(mean_correlation[i,j] - obs_correlation[j,1])/(std_correlation[i,j]*std_correlation[i,j] + obs_correlation[j,2]*obs_correlation[j,2])
    std_correlation[i,:]  = np.std(corr,axis=0)
    mean_correlation[i,:] = np.mean(corr,axis=0)

    chisquare=np.zeros([n_models,theta_bins])
    chi=np.zeros(n_models)
for i in range(n_models):
    correlation= mean_correlation[i,:]
    std= std_correlation[i,:]
    
    out = op.leastsq(errfunc, pinit,args=(theta, correlation,std), full_output=1)
    
    pfinal = out[0]
    covar = out[1]

    
   # print pfinal
    #print "covariance matriz \n\n", covar, "\n\n"
    
    ro[i] = float(pfinal[0])
    slope[i] = float(pfinal[1])
    
    ro_err[i]=np.sqrt(covar[0][0])
    slope_err[i]=np.sqrt(covar[1][1])
    
    print "ro=",ro[i],"+-" ,ro_err[i], "\n"
    print "slope=",slope[i],"+-" ,slope_err[i], "\n"
    #print powerlaw(theta,ro[i],slope[i]), "\n\n"
    for j in range(theta_bins):
        chisquare[i,j]= chisquare[i,j] + (mean_correlation[i,j] - obs_correlation[j,1])*(mean_correlation[i,j] - obs_correlation[j,1])/(std_correlation[i,j]*std_correlation[i,j] + obs_correlation[j,2]*obs_correlation[j,2])
        chi[i]= chi[i] + chisquare[i,j] 
corr_stat_file=pro_path + "data/mock_survey/" + "correlation_stat.dat"
corr_stat=np.empty( [ n_models , 4 ] )
chi=chi
sort_index=np.argsort(chi)
corr_stat[:,0]=m_min_m[sort_index]
corr_stat[:,1]=m_max_m[sort_index]
corr_stat[:,2]=f_occ_m[sort_index]
corr_stat[:,3]=chi[sort_index]
np.savetxt(corr_stat_file,corr_stat,fmt='%2.3lf')
fig2=P.figure()
ro_plot=fig2.add_subplot(111)
ro_plot.errorbar(ro,slope, yerr=slope_err, xerr=ro_err, label="Models",elinewidth=1.5,fmt="o")
ro_plot.errorbar(ro_obs,slope_obs, yerr=slope_err_obs, xerr=ro_err_obs, label="Hayashino et al 2004",elinewidth=4.5)
ro_plot.set_xlabel(r'$\theta_{0}$', fontsize=20)
ro_plot.set_ylabel(r"$\gamma$",fontsize=20)
ro_plot.set_title("Angular Correlation parameters. Match Surveys")
plot_name="power_law_correlation.png"
ro_plot.legend(shadow=False,loc=2)


fig2.savefig(plot_name)


fig3=P.figure()
mmin_plot=fig3.add_subplot(111)
mmax_index=np.where(m_max_m<12.0) #and ro<ro_max) #and slope_min<slope<slope_max)
mmin_plot.errorbar(m_min_m[mmax_index],ro[mmax_index], yerr=ro_err[mmax_index],label=r"$\log M_{max}<12.0$",elinewidth=1.5,fmt="o")
mmax_index=np.where(m_max_m>12.0) #and ro<ro_max) #and slope_min<slope<slope_max)
mmin_plot.errorbar(m_min_m[mmax_index]+0.01,ro[mmax_index], yerr=ro_err[mmax_index],label=r"$\log M_{max}>12.0$",elinewidth=1.5,fmt="o")
mmin_plot.set_ylabel(r'$\theta_{0}$', fontsize=20)
mmin_plot.set_xlabel(r"$M_{min}$",fontsize=20)
mmin_plot.set_title(r"$M_{min}$ vs $\theta_{0}$. Match Surveys")

plot_name="mmin_vs_correlation.png"
mmin_plot.legend(shadow=False,loc=2)
collection = collections.BrokenBarHCollection.span_where(m_min_m, ymin=15, ymax=23, where=ro>0,  facecolor='red', alpha=0.5)
mmin_plot.add_collection(collection)


fig3.savefig(plot_name)
corr_index=np.where(( (ro_max>ro) & (ro>ro_min) ) & ( (slope>slope_min) & (slope<slope_max) )== True) #and ro<ro_max) #and slope_min<slope<slope_max)

print corr_index
corr_stat=np.zeros( [ np.size(corr_index) , 5 ] )
corr_stat[:,0]=m_min_m[corr_index]
corr_stat[:,1]=m_max_m[corr_index]
corr_stat[:,2]=f_occ_m[corr_index]
corr_stat[:,3]=ro[corr_index]
corr_stat[:,4]=slope[corr_index]
corr_stat_file=pro_path + "data/mock_survey/" + "powerlaw_bestmodels_2sigma.dat"
np.savetxt(corr_stat_file,corr_stat,fmt='%2.3lf')
for i in corr_index:
    fig1=P.figure()
    ax_mean=fig1.add_subplot(111)
    model='(Mmin,Mmax,focc)=({0},{1},{2})'.format(m_min_m[i], m_max_m[i], f_occ_m[i])
    ax_mean.errorbar(obs_correlation[0:theta_bins,0],mean_correlation[i,:], std_correlation[i,:], label=model,elinewidth=1.5)

    ax_mean.set_xlabel(r'$\theta$', fontsize=16)
    ax_mean.set_ylabel(r"$\xi(\theta)$",fontsize=16)    
    ax_mean.errorbar(obs_correlation[0:theta_bins,0]-3.0, obs_correlation[0:theta_bins,1], obs_correlation[0:theta_bins,2],label="Hayashino et al 2004",elinewidth=3.0,fmt="o-")
    ax_mean.plot(obs_correlation[0:theta_bins,0]-3.0, powerlaw(theta,ro[i],slope[i]), label="power law fit" )
    ax_mean.set_title("highest density region")
    ax_mean.legend(shadow=False)
    plot_name="mean_correlation_" + model + ".png"
    fig1.savefig(plot_name)

               
        
