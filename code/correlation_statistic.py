import numpy as np
import pylab as P
import matplotlib.collections as collections
import matplotlib 
from scipy import stats
from scipy import optimize as op
import sys
#plot_dir="../paper/plots/"
plot_dir="../paper/plots/"

matplotlib.use('pdf')
#defining the power-law function to fit the correlation function
powerlaw = lambda x, amp, index: (x/amp)**(index)
fitfunc = lambda p, x: (x/p[0])**(p[1])
errfunc = lambda p, x, y,err: (y - fitfunc(p, x))/err

pinit = [2.0,-1.0] 

n_models=101
n_mocks=15
theta_bins=13
pro_path="/Users/jemejia/CosmisVarianceLAES/"

obs_correlation_file=pro_path + "data/obs/hayashino_whole_SSA22_field.txt"
obs=np.loadtxt(obs_correlation_file,skiprows=4)
theta1=obs[0:theta_bins,0]
if(sys.argv[1]=="maxden"):
    obs_correlation_file=pro_path + "data/obs/hayashino_whole_SSA22_field.txt"
if(sys.argv[1]=="meanden"):
    obs_correlation_file=pro_path + "data/obs/ouchi2010.txt"


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
slope_obs=-pfinal[1]
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
model_numbers=np.arange(len(f_occ_m))

mean_correlation= np.zeros( [n_models,theta_bins] )
std_correlation = np.zeros( [n_models,theta_bins] )
chi_individual  = np.zeros( [n_models,n_mocks] )

ro_ind  = np.zeros( [n_models,n_mocks] )
ro_err_ind  = np.zeros( [n_models,n_mocks] )
slope_ind  = np.zeros( [n_models,n_mocks] )
slope_err_ind  = np.zeros( [n_models,n_mocks] )



slope = np.zeros( n_models )
slope_err = np.zeros( n_models )
ro = np.zeros( n_models )
ro_err = np.zeros( n_models )

for i in range(n_models):
    corr=np.zeros( (n_mocks,theta_bins) )
    corr_std=np.zeros( (n_mocks,theta_bins) )
    
    for j in range(n_mocks):
        model_correlation_filename=pro_path +  "data/laes/correlation/" + sys.argv[1] + "_model_" + str(i) + "_mock_" + str(j) + ".txt"
        
        lae_correlation  = np.loadtxt(model_correlation_filename)
        corr[j,:] = lae_correlation[:,1]
        corr_std[j,:] = lae_correlation[:,2]
        #print lae_correlation[:,0],
        out = op.leastsq(errfunc, pinit,args=(theta1, corr[j,:],corr_std[j,:]), full_output=1)
        pfinal = out[0]
        #pinit=pfinal
        covar = out[1]
        
        if( covar is None ): 
            ro_ind[i,j]=ro_ind[i,j-1]
            slope_ind[i,j]=slope_ind[i,j-1]
            ro_err_ind[i,j]=ro_err_ind[i,j-1]
            slope_err_ind[i,j]=slope_err_ind[i,j-1]

        else:
        
            ro_ind[i,j] = float(pfinal[0])
            slope_ind[i,j] = -float(pfinal[1])
        
            ro_err_ind[i,j]=covar[0][0] # error squared
            slope_err_ind[i,j]=covar[1][1] # error squared
            
    ro[i]=np.mean(ro_ind[i,:]) 
    ro_err[i]=np.std(ro_ind[i,:]) + np.sqrt(np.sum( ro_err_ind[i,:] ))/n_mocks 
    slope[i]=np.mean(slope_ind[i,:])
    slope_err[i]=np.std(slope_ind[i,:]) + np.sqrt(np.sum(slope_err_ind[i,:]))/n_mocks
    std_correlation[i,:]  = np.std(corr,axis=0)
    mean_correlation[i,:] = np.mean(corr,axis=0)
    print "ro=",ro[i],"+-" ,ro_err[i], "\n"
    print "slope=",slope[i],"+-" ,slope_err[i], "\n"
chisquare=np.zeros([n_models,theta_bins])
chi=np.zeros(n_models)
"""
for i in range(n_models):
    
    correlation= mean_correlation[i,:]
    std= std_correlation[i,:]
    
    out = op.leastsq(errfunc, pinit,args=(lae_correlation[:,0], correlation,std), full_output=1)
    
    pfinal = out[0]
    covar = out[1]

        
    ro[i] = float(pfinal[0])
    slope[i] = float(pfinal[1])
    
    ro_err[i]=np.sqrt(covar[0][0])
    slope_err[i]=np.sqrt(covar[1][1])
    
    print "ro=",ro[i],"+-" ,ro_err[i], "\n"
    print "slope=",slope[i],"+-" ,slope_err[i], "\n"
    
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
"""
fig2=P.figure()
ro_plot=fig2.add_subplot(111)
ro_index=np.where(slope_err<2.0)
ro_plot.errorbar(ro[ro_index],slope[ro_index], yerr=slope_err[ro_index], xerr=ro_err[ro_index], label="Models",elinewidth=1.5,fmt="o")

if sys.argv[1]=="maxden" :
    ro_plot.errorbar(ro_obs,slope_obs, yerr=slope_err_obs, xerr=ro_err_obs, label="Hayashino et al 2004",elinewidth=4.5)
if sys.argv[1]=="meanden":
    ro_plot.errorbar(ro_obs,slope_obs, yerr=slope_err_obs, xerr=ro_err_obs, label="Ouchi et al 2010",elinewidth=4.5)
P.xlim(0,50)
P.ylim(0.5,3.5)
ro_plot.set_xlabel(r'$\theta_{0}(arcsec)$', fontsize=20)
ro_plot.set_ylabel(r"$\beta$",fontsize=20)
##ro_plot.set_title("Angular Correlation parameters. Match survey",fontsize=20)
plot_name=plot_dir + "power_law_correlation_" + sys.argv[1] + ".pdf"
ro_plot.legend(shadow=False,loc=2)

fig2.savefig(plot_name,format="pdf")

P.xlim()
P.ylim()
fig3=P.figure()
mmin_plot=fig3.add_subplot(111)
mmax_index=np.where(m_max_m<12.0) #and ro<ro_max) #and slope_min<slope<slope_max)
mmin_plot.errorbar(m_min_m[mmax_index],ro[mmax_index], yerr=ro_err[mmax_index],label=r"$\log M_{max}<12.0$",elinewidth=1.5,fmt="s")
mmax_index=np.where(m_max_m>12.0) #and ro<ro_max) #and slope_min<slope<slope_max)
mmin_plot.errorbar(m_min_m[mmax_index]+0.01,ro[mmax_index], yerr=ro_err[mmax_index],label=r"$\log M_{max}>12.0$",elinewidth=1.5,fmt="o")
mmin_plot.set_ylabel(r'$\theta_{0}$', fontsize=20)
mmin_plot.set_xlabel(r"$M_{min}$",fontsize=20)

if sys.argv[1]=="maxden" :
    mmin_plot.set_title(r"$M_{min}$ vs $\theta_{0}$. Match Survey. Hayashino 2004",fontsize=20)
if sys.argv[1]=="meanden":
    mmin_plot.set_title(r"$M_{min}$ vs $\theta_{0}$. Match Survey. Ouchi 2010",fontsize=20)

P.ylim(0,100)
plot_name=plot_dir +"mmin_vs_correlation_" +sys.argv[1] + ".pdf"
mmin_plot.legend(shadow=False,loc=2)
collection = collections.BrokenBarHCollection.span_where(m_min_m, ymin=ro_min, ymax=ro_max, where=ro>0,  facecolor='red', alpha=0.5)
mmin_plot.add_collection(collection)
fig3.savefig(plot_name,format="pdf")
P.ylim()

corr_index=np.where(( (ro_max>(ro-ro_err)) & ((ro+ro_err)>ro_min) ) & ( ((slope+slope_err)>slope_min) & ((slope-slope_err)<slope_max) )== True) #and ro<ro_max) #and slope_min<slope<slope_max)
#corr_index=np.where(( (ro_max>ro) & (ro>ro_min) ) == True) #and ro<ro_max) #and slope_min<slope<slope_max)
print corr_index

m_min=m_min_m[corr_index]
m_max=m_max_m[corr_index]
f_occ=f_occ_m[corr_index]
model_numbers=model_numbers[corr_index]
print model_numbers, np.size(model_numbers)
ro_b=ro[corr_index]
slope_b=slope[corr_index]


best_index=np.where(f_occ>=0.0)
corr_stat=np.zeros( [ np.size(best_index) , 5 ] )
corr_stat[:,0]=m_min[best_index]
corr_stat[:,1]=m_max[best_index]
corr_stat[:,2]=f_occ[best_index]
corr_stat[:,3]=ro[best_index]
corr_stat[:,4]=slope[best_index]
corr_stat_file=pro_path + "data/mock_survey/" + "powerlaw_bestmodels_" + sys.argv[1] + ".dat"
np.savetxt(corr_stat_file,corr_stat,fmt='%2.3lf')

fig4=P.figure()
mass_plot=fig4.add_subplot(111)
mass_index=np.where(f_occ<=0.2)

mass_plot.plot(m_min[mass_index],m_max[mass_index] - m_min[mass_index],'bs',label=r"$\log f_{occ}<=0.2$")
mass_index=np.where(f_occ>0.2)
mass_plot.plot(m_min[mass_index],m_max[mass_index] - m_min[mass_index],'kp',label=r"$\log f_{occ}>0.2$")
#mass_plot.plot(m_min[best_index], m_max[best_index],'ks')
mass_plot.set_xlabel( r'$M_{min}$', fontsize=20 )
mass_plot.set_ylabel( r"$M_{max}$",fontsize=20 )
mass_plot.set_title(r"Observationally consistent models. Full Field",fontsize=20)
if sys.argv[1]=="maxden" :
    mass_plot.set_title(r"$\log M_{min}$ vs $\Delta log M$. Observationally consistent models. Hayashino 2004",fontsize=20)
if sys.argv[1]=="meanden":
    mass_plot.set_title(r"$\log M_{min}$ vs $\Delta log M$. Observationally consistent models. Ouchi 2010",fontsize=20)
    

P.xlim(10.2,11.3)
#P.ylim(10.4,12.0)
plot_name=plot_dir +"mmin_vs_dm_" + sys.argv[1] + ""
mass_plot.legend(shadow=False,loc=2)
fig4.savefig(plot_name,format="pdf")


fig5=P.figure()
f_occ_plot=fig5.add_subplot(111)
#mass_index=np.where(f_occ=<0.2)
f_occ_plot.plot(m_min,f_occ, 'kp')
f_occ_plot.set_xlabel(r'$M_{min}$', fontsize=20)
f_occ_plot.set_ylabel(r"$f_{occ}$",fontsize=20)
if sys.argv[1]=="maxden" :
    f_occ_plot.set_title(r"$M_{min}$ vs $f_{occ}$. Observationally consistent models. Hayashino 2004",fontsize=20)
if sys.argv[1]=="meanden":
    f_occ_plot.set_title(r"$M_{min}$ vs $f_{occ}$. Observationally consistent models. Ouchi 2010",fontsize=20)



#P.xlim(10.2,13.0)
collection = collections.BrokenBarHCollection.span_where(m_min_m, ymin=0.1, ymax=0.2, where=f_occ>0,  facecolor='blue', alpha=0.5)
f_occ_plot.add_collection(collection)
#P.ylim(0.,12.0)
plot_name=plot_dir +"mmin_vs_focc_"  + sys.argv[1] + ".pdf"
f_occ_plot.legend(shadow=False,loc=2)
fig5.savefig(plot_name,format="pdf")




               
        
