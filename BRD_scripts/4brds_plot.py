
# coding: utf-8
#author: Xiaox
# In[1]:

#This notebook is created to look at RMSF differences b/w simulated and crystal structures.
get_ipython().magic(u'load_ext autoreload')
get_ipython().magic(u'autoreload 2')
get_ipython().magic(u'matplotlib inline')


# In[2]:

get_ipython().magic(u'run /Users/xiaox/PycharmProjects/BRD/ipn_env/ipn_head.py')


# In[6]:

os.chdir('/Users/xiaox/Documents/0_Projects/0_BRD/7_5BRDs/New/')


# In[5]:

import sys
sys.path.append('/Library/Python/2.7/site-packages/')
import MDAnalysis as mda


# In[47]:

#rmsf
os.chdir('/Users/xiaox/Documents/0_Projects/0_BRD/7_5BRDs/New/rmsf-rmsd/')
AT=np.loadtxt("ATAD2_run2_Ca_RMSF.dat",skiprows=29)
BA=np.loadtxt("BAZ2B_run2_Ca_RMSF.dat",skiprows=20)
B2=np.loadtxt("BRD2_1_run1_Ca_RMSF.dat",skiprows=25)
#B4=np.loadtxt("BRD4_1_run1_Ca_RMSF.dat",skiprows=40)
CBP=np.loadtxt("CREBBP_run1_Ca_RMSF.dat",skiprows=26)
#%%
print len(AT),len(BA),len(B2),len(CBP)#,len(B4)

# BA Wpf 20, AT 29, B2 25, B4 40, CBP 26
# BA Wpf 20, AT 29, B2 25, B4 40, CBP 26
font = {'family' : 'serif',
            'color'  : 'darkred',
            'weight' : 'normal',
            'size'   : 24,
            }
plt.gcf()
plt.rc('font', size='20')
plt.figure(figsize=(10, 5))
plt.xlabel('Residue', fontdict=font)
plt.ylabel('RMSF($\AA$)', fontdict=font)
plt.xlim(xmin=97, xmax=170)
#plt.grid()
plt.plot(range(97,97+len(AT[:73,1])),AT[:73,1],"#009999",label='ATAD2',linewidth=2) 
plt.plot(range(97,97+len(BA[:73,1])),BA[:73,1],"#A0A0A0",label='BAZ2B',linewidth=2)
plt.plot(range(97,97+len(B2[:73,1])),B2[:73,1],"#FF99CC",label='BRD2(1)',linewidth=2)
#plt.plot(range(97,97+len(B4[:73,1])),B4[:73,1],"#FF007F",label='BRD4(1)',linewidth=2)
plt.plot(range(97,97+len(CBP[:73,1])),CBP[:73,1],"#FFB266",label='CREBBP',linewidth=2)
plt.axvspan(0+97,2+97, edgecolor='#E5FFCC',facecolor='#E5FFCC', alpha=0.5)
plt.axvspan(5+97,7+97, edgecolor='#E5FFCC',facecolor='#E5FFCC', alpha=0.5)
plt.axvspan(14+97,16+97, edgecolor='#E5FFCC',facecolor='#E5FFCC', alpha=0.5)
plt.axvspan(23+97,25+97, edgecolor='#E5FFCC',facecolor='#E5FFCC', alpha=0.5)
plt.axvspan(58+97,59+97, edgecolor='#E5FFCC',facecolor='#E5FFCC', alpha=0.5)
plt.axvline(x=54+97,linewidth=2, color='#E5FFCC')
plt.axvline(x=51+97, linewidth=2, color='#E5FFCC')
#index = np.array( (1,6,15,24,51,54,58) )
#plt.xticks( index, ('WPF','PVD','PDY','PMD','V1','N1','YN'),)
plt.legend(fontsize=20)
plt.tight_layout()
#plt.show()
plt.savefig("../RMSF_4brds.png", dpi=400)


# In[21]:

#bottleneck
os.chdir('/Users/xiaox/Documents/0_Projects/0_BRD/7_5BRDs/New/caver/')

AT=np.loadtxt("Caver_ATAD2.txt",skiprows=1)
BA=np.loadtxt("Caver_BAZ2B.txt",skiprows=1)
B2=np.loadtxt("Caver_BRD2.txt",skiprows=1)
#B4=np.loadtxt("Caver_BRD4.txt",skiprows=1)
CBP=np.loadtxt("Caver_CREBBP.txt",skiprows=1)
print len(AT[0]),len(BA),len(B2),len(CBP)#len(B4),
for i in range(1,2001):
    if AT[i,0] != AT[i-1,0] +1:
        AT=np.insert(AT, i, np.array((AT[i-1,0] +1, 0)), 0)
    if BA[i,0] != BA[i-1,0] +1:
        BA=np.insert(BA, i, np.array((BA[i-1,0] +1, 0)), 0)
    if B2[i,0] != B2[i-1,0] +1:
        B2=np.insert(B2, i, np.array((B2[i-1,0] +1, 0)), 0)
    #if B4[i,0] != B4[i-1,0] +1:
        #B4=np.insert(B4, i, np.array((B4[i-1,0] +1, 0)), 0)
    if CBP[i,0] != CBP[i-1,0] +1:
        CBP=np.insert(CBP, i, np.array((CBP[i-1,0] +1, 0)), 0)
#np.savetxt("test.txt",B4)
#%%
font = {'family' : 'serif',
            'color'  : 'darkred',
            'weight' : 'normal',
            'size'   : 24,
            }

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 5))
plt.rc('font', size='15')
# generate some random test data
all_data=[]
all_data.append(AT[:,1])
all_data.append(BA[:,1])
all_data.append(B2[:,1])
#all_data.append(B4[:,1])
all_data.append(CBP[:,1])
# plot violin plot
ax.violinplot(all_data,
                   showmeans=False,
                   showmedians=True)
ax.set_title('Bottleneck')
ax.yaxis.grid(True)
ax.set_xticks([y+1 for y in range(len(all_data))])
ax.set_xlabel('Protein', fontdict=font)
ax.set_ylabel('Radius($\AA$)', fontdict=font)

# add x-tick labels
plt.gcf()
plt.setp(ax, xticks=[y+1 for y in range(len(all_data))],
         xticklabels=['ATAD2', 'BAZ2B', 'BRD2(1)', 'CREBBP'])
plt.legend(fontsize=10)
plt.tight_layout()
#plt.show()
plt.savefig("../4brds_Caver_violin.png", dpi=400)
#%%


# In[27]:

plt.rc('font', size='8')
plt.figure(figsize=(12, 10))
x=np.array(range(len(AT[:,1])))/2
#plt.grid()
fig,(ax0, ax1, ax2, ax3) = plt.subplots(nrows=4, sharex=True)
ax0.plot(x,AT[:,1],"#009999",label='ATAD2',marker='.',linestyle = 'None',markersize=1.8) 
ax0.legend(fontsize=8)
ax0.set_ylim(ymax=4.0)
ax1.plot(x,BA[:,1],"#A0A0A0",label='BAZ2B',marker='.',linestyle = 'None',markersize=1.8)
ax1.legend(fontsize=8)
ax1.set_ylim(ymax=4.0)
ax2.plot(x,B2[:,1],"#FF99CC",label='BRD2(1)',marker='.',linestyle = 'None',markersize=1.8)
ax2.legend(fontsize=8)
ax2.set_ylim(ymax=4.0)
#ax3.plot(x,B4[:,1],"#FF007F",label='BRD4(1)',marker='.',linestyle = 'None',markersize=1.8)
#ax3.legend(fontsize=8)
#ax3.set_ylim(ymax=4.0)
ax3.plot(x,CBP[:,1],"#FFB266",label='CREBBP',marker='.',linestyle = 'None',markersize=1.8)
ax3.legend(fontsize=8)
ax3.set_ylim(ymax=4.0)
#fig.legend(fontsize=8)
plt.tight_layout()
ax.set_title('Bottleneck')
plt.xlabel('Time(ns)', fontdict=font)
plt.ylabel('Radius($\AA$)', fontdict=font)
plt.savefig("../4brds_bottleneck.png", dpi=400)


# In[31]:

os.chdir('/Users/xiaox/Documents/0_Projects/0_BRD/7_5BRDs/New/energy/')
#%% 

AT=np.loadtxt("A2_ene.xvg",skiprows=21)
BA=np.loadtxt("BA2_ene.xvg",skiprows=21)
B2=np.loadtxt("B21_ene.xvg",skiprows=21)
#B4=np.loadtxt("B41_ene.xvg",skiprows=21)
CBP=np.loadtxt("C1_ene.xvg",skiprows=21)
font = {'family' : 'serif',
            'color'  : 'darkred',
            'weight' : 'normal',
            'size'   : 16,
            }

x=np.array(range(len(AT[:2001,3])))/2
#plt.grid()
plt.figure(figsize=(12, 12))
fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=4, sharex=True)
ax0.plot(x,AT[:2001,3]/100000,"#009999",label='ATAD2',marker='.',linestyle = 'None',markersize=1.8) 
#ax0.legend(fontsize=8)
#ax0.set_ylim(ymax=4.0)
ax1.plot(x,BA[:2001,3]/100000,"#A0A0A0",label='BAZ2B',marker='.',linestyle = 'None',markersize=1.8)
#ax1.legend(fontsize=8)
#ax1.set_ylim(ymax=4.0)
ax2.plot(x,B2[:2001,3]/100000,"#FF99CC",label='BRD2(1)',marker='.',linestyle = 'None',markersize=1.8)
#ax2.legend(fontsize=8)
#ax2.set_ylim(ymax=4.0)
#ax3.plot(x,B4[:2001,3]/100000,"#FF007F",label='BRD4(1)',marker='.',linestyle = 'None',markersize=1.8)
#ax3.legend(fontsize=8)
#ax3.set_ylim(ymax=4.0)
ax3.plot(x,CBP[:2001,3]/100000,"#FFB266",label='CREBBP',marker='.',linestyle = 'None',markersize=1.8)
#ax4.legend(fontsize=8)
#ax4.set_ylim(ymax=4.0)
#fig.legend(fontsize=8)
plt.rc('font', size='8')
plt.xlabel('Time(ns)', fontdict=font)
plt.ylabel('$(\degree)$', fontdict=font)
plt.tight_layout()
plt.savefig("../4brds_ene_totalenergy_ts.png", dpi=400)


# In[37]:

os.chdir('/Users/xiaox/Documents/0_Projects/0_BRD/7_5BRDs/New/dihedral/')

AT=np.loadtxt("ATAD2_run2_waterHB_Dih_V2.dat",skiprows=1)
BA=np.loadtxt("BAZ2B_run2_waterHB_Dih_V2.dat",skiprows=1)
B2=np.loadtxt("BRD2_1_run1_waterHB_Dih_V2.dat",skiprows=1)
#B4=np.loadtxt("BRD4_1_run1_waterHB_Dih_V2.dat",skiprows=1)
CBP=np.loadtxt("CREBBP_run1_waterHB_Dih_V2.dat",skiprows=1)
#%%
font = {'family' : 'serif',
            'color'  : 'darkred',
            'weight' : 'normal',
            'size'   : 16,
            }

x=np.array(range(len(AT[:2001,3])))/2
#plt.grid()
fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=4, sharex=True)
ax0.plot(x,AT[:2001,3],"#009999",label='ATAD2',marker='.',linestyle = 'None',markersize=1.8) 
#ax0.legend(fontsize=8)
ax0.set_ylim((-180,180))
ax1.plot(x,BA[:2001,3],"#A0A0A0",label='BAZ2B',marker='.',linestyle = 'None',markersize=1.8)
#ax1.legend(fontsize=8)
ax1.set_ylim((-180,180))
ax2.plot(x,B2[:2001,3],"#FF99CC",label='BRD2(1)',marker='.',linestyle = 'None',markersize=1.8)
#ax2.legend(fontsize=8)
ax2.set_ylim((-180,180))
#ax3.plot(x,B4[:2001,3],"#FF007F",label='BRD4(1)',marker='.',linestyle = 'None',markersize=1.8)
#ax3.legend(fontsize=8)
#ax3.set_ylim((-180,180))
ax3.plot(x,CBP[:2001,3],"#FFB266",label='CREBBP',marker='.',linestyle = 'None',markersize=1.8)
#ax4.legend(fontsize=8)
ax3.set_ylim((-180,180))
#fig.legend(fontsize=8)
plt.rc('font', size='8')
plt.xlabel('Time(ns)', fontdict=font)
plt.ylabel('$(\degree)$', fontdict=font)
plt.tight_layout()
plt.savefig("../4brds_dihedral.png", dpi=400)


# In[44]:

os.chdir('/Users/xiaox/Documents/0_Projects/0_BRD/7_5BRDs/New/rmsf-rmsd/')
#%% 

AT=np.loadtxt("ATAD2_run2_RMSD.dat",skiprows=1)
BA=np.loadtxt("BAZ2B_run2_RMSD.dat",skiprows=1)
B2=np.loadtxt("BRD2_1_run1_RMSD.dat",skiprows=1)
B4=np.loadtxt("BRD4_1_run1_RMSD.dat",skiprows=1)
CBP=np.loadtxt("CREBBP_run1_RMSD.dat",skiprows=1)
font = {'family' : 'serif',
            'color'  : 'darkred',
            'weight' : 'normal',
            'size'   : 16,
            }

plt.rc('font', size='6')

x=np.array(range(len(AT[:2001,1])))/2
#plt.grid()
fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=4, sharex=True)
ax0.plot(x,AT[:2001,1],"#009999",label='ATAD2',marker='.',linestyle = 'None',markersize=1.8) 
#ax0.legend(fontsize=8)
ax0.set_ylim(ymax=4.0)
ax1.plot(x,BA[:2001,1],"#A0A0A0",label='BAZ2B',marker='.',linestyle = 'None',markersize=1.8)
#ax1.legend(fontsize=8)
ax1.set_ylim(ymax=4.0)
ax2.plot(x,B2[:2001,1],"#FF99CC",label='BRD2(1)',marker='.',linestyle = 'None',markersize=1.8)
#ax2.legend(fontsize=8)
ax2.set_ylim(ymax=4.0)
#ax3.plot(x,B4[:2001,1],"#FF007F",label='BRD4(1)',marker='.',linestyle = 'None',markersize=1.8)
#ax3.legend(fontsize=8)
#ax3.set_ylim(ymax=4.0)
ax3.plot(x,CBP[:2001,1],"#FFB266",label='CREBBP',marker='.',linestyle = 'None',markersize=1.8)
#ax4.legend(fontsize=8)
ax3.set_ylim(ymax=4.0)
#fig.legend(fontsize=8)
plt.xlabel('Time(ns)', fontdict=font)
plt.tight_layout()
plt.savefig("../4brds_RMSD_Ca_ts.png", dpi=400)


# In[53]:

os.chdir('/Users/xiaox/Documents/0_Projects/0_BRD/7_5BRDs/New/distance/')

means_dis = (9.6,9.2,9.2,9.1, 11.8,9.5,11.5,10.4, 16.3,12.6,17.0,13.2, 12.8,12.3,14.5,10.6,17.3,13.6,16.7,13.7, 9.6,7.5,7.2,9.2, 14.0,13.1,11.7,15.0, 13.1,6.5,8.2,7.5)
std_dis = (0.8,1.8,0.7,0.4, 1.7,1.5,0.7,1.2, 2.8,1.8,1.2,2.1, 2.3,2.4,1.2,2.4,3.3,2.5,1.5,1.7, 1.9,1.5,0.9,1.3, 3.1,2.0,1.0,2.3, 1.5,0.5,0.6,0.8) # 1000ns trajectories


means_dis2 = (9.91,8.83,9.28,9.36, 10.94,10.73,10.16,10.33, 14.56,13.76,14.33,13.67, 12.93,14.21,10.97,11.07, 15.61,16.26,13.18,13.03, 8.37,8.33,7.99,7.98, 11.44,11.26,11.42,11.36, 12.08,6.42,7.96,7.00)
std_dis2 = (0.22,0.15,0.07,0.15, 0.21,0.19,0.31,0.26, 0.43,0.37,0.49,0.70, 0.54,0.47,0.32,0.79,0.72,0.59,0.40,0.58, 0.14,0.18,0.19,0.18, 0.16,0.19,0.10,0.16, 0.10,0.12,0.17,0.29) # crystal structure survey


n_groups = len( means_dis)

fig = plt.figure( figsize=(20,8)); 
ax = fig.add_subplot(111);

index = np.array( (1.2,1.4,1.6,1.8, 2.2,2.4,2.6,2.8, 3.2,3.4,3.6,3.8, 4.2,4.4,4.6,    4.8, 5.2,5.4,5.6,5.8,  6.2,6.4,6.6,6.8,  7.2,7.4,7.6,7.8, 8.2,8.4,8.6,8.8) )

bar_width = 0.2
ax.set_xlim( (min(index)-0.3, max(index)+bar_width+0.3 ));

opacity = 0.4
error_config = {'ecolor': '0.3'}

font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 29,
            }

rects1 = ax.bar(index, means_dis, bar_width,
                 alpha=opacity,
                 color='c',
                 yerr=std_dis,
                 capsize=10,
                 error_kw = {'ecolor':'c', 'elinewidth':'4', 'capthick':'3'})
rects1 = ax.bar(index, means_dis2, bar_width,
                 alpha=opacity+0.2,
                 color='k',
                 yerr=std_dis2,
                 capsize=10,
                 error_kw = {'ecolor':'k', 'elinewidth':'4', 'capthick':'3'})

plt.ylabel('Distance($\AA$)',fontdict=font)
#plt.xlabel( 'pDy-Yn    I94-Yn    I94-HC1    L92-wPf    L92-HC1    pVd-wPf    pVd-Yn    wPf-HC1', fontdict=font)
plt.xlabel( 'pDy-Yn    I94-Yn    I94-HC1    L92-wPf    L92-HC1    pVd-wPf    pVd-Yn    wPf-HC1', fontdict=font)
plt.xticks( index, ('ATAD2','BAZ2B','BRD2(1)','CREBBP','ATAD2','BAZ2B','BRD2(1)','CREBBP',    'ATAD2','BAZ2B','BRD2(1)','CREBBP','ATAD2','BAZ2B','BRD2(1)','CREBBP','ATAD2','BAZ2B','BRD2(1)','CREBBP',    'ATAD2','BAZ2B','BRD2(1)','CREBBP','ATAD2','BAZ2B','BRD2(1)','CREBBP','ATAD2','BAZ2B','BRD2(1)','CREBBP',    ), rotation=60)

plt.tick_params(axis='y', length=15, width=1.5, pad=6, labelsize=22 );
plt.tick_params(axis='x', length=0, width=1.5, pad=6, labelsize=22 );
for each_spine in ax.spines.values():
    each_spine.set_linewidth(1.5);
#plt.legend()
plt.ylim(ymin=4.5,ymax=24)
#plt.grid(which='minor')
plt.tight_layout()
#plt.show()
plt.savefig('../4brds_distance_errorbar.png',dpi=400)


# In[57]:

os.chdir('/Users/xiaox/Documents/0_Projects/0_BRD/7_5BRDs/New/distance/') 
#%%
means_dis = (8.06,8.48,9.17,7.38, 8.29,7.88,7.30,7.99, 9.34,8.24,8.57,9.08, 5.21,5.17,5.01,5.02,7.34,7.01,6.67,6.72, 11.70,11.07,9.92,12.21, 10.69,10.33,9.79,10.18)
std_dis = (1.94,2.11,0.81,0.86, 1.29,1.26,0.67,0.59, 1.33,0.84,0.69,0.69, 0.32,0.29,0.33,0.28,0.39,0.38,0.31,0.29, 2.06,0.98,0.89,1.07, 1.30,1.32,0.67,0.58) #1000ns_trajectories

means_dis2 = (7.96,7.50,8.33,7.84, 6.96,6.88,7.21,7.16, 8.38,7.80,8.64,8.53, 5.12,4.90,4.72,4.79,6.77,6.65,6.30,6.40, 9.68,9.56,9.44,9.72, 9.28,9.20,9.54,9.40)
std_dis2 = (0.12,0.22,0.13,0.26, 0.09,0.07,0.07,0.09, 0.10,0.13,0.13,0.12, 0.14,0.07,0.07,0.06,0.04,0.05,0.04,0.07, 0.07,0.07,0.08,0.14, 0.08,0.07,0.08,0.15) # Holo-crystal

n_groups = len( means_dis)

fig = plt.figure( figsize=(20,8)); 
ax = fig.add_subplot(111);

index = np.array( (1.2,1.4,1.6,1.8, 2.2,2.4,2.6,2.8, 3.2,3.4,3.6,3.8, 4.2,4.4,4.6,    4.8, 5.2,5.4,5.6,5.8,  6.2,6.4,6.6,6.8,  7.2,7.4,7.6,7.8) )

bar_width = 0.2
ax.set_xlim( (min(index)-0.3, max(index)+bar_width+0.3 ));

opacity = 0.4
error_config = {'ecolor': '0.3'}

font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 29,
            }

rects1 = ax.bar(index, means_dis, bar_width,
                 alpha=opacity,
                 color='c',
                 yerr=std_dis,
                 capsize=10,
                 error_kw = {'ecolor':'c', 'elinewidth':'4', 'capthick':'3'})
rects1 = ax.bar(index, means_dis2, bar_width,
                 alpha=opacity+0.2,
                 color='k',
                 yerr=std_dis2,
                 capsize=10,
                 error_kw = {'ecolor':'k', 'elinewidth':'4', 'capthick':'3'})

plt.ylabel('Distance($\AA$)',fontdict=font)
#plt.xlabel( 'pdY_Oh-wPf_O pdY_Oh-pMd_O pdY_Oh-V1_O N1_ND2-V1_O pMd_N-V1_O preP_O-V1_O pMd_N-wPf_O', fontdict=font)
plt.xticks( index, ('ATAD2','BAZ2B','BRD2(1)','CREBBP','ATAD2','BAZ2B','BRD2(1)','CREBBP',    'ATAD2','BAZ2B','BRD2(1)','CREBBP','ATAD2','BAZ2B','BRD2(1)','CREBBP','ATAD2','BAZ2B','BRD2(1)','CREBBP',    'ATAD2','BAZ2B','BRD2(1)','CREBBP','ATAD2','BAZ2B','BRD2(1)','CREBBP'), rotation=60)

plt.tick_params(axis='y', length=15, width=1.5, pad=6, labelsize=22 );
plt.tick_params(axis='x', length=0, width=1.5, pad=6, labelsize=22 );
for each_spine in ax.spines.values():
    each_spine.set_linewidth(1.5);
#plt.legend()
plt.ylim(ymin=4.5,ymax=14)
#plt.grid(which='minor')
plt.tight_layout()
#plt.show()
plt.savefig('../4brds_distance_Water.png',dpi=400)


# In[60]:

os.chdir('/Users/xiaox/Documents/0_Projects/0_BRD/7_5BRDs/New/dihedral/')

AT=np.loadtxt("ATAD2_run2_waterHB_Dih_V2.dat",skiprows=1)
BA=np.loadtxt("BAZ2B_run2_waterHB_Dih_V2.dat",skiprows=1)
B2=np.loadtxt("BRD2_1_run1_waterHB_Dih_V2.dat",skiprows=1)
#B4=np.loadtxt("BRD4_1_run1_waterHB_Dih_V2.dat",skiprows=1)
CBP=np.loadtxt("CREBBP_run1_waterHB_Dih_V2.dat",skiprows=1)
print len(AT),len(BA),len(B2),len(CBP)
font = {'family' : 'serif',
            'color'  : 'darkred',
            'weight' : 'normal',
            'size'   : 24,
            }

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 5))

plt.rc('font', size='20')
# generate some random test data
all_data=[]
all_data.append(AT[:2000,3]+ np.less(AT[:2000,3],0)*360)
all_data.append(BA[:2000,3]+ np.less(BA[:2000,3],0)*360)
all_data.append(B2[:2000,3]+ np.less(B2[:2000,3],0)*360)
#all_data.append(B4[:2000,3]+ np.less(B4[:2000,3],0)*360)
all_data.append(CBP[:2000,3]+ np.less(CBP[:2000,3],0)*360)

ax.violinplot(all_data,
                   showmeans=False,
                   showmedians=True)
#ax.set_title('pdY_O Dihedral')
ax.yaxis.grid(True)
ax.set_xticks([y+1 for y in range(len(all_data))])
ax.set_xlabel('Protein', fontdict=font)
ax.set_ylabel('$(\degree)$', fontdict=font)
plt.ylim(ymin=-1,ymax=361)
# add x-tick labels
plt.gcf()
plt.setp(ax, xticks=[y+1 for y in range(len(all_data))],
         xticklabels=['ATAD2', 'BAZ2B', 'BRD2(1)', 'CREBBP'])
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig("../4brds_pdY_O.png", dpi=400)


# In[71]:

# figure 1 rmsf bfactor pca 
# 1.1 rmsf
os.chdir('/Users/xiaox/Documents/0_Projects/0_BRD/7_5BRDs/otherBRDs/apo_crys/')
from glob import glob
def bfac2rmsf(crysfile):
    # crystal rmsf
    
    pdb=mda.Universe(crysfile,crysfile)
    Mc_sel=pdb.selectAtoms('name CA or name O or name N or name C')
    Bkb_sel=pdb.selectAtoms('name CA or name N or name C')
    Ca_sel=pdb.selectAtoms('name CA')
    
    Mc_Bfac=Mc_sel.atoms.bfactors # rmsf^2=3B/(8*pi^2)
    Mc_rmsf=np.sqrt(3. * Mc_Bfac / 8.) / np.pi /10 #convert Bfac into rmsf scale, unit: nm
    
    Bkb_Bfac=Bkb_sel.atoms.bfactors
    Bkb_rmsf=np.sqrt(3. * Bkb_Bfac / 8.) / np.pi /10
    
    Ca_Bfac=Ca_sel.atoms.bfactors
    Ca_rmsf=np.sqrt(3. * Ca_Bfac / 8.) / np.pi /10
    return Mc_rmsf, Bkb_rmsf, Ca_rmsf
from glob import glob

ofile_Mc=file( 'Mc_rmsf_crys.dat', 'w')
ofile_Bkb=file( 'Bkb_rmsf_crys.dat', 'w')
ofile_Ca=file( 'Ca_rmsf_crys.dat', 'w')

for pro in ['ATAD2','BAZ2B','BRD2(1)','CREBBP']:
    files = glob(str(pro)+'_crys.pdb')
    
    for pdbfile in files:
        Mc_pdb, Bkb_pdb, Ca_pdb = bfac2rmsf(pdbfile)
        Mc_str=''
        for each in Mc_pdb:
            Mc_str += str(each)
            Mc_str +=' '
        Bkb_str=''
        for each in Bkb_pdb:
            Bkb_str += str(each)
            Bkb_str +=' '
        Ca_str=''
        for each in Ca_pdb:
            Ca_str += str(each)
            Ca_str +=' '
        ofile_Mc.write(str(pdbfile)[:-4]+' '+str(Mc_str)+'\n')
        ofile_Bkb.write(str(pdbfile)[:-4]+' '+str(Bkb_str)+'\n')
        ofile_Ca.write(str(pdbfile)[:-4]+' '+str(Ca_str)+'\n')
ofile_Bkb.close()
ofile_Ca.close()
ofile_Mc.close()


# In[92]:

proall=['ATAD2','BAZ2B','BRD2(1)','CREBBP']
for pro in proall:
    Ca_MD=np.loadtxt('../rmsf_Ca/'+pro+'_Ca_rmsf.xvg', skiprows=12,usecols=(1,))
    Ca_crys=np.loadtxt(pro.split('_')[0]+'_Ca.dat')
    Ca_ensemble=np.load('../../New/ensemble/'+pro+'_pdb_rmsf.npy')/10
    plt.gcf()
    plt.rc('font', size='12')
    plt.title(pro,fontsize='18')
    plt.plot(Ca_MD,'b.-',label='MD')
    plt.plot(Ca_crys,'m.-',label='B-factor')
    plt.plot(Ca_ensemble,'y.-',label='Ensemble')
    plt.legend()
    plt.xlabel('Residue ID',fontsize=18)
    plt.ylabel('RMSF(nm)',fontsize=18)
    plt.ylim([0,0.6])
    #plt.show()
    plt.savefig('../../New/'+pro+'rmsf.png',dpi=400)
    plt.clf()


# In[ ]:



