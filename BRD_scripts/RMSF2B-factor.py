
# coding: utf-8

# # Preload

# In[1]:

#This notebook is created to look at RMSF differences b/w simulated and crystal structures.
get_ipython().magic(u'load_ext autoreload')
get_ipython().magic(u'autoreload 2')
get_ipython().magic(u'matplotlib inline')


# In[2]:

get_ipython().magic(u'run /Users/xiaox/PycharmProjects/BRD/ipn_env/ipn_head.py')


# In[3]:

os.chdir('/Users/xiaox/Documents/0_Projects/0_BRD/7_5BRDs/2000us/rmsf/')


# In[4]:

import sys
sys.path.append('/Library/Python/2.7/site-packages/')
import MDAnalysis as mda


# In[5]:

def BfComp(Mcfile,Bkbfile,Cafile,crysfile):
    import pandas as pd
    # simulation rmsf
    
    Mc=np.loadtxt(Mcfile,skiprows=12,usecols=(1,)) # load the rmsf.xvg file
    #print Mc
    Bkb=np.loadtxt(Bkbfile,skiprows=12,usecols=(1,))
    Ca=np.loadtxt(Cafile,skiprows=12,usecols=(1,))
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
    df_Mc = pd.DataFrame(np.array([Mc[:-1], Mc_rmsf]).T, columns=[ 'MchainS', 'MchainC'])
    df_Bkb = pd.DataFrame(np.array([Bkb, Bkb_rmsf]).T, columns=[ 'BkbS', 'BkbC'])
    df_Ca = pd.DataFrame(np.array([Ca, Ca_rmsf]).T, columns=[ 'CaS', 'CaC'])
#     df = pd.DataFrame(np.array([Mc, Mc_rmsf, Bkb, Bkb_rmsf, Ca, Ca_rmsf]).T, 
#                       columns=[ 'MchainS', 'MchainC', 'BkbS', 'BkbC', 'CaS', 'CaC'])
    return df_Mc,df_Bkb,df_Ca


# In[6]:

for pro in ["ATAD2","BAZ2B","BRD2_1","BRD4_1","CREBBP"]:
    
    mc=str(pro)+'-Mc_rmsf.xvg'
    bkb=str(pro)+'-Bkb_rmsf.xvg'
    ca=str(pro)+'-Ca_rmsf.xvg'
    crys=str(pro)+'_swiss.pdb'
    df_mc,df_bkb,df_ca=BfComp(mc,bkb,ca,crys)
    print str(pro)
    print df_bkb.corr()
    #print df_bkb.corr()
    #print df_ca.corr()
    df_bkb.plot()
    plt.gcf()
    plt.ylim((0,1))


# In[8]:

def BfComp_2(Mcfile,Bkbfile,Cafile,crysfile):
    
    # simulation rmsf
    
    Mc1=np.loadtxt(str(Mcfile)+'_1.xvg',skiprows=12,usecols=(1,)) # load the rmsf.xvg file
    Mc2=np.loadtxt(str(Mcfile)+'_2.xvg',skiprows=12,usecols=(1,))
    #print Mc
    Bkb1=np.loadtxt(str(Bkbfile)+'_1.xvg',skiprows=12,usecols=(1,))
    Bkb2=np.loadtxt(str(Bkbfile)+'_2.xvg',skiprows=12,usecols=(1,))
    Ca1=np.loadtxt(str(Cafile)+'_1.xvg',skiprows=12,usecols=(1,))
    Ca2=np.loadtxt(str(Cafile)+'_2.xvg',skiprows=12,usecols=(1,))
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
    df_Mc = pd.DataFrame(np.array([Mc1[:-1],Mc2[:-1], Mc_rmsf]).T, columns=[ 'MchainS1','MchainS2', 'MchainC'])
    df_Bkb = pd.DataFrame(np.array([Bkb1, Bkb2, Bkb_rmsf]).T, columns=[ 'BkbS1','BkbS2', 'BkbC'])
    df_Ca = pd.DataFrame(np.array([Ca1, Ca2, Ca_rmsf]).T, columns=[ 'CaS1', 'CaS2', 'CaC'])
#     df = pd.DataFrame(np.array([Mc, Mc_rmsf, Bkb, Bkb_rmsf, Ca, Ca_rmsf]).T, 
#                       columns=[ 'MchainS', 'MchainC', 'BkbS', 'BkbC', 'CaS', 'CaC'])
    return df_Mc,df_Bkb,df_Ca


# In[9]:

import pandas as pd
Df_mc=[]
for pro in ["ATAD2","BAZ2B","BRD2_1","BRD4_1","CREBBP"]:
    
    mc=str(pro)+'-Mc_rmsf'
    bkb=str(pro)+'-Bkb_rmsf'
    ca=str(pro)+'-Ca_rmsf'
    crys=str(pro)+'_swiss.pdb'
    df_mc,df_bkb,df_ca=BfComp_2(mc,bkb,ca,crys)
    Df_mc.append(df_mc)
    print str(pro), df_mc.shape
    print df_mc.corr()
    #print df_bkb.corr()
    #print df_ca.corr()
    #df_mc.plot()
    #plt.gcf()
    #plt.ylim((0,1))


# In[29]:

all_avg=np.loadtxt("Mc_rmsf_crys_avg.txt",skiprows=2)
ATAD2_avg=all_avg[:,0]
BAZ2B_avg=all_avg[:412,2]
BRD2_1_avg=all_avg[:440,4]
BRD4_1_avg=all_avg[:504,6]
CREBBP_avg=all_avg[:440,8]
ATAD2_sq=np.sqrt(np.loadtxt("ATAD2_pro_pca_sqfluct.txt"))/10
BAZ2B_sq=np.sqrt(np.loadtxt("BAZ2B_pro_pca_sqfluct.txt"))/10
BRD2_1_sq=np.sqrt(np.loadtxt("BRD2_1_pro_pca_sqfluct.txt"))/10
BRD4_1_sq=np.sqrt(np.loadtxt("BRD4_1_pro_pca_sqfluct.txt"))/10
CREBBP_sq=np.sqrt(np.loadtxt("CREBBP_pro_pca_sqfluct.txt"))/10

font = {'family' : 'serif',
            'color'  : 'darkred',
            'weight' : 'normal',
            'size'   : 18,
            }
print np.corrcoef(ATAD2_avg,Df_mc[0]['MchainS1'][:])
print np.corrcoef(ATAD2_avg,Df_mc[0]['MchainS2'][:])
plt.gcf()
plt.plot(np.array(range(len(ATAD2_avg))),ATAD2_avg,'k',label='ATAD2_avg')
plt.plot(np.array(range(len(ATAD2_avg))),Df_mc[0]['MchainC'][:],label='ATAD2_init')
plt.plot(np.array(range(len(ATAD2_avg))),Df_mc[0]['MchainS1'][:],label='ATAD2_run1')
plt.plot(np.array(range(len(ATAD2_avg))),Df_mc[0]['MchainS2'][:],label='ATAD2_run2')
plt.xlabel('Residue', fontdict=font)
#plt.ylabel('RMSF($\AA$)', fontdict=font)
plt.ylabel('RMSF(nm)', fontdict=font)
#plt.plot(ATAD2_sq,'y',label='ATAD2_sq')
plt.xlim((0,len(ATAD2_avg)))
plt.ylim((0,1))
plt.legend()
#plt.show()
plt.tight_layout()
plt.savefig('./atad2_rmsf.png', dpi=300)
plt.clf()


# In[17]:

ATAD2_avg


# In[43]:

print np.corrcoef(BAZ2B_avg,Df_mc[1]['MchainS1'][:])
print np.corrcoef(BAZ2B_avg,Df_mc[1]['MchainS2'][:])

plt.plot(BAZ2B_avg,'k',label='BAZ2B_avg')
plt.plot(Df_mc[1]['MchainC'][:],label='BAZ2B_init')
plt.plot(Df_mc[1]['MchainS1'][:],label='BAZ2B_run1')
plt.plot(Df_mc[1]['MchainS2'][:],label='BAZ2B_run2')
#plt.plot(BAZ2B_sq,'y',label='BAZ2B_sq')
plt.ylim((0,1))
plt.legend()
#plt.show()
plt.savefig('./5baz2b_rmsf.png', dpi=300)
plt.clf()
print np.corrcoef(BRD2_1_avg,Df_mc[2]['MchainS1'][16:])
print np.corrcoef(BRD2_1_avg,Df_mc[2]['MchainS2'][16:])

plt.plot(BRD2_1_avg,'k',label='BRD2(1)_avg')
plt.plot(Df_mc[2]['MchainC'][16:],label='BRD2(1)_init')
plt.plot(Df_mc[2]['MchainS1'][16:],label='BRD2(1)_run1')
plt.plot(Df_mc[2]['MchainS2'][16:],label='BRD2(1)_run2')
#plt.plot(BRD2_1_sq,'y',label='BRD2(1)_sq')
plt.ylim((0,1))
plt.legend()
#plt.show()
plt.savefig('./brd2_rmsf.png', dpi=300)
plt.clf()
print np.corrcoef(BRD4_1_avg,Df_mc[3]['MchainS1'][4:])
print np.corrcoef(BRD4_1_avg,Df_mc[3]['MchainS2'][4:])

plt.plot(BRD4_1_avg,'k',label='BRD4(1)_avg')
plt.plot(Df_mc[3]['MchainC'][4:],label='BRD4(1)_init')
plt.plot(Df_mc[3]['MchainS1'][4:],label='BRD4(1)_run1')
plt.plot(Df_mc[3]['MchainS2'][4:],label='BRD4(1)_run2')
#plt.plot(BRD4_1_sq,'y',label='BRD4(1)_sq')
plt.ylim((0,1))
plt.legend()
#plt.show()
plt.savefig('./brd4_rmsf.png', dpi=300)
plt.clf()
print np.corrcoef(CREBBP_avg,Df_mc[4]['MchainS1'][8:-8])
print np.corrcoef(CREBBP_avg,Df_mc[4]['MchainS2'][8:-8])

plt.plot(CREBBP_avg,'k',label='CREBBP_avg')
plt.plot(Df_mc[4]['MchainC'][8:-8],label='CREBBP_init')
plt.plot(Df_mc[4]['MchainS1'][8:-8],label='CREBBP_run1')
plt.plot(Df_mc[4]['MchainS2'][8:-8],label='CREBBP_run2')
#plt.plot(CREBBP_sq,'y',label='CREBBP_sq')
plt.ylim((0,1))
plt.legend()
#plt.show()
plt.ylim((0,1))
plt.savefig('./cbp_rmsf.png', dpi=300)
plt.clf()


# In[30]:


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


# In[29]:

os.chdir('/Users/xiaox/Documents/0_Projects/0_BRD/2_pdb_new_15.4.15/Holo/pro_ONLY_PDB/')


# In[31]:

#read crystal structure, save the B-factors to dat file, by processing with excel, get .csv file.

from glob import glob
ofile_Mc=file( 'Mc_rmsf_crys.dat', 'w')
ofile_Bkb=file( 'Bkb_rmsf_crys.dat', 'w')
ofile_Ca=file( 'Ca_rmsf_crys.dat', 'w')

for pro in ["ATAD2","BAZ2B","BRD2-1","BRD4-1","CREBBP"]:
    files = glob(str(pro)+'*')
    
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


# # Subfamilies and the initial structures

# In[23]:

os.chdir('/Users/xiaox/Documents/0_Projects/0_BRD/7_5BRDs/otherBRDs/apo_crys/')
from glob import glob


# In[39]:

ls *_Mc.dat


# In[27]:

ls *_Mc_rmsf.xvg


# In[64]:

#Note that this list is changed according to the crystal structure names.
# divide the prolists into subfamilies.
#Subfamily I
proall=[]
proall.append(['PCAF_run2','CECR2_run4','GCN5L2_run2','FALZ_run2'])
proall.append(['BRD2-1_run1_','BRD2-1_run2_','BRD2-2_run3','BRD3-1_run3','BRD3-2_run3',               'BRD4-1_run1_','BRD4-1_run3_','BRD4-2_run2','BRDT-1','BRDT-1_run2'])
proall.append(['WDR9-2_run2','CREBBP_run1_','CREBBP_run2_','EP300_run3','PHIP-2_run3'])
proall.append(['KIAA1240','BRD1','BRD9_run2','BRD7_run3','ATAD2_run2_','BRPF1B_run2','BRPF1B_run3']) #BRD7, BRPF1B, NMR
proall.append(['TIF1Apdb_run3','BAZ2B_run2_','BAZ2B_run3_'])
proall.append(['TAF1-2_run3','TAF1L_2_run3','TAF1-1'])
proall.append(['SMARCA4','PB1','PB2','PB3','PB4','PB5']) #,'SMARCA2_run2'])


# In[56]:

for pro in proall[0]:
    Mc_MD=np.loadtxt(pro+'_Mc_rmsf.xvg', skiprows=12,usecols=(1,))
    Mc_crys=np.loadtxt(pro.split('_')[0]+'_Mc.dat')
    plt.title(pro)
    plt.plot(Mc_MD,'b.-')
    plt.plot(Mc_crys,'m.-')
    print len(Mc_MD), len(Mc_crys)
    plt.show()


# In[57]:

for pro in proall[1]:
    Mc_MD=np.loadtxt(pro+'_Mc_rmsf.xvg', skiprows=12,usecols=(1,))
    Mc_crys=np.loadtxt(pro.split('_')[0]+'_Mc.dat')
    plt.title(pro)
    plt.plot(Mc_MD,'b.-')
    plt.plot(Mc_crys,'m.-')
    print len(Mc_MD), len(Mc_crys)
    plt.show()


# In[58]:

for pro in proall[2]:
    Mc_MD=np.loadtxt(pro+'_Mc_rmsf.xvg', skiprows=12,usecols=(1,))
    Mc_crys=np.loadtxt(pro.split('_')[0]+'_Mc.dat')
    plt.title(pro)
    plt.plot(Mc_MD,'b.-')
    plt.plot(Mc_crys,'m.-')
    print len(Mc_MD), len(Mc_crys)
    plt.show()


# In[59]:

for pro in proall[3]:
    Mc_MD=np.loadtxt(pro+'_Mc_rmsf.xvg', skiprows=12,usecols=(1,))
    Mc_crys=np.loadtxt(pro.split('_')[0]+'_Mc.dat')
    plt.title(pro)
    plt.plot(Mc_MD,'b.-')
    plt.plot(Mc_crys,'m.-')
    print len(Mc_MD), len(Mc_crys)
    plt.show()


# In[60]:

for pro in proall[4]:
    Mc_MD=np.loadtxt(pro+'_Mc_rmsf.xvg', skiprows=12,usecols=(1,))
    Mc_crys=np.loadtxt(pro.split('_')[0]+'_Mc.dat')
    plt.title(pro)
    plt.plot(Mc_MD,'b.-')
    plt.plot(Mc_crys,'m.-')
    print len(Mc_MD), len(Mc_crys)
    plt.show()


# In[62]:

for pro in proall[5]:
    Mc_MD=np.loadtxt(pro+'_Mc_rmsf.xvg', skiprows=12,usecols=(1,))
    Mc_crys=np.loadtxt(pro.split('_')[0]+'_Mc.dat')
    plt.title(pro)
    plt.plot(Mc_MD,'b.-')
    plt.plot(Mc_crys,'m.-')
    print len(Mc_MD), len(Mc_crys)
    plt.show()


# In[65]:

for pro in proall[6]:
    Mc_MD=np.loadtxt(pro+'_Mc_rmsf.xvg', skiprows=12,usecols=(1,))
    Mc_crys=np.loadtxt(pro.split('_')[0]+'_Mc.dat')
    plt.title(pro)
    plt.plot(Mc_MD,'b.-')
    plt.plot(Mc_crys,'m.-')
    print len(Mc_MD), len(Mc_crys)
    plt.show()


# In[ ]:



