import numpy as np
import pandas as pd
import subprocess
import seird
import tool 
import enkf
import param as par
import matplotlib.pyplot as plt
plt.switch_backend('agg')

# 
#exp = 'test'
#exp = 'hdx_enkf_nomu3'
exp = 'hdx_enkf_mu0.01'
#exp = 'hdx_enkf_mu0.1'
outDir = '/home/okazaki/tmp/covid19/out/'+exp+'/'
figDir = '/home/okazaki/tmp/covid19/fig/'+exp+'/'
dateStart = '2020-01-01'
dateEnd = '2020-06-11'
dates = np.arange(dateStart,dateEnd,dtype='datetime64[D]')
df1 = pd.read_csv('/home/okazaki/tmp/covid19/data/org/owid-covid-data.csv', index_col=0)
countryList = set(df1.index.values)
df2 = pd.read_csv('/home/okazaki/tmp/covid19/data/org/time_series_covid19_confirmed_global_iso3_regions.csv', index_col=146)
countryList2 = set(df2.index.values)
keyList = ['population','population_density','median_age','aged_65_older','aged_70_older','gdp_per_capita']
#countryList = ['USA']
#countryList = ['CHN']
#countryList = ['AUT']
#countryList = ['HKG']
#countryList = ['NZL']
#countryList = ['ITA','FRA','DEU','ESP','NZL','TWN','HKG','CHN','USA']
countryList = ['ITA','USA','TUR']


# mkdir
subprocess.call('mkdir -p '+outDir+'/npy/',shell=True)
subprocess.call('mkdir -p '+outDir+'/txt/',shell=True)
subprocess.call('mkdir -p '+figDir,shell=True)

# for text output
with open(outDir +'/txt/alldata.txt',mode='w') as f:
    outStr = 'country,lon,lat,'
    for key in keyList :
        outStr += key +','
    outStr += 'date_start,date_end,lambda,delta\n'
    f.write(outStr)


for country in countryList :

    # check name
    if pd.isnull(country) :
        continue
    if not country in countryList2 :
        print( country + ' not found --skipped')
        continue

    print('---------------------------')
    print(country)

    # initialize vars
    N = tool.getdata(country,'population')
    xf, xa, pa, pat = enkf.ensinit(N)
    sprdPat = np.std(pat,ddof=1,axis=1)
    #
    xf_save = np.zeros((len(dates),6,par.kmax),dtype=float)
    xa_save = np.zeros((len(dates),6,par.kmax),dtype=float)
    pa_save = np.full ((len(dates),5,par.kmax),-999.)
    dx_save = np.zeros((len(dates),5,par.kmax),dtype=float)
    yo_save = np.zeros((len(dates),par.nobs),dtype=float)
    omb_save = np.zeros((len(dates),par.nobs),dtype=float)
    Y       = np.empty((par.nobs,par.kmax),dtype=float)
    print('mean',np.mean(pa,axis=1))
    print('std',np.std(pa,axis=1,ddof=1))
    print('min',np.min(pa,axis=1))
    print('max',np.max(pa,axis=1))


    # before the first case
    for it,date in enumerate(dates) :
        yo = tool.setobs_hdx(country,date)

        if yo[0] == np.nan :
            continue
        elif yo[0] > 0 :
            break
        # save
        xf_save[it,:,:] = xf

    # set value
    xf[1,:] = yo[0] / pa[0,:]
    xf[2,:] = yo[0]
    xf[0,:] = xf[0,:] - xf[1,:] - xf[2,:]
    yo_save[it,:] = yo[:]
    dateStr = date
    itStr = it

    # after first case detected
    xa = xf
    for it,date in enumerate(dates) :
        if date <= dateStr :
            continue
            
        # forecast
        xf,dxf = enkf.ensfcst(xa,pa)
        #if date == dates[-1] :
        #    tmp = np.concatenate([xf,dxf,pa])
        #    cor = np.corrcoef(tmp)
        #    plt.imshow(cor,vmin=-1,vmax=1,cmap='bwr')
        #    plt.xlim(-0.5,15.5)
        #    plt.ylim(-0.5,15.5)
        #    plt.hlines(5.5, -0.5,15.5)
        #    plt.hlines(10.5,-0.5,15.5)
        #    plt.vlines(5.5, -0.5,15.5)
        #    plt.vlines(10.5,-0.5,15.5)
        #    
        #    plt.xticks(np.arange(16),['S','E','I','R','D','N','dE','dI','dE2R','dI2R','dD',r'$\sigma$',r'$\lambda$',r'$\gamma$',r'$\delta$',r'$\mu$'])
        #    plt.yticks(np.arange(16),['S','E','I','R','D','N','dE','dI','dE2R','dI2R','dD',r'$\sigma$',r'$\lambda$',r'$\gamma$',r'$\delta$',r'$\mu$'])
        #    ax = plt.gca()
        #    labels = ax.get_xticklabels()
        #    plt.setp(labels, rotation=90)
        #    plt.colorbar()
        #    #plt.scatter(pa[3],dxf[4])
        #    #plt.scatter(pa[3],xf[2])
        #    plt.show()
        # observation
        yo = tool.setobs_hdx(country,date)
        yo0 = tool.setobs_hdx(country,dates[it-1])
        yo[:3] = yo[:3] - yo0[:3]
        yo[3] = yo0[3]
        yo = np.where(yo < 0, 0, yo)
        if np.any(yo == par.undef) :
            # no analysis
            xa = xf
            yo = 0.
        else :
            # obsope
            Y[0,:] = dxf[1,:]
            Y[1,:] = dxf[3,:]
            Y[2,:] = dxf[4,:]
            Y[3,:] = xa [2,:]
            # analysis
            dza = enkf.LETKF( np.concatenate([dxf,pat]), Y, yo )
            dxa = dza[:5,:]
            pat = dza[5:,:]
            # relaxation
            pat = enkf.relax(sprdPat,pat)
            pa  = enkf.arctransParam(pat)
            # update
            xa = enkf.ensupdate(xa,dxa)


        # save
        xf_save[it,:,:] = xf
        xa_save[it,:,:] = xa
        pa_save[it,:,:] = pa
        dx_save[it,:,:] = dx_save[it-1,:] + dxa
        yo_save[it,:]   = yo_save[it-1,:] + yo 
        omb_save[it,0]  = yo[0] - np.mean(dxf[1,:])
        omb_save[it,1]  = yo[1] - np.mean(dxf[3,:])
        omb_save[it,2]  = yo[2] - np.mean(dxf[4,:])
        omb_save[it,3]  = yo[3] - np.mean(xf [2,:])

        # monitor
        np.set_printoptions(precision=2,suppress=True)
        #print(date,np.mean(pa,axis=1)[0:5])
        print(date,np.mean(xf,axis=1)[0:5])
        print(date,np.mean(xa,axis=1)[0:5])
        print(date,np.mean(dxf,axis=1),yo)
        print(date,np.mean(dxa,axis=1),yo)
        np.set_printoptions(precision=3,suppress=False)
        print(date,np.mean(pa,axis=1))
        #print(pa)

    #
    # save data
    #
    if par.outNorm :
        xf_save /= N
        xa_save /= N
        yo_save /= N
        pa_save[:,1,:] *= N
        dx_save /= N
    np.save( outDir +'/npy/'+ country +'_xf.npy', xf_save )
    np.save( outDir +'/npy/'+ country +'_xa.npy', xa_save )
    np.save( outDir +'/npy/'+ country +'_pa.npy', pa_save )
    np.save( outDir +'/npy/'+ country +'_yo.npy', yo_save )

    #
    # text
    #
    outFile = country +'_timeseries.txt'
    with open(outDir +'/txt/'+ outFile,mode='w') as f:
        # header
        outStr='date,S,E,I,R,D,N,sigma,lambda,gamma,delta,mu\n'
        #outStr='date,S,E,I,R,D,N,sigma,lambda,gamma,delta,mu,rmse_dI,rmse_dR,rmse_dD\n'
        f.write(outStr)

        for it in np.arange(itStr+1,len(dates)) :
            date = dates[it]
            outStr = str(date).replace('-','')
            # state
            for i in range(6) :
                data = np.mean(xa_save[it,i,:],axis=0)
                outStr += ','+ str(data)
            # param
            for i in range(5) :
                data = np.mean(pa_save[it,i,:],axis=0)
                outStr += ','+ str(data)
            ## skill
            #for i in range(3) :
            #    data = np.sqrt(np.mean(oma_save * oma_save,axis=0))[i]
            #    print(data)
            #    outStr += ','+ str(data)
            outStr += '\n'

            # write
            f.write(outStr)


    outFile = 'alldata.txt'
    with open(outDir +'/txt/'+ outFile,mode='a') as f:
        # initialize var
        outStr = country +','

        # read obs
        data = tool.getdata(country,'lon')
        outStr += str(data) +','
        data = tool.getdata(country,'lat')
        outStr += str(data) +','
        for k,key in enumerate(keyList) :
            data = df1.at[country,key][0] 
            outStr += str(data) +','

        dateStr = dates[itStr+14]
        dateEnd = dates[itStr+28]
        data = str(dateStr).replace('-','')
        outStr += str(data) +','
        data = str(dateEnd).replace('-','') 
        outStr += str(data) +','
        data = np.mean(np.mean(pa_save[itStr+14:itStr+29,:,:],axis=2),axis=0)[1] 
        outStr += str(data) +','
        data = np.mean(np.mean(pa_save[itStr+14:itStr+29,:,:],axis=2),axis=0)[3] 
        outStr += str(data) +'\n'

        # write
        f.write(outStr)

    #
    # plot
    #
    fig,ax = plt.subplots(2,1,sharex=True,figsize=(8,8))
    for k in range(par.kmax) :
        ax[0].plot(dates[itStr+1:],xa_save[itStr+1:,0,k],lw=0.5,c='lightgreen')
        ax[0].plot(dates[itStr+1:],xa_save[itStr+1:,1,k],lw=0.5,c='lightyellow')
        ax[0].plot(dates[itStr+1:],xa_save[itStr+1:,2,k],lw=0.5,c='mistyrose')
        ax[0].plot(dates[itStr+1:],xa_save[itStr+1:,3,k],lw=0.5,c='powderblue')
        ax[0].plot(dates[itStr+1:],xa_save[itStr+1:,4,k],lw=0.5,c='gainsboro')
    ax[0].plot(dates[itStr+1:],np.mean(xa_save[itStr+1:,0,:],axis=1),c='g',label='S')
    ax[0].plot(dates[itStr+1:],np.mean(xa_save[itStr+1:,1,:],axis=1),c='gold',label='E')
    ax[0].plot(dates[itStr+1:],np.mean(xa_save[itStr+1:,2,:],axis=1),c='r',label='I')
    ax[0].plot(dates[itStr+1:],np.mean(xa_save[itStr+1:,3,:],axis=1),c='b',label='R')
    ax[0].plot(dates[itStr+1:],np.mean(xa_save[itStr+1:,4,:],axis=1),c='k',label='D')
    ax[0].plot(dates[itStr+1:],np.mean(dx_save[itStr+1:,1,:],axis=1),c='r',ls='dashed',label='CI')
    ax[0].plot(dates[itStr+1:],yo_save[itStr+1:,0],c='r',marker='o',markerfacecolor='None',linewidth=0,label='CI_obs')
    ax[0].plot(dates[itStr+1:],yo_save[itStr+1:,2],c='k',marker='o',markerfacecolor='None',linewidth=0,label='D_obs')
    if par.mu != 0 :
        ax[0].plot(dates[itStr+1:],np.mean(dx_save[itStr+1:,3,:],axis=1),c='b',ls='dashed',label='I->R')
        ax[0].plot(dates[itStr+1:],yo_save[itStr+1:,1],c='b',marker='o',markerfacecolor='None',linewidth=0,label='I->R_obs')
    else :
        ax[0].plot(dates[itStr+1:],yo_save[itStr+1:,1],c='b',marker='o',markerfacecolor='None',linewidth=0,label='R_obs')
    ax[0].set_yscale('log')
    ax[0].set_title(country)
    ax[0].set_ylabel('Population')
    ax[0].legend()
    for k in range(par.kmax) :
        ax[1].plot(dates[itStr+1:],pa_save[itStr+1:,1,k],lw=0.5,c='mistyrose')
    ax2=ax[1].twinx()
    for k in range(par.kmax) :
        ax2.plot(dates[itStr+1:],pa_save[itStr+1:,3,k],lw=0.5,c='gainsboro')
    ax[1].plot(dates[itStr+1:],np.mean(pa_save[itStr+1:,1,:],axis=1),c='r',label=r'$\lambda$')
    ax2.plot(dates[itStr+1:],np.mean(pa_save[itStr+1:,3,:],axis=1),c='k',label=r'$\delta$')
    ax[1].set_ylabel(r'$\lambda$')
    ax[1].set_yscale('log')
    ax2.set_ylabel(r'$\delta$')
    ax[1].set_xlabel('Date')
    h1, l1 = ax[1].get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax2.legend(h1+h2, l1+l2)
    labels = ax[1].get_xticklabels()
    plt.setp(labels, rotation=45)
    plt.tight_layout()
    plt.show()
    plt.savefig( figDir +'/timeseries_'+ country +'.png')
    plt.close()


    plt.plot(dates[itStr+1:], np.sqrt(omb_save * omb_save)[itStr+1:,0],c='r',label='nI')
    plt.plot(dates[itStr+1:], np.sqrt(omb_save * omb_save)[itStr+1:,1],c='b',label='nI->R')
    plt.plot(dates[itStr+1:], np.sqrt(omb_save * omb_save)[itStr+1:,2],c='k',label='nD')
    #plt.plot(dates[itStr+1:], np.sqrt(omb_save * omb_save)[itStr+1:,3],c='r',ls='dashed',label='I')
    plt.ylabel('ens mean fcst RMSE')
    plt.xlabel('Date')
    ax = plt.gca()
    labels = ax.get_xticklabels()
    plt.setp(labels, rotation=45)
    plt.title(country)
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.savefig( figDir +'/rmsd_'+ country +'.png')
    plt.close()
