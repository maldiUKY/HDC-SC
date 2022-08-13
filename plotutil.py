#©2021 Qi Sun

import scanpy as sc
import matplotlib.pyplot as plt
from collections import Counter

diseasedict={'covid':['COVID_A40','COVID_A41','COVID_A42','COVID_A43','COVID_A44'],
             'fibB':['FibrosisB_A6','FibrosisB_A9','FibrosisB_A13','FibrosisB_B5','FibrosisB_B9'],
             'fibC':['FibrosisC_A4','FibrosisC_A5','FibrosisC_B2','FibrosisC_B5'],
             'flu':['Flu_A16','Flu_A18','Flu_A19','Flu_A20','Flu_A21'],
             'normal':['normal_A3','normal_B4','normal_A6']}


def changelabel(adata,color):
    labelcountdict={}
    labelcountliststr=[]
    labellist = adata.obs[color].tolist()
    for i in labellist:
        if i not in labelcountdict.keys():
            labelcountdict[i]=0
        labelcountdict[i]=labelcountdict[i]+1
    labelcountlist=sorted(labelcountdict.items(), key=lambda x: x[1], reverse=True)
    for (label,count) in labelcountlist:
        labelcountliststr.append(label+':'+str(count))
    return labelcountliststr


def plotsingleglycan(modeldir,fildid,figdir,glycan):
    color=glycan
    adata = sc.read(modeldir + fildid + 'DE.h5ad')  # note modelname
    ad=sc.read(modeldir + fildid + '.h5ad')
    adata.obsm=ad.obsm
    markerrange = adata[:, adata.var_names == color].X.tolist()
    minv = np.quantile(markerrange, 0.01)
    maxv = np.quantile(markerrange, 0.99)
    fig, axs = plt.subplots(1, 2, figsize=(12, 7))
    sc.pl.umap(
        adata, color=[color],
        vmin=minv,
        vmax=maxv,
        ax=axs[0],
        show=False,
        cmap='RdYlBu_r',
    )
    sc.pl.spatial(
        adata,
        img_key=None,
        vmin=minv,
        vmax=maxv,
        library_id=None,
        color=color,
        title=str(color),
        size=20,
        show=False,
        cmap = 'RdYlBu_r',
        ax=axs[1],
    )
    plt.savefig(figdir + 'DE_' + fildid + '_RdYlBu_' + color + '.png')



def plotsingle(modeldir,sampleid,figdir):
    color='leiden'
    adata = sc.read(modeldir+sampleid+'.h5ad')
    print(sampleid)

    print(adata)
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))

    colorlist=adata.obs[color].unique().tolist()

    fig, axs = plt.subplots(1, 2, figsize=(12, 7))
    sc.pl.umap(
        adata, legend_loc='on data',color=[color], palette=[
            v
            for k, v in clusters_colors.items()
            if k in colorlist],
        ax=axs[0],
        show=False,
    )

    labelcountlist=changelabel(adata,color)
    supstr=', '.join(labelcountlist)
    sc.pl.spatial(
        adata,
        img_key=None,
        #legend_loc='on data',
        #legend_fontsize=3,
        library_id=None,
        color=color,
        title=sampleid,
        size=20,
        palette=[
            v
            for k, v in clusters_colors.items()
            if k in colorlist
        ],
        ax=axs[1],
        show=False
    )
    print(supstr)
    fig.suptitle(supstr)
    plt.savefig(figdir+sampleid+'_'+color+'.png')


def plotsinglematrix(modeldir,sampleid,figdir,matrixlist):
    color='leiden'
    adata = sc.read(modeldir+sampleid+'.h5ad')
    print(sampleid)
    print(adata)

    allcolor=adata.obs[color].unique().tolist()

    tissuecolorlist=[]
    for i in allcolor:
        if i not in matrixlist:
            tissuecolorlist.append(i)

    indicetissue = []
    indice = []

    for feat in tissuecolorlist:
        ad = adata[adata.obs[color] == feat, :]
        indicetissue.extend(ad.obs.index)
    adtissue=adata[indicetissue, :]

    for feat in matrixlist:
        ad = adata[adata.obs[color] == feat, :]
        indice.extend(ad.obs.index)
    admatrix=adata[indice, :]
    admatrix.obs[color] = 'matrix'

    adata=adtissue.concatenate(admatrix)

    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    clusters_colors_grey={}
    for cluster,colorname in clusters_colors.items():
        if cluster in tissuecolorlist:
            clusters_colors_grey[cluster] = colorname
    clusters_colors_grey['matrix'] = 'lightgrey'

    fig, axs = plt.subplots(1, 2, figsize=(12, 7))
    sc.pl.umap(
        adata, legend_loc='on data',color=[color], palette=[
            v
            for k, v in clusters_colors_grey.items()
            if k in adata.obs[color].unique().tolist()],
        ax=axs[0],
        show=False
    )

    labelcountlist=changelabel(adata,color)
    supstr=', '.join(labelcountlist)
    sc.pl.spatial(
        adtissue,
        img_key=None,
        #legend_loc='on data',
        #legend_fontsize=3,
        library_id=None,
        color=color,
        title=sampleid,
        size=20,
        palette=[
            v
            for k, v in clusters_colors_grey.items()
            if k in adtissue.obs[color].unique().tolist()
        ],
        ax=axs[1],
        show=False
    )
    print(supstr)
    fig.suptitle(supstr)
    plt.savefig(figdir+sampleid+'_'+color+'_matrix.png')



def plotsinglegrey_name(modeldir,sampleid,figdir,clusterlist,clustername):
    color='leiden'
    adata = sc.read(modeldir+sampleid+'.h5ad')
    print(sampleid)
    print(adata)

    allcolor=adata.obs[color].unique().tolist()

    greyclusterlist=[]
    for clusterid in allcolor:
        if clusterid not in clusterlist:
            greyclusterlist.append(clusterid)

    adatacluster = adata[adata.obs['leiden'] == clusterlist[0], :]
    adatacluster.obs['leiden'] = clustername[clusterlist[0]]

    for i in range(1, len(clusterlist)):
        ad = adata[adata.obs['leiden'] == clusterlist[i], :]
        ad.obs['leiden'] = clustername[clusterlist[i]]
        adatacluster = adatacluster.concatenate(ad)

    indicegrey = []
    for clusterid in greyclusterlist:
        ad = adata[adata.obs[color] == clusterid, :]
        indicegrey.extend(ad.obs.index)
    adgrey = adata[indicegrey, :]
    adgrey.obs[color] = 'NA'

    adatacluster = adatacluster.concatenate(adgrey)
    adata=adatacluster

    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    clusters_colors_grey=[]
    for clusterid,colorname in clusters_colors.items():
        if clusterid in clusterlist:
            clusters_colors_grey.append((clustername[clusterid],colorname))
    clusters_colors_grey.append(('NA','lightgrey'))

    m=[v for k, v in clusters_colors_grey if k in adata.obs[color].unique().tolist()]

    fig, axs = plt.subplots(1, 2, figsize=(14, 7))
    sc.pl.umap(
        adata, legend_loc='on data',color=[color], palette=[
            v
            for k, v in clusters_colors_grey
            if k in adata.obs[color].unique().tolist()],
        ax=axs[0],
        show=False
    )

    labelcountlist=changelabel(adata,color)
    supstr=', '.join(labelcountlist)
    sc.pl.spatial(
        adata,
        img_key=None,
        #legend_loc='best',
        legend_fontsize=8,
        library_id=None,
        color=color,
        title=None,
        size=20,
        palette=[
            v
            for k, v in clusters_colors_grey
            if k in adata.obs[color].unique().tolist()
        ],
        ax=axs[1],
        show=False
    )
    print(supstr)
    #fig.suptitle(supstr)
    plt.savefig(figdir+sampleid+'_'+color+'_ROIclustername.png')

def plotsinglegrey(modeldir,sampleid,figdir,colorlist,appendix):
    color='leiden'
    adata = sc.read(modeldir+sampleid+'.h5ad')
    print(sampleid)
    print(adata)

    allcolor=adata.obs[color].unique().tolist()

    greycolorlist=[]
    for i in allcolor:
        if i not in colorlist:
            greycolorlist.append(i)

    indicegrey = []
    indice = []

    for feat in greycolorlist:
        ad = adata[adata.obs[color] == feat, :]
        indicegrey.extend(ad.obs.index)
    adgrey=adata[indicegrey, :]
    adgrey.obs[color]='NA'

    for feat in colorlist:
        ad = adata[adata.obs[color] == feat, :]
        indice.extend(ad.obs.index)
    ad=adata[indice, :]

    adata=adgrey.concatenate(ad)

    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    clusters_colors_grey={}
    for cluster,colorname in clusters_colors.items():
        if cluster in colorlist:
            clusters_colors_grey[cluster] = colorname
    clusters_colors_grey['NA'] = 'lightgrey'

    fig, axs = plt.subplots(1, 2, figsize=(12, 7))
    sc.pl.umap(
        adata, legend_loc='on data',color=[color], palette=[
            v
            for k, v in clusters_colors_grey.items()
            if k in adata.obs[color].unique().tolist()],
        ax=axs[0],
        show=False
    )

    labelcountlist=changelabel(adata,color)
    supstr=', '.join(labelcountlist)
    sc.pl.spatial(
        adata,
        img_key=None,
        #legend_loc='on data',
        #legend_fontsize=3,
        library_id=None,
        color=color,
        title=sampleid,
        size=20,
        palette=[
            v
            for k, v in clusters_colors_grey.items()
            if k in adata.obs[color].unique().tolist()
        ],
        ax=axs[1],
        show=False
    )
    print(supstr)
    fig.suptitle(supstr)
    #plt.savefig(figdir+sampleid+'_'+color+'greymatrix.png')
    plt.savefig(figdir + sampleid + '_' + color + appendix+'.png')

def violinplot(adata,color):
    sc.pl.violin(adata, color, rotation=30,groupby='spatial')


def umapplotDE(modeldir,fildid,color):
    adata = sc.read(modeldir + fildid + '.h5ad')  # note modelname
    print(adata)
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    #sc.pp.neighbors(adata)
    #sc.tl.umap(adata)

    sc.pl.umap(
        adata, legend_loc='on data', color=[color,'library_id'], palette=[
            v
            for k, v in clusters_colors.items()
            if k in adata.obs.leiden.unique().tolist()],
    )

import  numpy as np
def spatialplotDE(modeldir,fildid,color,figdir,sampleids,disease):
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    adata = sc.read(modeldir + fildid + 'DE.h5ad')  # note modelname
    #sampleids=diseasedict[disease]#note modification
    adata = adata[adata.obs['library_id'].isin(sampleids), :]
    libs=adata.obs.library_id.unique().tolist()

    markerrange = adata[:, adata.var_names == color].X.tolist()
    minv = np.quantile(markerrange,0.01)
    maxv = np.quantile(markerrange,0.99)

    fig, axs = plt.subplots(1, len(libs), figsize=(10, 5))
    for i, library in enumerate(libs):
        ad = adata[adata.obs.library_id == library, :].copy()
        #ad = ad[ad.obs['leiden']=='15', :]
        sc.pl.spatial(
            ad,
            img_key=None,
            vmin=minv,
            vmax=maxv,
            #legend_fontsize=0,
            library_id=None,
            color=color,
            title=libs[i],
            size=20,
            palette=[
                v
                for k, v in clusters_colors.items()
                if k in ad.obs.leiden.unique().tolist()
            ],
            show=False,
            ax=axs[i]
        )
    plt.savefig(figdir+'DE_'+disease+'_'+color+'.png')



def spatialsingleplotDE(modeldir,fildid,color,figdir):
    adata = sc.read(modeldir + fildid + 'DE.h5ad')  # note modelname
    adata = adata[adata.obs['leiden'].isin(['7','9']), :]
    markerrange = adata[:, adata.var_names == color].X.tolist()
    minv = np.quantile(markerrange,0.01)
    maxv = np.quantile(markerrange,0.99)
    fig, axs = plt.subplots(1, 1, figsize=(10, 5))
    sc.pl.spatial(
        adata,
        img_key=None,
        vmin=minv,
        vmax=maxv,
        library_id=None,
        color=color,
        title=fildid+'_'+str(color),
        size=20,
        show=False
    )
    plt.savefig(figdir+'DE_'+fildid+'_'+color+'.png')


def umapplot(modeldir,sampleid,figdir):
    #adata = sc.read(modeldir + sampleid + '_multimodels0.5.h5ad')#note modelname
    #adata = sc.read(modeldir + sampleid + '0.5.h5ad')  # note modelname
    adata = sc.read(modeldir + sampleid+'1.h5ad')
    print(sampleid)
    figurename=figdir+sampleid+'umap.png'
    print(figurename)
    print(adata)
    color='leiden'
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    sc.pl.umap(#legend_loc='on data',
        adata,  color=[color,'library_id'], palette=[
            v
            for k, v in clusters_colors.items()
            if k in adata.obs.leiden.unique().tolist()],

    )

def spatialplot(modeldir,fildid,figdir,sampleids,disease):
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    #adata = sc.read(modeldir + fildid + '_multimodels0.5.h5ad')#note modelname
    #adata = sc.read(modeldir + fildid + '0.5.h5ad')  # note modelname
    adata = sc.read(modeldir + fildid + '1.h5ad')  # note modelname
    #sampleids=diseasedict[disease]
    adata = adata[adata.obs['library_id'].isin(sampleids), :]

    libs=adata.obs.library_id.unique().tolist()

    fig, axs = plt.subplots(1, len(libs), figsize=(10, 5),constrained_layout=True)
    color = 'leiden'
    for i, library in enumerate(libs):
        ad = adata[adata.obs.library_id == library, :].copy()
        sc.pl.spatial(
            ad,
            img_key=None,
            #vmin=0,
            #vmax=0.14,
            #legend_fontsize=0,
            library_id=None,
            color=color,
            title=libs[i],
            size=20,
            palette=[
                v
                for k, v in clusters_colors.items()
                if k in ad.obs.leiden.unique().tolist()
            ],
            show=False,
            ax=axs[i]
        )
    #plt.savefig(figdir+fildid+'_'+disease+'_'+color+'0.5.png')
    #plt.savefig(figdir + fildid + '_' + disease + '_' + color + '0.5test.png')
    plt.savefig(figdir + fildid + '_' + disease + '_' + color + 'R1.png')


def spatialplotgrey(modeldir,fildid,figdir,sampleids,disease,colorlist):
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    clusters_colors_grey = {}
    for cluster, colorname in clusters_colors.items():
        if cluster not in colorlist:
            clusters_colors_grey[cluster] = 'lightgrey'
        else:
            clusters_colors_grey[cluster] = colorname
    adata = sc.read(modeldir + fildid + '0.5.h5ad')  # note modelname


    adata = adata[adata.obs['library_id'].isin(sampleids), :]

    libs=adata.obs.library_id.unique().tolist()

    fig, axs = plt.subplots(1, len(libs), figsize=(10, 5),constrained_layout=True)
    color = 'leiden'
    for i, library in enumerate(libs):
        ad = adata[adata.obs.library_id == library, :].copy()
        sc.pl.spatial(
            ad,
            img_key=None,
            #vmin=0,
            #vmax=0.14,
            #legend_fontsize=0,
            library_id=None,
            color=color,
            title=libs[i],
            size=20,
            palette=[
                v
                for k, v in clusters_colors_grey.items()
                if k in ad.obs.leiden.unique().tolist()
            ],
            show=False,
            ax=axs[i]
        )
    plt.savefig(figdir + fildid + '_' + disease + '_' + color + 'grey.png')



def getReverselist(labels,labellist):
    notlabels=[]
    for label in labels:
        if label not in labellist:
            notlabels.append(label)
    return notlabels

def spatialmatrixplot(modeldir,fildid,figdir,disease,labellist):
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    adata = sc.read(modeldir + fildid + '.h5ad')#note modelname
    color = 'leiden'
    if disease in ['covid','fibB','fibC','flu','normal']:
        sampleids=diseasedict[disease]
        adata = adata[adata.obs['library_id'].isin(sampleids), :]
        libs = adata.obs.library_id.unique().tolist()
        fig, axs = plt.subplots(1, len(libs), figsize=(10, 5))
        for i, library in enumerate(libs):
            allnode=adata.obs.size
            ad = adata[adata.obs.library_id == library, :].copy()
            ad = ad[ad.obs['leiden'].isin(labellist), :]
            selnode=ad.obs.size
            percent=round(float(selnode)/allnode,2)
            sc.pl.spatial(
                ad,
                img_key=None,
                vmin=0,
                vmax=0.14,
                palette=[
                    v
                    for k, v in clusters_colors.items()
                    if k in ad.obs.leiden.unique().tolist()],
                library_id=None,
                color=color,
                title=libs[i]+' matrix percent: '+str(percent),
                size=20,
                show=False,
                ax=axs[i]
            )
    else:
        fig, axs = plt.subplots(1, 1, figsize=(10, 5))
        allnode = adata.obs.size
        ad = adata[adata.obs['leiden'].isin(labellist), :]
        selnode = ad.obs.size
        percent = round(float(selnode) / allnode, 2)
        sc.pl.spatial(
            ad,
            img_key=None,
            vmin=0,
            vmax=0.14,
            palette=[
                v
                for k, v in clusters_colors.items()
                if k in ad.obs.leiden.unique().tolist()],
            library_id=None,
            color=color,
            title=disease+' matrix percent: '+str(percent)+ ' '+str(selnode)+'/'+str(allnode),
            size=20,
            show=False
        )

    plt.savefig(figdir+fildid+'_'+disease+'_'+color+'R1.png')


def spatialIntegratedplot(modeldir,fildid,figdir,sampleid):
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    adata = sc.read(modeldir + fildid + '.h5ad')
    adata = adata[adata.obs['library_id']==sampleid, :]

    fig, axs = plt.subplots(1, 2, figsize=(12, 7))
    sc.pl.umap(
        adata, legend_loc='on data', color=["leiden"], palette=[
            v
            for k, v in clusters_colors.items()
            if k in adata.obs.leiden.unique().tolist()],
        ax=axs[0],
        show=False,
    )

    labelcountlist = changelabel(adata)
    supstr = ', '.join(labelcountlist)
    sc.pl.spatial(
        adata,
        img_key=None,
        # legend_loc='on data',
        # legend_fontsize=3,
        library_id=None,
        color='leiden',
        title=sampleid,
        size=20,
        palette=[
            v
            for k, v in clusters_colors.items()
            if k in adata.obs.leiden.unique().tolist()
        ],
        ax=axs[1],
        show=False
    )
    print(supstr)
    fig.suptitle(supstr)
    # plt.show()
    plt.savefig(figdir + fildid+'_'+sampleid + '.png')
