import scanpy as sc
import pandas
import numpy as np

def getcolors(df,featcluster,clusters_colors):
    rowcolorlist=[]
    colcolorlist=[]
    for row in df.index:
        rowcolorlist.append(clusters_colors[row])
    for col in df.columns:
        colcolorlist.append(clusters_colors[featcluster[col]])
    return rowcolorlist,colcolorlist

def readmatrix():
    matrixlist=[]
    f=open('./data/matrix/matrix.txt')
    for line in f:
        matrixlist.append(line.strip())
    return matrixlist

def DEana(fildid,clusterid,topN,featlist):
    matrixlist=readmatrix()
    adata = sc.read(modeldir + fildid + 'DE.h5ad')
    adata = adata[adata.obs['leiden'].isin(featlist), :]
    METHOD = "wilcoxon"
    sc.tl.rank_genes_groups(adata, method=METHOD, groups=[clusterid], groupby='leiden', key_added='de')
    result = adata.uns['de']
    groups = [clusterid]
    sta = pandas.DataFrame(
        {group + '_' + key: result[key][group]
         for group in groups for key in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']})
    name=clusterid + '_names'
    inverse_boolean_series = ~sta[name].isin(matrixlist)
    sta = sta[inverse_boolean_series]
    tops=sta[name].tolist()[:topN]
    sta.to_csv('./heatmaptmp/' + fildid +'-'+clusterid+ '.csv', sep=',')
    return tops

def convertrowcolor(rowcolorlist):
    colorlist=[]
    for color in rowcolorlist:
        if color not in colorlist:
            colorlist.append(color)
    return colorlist

def DEana_foldchange(fildid,clusterid,topN,featlist):
    matrixlist=readmatrix()
    adata = sc.read(modeldir + fildid + 'DE.h5ad')
    adata = adata[adata.obs['leiden'].isin(featlist), :]
    METHOD = "wilcoxon"
    sc.tl.rank_genes_groups(adata, method=METHOD, tie_correct=True,groups=[clusterid], groupby='leiden', key_added='de')
    result = adata.uns['de']
    groups = [clusterid]
    sta = pandas.DataFrame(
        {group + '_' + key: result[key][group]
         for group in groups for key in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']})
    name=clusterid + '_names'
    foldchange=clusterid + '_scores'#zhuyi shi zhege!!!
    inverse_boolean_series = ~sta[name].isin(matrixlist)
    sta = sta[inverse_boolean_series]
    sta=sta.sort_values(by=foldchange, ascending=False).head(topN)
    #sta = sta[sta[clusterid+'_logfoldchanges'] > 0.5] #buyong,buyongha
    tops=sta[name].tolist()
    sta.to_csv('./heatmaptmp/' + fildid +'-'+clusterid+ '.csv', sep=',')
    return tops

def getOriData(modeldir,fildid):
    adata = sc.read(modeldir + fildid + 'DE.h5ad')
    print(adata)
    pd = pandas.DataFrame(adata.X, columns=adata.var_names.tolist(), index=adata.obs['leiden'].tolist())
    for x in pd.columns:
        pd[x] = pd[x] - pd[x].mean()
    #pd=pd.sample(frac=0.3)
    #pd=pd.sample(axis='columns',frac=0.7)
    #adata=sc.AnnData(pd)
    sc.pp.subsample(adata, fraction=0.05, random_state=1, copy=False)
    sc.pp.highly_variable_genes(adata,n_top_genes=40)
    adata = adata[:, adata.var.highly_variable]
    print(adata)
    markerrange = adata.X.tolist()
    from itertools import chain
    markerrange = list(chain.from_iterable(markerrange))
    minv = np.quantile(markerrange, 0.01)
    maxv = np.quantile(markerrange, 0.99)

    g = sc.pl.heatmap(adata, adata.var_names, groupby='spatial', dendrogram=False, swap_axes=True, vmin=minv, vmax=maxv,
                      show=False,
                      cmap='RdYlBu_r', save='_' + fildid + '_raw.png')

    import seaborn as sns
    import matplotlib.pyplot as plt
    #g = sns.heatmap(pd.transpose(), cmap='RdYlBu_r',vmin=minv, vmax=maxv)
    #plt.show()

    return adata

def genHeatmapdata(figdir,fildid,featlist,topN):
    tops = []
    topsdict = {}
    featcluster = {}

    clusters_colors_dict = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    clusters_colors = {}
    for cluster, colorname in clusters_colors_dict.items():
        if cluster in featlist:
            clusters_colors[cluster] = colorname
    clusters_colors['double'] = 'black'

    for clusterid in featlist:
        if fildid == 'brainfrontal' or fildid == 'brainhippo' or fildid == 'normal':
            seltops = DEana_foldchange('normal', clusterid, topN, featlist)
        else:
            seltops = DEana_foldchange(fildid, clusterid, topN, featlist)
        for feat in seltops:
            if feat not in featcluster.keys():
                featcluster[feat] = clusterid
            else:
                featcluster[feat] = 'double'
        topsdict[clusterid] = seltops

    for feat in featlist:
        seltops=topsdict[feat]
        for selfeat in seltops:
            if selfeat not in tops:
                tops.append(selfeat)

    if fildid == 'brainfrontal' or fildid == 'brainhippo' or fildid == 'normal':
        adatalabel = sc.read(modeldir + 'normal' + '0.5.h5ad')
        adata = sc.read(modeldir + 'normal' + 'DE.h5ad')
    else:
        adatalabel = sc.read(modeldir + fildid + '.h5ad')
        adata = sc.read(modeldir + fildid + 'DE.h5ad')
    adata.obs['leiden'] = adatalabel.obs['leiden']

    indice = []
    for feat in featlist:
        ad = adata[adata.obs['leiden'] == feat, :]
        #index = random.choices(ad.obs.index, k=100)
        index = ad.obs.index
        indice.extend(index)

    adata = adata[indice, :]

    colindice = []
    for topfeat in tops:
        ad = adata[:, adata.var_names == topfeat]
        colindice.extend(ad.var_names)
    adata = adata[:, colindice]

    pd = pandas.DataFrame(adata.X, columns=adata.var_names.tolist(), index=adata.obs['leiden'].tolist())
    for x in pd.columns:
        pd[x] = pd[x] - pd[x].mean()
    adata.X = pd
    rowcolorlist, colcolorlist = getcolors(pd, featcluster, clusters_colors)
    colorlsit=convertrowcolor(rowcolorlist)

    adata.uns['leiden_colors']=colorlsit
    return adata,pd,rowcolorlist, colcolorlist,clusters_colors



def genHeatmapdata_name(figdir,fildid,clusterlist,clustername,topN):
    tops = []
    topsdict = {}
    featcluster = {}

    clusters_colors_dict = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    clusters_colors = {}
    for clusterid, colorname in clusters_colors_dict.items():
        if clusterid in clusterlist:
            clusters_colors[clustername[clusterid]] = colorname

    for clusterid in clusterlist:
        if fildid == 'brainfrontal' or fildid == 'brainhippo' or fildid == 'normal':
            seltops = DEana_foldchange('normal', clusterid, topN, clusterlist)
        else:
            seltops = DEana_foldchange(fildid, clusterid, topN, clusterlist)
        for feat in seltops:
            if feat not in featcluster.keys():
                featcluster[feat] = clustername[clusterid]
            else:
                featcluster[feat] = 'double'
        topsdict[clusterid] = seltops

    for clusterid in clusterlist:
        seltops=topsdict[clusterid]
        for selfeat in seltops:
            if selfeat not in tops:
                tops.append(selfeat)

    if fildid == 'brainfrontal' or fildid == 'brainhippo' or fildid == 'normal':
        adatalabel = sc.read(modeldir + 'normal' + '0.5.h5ad')
        adata = sc.read(modeldir + 'normal' + 'DE.h5ad')
    else:
        adatalabel = sc.read(modeldir + fildid + '.h5ad')
        adata = sc.read(modeldir + fildid + 'DE.h5ad')
    adata.obs['leiden'] = adatalabel.obs['leiden']

    adatacluster = adata[adata.obs['leiden'] == clusterlist[0], :]
    adatacluster.obs['leiden'] = clustername[clusterlist[0]]

    for i in range(1,len(clusterlist)):
        ad = adata[adata.obs['leiden'] == clusterlist[i], :]
        ad.obs['leiden'] = clustername[clusterlist[i]]
        adatacluster = adatacluster.concatenate(ad)

    adata = adatacluster

    colindice = []
    for topfeat in tops:
        ad = adata[:, adata.var_names == topfeat]
        colindice.extend(ad.var_names)
    adata = adata[:, colindice]

    pd = pandas.DataFrame(adata.X, columns=adata.var_names.tolist(), index=adata.obs['leiden'].tolist())
    for x in pd.columns:
        pd[x] = pd[x] - pd[x].mean()
    adata.X = pd

    colorlsit=[]
    for clusterid in clusterlist:
        colorlsit.append(clusters_colors[clustername[clusterid]])
    adata.uns['leiden_colors']=colorlsit
    return adata


def genHeatmap(adata,figdir,fildid):
    markerrange = adata.X.tolist()
    from itertools import chain
    markerrange = list(chain.from_iterable(markerrange))
    minv = np.quantile(markerrange, 0.01)
    maxv = np.quantile(markerrange, 0.99)

    #g=sc.pl.heatmap(adata, adata.var_names, groupby='leiden', dendrogram=False, swap_axes=True, vmin=minv, vmax=maxv, show=False,
    #              cmap='RdYlBu_r',save='_'+fildid+'_name.png')

    g = sc.pl.heatmap(adata, adata.var_names, groupby='leiden', dendrogram=False, swap_axes=True, vmin=minv, vmax=maxv,
                      show=False,
                      cmap='RdYlBu_r', save='_' + fildid + '.png')

    #sc.pl.matrixplot(adata, adata.var_names, groupby='leiden', dendrogram=False, cmap='RdYlBu_r', show=False,swap_axes=True,save='_'+fildid+'_name.png')





def genclustermap(pd,rowcolorlist, colcolorlist,clusters_colors):
    import seaborn as sns


    g = sns.clustermap(pd, robust=True,
                       metric="correlation",
                       xticklabels=1,
                       figsize=(20, 20),
                       row_cluster=True, col_cluster=True,
                       row_colors=rowcolorlist,
                       col_colors=colcolorlist,
                       cmap='RdYlBu_r',
                       cbar_kws={'label': 'intensity'}, yticklabels=False)

    for label in clusters_colors.keys():
        g.ax_col_dendrogram.bar(0, 0, color=clusters_colors[label], label=label, linewidth=10)
    g.ax_col_dendrogram.legend(loc="best", ncol=5)
    g.ax_col_dendrogram.set_visible(True)
    g.ax_row_dendrogram.set_visible(True)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90,fontsize=15)
    # g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=15)

    # Adjust the postion of the main colorbar for the heatmap
    g.cax.set_position([1, .2, .03, .45])

    g.savefig(figdir + fildid + '_heatmap.png')


def genclustermap_transpose(pd,rowcolorlist, colcolorlist,clusters_colors):

    import seaborn as sns

    pd = pd.transpose()
    g = sns.clustermap(pd, robust=True,
                       metric="correlation",
                       #xticklabels=1,
                       figsize=(20, 20),
                       row_cluster=False, col_cluster=False,
                       #row_colors=colcolorlist,
                       col_colors=rowcolorlist,
                       cmap='RdYlBu_r',
                       cbar_kws={'label': 'intensity'}, xticklabels=False)

    for label in clusters_colors.keys():
        g.ax_col_dendrogram.bar(0, 0, color=clusters_colors[label], label=label, linewidth=10)
    g.ax_col_dendrogram.legend(loc="best", ncol=5)
    g.ax_col_dendrogram.set_visible(False)
    g.ax_row_dendrogram.set_visible(False)
    # g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90,fontsize=15)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=15)
    #g.yaxis.set_tick_params(horizontalalignment='left')

    # Adjust the postion of the main colorbar for the heatmap
    g.cax.set_position([1, .2, .03, .45])

    g.savefig(figdir + fildid + '_heatmap.png')

def biclustering(figdir,fildid,featlist,topN):
    tops = []
    topsdict={}
    featcluster = {}

    clusters_colors_dict = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    clusters_colors = {}
    for cluster, colorname in clusters_colors_dict.items():
        if cluster in featlist:
            clusters_colors[cluster] = colorname
    clusters_colors['double'] = 'black'

    for clusterid in featlist:
        if fildid == 'brainfrontal' or fildid == 'brainhippo' or fildid == 'normal':
            seltops=DEana_foldchange('normal',clusterid, topN,featlist)
        else:
            seltops = DEana_foldchange(fildid, clusterid, topN, featlist)
        for feat in seltops:
            if feat not in featcluster.keys():
                featcluster[feat]=clusterid
            else:
                featcluster[feat]='double'
        topsdict[clusterid]=seltops

    for feat in featlist:
        tops.extend(topsdict[feat])

    #tops=list(set(tops))

    if fildid=='brainfrontal' or fildid=='brainhippo' or fildid=='normal':
        adatalabel = sc.read(modeldir + 'normal' + '0.5.h5ad')
        adata = sc.read(modeldir + 'normal' + 'DE.h5ad')
    else:
        adatalabel = sc.read(modeldir + fildid + '.h5ad')
        adata = sc.read(modeldir + fildid + 'DE.h5ad')
    adata.obs['leiden']=adatalabel.obs['leiden']

    indice=[]
    for feat in featlist:
        ad = adata[adata.obs['leiden'] == feat, :]
        #index = random.choices(ad.obs.index, k=100)
        index=ad.obs.index
        indice.extend(index)

    adata = adata[indice, :]

    colindice=[]
    for topfeat in tops:
        ad = adata[:,adata.var_names==topfeat]
        colindice.extend(ad.var_names)
    adata = adata[:,colindice]

    pd = pandas.DataFrame(adata.X, columns=adata.var_names.tolist(), index=adata.obs['leiden'].tolist())
    for x in pd.columns:
        pd[x] = pd[x] - pd[x].mean()
    adata.X=pd

    sc.pl.heatmap(adata, adata.var_names, groupby='leiden', dendrogram=False, swap_axes=True,
                  cmap='RdYlBu_r')
    sc.pl.matrixplot(adata, adata.var_names, groupby='leiden', dendrogram=True,cmap='RdYlBu_r',swap_axes=True,)

    rowcolorlist,colcolorlist=getcolors(pd,featcluster,clusters_colors)

    print(adata)

    #pd=pd.transpose()


    import seaborn as sns
    import matplotlib
    from matplotlib.pyplot import gcf
    #cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["blue","yellow","red"])

    g = sns.clustermap(pd,robust=True,
                       metric="correlation",
                       #xticklabels=1,
                       figsize=(20,20),
                       row_cluster=False,col_cluster=False,
                       #row_colors=colcolorlist,
                       col_colors=rowcolorlist,
                       cmap='RdYlBu_r',
                       cbar_kws={'label': 'intensity'},yticklabels=True,xticklabels=False)

    for label in clusters_colors.keys():
        g.ax_col_dendrogram.bar(0, 0, color=clusters_colors[label],label=label, linewidth=10)
    g.ax_col_dendrogram.legend(loc="best",ncol=5)
    g.ax_col_dendrogram.set_visible(False)
    g.ax_row_dendrogram.set_visible(False)
    #g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90,fontsize=15)
    #g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=15)

    # Adjust the postion of the main colorbar for the heatmap
    g.cax.set_position([1, .2, .03, .45])


    g.savefig(figdir+fildid+'_heatmap.png')

fildidfeatlist={'liver':['0', '1', '6', '12'],'Flu_A19':['9','10'],'FibrosisB_A6':['0','2', '4','11'],
                'FibrosisB_B5':['0','3', '10'],'Flu_A20':['0','3'],'brainfrontal':['4','7','14'],
                'brainhippo':['6','11'],'normal':['4','6','7','11','14'],'normal_B4':['1','4','5','6'],
                'COVID_A42':['2','5','7','9'],'COVID_A43':['0','6','8','12'],'FibrosisC_A5':['5','6','10'],
                'FibrosisC_A4':['2','13']}

covidlist=['Flu_A19','Flu_A20','FibrosisB_A6','FibrosisB_B5','normal_B4','COVID_A42','COVID_A43','FibrosisC_A5','FibrosisC_A4']




clustername={'FibrosisC_A5':{'6':'DAD','10':'Mucin','5':'Fibrosis'},
             'liver':{'6':'6_Hepatocytes','0':'0_Blood vessels L','12':'12_Blood vessels C','1':'1_unknown'}}


if __name__ == "__main__":
    dir = ''
    fildid='liver'

    if fildid=='liver':
        modeldir = dir + 'models_liver/'
        figdir = dir + 'figures_liver/'

    elif fildid in covidlist:
        modeldir='../covid/'+'models_resolution1/'
        figdir = dir + 'figures_covid/'

    elif fildid=='brainfrontal' or fildid=='brainhippo' or fildid=='normal':
        modeldir = dir + 'models_AD/'
        figdir = dir + 'figures_AD/'

    clusterlist=fildidfeatlist[fildid]
    samples = []
    batch = True
    topN=10

    adata=getOriData(modeldir, fildid)

    adata, pd,rowcolorlist, colcolorlist,clusters_colors=genHeatmapdata(figdir,fildid, clusterlist,topN)

    genHeatmap(adata, figdir, fildid)
