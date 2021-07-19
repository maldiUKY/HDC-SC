import numpy as np
import scanpy as sc
import pandas as pd
from util import concatenateadata,genadata,readOrifile
from collections import Counter


diseasedict={'covid':['COVID_A40','COVID_A41','COVID_A42','COVID_A43','COVID_A44'],
             'fibB':['FibrosisB_A6','FibrosisB_A9','FibrosisB_A13','FibrosisB_B5','FibrosisB_B9'],
             'fibC':['FibrosisC_A4','FibrosisC_A5','FibrosisC_B2','FibrosisC_B5'],
             'flu':['Flu_A16','Flu_A18','Flu_A19','Flu_A20','Flu_A21'],
             'normal':['normal_A3','normal_B4','normal_A6']}


def getadata(dir,samples,samplefiles,modelname):
    num=0
    adatas=[]
    for sample in samples:
        infile=dir+samplefiles[sample]
        print(infile)
        xs,ys,nodefeatlabels,nodefeatlabelsdict,nodefeats=readOrifile(infile,num)
        varindex = pd.DataFrame(index=nodefeatlabels)
        adata=genadata(xs,ys,nodefeats,sample,varindex)
        sc.pp.normalize_total(adata, target_sum=1, inplace=True)
        sc.pp.log1p(adata, base=2)
        adatas.append(adata)
    adata=concatenateadata(adatas)
    adatamodel = sc.read(modelname)
    adata.obs['leiden']=adatamodel.obs['leiden']
    return adata


def countlabels(labels,labellist):
    labelcount=[]
    for label in labels:
        labelcount.append(labellist.count(label))
    return labelcount


def outputlabelcount(adata,samples,modeldir,fildid):
    labels = adata.obs.leiden.unique().tolist()
    labelcounts=[]
    for sample in samples:
        ad = adata[adata.obs['library_id']==sample, :]
        labellist = ad.obs['leiden'].tolist()
        labelcount=countlabels(labels,labellist)
        labelcounts.append(labelcount)
    df = pd.DataFrame(data=labelcounts,index=samples,columns=labels)
    df.to_csv(modeldir+fildid+'.csv')


def getLabellist(adata,samples):
    labellist=[]
    ad = adata[adata.obs['library_id'].isin(samples), :]
    labels=ad.obs['leiden'].tolist()
    for label, count in Counter(labels).items():
        if count > 2000:
            labellist.append(label)
    print(Counter(labels))
    return labellist



seldict= {'covidnormal':'normal','covidflu':'flu',
          'covidfibB':'fibB','covidfibC':'fibC'}


def getRepresentpix(adata):
    sels=['covid','flu','fibB','fibC','normal']
    adatas=[]
    for sel in sels:
        label=getLabellist(adata, diseasedict[sel])
        ad = adata[adata.obs['library_id'].isin(diseasedict[sel]), :]
        ad = ad[ad.obs['leiden'].isin(label), :]
        ad.uns['group'] = {sel: {}}
        adatas.append(ad)
    adata = adatas[0].concatenate(
        adatas[1],adatas[2],adatas[3],adatas[4],
        batch_key="group_id",
        uns_merge="unique",
        batch_categories=[
            k
            for d in [
                adatas[0].uns["group"],
                adatas[1].uns["group"],
                adatas[2].uns["group"],
                adatas[3].uns["group"],
                adatas[4].uns["group"]
            ]
            for k, v in d.items()
        ],
    )
    return adata


def findde(adata,modeldir,fildid):
    sel=seldict[fildid]

    covidlabel = getLabellist(adata, diseasedict['covid'])
    adcovid = adata[adata.obs['library_id'].isin(diseasedict['covid']), :]
    adcovid = adcovid[adcovid.obs['leiden'].isin(covidlabel), :]
    adcovid.uns['group'] = {'covid': {}}

    otherlabel = getLabellist(adata, diseasedict[sel])
    adother=adata[adata.obs['library_id'].isin(diseasedict[sel]), :]
    adother = adother[adother.obs['leiden'].isin(otherlabel), :]
    adother.uns['group'] = {sel: {}}

    print('covid:'+','.join(covidlabel))
    print(sel + ':'+','.join(otherlabel))

    adatas = [adcovid, adother]

    adatas = adatas[0].concatenate(
        adatas[1],
        batch_key="group_id",
        uns_merge="unique",
        batch_categories=[
            k
            for d in [
                adatas[0].uns["group"],
                adatas[1].uns["group"]
            ]
            for k, v in d.items()
        ],
    )

    METHOD = "wilcoxon"
    TOPN = 7
    sc.tl.rank_genes_groups(adatas, method=METHOD, groupby='group_id', key_added='group_covidother_results')
    sc.pl.rank_genes_groups_stacked_violin(adatas, key='group_covidother_results', n_genes=TOPN, swap_axes=True)

    result = adatas.uns['group_covidother_results']
    groups=['covid']
    sta = pd.DataFrame(
        {group + '_' + key: result[key][group]
         for group in groups for key in ['names', 'scores','pvals','pvals_adj', 'logfoldchanges']})
    sta.to_csv('./stanew/'+fildid+'.csv', sep=',')





def findmatrixDE(adata,labellist,fildid):
    METHOD = "wilcoxon"
    TOPN = 5
    sc.tl.rank_genes_groups(adata, method=METHOD, groups=labellist,groupby='leiden', key_added='matrix')
    sc.pl.rank_genes_groups_stacked_violin(adata, key='matrix', n_genes=TOPN, swap_axes=True,show=False,
                                           save=fildid+'_top5.png',title=fildid+' top 5 markers for matrix')

    result = adata.uns['matrix']
    groups=labellist
    sta = pd.DataFrame(
        {group + '_' + key: result[key][group]
         for group in groups for key in ['names', 'scores','pvals','pvals_adj', 'logfoldchanges']})
    sta.to_csv('./stamatrix/'+fildid+'.csv', sep=',')



def findnormalDE(adata,modeldir,fildid):
    print(fildid)
    sel = fildid.lstrip('normal')
    covidlabel = getLabellist(adata, diseasedict['normal'])
    adcovid = adata[adata.obs['library_id'].isin(diseasedict['normal']), :]
    adcovid = adcovid[adcovid.obs['leiden'].isin(covidlabel), :]
    adcovid.uns['group'] = {'normal': {}}

    otherlabel = getLabellist(adata, diseasedict[sel])
    adother=adata[adata.obs['library_id'].isin(diseasedict[sel]), :]
    adother = adother[adother.obs['leiden'].isin(otherlabel), :]
    adother.uns['group'] = {sel: {}}

    print('normal:'+','.join(covidlabel))
    print(sel + ':'+','.join(otherlabel))

    adatas = [adcovid, adother]

    adatas = adatas[0].concatenate(
        adatas[1],
        batch_key="group_id",
        uns_merge="unique",
        batch_categories=[
            k
            for d in [
                adatas[0].uns["group"],
                adatas[1].uns["group"]
            ]
            for k, v in d.items()
        ],
    )
    METHOD = "wilcoxon"
    TOPN = 7
    sc.tl.rank_genes_groups(adatas, method=METHOD, groupby='group_id', key_added='group_covidother_results')
    #sc.pl.rank_genes_groups_stacked_violin(adatas, key='group_covidother_results', n_genes=TOPN, swap_axes=True)

    result = adatas.uns['group_covidother_results']
    groups=[sel]
    sta = pd.DataFrame(
        {group + '_' + key: result[key][group]
         for group in groups for key in ['names', 'scores','pvals','pvals_adj', 'logfoldchanges']})
    sta.to_csv('./stanew/'+fildid+'.csv', sep=',')


def finddecovidDE(adata):
    METHOD = "wilcoxon"
    TOPN = 7
    sc.tl.rank_genes_groups(adata, groupby='leiden', method=METHOD, key_added='15toother', groups=['15'])
    result = adata.uns['15toother']
    groups = ['15']
    sta = pd.DataFrame(
        {group + '_' + key: result[key][group]
         for group in groups for key in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']})
    sta.to_csv('./sta/' + 'covid15tootherall' + '.csv', sep=',')

def finddeall(adata,modeldir,fildid):
    for sel in ['normal','flu','fibB','fibC']:
        covidlabel = getLabellist(adata, diseasedict['covid'])
        adcovid = adata[adata.obs['library_id'].isin(diseasedict['covid']), :]
        adcovid = adcovid[adcovid.obs['leiden'].isin(covidlabel), :]
        adcovid.uns['group'] = {'covid': {}}

        otherlabel = getLabellist(adata, diseasedict[sel])
        adother=adata[adata.obs['library_id'].isin(diseasedict[sel]), :]
        adother = adother[adother.obs['leiden'].isin(otherlabel), :]
        adother.uns['group'] = {sel: {}}

        print('covid:'+','.join(covidlabel))
        print(sel + ':'+','.join(otherlabel))

        adatas = [adcovid, adother]

        adatas = adatas[0].concatenate(
            adatas[1],
            batch_key="group_id",
            uns_merge="unique",
            batch_categories=[
                k
                for d in [
                    adatas[0].uns["group"],
                    adatas[1].uns["group"]
                ]
                for k, v in d.items()
            ],
        )

        METHOD = "wilcoxon"
        TOPN = 7
        sc.tl.rank_genes_groups(adatas, method=METHOD, groupby='group_id', key_added='group_covidother_results')
        sc.pl.rank_genes_groups_stacked_violin(adatas, key='group_covidother_results', n_genes=TOPN, swap_axes=True)

        result = adatas.uns['group_covidother_results']
        groups=['covid']
        sta = pd.DataFrame(
            {group + '_' + key: result[key][group]
             for group in groups for key in ['names', 'scores','pvals','pvals_adj', 'logfoldchanges']})
        sta.to_csv('./sta/'+fildid+'_'+sel+'.csv', sep=',')



def finddeallnormal(adata,modeldir,fildid):
    for sel in ['covid','flu','fibB','fibC']:
        covidlabel = getLabellist(adata, diseasedict['normal'])
        adcovid = adata[adata.obs['library_id'].isin(diseasedict['normal']), :]
        adcovid = adcovid[adcovid.obs['leiden'].isin(covidlabel), :]
        adcovid.uns['group'] = {'normal': {}}

        otherlabel = getLabellist(adata, diseasedict[sel])
        adother=adata[adata.obs['library_id'].isin(diseasedict[sel]), :]
        adother = adother[adother.obs['leiden'].isin(otherlabel), :]
        adother.uns['group'] = {sel: {}}

        print('normal:'+','.join(covidlabel))
        print(sel + ':'+','.join(otherlabel))

        adatas = [adcovid, adother]

        adatas = adatas[0].concatenate(
            adatas[1],
            batch_key="group_id",
            uns_merge="unique",
            batch_categories=[
                k
                for d in [
                    adatas[0].uns["group"],
                    adatas[1].uns["group"]
                ]
                for k, v in d.items()
            ],
        )

        METHOD = "wilcoxon"
        TOPN = 7
        sc.tl.rank_genes_groups(adatas, method=METHOD, groupby='group_id', key_added='group_covidother_results')
        #sc.pl.rank_genes_groups_stacked_violin(adatas, key='group_covidother_results', n_genes=TOPN, swap_axes=True)

        result = adatas.uns['group_covidother_results']
        groups=[sel]
        sta = pd.DataFrame(
            {group + '_' + key: result[key][group]
             for group in groups for key in ['names', 'scores','pvals','pvals_adj', 'logfoldchanges']})
        sta.to_csv('./stanew/'+fildid+'_'+sel+'.csv', sep=',')


def finddiseasenormalde(adata,fildid):
    for sel in ['covid','flu','fibB','fibC']:
        covidlabel = getLabellist(adata, diseasedict['normal'])
        adcovid = adata[adata.obs['library_id'].isin(diseasedict['normal']), :]
        adcovid = adcovid[adcovid.obs['leiden'].isin(covidlabel), :]
        adcovid.uns['group'] = {'normal': {}}

        otherlabel = getLabellist(adata, diseasedict[sel])
        adother=adata[adata.obs['library_id'].isin(diseasedict[sel]), :]
        adother = adother[adother.obs['leiden'].isin(otherlabel), :]
        adother.uns['group'] = {sel: {}}

        print('normal:'+','.join(covidlabel))
        print(sel + ':'+','.join(otherlabel))

        adatas = [adcovid, adother]

        adatas = adatas[0].concatenate(
            adatas[1],
            batch_key="group_id",
            uns_merge="unique",
            batch_categories=[
                k
                for d in [
                    adatas[0].uns["group"],
                    adatas[1].uns["group"]
                ]
                for k, v in d.items()
            ],
        )

        METHOD = "wilcoxon"
        TOPN = 7
        sc.tl.rank_genes_groups(adatas, method=METHOD, groupby='group_id', key_added='group_covidother_results')
        sc.pl.rank_genes_groups_stacked_violin(adatas, key='group_covidother_results', n_genes=TOPN, swap_axes=True)

        result = adatas.uns['group_covidother_results']
        groups=['normal']
        sta = pd.DataFrame(
            {group + '_' + key: result[key][group]
             for group in groups for key in ['names', 'scores','pvals','pvals_adj', 'logfoldchanges']})
        sta.to_csv('./sta/'+fildid+'_normal_'+sel+'.csv', sep=',')


def findDE(adata,labellist,ref,fildid):
    METHOD = "wilcoxon"
    TOPN = 5
    sc.tl.rank_genes_groups(adata, method=METHOD, groups=labellist,reference=ref,groupby='leiden', key_added='matrix')
    sc.pl.rank_genes_groups_stacked_violin(adata, key='matrix',n_genes=TOPN, swap_axes=True,show=False,
                                           save=fildid+'_top5.png',title=fildid+' top 5 markers for matrix')

    result = adata.uns['matrix']
    groups=labellist
    sta = pd.DataFrame(
        {group + '_' + key: result[key][group]
         for group in groups for key in ['names', 'scores','pvals','pvals_adj', 'logfoldchanges']})
    sta.to_csv('./sta_covid/'+fildid+'.csv', sep=',')