from util import scanpycluster
from plotutil import plotsingle,plotsinglegrey,plotsinglematrix,plotsingleglycan,plotsinglegrey_name
from deutil import getadata


fildidfeatlist={'liver':['0', '1', '6', '12'],'Flu_A19':['9','10'],'FibrosisB_A6':['0','2', '4','11'],
                'FibrosisB_B5':['0','3', '10'],'Flu_A20':['0','3'],'brainfrontal':['4','7','14'],
                'brainhippo':['6','11'],'normal':['4','6','7','11','14'],'normal_B4':['1','4','5','6'],
                'COVID_A42':['2','5','7','9'],'COVID_A43':['0','6','8','12'],'FibrosisC_A5':['5','6','10'],
                'FibrosisC_A4':['2','13']}

clustername={'FibrosisC_A5':{'6':'DAD','10':'Mucin','5':'Fibrosis'},
             'liver':{'6':'6_Hepatocytes','0':'0_Blood vessels L','12':'12_Blood vessels C','1':'1_unknown'}}

samplefiles={ 'liver':'BH3259 A5 normal liver PNG ISO tleh.txt'}

def clustersamples(fildid):
    modelname = fildid + '.h5ad'
    samples=[]
    samples.append(fildid)
    scanpycluster(datadir, fildid, samplefiles, batch, modeldir + modelname)


def desamples(modeldir,fildid):
    modelname=fildid+'.h5ad'
    samples=[fildid]
    adata=getadata(datadir, samples, samplefiles, modeldir+modelname)
    adata.write(modeldir+fildid+'DE.h5ad')


runloc='LOC'
if runloc=='LOC':
    dir='./'
datadir=dir+'data/'

modeldir=dir+'models_liver/'
figdir=dir+'figures_liver/'
samples=[]
batch=True
fildid='liver'


########cluster a single sample############
#clustersamples(fildid)


#generate adata for DE analysis
#desamples(modeldir,fildid)


#######umap plot and spatial plot##########
#plotsingle(modeldir,fildid,figdir)
#clusterlist=['0','1','6','12']
#plotsinglegrey(modeldir,fildid,figdir,clusterlist,'_ROIcluster')
#plotsinglegrey_name(modeldir,fildid,figdir,clusterlist,clustername[fildid])

#matrixlist=['16','14','13','11']
#plotsinglegrey(modeldir,fildid,figdir,matrixlist,'_matrix')
#plotsinglematrix(modeldir,fildid,figdir,matrixlist)


glycan=['689.1874','1079.5204','1663.5862']

#glycan=['1175.3761','1743.5917','1663.5862']

for feat in glycan:
    plotsingleglycan(modeldir,fildid,figdir,feat)



