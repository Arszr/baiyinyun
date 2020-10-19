import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
# print(robjects.__file__)
import sys
sys.path.append('./')
import importlib
import json
import os
from web_app.function.WordCould import word_img
# importlib.reload(sys)
# #sys.setdefaultencoding('gbk')


class Ubiquitination():

    def __init__(self,disease,path):

        self.disease=disease
        path=path
        self.num=11
        # print('泛素化',self.disease,path)
        self.datas_path=self.data_path(self.disease)
        self.save=self.save_path(path,self.disease)

        # self.analysis(self.datas_path,self.disease,self.save)

    def load_R(self):
        pass

    def data_path(self,name):

        exp_path='./web_app/data/disease/exp_data/{}.txt'.format(name)
        clinical_path='./web_app/data/disease/clinical/{}.txt'.format(name)
        ubiquitina_path='./web_app/data/data/ubiq/UbiqGene.txt'
        # print(exp_path)
        return (exp_path,clinical_path,ubiquitina_path)

    def save_path(self,path,disease):
        path=path
        disease=disease

        sp=path+'/Ubiquitination/'

        if not os.path.exists(sp):
            os.makedirs(sp)
        
        sp=sp+disease+'/'
        
        if not os.path.exists(sp):
            os.makedirs(sp)
        
        # print(sp)
        return sp

    def analysis(self,data,name,save_path):
        
        data_path=data
        name=name
        save_path=save_path
        # print(data_path[0],'TCGA-BRCA',save_path)

        lime_all='./web_app/data/Difference/{}/limma_DEG_all.csv'.format(name)
        lime_n='./web_app/data/Difference/{}/limma_DEG_0.05.csv'.format(name)
        ubiq='./web_app/data/data/ubiq/UbiqGene.txt'
        pheno='./web_app/data/Difference/{}/pheno.csv'.format(name)
        exp_data='./web_app/data/disease/exp_data/{}.csv'.format(name)
        cli='./web_app/data/disease/clinical/{}.csv'.format(name)

        return (lime_all,lime_n,ubiq,pheno,exp_data,cli,save_path)

        fig1_result=self.fig1(lime_all,lime_n,ubiq,pheno,exp_data,cli,usa,save_path)

        
        


        # print(multiple[0])
        # print(single[0],single[1])

    def fig1(self,lime_all,lime_n,ubiq,pheno,exp_data,cli,save_path):

        lime_all=lime_all
        lime_n=lime_n
        ubiq=ubiq
        pheno=pheno
        exp_data=exp_data
        cli=cli
        save_path=save_path+'Fig1/'

        if not os.path.exists(save_path):
            os.makedirs(save_path)

        r=robjects.r
        # 加载差异分析文件
        r.source('./web_app/script/Conversion_Difference.r')
        r.source('./web_app/script/Single.r')
        r.source('./web_app/script/Survival_.r')
        r.source('./web_app/script/RelatedBubbles.r')
        # 调用差异分析函数完成差异分析
        # 构建差异基因，绘制差异火山图，热图
        difference=r.Difference(lime_all,lime_n,ubiq,pheno,exp_data,save_path)
        # print(difference[0],difference[1])

        # 单多因素分析
        single=r.SingleFactor(cli,exp_data,difference[0],save_path)
        # # print([i for i in single])
        

        survival=r.Survival_(single[0],single[1],difference[0],pheno,cli,save_path)
        # survival=r.Survival_(single[0],single[1],difference[0],pheno,save_path)
        # # # print([i for i in survival])
        # # # # 相关性气泡图
        
        bubble=r.RelatedBubbles(survival[0],cli,save_path)

        word_img(single[1],save_path)
        # print([i for i in bubble])

        result={
            'code':1,
            'difference':[i for i in difference],
            'single':[i for i in single],
            'survival':[i for i in survival],
            'bubble':[i for i in bubble],
        }
        return result

    def fig2(self,save_path):
        # save_path=save_path+'Fig2/'

        # if not os.path.exists(save_path):
        #     os.makedirs(save_path)

        r=robjects.r
        # 加载差异分析文件
        r.source('./web_app/script/GeneSurvivalModel/Heatmap.r')

        result={
            'code':2,
            
        }
        return result
        # 条带热图
        # r.source('./web_app/script/Heatmap.r')
        # r.Heatmap('./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/train.csv',
        # "./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Cox_genes_OS_pValue.csv",
        # './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/')

        # Lasso折线图
        # r.source('./web_app/script/LineLasso.r')
        # r.LineLasso(
        #     './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/train.csv',
        #     './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/CoxSingle_train.csv',
        #     './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/'
        # )
        # # 发散曲线
        # r.source('./web_app/script/CurveLasso.r')
        # r.CurveLasso('./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/train.csv',
        # './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/signature_gene.txt',
        # './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/')

        # # 随机生存森林
        # r.source('./web_app/script/RandomSurvivalForest.r')
        # r.Rsf('./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/train.csv',
        # './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/signature_gene.txt',
        # './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/')

    def epath(self):
        return self.save
    def progress(self):

        return 1
if __name__ == "__main__":
    a=Ubiquitination('TCGA-LIHC','./web_app/temp/Arsz')
    x=a.analysis(a.datas_path,a.disease,a.save)
    f1=a.fig1(x[0],x[1],x[2],x[3],x[4],x[5],x[6])
