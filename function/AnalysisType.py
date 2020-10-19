import os,re
# os.environ["R_HOME"] = "C:\\Program Files\\R\\R-4.0.2"
import string
import pandas as pd
import numpy as np

import rpy2.robjects as robjects
# print(robjects.__file__)
import sys
import importlib
import json
import zipfile
import datetime
import web_app.db_open as db
import web_app.BaiEmailOS  as wemail
import threading
importlib.reload(sys)
#sys.setdefaultencoding('gbk')


class AnalysisType(threading.Thread):
    

    def __init__(self,user,gene_set,cancer_set,target=None,args=None,name=None):
        threading.Thread.__init__(self)
        self._target = target
        self._args = args
        self.setName(str(user))

        self.user=str(user)
        self.gene_set=str(gene_set)
        self.cancer_set=cancer_set
        
        self.process=0
        self.sumnum=2
        self.zipfile=None
        # self.process={}
        self.global_path=self.floter(self.user)
        # self.methods(self.gene_set,self.cancer_set)
        self.mes=self.methods



    def run(self):
        try:
            
            self.mes(self.gene_set,self.cancer_set)
        finally:
            # del self._target, self._args
            # del self.mes
            pass
            

    def floter(self,path):
        temp_path='./web_app/temp/'+path
        if not os.path.exists(temp_path):
            os.makedirs(temp_path)

        return temp_path


    def methods(self,gene_set,cancer_set):

        f = open("./web_app/json/gentSet_id.json", encoding='utf-8')  
        setting = json.load(f)
        methods_name=setting[gene_set]
        cancer_list=cancer_set

        # 导入模块
        moudel='web_app.function.'+methods_name
        
        metaclass=importlib.import_module(moudel)
        # 获取模块中的类
        met_class=getattr(metaclass,methods_name)
        self.sumnum=self.sumnum*len(cancer_list)
        for index,_ in enumerate(cancer_list):
            a=met_class(_,self.global_path)
            x=a.analysis(a.datas_path,a.disease,a.save)
            f1=a.fig1(x[0],x[1],x[2],x[3],x[4],x[5],x[6])
            # for i in range(100000000):
            #     i+=1
            # f1={'code':1}
            # print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
            self.process=f1['code']+(index+1)
            for i in range(10000000):
                i+=1
            f2=a.fig2(1)
            self.process=f2['code']+(index+1)
            
            
            if self.process==self.sumnum:
                # print(a.epath())
                path=a.epath().split('/')[:-2]
                path='/'.join(path)

                plint=zipPath=a.epath().split('/')
                zipPath[2]='download'
                zipPath='/'.join(zipPath)
                # print(zipPath)
                j='.'
                for i in range(1,len(plint)):
                    j=j+'/'+plint[i]
                    if not os.path.exists(j):
                        os.makedirs(j)
                # print(j)
                dtime=datetime.datetime.now()
                x=datetime.datetime.strftime(dtime,'%Y%m%d%H%M%S')
                # zipFilePath=path+"/%s.zip"%(x)
                zipfile=self.zip_file_path(path, j, x+'.zip')
                addmsg=db.add_zipfile(self.user,zipfile[1],zipfile[-1])
                if addmsg==200:
                    # print('添加成功')
                    self.zipfile=zipfile

                else:
                    # print('添加失败')
                    pass

                

        # print(type(met_class))
    
    '''
        def resultzip(self,path):
            path=path.split('/')[:-3]
            path='/'.join(path)
            dtime=datetime.datetime.now()
            x=datetime.datetime.strftime(dtime,'%Y%m%d%H%M%S')
            zipFilePath=path+"/%s.zip"%(x)

            zipFile=zipfile.ZipFile(zipFilePath,"w",zipfile.ZIP_DEFLATED) 

            absDir=path
            print(absDir,zipFilePath)
            
            self.writeAllFileToZip(absDir,zipFile) #开始压缩。如果当前工作目录跟脚本所在目录一样，直接运行这个函数。
            # #执行这条压缩命令前，要保证当前工作目录是脚本所在目录(absDir的父级目录)。否则会报找不到文件的错误。
            print("压缩成功")
            
            return None 
            
        def writeAllFileToZip(self,absDir,zipFile):
            # z = zipFile
            # startdir = absDir
            # for dirpath, dirnames, filenames in os.walk(startdir):
            #     for filename in filenames:
            #         z.write(os.path.join(dirpath, filename))
            # z.close()
            # return
            for f in os.listdir(absDir):
                absFile=os.path.join(absDir,f) #子文件的绝对路径
                print('-------',absFile)
                if os.path.isdir(absFile): #判断是文件夹，继续深度读取。
                    relFile=absFile #改成相对路径，否则解压zip是/User/xxx开头的文件。
                    zipFile.write(relFile) #在zip文件中创建文件夹
                    self.writeAllFileToZip(absFile,zipFile) #递归操作
                else: #判断是普通文件，直接写到zip文件中。
                    relFile=absFile #改成相对路径
                    zipFile.write(relFile)
    '''      
    
    def get_zip_file(self,input_path, result):
        """
        对目录进行深度优先遍历
        :param input_path:
        :param result:
        :return:
        """
        files = os.listdir(input_path)
        for file in files:
            if os.path.isdir(input_path + '/' + file):
                self.get_zip_file(input_path + '/' + file, result)
            else:
                result.append(input_path + '/' + file)


    def zip_file_path(self,input_path, output_path, output_name):
        """
        压缩文件
        :param input_path: 压缩的文件夹路径
        :param output_path: 解压（输出）的路径
        :param output_name: 压缩包名称
        :return:
        """
        f = zipfile.ZipFile(output_path + output_name, 'w', zipfile.ZIP_DEFLATED)
        filelists = []
        self.get_zip_file(input_path, filelists)
        for file in filelists:
            
            if file.endswith(".pdf") or file.endswith(".png") or file.endswith(".svg"):
                name=file.split('/')[-1]
                # f.write(filename=file,arcname=name)
                f.write(filename=file)
        # 调用了close方法才会保证完成压缩
        f.close()
        return (output_path,output_name,output_path+output_name)

    def pull_email(self,username,filepath):
        user=db.query_user(username)
        filepath=filepath
        email=user.email
        x=wemail.send_mail('你有一封新邮件,注意查收哦~',email,'\n\n所以爱会消失对吗？',filepath,filetype='zip')

        return 200

    def get(self):

        return {'process':self.process,'sumnum':self.sumnum,'zipfile':self.zipfile,'user':self.user}
    


if __name__ == "__main__":
    AnalysisType('Arsz','1',['TCGA-BRCA'])