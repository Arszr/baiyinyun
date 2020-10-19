# coding: utf-8
import os
import sys
sys.path.append("..")
from wordcloud import WordCloud,STOPWORDS
import itchat
import re
import pandas as pd
from wordcloud import WordCloud, ImageColorGenerator
import numpy as np
import PIL.Image as Image
import matplotlib.colors as colors
import matplotlib
matplotlib.rcParams['backend'] = 'SVG'
import matplotlib.pyplot as plt
from multiprocessing import Pool
import multiprocessing

def word_img(data,save,bg_path="./web_app/data/image/wordcloud_bg/wordcloud_bg3.jpg"):

    path=data
    text = pd.read_csv(path, sep=',')
    name=list(text.factors)
    value=text.iloc[:,5]*10
    dic = dict(zip(name, value))
    #colormaps=colors.ListedColormap(['#00cc03','#82d0fe','#3949fb','#36c80e'])
    colormaps=colors.ListedColormap(["#84aaff","#E78AC3","#52deb2","#d9754e"])
    # 对文本进行分词
    # cut_text =''.join(jieba.cut(text))
    # 读取图片
    color_mask = np.array(Image.open(bg_path))
    # 生成词云
    # 导入字体，不导入也并不影响，若是中文的或者有其他的字符需要自己选择合适的字体包
    cloud = WordCloud(
                    background_color="white",
                    #   width =1000,
                    #   height=1000,
                    collocations=True,
                    relative_scaling=.5,
                    scale=2,
                    mask=color_mask,
                    prefer_horizontal  =0.9,
                    max_words=500,
                    max_font_size=100,
                    # colormap=colormaps
                    )
    word_cloud = cloud.fit_words(dic)
    # 输出图片
    plt.axis( 'off' )
    # plt.figsize=(10.0, 10.0)
    plt.imshow(word_cloud)

    plt.savefig(save+"Single_element.svg",format = 'svg',dpi=500)

if __name__ == "__main__":
    word_img('./web_app/temp/admin/Ubiquitination/TCGA-BRCA/Fig1/CoxSingle_train.csv','./web_app/temp/admin/Ubiquitination/TCGA-BRCA/Fig1/')

    
