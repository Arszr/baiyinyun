"""
Routes and views for the flask application.
"""

from flask_wtf.csrf import CsrfProtect
from flask import Flask, request, session, redirect, url_for, render_template, make_response, jsonify
from flask_mail import Message
from . import mail

import sys
sys.path.append('./')
from web_app import app
from threading import Thread
import web_app.BaiEmailOS  as wemail
import datetime
import json

# from web_app.db import User
from web_app.function.AnalysisType import AnalysisType
import web_app.db_open as db
import os
from flask.helpers import send_from_directory

dict_stats = {}

@app.route('/',methods=['GET'])
def home():
    return '哈哈哈'

@app.route('/api/login',methods=['GET','POST'])
def login():

    if request.method == 'POST':

        # print(request.form)
        user = request.data.decode('utf-8')
        if user=='' or user==None:
            username=request.form.get('username')
            password=request.form.get('password')
        # #获取到POST过来的数据，因为我这里传过来的数据需要转换一下编码。根据晶具体情况而定
        else:
            user_json = json.loads(user)
            username=user_json['username']
            password=user_json['password']

        user_=db.query_user(username)

        if username==user_.userName and password== user_.passWord:
            # print(user_json)
            #把区获取到的数据转为JSON格式。
            name=user_.userName
            # print({name})
            return jsonify({'meta':{'msg':'登录成功','code':200},'data':{'name':name}})

        elif username!=user_.userName or password!= user_.passWord:
            return jsonify({'meta':{'msg':'账号或密码错误','code':201}})

    return '123'

@app.route('/api/registered',methods=['GET','POST'])
def registered():

    user='Arsz'
    if request.method=='POST':
        pass
    else:
        user_flag=db.query_user(user)
        if user_flag==None:
            addmsg=db.add(user,'123','12345678910','948563218@qq.com')
            if addmsg==200:
                addmsg='注册成功'
                code=200
        else:
            addmsg='用户名已存在！'
            code=201
            
    return jsonify({'msg':addmsg,'code':code})

@app.route('/api/set_title',methods=['GET'])
def set_title():

    f = open("./web_app/json/MultiTypePrognosticModel.json", encoding='utf-8')
    setting = json.load(f)
    return jsonify(setting)

@app.route('/api/Sidebar',methods=['GET'])
def Sidebar():

    f = open("./web_app/json/Sidebar.json", encoding='utf-8')
    Sidebar = json.load(f)
    return jsonify(Sidebar)


@app.route('/api/fansu',methods=['POST'])
def fansu():

    method = request.json
    if method==None:
        method = request.data.decode('utf-8')
        method=json.loads(method)
        method_json=method['data']
    # print(method['data'])
    #获取到POST过来的数据，因为我这里传过来的数据需要转换一下编码。根据晶具体情况而定
    # json.loads(method)
    else:
        method_json = json.loads(method['data'])

    gene_set=method_json['gene_set']['gid']

    ans_cancer=method_json['cancer']
    cancer_set=[_['name'] for _ in ans_cancer.values()]

    user=method_json['user']['name']
    # thr = Thread(target=_send_async_mail, args=[app, message])
    
    dict_stats[user]=AnalysisType(user,gene_set,cancer_set)
    dict_stats[user].start()
    # print(gene_set,cancer_set)
    print({'code':456,'msg':'请求成功'})
    return jsonify({'code':456,'msg':'请求成功'})

@app.route('/api/analprogress',methods=['GET','POST'])
def analprogress():
    method = request.data.decode('utf-8')
    method_json=json.loads(method)
    userx=method_json['name']
    x=dict_stats[userx].get()
    # print(x)
    process=x['process']
    sumnum=x['sumnum']
    zipfile=x['zipfile']
    username=x['user']

    if process==0 or process<sumnum:
        flag=0
    elif process==sumnum:
        flag=1
    name=None
    if zipfile!=None:
        name=zipfile[1]
    if username!= None and zipfile!=None:
        user=db.query_user(username)
        filepath=zipfile[-1]
        email=user.email
        x=wemail.send_mail('你有一封新邮件,注意查收哦~',email,'\n\n所以爱会消失对吗？',filepath,filetype='zip/zip')
    
    print({'code':process,'total':sumnum,'msg':flag,'name':name})

    return jsonify({'code':process,'total':sumnum,'msg':flag,'name':name})

@app.route('/api/email',methods=['POST'])
def email():

    method = request.data.decode('utf-8')

    email='17313150835@qq.com'
    
    #     message = Message('你有一封新邮件~',
    #                       body='\n\n所以爱会消失对吗？',
    #                       recipients=[email])
    #     mail.send(message)
    x=wemail.send_mail('你有一封新邮件,注意查收哦~',email,'\n\n所以爱会消失对吗？')
    print(x)
    return jsonify('发送成功')

@app.route('/api/download',methods=['POST'])
def api_download():
    # os.path.realpath(__file__)打印出来的是当前路径/soft/flask/day1/wenwen.py
    # os.path.dirname(os.path.realpath(__file__))打印出来的是/soft/flask/day1/
    if request.method == 'POST':

        # print(request.form)
        user = request.data.decode('utf-8')
        # if user=='' or user==None:
        #     username=request.form.get('username')
        # # #获取到POST过来的数据，因为我这里传过来的数据需要转换一下编码。根据晶具体情况而定
        # else:
        user_json = json.loads(user)
        #     username=user_json['username']
        #     password=user_json['password']
        username=user_json['username']
        filename=user_json['filename']
        user_=db.query_zipfile(username,filename)

        # print(user_.userName,user_.filename,user_.filepath)

        return jsonify({"url":'/download/{}/{}'.format(user_.userName,user_.filename)})

@app.route('/api/remove',methods=['POST'])
def apiremove():


    return jsonify({'code':-1})

# @app.route('/api/download/<path:username>/<path:filename>',methods=['GET','POST'])
@app.route('/download/<path:username>/<path:filename>',methods=['GET','POST'])
def download(username,filename):
    # os.path.realpath(__file__)打印出来的是当前路径/soft/flask/day1/wenwen.py
    # os.path.dirname(os.path.realpath(__file__))打印出来的是/soft/flask/day1/
    username=username
    filename=filename
    user_=db.query_zipfile(username,filename)
    
    # current_dir = os.path.dirname(os.path.realpath(__file__))
    path='E:/python_code/Baiyinyun'+user_.filepath[1:]
    path=path.split('/')[:-1]
    path='/'.join(path)
    # print(current_dir)
    # return jsonify({'code':1})
    # print(username,filename)

    return send_from_directory(path, user_.filename,as_attachment=True)

