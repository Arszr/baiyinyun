

import sys,string
sys.path.append('./')
from web_app import app
from flask import render_template
from threading import Thread
from flask_mail import Message
from . import mail

def _send_async_mail(app, message):
    with app.app_context():# 将Flask的app推入栈中
        mail.send(message)
        
def send_mail(subject, to, body,filepath=None,filetype=None):
    message = Message(subject, recipients=[to])
    message.body = render_template('mail.txt')
    message.html = render_template('mail.html')
    filepath=filepath
    filetype=filetype
    filepath=path='E:/python_code/Baiyinyun'+filepath[1:]

    with app.open_resource(filepath) as fp:
        message.attach(filepath.split('/')[-1], filetype, fp.read())

    thr = Thread(target=_send_async_mail, args=[app, message])
    thr.start()
    # print(thr)
    return thr