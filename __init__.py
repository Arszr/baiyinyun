"""
The flask application package.
"""

from flask import Flask
from flask_cors import CORS
from flask_restful import reqparse, Api, Resource
from flask_httpauth import HTTPTokenAuth
from flask_mail import Mail
from flask_sqlalchemy import SQLAlchemy

import web_app.config as config

app = Flask(__name__)
app.config.from_object(config)

db = SQLAlchemy(app, use_native_unicode='utf8')
# db.create_all()

mail = Mail()    #测试时可以直接在Mail()中写入app对象
api=Api()

db.init_app(app)
mail.init_app(app)    #这种方式是开发的时候常用的，因为我们要在其他模块中使用mail对象
api.init_app(app)
# set the secret key.  keep this really secret:
app.secret_key = 'A0Zr98j/3yX R~XHH!jmN]LWX/,?RT'
CORS(app)
