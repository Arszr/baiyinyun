#start.py
# from flask import Flask

# app = Flask(__name__) 

# @app.route('/')
# def hello():
#     return 'hello docker&flask'

# if __name__ == '__main__':
#     app.run()

from flask import Flask, g
from flask_restful import reqparse, Api, Resource
from flask_httpauth import HTTPTokenAuth

# Flask相关变量声明
app = Flask(__name__)
api = Api(app)




if __name__ == "__main__":
    app.run(debug=True)