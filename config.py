JSON_AS_ASCII=False

# SECRET_KEY = 'hard to guess string'
SQLALCHEMY_COMMIT_ON_TEARDOWN = True
# SQLALCHEMY_DATABASE_URI = 'mysql+pymysql://root:123456@192.168.110.177:80/phpMyAdmin4.8.5/baiyinyun?charset=utf8'
DIALECT = 'mysql'
DRIVER = 'pymysql'
USERNAME = 'root'
PASSWORD = 'rootx'
HOST = '192.168.110.241'
PORT = '3306'
DATABASE = 'baiyinyun'

SQLALCHEMY_DATABASE_URI = "{}+{}://{}:{}@{}:{}/{}?charset=utf8".format(DIALECT, DRIVER, USERNAME, PASSWORD, HOST, PORT,DATABASE)
SQLALCHEMY_TRACK_MODIFICATIONS = False
SQLALCHEMY_POOL_RECYCLE = 5

# 邮箱配置
MAIL_SERVER = 'smtp.163.com'  # 163smtp服务器
MAIL_PORT = 25  # 端口号

MAIL_USERNAME = "13635042214@163.com"
MAIL_PASSWORD = "abc123"   #填写客户端授权码

FLASKY_MAIL_SUBJECT_PREFIX = '[Flasky]'
FLASKY_MAIL_SENDER = '13635042214@163.com'
FLASKY_ADMIN = "13635042214@163.com"
MAIL_DEFAULT_SENDER = ("Arsz","13635042214@163.com")
MAIL_NAME= "Arsz"