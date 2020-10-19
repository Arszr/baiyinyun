
from web_app.DBBase import User,Session,DownloadFile
# import web_app.DBBase
session = Session()
# 操作数据库
# 增
def add(username,password,mobile,email):
    user = User(
        userName=username,
        passWord=password,
        mobile=mobile,
        email=email
    )
    session.add(user)
    session.commit()
    session.close()
    return 200

# 查 
def query_user(username):
    # pass
    user=session.query(User).filter(User.userName == username).first()
    session.close()
    # print(user)
    return user
    # user__all = session.query(User).all()
    # for user in user__all:
    #     print("用户名: {}      密码: {}".format(user.userName, user.passWord))

    # count = session.query(User).count()
    # print("总记录条数为: {}".format(count))

    # first_name_eq_jack = session.query(User).filter_by(userName='jack6').first()
    # print("查询用户名为 jack的用户的id为: {}".format(first_name_eq_jack.id))

    # user__slice = session.query(User).slice(0, 3)
    # for slice in user__slice:
    #     print("分页查询三条的结果为: id为; {}     用户名为: {}".format(slice.id, slice.userName))
# 删
def delete(id):
    user = session.query(User).filter(User.id == id).first()

    session.delete(user)
    session.commit()
    session.close()

# 改
def update():
    user = session.query(User).filter(User.id == 8).first()
    user.userName = 'update'
    user.passWord = '123321'

    session.add(user)
    session.commit()
    session.close()


def add_zipfile(username,filename,filepath):
    zipfile = DownloadFile(
        userName=username,
        filename=filename,
        filepath=filepath,
    )
    # session.ping(reconnect=True)
    session.add(zipfile)
    session.commit()
    session.close()
    return 200

def query_zipfile(username,zipfile):
    filedata=session.query(DownloadFile).filter(DownloadFile.userName == username ,DownloadFile.filename==zipfile).first()
    session.close()
    # print(filedata)
    return filedata