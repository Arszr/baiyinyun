from sqlalchemy import Column, String, INTEGER,create_engine,DateTime,ForeignKey
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
import datetime
from . import db

"""
    说明:
        @ 这里的sessionmaker是一个工厂类,按照我们创建的engine绑定了数据库引擎
        @ 当我们使用这个Session的时候都会创建一个绑定引擎的Session,通过Session操作数据库
"""
# engine = create_engine(db)
Session = db.session

"""
    说明:
        @ 以下代码会创建flask数据库中user表的映射关系
        通过继承父类,指定字段类型表名将会建立映射关系
        
"""
Base = declarative_base()

# class Student(Base):
#     __tablename__ = 'student'
#     id = Column(Integer, primary_key=True)
#     name = Column(String(32))
#     school_id = Column(Integer, ForeignKey('school.id'))
#     stu2sch = relationship('School', backref='sch2stu')

# class School(Base):
#     __tablename__ = 'school'
#     id = Column(Integer, primary_key=True)
#     name = Column(String(32))
    
# class BaseModel(db.Model):
#     __abstract__ = True # 抽象类，可以将其他数据表中的公共字段存放在这个类中，然后继承该类
#     id =db.Column(db.Integer,primary_key=True,autoincrement=True)

class User(db.Model):
    __tablename__ = 'user'  # 指定表名字为 user
    # 主键id, 自增, Integer类型
    bid=Column(INTEGER,primary_key=True, autoincrement=True, nullable=False)
    # userName字段 varchar类型 限制45 不为空
    userName = Column(String(45), nullable=False,unique=True)
    # passWord字段, varchar类型 限制45 不为空
    passWord = Column(String(45), nullable=False)
    mobile = Column(String(45), nullable=False)
    email = Column(String(45), nullable=False)

    create_time = Column(DateTime, default=datetime.datetime.now)
    # 更新时间
    update_time = Column(DateTime, default=datetime.datetime.now, onupdate=datetime.datetime.now)

class DownloadFile(db.Model):

    __tablename__ = 'zipfile'  # 指定表名字为 user
    # 主键id, 自增, Integer类型
    zid=Column(INTEGER,primary_key=True, autoincrement=True, nullable=False)
    # userName字段 varchar类型 限制45 不为空
    userName = Column(String(45), nullable=False)
    # passWord字段, varchar类型 限制45 不为空
    filename=Column(String(100),nullable=False)
    filepath=Column(String(200),nullable=False)
    create_time = Column(DateTime, default=datetime.datetime.now)

db.create_all()