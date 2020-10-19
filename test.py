import requests
import json,sys,io
import time
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding='utf-8')
# api路径
url="http://192.168.110.227:5555/api/fansu"
url1="http://192.168.110.227:5555/api/analprogress"
urld="http://192.168.110.227:5555/api/download"

url="http://192.168.110.227:5555/api/email"
# print(int(time.time()))
parms = {
    'data':{'user':{
        'name':'Arsz'
    },
    'gene_set': {
        'gid': 1,
        'name': '泛素化'
    },
    'cancer': {
        'TCGA-BRCA': {
            'cid': 1,
            'name': 'TCGA-BRCA'
        },
        'TCGA-LIHC': {
            'cid': 2,
            'name': 'TCGA-LIHC'
        }
    }}
}
dow = {
    'username':'Arsz',
    'filename':'20200925114903.zip'
}
xx={'name':'Arsz'}
headers = {
    'User-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/85.0.4183.102 Safari/537.36',
    'Spam': 'Eggs'
}

res = requests.post(url, data=json.dumps(parms),headers=headers)  # 发送请求

print(res.text,res.status_code)
x=0
while 1:
    rep = requests.post(url1, data=json.dumps(xx),headers=headers)  # 发送请求
    print(rep.text)
    x=json.loads(rep.text)['code']
    print(x)
    time.sleep(1)
    if x==4:
        break    
res1 = requests.post(urld, data=json.dumps(dow),headers=headers)  # 发送请求
print(res1.text)
# print(x)