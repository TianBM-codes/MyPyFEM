#!/usr/bin/env python3
"""
发送 POST 请求到 127.0.0.1:5000/api/startserver
携带参数:
  save_path="./test"
  structId=0
"""

import requests
import json

url = "http://127.0.0.1:5000/api/stopserver"
payload = {
    "save_path": "./test",
    "structId": 1
}

try:
    resp = requests.post(url, json=payload, timeout=5)
    resp.raise_for_status()
    print("返回状态码:", resp.status_code)
    print("返回内容:", resp.text)
except requests.RequestException as e:
    print("请求失败:", e)
