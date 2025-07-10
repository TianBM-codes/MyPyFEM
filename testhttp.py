import requests

url = "http://localhost:5000/api/re_calculate_stiff"
params = {"e_value": 2.4e11}  # 参数放在 URL
data = {}  # 请求体可以为空

response = requests.post(url, params=params, json=data)
print("Status Code:", response.status_code)
print("Response:", response.json())