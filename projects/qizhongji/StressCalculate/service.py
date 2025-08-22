from data_insert import *
from flask import Flask, request, jsonify


app = Flask(__name__)
@app.route('/start_cal', methods=['POST'])
def start_cal():
    """
    POST 请求调用 main() 方法
    """
    main()
    return jsonify({"code": 200, "msg": "计算完成"})


if __name__ == '__main__':
    # 启动 Flask 服务
    app.run(host='0.0.0.0', debug=False, port=52000)