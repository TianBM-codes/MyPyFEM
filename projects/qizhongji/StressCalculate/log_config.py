import logging

# 配置日志记录
def setup_logger():
    logger = logging.getLogger()  # 获取根日志记录器
    logger.handlers.clear()
    logger.setLevel(logging.DEBUG)  # 设置日志级别为DEBUG（记录所有级别的信息）

    # 创建文件处理器，将日志输出到文件
    file_handler = logging.FileHandler("rizhi.log", encoding="utf-8", mode="a")
    file_handler.setLevel(logging.INFO)
    file_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(file_formatter)

    # 控制台日志（DEBUG 及以上）
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    console_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console_handler.setFormatter(console_formatter)

    # 添加处理器
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    # 禁用第三方库的日志
    logging.getLogger("werkzeug").disabled = True  # 彻底禁用 Flask 请求日志
    logging.getLogger("pymysql").setLevel(logging.WARNING)  # 只记录警告及以上

# 调用日志配置
setup_logger()
