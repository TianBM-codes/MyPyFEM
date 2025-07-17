import numpy as np
from scipy.signal.windows import hamming
from scipy import signal
from scipy.integrate import cumulative_trapezoid
from scipy.signal import butter, filtfilt
import logging

# 获取与主模块相同的日志记录器
logger = logging.getLogger(__name__)

## 以下两个库用于测试
import pylab as pl
import pandas as pd


class VibAnalysis:
    data = ''
    fs = ''
    cutoff_fs = ''
    butter_a = ''
    butter_b = ''
    data_no_dc = ''
    freq = ''
    fft_y = ''

    def setFilter(self, fs, cutoff_fs=0.1):
        ## 设置滤波器参数，用于高通滤波
        ## fs：数据采样频率
        ## cutoff_fs： 滤波截止频率
        ## 说明，需要先设置滤波器参数，再设置待处理数据序列；
        try:
            self.fs = fs
            self.cutoff_fs = cutoff_fs
            cutoff_frequency = 1.5  # 截止频率，单位为 Hz
            # 设计高通滤波器
            self.butter_b, self.butter_a = butter(5, self.cutoff_fs / (self.fs / 2), btype='high')
        except Exception as e:
            print("in setFilter, Error: ", e)
            logger.error(f"Error in setFilter!  {e}")

    def setData(self, data, fs):
        ## 设置待分析数据序列
        ## data：待分析数据序列
        ## fs: 数据采样频率
        ## 说明，设置数据前，须先配置好滤波参数，否则不滤波
        try:
            self.fs = fs
            # 应用高通滤波器
            self.data_no_dc = filtfilt(self.butter_b, self.butter_a, data)
            return self.data_no_dc
        except Exception as e:
            print("in setData, Error: ", e)
            logger.error(f"Error in setData!  {e}")
            self.data_no_dc = data  ## 不滤波
            return data

    def FFT(self):
        ## 快速傅里叶变换分析，应用hamming窗
        ## 返回([频率]， [幅值])
        try:

            L = len(self.data_no_dc)  # 信号长度
            N = int(np.power(2, np.ceil(np.log2(L))))  # 下一个最近二次幂
            # N=L
            # 应用汉明窗
            window = hamming(L)
            signal_windowed = self.data_no_dc * window

            FFT_y1 = np.abs(np.fft.fft(signal_windowed, N)) / L * 2  # N点FFT 变化,但处于信号长度
            # 归一化处理
            window_sum = np.sum(window)

            FFT_y1 = FFT_y1  # / window_sum

            self.freq = np.arange(int(N / 2)) * self.fs / N  # 频率坐标
            self.fft_y = FFT_y1[range(int(N / 2))]  # 取一半

            return self.freq, self.fft_y
        except Exception as e:
            print("Error: ", e)
            logger.error(f"Error in FFT!  {e}")
            return [], []

    def findpeak(self, h, dis):
        ## 峰值查找函数
        ## h：幅值过滤阈值，即只考虑幅值大于h的峰值
        ## dis：峰值水平距离阈值，即极大值点的距离至少大于等于dis个水平单位
        ## return [峰值索引]
        try:
            num_peak_3 = signal.find_peaks(self.fft_y[:len(self.freq)], height=h, distance=dis)  # 低于指定height的信号都不考虑; distance表极大值点的距离至少大于等于10个水平单位
            # print(num_peak_3)
            id = num_peak_3[0]  ##峰值id
            if (len(id) > 1):
                # dpidx = peak[0]
                pkv = num_peak_3[1]['peak_heights']  ## 峰值
                fv = [id[i] * pkv[i] for i in range(len(pkv))]  ## 峰值与频率积，即面积
                dfv = [np.abs(np.round(100 * (fv[i + 1] - fv[i]) / fv[i + 1], 1)) for i in range(len(fv) - 1)]  ## 两个峰值之间的峰值频率面积差值
                dpkv = [np.abs(np.round(100 * (pkv[i + 1] - pkv[i]) / pkv[i + 1], 1)) for i in range(len(pkv) - 1)]  ## 两个峰值之间差值
                # print("fv: ", fv)
                # print("dfv: ", dfv)
                # print("dpkv: ", dpkv)
                idx = []
                val = 10  ## 判断阈值，为百分比，当两个峰值的差值百分比大于阈值时，人为是两个峰值，否则是一个
                for i, (k, v) in enumerate(zip(dfv, dpkv)):
                    if (k > val or v > val):
                        idx.append(i)
                pkid = [0]
                for i in range(len(idx) - 1):
                    pkid.append(np.where(pkv == np.max(pkv[idx[i] + 1:idx[i + 1] + 1]))[0][0])
                pkid.append(np.where(pkv == np.max(pkv[idx[-1] + 1:]))[0][0])

                return num_peak_3[0][pkid]
            elif (len(id) == 1):
                return num_peak_3[0]
            else:
                return []
        except Exception as e:
            print("Error: ", e)
            logger.error(f"Error in findpeak!  {e}")
            return []

    def integrate(self):
        ## 数值积分函数，返回待分析数据的一次积分和二次积分结果
        ## return [二次积分结果]，[一次积分结果]
        ## 说明，如果配置了滤波参数，则每次积分前和积分后均应用一次滤波操作；
        try:
            # 定义时间间隔
            dt = 1.0 / self.fs  # 假设每个数据点之间的时间间隔为 1 秒

            # 第一次积分：从加速度得到速度
            v = cumulative_trapezoid(self.data_no_dc, dx=dt, initial=0)

            try:
                vv = filtfilt(self.butter_b, self.butter_a, v)  ## 先滤波处理，去除直流分量，因为积分速度只考虑动态速度
            except Exception as ef:
                print("in integrate first filter info: ", ef)
                logger.info(f"in integrate first filter!  {ef}")
                vv = v  ## 不进行滤波

            # 第二次积分：从速度得到位移
            s = cumulative_trapezoid(vv, dx=dt, initial=0)
            try:
                ss = filtfilt(self.butter_b, self.butter_a, s)  ## 滤波处理，去除直流分量，因为积分位移只考虑动态位移
            except Exception as ef:
                print("in integrate second filter info: ", ef)
                logger.info(f"in integrate second filter!  {ef}")
                ss = s  ## 不进行滤波

            return ss, vv
        except Exception as e:
            print("Error: ", e)
            logger.error(f"Error in integrate!  {e}")
            return [], []


if __name__ == "__main__":
    ## 导入测试数据
    b = np.loadtxt('./data/jingli3-2.csv', dtype=str, delimiter=',', skiprows=0, encoding='utf-8')
    bb = pd.DataFrame(b)
    bb = bb.apply(pd.to_numeric, errors='coerce').fillna(0.0)
    data = bb[3]

    ## 选择测试数据区间
    startPoint = 0  # 44000 #4-43000 #3-36500# 2-44000
    ts = 110
    sampling_rate = 500  # 31.75
    fs = 500
    t = np.arange(0, ts, 1.0 / sampling_rate)
    fft_size = len(t)
    t0 = 47
    t1 = 55
    data = data[startPoint:startPoint + fft_size]
    data = data[sampling_rate * t0:sampling_rate * t1]
    t = t[sampling_rate * t0:sampling_rate * t1]

    vib = VibAnalysis()  ## 创建动态数据分析对象
    vib.setFilter(fs, 1.5)  ## 配置滤波器参数
    data_noDC = vib.setData(data, fs)  ## 设置待分析数据
    frq, fy = vib.FFT()  ## 快速傅里叶分析
    print(f"fy mean: {np.mean(fy)}, std: {np.std(fy)}, m+3std:{np.mean(fy) + 6 * np.std(fy)}")

    peak = vib.findpeak(np.mean(fy) + 6 * np.std(fy), 40)  ## 幅频曲线峰值提取
    s, v = vib.integrate()  ## 积分计算，得到位移和速度

    print("peak: ", peak)
    print("pkf: ", frq[peak])
    print("pkv: ", fy[peak])
    pl.subplot(4, 1, 1)
    pl.plot(t, data)
    pl.subplot(4, 1, 2)
    pl.plot(frq, fy)
    pl.plot(frq[peak], fy[peak], 'o')
    pl.subplot(4, 1, 3)
    pl.plot(t, v)
    pl.subplot(4, 1, 4)
    pl.plot(t, s)
    pl.show()
