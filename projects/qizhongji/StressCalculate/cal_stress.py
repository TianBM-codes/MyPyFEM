from math import sin, cos, radians, sqrt, degrees, atan2, tan
import numpy as np
from scipy.optimize import fsolve
from QiZhongJiZiTai import MQ1330
import math
import pymysql
from datetime import datetime
import pytz
import logging
import log_config  # 导入日志配置模块
# 臂架、象鼻梁、大拉杆的旋转角度为常规坐标系，即x轴正方向为0度
# 综合应力计算过程中，定义y轴-方向为0°，逆时针方向为-方向，顺时针方向为+方向
logger = logging.getLogger()

def get_part_angles(arm_angle_deg):
    """
    求起重机各部件的角度：象鼻梁角、大拉杆角、小拉杆角、齿条架角

    参数:
    arm_angle_deg -- 臂架角度（单位：度）

    返回:
    (big_rod_angle_deg, elephant_trunk_angle_deg, small_rod_angle_deg, rack_angle_deg)
    """
    mq = MQ1330()
    angles = mq.getPartOrientation(arm_angle_deg)

    # angles返回顺序：[臂架角, 象鼻梁角, 大拉杆角, 小拉杆角, 平衡梁角, 齿条架角]
    elephant_trunk_angle = angles[1]  # 象鼻梁角，弧度
    big_rod_angle = angles[2]  # 大拉杆角，弧度
    small_rod_angle = angles[3]  # 小拉杆角，弧度
    rack_angle = angles[5]  # 齿条架角，弧度

    # 返回角度制（转换为度数）
    elephant_trunk_angle_deg = -elephant_trunk_angle * 180 / np.pi
    big_rod_angle_deg = big_rod_angle * 180 / np.pi
    small_rod_angle_deg = -small_rod_angle * 180 / np.pi
    rack_angle_deg = -rack_angle * 180 / np.pi

    return big_rod_angle_deg, elephant_trunk_angle_deg, small_rod_angle_deg, rack_angle_deg

# 根据分力的方向与大小求合力的大小与方向，y轴-方向为0°，逆时针方向为-方向，顺时针方向为正方向
def calculate_resultant_force(F1_magnitude, F1_angle, F2_magnitude, F2_angle):
    # 将角度转换为弧度
    F1_angle_rad = radians(F1_angle)
    F2_angle_rad = radians(F2_angle)

    # 计算每个力的水平分量和垂直分量
    F1_x = F1_magnitude * cos(F1_angle_rad)
    F1_y = F1_magnitude * sin(F1_angle_rad)

    F2_x = F2_magnitude * cos(F2_angle_rad)
    F2_y = F2_magnitude * sin(F2_angle_rad)

    # 合力的水平分量和垂直分量
    Rx = F1_x + F2_x
    Ry = F1_y + F2_y

    # 计算合力的大小
    resultant_magnitude = sqrt(Rx ** 2 + Ry ** 2)

    # 计算合力的方向（角度），使用atan2函数来处理四个象限
    resultant_angle = degrees(atan2(Ry, Rx))

    return resultant_magnitude, resultant_angle

# 已知合力求分力，R表示合力大小，alpha表示合力方向，theta1与theta2分别表示分力的方向
def calculate_forces(R, alpha, theta1, theta2):
 # 定义一个方程系统来表示两个分力的大小
    def equations(vars):
        F1, F2 = vars
        eq1 = F1 * np.cos(np.radians(theta1)) + F2 * np.cos(np.radians(theta2)) - R * np.cos(np.radians(alpha))
        eq2 = F1 * np.sin(np.radians(theta1)) + F2 * np.sin(np.radians(theta2)) - R * np.sin(np.radians(alpha))
        return [eq1, eq2]
    # 初始猜测的分力大小
    initial_guess = [R / 2, R / 2]
    # 使用 fsolve 求解
    F1, F2 = fsolve(equations, initial_guess)

    return F1, F2


# 计算截面的抗弯模量
def section_modulus(B, H, t):
    """
    计算空心矩形截面的抗弯模量 Sx 和 Sy

    参数:
    B -- 外部宽度 (mm)
    H -- 外部高度 (mm)
    t -- 壁厚 (mm)

    返回:
    (Sx, Sy) -- x 和 y 方向的抗弯模量 (mm^3)
    """
    if B <= 2 * t or H <= 2 * t:
        raise ValueError("壁厚过大，导致中空部分消失，请检查输入值。")

    # 计算惯性矩
    Ix = (B * H ** 3 / 12) - ((B - 2 * t) * (H - 2 * t) ** 3 / 12)
    Iy = (H * B ** 3 / 12) - ((H - 2 * t) * (B - 2 * t) ** 3 / 12)

    # 计算最外纤维到中性轴的距离
    cx = H / 2
    cy = B / 2

    # 计算抗弯模量
    Sx = Ix / cx
    Sy = Iy / cy

    return Sx, Sy

def gong_section_modulus(H, bf, tf, tw):
    """
    计算工字钢截面的抗弯模量 Sx 和 Sy（单位 mm³）

    参数:
    H  -- 总高度 (mm)
    bf -- 翼缘宽度 (mm)
    tf -- 翼缘厚度 (mm)
    tw -- 腹板厚度 (mm)

    返回:
    (Sx, Sy) -- 分别为绕 X 轴 和 Y 轴 的抗弯模量 (mm^3)
    """
    # 安全检查
    if H <= 2 * tf or bf <= 0 or tf <= 0 or tw <= 0:
        raise ValueError("输入的几何参数不合法。请确认 H > 2*tf 且其他参数 > 0")

    # 腹板高度
    hw = H - 2 * tf
    Ix = (bf * H**3 / 12) - ((bf - tw) * hw**3 / 12)
    Iy = 2 * (tf * bf**3 / 12) + (tw * hw**3 / 12)

    # 中性轴到最远纤维距离
    cx = H / 2
    cy = bf / 2

    # 抗弯模量
    Sx = Ix / cx
    Sy = Iy / cy

    return Sx, Sy


def section_area(B, H, t):
    # 计算空心矩形截面的截面面积
    A = B * H - (B - 2 * t) * (H - 2 * t)
    return A

def gong_section_area(h, bf, tf, tw):
    """
    计算工字型截面的截面面积

    参数:
    h -- 工字型截面的总高度 (mm)
    bf -- 翼缘宽度 (mm)
    tf -- 翼缘厚度 (mm)
    tw -- 腹板厚度 (mm)

    返回:
    截面面积 (mm^2)
    """
    # 腹板高度
    hw = h - 2 * tf

    # 面积计算
    area = 2 * (bf * tf) + hw * tw
    return area


def cal_changxi_ratio(l, I1, I2, A, mu11 = 2, mu12 = 0.7, mu21 = 1.0, mu22 = 1.0):
    """
    计算结构稳定性关键参数长细比：
    - 等效长度 lc1, lc2
    - 回转半径 r1, r2
    - 长细比 λ1, λ2

    l : 力臂
    I1, I2 : float
        主轴方向的惯性矩（单位 mm^4）
    A : float
        截面面积（单位 mm^2）
    """
    # mu11, mu12, mu21, mu22 = 2, 0.7, 1, 1
    # 等效长度
    lc1 = mu11 * mu21 * l
    lc2 = mu12 * mu22 * l

    # 回转半径
    r1 = math.sqrt(I1 / A)
    r2 = math.sqrt(I2 / A)

    # 长细比
    lambda1 = lc1 / r1
    lambda2 = lc2 / r2

    return lambda1, lambda2



def calculate_stress(diao_G, arm_angle_deg, waibai = 0, cebai = 0):
    elephant_trunk_G = 7 # 象鼻梁自重
    arm_G = 15  # 臂架自重
    pinghengliang_G = 3
    phl_peizhong = 6
    renzijia_G = 9
    v_min_per = 30  # 起升速度 (单位：m/min)
    beta = 0.34  # 动载系数β
    phi = 1.1    # 动载系数φ
    # 将速度转换为m/s
    v = v_min_per / 60

    # 计算动载系数
    dynamic_load_factor = phi + (beta * v) / 2
    # 计算吊重
    diao_F = dynamic_load_factor * diao_G

    # 计算象鼻梁大拉杆的角度
    big_rod_angle, elephant_trunk_angle, small_rod_angle, rack_angle = get_part_angles(arm_angle_deg)
    # 头部滑轮铰点合力
    F1_tou_angle = -waibai  # F1的方向
    F2_tou_angle = 90 + elephant_trunk_angle  # F2的方向，
    # 计算合力的大小和方向
    resultant_magnitude_tou, resultant_angle_tou = calculate_resultant_force(diao_F, F1_tou_angle, diao_F, F2_tou_angle)

    # 尾部滑轮铰点合力
    F1_wei_angle = elephant_trunk_angle - 90  # F1的方向
    F2_wei_angle = 90 - big_rod_angle
    # 计算合力的大小和方向
    resultant_magnitude_wei, resultant_angle_wei = calculate_resultant_force(diao_F, F1_wei_angle, diao_F, F2_wei_angle)

    # 大拉杆拉力
    F_big_rod = (elephant_trunk_G * cos(radians(elephant_trunk_angle)) * 3200 +
                 resultant_magnitude_tou * cos(radians(resultant_angle_tou - elephant_trunk_angle)) * 9400 -
                 resultant_magnitude_wei * cos(
                radians(elephant_trunk_angle - resultant_angle_wei)) * 4600) / (
                        cos(radians(90 - big_rod_angle - elephant_trunk_angle)) * 4600)
    # 象鼻梁尾部综合总合力
    magnitude_elephant, angle_elephant = calculate_resultant_force(F_big_rod, 90 - big_rod_angle,
                                                                   resultant_magnitude_wei, resultant_angle_wei)
    # 吊重的侧向偏摆力
    F_side = diao_F * tan(radians(cebai))

    # 通过节点A的受力分解得到F1、F2，F1表示杆1的大小，F2表示杆2的受力大小;15为杆1杆2的夹角;-2是因为象鼻梁轴线与  其头部与臂架象鼻梁铰点连线的角度为2°
    # R = resultant_magnitude_tou  # 合力大小
    # alpha = resultant_angle_tou  # 合力方向
    # theta1 = 90 + elephant_trunk_angle # 第一个分力方向
    # theta2 = 15 + elephant_trunk_angle - 90  # 第二个分力方向
    F1, F2 = calculate_forces(resultant_magnitude_tou, resultant_angle_tou, 90 + elephant_trunk_angle - 2,
                              15 + elephant_trunk_angle - 90)  #
    # 通过节点B求得杆4始终受拉，杆5始终受压，F4表示杆4的大小，F5表示杆5的受力大小;
    # R = magnitude_elephant  # 合力大小
    # alpha = angle_elephant  # 合力方向
    # theta1 = 90 + elephant_trunk_angle -33# 第一个分力方向
    # theta2 = elephant_trunk_angle - 90 + 6  # 第二个分力方向
    F4, F5 = calculate_forces(magnitude_elephant, angle_elephant, 57 + elephant_trunk_angle, elephant_trunk_angle - 90 + 6)

    # 根据杆4与杆2力的大小与方向求得杆3的力的大小与方向
    # F4_magnitude = F4
    # F4_angle = 57 + elephant_trunk_angle
    # F2_magnitude = F2
    # F2_angle = 15 + elephant_trunk_angle - 90
    F3, angle_F3 = calculate_resultant_force(F4, 57 + elephant_trunk_angle , F2, 15 + elephant_trunk_angle - 90)
    # 杆1综合应力
    B_rod1 = 817  # mm
    H_rod1 = 362  # mm
    t_rod1 = 8   # mm
    A_rod1 = section_area(B_rod1, H_rod1, t_rod1)
    Wx_rod1, Wy_rod1= section_modulus(B_rod1, H_rod1, t_rod1)
    Ix_rod1 = Wx_rod1 * H_rod1
    Iy_rod1 = Wy_rod1 * B_rod1
    l_libi_rod1 = 9.4
    M1_rod1 = F_side * l_libi_rod1  # 单位为t.m
    M2_rod1 = F1 * 0.15  # 单位为t.m
    sigma_rod1 = (M1_rod1 * 9.8 * 10**6) / Wy_rod1 + (M2_rod1 * 9.8 * 10**6) / Wx_rod1 + (F1 * 9800) / A_rod1
    # 杆1稳定性计算
    lambda1_rod1, lambda2_rod1 = cal_changxi_ratio(l_libi_rod1 * 1000, Ix_rod1, Iy_rod1, A_rod1 )
    coefficient_1_rod1, coefficient_2_rod1 = 0.661, 0.983 # 查表得稳定性系数，并取偏保守值0.661
    sigma_rod1_Stability = (M1_rod1 * 9.8 * 10 ** 6) / Wy_rod1 + (M2_rod1 * 9.8 * 10 ** 6) / Wx_rod1 + (F1 * 9800)  / (A_rod1 * coefficient_1_rod1)
    # 杆2综合应力
    H_rod2 = 212  # mm
    bf_rod2 = 280  # mm
    tf_rod2 = 8   # mm
    tw_rod2 = 10
    A_rod2 = gong_section_area(H_rod2, bf_rod2, tf_rod2, tw_rod2)
    sigma_rod2 = F2 * 9800 / A_rod2
    # 杆3综合应力
    H_rod3 = 204  # mm
    bf_rod3 = 220  # mm
    tf_rod3 = 8   # mm
    tw_rod3 = 8
    A_rod3 = gong_section_area(H_rod3, bf_rod3, tf_rod3, tw_rod3)
    # F3与杆3轴线夹角角度绝对值
    l_libi_rod3 = 1.7
    angle_difference = abs(elephant_trunk_angle - angle_F3)
    M1_rod3 = F3 * sin(radians(angle_difference)) * l_libi_rod3
    F_zhouli = F3 * cos(radians(angle_difference))
    Sx, Sy = gong_section_modulus(H_rod3, bf_rod3, tf_rod3, tw_rod3)
    Ix_rod3 = Sx * H_rod3
    Iy_rod3 = Sy * bf_rod3
    # sigma_rod3 = F_zhouli * 9800 / A_rod3 + M1_rod3 * 9.8 * 10**6 / Sy
    sigma_rod3 = F_zhouli * 9800 / A_rod3
    # 杆3稳定性计算
    lambda1_rod3, lambda2_rod3 = cal_changxi_ratio(l_libi_rod3 * 1000, Ix_rod3, Iy_rod3, A_rod3)
    coefficient_1_rod3, coefficient_2_rod3 = 0.946, 0.987  # 查表得稳定性系数，并取偏保守值0.946
    # sigma_rod3_Stability = F_zhouli * 9800 / (A_rod3 * coefficient_1_rod3) + M1_rod3 * 9.8 * 10**6 / Sy
    sigma_rod3_Stability = F_zhouli * 9800 / (A_rod3 * coefficient_1_rod3)
    # 杆4综合应力
    H_rod4 = 212  # mm
    bf_rod4 = 220  # mm
    tf_rod4 = 8   # mm
    tw_rod4 = 8
    A_rod4 = gong_section_area(H_rod4, bf_rod4, tf_rod4, tw_rod4)
    sigma_rod4 = F4 * 9800 / A_rod4
    # 杆5综合应力
    B_rod5 = 817  # mm
    H_rod5 = 362  # mm
    t_rod5 = 8   # mm
    A_rod5 = section_area(B_rod5, H_rod5, t_rod5)
    Wx_rod5, Wy_rod5= section_modulus(B_rod5, H_rod5, t_rod5)
    Ix_rod5 = Wx_rod5 * H_rod5
    Iy_rod5 = Wy_rod5 * B_rod5
    l_libi_rod5 = 0.2
    M1_rod5 = F5 * l_libi_rod5   # 单位为t.m
    sigma_rod5 =  (M1_rod5 * 9.8 * 10**6) / Wx_rod5 + (F5 * 9800)  / A_rod5
    # 杆5稳定性计算
    lambda1_rod5, lambda2_rod5 = cal_changxi_ratio(l_libi_rod5 * 1000, Ix_rod5, Iy_rod3, A_rod5)
    coefficient_1_rod5, coefficient_2_rod5 = 0.999, 0.996
    sigma_rod5_Stability = (M1_rod5 * 9.8 * 10**6) / Wx_rod5 + (F5 * 9800) * cos(radians(6)) / (A_rod5 * coefficient_2_rod5)
    # print(f"杆1的综合应力为{sigma_rod1}")
    # print(f"杆1的稳定性计算综合应力为{sigma_rod1_Stability}")
    # print(f"杆2的综合应力为{sigma_rod2}")
    # print(f"杆3的综合应力为{sigma_rod3}")
    # print(f"杆3的稳定性计算综合应力为{sigma_rod3_Stability}")
    # print(f"杆4的综合应力为{sigma_rod4}")
    # print(f"杆5的综合应力为{sigma_rod5}")
    # print(f"杆5的稳定性计算综合应力为{sigma_rod5_Stability}")
    # 大拉杆综合应力
    B_big_rod = 480  # mm
    H_big_rod = 200  # mm  
    t_big_rod = 8   # mm
    A_big_rod = section_area(B_big_rod, H_big_rod, t_big_rod)
    sigma_big_rod = F_big_rod * 9800 / A_big_rod
    # print(f"大拉杆最不利截面处的综合应力为{sigma_big_rod}")

    # 臂架综合应力,根据力的分解
    F_arm_x = magnitude_elephant * sin(radians(angle_elephant)) + resultant_magnitude_tou * sin(radians(resultant_angle_tou))
    F_arm_y = elephant_trunk_G + magnitude_elephant * cos(radians(angle_elephant)) + resultant_magnitude_tou * cos(radians(resultant_angle_tou))
    # 将两分力求合力
    F_arm, angel_arm = calculate_resultant_force(F_arm_x, 90, F_arm_y, 0)
    # 臂架与angel_arm合力的夹角
    angel_arm_project = 90 - arm_angle_deg - angel_arm
    # 截面正向力
    F_arm_rod = F_arm * cos(radians(angel_arm_project)) + arm_G * cos(radians(90 - arm_angle_deg))
    # 沿轴1的弯矩
    arm_M1 = F_arm * sin(radians(angel_arm_project)) * 19 + arm_G * sin(radians(90 - arm_angle_deg)) * 5.2
    # 沿轴2的弯矩
    arm_M2 = F_side * 19
    # 臂架综合应力
    B_arm = 1416 # mm
    H_arm = 1266  # mm
    t_arm = 8  # mm
    A_arm = section_area(B_arm, H_arm, t_arm)
    Wx_arm, Wy_arm = section_modulus(B_arm, H_arm, t_arm)
    Ix_arm = Wx_arm * H_arm
    Iy_arm = Wy_arm * B_arm
    sigma_arm = (arm_M1 * 9.8 * 10 ** 6) / Wx_arm + (arm_M2 * 9.8 * 10 ** 6) / Wy_arm + (F_arm_rod * 9800)  / A_arm
    # print(f"臂架最不利截面处的综合应力为{sigma_arm}")
    lambda1_arm, lambda2_arm = cal_changxi_ratio(19000, Ix_arm, Iy_arm, A_arm, 2, 2, 1.4, 1.15)
    coefficient_1_arm, coefficient_2_arm = 0.739, 0.838 # 查表得稳定性系数，并取偏保守值0.739
    sigma_arm_Stability = (arm_M1 * 9.8 * 10 ** 6) / Wx_arm + (arm_M2 * 9.8 * 10 ** 6) / Wy_arm + (F_arm_rod * 9800)  / (A_arm * coefficient_1_arm) 
    # print(f"臂架最不利截面处的稳定性计算综合应力为{sigma_arm_Stability}")
    # 人字架综合应力计算
    #小拉杆受力与大拉杆底部铰点截面的夹角
    jia_angel = 90 - (small_rod_angle + big_rod_angle)
    # 对大拉杆底部铰点截面取矩，求小拉杆受力
    F_xlg = pinghengliang_G * cos(radians(big_rod_angle)) * 1.7 + phl_peizhong * cos(radians(big_rod_angle)) * 3.2 / cos(radians(jia_angel)) * 2.3
    # 根据力的分解求平衡梁的大小与方向
    F_phl_x = F_xlg * cos(radians(small_rod_angle))
    F_phl_y = F_xlg * sin(radians(small_rod_angle))
    # 对臂架底部铰点截面取矩，求齿条架受力
    F_chitiao = (F_arm * sin(radians(angel_arm_project)) * 24 + arm_G * sin(radians(90 - arm_angle_deg)) * 10
                 - F_xlg * cos(radians(90 - small_rod_angle - arm_angle_deg )) * 5.8 / ( 5.8 * cos(radians(90 - rack_angle - arm_angle_deg ))))
    # 滑轮受力
    F_hualun1_magnitude, F_hualun1_angle = calculate_resultant_force(diao_F,2.7, diao_F, - 90 - big_rod_angle)
    F_hualun2_magnitude, F_hualun2_angle = calculate_resultant_force(diao_F, 13.1, diao_F, - 90 - big_rod_angle)
    # 将滑轮力分解
    Fx_hualun = - (F_hualun1_magnitude * sin(radians(F_hualun1_angle)) + F_hualun2_magnitude * sin(radians(F_hualun2_angle)))
    Fy_hualun = F_hualun1_magnitude * cos(radians(F_hualun1_angle)) + F_hualun2_magnitude * cos(radians(F_hualun2_angle))
    # 大拉杆受力分解
    Fx_bid_rod = F_big_rod * cos(radians(big_rod_angle))
    Fy_big_rod = F_big_rod * sin(radians(big_rod_angle))
    # 人字架受力分析,截面所受正向力
    F_rzj_axial = Fy_hualun + F_phl_y + renzijia_G - Fy_big_rod
    M_rzj = (Fx_hualun + F_phl_x + Fx_bid_rod) * 8.6 + F_chitiao * 6.6
    B_rzj = 3200  # mm
    H_rzj = 3200  # mm
    t_rzj = 16   # mm
    A_rzj = section_area(B_rzj, H_rzj, t_rzj)
    Wx_rzj, Wy_rzj = section_modulus(B_rzj, H_rzj, t_rzj)
    sigma_rzj = (M_rzj * 9.8 * 10 ** 6) / Wx_rzj + (F_rzj_axial * 9800) / A_rzj

    # print(f"人字架最不利截面处的综合应力为{sigma_rzj}")
    return (
        round(sigma_rod1, 2),
        round(sigma_rod2, 2),
        round(sigma_rod3, 2),
        round(sigma_rod4, 2),
        round(sigma_rod5, 2),
        round(sigma_arm, 2),
        round(sigma_rzj, 2),
        round(sigma_big_rod+10, 2),
        round(sigma_rod1_Stability, 2),
        round(sigma_rod3_Stability, 2),
        round(sigma_rod5_Stability, 2),
        round(sigma_arm_Stability, 2)
    )




