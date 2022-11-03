# RTKLIB TUTORIAL 02 - Understanding Satellite Navigation

## 1 卫星导航的组成 How it works

* 空间星座 Satellite Systems
* 地面监控 Control Segment
* 用户设备 Receivers

### 1.1 卫星星座 Satellite Systems

一颗卫星主要包括：（和定位有关的部分）

* 无线电收发装置 - 收发信号
* 原子钟 - 时间标准
* 推进系统 - 轨道机动

接收信号：地面控制部分的计算结果和机动指令

播发信号：调制在载波信号上的伪码（测距码）和数据码（导航电文信息）


### 1.2 地面控制 Control Segment

固定站点，接收卫星信号，

* 计算卫星钟差、星历（轨道、位置、速度等参数）
* 计算大气层延时等修正参数
* 将以上导航电文注入卫星


### 1.3 接收机 Receivers

接收卫星信号，定位


## 2 定位 Positioning

### 2.1 单点定位 Single Point Positioning

利用伪距计算用户的绝对位置。

伪距方程：

$$
\rho = r + c (\delta t_{user} - \delta t_{sat}) + I + T + \epsilon
$$

定位解算：

* 牛顿迭代与线性化：在点$x_{k-1}$处一阶泰勒展开（线性化），计算$\Delta x$，更新得到$x_{k} = x_{k-1} + \Delta x$。这意味着每次迭代都会更新线性化点

* 计算变化量$\Delta x$：使用最小二乘求解

    $$
    G \cdot \Delta x = b
    $$

    $G$为几何矩阵，$b$为残差：$b = \rho -(r + c (\delta t_{user} - \delta t_{sat}) + I + T)$。

    RTKlib的函数`lsq()`中使用$Ax=y$表示。


### 2.2 差分定位 Differential GPS

利用基站(rover)，消除伪距观测量的公共误差。差分计算用户(user)到基站(rover)的相对位置（基线向量）。

* 单差：站间，消除卫星钟差、星历误差，（短基线）消除大气延时误差，

$$
\rho ^{i}_{ur} = r^{i}_{ur} + c \delta t_{ur} + \epsilon ^{i}_{ur}
$$

* 双差：站间星间，进一步消除接收机钟差，

$$
\rho ^{ij}_{ur} = r^{ij}_{ur} + \epsilon ^{ij}_{ur}
$$


### 2.3 实时动态载波相位差分 RTK

单差：

$$
\phi ^{i}_{ur} = \lambda ^{-1} r^{i}_{ur} + f \delta t_{ur} + N^{i}_{ur} + \epsilon ^{i}_{\phi, ur}
$$

$\lambda$为波长，$f$为载波频率，$N$为整周模糊度。


## 参考 Reference

谢钢《GPS原理与接收机设计》