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

### 伪距定位 Pseudorange Positioning

单点定位 Single Point Positioning

差分定位 Differential GPS


### 实时载波相位差分 RTK



## 参考 Reference

谢钢《GPS原理与接收机设计》