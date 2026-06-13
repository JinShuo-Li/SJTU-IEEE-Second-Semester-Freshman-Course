# 振动力学基础

振动是物理量在某一平衡值附近反复变化的过程. 位移的机械振动, 电路中的电磁振荡, 分子附近的小振动都可以用同一套数学结构处理. 本章重点是简谐振动, 微振动近似, 振动合成, 阻尼振动和受迫振动.

## 振动的分类

按照是否存在外界周期驱动, 振动可分为:

1. **自由振动**: 系统偏离平衡位置后, 只在自身恢复力作用下振动.
2. **受迫振动**: 系统受到外界周期性驱动力作用.

按照是否存在能量耗散, 可分为:

1. **无阻尼振动**: 总机械能守恒.
2. **阻尼振动**: 机械能逐渐耗散, 振幅随时间衰减.

按照运动方程是否线性, 可分为:

1. **谐振动**: 恢复力与位移成正比, 运动方程线性.
2. **非谐振动**: 恢复力与位移不严格成正比, 但在平衡位置附近常可近似为简谐振动.

一个周期性振动可以分解为一系列频率分立的简谐振动. 因而简谐振动不仅是最简单的振动, 也是分析复杂振动的基本单元.

## 简谐振动的运动学描述

简谐振动的标准表达式为:

$$
x(t)=A\cos(\omega t+\varphi).
$$

其中:

- $A$ 为振幅, 表示离开平衡位置的最大位移.
- $\omega$ 为角频率, 单位为 $\mathrm{rad/s}$.
- $T$ 为周期, $\nu$ 为频率.
- $\varphi$ 为初相位.
- $\omega t+\varphi$ 为 $t$ 时刻相位.

频率, 周期与角频率关系为:

$$
\nu=\frac{1}{T},\qquad \omega=2\pi \nu=\frac{2\pi}{T}.
$$

速度与加速度分别为:

$$
v=\frac{dx}{dt}=-A\omega\sin(\omega t+\varphi),
$$

$$
a=\frac{d^2x}{dt^2}=-A\omega^2\cos(\omega t+\varphi)=-\omega^2x.
$$

因此简谐振动满足微分方程:

$$
\frac{d^2x}{dt^2}+\omega^2x=0.
$$

这个方程的核心物理含义是: 加速度始终指向平衡位置, 且大小与位移成正比.

## 相位与相位差

对于两个简谐振动

$$
x_1=A_1\cos(\omega_1t+\varphi_1),\qquad
x_2=A_2\cos(\omega_2t+\varphi_2),
$$

相位差为:

$$
\Delta\varphi=(\omega_2t+\varphi_2)-(\omega_1t+\varphi_1).
$$

若频率相同, 则相位差不随时间变化:

$$
\Delta\varphi=\varphi_2-\varphi_1.
$$

- 当 $\Delta\varphi=2k\pi$ 时, 两振动同相.
- 当 $\Delta\varphi=(2k+1)\pi$ 时, 两振动反相.
- 若 $0<\Delta\varphi<\pi$, 则 $x_2$ 比 $x_1$ 超前 $\Delta\varphi$.

相位描述的是振动状态. 对简谐振动而言, 相位一旦确定, 位移和速度也随之确定.

## 简谐振动的表示方法

### 解析法

直接使用

$$
x=A\cos(\omega t+\varphi)
$$

来描述振动. 已知 $A,\omega,\varphi$ 即可写出完整运动方程; 反之由运动方程也能读出这些量.

### 曲线法

画出 $x-t$ 曲线. 从图像中可读出:

- 最大偏离量给出振幅 $A$.
- 相邻两个同相位点的时间间隔给出周期 $T$.
- $t=0$ 处的位移和速度方向共同确定初相 $\varphi$.

### 振幅矢量法

令长度为 $A$ 的矢量以角速度 $\omega$ 绕原点逆时针旋转, 其在 $x$ 轴上的投影就是:

$$
x=A\cos(\omega t+\varphi).
$$

振幅矢量法特别适合处理同频率简谐振动的合成.

## 简谐振动的动力学方程

若一维系统所受回复力满足

$$
F=-kx,
$$

由牛顿第二定律得:

$$
m\frac{d^2x}{dt^2}=-kx.
$$

于是:

$$
\frac{d^2x}{dt^2}+\frac{k}{m}x=0.
$$

所以角频率为:

$$
\omega=\sqrt{\frac{k}{m}},\qquad T=2\pi\sqrt{\frac{m}{k}}.
$$

简谐振动的判断标准不是轨迹形状像不像正弦曲线, 而是系统在平衡位置附近是否满足线性回复力:

$$
F=-kx.
$$

## 弹簧振子的能量

对水平弹簧振子:

$$
x=A\cos(\omega t+\varphi),\qquad v=-A\omega\sin(\omega t+\varphi).
$$

动能为:

$$
E_k=\frac{1}{2}mv^2=\frac{1}{2}mA^2\omega^2\sin^2(\omega t+\varphi).
$$

势能为:

$$
E_p=\frac{1}{2}kx^2=\frac{1}{2}kA^2\cos^2(\omega t+\varphi).
$$

利用 $k=m\omega^2$, 总能量为:

$$
E=E_k+E_p=\frac{1}{2}kA^2=\frac{1}{2}m\omega^2A^2.
$$

动能和势能随时间周期性交换, 但无阻尼情况下总机械能守恒. 能量与振幅平方成正比.

## 常见简谐系统

### 单摆的小角振动

单摆长度为 $l$, 摆角为 $\theta$. 切向运动方程为:

$$
ml\frac{d^2\theta}{dt^2}=-mg\sin\theta.
$$

当 $|\theta|\ll 1$ 时, $\sin\theta\approx\theta$, 得:

$$
\frac{d^2\theta}{dt^2}+\frac{g}{l}\theta=0.
$$

因此:

$$
\omega=\sqrt{\frac{g}{l}},\qquad T=2\pi\sqrt{\frac{l}{g}}.
$$

单摆周期在小角近似下与质量无关, 与振幅近似无关.

### 复摆

复摆是能绕固定水平轴转动的刚体. 设转轴到质心距离为 $r_c$, 刚体对转轴转动惯量为 $J$, 摆角为 $\theta$. 力矩方程为:

$$
J\frac{d^2\theta}{dt^2}=-mgr_c\sin\theta.
$$

小角近似下:

$$
J\frac{d^2\theta}{dt^2}+mgr_c\theta=0.
$$

所以:

$$
\omega=\sqrt{\frac{mgr_c}{J}},\qquad
T=2\pi\sqrt{\frac{J}{mgr_c}}.
$$

### 扭摆

若扭转回复力矩为

$$
M=-\kappa\theta,
$$

转动惯量为 $J$, 则:

$$
J\frac{d^2\theta}{dt^2}+\kappa\theta=0.
$$

于是:

$$
\omega=\sqrt{\frac{\kappa}{J}},\qquad T=2\pi\sqrt{\frac{J}{\kappa}}.
$$

## 微振动的简谐近似

许多系统在平衡位置附近虽然不是严格简谐, 但可近似为简谐振动. 设势能函数为 $U(x)$, 平衡位置为 $x_0$. 在稳定平衡位置有:

$$
\left.\frac{dU}{dx}\right|_{x_0}=0,
\qquad
\left.\frac{d^2U}{dx^2}\right|_{x_0}>0.
$$

在 $x_0$ 附近作泰勒展开:

$$
U(x)\approx U(x_0)+\frac{1}{2}U''(x_0)(x-x_0)^2.
$$

回复力为:

$$
F=-\frac{dU}{dx}\approx -U''(x_0)(x-x_0).
$$

所以等效劲度系数为:

$$
k_{\mathrm{eff}}=U''(x_0),
$$

小振动角频率为:

$$
\omega=\sqrt{\frac{k_{\mathrm{eff}}}{m}}.
$$

这说明简谐近似的本质是势能曲线在稳定平衡点附近近似为抛物线.

## 同方向同频率简谐振动的合成

设两个同方向, 同频率的简谐振动为:

$$
x_1=A_1\cos(\omega t+\varphi_1),\qquad
x_2=A_2\cos(\omega t+\varphi_2).
$$

合振动仍为简谐振动:

$$
x=x_1+x_2=A\cos(\omega t+\varphi).
$$

振幅满足:

$$
A=\sqrt{A_1^2+A_2^2+2A_1A_2\cos(\varphi_2-\varphi_1)}.
$$

初相满足:

$$
\tan\varphi=\frac{A_1\sin\varphi_1+A_2\sin\varphi_2}{A_1\cos\varphi_1+A_2\cos\varphi_2}.
$$

特殊地:

- 同相: $A=A_1+A_2$.
- 反相: $A=|A_1-A_2|$.
- 相位差为 $\pi/2$: $A=\sqrt{A_1^2+A_2^2}$.

## 同方向不同频率简谐振动的合成与拍

若

$$
x_1=A\cos\omega_1t,
\qquad
x_2=A\cos\omega_2t,
$$

则:

$$
x=x_1+x_2=2A\cos\left(\frac{\omega_1-\omega_2}{2}t\right)
\cos\left(\frac{\omega_1+\omega_2}{2}t\right).
$$

当 $\omega_1$ 与 $\omega_2$ 接近时, 合振动表现为高频振动被缓慢变化的包络调制, 这种现象称为拍. 拍频为:

$$
\nu_{\mathrm{beat}}=|\nu_1-\nu_2|.
$$

若用角频率表示, 包络角频率与 $|\omega_1-\omega_2|$ 有关, 但实际每秒强弱变化次数对应频率差 $|\nu_1-\nu_2|$.

## 垂直简谐振动的合成

设两个相互垂直的同频率简谐振动:

$$
x=A\cos\omega t,
\qquad
y=B\cos(\omega t+\delta).
$$

消去 $t$ 可得轨迹方程:

$$
\frac{x^2}{A^2}+\frac{y^2}{B^2}-\frac{2xy}{AB}\cos\delta=\sin^2\delta.
$$

这通常是椭圆. 特殊情形:

1. $\delta=0$ 或 $\pi$: 轨迹为直线.
2. $\delta=\pm\pi/2$ 且 $A=B$: 轨迹为圆.
3. $\delta=\pm\pi/2$ 且 $A\neq B$: 轨迹为正椭圆.

若两个垂直振动频率不同, 合成轨迹称为李萨如图形. 当频率比为有理数时, 轨迹闭合; 当频率比为无理数时, 轨迹不闭合.

## 阻尼振动

实际系统中常存在阻力. 若阻力与速度成正比:

$$
F_d=-b\frac{dx}{dt},
$$

则运动方程为:

$$
m\frac{d^2x}{dt^2}+b\frac{dx}{dt}+kx=0.
$$

令

$$
\omega_0=\sqrt{\frac{k}{m}},\qquad \beta=\frac{b}{2m},
$$

方程化为:

$$
\frac{d^2x}{dt^2}+2\beta\frac{dx}{dt}+\omega_0^2x=0.
$$

### 欠阻尼

当 $\beta<\omega_0$ 时, 系统仍振动, 但振幅指数衰减:

$$
x=A_0e^{-\beta t}\cos(\omega_d t+\varphi),
$$

其中阻尼角频率为:

$$
\omega_d=\sqrt{\omega_0^2-\beta^2}.
$$

### 临界阻尼

当 $\beta=\omega_0$ 时, 系统不振荡, 且以最快方式回到平衡位置:

$$
x=(C_1+C_2t)e^{-\beta t}.
$$

### 过阻尼

当 $\beta>\omega_0$ 时, 系统不振荡, 缓慢回到平衡位置:

$$
x=C_1e^{(-\beta+\sqrt{\beta^2-\omega_0^2})t}+C_2e^{(-\beta-\sqrt{\beta^2-\omega_0^2})t}.
$$

阻尼振动中, 振幅按 $e^{-\beta t}$ 衰减, 能量约按 $e^{-2\beta t}$ 衰减.

## 受迫振动

受迫振动方程为:

$$
m\frac{d^2x}{dt^2}+b\frac{dx}{dt}+kx=F_0\cos\omega t.
$$

写成标准形式:

$$
\frac{d^2x}{dt^2}+2\beta\frac{dx}{dt}+\omega_0^2x=f_0\cos\omega t,
\qquad f_0=\frac{F_0}{m}.
$$

解由暂态项和稳态项组成:

$$
x=x_{\mathrm{transient}}+x_{\mathrm{steady}}.
$$

暂态项随阻尼衰减, 足够长时间后只剩稳态振动:

$$
x=A\cos(\omega t-\delta).
$$

稳态振幅为:

$$
A=\frac{F_0/m}{\sqrt{(\omega_0^2-\omega^2)^2+(2\beta\omega)^2}}.
$$

相位滞后满足:

$$
\tan\delta=\frac{2\beta\omega}{\omega_0^2-\omega^2}.
$$

当驱动力频率较低时, 振动近似与驱动力同相; 当频率很高时, 振动近似反相; 当 $\omega=\omega_0$ 时, 相位滞后接近 $\pi/2$.

## 共振

受迫振动中, 稳态振幅随驱动力频率变化. 当阻尼较小时, 振幅在接近固有频率处达到极大, 称为共振.

位移振幅最大的频率为:

$$
\omega_r=\sqrt{\omega_0^2-2\beta^2}.
$$

当阻尼很小时:

$$
\omega_r\approx\omega_0.
$$

速度振幅

$$
v_{\max}=\omega A
$$

在 $\omega=\omega_0$ 附近最大. 工程中既可以利用共振, 如乐器共鸣和无线电调谐, 也需要避免共振, 如桥梁, 建筑物和机械结构的振动破坏.

## 本章公式总表

| 内容 | 公式 |
| --- | --- |
| 简谐振动 | $x=A\cos(\omega t+\varphi)$ |
| 速度 | $v=-A\omega\sin(\omega t+\varphi)$ |
| 加速度 | $a=-\omega^2x$ |
| 周期 | $T=2\pi/\omega$ |
| 弹簧振子 | $\omega=\sqrt{k/m}$ |
| 单摆小角振动 | $T=2\pi\sqrt{l/g}$ |
| 复摆 | $T=2\pi\sqrt{J/(mgr_c)}$ |
| 总能量 | $E=\frac12 kA^2$ |
| 同频合成振幅 | $A=\sqrt{A_1^2+A_2^2+2A_1A_2\cos\Delta\varphi}$ |
| 拍频 | $\nu_{\mathrm{beat}}=|\nu_1-\nu_2|$ |
| 阻尼频率 | $\omega_d=\sqrt{\omega_0^2-\beta^2}$ |
| 受迫振幅 | $A=\dfrac{F_0/m}{\sqrt{(\omega_0^2-\omega^2)^2+(2\beta\omega)^2}}$ |
