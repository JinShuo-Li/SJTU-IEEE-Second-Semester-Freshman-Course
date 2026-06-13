# 线性代数 II

线性空间、Jordan 标准型、矩阵分解与矩阵微积分

目录:

- [线性空间与线性变换](#线性空间与线性变换)
  - [线性空间](#线性空间)
    - [线性组合、张成空间与基](#线性组合-张成空间与基)
  - [子空间及其运算](#子空间及其运算)
    - [子空间的和与交](#子空间的和与交)
    - [直和](#直和)
    - [维数公式](#维数公式)
  - [补子空间、正交补与投影](#补子空间-正交补与投影)
  - [线性变换与矩阵](#线性变换与矩阵)
    - [矩阵表示](#矩阵表示)
    - [基变换与相似](#基变换与相似)
  - [伴随变换](#伴随变换)
  - [投影变换](#投影变换)
  - [矩阵的张量积](#矩阵的张量积)
    - [向量化公式](#向量化公式)
- [方阵的 Jordan 标准型与零化多项式](#方阵的-jordan-标准型与零化多项式)
  - [特征值、特征向量与对角化](#特征值-特征向量与对角化)
  - [广义特征向量与 Jordan 链](#广义特征向量与-jordan-链)
  - [幂零矩阵的结构](#幂零矩阵的结构)
  - [Jordan 标准型](#jordan-标准型)
    - [Jordan 标准型的直接意义](#jordan-标准型的直接意义)
  - [零化多项式与最小多项式](#零化多项式与最小多项式)
    - [零化多项式的用法](#零化多项式的用法)
  - [Jordan 标准型的求法](#jordan-标准型的求法)
    - [步骤一：求特征值与代数重数](#步骤一-求特征值与代数重数)
    - [步骤二：用核维数确定块数](#步骤二-用核维数确定块数)
    - [步骤三：构造 Jordan 链](#步骤三-构造-jordan-链)
- [矩阵分解及其应用](#矩阵分解及其应用)
  - [矩阵的满秩分解](#矩阵的满秩分解)
    - [用初等变换求满秩分解](#用初等变换求满秩分解)
    - [满秩分解的意义](#满秩分解的意义)
  - [LU 分解](#lu-分解)
  - [QR 分解](#qr-分解)
    - [QR 分解与最小二乘](#qr-分解与最小二乘)
  - [谱分解](#谱分解)
    - [谱分解的几何意义](#谱分解的几何意义)
  - [Schur 三角分解](#schur-三角分解)
  - [奇异值分解](#奇异值分解)
    - [SVD 的几何意义](#svd-的几何意义)
    - [Moore--Penrose 伪逆](#moore--penrose-伪逆)
    - [低秩逼近](#低秩逼近)
- [向量范数](#向量范数)
  - [范数的定义](#范数的定义)
  - [范数的等价性](#范数的等价性)
    - [常见范数之间的估计](#常见范数之间的估计)
  - [由范数诱导的距离与收敛](#由范数诱导的距离与收敛)
  - [矩阵范数与诱导范数](#矩阵范数与诱导范数)
  - [范数与条件数](#范数与条件数)
- [梯度矩阵与矩阵微积分](#梯度矩阵与矩阵微积分)
  - [标量函数的梯度](#标量函数的梯度)
  - [Jacobian 矩阵](#jacobian-矩阵)
  - [Hessian 矩阵](#hessian-矩阵)
  - [矩阵变量函数的微分](#矩阵变量函数的微分)
    - [迹技巧](#迹技巧)
  - [常用矩阵求导公式](#常用矩阵求导公式)
  - [链式法则](#链式法则)
  - [二次优化中的梯度与 Hessian](#二次优化中的梯度与-hessian)

## 线性空间与线性变换

线性代数的对象并不只是列向量和矩阵。更本质的对象是线性空间以及线性空间之间保持加法和数乘结构的映射。矩阵是线性变换在选定基以后得到的坐标表达；当基改变时，矩阵会改变，但线性变换自身不变。因此，本章先从抽象线性空间出发，再解释矩阵如何作为线性变换的具体表示。

### 线性空间

**定义（线性空间）**:

设 $\mathbb F$ 是一个数域，$V$ 是一个非空集合。如果在 $V$ 上定义了向量加法
$$
+:V\times V\to V,
$$
以及数乘
$$
\mathbb F\times V\to V,
$$
并且对任意 $\alpha,\beta,\gamma\in V$ 与任意 $a,b\in\mathbb F$ 满足：
1. $\alpha+\beta=\beta+\alpha$；
2. $(\alpha+\beta)+\gamma=\alpha+(\beta+\gamma)$；
3. 存在零向量 $0\in V$，使得 $\alpha+0=\alpha$；
4. 对每个 $\alpha\in V$，存在 $-\alpha\in V$，使得 $\alpha+(-\alpha)=0$；
5. $a(b\alpha)=(ab)\alpha$；
6. $1\alpha=\alpha$；
7. $(a+b)\alpha=a\alpha+b\alpha$；
8. $a(\alpha+\beta)=a\alpha+a\beta$，

则称 $V$ 是数域 $\mathbb F$ 上的线性空间。$V$ 中的元素称为向量，$\mathbb F$ 中的元素称为纯量或数。

这些公理的意义是：在 $V$ 中可以做线性组合，并且线性组合的代数规则与普通向量一致。只要可以稳定地讨论
$$
a_1\alpha_1+\cdots+a_k\alpha_k,
$$
就可以谈论线性相关、基、维数、线性变换等概念。

**例**:

以下对象都是线性空间：
1. $\mathbb F^n$，即所有 $n$ 维列向量构成的集合；
2. $\mathbb F^{m\times n}$，即所有 $m\times n$ 矩阵构成的集合；
3. $\mathbb F[x]_{\le n}$，即次数不超过 $n$ 的多项式构成的集合；
4. $C[a,b]$，即区间 $[a,b]$ 上连续函数构成的实线性空间；
5. 微分方程 $y''+p(x)y'+q(x)y=0$ 的所有解构成的集合。

其中最后一个例子体现了线性空间概念的价值：向量不一定是几何箭头，也可以是函数或方程的解。

#### 线性组合、张成空间与基

**定义（线性组合与张成空间）**:

设 $\alpha_1,\ldots,\alpha_k\in V$。形如
$$
a_1\alpha_1+\cdots+a_k\alpha_k,\qquad a_i\in\mathbb F,
$$
的向量称为 $\alpha_1,\ldots,\alpha_k$ 的一个线性组合。所有这样的线性组合构成的集合记为
$$
\operatorname{span}(\alpha_1,
\ldots,\alpha_k).
$$

**定义（线性相关与线性无关）**:

向量组 $\alpha_1,\ldots,\alpha_k$ 称为线性相关，如果存在不全为零的数 $a_1,
\ldots,a_k\in\mathbb F$，使得
$$
a_1\alpha_1+\cdots+a_k\alpha_k=0.
$$
否则称为线性无关。

线性无关意味着每个向量都不能由其余向量线性表出。线性相关意味着向量组中有冗余信息。基就是既无冗余又足够生成整个空间的向量组。

**定义（基与维数）**:

若向量组 $\alpha_1,\ldots,\alpha_n$ 线性无关，且
$$
V=\operatorname{span}(\alpha_1,
\ldots,\alpha_n),
$$
则称它是 $V$ 的一组基。有限维线性空间的任意两组基所含向量个数相同，这个共同的数称为 $V$ 的维数，记为 $\dim V$。

**定理（基扩充定理）**:

设 $V$ 是有限维线性空间。若 $\alpha_1,
\ldots,\alpha_r$ 是 $V$ 中的线性无关组，则可以添加若干向量，使其成为 $V$ 的一组基。

**证明**:

若 $\operatorname{span}(\alpha_1,
\ldots,\alpha_r)=V$，则已经是一组基。否则，存在 $\beta_1\in V$ 不属于该张成空间。此时 $\alpha_1,
\ldots,\alpha_r,\beta_1$ 仍线性无关；否则 $\beta_1$ 可由 $\alpha_i$ 线性表出，矛盾。若新的张成空间仍非 $V$，继续取 $\beta_2$。由于 $V$ 有有限基，线性无关组的长度不可能无限增加，过程必在有限步后停止，得到一组基。

证明完毕.

**定理（替换定理）**:

设 $\alpha_1,\ldots,\alpha_m$ 线性无关，$\beta_1,\ldots,\beta_n$ 张成 $V$。则 $m\le n$，并且可以从 $\beta_1,
\ldots,\beta_n$ 中删去 $m$ 个向量，把 $\alpha_1,\ldots,\alpha_m$ 换进去后仍张成 $V$。

**证明**:

先把 $\alpha_1$ 用 $\beta_1,
\ldots,\beta_n$ 线性表示。因为 $\alpha_1\ne 0$，至少有一个系数非零，不妨为 $\beta_1$ 的系数非零，于是 $\beta_1$ 可由 $\alpha_1,\beta_2,
\ldots,\beta_n$ 线性表出，所以后者仍张成 $V$。再把 $\alpha_2$ 用新的张成组表示。由于 $\alpha_1,\alpha_2$ 线性无关，$\alpha_2$ 的表达式中至少有一个来自剩余 $\beta$ 的系数非零，否则 $\alpha_2\in\operatorname{span}(\alpha_1)$。于是可继续替换。重复 $m$ 次即可。若 $m>n$，第 $n+1$ 步将无可替换的 $\beta$，导致 $\alpha_{n+1}\in\operatorname{span}(\alpha_1,
\ldots,\alpha_n)$，矛盾。所以 $m\le n$。

证明完毕.

### 子空间及其运算

**定义（子空间）**:

设 $V$ 是 $\mathbb F$ 上的线性空间。若非空子集 $W\subseteq V$ 对 $V$ 中的加法和数乘也构成一个线性空间，则称 $W$ 是 $V$ 的子空间。

**定理（子空间判别法）**:

设 $W\subseteq V$ 非空。则 $W$ 是 $V$ 的子空间，当且仅当对任意 $u,v\in W$ 与任意 $a,b\in\mathbb F$，都有
$$
au+bv\in W.
$$
等价地，只需验证 $W$ 对加法和数乘封闭。

**证明**:

若 $W$ 是子空间，则线性组合封闭显然成立。反过来，若 $W$ 对任意二元线性组合封闭，由非空性取 $w\in W$，令 $a=0$ 得 $0\in W$；令 $a=-1$ 得 $-w\in W$；令 $a=b=1$ 得加法封闭；令 $b=0$ 得数乘封闭。其余线性空间公理继承自 $V$，故 $W$ 为子空间。

证明完毕.

#### 子空间的和与交

**定义（子空间的和与交）**:

设 $U,W$ 是 $V$ 的子空间。定义
$$
U+W=\{u+w:u\in U,
 w\in W\},
$$
以及
$$
U\cap W=\{v\in V:v\in U \text{ 且 } v\in W\}.
$$

**定理**:

若 $U,W$ 是 $V$ 的子空间，则 $U+W$ 与 $U\cap W$ 都是 $V$ 的子空间。

**证明**:

先看交。显然 $0\in U\cap W$。若 $x,y\in U\cap W$，$a,b\in\mathbb F$，则 $x,y\in U$ 且 $x,y\in W$，由 $U,W$ 的封闭性得 $ax+by\in U$ 且 $ax+by\in W$，所以 $ax+by\in U\cap W$。

再看和。$0=0+0\in U+W$。若 $x=u_1+w_1$，$y=u_2+w_2$，其中 $u_i\in U,w_i\in W$，则
$$
ax+by=(au_1+bu_2)+(aw_1+bw_2).
$$
由封闭性 $au_1+bu_2\in U$，$aw_1+bw_2\in W$，故 $ax+by\in U+W$。

证明完毕.

**注**:

两个子空间的并通常不是子空间。比如 $\mathbb R^2$ 中的 $x$ 轴与 $y$ 轴都是子空间，但它们的并不是子空间，因为 $(1,0)+(0,1)=(1,1)$ 不属于二者的并。

#### 直和

**定义（直和）**:

设 $U,W$ 是 $V$ 的两个子空间。如果 $U\cap W=\{0\}$，则称 $U+W$ 是 $U$ 与 $W$ 的直和，记为
$$
U\oplus W.
$$
若进一步 $V=U+W$，则称 $V$ 是 $U$ 与 $W$ 的直和，记为
$$
V=U\oplus W.
$$

**定理（直和的等价刻画）**:

设 $U,W$ 是 $V$ 的子空间。以下命题等价：
1. $U+W$ 是直和，即 $U\cap W=\{0\}$；
2. $U+W$ 中每个向量 $v$ 都存在唯一分解
$$
  v=u+w,
  \qquad u\in U,
  \ w\in W;
$$
3. 零向量在 $U+W$ 中的分解唯一；
4. $U$ 的一组基与 $W$ 的一组基合起来构成 $U+W$ 的一组基。

**证明**:

$(1)\Rightarrow(2)$：存在性来自 $v\in U+W$ 的定义。若
$$
v=u+w=u'+w',
$$
则
$$
u-u'=w'-w.
$$
左边属于 $U$，右边属于 $W$，故它属于 $U\cap W=\{0\}$，所以 $u=u'$ 且 $w=w'$。

$(2)\Rightarrow(3)$：零向量作为 $0=0+0$ 有一组分解。由唯一性，它没有其他分解。

$(3)\Rightarrow(1)$：若 $x\in U\cap W$，则
$$
0=0+0=x+(-x),
$$
其中 $x\in U$ 且 $-x\in W$。由零向量分解唯一，得 $x=0$。因此 $U\cap W=\{0\}$。

$(1)\Rightarrow(4)$：设 $\mathcal B_U=(u_1,
\ldots,u_r)$ 是 $U$ 的基，$\mathcal B_W=(w_1,
\ldots,w_s)$ 是 $W$ 的基。任取 $v\in U+W$，有 $v=u+w$，而 $u,w$ 分别可由对应基线性表出，故合并向量组生成 $U+W$。若
$$
a_1u_1+\cdots+a_ru_r+b_1w_1+\cdots+b_sw_s=0,
$$
则
$$
a_1u_1+\cdots+a_ru_r=-(b_1w_1+\cdots+b_sw_s).
$$
左边属于 $U$，右边属于 $W$，故共同向量属于 $U\cap W$，只能为 $0$。于是由两组基的线性无关性得到全部系数为 $0$，合并向量组线性无关。因此它是 $U+W$ 的基。

$(4)\Rightarrow(1)$：若 $x\in U\cap W$，则 $x$ 可以同时用 $U$ 的基和 $W$ 的基表出。于是零向量有表达
$$
x+(-x)=0.
$$
把 $x$ 和 $-x$ 分别展开到两组基上，就得到合并基的一组线性组合为零。由于合并向量组线性无关，所有系数为零，故 $x=0$。

证明完毕.

**定义（多个子空间的直和）**:

设 $W_1,
\ldots,W_k$ 是 $V$ 的子空间。如果每个
$$
v\in W_1+\cdots+W_k
$$
均可唯一写成
$$
v=w_1+\cdots+w_k,
\qquad w_i\in W_i,
$$
则称该和为直和，记为
$$
W_1\oplus\cdots\oplus W_k.
$$

**命题**:

$W_1+\cdots+W_k$ 是直和，当且仅当对每个 $i$，
$$
W_i\cap (W_1+\cdots+W_{i-1}+W_{i+1}+\cdots+W_k)=\{0\}.
$$

**证明**:

若和为直和，取 $x\in W_i\cap\sum_{j\ne i}W_j$，则 $x$ 一方面可以作为第 $i$ 个分量表示，另一方面可以由其他分量表示，所以零向量有两种分解，故 $x=0$。

反过来，若上述交均为零，假设
$$
w_1+\cdots+w_k=w'_1+\cdots+w'_k.
$$
则对任意 $i$，有
$$
w_i-w_i'=\sum_{j\ne i}(w'_j-w_j).
$$
左边属于 $W_i$，右边属于其他子空间之和，所以 $w_i-w_i'=0$。因此分解唯一。

证明完毕.

**注**:

仅有两两相交为零并不足以推出多个子空间的直和。例如在 $\mathbb R^2$ 中，令
$$
W_1=\operatorname{span}(1,0),
\quad W_2=\operatorname{span}(0,1),
\quad W_3=\operatorname{span}(1,1).
$$
任意两个一维子空间交为 $\{0\}$，但
$$
(1,0)+(0,1)-(1,1)=0
$$
给出了零向量的非平凡分解，所以 $W_1+W_2+W_3$ 不是直和。

#### 维数公式

**定理（Grassmann 维数公式）**:

设 $U,W$ 是有限维线性空间 $V$ 的子空间，则
$$
\dim(U+W)=\dim U+\dim W-\dim(U\cap W).
$$

**证明**:

设 $\dim(U\cap W)=r$，取 $U\cap W$ 的一组基
$$
\alpha_1,
\ldots,\alpha_r.
$$
将其扩充为 $U$ 的一组基：
$$
\alpha_1,
\ldots,\alpha_r,
\beta_1,
\ldots,\beta_p,
$$
其中 $p=\dim U-r$；同时扩充为 $W$ 的一组基：
$$
\alpha_1,
\ldots,\alpha_r,
\gamma_1,
\ldots,\gamma_q,
$$
其中 $q=\dim W-r$。

证明
$$
\alpha_1,
\ldots,\alpha_r,
\beta_1,
\ldots,\beta_p,
\gamma_1,
\ldots,\gamma_q
$$
是 $U+W$ 的一组基。生成性显然：任意 $v=u+w\in U+W$，分别把 $u,w$ 展开即可。

再证线性无关。若
$$
\sum_{i=1}^r a_i\alpha_i+
\sum_{j=1}^p b_j\beta_j+
\sum_{l=1}^q c_l\gamma_l=0,
$$
则
$$
\sum_{l=1}^q c_l\gamma_l
=-\left(\sum_{i=1}^r a_i\alpha_i+
\sum_{j=1}^p b_j\beta_j\right).
$$
左边属于 $W$，右边属于 $U$，故该向量属于 $U\cap W$，可写为 $\sum d_i\alpha_i$。于是
$$
\sum_{l=1}^q c_l\gamma_l-
\sum_{i=1}^r d_i\alpha_i=0.
$$
由于 $\alpha_1,
\ldots,\alpha_r,
\gamma_1,
\ldots,\gamma_q$ 是 $W$ 的基，线性无关，所以 $c_l=0$。代回原式，由 $U$ 的基的线性无关性可得 $a_i=b_j=0$。故合并向量组线性无关。

于是
$$
\dim(U+W)=r+p+q
=r+(\dim U-r)+(\dim W-r),
$$
即
$$
\dim(U+W)=\dim U+\dim W-\dim(U\cap W).
$$

证明完毕.

**推论**:

若 $V=U\oplus W$，则
$$
\dim V=\dim U+\dim W.
$$
更一般地，若
$$
V=W_1\oplus\cdots\oplus W_k,
$$
则
$$
\dim V=\sum_{i=1}^k\dim W_i.
$$

### 补子空间、正交补与投影

**定义（补子空间）**:

设 $U$ 是 $V$ 的子空间。若存在子空间 $W\subseteq V$ 使得
$$
V=U\oplus W,
$$
则称 $W$ 是 $U$ 在 $V$ 中的一个补子空间。

**定理（有限维补空间存在性）**:

若 $V$ 有限维，$U$ 是 $V$ 的子空间，则 $U$ 在 $V$ 中至少存在一个补子空间。

**证明**:

取 $U$ 的一组基 $u_1,
\ldots,u_r$，由基扩充定理将其扩充为 $V$ 的一组基
$$
u_1,
\ldots,u_r,w_1,
\ldots,w_s.
$$
令
$$
W=\operatorname{span}(w_1,
\ldots,w_s).
$$
任意 $v\in V$ 都可唯一写成这些基向量的线性组合，故可唯一分解为 $u+w$，其中 $u\in U,w\in W$。所以 $V=U\oplus W$。

证明完毕.

补子空间通常不唯一。例如 $\mathbb R^2$ 中 $x$ 轴的补子空间可以是 $y$ 轴，也可以是任意不等于 $x$ 轴的一维直线。若空间带有内积，则可以用正交性选出一个自然的补空间。

**定义（内积与正交补）**:

实线性空间 $V$ 上的内积是映射 $\left\langle \cdot,\cdot\right\rangle:V\times V\to\mathbb R$，满足对第一变量线性、对称性与正定性。复线性空间上通常取第一变量线性、第二变量共轭线性，且满足
$$
\left\langle x,y\right\rangle=\overline{\left\langle y,x\right\rangle},
\qquad \left\langle x,x\right\rangle>0\ (x\ne0).
$$
若 $U\subseteq V$ 是子空间，定义
$$
U^\perp=\{v\in V:\left\langle u,v\right\rangle=0,
\ \forall u\in U\}.
$$

**定理（正交分解）**:

设 $V$ 是有限维内积空间，$U$ 是 $V$ 的子空间，则
$$
V=U\oplus U^\perp.
$$

**证明**:

先证 $U\cap U^\perp=\{0\}$。若 $v\in U\cap U^\perp$，则 $v\in U$ 且 $v$ 与 $U$ 中所有向量正交，特别地 $\left\langle v,v\right\rangle=0$，由正定性得 $v=0$。

再证 $U+U^\perp=V$。取 $U$ 的一组标准正交基 $e_1,
\ldots,e_k$。任意 $v\in V$，令
$$
u=\sum_{i=1}^k\left\langle v,e_i\right\ranglee_i.
$$
则 $u\in U$。令 $w=v-u$。对任意 $j$，
$$
\left\langle w,e_j\right\rangle
=\left\langle v,e_j\right\rangle-\sum_{i=1}^k\left\langle v,e_i\right\rangle\left\langle e_i,e_j\right\rangle
=\left\langle v,e_j\right\rangle-\left\langle v,e_j\right\rangle=0.
$$
因为 $e_1,
\ldots,e_k$ 张成 $U$，所以 $w\in U^\perp$。于是 $v=u+w\in U+U^\perp$。综上，$V=U\oplus U^\perp$。

证明完毕.

**定义（正交投影）**:

设 $V=U\oplus U^\perp$。对任意 $v\in V$，唯一分解
$$
v=u+w,
\qquad u\in U,
\ w\in U^\perp.
$$
称 $u$ 为 $v$ 到 $U$ 上的正交投影，记为 $P_Uv$。

若 $e_1,
\ldots,e_k$ 是 $U$ 的标准正交基，则
$$
P_Uv=\sum_{i=1}^k\left\langle v,e_i\right\ranglee_i.
$$
如果 $A=[a_1\ \cdots\ a_k]$ 的列向量构成 $U\subseteq\mathbb R^n$ 的一组基，但不一定正交，则 $v$ 到 $U$ 的投影为
$$
P_Uv=A(A^TA)^{-1}A^Tv.
$$
这是因为投影点 $Ax$ 满足残差 $v-Ax$ 与 $U$ 正交，即
$$
A^T(v-Ax)=0.
$$
由此得到正规方程 $A^TAx=A^Tv$。

### 线性变换与矩阵

**定义（线性变换）**:

设 $V,W$ 是 $\mathbb F$ 上的线性空间。映射 $T:V\to W$ 称为线性变换，如果对任意 $\alpha,
\beta\in V$ 与任意 $c\in\mathbb F$，
$$
T(\alpha+\beta)=T\alpha+T\beta,
\qquad T(c\alpha)=cT\alpha.
$$
等价地，对任意 $a,b\in\mathbb F$，
$$
T(a\alpha+b\beta)=aT\alpha+bT\beta.
$$

**定义（核与像）**:

线性变换 $T:V\to W$ 的核与像分别定义为
$$
\operatorname{Ker} T=\{v\in V:T v=0\},
\qquad
\operatorname{Im} T=\{T v:v\in V\}.
$$

**定理**:

$\operatorname{Ker} T$ 是 $V$ 的子空间，$\operatorname{Im} T$ 是 $W$ 的子空间。

**证明**:

若 $x,y\in\operatorname{Ker} T$，$a,b\in\mathbb F$，则
$$
T(ax+by)=aTx+bTy=0,
$$
所以 $ax+by\in\operatorname{Ker} T$。若 $u=Tx$，$v=Ty$ 属于 $\operatorname{Im} T$，则
$$
au+bv=T(ax+by),
$$
所以 $au+bv\in\operatorname{Im} T$。

证明完毕.

**定义（秩与零度）**:

若 $V,W$ 有限维，定义
$$
\operatorname{rank} T=\dim\operatorname{Im} T,
\qquad
\operatorname{nullity} T=\dim\operatorname{Ker} T.
$$

**定理（秩-零度定理）**:

设 $T:V\to W$ 是有限维线性空间之间的线性变换，则
$$
\operatorname{rank} T+\operatorname{nullity} T=\dim V.
$$

**证明**:

取 $\operatorname{Ker} T$ 的一组基 $u_1,
\ldots,u_k$，将它扩充为 $V$ 的一组基
$$
u_1,
\ldots,u_k,v_1,
\ldots,v_r.
$$
证明 $Tv_1,
\ldots,Tv_r$ 是 $\operatorname{Im} T$ 的一组基。

生成性：任意 $y\in\operatorname{Im} T$，有 $y=Tx$。把 $x$ 展开：
$$
x=\sum_{i=1}^k a_i u_i+
\sum_{j=1}^r b_j v_j.
$$
于是
$$
y=Tx=\sum_{j=1}^r b_jTv_j,
$$
故 $Tv_j$ 张成 $\operatorname{Im} T$。

线性无关性：若
$$
\sum_{j=1}^r c_jTv_j=0,
$$
则
$$
T\left(\sum_{j=1}^r c_jv_j\right)=0,
$$
即 $\sum c_jv_j\in\operatorname{Ker} T$，可写为 $\sum d_i u_i$。于是
$$
\sum_{j=1}^r c_jv_j-
\sum_{i=1}^k d_iu_i=0.
$$
由于 $u_i,v_j$ 构成 $V$ 的基，全部系数为零，特别地 $c_j=0$。故 $Tv_1,
\ldots,Tv_r$ 线性无关。

所以 $\operatorname{rank} T=r$，$\operatorname{nullity} T=k$，且 $\dim V=k+r$。

证明完毕.

#### 矩阵表示

设 $\mathcal B=(v_1,
\ldots,v_n)$ 是 $V$ 的一组基，$\mathcal C=(w_1,
\ldots,w_m)$ 是 $W$ 的一组基。对线性变换 $T:V\to W$，每个 $Tv_j$ 可唯一写成
$$
Tv_j=a_{1j}w_1+
\cdots+a_{mj}w_m.
$$
把这些坐标作为矩阵的第 $j$ 列，得到
$$
[T]_{\mathcal C\leftarrow\mathcal B}
=\begin{bmatrix}
a_{11}&a_{12}&\cdots&a_{1n}\\
a_{21}&a_{22}&\cdots&a_{2n}\\
\vdots&\vdots&\ddots&\vdots\\
a_{m1}&a_{m2}&\cdots&a_{mn}
\end{bmatrix}.
$$
若 $x=[v]_{\mathcal B}$ 是 $v$ 在基 $\mathcal B$ 下的坐标，则
$$
[Tv]_{\mathcal C}=[T]_{\mathcal C\leftarrow\mathcal B}[v]_{\mathcal B}.
$$
因此矩阵乘法本质上就是线性变换作用在坐标上的规则。

**命题（复合变换的矩阵）**:

若 $T:U\to V$，$S:V\to W$，并且各空间选定基，则
$$
[S\circ T]=[S][T].
$$

**证明**:

对任意 $x\in U$，
$$
[(S\circ T)x]=[S(Tx)]=[S][Tx]=[S][T][x].
$$
由坐标表示唯一性，得到结论。

证明完毕.

#### 基变换与相似

若 $V$ 的两组基为 $\mathcal B=(v_1,
\ldots,v_n)$ 与 $\mathcal B'=(v'_1,
\ldots,v'_n)$。定义过渡矩阵 $P$ 为
$$
(v'_1,
\ldots,v'_n)=(v_1,
\ldots,v_n)P.
$$
于是任意向量 $x\in V$ 的坐标满足
$$
[x]_{\mathcal B}=P[x]_{\mathcal B'}.
$$
如果 $T:V\to V$ 在两组基下的矩阵分别为 $A$ 与 $A'$，则
$$
A'=P^{-1}AP.
$$

**定义（相似矩阵）**:

若存在可逆矩阵 $P$ 使得
$$
B=P^{-1}AP,
$$
则称 $A$ 与 $B$ 相似。

相似矩阵表示同一个线性变换在不同基下的坐标。因此，相似变换不改变线性变换的本质性质，例如特征值、迹、行列式、秩、最小多项式和 Jordan 块结构。

**命题（相似不变量）**:

若 $B=P^{-1}AP$，则
$$
\operatorname{tr} B=\operatorname{tr} A,
\qquad
\det B=\det A,
\qquad
\operatorname{rank} B=\operatorname{rank} A.
$$
进一步，对任意多项式 $f$，有
$$
f(B)=P^{-1}f(A)P.
$$

**证明**:

迹与行列式由
$$
\operatorname{tr}(P^{-1}AP)=\operatorname{tr}(APP^{-1})=\operatorname{tr} A,
$$
以及
$$
\det(P^{-1}AP)=\det(P^{-1})\det(A)\det(P)=\det A
$$
得到。秩由可逆矩阵左乘右乘不改变秩得到。多项式结论由
$$
B^k=P^{-1}A^kP
$$
逐项相加得到。

证明完毕.

### 伴随变换

**定义（伴随变换）**:

设 $V,W$ 是内积空间，$T:V\to W$ 为线性变换。如果线性变换 $T^*:W\to V$ 满足对任意 $v\in V,w\in W$，
$$
\left\langle Tv,w\right\rangle=\left\langle v,T^*w\right\rangle,
$$
则称 $T^*$ 为 $T$ 的伴随变换。

在有限维内积空间中，伴随变换总存在且唯一。

**定理（伴随的矩阵表示）**:

若在 $V,W$ 中分别取标准正交基，$T$ 的矩阵为 $A$，则 $T^*$ 的矩阵为
$$
A^*=\overline{A}^{\mathsf T}.
$$
在实内积空间中，$A^*=A^{\mathsf T}$。

**证明**:

设 $e_1,
\ldots,e_n$ 是 $V$ 的标准正交基，$f_1,
\ldots,f_m$ 是 $W$ 的标准正交基。若
$$
Te_j=\sum_{i=1}^m a_{ij}f_i,
\qquad
T^*f_i=\sum_{j=1}^n b_{ji}e_j,
$$
则
$$
a_{ij}=\left\langle Te_j,f_i\right\rangle=\left\langle e_j,T^*f_i\right\rangle=\overline{b_{ji}}
$$
在通常的复内积约定下等价于 $b_{ji}=\overline{a_{ij}}$。因此 $T^*$ 的矩阵为 $\overline{A}^{\mathsf T}$。

证明完毕.

**命题（伴随的基本性质）**:

若 $S,T$ 为线性变换，$a,b\in\mathbb F$，则
1. $(aS+bT)^*=\overline a S^*+\overline b T^*$；
2. $(ST)^*=T^*S^*$；
3. $(T^*)^*=T$；
4. $\operatorname{Ker} T=(\operatorname{Im} T^*)^\perp$，$\operatorname{Ker} T^*=(\operatorname{Im} T)^\perp$；
5. 若 $T$ 是单射，则 $T^*$ 是满射；若 $T$ 是满射，则 $T^*$ 是单射。

**证明**:

前三条由伴随定义直接验证。以第二条为例，
$$
\left\langle STx,y\right\rangle=\left\langle Tx,S^*y\right\rangle=\left\langle x,T^*S^*y\right\rangle,
$$
所以 $(ST)^*=T^*S^*$。

第四条中，$x\in\operatorname{Ker} T$ 当且仅当 $Tx=0$，当且仅当对任意 $y$ 有 $\left\langle Tx,y\right\rangle=0$，当且仅当对任意 $y$ 有 $\left\langle x,T^*y\right\rangle=0$，即 $x\perp\operatorname{Im} T^*$。另一式同理。

第五条由秩-零度定理或第四条可得。例如 $T$ 单射意味着 $\operatorname{Ker} T=\{0\}$，所以 $(\operatorname{Im} T^*)^\perp=\{0\}$，有限维下推出 $\operatorname{Im} T^*$ 为全空间。

证明完毕.

### 投影变换

**定义（投影变换）**:

线性变换 $P:V\to V$ 称为投影变换，如果存在直和分解
$$
V=U\oplus W
$$
使得对任意 $v=u+w$，其中 $u\in U,w\in W$，有
$$
Pv=u.
$$
此时称 $P$ 是沿 $W$ 到 $U$ 上的投影。

**定理（投影与幂等）**:

线性变换 $P:V\to V$ 是投影变换，当且仅当
$$
P^2=P.
$$

**证明**:

若 $P$ 是沿 $W$ 到 $U$ 上的投影，则对 $v=u+w$，有 $Pv=u$，而 $u=u+0$，故 $P(Pv)=Pu=u=Pv$，即 $P^2=P$。

反过来，设 $P^2=P$。令
$$
U=\operatorname{Im} P,
\qquad W=\operatorname{Ker} P.
$$
对任意 $v\in V$，
$$
v=Pv+(v-Pv).
$$
其中 $Pv\in U$，并且
$$
P(v-Pv)=Pv-P^2v=0,
$$
故 $v-Pv\in W$。所以 $V=U+W$。若 $x\in U\cap W$，则 $x=Py$ 且 $Px=0$，于是
$$
x=Py=P^2y=P(Py)=Px=0.
$$
故 $V=U\oplus W$，并且 $P$ 正是沿 $W$ 到 $U$ 上的投影。

证明完毕.

**命题（正交投影的矩阵性质）**:

内积空间中的线性变换 $P$ 是某个子空间上的正交投影，当且仅当
$$
P^2=P,
\qquad P^*=P.
$$

**证明**:

若 $P$ 是到 $U$ 上的正交投影，则 $P^2=P$。并且对任意 $x=u+u^\perp$，$y=v+v^\perp$，其中 $u,v\in U$，$u^\perp,v^\perp\in U^\perp$，有
$$
\left\langle Px,y\right\rangle=\left\langle u,v+v^\perp\right\rangle=\left\langle u,v\right\rangle
=\left\langle u+u^\perp,v\right\rangle=\left\langle x,Py\right\rangle,
$$
故 $P^*=P$。

反过来，若 $P^2=P$，则 $V=\operatorname{Im} P\oplus\operatorname{Ker} P$。由 $P^*=P$，若 $u=Px\in\operatorname{Im} P$，$w\in\operatorname{Ker} P$，则
$$
\left\langle u,w\right\rangle=\left\langle Px,w\right\rangle=\left\langle x,P^*w\right\rangle=\left\langle x,Pw\right\rangle=0.
$$
故 $\operatorname{Ker} P=(\operatorname{Im} P)^\perp$，所以 $P$ 是正交投影。

证明完毕.

### 矩阵的张量积

张量积在矩阵论中最常见的具体形式是 Kronecker 积。它用来描述分块结构、双线性映射、矩阵方程、量子态空间以及多维离散线性算子。

**定义（矩阵的张量积或 Kronecker 积）**:

设 $A=(a_{ij})\in\mathbb F^{m\times n}$，$B\in\mathbb F^{p\times q}$。定义
$$
A\otimes B=
\begin{bmatrix}
a_{11}B&a_{12}B&\cdots&a_{1n}B\\
a_{21}B&a_{22}B&\cdots&a_{2n}B\\
\vdots&\vdots&\ddots&\vdots\\
a_{m1}B&a_{m2}B&\cdots&a_{mn}B
\end{bmatrix}
\in\mathbb F^{mp\times nq}.
$$

**例**:

若
$$
A=\begin{bmatrix}1&2\\3&4\end{bmatrix},
\qquad
B=\begin{bmatrix}0&1\\1&0\end{bmatrix},
$$
则
$$
A\otimes B=
\begin{bmatrix}
0&1&0&2\\
1&0&2&0\\
0&3&0&4\\
3&0&4&0
\end{bmatrix}.
$$

**命题（Kronecker 积的基本运算）**:

只要矩阵阶数相容，就有
1. $(A+C)\otimes B=A\otimes B+C\otimes B$；
2. $A\otimes(B+D)=A\otimes B+A\otimes D$；
3. $(\lambda A)\otimes B=A\otimes(\lambda B)=\lambda(A\otimes B)$；
4. $(A\otimes B)^{\mathsf T}=A^{\mathsf T}\otimes B^{\mathsf T}$；
5. $(A\otimes B)^*=A^*\otimes B^*$；
6. $(A\otimes B)(C\otimes D)=(AC)\otimes(BD)$。

**证明**:

前五条直接由分块定义验证。第六条是最重要的性质。设 $A=(a_{ij})$，$C=(c_{jk})$。则 $(A\otimes B)(C\otimes D)$ 的第 $(i,k)$ 个分块为
$$
\sum_j (a_{ij}B)(c_{jk}D)=\sum_j a_{ij}c_{jk}BD=(AC)_{ik}BD,
$$
这正是 $(AC)\otimes(BD)$ 的第 $(i,k)$ 个分块。

证明完毕.

**推论**:

若 $A,B$ 均可逆，则
$$
(A\otimes B)^{-1}=A^{-1}\otimes B^{-1}.
$$

**证明**:

由乘法性质，
$$
(A\otimes B)(A^{-1}\otimes B^{-1})=I\otimes I=I.
$$
右逆同理。

证明完毕.

**定理（Kronecker 积的秩、迹与行列式）**:

设 $A\in\mathbb F^{m\times n}$，$B\in\mathbb F^{p\times q}$，则
$$
\operatorname{rank}(A\otimes B)=\operatorname{rank} A\cdot\operatorname{rank} B.
$$
若 $A\in\mathbb F^{m\times m}$，$B\in\mathbb F^{p\times p}$，则
$$
\operatorname{tr}(A\otimes B)=\operatorname{tr}(A)\operatorname{tr}(B),
$$
且
$$
\det(A\otimes B)=\det(A)^p\det(B)^m.
$$

**证明**:

秩公式可用等价标准形证明。存在可逆矩阵 $P,Q,R,S$ 使
$$
PAQ=\begin{bmatrix}I_r&0\\0&0\end{bmatrix},
\qquad
RBS=\begin{bmatrix}I_s&0\\0&0\end{bmatrix},
$$
其中 $r=\operatorname{rank} A,s=\operatorname{rank} B$。由 Kronecker 积乘法性质，
$$
(P\otimes R)(A\otimes B)(Q\otimes S)=(PAQ)\otimes(RBS).
$$
左右乘可逆矩阵不改变秩，而右端显然有 $rs$ 个主对角单位元，所以秩为 $rs$。

迹公式由分块主对角线求和可得：$A\otimes B$ 的主对角分块为 $a_{11}B,
\ldots,a_{mm}B$，其迹之和为 $\sum_i a_{ii}\operatorname{tr} B=\operatorname{tr} A\operatorname{tr} B$。

行列式公式可先在 $A,B$ 可对角化时由特征值乘积得到，再由多项式恒等原理推广到一般矩阵；也可用 Schur 分解证明。若 $A$ 的特征值为 $\lambda_1,
\ldots,\lambda_m$，$B$ 的特征值为 $\mu_1,
\ldots,\mu_p$，则 $A\otimes B$ 的特征值为所有 $\lambda_i\mu_j$，其乘积为
$$
\prod_{i=1}^m\prod_{j=1}^p\lambda_i\mu_j
=\left(\prod_i\lambda_i\right)^p
\left(\prod_j\mu_j\right)^m.
$$

证明完毕.

**定理（特征值）**:

若 $Ax=\lambda x$，$By=\mu y$，其中 $x\ne0,y\ne0$，则
$$
(A\otimes B)(x\otimes y)=\lambda\mu(x\otimes y).
$$
因此，$A\otimes B$ 的特征值由 $A$ 与 $B$ 的特征值两两相乘得到，并按代数重数计数。

**证明**:

由 Kronecker 积的乘法性质，
$$
(A\otimes B)(x\otimes y)=Ax\otimes By
=\lambda x\otimes \mu y
=\lambda\mu(x\otimes y).
$$
若采用 Schur 分解或 Jordan 标准型，可得到完整的重数结论。

证明完毕.

#### 向量化公式

对矩阵 $X=(x_{ij})\in\mathbb F^{m\times n}$，定义 $\operatorname{vec}(X)$ 为把 $X$ 的列依次堆叠得到的 $mn$ 维列向量：
$$
\operatorname{vec}(X)=(x_{11},\ldots,x_{m1},x_{12},\ldots,x_{m2},\ldots,x_{1n},\ldots,x_{mn})^{\mathsf T}.
$$

**定理（vec 公式）**:

若矩阵阶数相容，则
$$
\operatorname{vec}(AXB)=(B^{\mathsf T}\otimes A)\operatorname{vec}(X).
$$
特别地，
$$
\operatorname{vec}(AX)=(I\otimes A)\operatorname{vec}(X),
\qquad
\operatorname{vec}(XB)=(B^{\mathsf T}\otimes I)\operatorname{vec}(X).
$$

**证明**:

设 $X=[x_1\ \cdots\ x_n]$ 按列分块。左乘 $A$ 时，$AX=[Ax_1\ \cdots\ Ax_n]$，所以向量化等于块对角矩阵 $I\otimes A$ 作用于 $\operatorname{vec}(X)$。右乘 $B$ 时，$XB$ 的第 $j$ 列为 $\sum_i b_{ij}x_i$，这正对应 $B^{\mathsf T}\otimes I$ 对列堆叠向量的作用。两式合成即
$$
\operatorname{vec}(AXB)=\operatorname{vec}((AX)B)=(B^{\mathsf T}\otimes I)\operatorname{vec}(AX)=(B^{\mathsf T}\otimes I)(I\otimes A)\operatorname{vec}(X)=(B^{\mathsf T}\otimes A)\operatorname{vec}(X).
$$

证明完毕.

**例（矩阵方程）**:

矩阵方程
$$
AXB=C
$$
可化为普通线性方程
$$
(B^{\mathsf T}\otimes A)\operatorname{vec}(X)=\operatorname{vec}(C).
$$
这说明 Kronecker 积是处理矩阵未知量方程的标准工具。

## 方阵的 Jordan 标准型与零化多项式

方阵相似分类的核心问题是：给定一个线性变换，能不能选择一组好的基，使得它的矩阵尽可能简单？对角化是最理想的情况；当矩阵不能对角化时，Jordan 标准型给出次优但仍十分清晰的结构。Jordan 标准型不仅解释了特征值和特征向量的不足，还把矩阵函数、矩阵幂、微分方程组、零化多项式等问题统一起来。

本章默认数域为 $\mathbb C$，除非特别说明。因为复数域上每个非常数多项式都能分解为一次因子，从而每个复方阵都有 Jordan 标准型。实矩阵也可以在复数域上讨论其 Jordan 标准型；若坚持实数域，则会出现 $2\times2$ 实块，对理论表述不如复数域简洁。

### 特征值、特征向量与对角化

**定义（特征值与特征向量）**:

设 $A\in\mathbb C^{n\times n}$。若存在 $0\ne x\in\mathbb C^n$ 与 $\lambda\in\mathbb C$，使得
$$
Ax=\lambda x,
$$
则称 $\lambda$ 是 $A$ 的特征值，$x$ 是对应于 $\lambda$ 的特征向量。

特征方程为
$$
\det(\lambda I-A)=0.
$$
其左边
$$
\chi_A(\lambda)=\det(\lambda I-A)
$$
称为 $A$ 的特征多项式。

**定义（代数重数与几何重数）**:

设 $\lambda$ 是 $A$ 的特征值。$\lambda$ 作为特征多项式 $\chi_A(t)$ 根的重数称为代数重数，记为 $m_a(\lambda)$。特征子空间
$$
E_\lambda=\operatorname{Ker}(A-\lambda I)
$$
的维数称为几何重数，记为
$$
m_g(\lambda)=\dim\operatorname{Ker}(A-\lambda I).
$$

**命题**:

对任意特征值 $\lambda$，有
$$
1\le m_g(\lambda)\le m_a(\lambda).
$$

**证明**:

左边来自特征值定义。证明右边可取 $E_\lambda$ 的一组基并扩充为 $\mathbb C^n$ 的基。在这组基下，$A$ 的矩阵具有分块形式
$$
\begin{bmatrix}
\lambda I_{m_g}& *\\
0&B
\end{bmatrix}.
$$
因此
$$
\chi_A(t)=\det(tI-A)=(t-\lambda)^{m_g}\det(tI-B),
$$
故 $m_a(\lambda)\ge m_g(\lambda)$。

证明完毕.

**定理（对角化判别）**:

矩阵 $A\in\mathbb C^{n\times n}$ 可对角化，当且仅当 $\mathbb C^n$ 有一组由 $A$ 的特征向量构成的基。等价地，对每个特征值 $\lambda$，
$$
m_g(\lambda)=m_a(\lambda),
$$
且所有特征值代数重数之和为 $n$。

**证明**:

若 $A=PDP^{-1}$，其中 $D$ 为对角矩阵，则 $P$ 的列向量就是 $A$ 的特征向量，且构成基。反之，若有特征向量基 $v_1,
\ldots,v_n$，令 $P=[v_1\ \cdots\ v_n]$，则
$$
AP=P\operatorname{diag}(\lambda_1,
\ldots,\lambda_n),
$$
所以 $P^{-1}AP$ 为对角矩阵。

不同特征值对应的特征子空间线性无关，因此所有特征向量能构成基，当且仅当各特征子空间维数之和为 $n$。而特征多项式在 $\mathbb C$ 上完全分解，所有代数重数之和为 $n$，故等价于每个 $m_g(\lambda)=m_a(\lambda)$。

证明完毕.

### 广义特征向量与 Jordan 链

当 $m_g(\lambda)<m_a(\lambda)$ 时，普通特征向量不够多。解决办法是引入广义特征向量。

**定义（广义特征向量）**:

设 $\lambda$ 是 $A$ 的特征值。若非零向量 $x$ 满足
$$
(A-\lambda I)^k x=0
$$
对某个正整数 $k$ 成立，则称 $x$ 是 $A$ 关于 $\lambda$ 的广义特征向量。

对应的广义特征子空间定义为
$$
G_\lambda=\bigcup_{k\ge1}\operatorname{Ker}(A-\lambda I)^k.
$$
在有限维空间中，当 $k$ 足够大时该升链稳定：
$$
\operatorname{Ker}(A-\lambda I)
\subseteq \operatorname{Ker}(A-\lambda I)^2
\subseteq\cdots.
$$
若 $m_a(\lambda)=m$，则
$$
G_\lambda=\operatorname{Ker}(A-\lambda I)^m.
$$

**定义（Jordan 块）**:

大小为 $s$、特征值为 $\lambda$ 的 Jordan 块定义为
$$
J_s(\lambda)=
\begin{bmatrix}
\lambda&1&0&\cdots&0\\
0&\lambda&1&\cdots&0\\
\vdots&\vdots&\ddots&\ddots&\vdots\\
0&0&\cdots&\lambda&1\\
0&0&\cdots&0&\lambda
\end{bmatrix}
=\lambda I_s+N_s,
$$
其中 $N_s$ 是上超对角线为 $1$、其余为 $0$ 的幂零矩阵。

**定义（Jordan 链）**:

向量组
$$
v_1,v_2,
\ldots,v_s
$$
称为关于特征值 $\lambda$ 的一条 Jordan 链，如果
$$
(A-\lambda I)v_1=0,
$$
且
$$
(A-\lambda I)v_2=v_1,
\quad
(A-\lambda I)v_3=v_2,
\quad \ldots,\quad
(A-\lambda I)v_s=v_{s-1}.
$$
此时在有序基 $v_1,
\ldots,v_s$ 下，$A$ 在该链张成空间上的矩阵为
$$
\begin{bmatrix}
\lambda&1&0&\cdots&0\\
0&\lambda&1&\cdots&0\\
\vdots&\vdots&\ddots&\ddots&\vdots\\
0&0&\cdots&\lambda&1\\
0&0&\cdots&0&\lambda
\end{bmatrix}.
$$

**证明（说明）**:

由链条件，
$$
Av_1=\lambda v_1,
\qquad
Av_j=v_{j-1}+\lambda v_j\quad (j\ge2).
$$
因此第 $j$ 列坐标正好为第 $j$ 个对角元 $\lambda$ 与第 $j-1$ 行的 $1$。

证明完毕.

### 幂零矩阵的结构

Jordan 标准型的关键先理解幂零矩阵。矩阵 $N$ 称为幂零，如果存在 $k$ 使 $N^k=0$。Jordan 块 $J_s(\lambda)$ 减去 $\lambda I$ 后就是幂零块 $N_s$。

**定义（幂零指数）**:

若 $N^r=0$ 且 $N^{r-1}\ne0$，则称 $r$ 为幂零矩阵 $N$ 的幂零指数。

**引理**:

若 $N$ 是幂零矩阵，则所有特征值均为 $0$。

**证明**:

若 $Nx=\lambda x$，$x\ne0$，则 $N^kx=\lambda^kx$。当 $N^k=0$ 时得到 $\lambda^kx=0$，故 $\lambda=0$。

证明完毕.

对于一个幂零 Jordan 块 $N_s$，容易计算
$$
\dim\operatorname{Ker} N_s^k=\min(k,s).
$$
若幂零矩阵由若干幂零 Jordan 块组成，块大小为 $s_1,
\ldots,s_r$，则
$$
\dim\operatorname{Ker} N^k=\sum_{i=1}^r \min(k,s_i).
$$
这给出了由核维数恢复 Jordan 块大小的公式。

**定理（核维数与 Jordan 块数）**:

设 $A$ 关于特征值 $\lambda$ 的 Jordan 块大小为 $s_1,
\ldots,s_r$。令
$$
d_k=\dim\operatorname{Ker}(A-\lambda I)^k,
\qquad d_0=0.
$$
则：
1. $d_k-d_{k-1}$ 等于大小至少为 $k$ 的 Jordan 块个数；
2. 大小恰为 $k$ 的 Jordan 块个数为
$$
  (d_k-d_{k-1})-(d_{k+1}-d_k)=2d_k-d_{k-1}-d_{k+1}.
$$

**证明**:

对一个大小为 $s$ 的块，$\dim\operatorname{Ker} N_s^k=\min(k,s)$，于是
$$
\min(k,s)-\min(k-1,s)=
\begin{cases}
1,&s\ge k,\\
0,&s<k.
\end{cases}
$$
对所有 Jordan 块求和得到第一条。第二条是“至少为 $k$”的块数减去“至少为 $k+1$”的块数。

证明完毕.

### Jordan 标准型

**定理（Jordan 标准型）**:

任意复方阵 $A\in\mathbb C^{n\times n}$ 都相似于一个分块对角矩阵
$$
J=\operatorname{diag}(J_{s_{11}}(\lambda_1),
\ldots,J_{s_{1r_1}}(\lambda_1),
\ldots,
J_{s_{p1}}(\lambda_p),
\ldots,J_{s_{pr_p}}(\lambda_p)),
$$
其中 $\lambda_1,
\ldots,\lambda_p$ 是 $A$ 的不同特征值。该形式除 Jordan 块排列次序外唯一。

**证明（证明思路）**:

完整证明可分为三步。

第一步，利用特征多项式的分解和 Bezout 恒等式证明广义特征子空间分解：
$$
\mathbb C^n=G_{\lambda_1}\oplus\cdots\oplus G_{\lambda_p},
$$
其中
$$
G_{\lambda_i}=\operatorname{Ker}(A-\lambda_i I)^{m_i},
$$
$m_i$ 是 $\lambda_i$ 的代数重数。该分解称为主分解或初等因子分解。

第二步，在每个 $G_\lambda$ 上研究
$$
N=A-\lambda I.
$$
此时 $N$ 是幂零变换。通过逐层考察
$$
\operatorname{Ker} N\subseteq\operatorname{Ker} N^2\subseteq\cdots
$$
可以选择若干条 Jordan 链，使其合起来构成 $G_\lambda$ 的一组基。

第三步，把所有特征值对应的 Jordan 链基合并，就得到整个空间的一组基。在这组基下，$A$ 的矩阵就是 Jordan 块的分块对角阵。

唯一性来自上一节的核维数公式。因为 $\dim\operatorname{Ker}(A-\lambda I)^k$ 是相似不变量，所以每个特征值下各种大小的 Jordan 块个数都被唯一确定。

证明完毕.

**注**:

Jordan 标准型不是数值计算中稳定的工具，因为矩阵轻微扰动可能显著改变 Jordan 块结构。但在理论分析中，它极其重要，因为它提供了方阵相似分类的精确答案。

#### Jordan 标准型的直接意义

若
$$
A=PJP^{-1},
$$
则许多关于 $A$ 的问题转化为关于 Jordan 块的问题。例如：
$$
A^k=PJ^kP^{-1},
\qquad
f(A)=Pf(J)P^{-1},
\qquad
\exp(A)=P\exp(J)P^{-1}.
$$
对单个块
$$
J_s(\lambda)=\lambda I+N,
\quad N^s=0,
$$
由二项式公式
$$
J_s(\lambda)^k=(\lambda I+N)^k
=\sum_{r=0}^{s-1}\binom{k}{r}\lambda^{k-r}N^r.
$$
若 $f$ 在 $\lambda$ 附近可展开，则
$$
f(J_s(\lambda))=
\sum_{r=0}^{s-1}\frac{f^{(r)}(\lambda)}{r!}N^r.
$$
这解释了为什么非对角 Jordan 块会引入导数项。

### 零化多项式与最小多项式

**定义（零化多项式）**:

设 $A\in\mathbb F^{n\times n}$，若多项式 $p(t)\in\mathbb F[t]$ 满足
$$
p(A)=0,
$$
则称 $p$ 是 $A$ 的零化多项式。

零化多项式反映了矩阵满足的代数方程。最基本的零化多项式来自 Cayley--Hamilton 定理。

**定理（Cayley--Hamilton 定理）**:

设
$$
\chi_A(t)=\det(tI-A)
$$
是 $A$ 的特征多项式，则
$$
\chi_A(A)=0.
$$

**证明**:

记 $\operatorname{adj}(tI-A)$ 为 $tI-A$ 的伴随矩阵。由伴随矩阵恒等式
$$
(tI-A)\operatorname{adj}(tI-A)=\det(tI-A)I=\chi_A(t)I.
$$
注意这是一条关于 $t$ 的矩阵多项式恒等式。由于 $A$ 与 $A$ 的多项式可交换，可以把 $t$ 替换为矩阵 $A$，得到
$$
(AI-A)\operatorname{adj}(AI-A)=0=\chi_A(A)I.
$$
故 $\chi_A(A)=0$。

证明完毕.

**定义（最小多项式）**:

所有使 $p(A)=0$ 的非零多项式中，存在唯一的首一且次数最低的多项式，称为 $A$ 的最小多项式，记为 $m_A(t)$。

**定理（最小多项式的整除性质）**:

设 $p(t)$ 是任意多项式，则
$$
p(A)=0
$$
当且仅当
$$
m_A(t)\mid p(t).
$$
特别地，
$$
m_A(t)\mid \chi_A(t).
$$

**证明**:

由多项式除法，存在 $q(t),r(t)$，使
$$
p(t)=q(t)m_A(t)+r(t),
\qquad \deg r<\deg m_A.
$$
代入 $A$，若 $p(A)=0$，则
$$
r(A)=p(A)-q(A)m_A(A)=0.
$$
由 $m_A$ 的最低次数性，只能有 $r=0$，故 $m_A\mid p$。反向显然。

Cayley--Hamilton 定理给出 $\chi_A(A)=0$，于是 $m_A\mid\chi_A$。

证明完毕.

**定理（最小多项式与 Jordan 块）**:

设 $A$ 的不同特征值为 $\lambda_1,
\ldots,\lambda_p$。对每个 $\lambda_i$，令 $s_i$ 为 $A$ 关于 $\lambda_i$ 的最大 Jordan 块大小。则
$$
m_A(t)=\prod_{i=1}^p (t-\lambda_i)^{s_i}.
$$

**证明**:

相似不改变最小多项式，因此可设 $A=J$ 为 Jordan 标准型。对单个块
$$
J_s(\lambda)-\lambda I=N_s,
$$
且 $N_s^s=0$，$N_s^{s-1}\ne0$。所以要使该块被零化，$p(t)$ 必须含有因子 $(t-\lambda)^s$。对同一特征值的多个块，只需要最高阶那个块决定指数。对不同特征值，由互素因子相乘得到整体最小多项式。

证明完毕.

**推论（对角化判别）**:

矩阵 $A$ 可对角化，当且仅当其最小多项式没有重根，即
$$
m_A(t)=\prod_{i=1}^p(t-\lambda_i).
$$

**证明**:

$A$ 可对角化当且仅当所有 Jordan 块大小均为 $1$。由最小多项式与最大 Jordan 块大小的关系，这等价于每个因子的指数均为 $1$。

证明完毕.

#### 零化多项式的用法

若知道 $A$ 满足一个低次数方程，就可以化简高次幂。例如
$$
A^2-3A+2I=0
$$
意味着
$$
A^2=3A-2I.
$$
于是 $A^k$ 都可以写成 $A$ 与 $I$ 的线性组合。一般地，若 $m_A$ 的次数为 $r$，则所有矩阵多项式 $p(A)$ 都可以化成
$$
c_0I+c_1A+\cdots+c_{r-1}A^{r-1}.
$$

### Jordan 标准型的求法

实际求 Jordan 标准型时，核心不是先找所有广义特征向量，而是先确定每个特征值的 Jordan 块大小。

#### 步骤一：求特征值与代数重数

计算
$$
\chi_A(t)=\det(tI-A).
$$
将其分解为
$$
\chi_A(t)=\prod_{i=1}^p(t-\lambda_i)^{m_i}.
$$
其中 $m_i$ 是 $\lambda_i$ 的代数重数。所有关于 $\lambda_i$ 的 Jordan 块大小之和为 $m_i$。

#### 步骤二：用核维数确定块数

对每个特征值 $\lambda$，计算
$$
d_k=\dim\operatorname{Ker}(A-\lambda I)^k,
\qquad k=1,2,
\ldots,m_a(\lambda).
$$
则
$$
d_1
$$
就是该特征值对应的 Jordan 块总数，因为每个 Jordan 块贡献一个特征向量。

大小至少为 $k$ 的 Jordan 块个数为
$$
d_k-d_{k-1}.
$$
大小恰为 $k$ 的 Jordan 块个数为
$$
2d_k-d_{k-1}-d_{k+1}.
$$
通常计算到 $d_k=m_a(\lambda)$ 就可以停止。

#### 步骤三：构造 Jordan 链

确定块大小后，需要找相似变换矩阵 $P$ 时，才需要构造基。对每个块大小 $s$，找向量 $v_s$ 满足
$$
(A-\lambda I)^s v_s=0,
\qquad
(A-\lambda I)^{s-1}v_s\ne0.
$$
然后定义
$$
v_{s-1}=(A-\lambda I)v_s,
\quad
v_{s-2}=(A-\lambda I)^2v_s,
\quad \ldots,
\quad
v_1=(A-\lambda I)^{s-1}v_s.
$$
这样得到一条 Jordan 链
$$
v_1,
\ldots,v_s.
$$
把所有 Jordan 链合并为矩阵 $P$ 的列，就有
$$
P^{-1}AP=J.
$$

**例**:

设
$$
A=
\begin{bmatrix}
2&1&0\\
0&2&0\\
0&0&2
\end{bmatrix}.
$$
唯一特征值为 $2$，代数重数为 $3$。计算
$$
A-2I=
\begin{bmatrix}
0&1&0\\
0&0&0\\
0&0&0
\end{bmatrix}.
$$
于是
$$
d_1=\dim\operatorname{Ker}(A-2I)=2,
\qquad
(A-2I)^2=0,
\qquad d_2=3.
$$
大小至少为 $1$ 的块有 $d_1-d_0=2$ 个，大小至少为 $2$ 的块有 $d_2-d_1=1$ 个。因此 Jordan 块大小为 $2$ 和 $1$，即
$$
J=\operatorname{diag}(J_2(2),J_1(2)).
$$
事实上 $A$ 已经是该 Jordan 标准型。

**例**:

设 $A$ 的特征多项式为
$$
\chi_A(t)=(t-1)^5(t+2)^3.
$$
若
$$
\dim\operatorname{Ker}(A-I)=2,
\quad
\dim\operatorname{Ker}(A-I)^2=4,
\quad
\dim\operatorname{Ker}(A-I)^3=5,
$$
则关于 $\lambda=1$ 的块数如下：大小至少 $1$ 的块有 $2$ 个，大小至少 $2$ 的块有 $4-2=2$ 个，大小至少 $3$ 的块有 $5-4=1$ 个。因此块大小为 $3$ 和 $2$。

若关于 $\lambda=-2$，有
$$
\dim\operatorname{Ker}(A+2I)=1,
\quad
\dim\operatorname{Ker}(A+2I)^2=2,
\quad
\dim\operatorname{Ker}(A+2I)^3=3,
$$
则只有一个大小为 $3$ 的块。于是 Jordan 型为
$$
\operatorname{diag}(J_3(1),J_2(1),J_3(-2)).
$$
最小多项式为
$$
m_A(t)=(t-1)^3(t+2)^3.
$$

## 矩阵分解及其应用

矩阵分解的思想是把复杂矩阵拆成若干结构更简单的矩阵。不同分解服务于不同目标：满秩分解揭示秩结构，LU 分解服务线性方程组求解，QR 分解服务正交化与最小二乘，谱分解揭示自伴随算子的本征结构，Schur 分解提供一般复矩阵的稳定三角化，奇异值分解则统一处理任意矩阵的几何伸缩、低秩逼近和伪逆。

### 矩阵的满秩分解

**定义（满秩分解）**:

设 $A\in\mathbb F^{m\times n}$，$\operatorname{rank} A=r$。若存在 $F\in\mathbb F^{m\times r}$ 与 $G\in\mathbb F^{r\times n}$，使得
$$
A=FG,
$$
且
$$
\operatorname{rank} F=\operatorname{rank} G=r,
$$
则称 $A=FG$ 是 $A$ 的满秩分解。

**定理（满秩分解存在性）**:

任意秩为 $r$ 的矩阵 $A\in\mathbb F^{m\times n}$ 都存在满秩分解。

**证明**:

取 $A$ 的列空间的一组基，由这些列组成矩阵 $F\in\mathbb F^{m\times r}$。由于 $A$ 的每一列都属于 $\operatorname{Im} A=\operatorname{Im} F$，故每一列都可唯一写成 $F$ 的列向量线性组合。把这些坐标作为矩阵 $G$ 的列，就有 $A=FG$。显然 $F$ 列满秩，即 $\operatorname{rank} F=r$。又因为 $\operatorname{rank} A\le \operatorname{rank} G\le r$ 且 $\operatorname{rank} A=r$，所以 $\operatorname{rank} G=r$。

证明完毕.

#### 用初等变换求满秩分解

若通过行初等变换把 $A$ 化为行阶梯形矩阵 $R$：
$$
EA=R,
$$
其中 $E$ 可逆。设 $R$ 的非零行组成 $G\in\mathbb F^{r\times n}$。若在原矩阵 $A$ 中取与主元列对应的 $r$ 列组成 $F$，则通常可得到
$$
A=FG.
$$
更直接的方法是选取 $A$ 的列空间基作为 $F$，再解线性方程求每列坐标。

**例**:

设
$$
A=\begin{bmatrix}
1&2&3\\
2&4&6\\
1&1&1
\end{bmatrix}.
$$
记三列分别为 $a_1,a_2,a_3$。可以验证 $a_1,a_2$ 线性无关，且
$$
a_3=-a_1+2a_2.
$$
所以 $r=2$。取列空间基组成
$$
F=\begin{bmatrix}
1&2\\
2&4\\
1&1
\end{bmatrix},
$$
再把每一列在 $a_1,a_2$ 下的坐标写成 $G$ 的对应列，即
$$
G=\begin{bmatrix}
1&0&-1\\
0&1&2
\end{bmatrix}.
$$
于是
$$
FG=\begin{bmatrix}
1&2&3\\
2&4&6\\
1&1&1
\end{bmatrix}=A,
$$
这就是一个满秩分解。

#### 满秩分解的意义

满秩分解把 $A$ 分成
$$
\mathbb F^n \xrightarrow{G} \mathbb F^r \xrightarrow{F} \mathbb F^m.
$$
其中 $G$ 将输入压缩到秩为 $r$ 的有效坐标，$F$ 再把有效坐标嵌入到列空间中。因此满秩分解清楚地揭示：秩为 $r$ 的线性变换本质上只携带 $r$ 维信息。

### LU 分解

**定义（LU 分解）**:

设 $A\in\mathbb F^{n\times n}$。若
$$
A=LU,
$$
其中 $L$ 为下三角矩阵，$U$ 为上三角矩阵，则称其为 LU 分解。若 $L$ 的对角元全为 $1$，称为 Doolittle 型 LU 分解。

LU 分解来自 Gaussian 消元。它的意义是把求解
$$
Ax=b
$$
拆成两步：
$$
Ly=b,
\qquad
Ux=y.
$$
前者前代，后者回代，计算量远小于直接求逆。

**定理（无需换行的 LU 分解判别）**:

若 $A\in\mathbb F^{n\times n}$ 的所有顺序主子式均非零：
$$
\det A_k\ne0,
\qquad k=1,
\ldots,n,
$$
其中 $A_k$ 是 $A$ 左上角 $k\times k$ 子矩阵，则 $A$ 存在唯一的分解
$$
A=LU,
$$
其中 $L$ 为单位下三角矩阵，$U$ 为上三角矩阵。

**证明**:

Gaussian 消元第 $k$ 步需要第 $k$ 个主元非零。顺序主子式非零保证每一步主元非零，因此不需要行交换即可完成消元。把消元乘子放入 $L$ 的相应位置，把消元后的上三角矩阵记为 $U$，得 $A=LU$。

唯一性：若
$$
A=L_1U_1=L_2U_2,
$$
则
$$
L_2^{-1}L_1=U_2U_1^{-1}.
$$
左边是单位下三角矩阵，右边是上三角矩阵。一个矩阵若既是单位下三角又是上三角，只能是单位矩阵。因此 $L_1=L_2$，$U_1=U_2$。

证明完毕.

**注**:

一般矩阵未必能不换行地做 LU 分解。实际计算中常用带置换的分解
$$
PA=LU,
$$
其中 $P$ 是置换矩阵。这称为部分主元 LU 分解，具有更好的数值稳定性。

**例**:

设
$$
A=\begin{bmatrix}
2&1&1\\
4&-6&0\\
-2&7&2
\end{bmatrix}.
$$
消去第一列：乘子为 $2,-1$。消去后第二、三行变为
$$
(0,-8,-2),
\quad
(0,8,3).
$$
再消去第二列，乘子为 $-1$，得到
$$
U=\begin{bmatrix}
2&1&1\\
0&-8&-2\\
0&0&1
\end{bmatrix},
\qquad
L=\begin{bmatrix}
1&0&0\\
2&1&0\\
-1&-1&1
\end{bmatrix}.
$$
可验证 $A=LU$。

### QR 分解

**定义（QR 分解）**:

设 $A\in\mathbb R^{m\times n}$，$m\ge n$。若
$$
A=QR,
$$
其中 $Q\in\mathbb R^{m\times n}$ 的列向量标准正交，即 $Q^TQ=I_n$，$R\in\mathbb R^{n\times n}$ 为上三角矩阵，则称为约化 QR 分解。若 $Q\in\mathbb R^{m\times m}$ 是正交矩阵，$R\in\mathbb R^{m\times n}$ 为上梯形矩阵，则称为完全 QR 分解。

**定理（满列秩矩阵的 QR 分解）**:

若 $A=[a_1\ \cdots\ a_n]\in\mathbb R^{m\times n}$ 满列秩，则存在唯一分解
$$
A=QR,
$$
其中 $Q^TQ=I$，$R$ 为对角元为正的上三角矩阵。

**证明**:

对 $a_1,
\ldots,a_n$ 进行 Gram--Schmidt 正交化。令
$$
u_1=a_1,
\qquad
q_1=\frac{u_1}{\left\lVert u_1\right\rVert}.
$$
对 $k\ge2$，令
$$
u_k=a_k-\sum_{j=1}^{k-1}(q_j^Ta_k)q_j,
\qquad
q_k=\frac{u_k}{\left\lVert u_k\right\rVert}.
$$
由于 $A$ 满列秩，每一步 $u_k\ne0$。于是
$$
a_k=\sum_{j=1}^{k-1}(q_j^Ta_k)q_j+\left\lVert u_k\right\rVertq_k.
$$
令
$$
r_{jk}=q_j^Ta_k\quad(j<k),
\qquad
r_{kk}=\left\lVert u_k\right\rVert>0,
\qquad
r_{jk}=0\quad(j>k),
$$
得到 $A=QR$。

唯一性：若 $A=Q_1R_1=Q_2R_2$，则
$$
Q_2^TQ_1=R_2R_1^{-1}.
$$
左边是正交矩阵，右边是上三角矩阵且对角元为正。一个上三角正交矩阵若对角元为正，只能是 $I$。故 $Q_1=Q_2,R_1=R_2$。

证明完毕.

#### QR 分解与最小二乘

考虑超定方程
$$
Ax\approx b,
\qquad A\in\mathbb R^{m\times n},
\ m\ge n,
$$
最小二乘问题为
$$
\min_x\left\lVert Ax-b\right\rVert_2.
$$
若 $A=QR$ 为约化 QR 分解，则
$$
\left\lVert Ax-b\right\rVert_2^2
=\left\lVert QRx-b\right\rVert_2^2.
$$
将 $Q$ 扩充为正交矩阵 $[Q\ Q_\perp]$，由正交变换不改变二范数，
$$
\left\lVert QRx-b\right\rVert_2^2
=\left\lVert Rx-Q^Tb\right\rVert_2^2+
orm{Q_\perp^Tb}_2^2.
$$
第二项与 $x$ 无关，所以最优解满足
$$
Rx=Q^Tb.
$$
这比正规方程
$$
A^TAx=A^Tb
$$
数值上更稳定，因为正规方程会把条件数平方。

### 谱分解

**定义（Hermite 矩阵与正规矩阵）**:

复矩阵 $A$ 若满足
$$
A^*=A,
$$
称为 Hermite 矩阵；实情形即对称矩阵 $A^T=A$。若
$$
A^*A=AA^*,
$$
称为正规矩阵。

**定理（Hermite 矩阵的谱定理）**:

若 $A=A^*\in\mathbb C^{n\times n}$，则存在酉矩阵 $U$ 与实对角矩阵 $\Lambda$，使得
$$
A=U\Lambda U^*.
$$
也就是说，Hermite 矩阵存在一组标准正交特征向量，且所有特征值为实数。

**证明**:

先证特征值为实数。若 $Ax=\lambda x$，$x\ne0$，则
$$
\lambda\left\langle x,x\right\rangle=\left\langle Ax,x\right\rangle=\left\langle x,A^*x\right\rangle=\left\langle x,Ax\right\rangle=\overline{\lambda}\left\langle x,x\right\rangle,
$$
故 $\lambda=\overline{\lambda}$。

再取一个单位特征向量 $u_1$。其正交补 $u_1^\perp$ 在 $A$ 下不变：若 $y\perp u_1$，则
$$
\left\langle Ay,u_1\right\rangle=\left\langle y,Au_1\right\rangle=\left\langle y,\lambda u_1\right\rangle=0.
$$
因此可在 $u_1^\perp$ 上归纳，得到标准正交特征向量基。

证明完毕.

**推论（实对称矩阵谱分解）**:

若 $A=A^T\in\mathbb R^{n\times n}$，则存在正交矩阵 $Q$ 与实对角矩阵 $\Lambda$，使得
$$
A=Q\Lambda Q^T.
$$

**定理（正规矩阵的谱定理）**:

复矩阵 $A$ 酉对角化，当且仅当 $A$ 正规。即存在酉矩阵 $U$ 使
$$
A=U\Lambda U^*,
$$
当且仅当
$$
A^*A=AA^*.
$$

**证明**:

若 $A=U\Lambda U^*$，则
$$
A^*A=U\Lambda^*\Lambda U^*,
\qquad
AA^*=U\Lambda\Lambda^*U^*,
$$
二者相等，因为对角矩阵与其共轭转置可交换。

反过来，可用 Schur 分解。任意复矩阵都酉相似于上三角矩阵：
$$
A=UTU^*.
$$
若 $A$ 正规，则 $T$ 正规。上三角正规矩阵必为对角矩阵。证明如下：比较 $T^*T$ 与 $TT^*$ 的第一个对角元，可知第一行非对角元为零；再对剩余主子矩阵归纳，得 $T$ 对角。因此 $A$ 酉对角化。

证明完毕.

#### 谱分解的几何意义

若
$$
A=Q\Lambda Q^T
$$
为实对称矩阵的谱分解，则
$$
Ax=\sum_{i=1}^n\lambda_i(q_i^Tx)q_i.
$$
即 $A$ 在每个正交特征方向 $q_i$ 上只做伸缩，伸缩倍数为 $\lambda_i$。因此二次型
$$
x^TAx
$$
可化为
$$
x^TAx=\sum_{i=1}^n\lambda_i y_i^2,
\qquad y=Q^Tx.
$$
由此立即得到正定判别：
$$
A\text{ 正定}\quad\Longleftrightarrow\quad \lambda_i>0\ \text{对所有 }i.
$$

### Schur 三角分解

**定理（Schur 分解）**:

任意复矩阵 $A\in\mathbb C^{n\times n}$ 都存在酉矩阵 $U$，使得
$$
U^*AU=T
$$
为上三角矩阵。$T$ 的对角元正是 $A$ 的特征值。

**证明**:

对 $n$ 归纳。$n=1$ 显然。设 $\lambda$ 是 $A$ 的一个特征值，取单位特征向量 $u_1$。将 $u_1$ 扩充为 $\mathbb C^n$ 的标准正交基，形成酉矩阵 $U_1=[u_1\ U_2]$。则
$$
U_1^*AU_1=
\begin{bmatrix}
\lambda&*\\
0&A_1
\end{bmatrix}.
$$
由归纳假设，存在酉矩阵 $V$ 使 $V^*A_1V$ 上三角。令
$$
U=U_1\begin{bmatrix}1&0\\0&V\end{bmatrix},
$$
即可得到上三角化。

证明完毕.

Schur 分解比 Jordan 分解更适合数值计算，因为酉变换保持二范数，不会显著放大误差。它也是 QR 算法求特征值的理论基础。

### 奇异值分解

**定义（奇异值）**:

设 $A\in\mathbb C^{m\times n}$。矩阵 $A^*A$ 是半正定 Hermite 矩阵，其特征值非负。设这些特征值为
$$
\sigma_1^2\ge\sigma_2^2\ge\cdots\ge\sigma_n^2\ge0.
$$
非负数 $\sigma_i$ 称为 $A$ 的奇异值。

**定理（奇异值分解）**:

任意矩阵 $A\in\mathbb C^{m\times n}$ 都存在酉矩阵 $U\in\mathbb C^{m\times m}$、$V\in\mathbb C^{n\times n}$，以及矩形对角矩阵
$$
\Sigma=\begin{bmatrix}
\operatorname{diag}(\sigma_1,
\ldots,\sigma_r)&0\\
0&0
\end{bmatrix}\in\mathbb R^{m\times n},
$$
使得
$$
A=U\Sigma V^*.
$$
其中 $\sigma_1\ge\cdots\ge\sigma_r>0$，$r=\operatorname{rank} A$。

**证明**:

对 Hermite 半正定矩阵 $A^*A$ 作谱分解。存在酉矩阵 $V=[v_1\ \cdots\ v_n]$，使得
$$
A^*Av_i=\sigma_i^2v_i.
$$
对 $\sigma_i>0$ 的部分，定义
$$
u_i=\frac{Av_i}{\sigma_i}.
$$
则
$$
\left\langle u_i,u_j\right\rangle
=\frac{1}{\sigma_i\sigma_j}\left\langle Av_i,Av_j\right\rangle
=\frac{1}{\sigma_i\sigma_j}\left\langle v_i,A^*Av_j\right\rangle
=\frac{\sigma_j^2}{\sigma_i\sigma_j}\left\langle v_i,v_j\right\rangle
=\delta_{ij}.
$$
所以 $u_1,
\ldots,u_r$ 标准正交。将其扩充为 $\mathbb C^m$ 的标准正交基，得到酉矩阵 $U$。由构造，
$$
Av_i=\sigma_i u_i\quad(i\le r),
\qquad
Av_i=0\quad(i>r),
$$
这正是矩阵等式 $AV=U\Sigma$，即 $A=U\Sigma V^*$。

证明完毕.

#### SVD 的几何意义

分解
$$
A=U\Sigma V^*
$$
表示线性变换 $A$ 可以分成三步：
1. $V^*$：对输入空间做正交坐标旋转；
2. $\Sigma$：沿坐标轴做伸缩，并可能降维；
3. $U$：对输出空间做正交坐标旋转。

所以奇异值就是矩阵沿某些正交方向的伸缩倍数。

**定理（四个基本子空间的 SVD 描述）**:

设
$$
A=U\Sigma V^*,
$$
其中正奇异值个数为 $r$，$U=[u_1\ \cdots\ u_m]$，$V=[v_1\ \cdots\ v_n]$。则
$$
\operatorname{Im} A=\operatorname{span}(u_1,
\ldots,u_r),
$$
$$
\operatorname{Ker} A=\operatorname{span}(v_{r+1},
\ldots,v_n),
$$
$$
\operatorname{Im} A^*=\operatorname{span}(v_1,
\ldots,v_r),
$$
$$
\operatorname{Ker} A^*=\operatorname{span}(u_{r+1},
\ldots,u_m).
$$

**证明**:

由
$$
Av_i=\sigma_i u_i\quad(i\le r),
\qquad
Av_i=0\quad(i>r)
$$
直接得到列空间与核空间的描述。对 $A^*=V\Sigma^*U^*$ 同理。

证明完毕.

#### Moore--Penrose 伪逆

**定义（伪逆）**:

若 $A=U\Sigma V^*$，定义
$$
A^+=V\Sigma^+U^*,
$$
其中 $\Sigma^+$ 把每个正奇异值 $\sigma_i$ 替换为 $1/\sigma_i$，并转置矩形对角位置。

伪逆给出最小二乘问题的最小范数解。对任意 $b$，
$$
x=A^+b
$$
是使 $\left\lVert Ax-b\right\rVert_2$ 最小的解中欧氏范数最小者。

**证明（证明）**:

把变量换成 $y=V^*x$，把右端换成 $c=U^*b$。由于 $U,V$ 酉，最小化
$$
\left\lVert Ax-b\right\rVert_2=\left\lVert \Sigma y-c\right\rVert_2.
$$
若正奇异值为 $\sigma_1,
\ldots,\sigma_r$，则目标函数为
$$
\sum_{i=1}^r(\sigma_i y_i-c_i)^2+
\sum_{i=r+1}^m c_i^2.
$$
最小化要求
$$
y_i=\frac{c_i}{\sigma_i}\quad(i\le r),
$$
而 $i>r$ 的分量不影响残差。为了使 $\left\lVert x\right\rVert_2=\left\lVert y\right\rVert_2$ 最小，取 $y_i=0$ 对 $i>r$。这正是 $x=V\Sigma^+U^*b=A^+b$。

证明完毕.

#### 低秩逼近

**定理（Eckart--Young 定理）**:

设
$$
A=\sum_{i=1}^r\sigma_i u_iv_i^*
$$
是 $A$ 的 SVD 展开，且
$$
\sigma_1\ge\cdots\ge\sigma_r>0.
$$
在所有秩不超过 $k$ 的矩阵 $B$ 中，
$$
A_k=\sum_{i=1}^k\sigma_i u_iv_i^*
$$
是 $A$ 在二范数和 Frobenius 范数下的最佳秩 $k$ 逼近，并且
$$
\left\lVert A-A_k\right\rVert_2=\sigma_{k+1},
$$
$$
\left\lVert A-A_k\right\rVert_F=\left(\sum_{i=k+1}^r\sigma_i^2\right)^{1/2}.
$$

**证明（证明思路）**:

由于酉变换不改变二范数和 Frobenius 范数，可以把问题化为逼近对角矩阵 $\Sigma$。保留最大的 $k$ 个对角元显然使剩余最大对角元最小，也使剩余平方和最小。严格证明可用奇异值的 Courant--Fischer 极值刻画。

证明完毕.

## 向量范数

范数是“长度”的抽象。在线性空间中，范数使得可以讨论收敛、误差、稳定性与最优化。有限维空间中所有范数诱导的拓扑相同，但不同范数对具体估计和计算有不同意义。

### 范数的定义

**定义（向量范数）**:

设 $V$ 是 $\mathbb F$ 上的线性空间。映射 $\left\lVert \cdot\right\rVert:V\to\mathbb R$ 称为范数，如果对任意 $x,y\in V$ 与任意 $a\in\mathbb F$ 满足：
1. 非负性：$\left\lVert x\right\rVert\ge0$；
2. 正定性：$\left\lVert x\right\rVert=0$ 当且仅当 $x=0$；
3. 齐次性：$\left\lVert ax\right\rVert=\left\lvert a\right\rvert\left\lVert x\right\rVert$；
4. 三角不等式：$\left\lVert x+y\right\rVert\le\left\lVert x\right\rVert+\left\lVert y\right\rVert$。

**例（常见 $p$ 范数）**:

对 $x=(x_1,
\ldots,x_n)^T\in\mathbb F^n$，定义
$$
\left\lVert x\right\rVert_p=\left(\sum_{i=1}^n\left\lvert x_i\right\rvert^p\right)^{1/p},
\qquad 1\le p<\infty.
$$
当 $p=2$ 时为 Euclidean 范数；当 $p=1$ 时为绝对值和范数。再定义
$$
\left\lVert x\right\rVert_\infty=\max_{1\le i\le n}\left\lvert x_i\right\rvert.
$$
这些都是范数。

$p$ 范数满足三角不等式的关键是 Minkowski 不等式。证明它通常需要 Holder 不等式。

**定理（Holder 不等式）**:

设 $1<p<\infty$，且 $q$ 满足
$$
\frac1p+\frac1q=1.
$$
则对任意 $x,y\in\mathbb F^n$，
$$
\sum_{i=1}^n\left\lvert x_i y_i\right\rvert
\le
\left(\sum_{i=1}^n\left\lvert x_i\right\rvert^p\right)^{1/p}
\left(\sum_{i=1}^n\left\lvert y_i\right\rvert^q\right)^{1/q}.
$$

**证明**:

若 $x=0$ 或 $y=0$，结论显然。令
$$
a_i=\frac{\left\lvert x_i\right\rvert}{\left\lVert x\right\rVert_p},
\qquad
b_i=\frac{\left\lvert y_i\right\rvert}{\left\lVert y\right\rVert_q}.
$$
则 $\sum a_i^p=1$，$\sum b_i^q=1$。由 Young 不等式
$$
ab\le \frac{a^p}{p}+\frac{b^q}{q}
$$
可得
$$
\sum_i a_i b_i\le \frac1p\sum_i a_i^p+\frac1q\sum_i b_i^q=1.
$$
乘回 $\left\lVert x\right\rVert_p\left\lVert y\right\rVert_q$ 即得结论。

证明完毕.

**定理（Minkowski 不等式）**:

对 $1\le p<\infty$，有
$$
\left\lVert x+y\right\rVert_p\le\left\lVert x\right\rVert_p+
\left\lVert y\right\rVert_p.
$$

**证明**:

$p=1$ 时由绝对值三角不等式逐项相加得到。设 $1<p<\infty$，令 $q$ 为共轭指数。若 $x+y=0$ 显然。否则：
$$
\left\lVert x+y\right\rVert_p^p
=\sum_i \left\lvert x_i+y_i\right\rvert^p
\le \sum_i (\left\lvert x_i\right\rvert+\left\lvert y_i\right\rvert)\left\lvert x_i+y_i\right\rvert^{p-1}.
$$
对两项分别用 Holder 不等式：
$$
\sum_i \left\lvert x_i\right\rvert\left\lvert x_i+y_i\right\rvert^{p-1}
\le \left\lVert x\right\rVert_p
\left(\sum_i \left\lvert x_i+y_i\right\rvert^{(p-1)q}\right)^{1/q}.
$$
因为 $(p-1)q=p$，所以上式为
$$
\left\lVert x\right\rVert_p\left\lVert x+y\right\rVert_p^{p/q}.
$$
同理得到 $y$ 项。因此
$$
\left\lVert x+y\right\rVert_p^p
\le (\left\lVert x\right\rVert_p+\left\lVert y\right\rVert_p)\left\lVert x+y\right\rVert_p^{p/q}.
$$
由于 $p-p/q=1$，两边除以 $\left\lVert x+y\right\rVert_p^{p/q}$，得到结论。

证明完毕.

### 范数的等价性

**定义（范数等价）**:

同一线性空间 $V$ 上的两个范数 $\left\lVert \cdot\right\rVert_a$ 与 $\left\lVert \cdot\right\rVert_b$ 称为等价，如果存在常数 $c,C>0$，使对所有 $x\in V$，
$$
c\left\lVert x\right\rVert_a\le\left\lVert x\right\rVert_b\le C\left\lVert x\right\rVert_a.
$$

**定理（有限维范数等价）**:

有限维线性空间上的任意两个范数都等价。

**证明**:

只需证明任意范数 $\left\lVert \cdot\right\rVert$ 与 $\left\lVert \cdot\right\rVert_2$ 等价。取一组基，把空间识别为 $\mathbb F^n$。设标准基为 $e_1,
\ldots,e_n$。若 $x=\sum x_ie_i$，则
$$
\left\lVert x\right\rVert\le\sum_i\left\lvert x_i\right\rvert\left\lVert e_i\right\rVert
\le\left(\sum_i\left\lVert e_i\right\rVert^2\right)^{1/2}\left\lVert x\right\rVert_2.
$$
故有上界。

函数 $f(x)=\left\lVert x\right\rVert$ 在 Euclidean 单位球面
$$
S=\{x:\left\lVert x\right\rVert_2=1\}
$$
上连续，而 $S$ 紧。由正定性，$f(x)>0$ 在 $S$ 上成立，故存在最小值 $m>0$。于是当 $x\ne0$ 时，
$$
\left\lVert x\right\rVert=\left\lVert x\right\rVert_2\left\lVert \frac{x}{\left\lVert x\right\rVert_2}\right\rVert
\ge m\left\lVert x\right\rVert_2.
$$
综上得到等价性。

证明完毕.

#### 常见范数之间的估计

对 $x\in\mathbb F^n$，有
$$
\left\lVert x\right\rVert_\infty\le\left\lVert x\right\rVert_2\le\left\lVert x\right\rVert_1\le\sqrt n\left\lVert x\right\rVert_2\le n\left\lVert x\right\rVert_\infty.
$$
更一般地，若 $1\le p\le q\le\infty$，则
$$
\left\lVert x\right\rVert_q\le\left\lVert x\right\rVert_p\le n^{1/p-1/q}\left\lVert x\right\rVert_q.
$$
这些估计常用于误差控制：当维数固定时，选择不同范数只会改变常数；当维数增长时，常数可能依赖维数，需要额外注意。

### 由范数诱导的距离与收敛

范数自然诱导距离
$$
d(x,y)=\left\lVert x-y\right\rVert.
$$
于是可以定义向量序列收敛：
$$
x_k\to x
\quad\Longleftrightarrow\quad
\left\lVert x_k-x\right\rVert\to0.
$$
在有限维空间中，由范数等价性可知：若某个范数下收敛，则任意范数下都收敛。因此有限维线性代数中，不同范数主要影响估计常数，而不改变极限本身。

### 矩阵范数与诱导范数

虽然本章标题为向量范数，但矩阵范数与向量范数密不可分。

**定义（诱导矩阵范数）**:

给定 $\mathbb F^n$ 与 $\mathbb F^m$ 上的向量范数，矩阵 $A\in\mathbb F^{m\times n}$ 的诱导范数定义为
$$
\left\lVert A\right\rVert=\sup_{x\ne0}\frac{\left\lVert Ax\right\rVert}{\left\lVert x\right\rVert}
=\sup_{\left\lVert x\right\rVert=1}\left\lVert Ax\right\rVert.
$$

由定义立即有
$$
\left\lVert Ax\right\rVert\le\left\lVert A\right\rVert\left\lVert x\right\rVert.
$$
同时诱导范数满足次乘性：
$$
\left\lVert AB\right\rVert\le\left\lVert A\right\rVert\left\lVert B\right\rVert.
$$

**命题（常见诱导范数）**:

对矩阵 $A=(a_{ij})$，有
$$
\left\lVert A\right\rVert_1=\max_j\sum_i\left\lvert a_{ij}\right\rvert,
$$
即最大列和；
$$
\left\lVert A\right\rVert_\infty=\max_i\sum_j\left\lvert a_{ij}\right\rvert,
$$
即最大行和；
$$
\left\lVert A\right\rVert_2=\sigma_{\max}(A),
$$
即最大奇异值。

**证明**:

$1$ 范数情形：
$$
\left\lVert Ax\right\rVert_1=\sum_i\left\lvert \sum_j a_{ij}x_j\right\rvert
\le\sum_j\left(\sum_i\left\lvert a_{ij}\right\rvert\right)\left\lvert x_j\right\rvert
\le\left(\max_j\sum_i\left\lvert a_{ij}\right\rvert\right)\left\lVert x\right\rVert_1.
$$
取达到最大列和的列对应的单位向量，可取等号。

$\infty$ 范数类似。二范数情形：
$$
\left\lVert Ax\right\rVert_2^2=x^*A^*Ax\le\lambda_{\max}(A^*A)x^*x,
$$
故 $\left\lVert A\right\rVert_2\le\sqrt{\lambda_{\max}(A^*A)}=\sigma_{\max}$。取最大特征值对应的单位特征向量可取等号。

证明完毕.

**定义（Frobenius 范数）**:

矩阵 $A=(a_{ij})$ 的 Frobenius 范数定义为
$$
\left\lVert A\right\rVert_F=\left(\sum_{i,j}\left\lvert a_{ij}\right\rvert^2\right)^{1/2}
=\sqrt{\operatorname{tr}(A^*A)}.
$$
若 $\sigma_1,
\ldots,\sigma_r$ 是 $A$ 的奇异值，则
$$
\left\lVert A\right\rVert_F=\left(\sum_{i=1}^r\sigma_i^2\right)^{1/2}.
$$

### 范数与条件数

若 $A$ 可逆，线性方程 $Ax=b$ 的解为 $x=A^{-1}b$。若右端扰动为 $b+\delta b$，解扰动为 $x+\delta x$，则
$$
A\delta x=\delta b.
$$
于是
$$
\frac{\left\lVert \delta x\right\rVert}{\left\lVert x\right\rVert}
\le
\left\lVert A^{-1}\right\rVert\left\lVert \delta b\right\rVert\frac{1}{\left\lVert x\right\rVert}.
$$
又因为 $\left\lVert b\right\rVert=\left\lVert Ax\right\rVert\le\left\lVert A\right\rVert\left\lVert x\right\rVert$，所以
$$
\frac{\left\lVert \delta x\right\rVert}{\left\lVert x\right\rVert}
\le
\left\lVert A\right\rVert\left\lVert A^{-1}\right\rVert
\frac{\left\lVert \delta b\right\rVert}{\left\lVert b\right\rVert}.
$$
定义条件数
$$
\operatorname{cond}(A)=\left\lVert A\right\rVert\left\lVert A^{-1}\right\rVert.
$$
在二范数下，若 $A$ 可逆，则
$$
\operatorname{cond}_2(A)=\frac{\sigma_{\max}(A)}{\sigma_{\min}(A)}.
$$
条件数越大，线性方程组对扰动越敏感。

## 梯度矩阵与矩阵微积分

矩阵微积分的核心是把多变量函数的微分写成矩阵形式。它在最小二乘、优化、机器学习、控制理论和统计中非常常见。本章采用“微分优先”的记法：先写出 $df$，再从内积形式中读出梯度。这样比逐元素求导更清晰，也更不容易出错。

### 标量函数的梯度

设 $f:\mathbb R^n\to\mathbb R$ 可微，$x=(x_1,
\ldots,x_n)^T$。梯度定义为
$$
\nabla f(x)=
\begin{bmatrix}
\frac{\partial f}{\partial x_1}\\
\vdots\\
\frac{\partial f}{\partial x_n}
\end{bmatrix}.
$$
微分满足
$$
df=(\nabla f(x))^Tdx.
$$
也就是说，梯度是使一阶变化写成内积的那个向量。

**例**:

若
$$
f(x)=a^Tx,
$$
则
$$
df=a^Tdx,
$$
所以
$$
\nabla f(x)=a.
$$
若
$$
f(x)=x^Ta,
$$
结论相同。

**例**:

若
$$
f(x)=x^TAx,
$$
其中 $A$ 为常矩阵，则
$$
df=d(x^TAx)=(dx)^TAx+x^TA(dx).
$$
把两项都写成 $\cdot^Tdx$ 的形式：
$$
(dx)^TAx=(A^Tx)^Tdx,
\qquad
x^TA(dx)=(A x)^Tdx.
$$
因此
$$
df=((A^T+A)x)^Tdx,
$$
所以
$$
\nabla f(x)=(A+A^T)x.
$$
若 $A$ 对称，则
$$
\nabla(x^TAx)=2Ax.
$$

### Jacobian 矩阵

**定义（Jacobian）**:

设 $F:\mathbb R^n\to\mathbb R^m$，
$$
F(x)=
\begin{bmatrix}
f_1(x)\\
\vdots\\
f_m(x)
\end{bmatrix}.
$$
其 Jacobian 矩阵定义为
$$
J_F(x)=
\begin{bmatrix}
\frac{\partial f_1}{\partial x_1}&\cdots&\frac{\partial f_1}{\partial x_n}\\
\vdots&\ddots&\vdots\\
\frac{\partial f_m}{\partial x_1}&\cdots&\frac{\partial f_m}{\partial x_n}
\end{bmatrix}
\in\mathbb R^{m\times n}.
$$

Jacobian 的意义是给出一阶线性近似：
$$
F(x+h)=F(x)+J_F(x)h+o(\left\lVert h\right\rVert).
$$
若 $m=1$，则 $J_f(x)$ 是行向量
$$
J_f(x)=(\nabla f(x))^T.
$$

**例**:

若
$$
F(x)=Ax+b,
$$
其中 $A\in\mathbb R^{m\times n}$，则
$$
J_F(x)=A.
$$
若
$$
F(x)=
\begin{bmatrix}
x_1^2+x_2\\
\sin x_1+x_2^3
\end{bmatrix},
$$
则
$$
J_F(x)=
\begin{bmatrix}
2x_1&1\\
\cos x_1&3x_2^2
\end{bmatrix}.
$$

### Hessian 矩阵

若 $f:\mathbb R^n\to\mathbb R$ 二阶可微，定义 Hessian 矩阵
$$
\nabla^2 f(x)=
\begin{bmatrix}
\frac{\partial^2 f}{\partial x_1\partial x_1}&\cdots&\frac{\partial^2 f}{\partial x_1\partial x_n}\\
\vdots&\ddots&\vdots\\
\frac{\partial^2 f}{\partial x_n\partial x_1}&\cdots&\frac{\partial^2 f}{\partial x_n\partial x_n}
\end{bmatrix}.
$$
在二阶偏导连续时，Hessian 对称。

Taylor 展开可写为
$$
f(x+h)=f(x)+\nabla f(x)^Th+\frac12 h^T\nabla^2 f(x)h+o(\left\lVert h\right\rVert^2).
$$
Hessian 决定函数的局部曲率。若 $\nabla^2 f(x)$ 正定，则 $x$ 附近函数向各方向上弯；若负定，则向各方向下弯。

### 矩阵变量函数的微分

为了对矩阵变量求导，需要在矩阵空间上选择内积。最常用的是 Frobenius 内积：
$$
\left\langle A,B\right\rangle_F=\operatorname{tr}(A^TB).
$$
若 $f:\mathbb R^{m\times n}\to\mathbb R$ 可微，则梯度矩阵 $\nabla_X f$ 定义为满足
$$
df=\operatorname{tr}((\nabla_X f)^T dX).
$$
这与向量情形的
$$
df=(\nabla f)^Tdx
$$
完全一致。

#### 迹技巧

矩阵求导中经常用以下恒等式：
$$
\operatorname{tr}(AB)=\operatorname{tr}(BA),
$$
只要乘法阶数相容；更一般地，
$$
\operatorname{tr}(A_1A_2\cdots A_k)=\operatorname{tr}(A_kA_1\cdots A_{k-1}).
$$
此外，
$$
\operatorname{tr}(A^T dX)=\left\langle A,dX\right\rangle_F.
$$
因此只要把 $df$ 化成 $\operatorname{tr}(G^T dX)$，就可以读出
$$
\nabla_X f=G.
$$

**例**:

设
$$
f(X)=\operatorname{tr}(A^TX).
$$
则
$$
df=\operatorname{tr}(A^T dX),
$$
所以
$$
\nabla_X f=A.
$$

**例**:

设
$$
f(X)=\operatorname{tr}(AXB),
$$
其中 $A,B$ 为常矩阵。则
$$
df=\operatorname{tr}(A(dX)B)=\operatorname{tr}(BA dX).
$$
要写成 $\operatorname{tr}(G^T dX)$，需有
$$
G^T=BA,
$$
所以
$$
\nabla_X f=(BA)^T=A^TB^T.
$$

**例（Frobenius 范数平方）**:

设
$$
f(X)=\frac12\left\lVert X\right\rVert_F^2=\frac12\operatorname{tr}(X^TX).
$$
则
$$
df=\frac12\operatorname{tr}((dX)^TX+X^TdX)
=\operatorname{tr}(X^TdX),
$$
因此
$$
\nabla_X f=X.
$$

**例（矩阵最小二乘）**:

设
$$
f(X)=\frac12\left\lVert AX-B\right\rVert_F^2.
$$
令残差
$$
R=AX-B.
$$
则
$$
f=\frac12\operatorname{tr}(R^TR),
\qquad dR=A dX.
$$
于是
$$
df=\operatorname{tr}(R^T dR)=\operatorname{tr}(R^T A dX)=\operatorname{tr}((A^TR)^T dX).
$$
所以
$$
\nabla_X f=A^T(AX-B).
$$
令梯度为零，得到正规方程
$$
A^TAX=A^TB.
$$

**例（双侧最小二乘）**:

设
$$
f(X)=\frac12\left\lVert AXB-C\right\rVert_F^2.
$$
令 $R=AXB-C$，则
$$
dR=A(dX)B.
$$
因此
$$
df=\operatorname{tr}(R^T A(dX)B)=\operatorname{tr}(BR^TA dX).
$$
所以
$$
\nabla_X f=A^T(AXB-C)B^T.
$$

### 常用矩阵求导公式

以下公式均假设矩阵阶数相容。

| 函数 | 梯度 |
| --- | --- |
| $f(x)=a^Tx$ | $\nabla f=a$ |
| $f(x)=x^TAx$ | $\nabla f=(A+A^T)x$ |
| $f(x)=\frac12\left\lVert Ax-b\right\rVert_2^2$ | $\nabla f=A^T(Ax-b)$ |
| $f(X)=\operatorname{tr}(A^TX)$ | $\nabla_X f=A$ |
| $f(X)=\operatorname{tr}(AXB)$ | $\nabla_X f=A^TB^T$ |
| $f(X)=\frac12\left\lVert X\right\rVert_F^2$ | $\nabla_X f=X$ |
| $f(X)=\frac12\left\lVert AX-B\right\rVert_F^2$ | $\nabla_X f=A^T(AX-B)$ |
| $f(X)=\frac12\left\lVert AXB-C\right\rVert_F^2$ | $\nabla_X f=A^T(AXB-C)B^T$ |
| $f(X)=\log\det X$ | $\nabla_X f=X^{-T}$ |
| $f(X)=\operatorname{tr}(X^{-1}A)$ | $\nabla_X f=-X^{-T}A^TX^{-T}$ |

**证明（关于 $\log\det X$ 的证明）**:

利用行列式微分公式
$$
d\det X=\det X \operatorname{tr}(X^{-1}dX),
$$
可得
$$
d\log\det X=\frac{1}{\det X}d\det X
=\operatorname{tr}(X^{-1}dX)=\operatorname{tr}((X^{-T})^TdX).
$$
因此
$$
\nabla_X\log\det X=X^{-T}.
$$

证明完毕.

**证明（关于 $\operatorname{tr}(X^{-1}A)$ 的证明）**:

由
$$
X X^{-1}=I
$$
微分得
$$
dX X^{-1}+X d(X^{-1})=0,
$$
所以
$$
d(X^{-1})=-X^{-1}(dX)X^{-1}.
$$
于是
$$
d\operatorname{tr}(X^{-1}A)
=\operatorname{tr}(d(X^{-1})A)
=-\operatorname{tr}(X^{-1}(dX)X^{-1}A)
=-\operatorname{tr}(X^{-1}AX^{-1}dX).
$$
因此梯度 $G$ 满足
$$
G^T=-X^{-1}AX^{-1},
$$
即
$$
G=-X^{-T}A^TX^{-T}.
$$

证明完毕.

### 链式法则

若 $y=g(x)\in\mathbb R^m$，$f=f(y)\in\mathbb R$，则
$$
\nabla_x(f\circ g)(x)=J_g(x)^T\nabla_y f(y).
$$
这是多元链式法则的矩阵形式。

**例**:

设
$$
f(x)=\frac12\left\lVert g(x)\right\rVert_2^2.
$$
则
$$
\nabla f(x)=J_g(x)^Tg(x).
$$
特别地，当 $g(x)=Ax-b$ 时，
$$
\nabla f(x)=A^T(Ax-b).
$$

### 二次优化中的梯度与 Hessian

考虑二次函数
$$
f(x)=\frac12x^TAx-b^Tx+c,
$$
其中 $A=A^T$。则
$$
\nabla f(x)=Ax-b,
\qquad
\nabla^2 f(x)=A.
$$
若 $A$ 正定，则 $f$ 严格凸，唯一极小点由
$$
Ax=b
$$
给出。若 $A$ 半正定，则 $f$ 凸，但极小点可能不唯一；若 $A$ 有负特征值，则 $f$ 沿对应方向向下无界。

**定理（最小二乘的正规方程）**:

对
$$
f(x)=\frac12\left\lVert Ax-b\right\rVert_2^2,
$$
有
$$
\nabla f(x)=A^T(Ax-b).
$$
因此极小点满足
$$
A^TAx=A^Tb.
$$
若 $A$ 满列秩，则 $A^TA$ 正定，极小点唯一：
$$
x=(A^TA)^{-1}A^Tb.
$$

**证明**:

梯度公式由前面的矩阵微分得到。若 $A$ 满列秩，则对任意 $x\ne0$，
$$
x^TA^TAx=\left\lVert Ax\right\rVert_2^2>0,
$$
所以 $A^TA$ 正定，正规方程有唯一解。凸二次函数的驻点即全局极小点。

证明完毕.

\backmatter

\chapter*{常用结论索引}
\addcontentsline{toc}{chapter}{常用结论索引}

1. 子空间判别：非空子集对线性组合封闭即为子空间。
2. 维数公式：$\dim(U+W)=\dim U+\dim W-\dim(U\cap W)$。
3. 直和判别：分解唯一等价于交为零。
4. 秩-零度定理：$\operatorname{rank} T+\operatorname{nullity} T=\dim V$。
5. 基变换：同一线性变换在不同基下的矩阵相似，$A'=P^{-1}AP$。
6. 伴随矩阵：标准正交基下，伴随对应共轭转置 $A^*$。
7. 投影判别：$P$ 是投影当且仅当 $P^2=P$；正交投影还满足 $P^*=P$。
8. Kronecker 积：$(A\otimes B)(C\otimes D)=(AC)\otimes(BD)$，$\operatorname{vec}(AXB)=(B^T\otimes A)\operatorname{vec}(X)$。
9. Jordan 块大小：$d_k-d_{k-1}$ 等于大小至少为 $k$ 的 Jordan 块个数。
10. 最小多项式：指数等于相应特征值的最大 Jordan 块大小。
11. 对角化判别：最小多项式无重根等价于可对角化。
12. QR 分解：满列秩矩阵 $A=QR$，最小二乘化为 $Rx=Q^Tb$。
13. 谱分解：Hermite 矩阵可酉对角化，实对称矩阵可正交对角化。
14. SVD：任意矩阵 $A=U\Sigma V^*$，奇异值描述正交方向伸缩。
15. 伪逆：$A^+=V\Sigma^+U^*$ 给出最小二乘的最小范数解。
16. 有限维范数等价：所有范数定义相同的收敛概念。
17. 矩阵梯度：把 $df$ 化为 $\operatorname{tr}(G^T dX)$，则 $\nabla_X f=G$。
