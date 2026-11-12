# 概率论与统计推断

设有限总体 $U=\{1,\ldots,N\}$，单位 $i$ 的观测值为 $y_i$，目标参数为总体总量 $T=\sum_{i=1}^N y_i$。抽样设计给出每个可能样本 $s\subseteq U$ 的概率 $p(s)$。在设计推断中，$U$ 和 $(y_1,\ldots,y_N)$ 固定，估计量的分布由 $p(s)$ 决定。

## 1.1 样本空间、边际分布与条件分布
**从抽样规则构造概率。** 从 $N=20$ 个街区中简单随机不放回抽取两个。有序的抽取路径共 $20\times19$ 条，每个无序对 $\{i,j\}$ 对应 $(i,j)$ 与 $(j,i)$ 两条路径。因此无序样本数是 $20\times19/2=190$，每个样本的概率是 $2/(20\times19)=1/190$。

令 $I_i=\mathbf 1\{i\in s\}$。固定街区 $i$，另一街区有 19 种选择。因此 $P(I_i=1)=19/190=1/10$。再固定另一个街区 $j$，同时包含 $i,j$ 的样本只有一个。两街区的边际概率、联合概率和条件概率由此给出：
$$
\begin{aligned}
E[I_i]&=1\cdot P(I_i=1)+0\cdot P(I_i=0)=1/10,\\
E[I_iI_j]&=P(I_i=1,I_j=1)=1/190,\\
P(I_j=1\mid I_i=1)&=\frac{P(I_i=1,I_j=1)}{P(I_i=1)}=\frac{1/190}{1/10}=1/19 .
\end{aligned}
$$
固定 $i$ 入样后，只余一个名额分配给其他 19 个街区。由协方差定义，
$$
\operatorname{Cov}(I_i,I_j)=E[I_iI_j]-E[I_i]E[I_j]
=\frac1{190}-\frac1{100}=-\frac9{1900}.
$$

**独立性与正交性。** 若 $X,Y$ 独立且二阶矩存在，离散情形下
$$
E[XY]=\sum_{x,y}xyP(X=x)P(Y=y)
=\left(\sum_x xP(X=x)\right)\left(\sum_y yP(Y=y)\right)=E[X]E[Y].
$$
反向命题可用三点分布检验。设 $X$ 在 $\{-1,0,1\}$ 上等概率取值，$Y=X^2$。逐项求和得 $E[X]=0$、$E[Y]=2/3$、$E[XY]=(-1+0+1)/3=0$，故协方差为零。但 $P(Y=0\mid X=0)=1$，而 $P(Y=0)=1/3$，联合分布不满足独立性。

## 1.2 期望、条件期望和方差分解
**全期望公式。** 设 $E|Y|<\infty$。在离散情形，先以 $X$ 的取值划分样本空间，再使用条件概率的定义。对 $P(X=x)>0$ 的层，
$$
\begin{aligned}
E[Y]&=\sum_y yP(Y=y)\\
&=\sum_y y\sum_xP(Y=y,X=x)\\
&=\sum_xP(X=x)\sum_y yP(Y=y\mid X=x)\\
&=\sum_xP(X=x)E[Y\mid X=x].
\end{aligned}
$$
绝对可积性允许交换求和顺序。连续情形写作 $E[Y]=\int m(x)\,dF_X(x)$。设两层人口占比分别为 $0.75,0.25$，层均值为 $2,10$。总体均值为 $0.75\times2+0.25\times10=4$。若每层简单随机抽取同样多的人，合并样本均值为 $(\bar Y_1+\bar Y_2)/2$，其期望为 $(2+10)/2=6$。将层均值按人口占比加权，期望恢复为 4。

**全方差公式。** 设 $E[Y^2]<\infty$。令 $\mu=E[Y]$、$m=E[Y\mid X]$、$A=Y-m$、$B=m-\mu$。由于 $B$ 是 $X$ 的函数，给定 $X$ 时可将其提出期望：
$$
E[AB]=E\{E[AB\mid X]\}=E\{B(E[Y\mid X]-m)\}=0.
$$
再分别计算平方项。由条件方差定义，$E[A^2]=E\{\operatorname{Var}(Y\mid X)\}$。由全期望公式，$E[m]=\mu$，故 $E[B^2]=\operatorname{Var}(m)$。将三项代回：
$$
\operatorname{Var}(Y)=E[(A+B)^2]=E[A^2]+2E[AB]+E[B^2]
=E[\operatorname{Var}(Y\mid X)]+\operatorname{Var}(E[Y\mid X]).
$$
若上述两层的条件方差为 4 与 16，则组内项为 $0.75\times4+0.25\times16=7$；组间项为 $0.75(2-4)^2+0.25(10-4)^2=12$。总方差为 19。

## 1.3 估计量、标准误和一致性
**偏差与均方误差。** 对目标 $\theta$，令 $b=E[\hat\theta]-\theta$，$v=\hat\theta-E[\hat\theta]$。由定义 $E[v]=0$，而 $b$ 是常数。因此
$$
E[(\hat\theta-\theta)^2]=E[(v+b)^2]=E[v^2]+2bE[v]+b^2
=\operatorname{Var}(\hat\theta)+b^2.
$$
标准误是 $\sqrt{\operatorname{Var}(\hat\theta)}$。设 $Y_i$ 具有共同均值 $\mu$。样本均值 $\bar Y=n^{-1}\sum_iY_i$ 的期望为 $n^{-1}\sum_i E[Y_i]=n\mu/n=\mu$。计算方差时，求和中的相关性进入二阶矩。写 $U_i=Y_i-\mu$，则
$$
\begin{aligned}
\operatorname{Var}(\bar Y)
&=E\left[\left(\frac1n\sum_i U_i\right)^2\right]\\
&=\frac1{n^2}\sum_i\sum_j E[U_iU_j]\\
&=\frac1{n^2}\left(\sum_i\operatorname{Var}(Y_i)+2\sum_{i<j}\operatorname{Cov}(Y_i,Y_j)\right).
\end{aligned}
$$
最后一行将 $i=j$ 的对角项与 $(i,j),(j,i)$ 成对的非对角项分开。独立同分布时，仅余 $n$ 个方差 $\sigma^2$，故 $\operatorname{Var}(\bar Y)=n\sigma^2/n^2=\sigma^2/n$。样本量从 $n$ 增至 $4n$，标准误的比值是 $\sqrt{n/(4n)}=1/2$。

若任意两项之间的相关系数均为 $\rho>0$ 且边际方差相同，则每个非对角项为 $\rho\sigma^2$。由有序非对角项计数 $n(n-1)$，
$$
\operatorname{Var}(\bar Y)=\frac{n\sigma^2+n(n-1)\rho\sigma^2}{n^2}
=\sigma^2\left(\rho+\frac{1-\rho}{n}\right).
$$
当 $\rho$ 固定时，方差极限为 $\rho\sigma^2$。增加相关观测与增加独立观测对精度的影响由此不同。

**从方差到一致性。** 在独立同分布、有限正方差情形，对任意 $\epsilon>0$ 有 $\epsilon^2\mathbf1\{|\bar Y-\mu|>\epsilon\}\le(\bar Y-\mu)^2$。两侧取期望，再除以 $\epsilon^2$，得
$$
P(|\bar Y-\mu|>\epsilon)\le\frac{E[(\bar Y-\mu)^2]}{\epsilon^2}=\frac{\sigma^2}{n\epsilon^2}\longrightarrow0.
$$
这就是依概率收敛的定义。对标准化和应用独立同分布中心极限定理，得到 $\sum_i(Y_i-\mu)/(\sigma\sqrt n)=\sqrt n(\bar Y-\mu)/\sigma\Rightarrow N(0,1)$。

**标准正态极限的推导。** 将单个观测标准化为 $W=(Y-\mu)/\sigma$，使 $E[W]=0$、$E[W^2]=1$。记特征函数 $\phi(t)=E[e^{itW}]$。有限二阶矩允许在期望内求前两阶导数，$\phi(0)=1$、$\phi'(0)=iE[W]=0$、$\phi''(0)=-E[W^2]=-1$。在零点作二阶展开得 $\phi(t)=1-t^2/2+o(t^2)$。独立性使和的特征函数为各项特征函数之积，因此对固定 $t$，
$$
\begin{aligned}
E\!\left[\exp\left(it\frac{\sum_j W_j}{\sqrt n}\right)\right]
&=\prod_{j=1}^n E[\exp(itW_j/\sqrt n)]\\
&=\{\phi(t/\sqrt n)\}^n
=\left(1-\frac{t^2}{2n}+o(n^{-1})\right)^n\longrightarrow e^{-t^2/2}.
\end{aligned}
$$
最后一步用 $n\log(1+a_n/n)\to a$（$a_n\to a$）。右侧是标准正态特征函数，由特征函数连续性定理得到上述分布收敛。

## 1.4 检验、区间与功效
**估计标准误。** 将平方和展开并利用 $\sum_iY_i=n\bar Y$：
$$
s^2=\frac{\sum_i(Y_i-\bar Y)^2}{n-1}
=\frac{\sum_i Y_i^2-2\bar Y\sum_iY_i+n\bar Y^2}{n-1}
=\frac{n}{n-1}\left(\frac1n\sum_iY_i^2-\bar Y^2\right).
$$
由大数定律，两项分别趋于 $E[Y^2]$ 与 $\mu^2$，故 $s^2\xrightarrow p \sigma^2$。当 $\sigma>0$，$\sigma/s\xrightarrow p1$。把未知方差的统计量写成两个因子的乘积，由 Slutsky 定理，
$$
\frac{\bar Y-\mu}{s/\sqrt n}
=\frac{\sqrt n(\bar Y-\mu)}{\sigma}\frac{\sigma}{s}\Rightarrow N(0,1).
$$

**反演检验。** 令 $z=\Phi^{-1}(0.975)\simeq1.96$。双侧检验接受域的不等式可逐步移项：
$$
-z\le\frac{\bar Y-\mu}{s/\sqrt n}\le z
\iff -z\,s/\sqrt n\le\bar Y-\mu\le z\,s/\sqrt n
\iff \bar Y-zs/\sqrt n\le\mu\le\bar Y+zs/\sqrt n .
$$
该随机区间的覆盖概率趋于 $0.95$。若 $Y_i$ 独立正态，则 $Z=\sqrt n(\bar Y-\mu)/\sigma$ 与 $U=(n-1)s^2/\sigma^2$ 独立，分别服从标准正态和自由度 $n-1$ 的卡方分布。由于 $Z/\sqrt{U/(n-1)}=(\bar Y-\mu)/(s/\sqrt n)$，按 $t$ 分布的定义，使用 $t_{n-1,0.975}$ 得到精确区间。

正态样本中，独立性可从正交分解直接得到。令 $w_i=(Y_i-\mu)/\sigma$，取正交矩阵 $H$ 的第一行为 $\mathbf1'/\sqrt n$，令 $v=Hw$。由于 $w\sim N(0,I_n)$，有 $v\sim N(0,HI_nH')=N(0,I_n)$，故 $v_1,\ldots,v_n$ 为独立标准正态。第一坐标 $v_1=\sqrt n(\bar Y-\mu)/\sigma$。由正交变换保持平方和，
$$
\sum_{j=2}^n v_j^2=\sum_i w_i^2-v_1^2
=\frac{\sum_i(Y_i-\mu)^2-n(\bar Y-\mu)^2}{\sigma^2}
=\frac{\sum_i(Y_i-\bar Y)^2}{\sigma^2}.
$$
故分母所用的平方和是 $n-1$ 个独立标准正态的平方和，与第一坐标独立。这给出了上式中卡方分布与自由度的来源。

**功效。** 在已知 $SE=\sigma/\sqrt n$ 的正态模型下检验 $\mu=\mu_0$。若真实均值为 $\mu_1$，令 $\delta=(\mu_1-\mu_0)/SE$，则检验统计量为 $(\bar Y-\mu_1)/SE+\delta=Z+\delta$。左尾与右尾互斥，因此
$$
\begin{aligned}
P_{\mu_1}(\text{拒绝})&=P(Z+\delta<-z)+P(Z+\delta>z)\\
&=\Phi(-z-\delta)+1-\Phi(z-\delta).
\end{aligned}
$$
取 $\delta=1$，两尾分别约为 $0.0015$ 和 $0.1685$，合计约为 $0.17$。

## 1.5 Horvitz–Thompson 总量估计与街区抽样
Horvitz 与 Thompson（1952）表 2 给出 Ames 市 20 个街区的住户数 $y_i$ 与抽样前目测数 $s_i$，合计分别为 434 和 394。目标是估计 $T=\sum_i y_i$。下列概率均对抽样设计取值，街区的住户数固定。

**构造权重。** 设 $\pi_i=P(I_i=1)>0$，考虑样本总量 $\sum_i I_i w_i y_i$。其期望为 $\sum_i w_i y_i E[I_i]=\sum_i w_i\pi_i y_i$。要使这一式子对任意一组 $y_i$ 都等于 $\sum_i y_i$，每个 $y_i$ 的系数必须为一。因此 $w_i=1/\pi_i$，得到
$$
\hat T=\sum_{i\in s}\frac{y_i}{\pi_i},\qquad E[\hat T]=T,\qquad E[\hat T/N]=T/N .
$$

**从联合入样概率到方差。** 令 $a_i=y_i/\pi_i$，$\pi_{ij}=P(I_i=1,I_j=1)$，$\pi_{ii}=\pi_i$。先将估计误差中心化，再展开平方：
$$
\begin{aligned}
\hat T-T&=\sum_i a_i(I_i-\pi_i),\\
\operatorname{Var}(\hat T)&=\sum_{i,j}a_i a_j E[(I_i-\pi_i)(I_j-\pi_j)],\\
E[(I_i-\pi_i)(I_j-\pi_j)]&=\pi_{ij}-\pi_j\pi_i-\pi_i\pi_j+\pi_i\pi_j
=\pi_{ij}-\pi_i\pi_j .
\end{aligned}
$$
若每一对的 $\pi_{ij}>0$，样本方差估计量为 $\hat V=\sum_{i,j\in s} a_i a_j(\pi_{ij}-\pi_i\pi_j)/\pi_{ij}$。将样本求和改写为总体求和，每个加数前乘 $I_iI_j$。由于 $E[I_iI_j/\pi_{ij}]=1$，逐项取期望恢复上述设计方差。

**复算。** 令 $S=\sum_i s_i$。第一次按 $s_i/S$ 抽取；抽中 $i$ 后，其余街区目测数合计为 $S-s_i$。按乘法公式，路径 $(i,j)$ 的概率为 $(s_i/S)s_j/(S-s_i)$。两种顺序互斥，故无序样本及街区的入样概率依次为
$$
q_{ij}=\frac{s_i}{S}\frac{s_j}{S-s_i}+\frac{s_j}{S}\frac{s_i}{S-s_j},
\qquad \pi_i=\sum_{j\ne i}q_{ij}.
$$
对每个 $i<j$，计算 $t_{ij}=y_i/\pi_i+y_j/\pi_j$。设计期望为 $m=\sum_{i<j}q_{ij}t_{ij}$，标准差为 $\sqrt{\sum_{i<j}q_{ij}(t_{ij}-m)^2}$。等概率设计以 $q_{ij}=1/190$ 代入。下表枚举同一组原始数据的 190 个样本。

| 设计与估计量 | 期望（户） | 偏差（户） | 标准差（户） |
|---|---:|---:|---:|
| 等概率抽样，样本均值乘 20 | 434.000 | 0.000 | 127.361 |
| 按目测规模抽样，样本均值乘 20 | 494.461 | 60.461 | 132.984 |
| 按目测规模抽样，HT 估计量 | 434.000 | 0.000 | 54.780 |

第二行的总量估计为 $(20/2)\sum_i I_i y_i$。其期望减去 $T$ 为 $\sum_i(10\pi_i-1)y_i=60.461$ 户。由 $\sum_i\pi_i=E[\sum_i I_i]=2$，可以核对程序所用权重确为两次抽取后的入样概率。若某街区未列入抽样框，其入样概率为零；任取两个只在该街区住户数上不同的总体，所有可能样本均相同，总量却不同。

[论文原文，表 2](https://www.stat.cmu.edu/~brian/905-2008/papers/Horvitz-Thompson-1952-jasa.pdf)
