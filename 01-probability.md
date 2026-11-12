# 概率论与统计推断

设有限总体 $U=\{1,\ldots,N\}$，单位 $i$ 的观测值为 $y_i$，目标参数为总体总量 $T=\sum_{i=1}^N y_i$。抽样设计给出每个可能样本 $s\subseteq U$ 的概率 $p(s)$。在设计推断中，$U$ 和 $(y_1,\ldots,y_N)$ 固定，估计量的分布由 $p(s)$ 决定。

## 1.1 样本空间、边际分布与条件分布

从 $N=20$ 个街区中简单随机不放回抽取两个街区。样本空间包含 $\binom{20}{2}=190$ 个无序对，概率均为 $1/190$。令 $I_i=\mathbf{1}\{i\in s\}$。包含 $i$ 的无序对共有 19 个，故 $P(I_i=1)=19/190=1/10$。对任意 $i\ne j$，联合入样概率为 $P(I_i=1,I_j=1)=1/190$。

由条件概率的定义，
$$
P(I_j=1\mid I_i=1)=\frac{P(I_i=1,I_j=1)}{P(I_i=1)}=\frac{1/190}{1/10}=\frac1{19}.
$$
由于 $1/19<1/10$，两个街区的入样指示变量负相关。独立性要求联合分布等于边际分布的乘积。零协方差只要求 $E[XY]=E[X]E[Y]$。例如，若 $X$ 在 $\{-1,0,1\}$ 上等概率取值，且 $Y=X^2$，则 $E[XY]=E[X^3]=0$，而 $P(Y=0\mid X=0)=1$，$P(Y=0)=1/3$。

## 1.2 期望、条件期望和方差分解

设 $X$ 为类别，$Y$ 为年支出，且 $E[|Y|]<\infty$。在离散情形下，对联合分布按 $X$ 分组可得
$$
E[Y]=\sum_x\sum_y yP(Y=y\mid X=x)P(X=x)
=\sum_x E[Y\mid X=x]P(X=x).
$$
连续情形以积分代替求和。若两类人口占比分别为 0.75、0.25，条件均值分别为 2、10，则总体均值为 $0.75(2)+0.25(10)=4$。等额分层抽样下，直接合并两个层的样本得到的均值期望为 $(2+10)/2=6$。

设 $E[Y^2]<\infty$，记 $m(X)=E[Y\mid X]$。分解 $Y-E[Y]=Y-m(X)+m(X)-E[Y]$。由于 $E[Y-m(X)\mid X]=0$，两项乘积的条件期望为零。对两边平方并取期望，得到
$$
\operatorname{Var}(Y)=E[\operatorname{Var}(Y\mid X)]+\operatorname{Var}(m(X)).
$$
第一项是组内方差的总体加权平均，第二项是组均值的方差。若上述两组的条件方差分别为 4 和 16，则两项分别为 $0.75(4)+0.25(16)=7$ 和 $0.75(2-4)^2+0.25(10-4)^2=12$，总方差为 19。

## 1.3 估计量、标准误和一致性

设 $Y_1,\ldots,Y_n$ 具有共同均值 $\mu$，$\bar Y=n^{-1}\sum_iY_i$。由期望的线性性，$E[\bar Y]=\mu$。估计量 $\hat\theta$ 的标准误定义为其抽样分布的标准差。将 $\hat\theta-\theta$ 分解为 $\hat\theta-E[\hat\theta]$ 与 $E[\hat\theta]-\theta$ 两项，平方取期望后得 $E[(\hat\theta-\theta)^2]=\operatorname{Var}(\hat\theta)+(E[\hat\theta]-\theta)^2$。

对任何具有二阶矩的样本，
$$
\operatorname{Var}(\bar Y)=\frac1{n^2}
\left[\sum_i\operatorname{Var}(Y_i)+2\sum_{i<j}\operatorname{Cov}(Y_i,Y_j)\right].
$$
独立同分布时，非对角协方差均为零，故 $\operatorname{Var}(\bar Y)=\sigma^2/n$。样本量增至四倍，标准误降至原值的一半。若边际方差均为 $\sigma^2$，且不同观测的相关系数均为 $\rho>0$，则 $n(n-1)$ 个非对角项各贡献 $\rho\sigma^2/n^2$，从而 $\operatorname{Var}(\bar Y)=\sigma^2[1+(n-1)\rho]/n$。

独立同分布且方差有限时，对任意 $\epsilon>0$，切比雪夫不等式给出 $P(|\bar Y-\mu|>\epsilon)\le \sigma^2/(n\epsilon^2)\to0$，故 $\bar Y\xrightarrow{p}\mu$。若 $\sigma^2>0$，中心极限定理进一步给出 $\sqrt n(\bar Y-\mu)/\sigma\Rightarrow N(0,1)$。

## 1.4 检验、区间与功效

令 $s^2=(n-1)^{-1}\sum_i(Y_i-\bar Y)^2$。在上述独立同分布条件下，$s^2\xrightarrow{p}\sigma^2$。由 Slutsky 定理，$(\bar Y-\mu)/(s/\sqrt n)\Rightarrow N(0,1)$，反演该统计量得到渐近 95% 置信区间 $\bar Y\pm1.96s/\sqrt n$。其覆盖率为重复抽样时区间覆盖固定参数 $\mu$ 的概率。若样本来自正态总体，该统计量精确服从 $t_{n-1}$ 分布。

对 $H_0:\mu=\mu_0$ 的双侧检验，拒绝域为 $|(\bar Y-\mu_0)/SE|>1.96$。在已知标准误的正态模型下，令 $\delta=(\mu_1-\mu_0)/SE$，则真实参数为 $\mu_1$ 时的拒绝概率为 $1-\Phi(1.96-\delta)+\Phi(-1.96-\delta)$。$\delta=1$ 时，检验功效约为 17%。

## 1.5 Horvitz–Thompson 总量估计与街区抽样

Horvitz 与 Thompson（1952）表 2 报告了艾奥瓦州 Ames 市 20 个街区的真实住户数 $y_i$ 和目测住户数 $s_i$，两列合计分别为 434 和 394。$T=\sum_i y_i$ 为目标总量，$s_i$ 为抽样前已知的辅助变量。下文的期望与方差均对抽样设计计算。 [论文原文，表 2](https://www.stat.cmu.edu/~brian/905-2008/papers/Horvitz-Thompson-1952-jasa.pdf)

记一阶入样概率 $\pi_i=P(I_i=1)>0$。已知总体规模 $N$，总量与均值估计分别为
$$
\hat T_{HT}=\sum_{i=1}^{N}I_i\frac{y_i}{\pi_i},\qquad
\hat\mu_{HT}=\frac{\hat T_{HT}}N.
$$
由 $E[I_i]=\pi_i$，有 $E[\hat T_{HT}]=\sum_i E[I_i]y_i/\pi_i=T$。令 $\pi_{ij}=E[I_iI_j]$，其中 $\pi_{ii}=\pi_i$。将 $\operatorname{Cov}(I_i,I_j)=\pi_{ij}-\pi_i\pi_j$ 代入加权和的方差，得
$$
\operatorname{Var}(\hat T_{HT})=\sum_i\sum_j
(\pi_{ij}-\pi_i\pi_j)\frac{y_i}{\pi_i}\frac{y_j}{\pi_j}.
$$
当联合入样概率均为正时，在上式各项乘以 $I_iI_j/\pi_{ij}$ 即得仅依赖入样单位的设计无偏方差估计量。其无偏性由 $E[I_iI_j/\pi_{ij}]=1$ 给出。

**数值计算。** 取原文表 2 的 $y_i$ 和 $s_i$，在每次抽取时按剩余街区的目测住户数比例分配抽取概率。记 $S=\sum_i s_i$。样本 $\{i,j\}$ 由先抽 $i$ 或先抽 $j$ 两个互斥事件构成，故
$$
q_{ij}=\frac{s_i}{S}\frac{s_j}{S-s_i}+\frac{s_j}{S}\frac{s_i}{S-s_j}.
$$
一阶入样概率为 $\pi_i=\sum_{j\ne i}q_{ij}$。下表按上述两步设计枚举全部 190 个无序对，设计均值和方差的权重均为 $q_{ij}$。

| 设计与估计量 | 设计期望（户） | 偏差（户） | 设计标准差（户） |
|---|---:|---:|---:|
| 等概率抽样，样本平均乘 20 | 434.000 | 0.000 | 127.361 |
| 按目测规模抽样，样本平均乘 20 | 494.461 | 60.461 | 132.984 |
| 按目测规模抽样，HT 估计量 | 434.000 | 0.000 | 54.780 |

在该设计下，以 $N\bar y_s$ 估计总量的偏差为 $\sum_i(N\pi_i/2-1)y_i=60.461$ 户。HT 权重使每个街区对总量估计的期望贡献恢复为 $y_i$。若某街区未列入抽样框，$\pi_i=0$，其住户数对目标总量的贡献在该设计下不可识别。

第一次抽取概率 $s_i/S$ 满足 $\sum_i s_i/S=1$。最终入样概率则满足 $\sum_i\pi_i=E[\sum_iI_i]=2$，其中包括第二次抽取的概率。
