# 随机对照试验 RCT

在有限总体随机化模型中，$N$ 个单位及其潜在结果固定。以下期望 $E_Z$ 和方差 $\operatorname{Var}_Z$ 均由处理分配向量 $Z$ 的联合分布定义。

## 3.1 完全随机分配下的设计无偏性
设 $n_1+n_0=N$ 且两组均非空。完全随机分配规定每一种含 $n_1$ 个处理单位的名单概率相同。全部名单有 $\binom{N}{n_1}$ 种。若固定单位 $i$ 入选，则余下 $n_1-1$ 人从 $N-1$ 人中选出。将组合数展开，有
$$
P(Z_i=1)=\frac{\binom{N-1}{n_1-1}}{\binom{N}{n_1}}
=\frac{(N-1)!}{(n_1-1)!\,n_0!}\frac{n_1!\,n_0!}{N!}
=\frac{n_1}{N}.
$$
设实际处理等于分配、无单位间干扰，且处理定义一致。由 $Z_i^2=Z_i$ 和 $Z_i(1-Z_i)=0$，观测结果 $Y_i=Z_iY_i(1)+(1-Z_i)Y_i(0)$ 满足 $Z_iY_i=Z_iY_i(1)$ 与 $(1-Z_i)Y_i=(1-Z_i)Y_i(0)$。因此两组观测均值之差可写成固定潜在结果的随机加权和：
$$
\hat\tau=\frac1{n_1}\sum_i Z_iY_i(1)-\frac1{n_0}\sum_i(1-Z_i)Y_i(0).
$$
组人数由设计固定。用期望的线性性和 $E_Z[1-Z_i]=1-n_1/N=n_0/N$，逐项计算：
$$
\begin{aligned}
E_Z[\hat\tau]&=\frac1{n_1}\sum_iY_i(1)\frac{n_1}{N}
-\frac1{n_0}\sum_iY_i(0)\frac{n_0}{N}\\
&=\frac1N\sum_i [Y_i(1)-Y_i(0)]\equiv\tau_S .
\end{aligned}
$$
$\tau_S$ 对试验名单逐人等权。若研究目标是名单以外人群的平均效应，还需指定名单相对于该人群的抽取机制。

## 3.2 均值差的有限总体方差
**处理指示量的协方差。** 令 $p=n_1/N$。对 $i\ne j$，在 $i$ 已入选的条件下还剩 $n_1-1$ 个处理名额。因此
$$
E_Z[Z_iZ_j]=p\frac{n_1-1}{N-1},\qquad
\operatorname{Cov}_Z(Z_i,Z_j)=p\left(\frac{n_1-1}{N-1}-\frac{n_1}{N}\right)
=-\frac{p(1-p)}{N-1}.
$$
对 $i=j$，由 $Z_i^2=Z_i$ 得 $\operatorname{Var}_Z(Z_i)=p-p^2$。

**两个处理组样本均值的协方差。** 对固定序列 $a_i,b_i$，令 $a_i^*=a_i-\bar a$、$b_i^*=b_i-\bar b$，$S_{ab}=\sum_i a_i^*b_i^*/(N-1)$。由于 $\sum_iZ_i=n_1$，处理组均值 $\bar a_T-\bar a=n_1^{-1}\sum_i Z_i a_i^*$。又因中心化序列的和均为零，
$$
\sum_{i\ne j}a_i^*b_j^*=(\sum_i a_i^*)(\sum_j b_j^*)-\sum_i a_i^*b_i^*=-\sum_i a_i^*b_i^*.
$$
据此把均值协方差拆成同一人和不同人两类项：
$$
\begin{aligned}
\operatorname{Cov}_Z(\bar a_T,\bar b_T)
&=\frac1{n_1^2}\left[p(1-p)\sum_i a_i^*b_i^*
-\frac{p(1-p)}{N-1}\sum_{i\ne j} a_i^*b_j^*\right]\\
&=\frac{p(1-p)}{n_1^2}\frac{N}{N-1}(N-1)S_{ab}
=\left(\frac1{n_1}-\frac1N\right)S_{ab}.
\end{aligned}
$$
令 $a_i=b_i=Y_i(1)$ 得处理组均值方差；交换两组，再令 $a_i=b_i=Y_i(0)$，得对照组均值方差 $(1/n_0-1/N)S_0^2$。

**两组之间的协方差。** 固定的未处理总量可写作 $N\bar Y(0)=n_1\bar Y_T(0)+n_0\bar Y_C(0)$。解出 $\bar Y_C(0)$，代入协方差，常数项消失：
$$
\operatorname{Cov}_Z(\bar Y_T(1),\bar Y_C(0))
=-\frac{n_1}{n_0}\left(\frac1{n_1}-\frac1N\right)S_{10}
=-\frac{S_{10}}N .
$$
这里 $S_{10}$ 是两个潜在结果的有限总体协方差。对均值差使用 $\operatorname{Var}(A-B)=\operatorname{Var}(A)+\operatorname{Var}(B)-2\operatorname{Cov}(A,B)$ 并合并 $1/N$ 项：
$$
\begin{aligned}
\operatorname{Var}_Z(\hat\tau)
&=\left(\frac1{n_1}-\frac1N\right)S_1^2+
\left(\frac1{n_0}-\frac1N\right)S_0^2+\frac{2S_{10}}N\\
&=\frac{S_1^2}{n_1}+\frac{S_0^2}{n_0}-\frac{S_1^2+S_0^2-2S_{10}}N .
\end{aligned}
$$
括号是 $\tau_i=Y_i(1)-Y_i(0)$ 的有限总体方差 $S_\tau^2$，由对 $\tau_i-\bar\tau$ 平方得到。对每臂 $n_d\ge2$，组内常规样本方差 $s_d^2$ 设计无偏：将其写成组内所有无序数对平方差之和除以 $n_d(n_d-1)$，每对入组的概率是 $n_d(n_d-1)/[N(N-1)]$，故期望等于全体数对平方差之和除以 $N(N-1)$，即 $S_d^2$。因此 $\hat V=s_1^2/n_1+s_0^2/n_0$ 的期望比真实方差大 $S_\tau^2/N$。恒定效应下差额为零。

[Neyman（1923，1990 年英译）](https://doi.org/10.1214/ss/1177012031)

## 3.3 分层、集群与协变量调整
**分层分配。** 令 $N_h$ 为层规模，$w_h=N_h/N$。将 $\tau_S=N^{-1}\sum_i\tau_i$ 按层重排，得到 $\sum_h(N_h/N)(\sum_{i\in h}\tau_i/N_h)=\sum_h w_h\tau_h$。在每层内应用第 3.1 节的证明，$E_Z[\hat\tau_h]=\tau_h$。由期望线性性及各层独立分配，
$$
E_Z[\sum_h w_h\hat\tau_h]=\sum_h w_h\tau_h,\qquad
\operatorname{Var}_Z(\sum_h w_h\hat\tau_h)=\sum_h w_h^2\operatorname{Var}_Z(\hat\tau_h).
$$
第二式展开时所有跨层协方差为零。直接合并处理组时，层均值的系数是 $n_{1h}/n_1$，合并对照组时则是 $n_{0h}/n_0$。两臂采用相同人口权重，才对应上式的目标。

**集群分配。** 社区 $c$ 有 $N_c$ 人，中签概率为 $0<p_c<1$。记 $T_c(z)=\sum_{i\in c}Y_i(z)$，假定无跨社区干扰。实际社区总量 $T_c^{obs}=Z_cT_c(1)+(1-Z_c)T_c(0)$。于是
$$
E_Z[Z_cT_c^{obs}/p_c]=T_c(1)E_Z[Z_c]/p_c=T_c(1),
\quad E_Z[(1-Z_c)T_c^{obs}/(1-p_c)]=T_c(0).
$$
相减、对社区求和、除以总人数 $N$，得设计无偏的个人等权效应。令 $\tau_c=[T_c(1)-T_c(0)]/N_c$，这个目标是 $\sum_cN_c\tau_c/N$。若直接比较等权社区均值，目标是 $\sum_c\tau_c/C$。两个社区人数为 10 与 90，效应分别为 10 与 0 时，两种目标各为 $(10\times10+90\times0)/100=1$ 和 $(10+0)/2=5$。

**协变量调整。** 设 $X$ 为基线协变量，分别在每个试验臂估计回归 $Y_i=a_z+(X_i-\bar X)'b_z+u_i$。截距一阶条件 $\sum_{Z_i=z}\hat u_i=0$ 意味着
$$
\hat a_z=\bar Y_z-(\bar X_z-\bar X)'\hat b_z,\qquad
\hat a_1-\hat a_0=\bar Y_1-\bar Y_0-(\bar X_1-\bar X)'\hat b_1+(\bar X_0-\bar X)'\hat b_0 .
$$
这是在同一协变量均值 $\bar X$ 处比较两臂的预测值。合并回归取截距 $a_0$、处理系数 $a_1-a_0$、基线斜率 $b_0$、交互斜率 $b_1-b_0$，代入 $Z=0,1$ 即得两条拟合直线。Lin 的渐近精度结果针对这类完全交互调整；有限样本中协变量失衡与估计斜率之积为随机量。

[Su 与 Ding（2021）](https://doi.org/10.1111/rssb.12468) · [Lin（2013）](https://doi.org/10.1214/12-AOAS583)

## 3.4 不完全依从、失访与干扰
**分配与实际处理。** 令 $D(z)$ 为给定分配下的实际处理，$Y^Z(z)$ 为分配潜在结果。随机分配与一致性给出 $E[Y\mid Z=z]=E[Y^Z(z)\mid Z=z]=E[Y^Z(z)]$。两臂相减得到 $ITT=E[Y^Z(1)-Y^Z(0)]$。若分配只通过实际处理影响结果且 $D$ 为二元，写 $Y^Z(z)=Y(0)+D(z)[Y(1)-Y(0)]$，则
$$
ITT=E[(D(1)-D(0))(Y(1)-Y(0))].
$$
当 $D(1)\ge D(0)$ 时，括号 $D(1)-D(0)$ 等于依从者指示量 $\mathbf1_C$。因此 $ITT=P(C)E[Y(1)-Y(0)\mid C]$，同时第一阶段 $E[D(1)-D(0)]=P(C)$。比值的分母需为正。

**失访。** 令 $R(z)$ 为结果可观测指示，$r_z=P(R(z)=1)$，$m_z=E[Y^Z(z)\mid R(z)=1]$。对失访者记 $q_z=E[Y^Z(z)\mid R(z)=0]$。由全期望公式 $\mu_z=r_zm_z+(1-r_z)q_z$。若结果在 $[L,U]$，因 $1-r_z\ge0$，替换 $q_z$ 为 $L,U$ 分别得到 $\ell_z=r_zm_z+(1-r_z)L$、$u_z=r_zm_z+(1-r_z)U$。从 $\mu_1\ge\ell_1$、$\mu_0\le u_0$ 得 $\mu_1-\mu_0\ge\ell_1-u_0$，反向取端点得上界 $u_1-\ell_0$。例如二元结果下，若 $(r_1,m_1)=(0.8,0.5)$、$(r_0,m_0)=(0.9,0.4)$，两臂均值区间为 $[0.40,0.60]$、$[0.36,0.46]$，故 ITT 在 $[-0.06,0.24]$。完整观测均值差 $0.10$ 只比较各臂的留存者。

**干扰与随机化检验。** 令邻域暴露为 $s$，结果为 $Y_i(z,s)$。相同 $s$ 下比较 $Y_i(1,s)-Y_i(0,s)$ 得直接效应；相同 $z$ 下比较 $Y_i(z,s_1)-Y_i(z,s_0)$ 得溢出效应。两阶段随机分配可产生这些暴露单元。在所有允许分配的潜在结果均相同这一锐零假设下，任一分配 $a$ 下的结果都等于观测向量 $Y^{obs}$。若 $p(a)$ 是已知分配概率，双侧随机化检验按
$$
p_F=\sum_a p(a)\mathbf1\{|T(a,Y^{obs})|\ge |T(a^{obs},Y^{obs})|\}
$$
累加至少同样极端的分配概率。等概率设计下，该式化为极端名单数除以所有允许名单数。

## 3.5 俄勒冈 Medicaid 抽签的分配效应
Finkelstein 等（2012）研究俄勒冈 Medicaid 名额抽签。中签为家庭提供申请机会，参保另需提交申请并通过资格审查；多人登记的家庭获得机会的概率较高。原文控制登记家庭人数，并按家庭聚类估计标准误。

设 $H=h$ 表示登记人数层，$Z$ 表示中签，$D$ 为实际参保，$Y$ 为就医或财务结果。定义 $\Delta_V(h)=E[V\mid Z=1,h]-E[V\mid Z=0,h]$。在层内随机分配、排除限制和单调性下，将第 3.4 节的公式分别用于各层：
$$
\Delta_D(h)=P(C\mid h)=p_{Ch},\qquad
\Delta_Y(h)=p_{Ch}\tau_{Ch}.
$$
按名单中目标个体的层分布 $\omega_h=P(H=h)$ 加权，可得标准化申请机会效应 $\sum_h\omega_h\Delta_Y(h)$。再除以同权重的参保第一阶段：
$$
\frac{\sum_h\omega_h\Delta_Y(h)}{\sum_h\omega_h\Delta_D(h)}
=\sum_h\frac{\omega_hp_{Ch}}{\sum_j\omega_jp_{Cj}}\tau_{Ch}.
$$
由贝叶斯公式，最后一个和式中的权重为 $P(H=h\mid C)$；因此比率对因中签而参保的人按其实际层分布加权。

**排除限制的数值例。** 设标准化的参保第一阶段为 $0.25$，就医概率的简约式为 $0.05$。若申请机会仅通过参保改变就医，$0.05=0.25\tau_C$，故 $\tau_C=0.20$。若申请过程另使总体就医概率增加 $0.01$，同一简约式变为 $0.05=0.01+0.25\tau_C$，解得 $\tau_C=0.16$。这里的数字用于比较两种分配机制；随机化本身识别的申请机会效应都是 $0.05$。

[原文与第一年结果](https://economics.mit.edu/sites/default/files/2022-08/The%20Oregon%20Health%20Insurance%20Experiment%20%20Evidence%20f.pdf)
