# 断点回归 RDD

设资格由运行变量 $X$ 跨越阈值 $c$ 决定。断点回归的总体参数由条件均值在 $c$ 处的两个单侧极限定义。清晰断点下，跨越阈值使处理状态从 0 变为 1；模糊断点下，跨越阈值引起处理概率的非零跳跃。

## 4.1 清晰断点的局部处理效应
**门槛两侧分别观察到什么。** 设 $D=\mathbf1\{X\ge c\}$，潜在结果为 $Y(d)$，定义 $m_d(x)=E[Y(d)\mid X=x]$。由一致性 $Y=DY(1)+(1-D)Y(0)$。给定 $X=x$，$D$ 是已知常数，可以提出条件期望：
$$
m(x)=E[Y\mid X=x]=\mathbf1\{x\ge c\}m_1(x)+\mathbf1\{x<c\}m_0(x).
$$
对任何 $x>c$ 有 $m(x)=m_1(x)$，对任何 $x<c$ 有 $m(x)=m_0(x)$。假定两侧均有任意接近门槛的支持，且潜在均值在 $c$ 连续。分别沿两侧取极限：
$$
\begin{aligned}
m_+(c)&=\lim_{x\downarrow c}m(x)=\lim_{x\downarrow c}m_1(x)=m_1(c),\\
m_-(c)&=\lim_{x\uparrow c}m(x)=\lim_{x\uparrow c}m_0(x)=m_0(c),\\
m_+(c)-m_-(c)&=E[Y(1)\mid X=c]-E[Y(0)\mid X=c]=E[Y(1)-Y(0)\mid X=c].
\end{aligned}
$$
连续运行变量下，$P(X=c)$ 通常为零；这里的条件均值由连续延拓确定。该参数对门槛处的人群取平均。

**连续性偏离与捆绑政策。** 若未处理均值自身在门槛处跳跃，记 $m_d^\pm$ 为单侧极限。在观测跳跃中加减右侧未处理均值：
$$
m_1^+-m_0^-=(m_1^+-m_0^+)+(m_0^+-m_0^-).
$$
第一项是右侧极限人群的处理效应，第二项是未处理结果的结构性跳跃。若两项政策同时启动，观测到的是 $Y(1,1)$ 与 $Y(0,0)$。加减 $Y(0,1)$ 后，联合效应分解为 $E[Y(1,1)-Y(0,1)\mid c]+E[Y(0,1)-Y(0,0)\mid c]$。其中 $Y(0,1)$ 对应的制度组合未在断点两侧出现。

[Hahn、Todd 与 van der Klaauw（2001）](https://doi.org/10.1111/1468-0262.00183)

## 4.2 局部线性估计、带宽与偏差校正
**求出断点两侧的截距。** 设 $u_i=X_i-c$，带宽 $h>0$，单侧样本权重 $k_i=K(u_i/h)\ge0$。固定一侧，在该侧取和，最小化 $Q(a,b)=\sum_i k_i(Y_i-a-bu_i)^2$。分别求偏导、除以 $-2$ 并移项：
$$
\sum_i k_iY_i=a\sum_i k_i+b\sum_i k_i u_i,\qquad
\sum_i k_i u_iY_i=a\sum_i k_i u_i+b\sum_i k_i u_i^2.
$$
记 $S_j=\sum_i k_i u_i^j$、$T_j=\sum_i k_i u_i^jY_i$。将第一式乘 $S_2$、第二式乘 $S_1$ 后相减，消去斜率；将第二式乘 $S_0$、第一式乘 $S_1$ 后相减，消去截距：
$$
\hat a=\frac{S_2T_0-S_1T_1}{S_0S_2-S_1^2},\qquad
\hat b=\frac{S_0T_1-S_1T_0}{S_0S_2-S_1^2}.
$$
分母为 $S_0\sum_i k_i(u_i-S_1/S_0)^2$，故一侧至少需要有正权重且相异的运行变量取值。门槛处 $u=0$，单侧预测是 $\hat a$，跳跃估计为 $\hat a_+-\hat a_-$。

**主导偏差的来源。** 展开 $T_j$ 的定义，将截距写为 $\hat a=\sum_i\ell_iY_i$，其中 $\ell_i=k_i(S_2-S_1u_i)/(S_0S_2-S_1^2)$。分别对权重及其一次矩求和：
$$
\sum_i\ell_i=\frac{S_2S_0-S_1^2}{S_0S_2-S_1^2}=1,\qquad
\sum_i\ell_i u_i=\frac{S_2S_1-S_1S_2}{S_0S_2-S_1^2}=0 .
$$
假定均值函数在该侧有连续二阶导数。给定运行变量，代入 $m(c+u_i)=m(c)+m'(c)u_i+m''(c)u_i^2/2+o(u_i^2)$。常数项原样保留，一次项消失，于是
$$
E[\hat a\mid X]-m(c)=\tfrac12 m''(c)\sum_i\ell_i u_i^2+\sum_i\ell_i o(u_i^2).
$$
在有界支持核、门槛处正且连续的单侧密度和非奇异核矩条件下，窗口内 $u_i=O(h)$，有效样本量阶数为 $nh$，单个权重为 $O_p(1/(nh))$。故二阶矩总和为 $O_p(h^2)$。独立观测且条件方差在门槛处有限、连续时，$\operatorname{Var}(\hat a\mid X)=\sum_i\ell_i^2\sigma_i^2$ 为 $1/(nh)$ 阶。在 $h\to0$、$nh\to\infty$ 下，合并两侧得到跳跃估计的主导偏差 $Bh^2$ 和主导方差 $V/(nh)$。

具体看权重的阶数，令 $v_i=u_i/h$、$L_j=(nh)^{-1}\sum_i k_i v_i^j$，则 $S_j=nh^{j+1}L_j$。若 $L_j$ 有有限极限且 $L_0L_2-L_1^2$ 的极限为正，把三个 $S_j$ 代回截距的权重式：
$$
\ell_i=\frac{k_i}{nh}\frac{L_2-L_1v_i}{L_0L_2-L_1^2}.
$$
紧支持核使窗口内 $v_i$ 有界。因此 $O_p(nh)$ 个非零权重各为 $O_p(1/(nh))$，其平方和为 $O_p(nh)\times O_p((nh)^{-2})=O_p((nh)^{-1})$。同理，$\sum_i\ell_i u_i^2=h^2\sum_i\ell_i v_i^2$，最后的加权和有有限极限，从而出现二阶偏差。

**带宽与区间。** 若 $B\ne0$、$V>0$，主导均方误差为 $M(h)=B^2h^4+V/(nh)$。求导并令其为零：
$$
4B^2h^3-\frac{V}{nh^2}=0
\iff 4nB^2h^5=V
\iff h_*=\left(\frac{V}{4nB^2}\right)^{1/5}.
$$
二阶导数 $12B^2h^2+2V/(nh^3)>0$。代回得偏差阶数 $h_*^2\asymp n^{-2/5}$，标准误阶数 $(nh_*)^{-1/2}\asymp n^{-2/5}$，故以零偏差为中心的正态近似保留了一阶中心偏移。稳健偏差校正以局部二次拟合估计曲率，得到偏差估计 $\hat b_h$。对 $\hat\tau_{bc}=\hat\tau-\hat b_h$ 应展开
$$
\operatorname{Var}(\hat\tau_{bc})=\operatorname{Var}(\hat\tau)+\operatorname{Var}(\hat b_h)-2\operatorname{Cov}(\hat\tau,\hat b_h).
$$
两次拟合使用同一批观测，最后的协方差由两组等价权重的乘积进入方差估计。

[Calonico、Cattaneo 与 Titiunik（2014）](https://doi.org/10.3982/ECTA11757)

## 4.3 模糊断点、离散评分与多门槛
**模糊断点。** 将阈值资格记作 $Z$，实际处理记作 $D$。令 $D(z)$ 为在资格 $z$ 下会接受的处理，$\tau_i=Y_i(1)-Y_i(0)$。排除限制规定资格仅通过处理影响结果。于是潜在结果在资格 $z$ 下为
$$
Y_i(D_i(z))=D_i(z)Y_i(1)+(1-D_i(z))Y_i(0)=Y_i(0)+D_i(z)\tau_i .
$$
记 $\Delta_V$ 为 $E[V\mid X=x]$ 的右极限减左极限。假定 $E[D(z)\mid x]$ 与 $E[Y(D(z))\mid x]$ 在阈值连续，则两侧极限可分别替换为 $c$ 处两种资格状态的潜在均值。对 $D$ 与 $Y$ 各做一次减法：
$$
\begin{aligned}
\Delta_D &= E[D(1)-D(0)\mid c],\\
\Delta_Y &= E[(Y(0)+D(1)\tau_i)-(Y(0)+D(0)\tau_i)\mid c]\\
         &= E[(D(1)-D(0))\tau_i\mid c].
\end{aligned}
$$
在单调性下，$D(1)-D(0)=\mathbf1_C$，其中 $C$ 为资格改变实际处理的个体。故 $\Delta_D=P(C\mid c)$，$\Delta_Y=P(C\mid c)E[\tau_i\mid C,c]$。若 $\Delta_D>0$，相除得局部依从者平均效应。若右侧资格对结果还有加性直接效应 $\delta$，则 $\Delta_Y=\delta+\Delta_D\tau_C$，观测比率比 $\tau_C$ 多出 $\delta/\Delta_D$。

**离散评分。** 设最近的左右支持点为 $x_L<c<x_R$。相邻均值差可以写成
$$
m_1(x_R)-m_0(x_L)=\tau(c)+[m_1(x_R)-m_1(c)]+[m_0(c)-m_0(x_L)].
$$
若两侧满足 $|m'_d(x)|\le L_d$，由积分公式，第一括号绝对值至多为 $L_1(x_R-c)$，第二括号至多为 $L_0(c-x_L)$。故相邻差偏离阈值效应的绝对值不超过两上界之和。对固定的整数评分，增加同一评分上的人数可改善均值精度，支持点之间的距离仍然存在。局部随机化则指定一个有限窗口和窗口内的分配机制，在该窗口定义参数。多个阈值 $c_j$ 分别识别 $\tau(c_j)$ 时，给定目标权重 $w_j$，总体目标为 $\sum_j w_j\tau(c_j)$，其中 $\sum_j w_j=1$。

[Cattaneo、Idrobo 与 Titiunik（2024）](https://www.cambridge.org/core/elements/practical-introduction-to-regression-discontinuity-designs/C6A70A32359115510AAC370A7869AE2F)

## 4.4 Medicare 资格门槛与急诊住院样本
Card、Dobkin 与 Maestas（2009）利用美国 65 岁 Medicare 资格门槛研究急诊住院患者的治疗与死亡率。样本由已入院者构成。设 $A(z)$ 为资格状态 $z$ 下是否入院，$Y(z)$ 为同一源人群在该资格下的结果。以下均以 $X=c$ 为条件。

**入院选择如何进入断点。** 设 $P(A(z)=1\mid c)>0$，且分子、分母所需潜在矩均在门槛处连续。对二元 $A$，有 $E[Y\mid A=1,x]=E[AY\mid x]/P(A=1\mid x)$。分别沿两侧取极限，由一致性将分子换为 $E[A(z)Y(z)\mid c]$，分母换为 $P(A(z)=1\mid c)$。再用同一条件期望恒等式，得到住院样本的结果跳跃
$$
\theta_A=E[Y(1)\mid A(1)=1,c]-E[Y(0)\mid A(0)=1,c].
$$
加减 $E[Y(0)\mid A(1)=1,c]$，把同一群体的结果差与两个群体的构成差分开：
$$
\begin{aligned}
\theta_A&=E[Y(1)-Y(0)\mid A(1)=1,c]\\
&\quad+E[Y(0)\mid A(1)=1,c]-E[Y(0)\mid A(0)=1,c].
\end{aligned}
$$
若 $A(1)=A(0)$ 几乎处处成立，第二行两个条件事件相同，故整行为零。仅有入院概率相等时，设两侧入院率都是 $q$，第二行可重写为 $E[(A(1)-A(0))Y(0)\mid c]/q$；即使 $E[A(1)-A(0)\mid c]=0$，该结果加权差仍可非零。原文聚焦难以推迟的急症，并比较门槛附近入院频率和病情特征，为入院构成的稳定性提供证据。

**从资格效应换算为参保效应。** 在人群构成稳定的数值例中，若结果跳跃为 4，参保概率跳跃为 0.4，单调性和排除限制给出 $4=0.4\tau_C$，故 $\tau_C=10$。若另有大小为 1 的资格直接效应，则 $4=1+0.4\tau_C$，解得 $\tau_C=7.5$。资格跳跃识别的是资格的总效应；换算实际参保效应还用到了第一阶段和直接渠道的限制。

[原论文](https://doi.org/10.1162/qjec.2009.124.2.597)
