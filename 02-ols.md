# 从 OLS 到因果识别

设 $D\in\{0,1\}$ 为处理状态，$Y$ 为观测结果。二元处理的总体回归系数等于处理组与对照组的结果均值差。其因果解释取决于两组在同一处理状态下的潜在结果分布。Angrist 与 Pischke 以理想随机试验为基准，分析观察研究中得到可比反事实的条件。 [教材](https://press.princeton.edu/books/paperback/9780691120355/mostly-harmless-econometrics)

## 2.1 条件期望与总体线性投影
**平方损失与条件均值。** 设 $E[Y^2]<\infty$，预测函数 $g(X)$ 平方可积，记 $m(X)=E[Y\mid X]$。对同一个观测加减 $m(X)$，令 $A=Y-m(X)$、$B=m(X)-g(X)$。由于 $B$ 由 $X$ 决定，
$$
E[AB\mid X]=B\,E[Y-m(X)\mid X]=B(m(X)-m(X))=0 .
$$
对 $Y-g(X)=A+B$ 平方，再取期望，用 $E[AB]=E[E[AB\mid X]]=0$ 消去交叉项，
$$
E[(Y-g(X))^2]=E[A^2]+2E[AB]+E[B^2]
=E[(Y-m(X))^2]+E[(m(X)-g(X))^2].
$$
第一项与 $g$ 无关。第二项非负，且仅当 $g(X)=m(X)$ 几乎处处成立时等于零。这给出平方损失意义下的最优预测函数。

**将函数空间限制为线性空间。** 设 $X$ 为含常数项的列向量，分量均有二阶矩，$Q=E[XX']$ 正定。对 $g(X)=X'b$ 展开总体目标：
$$
L(b)=E[Y^2]-2b'E[XY]+b'Qb .
$$
对第 $j$ 个坐标求导。线性项贡献 $-2E[X_jY]$。二次项为 $\sum_{k,l}b_kQ_{kl}b_l$，分别对左、右两个系数求导得 $\sum_l Q_{jl}b_l+\sum_k b_kQ_{kj}=2(Qb)_j$，最后一步用 $Q'=Q$。合并所有坐标，
$$
\nabla_b L=-2E[XY]+2Qb=0
\quad\Longrightarrow\quad Qb=E[XY]
\quad\Longrightarrow\quad \beta^{LP}=Q^{-1}E[XY].
$$
Hessian 为 $2Q$，正定性保证解唯一。令 $u=Y-X'\beta^{LP}$，则 $E[Xu]=E[XY]-Q\beta^{LP}=0$，而 $E[u\mid X]=m(X)-X'\beta^{LP}$。后者几乎处处为零的充要条件是条件均值属于给定线性空间。

## 2.2 样本最小二乘与二元回归
**先求截距，再求斜率。** 在样本中以 $a+bD_i$ 拟合 $Y_i$。记 $r_i=Y_i-a-bD_i$。对 $Q(a,b)=\sum_i r_i^2$ 使用链式法则：
$$
\frac{\partial Q}{\partial a}=\sum_i 2r_i(-1),\qquad
\frac{\partial Q}{\partial b}=\sum_i 2r_i(-D_i).
$$
第一式等于零给出 $\sum_iY_i-na-b\sum_iD_i=0$。除以 $n$ 得 $a=\bar Y-b\bar D$。将其代入残差，有 $r_i=(Y_i-\bar Y)-b(D_i-\bar D)$。第二式等于零是 $\sum_iD_i r_i=0$。减去 $\bar D\sum_i r_i=0$ 后依次得到
$$
\begin{aligned}
0&=\sum_i(D_i-\bar D)r_i\\
 &=\sum_i(D_i-\bar D)(Y_i-\bar Y)-b\sum_i(D_i-\bar D)^2,\\
\hat b&=\frac{\sum_i(D_i-\bar D)(Y_i-\bar Y)}{\sum_i(D_i-\bar D)^2},
\qquad \hat a=\bar Y-\hat b\bar D .
\end{aligned}
$$
最后一步要求分母为正。若所有人处理状态相同，设计矩阵的两列线性相关，截距与斜率无法分别确定。

**二元处理。** 令 $n_d$ 为 $D_i=d$ 的人数，假定 $n_0,n_1>0$。对照组拟合值 $m_0=a$，处理组拟合值 $m_1=a+b$。从 $(m_0,m_1)$ 到 $(a,b)$ 的逆变换为 $a=m_0,\ b=m_1-m_0$，所以可分别最小化两组平方和。对第 $d$ 组，
$$
0=-2\sum_{D_i=d}(Y_i-m_d)
\iff n_dm_d=\sum_{D_i=d}Y_i
\iff \hat m_d=\bar Y_d .
$$
代回逆变换得 $\hat a=\bar Y_0$、$\hat b=\bar Y_1-\bar Y_0$。总体对应式用 $E[(Y-m_d)^2\mid D=d]$ 取代组内平均平方损失，其最优常数为 $E[Y\mid D=d]$，因此总体系数同样等于两组条件均值之差。

## 2.3 潜在结果与选择偏误
**从观测差到反事实差。** 令 $Y(1),Y(0)$ 为同一人在参训与未参训状态下的收入。处理定义一致、单位间无干扰时，$Y=DY(1)+(1-D)Y(0)$。在 $D=1$ 的样本中第二项为零，在 $D=0$ 的样本中第一项为零。因此二元回归的总体系数可按下列次序改写：
$$
\begin{aligned}
b&=E[Y(1)\mid D=1]-E[Y(0)\mid D=0]\\
 &=E[Y(1)\mid D=1]-E[Y(0)\mid D=1]\\
 &\quad+E[Y(0)\mid D=1]-E[Y(0)\mid D=0]\\
 &=\underbrace{E[Y(1)-Y(0)\mid D=1]}_{ATT}
   +\underbrace{E[Y(0)\mid D=1]-E[Y(0)\mid D=0]}_{B_0}.
\end{aligned}
$$
第二行加减的是参训者若未参训的收入均值；第三行用同一条件下期望的线性性。$b=ATT$ 当且仅当 $B_0=0$。设处理组实际均值为 22，对照组实际均值为 25，处理组反事实均值为 18。观测差为 $22-25=-3$，而 $ATT=22-18=4$。其差来自 $B_0=18-25=-7$。

**从 ATT 到 ATE。** 写 $p=P(D=1)$、$\tau=Y(1)-Y(0)$、$ATU=E[\tau\mid D=0]$。由全期望公式，
$$
ATE=pATT+(1-p)ATU,\qquad ATT-ATE=(1-p)(ATT-ATU).
$$
将 $b-ATE$ 写作 $(b-ATT)+(ATT-ATE)$，便得到 $B_0+(1-p)(ATT-ATU)$。第一项描述未处理水平上的选择，第二项描述处理收益上的选择。若 $D$ 与潜在结果对独立，则 $E[Y(0)\mid D=d]=E[Y(0)]$ 且 $E[\tau\mid D=d]=E[\tau]$，两项同时为零。

## 2.4 条件可忽略性、重叠与标准化
**标准化。** 令 $X$ 为处理前协变量，$e(x)=P(D=1\mid X=x)$，$m_d(x)=E[Y\mid D=d,X=x]$。假定一致性、$(Y(0),Y(1))\perp D\mid X$ 和 $0<e(X)<1$ 几乎处处成立。条件独立给出 $E[Y(d)\mid X=x]=E[Y(d)\mid D=d,X=x]$；一致性使右边等于 $m_d(x)$；重叠使两个条件均值均可从观测分布确定。逐臂积分再相减：
$$
\begin{aligned}
E[Y(d)]&=\int E[Y(d)\mid X=x]\,dF_X(x)=\int m_d(x)\,dF_X(x),\\
ATE&=\int [m_1(x)-m_0(x)]\,dF_X(x).
\end{aligned}
$$
对 ATT，处理组收入已观测，需补出的量是 $E[Y(0)\mid D=1]$。先在处理组内按 $X$ 分层，再用对照组估计反事实：
$$
E[Y(0)\mid D=1]=\int E[Y(0)\mid D=1,X=x]\,dF_{X\mid D=1}(x)
=\int m_0(x)\,dF_{X\mid D=1}(x).
$$
第二个等号只需未处理潜在结果的条件均值可比，以及处理组的支持被对照组覆盖。故 ATT 的积分权重为处理组的 $X$ 分布。

**逆概率加权。** 固定 $X=x$，按 $D=1,0$ 两种取值展开：
$$
\begin{aligned}
E[DY\mid x]&=e(x)\cdot 1\cdot m_1(x)+(1-e(x))\cdot0=e(x)m_1(x),\\
E[(1-D)Y\mid x]&=e(x)\cdot0+(1-e(x))\cdot1\cdot m_0(x).
\end{aligned}
$$
分别除以各自非零的入组概率，得到 $m_1(x),m_0(x)$。减去第二项、对 $X$ 取期望，并应用上一段的标准化等式：
$$
ATE=E\left[\frac{DY}{e(X)}-\frac{(1-D)Y}{1-e(X)}\right].
$$

**增广估计式的双重稳健性。** 令 $a_d(x)$ 为候选结果均值，$0<q(x)<1$ 为候选处理概率。省略 $x$ 记号，定义 $\psi_1=a_1+D(Y-a_1)/q$。使用上面已推得的 $E[DY\mid x]=em_1$ 与 $E[D\mid x]=e$：
$$
\begin{aligned}
E[\psi_1\mid x]-m_1
&=a_1-m_1+(em_1-ea_1)/q\\
&=(a_1-m_1)-(e/q)(a_1-m_1)=(a_1-m_1)(1-e/q).
\end{aligned}
$$
对照臂令 $\psi_0=a_0+(1-D)(Y-a_0)/(1-q)$，以 $1-e$ 代入同一展开，得 $(a_0-m_0)[1-(1-e)/(1-q)]$。若两臂结果模型均正确，或处理概率模型正确，则两臂的条件偏差均为零。取期望并作差得 $E[\psi_1-\psi_0]=ATE$。这是总体矩恒等式；将函数替换为拟合值时，估计误差由两个误差因子的乘积进入余项。

[Rosenbaum 与 Rubin（1983）](https://doi.org/10.1093/biomet/70.1.41)

## 2.5 LaLonde 职业培训研究中的外部对照组
LaLonde（1986）以 NSW 职业培训随机试验为基准，保留项目处理组，以 PSID 或 CPS–SSA 中的非实验样本替换实验对照组。原文表 5 对男性参加者比较多组估计。设 $T$ 是项目处理组，$R$ 是实验对照组，$C$ 是外部对照组；$F_T$ 表示 $T$ 中培训前协变量的分布。

**比较组替换引入的偏差。** 固定项目处理组，目标和替换后的总体比较量分别是
$$
\tau_T=E[Y(1)\mid T]-E[Y(0)\mid T],\qquad
\theta_C=E[Y\mid T]-\int E[Y\mid X=x,C]\,dF_T(x).
$$
在 $T$ 观察到 $Y(1)$，在 $C$ 观察到 $Y(0)$。两式相减，处理组实际收入抵消：
$$
\begin{aligned}
\theta_C-\tau_T
&= E[Y(0)\mid T]-\int E[Y(0)\mid x,C]\,dF_T(x)\\
&=\int\{E[Y(0)\mid x,T]-E[Y(0)\mid x,C]\}\,dF_T(x).
\end{aligned}
$$
最后一步对第一项在 $T$ 内使用全期望公式。重加权只替换 $X$ 的分布。若项目招收者在培训前遭遇暂时性收入下降，即使调整相同的收入水平，两组此后的自然恢复速度仍可能不同；该差异进入被积函数。培训前收入的时间序列因而能检验外部对照组是否重现处理前的收入路径。

**与试验基准的精度比较。** 简单均值规格给出 $\hat\theta_C=\bar Y_T-\bar Y_C$、$\hat\tau_R=\bar Y_T-\bar Y_R$。两者相减等于 $\bar Y_R-\bar Y_C$，共享的处理组均值完全抵消。一般调整规格下，
$$
\operatorname{Var}(\hat\theta_C-\hat\tau_R)
=\operatorname{Var}(\hat\theta_C)+\operatorname{Var}(\hat\tau_R)-2\operatorname{Cov}(\hat\theta_C,\hat\tau_R).
$$
这一区分决定比较两种方法时采用哪一项标准误。

[原论文与发表信息](https://ideas.repec.org/a/aea/aecrev/v76y1986i4p604-20.html) · [论文表 5](https://business.baylor.edu/scott_cunningham/teaching/lalonde-1986.pdf)
