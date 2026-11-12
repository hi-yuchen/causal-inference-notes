# 控制变量、固定效应与标准误

协变量和固定效应共同决定处理变量的残差及其在投影系数中的权重。给定回归规格，系数的抽样方差取决于残差化处理与误差乘积的协方差结构。以下先推导点估计，再讨论方差估计。

## 7.1 遗漏变量偏误与 Frisch–Waugh–Lovell 定理
**短回归的总体系数。** 设 $Y=\alpha+\tau D+\gamma X+u$ 且 $E[u\mid D,X]=0$。由全期望公式，$E[u]=0$ 且 $E[Du]=E[D E[u\mid D,X]]=0$，因此 $\operatorname{Cov}(D,u)=0$。对 $D$ 与结果取协方差，常数项消失：
$$
\begin{aligned}
\operatorname{Cov}(D,Y)&=\tau\operatorname{Var}(D)+\gamma\operatorname{Cov}(D,X)+\operatorname{Cov}(D,u),\\
b_{\mathrm{short}}&=\frac{\operatorname{Cov}(D,Y)}{\operatorname{Var}(D)}
=\tau+\gamma\frac{\operatorname{Cov}(D,X)}{\operatorname{Var}(D)} .
\end{aligned}
$$
第二行使用含截距一元回归的斜率公式，要求 $\operatorname{Var}(D)>0$。若 $D$ 为二元，记 $p=P(D=1)$。第 6.1 节的二元协方差恒等式给出分子 $p(1-p)(E[X\mid1]-E[X\mid0])$，分母为 $p(1-p)$，相除即为两组基线差。若基线差为 2、$\gamma=3$，省略 $X$ 使系数增加 $3\times2=6$。

**FWL 定理。** 以下 $X$ 表示含截距的满列秩样本矩阵，$D$ 为处理列。给定候选斜率 $t$，先拟合控制系数 $b$。正规方程及解为
$$
X'(Y-Dt-Xb)=0 \ \Longrightarrow\ X'X b=X'(Y-Dt)
\ \Longrightarrow\ b(t)=(X'X)^{-1}X'(Y-Dt).
$$
令 $P_X=X(X'X)^{-1}X'$、$M=I-P_X$，则拟合后的残差为 $(Y-Dt)-Xb(t)=M(Y-Dt)$。由 $P_X^2=P_X=P_X'$，有 $M'=M$、$M^2=I-2P_X+P_X^2=M$。因此剖面平方和为
$$
\begin{aligned}
Q(t)&=(Y-Dt)'M'M(Y-Dt)\\
&=Y'MY-2tD'MY+t^2D'MD .
\end{aligned}
$$
对 $t$ 求导，$-2D'MY+2tD'MD=0$。只要 $D'MD>0$，二阶导数为正，解为
$$
\hat\tau=\frac{D'MY}{D'MD}=\frac{(MD)'(MY)}{(MD)'(MD)} .
$$
最后一个等号再次使用 $M'M=M$。因此先对 $D,Y$ 各自剔除控制变量，再以两个残差回归，可恢复长回归的处理系数。

**剩余处理变异与精度。** 因 $MX=0$，代入结构方程得 $\hat\tau-\tau=D'Mu/(D'MD)$。若 $\operatorname{Var}(u\mid D,X)=\sigma^2I$，给定设计矩阵，分母为常数。于是
$$
\operatorname{Var}(\hat\tau\mid D,X)
=\frac{D'M(\sigma^2 I)MD}{(D'MD)^2}
=\frac{\sigma^2 D'MD}{(D'MD)^2}
=\frac{\sigma^2}{D'MD}.
$$
加入控制可改变残差方差及剩余处理变异，两者分别出现在分子和分母。

## 7.2 混杂、中介和碰撞变量
**中介。** 设就业为 $M=aD+v$，收入为 $Y=bD+cM+u$，方程可作结构解释。按自然就业路径代入得 $Y(d)=(b+ac)d+cv+u$，所以 $Y(1)-Y(0)=b+ac$。若把就业固定在同一值 $m$，则 $Y(1,m)-Y(0,m)=(b+cm+u)-(cm+u)=b$。前者含经就业传导的 $ac$，后者是受控直接效应。要以包含 $D,M$ 的回归估计 $b$，还需相应外生性，例如 $E[u\mid D,M]=0$，以及设计矩阵满秩。

**碰撞变量。** 若 $D$ 与能力 $U$ 独立且方差均为正，$S=D+U$、$Y=U$，干预 $D$ 对 $Y$ 的效应为零。边际回归斜率为 $\operatorname{Cov}(D,U)/\operatorname{Var}(D)=0$。给定 $S$，由定义解出 $U=S-D$，故 $Y=1\cdot S-1\cdot D$。两斜率使拟合误差逐个观测等于零，因此是最小二乘解；$S$ 中含有独立、正方差的 $U$，所以解唯一。控制共同结果 $S$ 在原本独立的 $D,U$ 间引入了条件关联。

**控制分层截距后的权重。** 令 $p_x=P(X=x)$、$e_x=P(D=1\mid X=x)$、$m_d(x)=E[Y\mid D=d,x]$。在全套层截距上投影 $D$，拟合值是各层均值 $e_x$，故 $\tilde D=D-e_x$。在层 $x$ 内枚举 $D$ 的两个取值：
$$
\begin{aligned}
E[\tilde D Y\mid x]&=e_x(1-e_x)m_1(x)+(1-e_x)(-e_x)m_0(x)\\
&=e_x(1-e_x)(m_1(x)-m_0(x)),\\
E[\tilde D^2\mid x]&=e_x(1-e_x)^2+(1-e_x)e_x^2=e_x(1-e_x).
\end{aligned}
$$
在条件可交换性、一致性与重叠下，$m_1(x)-m_0(x)=\tau_x$ 为层内因果效应。令 $v_x=e_x(1-e_x)$，总体 FWL 公式给出 $b=\sum_x p_x v_x\tau_x/\sum_xp_xv_x$。总体 ATE 等于 $\sum_x p_x\tau_x$，相减并通分：
$$
b-ATE=\frac{E[v_X\tau_X]-E[v_X]E[\tau_X]}{E[v_X]}
=\frac{\operatorname{Cov}(v_X,\tau_X)}{E[v_X]} .
$$
因此即使所有混杂都已由分层控制，处理概率更接近二分之一的层仍会在共同斜率中占更大权重。

[Cinelli、Forney 与 Pearl](https://doi.org/10.1177/00491241221099552)

## 7.3 加性固定效应、类别交互与交互式因子模型
**个体固定效应。** 在 $Y_{it}=\alpha_i+\tau D_{it}+u_{it}$ 中，给定候选 $b$，对 $\alpha_i$ 最小化 $\sum_t(Y_{it}-\alpha_i-bD_{it})^2$。一阶条件是 $\sum_t(Y_{it}-\alpha_i-bD_{it})=0$，解为 $\alpha_i(b)=\bar Y_i-b\bar D_i$。代回后第 $it$ 项残差为 $(Y_{it}-\bar Y_i)-b(D_{it}-\bar D_i)$。再求 $b$ 的一阶条件并整理：
$$
\hat\tau=\frac{\sum_{it}(D_{it}-\bar D_i)(Y_{it}-\bar Y_i)}
{\sum_{it}(D_{it}-\bar D_i)^2}
=\tau+\frac{\sum_{it}(D_{it}-\bar D_i)(u_{it}-\bar u_i)}
{\sum_{it}(D_{it}-\bar D_i)^2}.
$$
最后一步代入去均值后的结构方程。若健康冲击同时改变个人参保和收入，该冲击进入分子的组内乘积。时间不变的参保状态则使相应个体对分母的贡献为零。

**双向固定效应。** 在平衡等权面板中构造 $r_{it}=D_{it}-\bar D_i-\bar D_t+\bar D$。对 $t$ 取均值，前两项相消，后两项为 $-T^{-1}\sum_t\bar D_t+\bar D=0$。对 $i$ 取均值同理为零。因此 $r$ 与所有单位、时期虚拟变量正交。同时 $D_{it}-r_{it}=\bar D_i+\bar D_t-\bar D$ 属于这两个虚拟变量空间之和。由正交投影的唯一性，$r$ 就是双向残差。不平衡面板中，各单位观测到不同时期，个体内的 $\bar D_t$ 平均未必等于全局均值，故需对实际虚拟变量矩阵联合残差化。

**类别交互固定效应。** 设城市试点指标为 $A_{ct}$、个人资格为 $E_{ict}$。在城市×年份单元中，$A_{ct}$ 是常数，其单元均值等于自身，处理残差为零。若实际处理是 $D_{ict}=A_{ct}E_{ict}$，则单元均值为 $A_{ct}\bar E_{ct}$，相减得 $\tilde D_{ict}=A_{ct}(E_{ict}-\bar E_{ct})$。若再控制资格主效应，剔除单元均值后资格变为 $\tilde E=E-\bar E_{ct}$。对它再投影，斜率为 $\kappa=\sum A_{ct}\tilde E^2/\sum\tilde E^2$，最终处理残差为 $(A_{ct}-\kappa)\tilde E$。这要求试点与非试点单元均有足够的资格变异才能识别交互效应。

**交互式因子模型。** Bai 的设定为 $Y_{it}=X_{it}'\beta+\lambda_i'f_t+u_{it}$，共同冲击 $f_t$ 对不同单位有不同载荷。两地区、两时期的共同分量作双重差分得到
$$
(\lambda_A'f_1-\lambda_A'f_0)-(\lambda_B'f_1-\lambda_B'f_0)
=(\lambda_A-\lambda_B)'(f_1-f_0).
$$
若载荷不同且因子变化，该项仍可非零。将残差排成 $N\times T$ 矩阵 $R(\beta)$。给定斜率和 $r$ 个因子，目标为 $\min_{\Lambda,F}\|R-\Lambda F'\|_F^2$。采用 $F'F=I_r$ 的归一化，对 $\Lambda$ 的正规方程为 $\Lambda F'F=RF$，得 $\hat\Lambda=RF$。代回拟合为 $RFF'$，目标化为 $\operatorname{tr}(R'R)-\operatorname{tr}(F'R'RF)$。因此选择 $R'R$ 最大的 $r$ 个特征值所对应的正交特征向量使目标最小；给定因子分量，减去它后再用最小二乘更新 $\beta$。因子与载荷存在旋转不确定性，因为 $(\Lambda H)(FH^{-T})'=\Lambda HH^{-1}F'=\Lambda F'$。斜率识别还要求解释变量含有可与低秩成分分离的变异。

为写出因子更新的极值条件，令 $R'R=U\operatorname{diag}(s_1^2,\ldots,s_T^2)U'$，特征值由大到小排列，令 $V=U'F$。则 $V'V=I_r$。记 $a_j=\sum_{k=1}^r V_{jk}^2$，它是投影矩阵 $VV'$ 的第 $j$ 个对角元，满足 $0\le a_j\le1$、$\sum_j a_j=r$。因此
$$
\operatorname{tr}(F'R'RF)=\sum_{j=1}^T s_j^2 a_j\le\sum_{j=1}^r s_j^2 .
$$
上界由前 $r$ 个 $a_j=1$、其余为零达到，对应 $F$ 取 $U$ 的前 $r$ 列。给定拟合共同分量 $\hat L_{it}$，再将 $Y_{it}-\hat L_{it}$ 对 $X_{it}$ 回归，正规方程为 $(\sum X_{it}X_{it}')\hat\beta=\sum X_{it}(Y_{it}-\hat L_{it})$。这两组条件分别说明低秩分量与回归斜率如何由另一组参数的当前值更新。

[Bai（2009）](https://doi.org/10.3982/ECTA6135)

## 7.4 同方差、异方差稳健与聚类稳健方差
**从系数误差到方差矩阵。** 设 $Y=Wb+u$，$W$ 满列秩，$\Omega=\operatorname{Var}(u\mid W)$。把结构方程代入最小二乘解：
$$
\hat b=(W'W)^{-1}W'(Wb+u)=b+(W'W)^{-1}W'u .
$$
记 $A=(W'W)^{-1}W'$。给定 $W$ 时 $A$ 为常数，故 $\operatorname{Var}(Au\mid W)=E[A(u-E[u\mid W])(u-E[u\mid W])'A'\mid W]=A\Omega A'$。用 $A'=W(W'W)^{-1}$ 代入。若 $\Omega=\sigma^2I$，中间矩阵为 $W'\Omega W=\sigma^2 W'W$，与一侧的逆矩阵抵消，得 $\sigma^2(W'W)^{-1}$。

**异方差稳健方差。** 令 $w_i'$ 为第 $i$ 行，$\omega_{ij}$ 为 $\Omega$ 的元素。按矩阵乘法逐项展开 $W'\Omega W=\sum_i\sum_jw_i w_j'\omega_{ij}$。若不同观测条件不相关，非对角项消失，只余 $\sum_i w_iw_i'\sigma_i^2$。以 $\hat u_i^2$ 代替各未知条件方差得到 HC0 中间项，外侧仍各乘一次 $(W'W)^{-1}$。在同方差条件下，$\hat u=M_W u$，$M_W=I-W(W'W)^{-1}W'$，所以 $\operatorname{Var}(\hat u\mid W)=\sigma^2M_W^2=\sigma^2M_W$。第 $i$ 项为 $\sigma^2(1-h_{ii})$，其中 $h_{ii}$ 是投影矩阵的对角元。这给出 HC2 除以 $1-h_{ii}$ 的杠杆修正动机；HC3 使用 $\hat u_i^2/(1-h_{ii})^2$。

**聚类。** 若簇间独立，则 $\omega_{ij}=0$ 只对来自不同簇的观测成立。保留同簇全部项，以残差乘积估计后得
$$
\begin{aligned}
\hat B&=\sum_g\sum_{i,j\in g}w_iw_j'\hat u_i\hat u_j\\
&=\sum_g\left(\sum_{i\in g}w_i\hat u_i\right)
          \left(\sum_{j\in g}w_j\hat u_j\right)'=\sum_g s_gs_g',\\
\hat V_{\mathrm{CR0}}&=(W'W)^{-1}\hat B(W'W)^{-1}.
\end{aligned}
$$
对角元开方为各系数的聚类标准误。两维聚类中，若相关性可来自同一地区或同一年，观测对属于二者并集。由 $\mathbf1_{A\cup B}=\mathbf1_A+\mathbf1_B-\mathbf1_{A\cap B}$，对每对得分乘积应用同一恒等式，得到 $\hat V_A+\hat V_B-\hat V_{A\cap B}$（未加有限样本修正时）。

**固定效应与聚类的不同运算。** 某城市 $T_c$ 期扰动向量为 $v_c$。去均值等于左乘 $M_c=I-\mathbf1\mathbf1'/T_c$。其协方差变为 $M_c\Omega_c M_c$。即使原协方差为 $\sigma^2I$，该式也等于 $\sigma^2M_c$，非对角元为 $-\sigma^2/T_c$。这项变换自身不决定推断层级。若政策由城市统一分配，同城增加人数不增加城市层面的分配次数；簇的抽样和处理分配机制决定哪些随机性进入研究目标的方差。

[Abadie 等（2023）](https://doi.org/10.1093/qje/qjac038)

## 7.5 州级安慰剂法律与 DID 检验的拒绝率
Bertrand、Duflo 与 Mullainathan（2004）向美国州级工资数据反复指定安慰剂法律。由于虚拟法律由研究者生成，真实效应为零。若共有 $B$ 次分配，第 $b$ 次报告 $p_b$，名义水平 $\alpha$ 下的经验拒绝率为 $B^{-1}\sum_{b=1}^B\mathbf1\{p_b\le\alpha\}$。借助同一数据和同一分配规则，可以直接比较不同方差估计的实际拒绝率。

**序列相关的来源。** 设州内误差服从平稳过程 $u_{st}=\rho u_{s,t-1}+\eta_{st}$，$|\rho|<1$，创新独立于过去、跨州独立、均值零且方差为 $\sigma_\eta^2$。记平稳方差为 $v$，用创新与过去独立，得 $v=\rho^2 v+\sigma_\eta^2$，故 $v=\sigma_\eta^2/(1-\rho^2)$。对 $t>r$，两侧与 $u_{sr}$ 取协方差，创新项为零，故 $\operatorname{Cov}(u_{st},u_{sr})=\rho\operatorname{Cov}(u_{s,t-1},u_{sr})$。迭代 $t-r$ 次得 $\rho^{t-r}v$；利用协方差对称性，任意时期配对的表达式为 $\rho^{|t-r|}v$。

**序列相关如何进入政策系数。** 记 $\tilde D$ 为双向残差化的政策变量，$Q=\sum_{st}\tilde D_{st}^2$。设政策路径与误差独立。在零效应下，FWL 给出 $\hat\tau=Q^{-1}\sum_{st}\tilde D_{st}u_{st}$。条件于政策路径，$Q$ 固定。将求和的方差展开、删去跨州项、代入上段协方差：
$$
\begin{aligned}
\operatorname{Var}(\hat\tau\mid D)
&=Q^{-2}\sum_s\sum_{t,r}\tilde D_{st}\tilde D_{sr}\rho^{|t-r|}v\\
&=Q^{-2}v\left(\sum_{st}\tilde D_{st}^2+
2\sum_s\sum_{r<t}\tilde D_{st}\tilde D_{sr}\rho^{t-r}\right).
\end{aligned}
$$
把州—年当作独立观测只保留第一项。当第二项之和为正时，计算的标准误过小，零效应下也可能过度拒绝。按州聚类时，先把各期得分相加再平方，因此同州跨期乘积进入方差。其渐近依据是独立簇数增长，且单一州的得分不支配总方差。若使用原文的安慰剂分配机制，还可以固定工资路径、重新指定实施州和时期，按第 3.4 节的随机化检验公式求尾概率。

[原文](https://doi.org/10.1162/003355304772839588)
