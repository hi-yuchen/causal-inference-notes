# 控制变量、固定效应与标准误

协变量和固定效应共同决定处理变量的残差及其在投影系数中的权重。给定回归规格，系数的抽样方差取决于残差化处理与误差乘积的协方差结构。以下先推导点估计，再讨论方差估计。

## 7.1 遗漏变量偏误与 Frisch–Waugh–Lovell 定理

设总体机制为 $Y=\alpha+\tau D+\gamma X+u$，$E[u\mid D,X]=0$。只回归 $Y$ 对 $D$，其总体斜率为
$$
b=\frac{\operatorname{Cov}(D,Y)}{\operatorname{Var}(D)}
=\tau+\gamma\frac{\operatorname{Cov}(D,X)}{\operatorname{Var}(D)}.
$$
对二元处理，$\operatorname{Cov}(D,X)/\operatorname{Var}(D)=E[X\mid D=1]-E[X\mid D=0]$。若两组的处理前收入相差 2 个单位，且 $\gamma=3$，则 $b-\tau=6$。在该模型下，短回归与长回归的总体系数相等，当且仅当 $\gamma\operatorname{Cov}(D,X)=0$。

现在将控制变量及截距排入矩阵 $X$，假设满列秩，令 $M_X=I-X(X'X)^{-1}X'$。最小化 $\|Y-D\tau-X\gamma\|^2$，固定 $\tau$ 对 $\gamma$ 求解得 $\hat\gamma(\tau)=(X'X)^{-1}X'(Y-D\tau)$。代回目标函数得到 $\|M_X(Y-D\tau)\|^2$。再对 $\tau$ 求一阶条件，用 $M_X'=M_X=M_X^2$：
$$
\hat\tau=\frac{D'M_XY}{D'M_XD}
=\frac{\sum_i\tilde D_i\tilde Y_i}{\sum_i\tilde D_i^2}.
$$
其中 $\tilde D=M_XD$、$\tilde Y=M_XY$。该系数唯一存在的条件为 $D'M_XD>0$。若 $\operatorname{Var}(u\mid D,X)=\sigma^2 I$，则 $\operatorname{Var}(\hat\tau\mid D,X)=\sigma^2/(D'M_XD)$。给定误差方差，控制后剩余的处理变异越小，估计方差越大。

## 7.2 混杂、中介和碰撞变量

若培训前收入同时影响报名和未参训状态下的后续收入，则它构成混杂因素。若培训影响就业，就业继而影响收入，则就业构成中介。设结构方程为 $M=aD+v$、$Y=bD+cM+u$。固定外生扰动，将 $D$ 增加一单位使 $M$ 增加 $a$，使 $Y$ 增加 $b+ac$。固定 $M$ 的直接效应为 $b$；在线性模型中，$E[u\mid D,M]=0$ 是将控制中介后的回归系数解释为 $b$ 的充分条件。

设 $D$ 与能力 $U$ 独立且方差均为正，入选指标为 $S=D+U$，结果为 $Y=U$。$S$ 是路径 $D\to S\leftarrow U\to Y$ 上的碰撞变量。边际上 $\operatorname{Cov}(D,Y)=0$。给定 $S$ 有 $Y=S-D$，故 $Y$ 对 $(S,D)$ 的回归中，$D$ 的系数为 $-1$，而干预 $D$ 对 $Y$ 的效应为零。 [Cinelli、Forney 与 Pearl](https://doi.org/10.1177/00491241221099552)

对离散协变量 $X$ 引入全套分层截距。设处理为二元，$e_x=P(D=1\mid X=x)$、$p_x=P(X=x)$，$m_d(x)=E[Y\mid D=d,X=x]$。由 FWL，$\tilde D=D-e_X$。层内有 $E[\tilde D Y\mid X=x]=e_x(1-e_x)[m_1(x)-m_0(x)]$ 及 $E[\tilde D^2\mid X=x]=e_x(1-e_x)$。在条件可交换性、一致性和重叠条件下，记 $\tau_x=m_1(x)-m_0(x)$，则共同斜率为 $\sum_x p_xe_x(1-e_x)\tau_x/\sum_x p_xe_x(1-e_x)$。同一总体的 ATE 为 $\sum_x p_x\tau_x$；两者之差等于 $\operatorname{Cov}(e_X(1-e_X),\tau_X)/E[e_X(1-e_X)]$。

## 7.3 加性固定效应、类别交互与交互式因子模型

设 $Y_{it}=\alpha_i+\tau D_{it}+u_{it}$。对任意变量 $V$，其在个体虚拟变量空间上的样本投影为个体均值 $\bar V_i$。将投影从原变量中扣除，得到
$$Y_{it}-\bar Y_i=\tau(D_{it}-\bar D_i)+(u_{it}-\bar u_i).$$
组内估计量的分母为 $\sum_{i,t}(D_{it}-\bar D_i)^2$。固定效应吸收时间不变的加性项 $\alpha_i$。若时变健康冲击同时引起参保状态和结果变化，残差化后仍可有 $E[(D_{it}-\bar D_i)(u_{it}-\bar u_i)]\ne0$。平衡等权面板中，单位和时期双向固定效应对应的处理残差为 $D_{it}-\bar D_i-\bar D_t+\bar D$。不平衡面板中的精确残差由对两个虚拟变量空间之和的联合投影确定。

设 $A_{ct}$ 为城市 $c$ 在年份 $t$ 的试点状态。城市×年份固定效应吸收 $A_{ct}$。若处理定义为 $D_{ict}=A_{ct}E_{ict}$，其中 $E_{ict}$ 为个人资格，则仅控制城市×年份固定效应后有 $\tilde D_{ict}=A_{ct}(E_{ict}-\bar E_{ct})$。此时处理变异来自试点城市当年的资格差异。若另加资格主效应，则系数依赖试点和非试点单元之间资格差异的比较。

Bai 的 interactive fixed effects 写作
$$Y_{it}=X_{it}'\beta+\lambda_i'f_t+u_{it}.$$
设 $\Lambda$ 与 $F$ 的各行为 $\lambda_i'$ 与 $f_t'$，共同分量为 $\Lambda F'$，其秩不超过 $r$。若 $f_t$ 表示全国性需求冲击，地区 $A,B$ 的未处理趋势差含有 $(\lambda_A-\lambda_B)'(f_1-f_0)$；共同年份截距仅吸收载荷相同的部分。给定 $\beta$，最小二乘因子拟合取残差矩阵的秩 $r$ 截断奇异值分解。对可逆矩阵 $H$，载荷和因子的变换 $(\Lambda H,FH^{-T})$ 保持拟合值不变。因子归一化固定这一不确定性。斜率识别还要求解释变量含有可与共同因子分离的变异。 [Bai（2009）](https://doi.org/10.3982/ECTA6135)

## 7.4 同方差、异方差稳健与聚类稳健方差

这里以 $W$ 表示含全部回归变量的满秩矩阵。由 $\hat b-b=(W'W)^{-1}W'u$，设给定 $W$ 的误差协方差矩阵为 $\Omega$，可得
$$V(\hat b\mid W)=(W'W)^{-1} W'\Omega W (W'W)^{-1}.$$
若 $\Omega=\sigma^2I$，则上式化为 $\sigma^2(W'W)^{-1}$。若误差条件不相关但方差可异，HC0 以 $\sum_i w_iw_i'\hat u_i^2$ 估计中间项，其中 $w_i'$ 为 $W$ 的第 $i$ 行。令 $h_{ii}$ 为投影矩阵的第 $i$ 个对角元。HC2 和 HC3 分别以 $\hat u_i^2/(1-h_{ii})$ 和 $\hat u_i^2/(1-h_{ii})^2$ 调整残差平方。

设不同簇相互独立，簇内协方差任意。令 $s_g=\sum_{i\in g}w_i\hat u_i$，聚类稳健中间项为 $\sum_g s_gs_g'$。每个外积保留同簇全部 $i,j$ 的交叉项。若相关性来自两个非嵌套维度 $A,B$，两维聚类估计为 $\hat V_A+\hat V_B-\hat V_{A\cap B}$，交集项扣除重复计入的观测对。该近似依赖各维度的有效独立簇数。

城市固定效应消除加性的城市常量。设去除城市常量前的时变扰动为 $v_{ct}$，其协方差矩阵为 $\Omega_c$。令 $M_c=I-\mathbf 1\mathbf 1'/T_c$，去均值后的协方差矩阵为 $M_c\Omega_c M_c$，通常含有非零的非对角元。设计型推断还取决于城市的抽取方式和政策分配层级。同城增加个体观测时，城市层面的独立分配次数保持不变。 [Abadie 等（2023）](https://doi.org/10.1093/qje/qjac038)

## 7.5 州级安慰剂法律与 DID 检验的拒绝率

Bertrand、Duflo 与 Mullainathan（2004，QJE）对美国女性工资数据反复指定虚拟法律。每次指定给出一组实施州和实施年份；由于法律由研究者生成，其真实效应为零。模拟拒绝率由拒绝零假设的次数除以重复次数得到，因而可以与名义检验水平直接比较。 [原文](https://doi.org/10.1162/003355304772839588)

设州内误差满足平稳过程 $u_{st}=\rho u_{s,t-1}+\eta_{st}$，$0<\rho<1$，创新跨州、跨期独立，方差为 $\sigma_\eta^2$，且与政策分配独立。记 $Q=\sum_{s,t}\tilde D_{st}^2$，其中 $\tilde D$ 为剔除州和年份固定效应后的处理变量。条件于处理矩阵，政策系数的方差为 $Q^{-2}\sum_s\sum_{t,r}\tilde D_{st}\tilde D_{sr}\sigma_\eta^2\rho^{|t-r|}/(1-\rho^2)$。将州—年视为独立仅保留 $t=r$ 项。当遗漏交叉项之和为正时，标准误被低估，点估计仍可保持无偏。

按州聚类将每州各期的得分先求和再取外积，保留上式的州内跨期项。其大样本近似以独立州的数量增长为基础，且要求单个州的得分不支配总体方差。在原论文的安慰剂设计中，还可依据虚拟法律的已知分配规则重新指定实施州和时期，得到零效应下检验统计量的随机化分布。
