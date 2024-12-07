#set page("a4", columns: 3, flipped: true, margin: (x: 3%, y: 4%), header: [
  _MH3510 Cheat Sheet_
  #h(1fr)
  _Pu Fanyi_
])

#set text(font: ("Roboto", "Noto Sans CJK SC"), lang: "zh", region: "cn", size: 9pt)

#set heading(numbering: "1.")

#let SSE = math.class(
  "normal",
  math.op("SSE")
)

#let SSR = math.class(
  "normal",
  math.op("SSR")
)

#let syy = math.class(
  "normal",
  $S_(y y)$
)

#let SST = math.class(
  "normal",
  math.op("SST")
)

#let msreg = math.class(
  "normal",
  $op("MS")_op("Reg")$
)

#let df = math.class(
  "normal",
  math.op("df")
)

#let SSEXT = math.class(
  "normal",
  math.op("SSEXT")
)

#let sumall = math.class(
  "normal",
  $y_(dot dot)$
)

#let cf = math.class(
  "normal",
  math.op("c.f.")
)

#let MST = math.class(
  "normal",
  math.op("MST")
)

#let MSE = math.class(
  "normal",
  math.op("MSE")
)

#let ss = math.class(
  "normal",
  math.op("SS")
)

= Regression

== Simple Linear Regression

$cases(S_(x x) = sum_(i=1)^n (x_i - overline(x))^2,
S_(x y) = sum_(i=1)^n (x_i-overline(x))(y_i - overline(y)),
S_(y y) = sum_(i=1)^n (y_i - overline(y))^2) =>hat(beta)_0 = overline(y) - hat(beta)_1overline(x), hat(beta)_1 = S_(x y) / S_(x x)$

$(hat(beta)_1-beta_1)/(sqrt(sigma^2\/S_(x x))) tilde cal(N)(0, 1), (hat(beta)_0-beta_0)/sqrt((sigma^2 sum_(i=1)^n x_i^2)\/ (n S_(x x))) tilde cal(N)(0, 1)$

$(hat(beta)_1-beta_1)/(sqrt(s^2\/S_(x x))) tilde t_(n-2), (hat(beta)_0-beta_0)/sqrt((s^2 sum_(i=1)^n x_i^2)\/ (n S_(x x))) tilde t_(n-2)$

$(hat(y)_0 - bb(E)[y_0]) / (s sqrt(1/n + (x_0 - overline(x))^2 / (S_(x x)))) tilde t_(n-2), (y_text("new") - hat(y)_text("new"))/(s sqrt(1 + 1/n + (x_0 - overline(x))^2 / (S_(x x)))) tilde t_(n-2)$

$cal(R)^2 = SSR / syy = S^2_(x y) / (S_(x x) syy) = r^2_(x y)$

$cal(R)^2$ of the observed variation in the response variable can be explained
by the regression model.

$bb(E)[msreg] = sigma^2 + beta_1 S_(x x)$

$hat(beta)_0 = l^top y, hat(beta)_1 = k^top y, k_i = (x_i - overline(x)) / S_(x x), l_i = 1/n - k_i overline(x)$


// $\sim\mathcal{N}(0, 1), \frac{\hat{\beta}_0-\beta_0}{\sqrt{\left(\sigma^2\sum_{i=1}^{n}x_i^2\right)/\left(nS_{xx}\right)}}\sim\mathcal{N}(0, 1)$

// $\frac{\hat{\beta}_1-\beta_1}{\sqrt{s^2/S_{xx}}}\sim t_{n-2}, \frac{\hat{\beta}_0-\beta_0}{\sqrt{\left(s^2\sum_{i=1}^{n}x_i^2\right)/\left(nS_{xx}\right)}}\sim t_{n-2}$



== Multiple Linear Regression

=== Useful Formula

$bb(E)[y y^top]=bb(E)[y] bb(E)[y]^top + Sigma_y$

$bb(E)[y^top A y] = bb(E)[y]^top A bb(E)[y] + tr(A Sigma_y)$

$op("Cov") (A x + b) = A op("Cov") (x) A^top$

$x bot y arrow.r.double Sigma_(x+y) = Sigma_x + Sigma_y$

Outliers: $|e^*_i| = abs(e_i / sqrt(MSE))>3$

QQ-Plot: 左边：往下尾巴大；右边：往上尾巴大

=== Normal Distributions

$1/((2 pi)^n abs(Sigma)^(1/2)) exp (-1/2 (y - mu)^top Sigma^(-1) (y-mu) )$

$a^top y tilde cal(N)(a^top mu, a^top Sigma a)$

=== Definition

$y = X beta + epsilon, space epsilon tilde cal(N)(0, sigma^2 I)$

$y, beta, epsilon in bb(R)^(n times 1), X in bb(R)^(n times (p + 1))$

$hat(beta) = (X^top X)^(-1) X^top y, e = y - X hat(beta)$

=== Hat Matrix

$H = X (X^top X)^(-1) X^top in bb(R)^(n times n)$

$hat(y) = H y, e = (I - H) y$

$H$ 和 $I-H$ 是投影矩阵，$H y bot (I-H) y$：

$cases(
  H X = X,
  H^2 = H,
  (I-H)^2=I-H
)$

$H$ and $I-H$ are orthogonal matrices, so

$cases(
  H^top = H,
  (I-H)^top = I-H
)$

=== Sum of Squares

$A <- 1/n 1 1^top$ （全 $1/n$ 矩阵，$bb(R)^(n times n)$）

$syy = sum_(i=1)^n (y_i - overline(y))^2 = y^top (I-A) y$

$SSE = sum_(i=1)^n (y_i - hat(y)_i)^2 = e^top e = y^top (I - H) y$

$SSR = sum_(i=1)^n (hat(y)_i - overline(y))^2 = y^top (H-A) y$

$SSR tilde sigma^2 chi^2_p, SSE tilde sigma^2 chi^2_(n-p-1)$

=== ANOVA Table

#figure(
  table(
    align: horizon,
    inset: 5pt,
    columns: (auto, auto, auto, auto, auto),
    table.header([Source], [DF], [SS], [MS], [$cal(F)$]),
    [Regression], [$p$], [$SSR$], [$msreg$], [$F$],
    [Residual], [$n-p-1$], [$SSE$], [$s^2$], [],
    [Total], [$n-1$], [$syy$], [], []
  )
)

$bb(E)[s^2]=sigma^2, F = msreg / s^2 tilde cal(F)(p, n - p - 1)$

$cal(H)_0: beta_1 = dots.c = beta_p = 0 ==> p = bb(P)(cal(F)_(p, n-p-1)>F)$

=== Cochran's Theorem

$y tilde cal(N)_(n) (0, sigma^2 I)$, $A_i$ 满足：

1. $sum_(i=1)^m A_i = I$
2. $A_i^2 = A_i$（投影矩阵）
3. $A_i^top = A$

我们有：

1. $y A_i y$ 是互相 independent 的
2. $y A_i y tilde sigma^2 chi^2_(r_i)$，其中 $r_i = tr(A_i)$

=== $cal(R)^2$ Statistic

$cal(R)^2 = SSR / syy$

$p arrow.t$ 会导致 $cal(R)$ 接近 $1$，因此我们引入 adjusted $cal(R)^2$：

$cal(R)^2_a = 1 - (SSE \/ (n-p-1)) / (syy \/ (n - 1))$

$cal(R)^2$ too small $->$ cannot predict precisely

$cal(R)^2$ too large $->$ overfit

=== $t$ Statistic

$hat(beta) tilde cal(N)_(p+1) (beta, sigma^2 (X^top X)^(-1))$

$hat(beta_i) tilde cal(N)(beta_i, sigma^2(X^top X)^(-1)_(i+1, i+1))$

$cal(H)_0: beta_i=0 => t^* = hat(beta_i) / (op("SE")(hat(beta_i))) = hat(beta_i) / sqrt(s^2 (X^top X)^(-1)_(i + 1, i + 1))$

=== Extra Sum of Squares

$SSR(X_2 | X_1) &= SSR(X_1, X_2) - SSR(X_1)\ &= SSE(X_1) - SSE(X_1, X_2)$

$op("SSEXT") = SSE_omega - SSE_Omega = SSR_Omega - SSR_omega$

#figure(
  table(
    align: horizon,
    inset: 5pt,
    columns: (auto, auto, auto, auto, auto),
    table.header([Source], [DF], [SS], [MS], [$cal(F)$]),
    [$omega$], [$q$], [$SSR_omega$], [$SSR_omega / q$], [],
    [Extra], [$p-q$], [$SSEXT$], [$SSEXT / (p-q)$], [$F$],
    [Residual], [$n-p-1$], [$SSE$], [$SSE_Omega / (n-p-1)$], [],
    [Total], [$n-1$], [$syy$], [], []
  ),
)

$F = (SSEXT \/ (p-q))/(SSE_Omega \/ (n-p-1)) tilde cal(F)(p-q, n-p-1)$

$cal(H)_0$: We should include the extra variables

$p = bb(P)(cal(F)_(p-q, n-p-1) > F)$.

= Classification

== One-Way Classification Model

=== Definition

- *Factor*: Gender | Race
- *Levels*: M, F | CN, JA, RU
- *No. of levels*: $2$ | $3$
- *Treatments*: (M, CN), (F, CN), ...
- *No. of treatments*: $2 times 3 = 6$

$y_(i j) = theta_i + epsilon_(i j), space i in [1, k], j in [1, n_i], space epsilon stretch(tilde)^"iid" cal(N)(0, sigma^2)$

*Re-parameterization*: $theta_i = mu + tau_i$: $y_(i j) = mu + tau_i + epsilon_(i j)$

- Overall sample sum: $sumall = sum_(i=1)^k sum_(j=1)^(n_i) y_(i j)$
- $i$-th sample sum: $y_(i dot) = sum_(j=1)^(n_i) y_(i j)$
- Overall sample mean: $overline(y)_(dot dot) = y_(dot dot) / n$
- $i$-th sample mean: $overline(y)_(i dot) = y_(i dot) / n_i$

=== Sum of Squares

Correction factor: $cf = n (overline(y)_(dot dot))^2 = (y_(dot dot)^2) / n $

Total corrected SS: $syy = sum_(i=1)^k sum_(j=1)^(n_i) (y_(i j) - overline(y)_(dot dot))^2 = sum_(i=1)^k sum_(j=1)^(n_i) y_(i j)^2 - cf $

Between-group SS: $SST = sum_(i=1)^k n_i (overline(y)_(i dot) - overline(y)_(dot dot))^2 = sum_(i=1)^k y_i / n - cf$

Within-group SS: $SSE = sum_(i=1)^k sum_(j=1)^(n_i) (y_(i j) - overline(y)_(i dot))^2 = syy - SST$

=== One-way ANOVA Table

#figure(
  table(
    align: horizon,
    inset: 5pt,
    columns: (auto, auto, auto, auto, auto),
    table.header([Source], [DF], [SS], [MS], [$cal(F)$]),
    [Between groups], [$k-1$], [$SST$], [$MST$], [$F$],
    [Within groups], [$n-k$], [$SSE$], [$MSE$], [],
    [Total], [$n-1$], [$syy$], [], []
  ),
)

$MST = SST / k, MSE = SSE / (n-k), bb(E)[MSE] = sigma^2$

$F = MST / MSE tilde cal(F)_(k-1, n-k)$

=== Analysis of Parameters

我们很多时候喜欢研究关于 $theta$ 的一些线性式子，比如比较两个类别（$theta_i-theta_j$）等等，一般地：

$L=alpha^top theta, hat(L)=alpha^top overline(y)$

$hat(L) tilde cal(N)(L, sigma^2 sum_(i=1)^k alpha_i^2 / n_i)$

$hat(L) tilde t_(n-k)(L, s^2 sum_(i=1)^k alpha_i^2 / n_i)$

== Two-Way ANOVA

=== Definition

$y_(i j k) = theta_(i j) + epsilon_(i j k), space epsilon_(i j k) tilde cal(N)(0, sigma^2)$

$i = 1, dots, I, space j=1, dots, J, space k = 1,dots, r$

=== Re-parameterization

$theta_(i j) = mu + alpha_i + beta_j + (alpha beta)_(i j)$

$sum alpha_i = sum beta_j = sum_(i=1)^I (alpha beta)_(i j) = sum_(j=1)^J (alpha beta)_(i j) = 0$

=== Model Fitting

$cases(
  hat(mu) = overline(y)_(dot dot dot)\
  hat(alpha)_i = overline(y)_(i dot dot) - hat(mu)\
  hat(beta)_i = overline(y)_(dot j dot) - hat(mu)\
  hat(alpha beta)_(i j) = overline(y)_(i j) - hat(mu) - hat(alpha)_i - hat(beta)_j
)$

=== Sum of Squares

Correction factor: $cf = y^2_(dot dot dot) / (I J r)$

$syy = sum_i sum_j sum_k (y_(i j k) - overline(y)_(dot dot dot))^2$

$SST = r sum_i sum_j (overline(y)_(i j dot) - overline(y)_(dot dot dot))^2$

$SSE = sum_i sum_j sum_k (y_(i j k) - overline(y)_(i j dot))^2$

$ss_A = r J sum_i (overline(y)_(i dot dot) - overline(y)_(dot dot dot))^2$

$ss_B = r I sum_j (overline(y)_(dot j dot) - overline(y)_(dot dot dot))^2$

$ss_(A B) = r sum_i sum_j (overline(y)_(i j dot) - overline(y)_(i dot dot) - overline(y)_(dot j dot) + overline(y)_(dot dot dot))^2$

$syy = SST + SSE$

=== Two-way ANOVA Table

#figure(
  table(
    align: horizon,
    inset: 5pt,
    columns: (auto, auto, auto, auto, auto),
    table.header([Source], [DF], [SS], [MS], [$cal(F)$]),
    [Main eff. $A$], [$I - 1$], [$ss_A$], [$op("MS")_A$], [$F_A$],
    [Main eff. $B$], [$J - 1$], [$ss_B$], [$op("MS")_B$], [$F_B$],
    [Interaction], [$(I - 1)(J - 1)$], [$ss_(A B)$], [$op("MS")_(A B)$], [$F_(A B)$],
    [Error], [$I J (r - 1)$], [$SSE$], [$MSE$], [],
    [Total], [$I J r - 1$], [$syy$], [], [],
  ),
)

$op("MS")_A = ss_A / (I - 1), op("MS")_B = ss_B / (J - 1), op("MS")_(A B) = ss_(A B) / ((I - 1)(J - 1))$

$bb(E)[op("MS")_A] = sigma^2 + r J sum_i (overline(theta)_(i dot) - mu)^2 / (I - 1)$

$bb(E)[op("MS")_B] = sigma^2 + r I sum_j (overline(theta)_(dot j) - mu)^2 / (J - 1)$

$bb(E)[op("MS")_(A B)] = sigma^2 + r sum_i sum_j (theta_(i j) - overline(theta)_(i dot) - overline(theta)_(dot j) + mu)^2 / ((I - 1)(J - 1))$

$MSE = SSE / (I J (r - 1)), bb(E)[MSE] = sigma^2$

=== Analysis of Parameters

Parameters of interest: $L = sum_(i=1)^I sum_(j = 1)^J alpha_(i j) theta_(i j)$

$hat(L) = sum_(i=1)^I sum_(j = 1)^J alpha_(i j) overline(y)_(i j)$

$hat(L) tilde cal(N)(L, sigma^2 sum_(i=1)^I sum_(j=1)^J alpha_(i j)^2/r)$

$hat(L) tilde t_(I J (r - 1))(L, s^2 sum_(i=1)^I sum_(j=1)^J alpha_(i j)^2/r)$

= Universal Linear Modelling

== Dummy Variable 

Similar to one-hot: $k$ levels $=>$ $I_1, I_2, dots, I_(k-1)$

== ANCOVA

=== One Factor

$Y = mu + tau_1 T_1 + dots.c + tau_(k-1) T_(k-1) + gamma_1 X_1 + dots.c + gamma_p X_p + epsilon$

#figure(
  table(
    align: horizon,
    inset: 5pt,
    columns: (auto, auto, auto, auto, auto),
    table.header([Source], [DF], [SS], [MS], [$cal(F)$]),
    [Reg. fitting], [$p$], [$SSR_omega$], [$op("MSR")_omega$], [],
    [Extra], [$k - 1$], [$SSEXT$], [$op("MSEXT")$], [$F$],
    [Residual], [$n - p - k$], [$SSE_Omega$], [$MSE_Omega$], [],
    [Total], [$n - 1$], [$syy$], [], [],
  ),
)

$op("MSR")_omega = SSR_omega / p, op("MSEXT") = SSEXT / (k - 1), MSE_Omega = SSE_Omega / (n - p - k)$

$F = op("MSEXT") / MSE_Omega tilde cal(F)(k - 1, n - p - k)$

=== Two Factors

$Y = &mu + sum_(i=1)^(I - 1) alpha_i A_i
 + sum_(j=1)^(J - 1) beta_j B_j\
& + sum_(i=1)^(I - 1) sum_(j = 1)^(J - 1) (alpha beta)_(i j) A_i B_j + sum_(k=1)^p gamma_p X_p + epsilon$

Reduce both:

#figure(
  table(
    align: horizon,
    inset: 5pt,
    columns: (auto, auto, auto, auto, auto),
    table.header([Source], [DF], [SS], [MS], [$cal(F)$]),
    [Reg. fitting], [$p + I + J - 2$], [$SSR_omega$], [$op("MSR")_omega$], [],
    [Extra], [$(I - 1)(J - 1)$], [$SSEXT$], [$op("MSEXT")$], [$F$],
    [Residual], [$n - p - I J$], [$SSE_Omega$], [$MSE_Omega$], [],
    [Total], [$n - 1$], [$syy$], [], [],
  ),
)

$op("MSR")_omega = SSR_omega / (p + I + J - 2), op("MSEXT") = SSEXT / ((I - 1) (J - 1)), MSE_Omega = SSE_Omega / (n - p - I J)$

$F = op("MSEXT") / MSE_Omega tilde cal(F)((I - 1)(J - 1), n - p - I J)$

Reduct one:

#figure(
  table(
    align: horizon,
    inset: 5pt,
    columns: (auto, auto, auto, auto, auto),
    table.header([Source], [DF], [SS], [MS], [$cal(F)$]),
    [Concomitant], [$p$], [$SSR''$], [], [],
    [Main eff. $B$], [$J - 1$], [$ss_B$], [$op("MS")_B$], [$F_B$],
    [Main eff. $A$], [$I - 1$], [$ss_A$], [$op("MS")_A$], [$F_A$],
    [Interaction], [$(I - 1)(J - 1)$], [$SSEXT$], [$op("MSEXT")$], [$F_(A B)$],
    [Error], [$n - p - I J$], [$SSE_Omega$], [$MSE_Omega$], [],
    [Total], [$n - 1$], [$syy$], [], [],
  ),
)

$F_A = op("MS")_A / (MSE_Omega), F_B = op("MS")_B / (MSE_Omega), F_(A B) = op("MSEXT") / (MSE_Omega)$

Interaction between factors and concomitant variables:

$Y = &mu + sum_(i=1)^(I - 1) alpha_i A_i + sum_(j=1)^(J - 1) beta_j B_j + sum_(i=1)^(I - 1)sum_(j=1)^(J - 1)(alpha beta)_(i j) A_i B_j\
& + [
  sum_(i=1)^(I - 1) delta_i A_i+sum_(j=1)^(J - 1) zeta_j B_j + sum_(i=1)^(I - 1)sum_(j=1)^(J - 1)(delta zeta)_(i j) A_i B_j
]X\
&+ gamma X+ epsilon $
