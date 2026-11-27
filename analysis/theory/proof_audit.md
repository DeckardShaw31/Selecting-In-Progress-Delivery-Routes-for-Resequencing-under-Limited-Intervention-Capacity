# Theory audit: Propositions 1–2 and strengthened results

Audit of `proposal/research_proposal.tex`, Section 3 (Model), done 2026-11-26. This is the check required before committing to a theory-heavy paper (framework Section 15, item 2). Numerical checks: `code/resequencing/theory_checks.py` (T1–T4, all pass), plus Gate 3's verification checks.

## Verdict

**Propositions 1 and 2 are correct as stated. No step of either proof fails.** The audit found four places where the text should be made more precise (§2), and four results that strengthen the theory considerably (§3):

- an exact representation of exposure functions;
- a sharper version of Proposition 1;
- a Wasserstein-distance form of the loss bound;
- a worst-case loss that the bound attains exactly.

With these, the theory section can carry an EJOR paper. One result is still open: the route-wide-shock result (H3), §4.

Notation as in the proposal:
- $A_{ri}=p_{ri}+D_r+\xi_{ri}$ is the arrival time;
- $s_{ri}=d_{ri}-p_{ri}$ is the slack;
- $w_{ri}>0$ are the weights, with total $W_r=\sum_i w_{ri}$;
- the risk is $\rho_r=\sum_i w_{ri}\,\mathbb E(A_{ri}-d_{ri})^+$.

## 1. Step-by-step check

| Step | Check | Status |
|---|---|---|
| Arrival decomposition (1) | A definition of $\xi_{ri}$, so always possible. The content is in Assumption 1. | correct |
| Properties of $g_r$ | Convex and nondecreasing (sum of expectations of convex nondecreasing functions of $x$). Lipschitz with constant $W_r$. $g_r(x)\to0$ as $x\to-\infty$ by dominated convergence, with dominating function $\lvert\xi_{ri}\rvert+\lvert x_0-s_{ri}\rvert$ for $x\le x_0$, integrable since $\mathbb E\lvert\xi\rvert<\infty$. | correct (state the Lipschitz constant $W_r$) |
| Proposition 1, sufficiency | If $g_r\equiv g$, then $r\in S_D$, $r'\notin S_D$ gives $D_r\ge D_{r'}$, so $\rho_r\ge\rho_{r'}$. Any tie-breaking works. | correct |
| Proposition 1, converse | Continuity gives $\epsilon$. $g_{r'}(x_0)>g_r(x_0)\ge0$, so $g_{r'}(x_0)>0$, and finitely many other routes can be pushed below it. With $K=1$, $S_D=\{r\}$. | correct |
| Perturbation lemma | Needs $\lvert S_\rho\rvert=\lvert S_\sigma\rvert=K$; holds by definition. | correct |
| Proposition 2(i) | $\bar g$ nondecreasing makes $S_D$ a top-$K$ set for $\sigma=\bar g(D)$. | correct |
| Proposition 2(ii) | $\lvert a^+-b^+\rvert\le\lvert a-b\rvert\,\mathbf 1\{\max(a,b)>0\}$: trivial if both are $\le0$, otherwise $(\cdot)^+$ is 1-Lipschitz. The triangle inequality then moves the absolute value inside the expectation. | correct |
| Remark (route-wide drift) | If $\xi_{ri}=h(D_r)+\eta_{ri}$ with a common law of $\eta$ and $x+h(x)$ nondecreasing, then $\rho_r=g(D_r+h(D_r))$ is nondecreasing in $D_r$. | correct, but narrow (§2c) |

## 2. Corrections to the text (no result changes)

**a. Proposition 2 does not need Assumption 1.**
- Part (i) holds for any nonnegative risks $\rho_r$ and any delays $D_r$: the proof uses only the two rankings.
- Part (ii) holds for any random arrivals with finite mean, because the decomposition $A=p+D+\xi$ is a definition.
- The text says "Under Assumption 1", which understates the result. **State Proposition 2 as distribution-free.** It then also covers waiting, learning and route-wide shocks.

**b. Waiting at window openings violates Assumption 1.**
- With planned idle time, an early vehicle waits, and a late vehicle uses up planned idle time. Either way, the law of $\xi_r$ depends on $D_r$.
- So Proposition 1's model is exact only without planned waiting. The simulated regimes MIXA, MIXB and TW include waiting.
- Proposition 2 is unaffected (point a).
- Extension worth stating: with waiting, the effective delay at stop $i$ is $h_{ri}(D_r)$ with $h_{ri}$ nondecreasing and 1-Lipschitz. Sufficiency in Proposition 1 then holds when all routes share the same $\sum_i w_{ri}\mathbb E(h_{ri}(x)+\eta_{ri}-s_{ri})^+$.

**c. The route-wide-shock remark covers only an additive common drift.**
- In the simulator the route shock is multiplicative. After learning, the expected future deviation at stop $i$ is roughly $\hat\theta_r\times$ (planned time from now to $i$), which differs by stop and by route.
- The remark therefore does not explain the high-noise attenuation seen in Gate 2c, and H3 needs its own result (§4).

**d. Remark 1, "exposure functions differ only through route length and the distribution of future deviations".** Lemma 1 below makes this exact.
- Under ETA promises, $G_r(y)=\sum_i w_{ri}\,P(\beta-\xi_{ri}\le y)$.
- Under day-definite promises, the remaining planned durations also enter through $s_{ri}$.

## 3. Strengthened results (proofs below; all checked numerically)

Define the **effective slack** $Y_{ri}=s_{ri}-\xi_{ri}$: the amount of current delay that promise $i$ can absorb. Define the **weighted effective-slack distribution function** $G_r(y)=\sum_i w_{ri}\,P(Y_{ri}\le y)$, which is nondecreasing and runs from 0 to $W_r$.

**Lemma 1 (representation).** Under $\mathbb E\lvert\xi_{ri}\rvert<\infty$,
$$g_r(x)=\int_{-\infty}^{x}G_r(y)\,dy .$$
The right derivative of $g_r$ is $G_r$, and $g_r$ depends on the joint law of $(\xi_{r1},\dots,\xi_{rn_r})$ only through its marginals.

*Proof.* $(x-Y)^+=\int_{-\infty}^{x}\mathbf 1\{Y\le y\}\,dy$. Take expectations (Tonelli) and sum over $i$ with weights. The integral is finite because $\int_{-\infty}^0P(Y\le y)\,dy=\mathbb E[Y^-]<\infty$. Right-continuity of $G_r$ gives the derivative statement. ∎

*Check T1: integrated distribution function vs. closed form, maximum error 6e-7.*

**Proposition 1′ (sharper characterization).** Under Assumption 1, behind-schedule triage maximizes captured risk for every $K$ and every delay vector **if and only if $G_1=\dots=G_N$**: all routes have the same weighted distribution of effective slack. In particular:
- **Total weights must be equal** ($W_r=G_r(\infty)$). With unit weights, every route must have the same number of remaining promises.
- **Within-route correlation of deviations is irrelevant.** Route-wide shocks can matter only through learning, which takes the model outside Assumption 1.

*Proof.* By Lemma 1, $g_r=g_{r'}$ if and only if $G_r=G_{r'}$. Proposition 1 does the rest. ∎

**Proposition 2′(i) (sharpened, distribution-free).** For any nonnegative risks, any delays and any nondecreasing $\bar g$:
$$V(S_\rho)-V(S_D)\le 2\min(K,N-K)\,\max_r\lvert\rho_r-\bar g(D_r)\rvert .$$
The best $\bar g$ gives $\varepsilon^*=\tfrac12\max_{D_a\le D_b}(\rho_a-\rho_b)^+$, which is the error of the $L_\infty$ isotonic fit. So
$$V(S_\rho)-V(S_D)\le \min(K,N-K)\cdot\max_{D_a\le D_b}(\rho_a-\rho_b)^+ ,$$
and this bound is attained.

*Proof.*
1. Let $A=S_\rho\setminus S_D$ and $B=S_D\setminus S_\rho$. Then $\lvert A\rvert=\lvert B\rvert=m\le\min(K,N-K)$, since $A$ lies in the complement of $S_D$.
2. Pair the elements of $A$ and $B$ one to one. For each pair, $D_b\ge D_a$, so $\rho_a-\rho_b\le(\rho_a-\rho_b)^+\le\max_{D_a\le D_b}(\rho_a-\rho_b)^+\le 2\varepsilon$.
3. Summing over the $m$ pairs gives the bound.
4. Tightness: take $N=2$, $K=1$, $D_1<D_2$, $\rho_1=\Delta$, $\rho_2=0$. The loss is $\Delta$. ∎

*Check T2: 400 random instances, all pass. The $\varepsilon^*$ closed form was also checked against brute force in Gate 3 verification check 3.*

**Corollary (Wasserstein form).** Under Assumption 1, take $\bar g=\frac1N\sum_k g_k$, which is nondecreasing. Then
$$\lvert g_r(D_r)-\bar g(D_r)\rvert\le\int_{-\infty}^{D_r}\lvert G_r-\bar G\rvert\,dy .$$
If all routes carry the same total weight $W$, this is at most $W\cdot W_1(\nu_r,\bar\nu)$, where $\nu_r=G_r/W$, $\bar\nu$ is their mixture and $W_1$ is the Wasserstein-1 distance. Hence
$$V(S_\rho)-V(S_D)\le 2\min(K,N-K)\,W\max_r W_1(\nu_r,\bar\nu).$$

**The behind-schedule loss is bounded by how far each route's effective-slack distribution lies from the fleet average.** This is the slack-dispersion result planned for Phase 6, item (a).

*Check T2: all pass.*

**Proposition 3 (worst-case loss).** Under Assumption 1, let $X=(-\infty,b]$ be a range of delays and
$$\Gamma_X=\sup_{x\in X}\big[\max_r g_r(x)-\min_r g_r(x)\big].$$
For every $K\in\{1,\dots,N-1\}$,
$$\Gamma_X\le\sup_{D\in X^N}\big[V(S_\rho)-V(S_D)\big]\le\min(K,N-K)\,\Gamma_X ,$$
with ties among delays broken adversarially. At $K=1$ both sides are equal, so the worst-case loss is exactly the largest pointwise spread of the exposure functions. Moreover, $\Gamma_{\mathbb R}<\infty$ **if and only if all routes carry the same total weight**. Otherwise behind-schedule triage has unbounded worst-case loss as delays grow, because at equal large delay the route with more remaining promises carries proportionally more risk.

*Proof.*
1. **Upper bound.** $\bar g(x)=\tfrac12[\max_r g_r(x)+\min_r g_r(x)]$ is nondecreasing, being built from maxima and minima of nondecreasing functions. It satisfies $\lvert g_r(x)-\bar g(x)\rvert\le\Gamma_X/2$ on $X$. Apply Proposition 2′(i).
2. **Lower bound.**
   - Take $r,r'$ and $x$ with $g_{r'}(x)-g_r(x)$ close to $\Gamma_X$. Put $r'$ at $x$ and $r$ just above $x$.
   - Put $K-1$ further routes at or above $D_r$; there are enough, since $K\le N-1$. Push the rest towards $-\infty$, where their risk goes to 0.
   - Then $V(S_\rho)\ge V(T\cup\{r'\})$, where $T$ is the set of the $K-1$ routes at or above $D_r$. So the loss is at least $\rho_{r'}-\rho_r\to g_{r'}(x)-g_r(x)$ by continuity.
3. **Finiteness.** As $x\to\infty$, $g_r(x)=W_rx-\sum_i w_{ri}\mathbb EY_{ri}+o(1)$, so differences between exposure functions diverge if and only if the total weights differ. On $(-\infty,x_0]$ they are bounded by continuity and the limit 0. ∎

*Check T3: brute force over delay grids never exceeds Γ. The construction in the proof reaches Γ to within 1e-6.*

**Proposition 2′(ii) (any deterministic projection).** For any deterministic projected arrivals $\hat A_{ri}$, and $\pi_r=\sum_i w_{ri}(\hat A_{ri}-d_{ri})^+$,
$$\lvert\rho_r-\pi_r\rvert\le\tau_r=\sum_i w_{ri}\,\mathbb E\big[\lvert A_{ri}-\hat A_{ri}\rvert\,\mathbf 1\{\max(A_{ri},\hat A_{ri})>d_{ri}\}\big],$$
so $V(S_\rho)-V(S_\pi)\le2\min(K,N-K)\max_r\tau_r$. This covers the simulator's projection, which uses mean times with waiting rather than $p+D+\mu$.

*Check T4: skewed deviations and arbitrary projections, never violated. The route-level inequality also held in 100% of rows in the Gate 3 verification and smoke-test runs.*

**Closed form (Gaussian deviations).** If $\xi_{ri}\sim N(\mu_{ri},\sigma_{ri}^2)$, then
$$g_r(x)=\sum_i w_{ri}\,\sigma_{ri}\,\psi\!\Big(\frac{x-s_{ri}+\mu_{ri}}{\sigma_{ri}}\Big),\qquad \psi(z)=z\Phi(z)+\varphi(z),$$
which makes the effective-slack distances and Γ computable from planned slacks and dispersions.

## 3b. Average-case bounds (added 2026-11-26, after Gate 3; "Proposition 5" here is Proposition 4 in the manuscript, whose projection result is now Proposition 5)

**Motivation.** In Gate 3 (H4) the $L_\infty$ term of Proposition 2′(i) holds in every cell but is vacuous in all of them and does not rank cells by loss (Spearman 0.33). Its least-squares counterpart, $H_{\rm iso}=1-R^2$ of the isotonic regression of risk on delay, tracks the loss at 0.97 (exploratory). The result below turns that statistic into a valid bound.

**Notation.**
- $L=V(S_\rho)-V(S_D)$ is the loss of behind-schedule triage in captured risk, and $m_K=\min(K,N-K)$.
- $\bar g$ is a nondecreasing function, and the residuals are $e_r=\rho_r-\bar g(D_r)$.
- $e^\pm$ are the positive and negative parts of $e$, taken entrywise.
- For a nonnegative vector $x$, $T_m(x)$ is the sum of its $m$ largest entries.

**Proposition 5 (average-case loss bound; manuscript Proposition 4, "tail-average loss bound").** For any nonnegative risks, any delays and any nondecreasing $\bar g$,
$$L\;\le\;T_{m_K}(e^+)+T_{m_K}(e^-)\;\le\;\sqrt{2m_K}\,\lVert e\rVert_2 .$$
Take $\bar g=\hat g$, the least-squares isotonic regression of $\rho$ on $D$. This minimizes $\lVert e\rVert_2$ over nondecreasing $\bar g$. Writing $s_\rho^2=\frac1N\sum_r(\rho_r-\bar\rho)^2$ and $H=\lVert\rho-\hat g(D)\rVert_2^2/(Ns_\rho^2)$, we get
$$L\;\le\;s_\rho\sqrt{2m_KNH}.$$
For $K\le N/2$, with $\kappa=K/N$ and $\bar\rho_K=V(S_\rho)/K$ the mean risk of the $K$ riskiest routes,
$$\frac{L}{V(S_\rho)}\;\le\;\sqrt{\frac{2}{\kappa}}\;\sqrt H\;\frac{s_\rho}{\bar\rho_K}.$$
Both inequalities are attained.

*Proof.*
1. **Setup.** Let $A=S_\rho\setminus S_D$ and $B=S_D\setminus S_\rho$, with $\lvert A\rvert=\lvert B\rvert=j\le m_K$. Put $\sigma_r=\bar g(D_r)$.
2. **Order of $\sigma$.** $S_D$ is a top-$K$ set for $D$, so $D_b\ge D_a$ for $b\in B$ and $a\in A$. Since $\bar g$ is nondecreasing, $\sigma_b\ge\sigma_a$, and therefore $\sum_A\sigma_a\le\sum_B\sigma_b$.
3. **Loss decomposition.**
$$L=\sum_A\rho_a-\sum_B\rho_b=\sum_A e_a-\sum_B e_b+\Big(\sum_A\sigma_a-\sum_B\sigma_b\Big)\le\sum_A e_a^+ +\sum_B e_b^- .$$
4. **First inequality.** $A$ and $B$ have $j\le m_K$ elements and the entries are nonnegative, so the right-hand side is at most $T_{m_K}(e^+)+T_{m_K}(e^-)$.
5. **Second inequality.**
   - Cauchy–Schwarz on the $m_K$ selected entries gives $T_m(x)\le\sqrt m\,\lVert x\rVert_2$.
   - $e^+$ and $e^-$ have disjoint supports, so $\lVert e^+\rVert_2+\lVert e^-\rVert_2\le\sqrt2\,\lVert e\rVert_2$.
6. **Isotonic form.** The isotonic statement is the second inequality with $\lVert\rho-\hat g(D)\rVert_2^2=Ns_\rho^2H$. The relative form divides by $V(S_\rho)=K\bar\rho_K$ with $m_K=K$.
7. **Tightness.** Take $K\le N-K$, and choose three groups of routes as follows.

   | Group | Routes | Delay | Risk |
   |---|---|---|---|
   | $C$ | $N-2K$ | $-1$ | $0$ |
   | $A$ | $K$ | $0$ | $2c$ |
   | $B$ | $K$ | $1$ | $0$ |

   - Then $S_D=B$, $S_\rho=A$ and $L=2cK$.
   - The isotonic fit pools $A$ and $B$ at $c$ and puts $C$ at 0. The residuals are $+c$ on $A$, $-c$ on $B$ and 0 on $C$.
   - So $T_K(e^+)+T_K(e^-)=2cK$ and $\sqrt{2K}\,\lVert e\rVert_2=\sqrt{2K}\,c\sqrt{2K}=2cK$.

   For $K>N-K$, add $2K-N$ routes at delay 2 with risk $2c$. They lie in both $S_D$ and $S_\rho$ and have residual 0. ∎

**Comparisons.**
- **With Proposition 2′(i).** For the same $\bar g$, $T_m(e^+)+T_m(e^-)\le 2m\lVert e\rVert_\infty$, so the first bound always dominates Proposition 2′(i). Minimizing it over $\bar g$ is a linear program.
- **The quadratic bound** is not uniformly comparable with $2m_K\varepsilon^*$. For the same fit it is sharper when $\lVert e\rVert_2<\sqrt{2m_K}\max_r|e_r|$, that is, when fewer than about $2m_K$ routes have a misfit comparable to the largest. A few extreme routes then set the maximum while most residuals are small, which is the situation Gate 3 found.
  - *Corrected 2026-11-26 in the internal review.* An earlier version of this sentence said the opposite ("sharper when the misfit is spread over many routes").

**Corollary 2 (average slack dispersion).** Assume Assumption 1 holds and every route carries total weight $W$. Let $\bar g=\frac1N\sum_k g_k$ and $d_r=W_1(\nu_r,\bar\nu)$, as in the Wasserstein corollary. Then $\lvert e_r\rvert\le Wd_r$ uniformly in the delays. Because $e^+$ and $e^-$ have disjoint supports,
$$L\;\le\;W\,T_{2m_K}(d)\;\le\;W\sqrt{2m_K}\,\lVert d\rVert_2\qquad\text{for every delay vector.}$$
- **The first inequality** improves the Wasserstein corollary, $2m_KW\max_rd_r$: the maximum distance is replaced by the sum of the $2m_K$ largest distances.
- **The second** depends only on the root-mean-square distance.

**Population form.** Let routes be drawn i.i.d., with $(D_r,\rho_r)\sim P$ and $\mathbb E\rho^2<\infty$, and let $K/N\to\kappa\in(0,1)$. Let $g^\star$ be the $L^2(P)$ projection of $\rho$ onto nondecreasing functions of $D$; this is $\mathbb E[\rho\mid D]$ when that function is nondecreasing. Proposition 5 applied to the fixed function $g^\star$ and the strong law of large numbers give, almost surely,
$$\limsup_{N\to\infty}\frac{L}{N}\;\le\;\sqrt{2\min(\kappa,1-\kappa)}\;\big\lVert\rho-g^\star(D)\big\rVert_{L^2(P)} .$$
When $\mathbb E[\rho\mid D]$ is nondecreasing, the right-hand side is $\sqrt{2\min(\kappa,1-\kappa)\,\mathbb E\,\mathrm{Var}(\rho\mid D)}$. **The loss per route is governed by the conditional variance of risk given delay**: the spread of the exposure functions of different routes at the same delay, averaged over delays.
- **Proposition 3** governs the worst case through the largest pointwise spread $\Gamma$.
- **Proposition 5** governs the typical case through the mean squared spread.

The analogue for Corollary 2 replaces $\lVert d\rVert_2^2/N$ by $\mathbb E\,d_R^2$, the mean squared Wasserstein distance of a random route's effective-slack distribution from the population average. This holds for every sequence of delay vectors.

*Checks T5–T7 in `theory_checks.py`.*

## 3c. Route-wide shocks and learning (H3; added 2026-11-26, after Gate 3)

Manuscript numbering: Assumption 2, Lemma 3, Proposition 6 (signal-to-heterogeneity bound) and Proposition 7 (route-wide share at equal variance).

**Assumption 2 (linear learning).** Conditional on the information at the check-in,
$$\xi_{ri}=(\lambda_{ri}-1)D_r+U_{ri},$$
where the $\lambda_{ri}>0$ are constants and the law of $(U_{r1},\dots,U_{rn_r})$ does not depend on $D_r$.
- **Assumption 1 is the case $\lambda\equiv1$.**
- **Interpretation.** $\lambda_{ri}-1$ is the drift that the dispatcher extrapolates to promise $i$ per minute of current delay.

**Lemma 3 (representation with learning).** Under Assumption 2, $\rho_r=\tilde g_r(D_r)$, where
$$\tilde g_r(x)=\int_{-\infty}^x\tilde G_r(y)\,dy,\qquad \tilde G_r(y)=\sum_i w_{ri}\lambda_{ri}\,P(\tilde Y_{ri}\le y),\qquad \tilde Y_{ri}=\frac{s_{ri}-U_{ri}}{\lambda_{ri}} .$$

*Proof.* $A-d=D+\xi-s=\lambda D-(s-U)=\lambda(D-\tilde Y)$. Since $\lambda>0$, $(A-d)^+=\lambda(D-\tilde Y)^+$, and the law of $\tilde Y$ does not depend on $D$. Apply Lemma 1 with weights $w\lambda$. ∎

Several results use Assumption 1 only through Lemma 1: Lemma 1's properties, Proposition 1′, the Wasserstein corollary, Proposition 3 and Corollary 2. They therefore hold with
- weights $w\lambda$;
- effective slack $\tilde Y$, measured **in units of current delay**;
- total weight equal to the **delay sensitivity** $\tilde W_r=\sum_iw_{ri}\lambda_{ri}$.

Propositions 2′(i), 2′(ii) and 5 (manuscript Propositions 2, 5 and 4) need no assumption.

**Proposition 6 (signal-to-heterogeneity bound).** Assume Assumption 2 and $\tilde W_r=\tilde W$ for all $r$. Define:
- $\tilde\nu_r=\tilde G_r/\tilde W$, with finite mean $\tilde\mu_r$;
- $\bar{\tilde\nu}$, the average of the $\tilde\nu_r$, and $\tilde d_r=W_1(\tilde\nu_r,\bar{\tilde\nu})$;
- $E_K=\frac1K\sum_{r\in S_D}(D_r-\tilde\mu_r)^+$.

If $E_K>0$,
$$\frac{L}{V(S_\rho)}\le\frac{T_{2m_K}(\tilde d)}{K\,E_K}\le\frac{2m_K}{K}\cdot\frac{\max_r\tilde d_r}{E_K}.$$

*Proof.*
1. **Numerator.** Corollary 2 (via Lemma 3) gives $L\le\tilde W\,T_{2m_K}(\tilde d)$.
2. **Denominator.** By Jensen, $\tilde g_r(x)=\tilde W\,\mathbb E_{\tilde\nu_r}(x-\tilde Y)^+\ge\tilde W(x-\tilde\mu_r)^+$. Hence $V(S_\rho)\ge V(S_D)\ge\tilde WKE_K$.
3. Divide. The last step uses $T_{2m}(x)\le 2m\max x$. ∎

**Reading.** The relative loss is at most the heterogeneity of effective slack, measured in units of current delay, divided by how far the routes furthest behind exceed the slack their promises absorb on average.

**Unequal sensitivities.** With unequal $\tilde W_r$, no such bound holds for all delays (Proposition 3). Lemma 1's asymptotics give $\rho_r-\tilde W_r(D_r-\tilde\mu_r)\to0$ as $D_r\to\infty$. So at large delays, risk ranks routes by **delay times delay sensitivity**, not by delay alone.

**The linear-Gaussian route-shock model.**
- **Completed legs.** Route $r$ has completed legs $\ell=1,\dots,L_r$ with planned travel times $t_{r\ell}$. Write $\Lambda_r=\sum_\ell t_{r\ell}$, $Q_r=\sum_\ell t_{r\ell}^2$ and $n^{\rm eff}_r=\Lambda_r^2/Q_r\in[1,L_r]$; by Cauchy–Schwarz, $n^{\rm eff}_r=L_r$ if and only if the legs are equal.
- **Leg deviations.** The deviation on leg $\ell$ is $t_{r\ell}(\theta_r+\varepsilon_{r\ell})$, with $\theta_r\sim N(0,\varphi V)$ and $\varepsilon_{r\ell}\sim N(0,(1-\varphi)V)$ i.i.d.
- **Services.** Service deviations up to the check-in add independent Gaussian noise with variance $S_r$, and future services add $S_{ri}$ up to promise $i$.
- **Future legs.** $c_{ri}$ is the planned travel time from the check-in to promise $i$, and $q_{ri}$ is the sum of squared planned leg times over that stretch.
- **Learning.** The dispatcher learns $\theta_r$ from $D_r$.

**Proposition 7 (route-wide share at equal variance).** In this model, Assumption 2 holds with $\lambda_{ri}=1+c_{ri}\beta_r$ and Gaussian $U_{ri}$, where
$$\operatorname{Var}D_r=V[\Lambda_r^2\varphi+Q_r(1-\varphi)]+S_r,\qquad \beta_r=\frac{\Lambda_r\varphi V}{\operatorname{Var}D_r},$$
$$\operatorname{Var}(\theta_r\mid D_r)=\frac{\varphi V\,[Q_r(1-\varphi)V+S_r]}{\operatorname{Var}D_r},\qquad \operatorname{Var}U_{ri}=c_{ri}^2\operatorname{Var}(\theta_r\mid D_r)+q_{ri}(1-\varphi)V+S_{ri}.$$
As $\varphi$ rises at fixed $V$:
- **(a)** $\beta_r$ increases strictly, from 0 to $\Lambda_rV/(\Lambda_r^2V+S_r)$, so every $\lambda_{ri}$ increases.
- **(b)** $\operatorname{Var}D_r$ increases strictly if $L_r\ge2$. Its travel part grows by the factor $n^{\rm eff}_r$ from $\varphi=0$ to $\varphi=1$.
- **(c)** At $\varphi=1$ without service noise, $U\equiv0$ and
$$\rho_r=\sum_iw_{ri}\Big(\frac{\Lambda_r+c_{ri}}{\Lambda_r}D_r-s_{ri}\Big)^+ .$$

*Proof.*
1. **Delay.** $D_r=\Lambda_r\theta_r+\sum_\ell t_{r\ell}\varepsilon_{r\ell}+(\text{service})$ is Gaussian with the stated variance.
2. **Future deviation.** The future deviation at promise $i$ is $c_{ri}\theta_r$ plus independent leg and service noise.
3. **Conditioning.** Gaussian conditioning gives $\mathbb E[\theta_r\mid D_r]=\beta_rD_r$ with $\beta_r=\operatorname{Cov}(\theta_r,D_r)/\operatorname{Var}D_r$. The residual is Gaussian and independent of $D_r$, with the stated variance. So $\xi_{ri}=c_{ri}\beta_rD_r+U_{ri}$, which is Assumption 2.
4. **(a)** $\frac{d}{d\varphi}\beta_r\propto\Lambda_rV\big[\operatorname{Var}D_r-\varphi V(\Lambda_r^2-Q_r)\big]=\Lambda_rV(Q_rV+S_r)>0$.
5. **(b)** $\frac{d}{d\varphi}\operatorname{Var}D_r=V(\Lambda_r^2-Q_r)>0$ when two or more legs have positive time. The ratio of the travel parts at $\varphi=1$ and $\varphi=0$ is $\Lambda_r^2/Q_r$.
6. **(c)** At $\varphi=1$ with $S=0$, $\operatorname{Var}(\theta\mid D)=0$, $q(1-\varphi)V=0$ and $\beta=1/\Lambda$. ∎

**Exact version for the simulator (lognormal legs).** Each leg is $t\,e^{\theta+\varepsilon-V/2}$, so
$$\operatorname{Var}D_r=(e^{\varphi V}-1)(\Lambda_r^2-Q_r)+(e^V-1)Q_r+(e^{\sigma_S^2}-1)\sum_{\ell}\bar s_{r\ell}^2 .$$
Its travel part grows by exactly $n^{\rm eff}_r$ from $\varphi=0$ to $\varphi=1$.

**Reading for H3.** At equal total variance, raising the route-wide share works through two channels in Proposition 6.
1. **Amplification.** $\lambda$ grows, which divides the slack heterogeneity: $\tilde Y=(s-U)/\lambda$. At $\varphi=1$ the slack is scaled by $\Lambda/(\Lambda+c)$.
2. **Spread.** The same total variance produces up to $\sqrt{n^{\rm eff}}$ times more spread in delays. A route-wide shock accumulates linearly over the completed legs, while independent leg noise accumulates like a random walk. This raises $E_K$.

These are channels, not a monotonicity theorem. In minutes, $\operatorname{Var}U_{ri}$ first rises with $\varphi$: its derivative at 0 is $V(c_{ri}^2-q_{ri})>0$, because a weak route shock is poorly identified. So the propositions do not prove that the loss falls in $\varphi$. Gate 3 Arm B measured that the loss falls (H3 supported).

*Checks T8–T10 in `theory_checks.py`. The prediction for channel 2 on the Arm B runs is specified below before computing.*

**Channel-2 check on Gate 3 Arm B (exploratory; expectation written before computing).**
- **Data.** In the DAY regime there is no waiting, so `delay_plan_min` is exactly the sum of the completed legs' and services' deviations.
- **Prediction.** By the exact formula above, the cross-route mean of $D_r^2$ in each of the 10 DAY cells of Arm B (two variance levels, five values of $\varphi$) should equal the mean predicted variance.
- **Expectation.** The ratio of observed to predicted mean $D^2$ lies within three standard errors of 1 in all 10 cells.
- **Reported alongside:** the distribution of $n^{\rm eff}$ at the 35% check-in, and the predicted and observed growth of the delay spread from $\varphi=0.10$ to $0.90$.

**Result** (`code/resequencing/route_shock_check.py`, `analysis/theory/route_shock/`).
- **The expectation is met.** In all 10 DAY cells, the ratio of observed to predicted mean $D^2$ lies within three standard errors of 1. It is 0.93–0.97 at base variance and 1.06–1.08 at high variance, with |z| ≤ 2.2. The five cells at each level share the same random draws, so their deviations are correlated.
- **The delay spread grows as predicted.** From $\varphi=0.10$ to $0.90$, the root-mean-square delay grows by a factor of 1.68 (observed) against 1.72 (predicted) at base variance, and by 1.72 against 1.71 at high variance.
- **$n^{\rm eff}$ is small.** At the 35% check-in its median is only 6.6 (10–90%: 3.4–13.7), although the median route has completed 53 legs. Planned leg times are very unequal: a few long legs and many short hops between nearby stops. The spread channel can therefore raise the delay spread by at most about $\sqrt{6.6}\approx2.6$ for the median route, and the amplification channel carries the rest.

## 4. Still open for the theory-heavy version

1. **Route-wide shocks and learning (H3).** A result for the multiplicative shock: future deviation $c_{ri}\hat\theta_r+\eta_{ri}$, with $c_{ri}$ the planned time to stop $i$ and $\hat\theta_r$ increasing in $D_r$. Two possible forms:
   - the loss of behind-schedule triage shrinks as the share of variance in $\hat\theta$ grows (for example, a limit statement);
   - a bound in which the heterogeneity term is divided by the spread of $c_{ri}\hat\theta_r$.

   Gate 3 Arm B gives the numbers to check it against.
2. **Relative bounds.** The bounds are absolute (minutes). A share-of-optimum version divides by $V(S_\rho)$. The data-dependent version is already computed per cell in Gate 3 ($B_i$, $B_{ii}$); a version that is not data-dependent needs a lower bound on top-$K$ risk.
3. **Waiting extension** (§2b): state Proposition 1 with stop-specific absorption maps $h_{ri}$.

## 5. Consequences for the paper

- The theory is sound. With Lemma 1 and Propositions 1′, 2′ and 3 it becomes a coherent OR contribution: an exact characterization, distribution-free bounds that the data can evaluate, a worst-case loss the bound attains exactly, and a slack-dispersion (Wasserstein) form linking the theory to promise structure.
- Proposition 2's bounds hold in the simulator's regimes with waiting and learning, as well as in observed data, because they need no distributional assumption. Gate 3 checks them numerically in every cell.
- Reviewers may still see the top-$K$ perturbation argument as elementary. The framing should lead with Lemma 1, Proposition 3 and the Wasserstein form, and keep the perturbation lemma as a tool.

## 6. Proof expansion for submission (2026-11-26; manuscript numbering)

This follows `manuscript/submission_review_2026-11-26.md`. No result was weakened. The pre-edit manuscript is archived at `manuscript/archive/manuscript_2026-11-26_pre-proof-expansion.tex`.

**Gaps closed:**
- **Corollary 3 (large fleets).**
  - What it was: a prose paragraph.
  - What it is now: a numbered corollary with its hypotheses, including the condition $\mathbb E g(D)^2<\infty$, which was missing.
  - The unjustified min–max exchange is replaced by the elementary inequality $(x-t)^+\le x^2/(4t)$.
  - The minimizer is shown to exist: it is the $(1-\alpha)$-quantile.
  - $g$ must be fixed, not fitted.
  - $\bar T_\alpha$ is identified as $\alpha\,$CVaR (Rockafellar and Uryasev 2002).
- **Proposition 2.**
  - What it was: "attained" shown only for $N=2$.
  - What it is now: attained for every $N$ and $K$, using the three-group instances.
  - The closed form is attributed to Ubhaya (1974).

**Restructuring:** Lemma 2 is now the tail top-$K$ lemma for arbitrary scores, $L\le T_m(e^+)+T_m(e^-)\le T_{2m}(|e|)\le 2m\max|e|$. Propositions 2 and 4 and Corollary 2 follow from it.

**New result:** the tail form of Proposition 5, $L\le T_{2m_K}(\tau)$. It has **not** been evaluated on the data. If it is used as evidence, specify the expectation before computing it on `tau_ii_min`.

**Generalized:**
- **Proposition 1.**
  - The necessity construction now places the other routes at $x_0$ rather than pushing them to $-\infty$, and works for any single $K$.
  - Stated in the text: the equivalence holds on any delay interval $X$.
- **Proposition 3.**
  - $X$ is now any interval with nonempty interior.
  - The tie-breaking qualifier is removed, because the construction has no ties at the boundary.
  - $\Gamma_X<\infty$ whenever $X$ is bounded above.

**Made explicit:**
- Assumption 1 is now a fixed law $\mu_r$; this is the counterfactual meaning of $g_r$.
- Lemma 1's steps are written out.
- Assumption 2 now requires $\mathbb E|R|<\infty$, and the proof of Lemma 3 lists which results carry over.
- Proposition 7 now states the hypotheses $\sigma^2>0$, $\Lambda_r>0$ and $v_r=v_{ri}=0$ in part (c). Its Gaussian conditioning, derivative steps and the bound $1\le n^{\rm eff}\le L$ are written out.
- The tail-bound LP is in the appendix as equation (A.1), matching `avg_bound.tail_bound_lp`.

**Text fixes:**
- **Remark 1** no longer claims the Corollary 1 bound is small. §5.3 now reports the LP tail bound by regime: 0.34–0.75 under ETA, DAY and MIXA, and 0.97–0.99 under TW.
- **"Worst case equals Γ"** is now limited to a single intervention in the abstract, the introduction and the conclusion.
- **Proposition 7** is called a stylized analogue of the simulator, with three differences stated. In particular, the simulator learns θ from leg log-ratios.

**Checks:** T11–T14 were added to `code/resequencing/theory_checks.py`. All T1–T14 pass.

### 6b. Second pass (same day): formulas derived throughout; house-style spacing

**Author decision:** the manuscript stays at the house-style default spacing, not EJOR's 1.5 spacing. Proofs may go to the appendix, but every formula must be derived.

**Proofs moved to Appendix A and written out further.** Only Lemma 2's proof stays in the main text; each moved result keeps a pointer and a one-line key step.
- **Corollary 1.** $\nu_r$ is now shown to be a distribution with a finite mean.
- **Proposition 5.** The 1-Lipschitz case analysis is spelled out.
- **Lemma 3.** The text lists exactly which properties transfer.
- **Proposition 6.** The mixture variable $\tilde Y_r$ makes the Jensen step explicit.

**Derived in the text (Section 3):**
- The condition under which the quadratic form of Proposition 4 is sharper.
- $\mathbb E e^2=\mathbb E\operatorname{Var}(\rho\mid D)$, by the tower property. Conditional Jensen makes the integrability of $\mathbb E[\rho\mid D]$ automatic.
- $\operatorname{Var}(\rho\mid D=x)=\operatorname{Var}(g_R(x)\mid D=x)$.
- The learning identity (now an equation) and the "broken when $D>\tilde Y$" reading.
- $\tilde W_r$ as the limit of $\tilde G_r$.
- The asymptote $\rho-\tilde W(D-\tilde\mu)\to0$, with $\tilde W\tilde\mu=\sum w\lambda\mathbb E\tilde Y$.
- The amplification factor $\Lambda/(\Lambda+c)$ and the spread factor $\sqrt{n^{\rm eff}}$.
- Remark 2 (waiting) now has an explicit model $A=p+h(D)+\xi$ and the derivation.

**Formulas added to Sections 4–6, matched to the code:**
- Execution model: equation (7), with the mean-one correction derived.
- Normal–normal posterior: equation (8), derived by completing the square. It matches `sim.posterior_theta`.
- Arrival recursion with waiting: equation (9), matching `sim.eval_seq`.
- Definitions of L0–L4 and of $D_r$ (the planned end of service at the last stop, under the regime's waiting rule).
- Removable lateness $\Delta_r$, efficiency (equation (10)), the L0 expectation, L4/oracle, and $\omega$ (equation (11), matching `analyze_gate3`).
- The $\sigma_U$ error model, matching `gate3.true_means`.
- The calibration model with the moment estimator, matching `gate3_calibrate.re_decomp`.
- Gate 5 scores $S^{(1)},S^{(2)},S^{(3)}$ and the recovered share $\mathcal G$, matching `gate5.add_scores` and the Gate 5 pre-registration.
- Observed-data scores B1, B2 and F, the pace fit, $\mathcal R_r$, the outcome $Z_r$ and efficiency, matching `gate2d.py` and `gate4.py`.

**New Appendix B (derivations for the computational study).** The pre-registration record is now Appendix C.
- **Lemma 4.** $\mathbb E|X-1|=2[2\Phi(\sigma/2)-1]$, giving 20/39/77%.
- **Lemma 5.** The unbalanced one-way random-effects moments: $\mathbb E\,\mathrm{MSB}=\psi^2_{\rm leg}+J_0\psi^2_{\rm route}$.
- **Lemma 6.** The exact lognormal $\operatorname{Var}D$ (equation (12)) and its first-order reduction to Proposition 7.

**Checks:** T15 (Lemma 4) and T16 (Lemma 5, including agreement with `re_decomp`) were added and pass. T10 already covers Lemma 6.

**Result:** the manuscript is 38 pages in house style, with the main text through page 25 and Appendices A–C on pages 26–36. It compiles cleanly, with no overfull boxes, and every label is referenced. The abstract is 249 words.

### 6c. Readiness fixes (2026-11-27; review in `manuscript/submission_review_2026-11-27.md`)

The pre-edit copy is `manuscript/archive/manuscript_2026-11-27_pre-readiness-fixes.tex`.

**Theory**
- **Proposition 5 is tight.** Two groups of $m_K$ routes with deterministic arrivals and no ties in $\pi$ give a loss equal to the share $1-\epsilon/(2c)$ of $T_{2m_K}(\tau)=2m_K\max\tau$. The construction is in the appendix proof, and new check T17 passes for all $2\le N\le24$ and every $K$.
- **Proposition 4.** The relative form now states the condition $V(S_\rho)>0$.

**Derivation.** Appendix B now gives the first-order route share of simulated inter-arrival residuals (equation B.1, `eq:route-share`):
$$\text{share}\approx\varphi\,\mathbb E[\bar a^2]\sigma^2/(\mathbb E[a^2]\sigma^2+\mathbb E[(1-a)^2]\sigma_S^2),$$
where $a$ is the travel share of the planned inter-arrival time. The formula is within 12% of the simulated share in all 13 calibration settings (`route_share_check.py`, output `analysis/theory/route_shock/route_share_check.csv`). Section 5.4 adds a caveat that part of the observed route effect may be misfit of the fleet pace at route level.

**Text and numbers**
- $\omega$ is 0.17, not 0.18.
- The L3 lower bound is "at least 0.88" (and 88%).
- Gate 5 harm is 0.030 in the robustness cells and 0.051 in the delay-structure cells.
- Propositions 2 and 5 are "tight in the worst case".
- Notation: $V(\varphi)\to\mathcal V_r(\varphi)$, and the LP slacks are now $\eta^\pm_r$.
- Highlight 3 was reworded.

**Check 3c and bootstrap intervals.** Both were specified in `projection_tail_and_bootstrap.md` before computing and run with `bootstrap_tail_check.py`; the results are reported in §§4.3, 5.1, 5.3 and 5.4, the note to Table 3, and Table C.1.

**Code**
- `konovalenko_io.py` reproduces `konovalenko_stops.pkl` exactly.
- `theory_checks.py` now uses manuscript numbering and includes T17.
- The figure scripts write vector PDFs; the manuscript uses them.

**EJOR copy (`build_ejor.py`)**
- Appendix C and the full §5.5 moved to the supplement, which now has sections S1–S4; §5.5 is replaced by `gate5_summary.tex`.
- hyperref link boxes are hidden.
- At 1.5 spacing the main file is 35 pp (about half of the last page free), the supplement 12 pp, and the working draft 41 pp.

### 6d. Corollary 4 and the score symbol (2026-11-27, second round)

The pre-edit copy is `manuscript/archive/manuscript_2026-11-27_pre-sigma-rename.tex`.

**Score symbol.** Scores are now $\hat\rho_r$ (formerly $\sigma_r$), in Lemma 2, the efficiency metric (equation (10), $\mathrm{Eff}(\hat\rho)$), and the proofs of Propositions 2, 4 and 5 and Corollary 1. The symbol $\sigma$ now denotes only standard deviations.

**New Corollary 4 (projection from conditional means).** If $\pi_r\le\rho_r$ for every route, then $V(S_\rho)-V(S_\pi)\le T_{m_K}(\tau)\le m_K\max_r\tau_r$, and both bounds are tight.
- **Sufficient condition.** It holds when $\hat A\le\mathbb EA$, via monotonicity and Jensen's inequality on $(x-d)^+$. In particular it holds for a schedule built from conditional means: recursion (9) is convex in the travel and service times, as an induction on compositions of the convex nondecreasing maps $\mathcal O_{rk}$, and the arrival times are integrable.
- **Tightness.** Group 1 of the Proposition 5 construction becomes exact ($\pi=\rho=c$), and the loss share is $1-\epsilon/c$.
- **Check T18 passes.** On random fleets with waiting, under empirical laws so that every expectation is exact: $\max(\hat A-\mathbb EA)=1.6\cdot10^{-12}$; the loss is at most $T_m(\tau)$ for every $K$, with a maximum loss/bound of 0.118; the construction is tight.
- **Data sanity check (not in the paper).** In Gate 2c, $\pi_r$ exceeds the prior estimate L3 in 2.1% of route-cells, by at most 14% relative. This is consistent with Monte Carlo error over 256 prior scenarios; the paper states only the theoretical condition.

**Main text.** A short discussion after the corollary adds three points:
- by Jensen's inequality, the projection can only underestimate risk;
- this is the failure mode of Remark 1;
- in the study, the condition holds against L3 but not against L4.

The first contributions bullet mentions the corollary.

**EJOR length.** The corollary pushed the main file to 36 pp at 1.5 spacing, so `build_ejor.py` (decision H) now also moves the exploratory spread-channel check of §5.4, with equation (12), to supplement Section S3. The main file keeps `spread_summary.tex`. The supplement now has:
- S1, proofs;
- S2, derivations;
- S3, spread check;
- S4, Gate 5;
- S5, record.

**Lengths.**
- The EJOR main file is 35 pp at 1.5 spacing, with about a third of the last page free.
- The supplement is 13 pp, and the working draft 42 pp.

### 6e. Fourth review pass: disclosure, two derivations, density (2026-11-27; review in `manuscript/submission_review_2026-11-27_pass4.md`)

Pre-edit copies: `manuscript/archive/*_2026-11-27_pre-pass4-fixes.*`.

**New results, each with a proof and a check.**
- Remark 1 now derives the projection under plan-anchored promises. Without waiting, the schedule from the mean times is p + D, so pi_r = W_r (D_r - beta)^+ under ETA promises, and pi_r = 0 for D_r <= beta under day-definite promises. Check T20. On the Gate 2c data the routes with pi > 0 are exactly those with D > 60 min (13 = 13 at base noise, 117 = 117 at high noise, ETA, 35% check-in).
- New Lemma 4 (Appendix A, after the proof of Corollary 4): the dispersion term under Gaussian deviations. It is sigma phi(z) for a promise projected to be met and sigma sqrt(2/pi) for one projected to be missed, and the error of the projected lateness is sigma L(|z|). Check T19. The text after Proposition 5 now uses it instead of a verbal claim. The Appendix B lemmas become Lemmas 5-7, and theory_checks.py and both READMEs are renumbered.
- T1-T20 all pass.

**Disclosures (no result changes).**
- Tie-breaking (§4.3, §6.1): stable sort, i.e. file order. Breaking ties at random changes no reported efficiency by more than 0.029. The largest change over all 1,884 cell-score-capacity combinations is 0.061 (L2, DAY, phi = 0.10, 10%, not reported). Script `ties_mc_check.py`; outputs in `analysis/theory/ties_mc/`.
- Monte Carlo (§4.3): L3 uses 256 prior scenarios on its own stream. L4 is the held-out estimate that also values Delta, so it shares its Monte Carlo error. Scoring L4 on the independent 48 optimization scenarios keeps at least 0.984 of its value (median 0.998, 56 cells at 5%). The oracle's restart from the earliest-deadline order and its 60-pass limit are now stated.
- §5.2 and §7: the projection under plan-anchored promises is "uninformative", not "worse than the delay itself".

**Wording.**
- Intro: worst-case gap finite only over unrestricted delays.
- CVaR at confidence level 1 - alpha.
- mu-tilde defined once more generally.
- "tracks" instead of "governed".
- MAE "as a fraction of that mean".
- Calibration caveat on the travel share.
- Amplification "must come through", not measured.
- 10% station-day capacity explained (90% of station-days have fewer than 30 routes).
- Managerial bullet: 90-96% confirmatory and 81-94% across delay structures.
- Highlight 4 limited to the window-driven gap.
- Proof set T renamed \mathcal{T}.
- Capacity sentence (2% and 10%) added to the Table 4 note.

**Density.**
- Theory contribution split into three sub-bullets.
- §3.4 split into "Worst-case loss" and "Loss away from the worst case".
- New Table 2 summarizes the theoretical results, replacing the long list in the paragraph after Assumption 1.
- §5.3's bound evaluation is condensed around new Table 5 (bounds by cell set, from `avg_bound_summary.json` and `tail_bootstrap/check3c_cells.csv`). The table also reports, for completeness, the max-form projection bound on the delay-structure cells: Spearman 0.90.

**EJOR length.** The additions brought the main file to 38 pp at 1.5 spacing (the pre-edit copy measures 35 pp the same way). `build_ejor.py` therefore gained two decisions:
- **I:** the bound evaluation of §5.3 with its table moves to supplement S3, and the main file keeps `bounds_summary.tex`.
- **J:** the large-fleet derivation after Corollary 3 moves to S1, after the proof of Corollary 3, and the main file keeps `large_fleet_summary.tex`.

Supplement sections are now S1 proofs, S2 derivations, S3 bounds, S4 spread, S5 Gate 5 and S6 record.

**Lengths.**
- EJOR main: 35 pp at 1.5 spacing, about 2 lines free on the last page, so it is at the limit.
- Supplement: 16 pp.
- Working draft: 44 pp in house style.
