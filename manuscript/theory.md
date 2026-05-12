# SwarmLineage-OT Theory Sketch

Let each virtual cell agent carry latent state \(z_i(t)\), fate belief \(q_i(t)\), age or cell-cycle state \(a_i(t)\), density/resource score \(r_i(t)\), optional spatial coordinate \(s_i(t)\), and memory exposure \(m_i(t)\). The finite-agent dynamics are modelled as a birth-death-interacting SDE:

\[
dz_i =
\left[
v_{\theta}^{\mathrm{intrinsic}}(z_i,t)
+ v^{\mathrm{OT}}(z_i,t)
+ v_{\theta}^{\mathrm{swarm}}(z_i,\mathcal{N}_i,r_i,q_i)
+ v_{\theta}^{\mathrm{CCI}}(z_i,\mathcal{C}_i)
- \nabla U_{\mathrm{manifold}}(z_i)
\right]dt
+ \sigma_{\theta}(z_i,t,H_i,r_i,a_i)dW_i .
\]

Birth and death hazards are
\[
\lambda_i^{b}=\mathrm{softplus}(f_\theta^b(z_i,t,a_i,r_i,g_i,H_i,\mathcal{C}_i)),
\quad
\lambda_i^{d}=\mathrm{softplus}(f_\theta^d(z_i,t,r_i,\mathcal{C}_i)).
\]

The phenomenological memory field for fate \(k\) is
\[
\partial_t h_k = \mathrm{deposit}_k(\rho_t,q)-\eta h_k + D_h \Delta h_k .
\]

Here “pheromone” is a complex-systems memory variable, not a literal biological pheromone. In biological interpretation it can represent niche signal, secreted factor, ECM remodelling or historical flow.

## Mean-Field Limit

Under exchangeability, bounded hazards and Lipschitz drift/interactions, the empirical measure \(\rho_t^N=N^{-1}\sum_i\delta_{z_i(t)}\) formally converges to a density \(\rho_t(z)\) satisfying

\[
\partial_t \rho_t(z)
=
-\nabla\cdot\left(\rho_t v_\theta(z,t,\rho_t,h_t)\right)
+ \nabla\cdot\left(D_\theta(z,t,\rho_t)\nabla \rho_t\right)
+ [\beta_\theta(z,t,\rho_t)-\delta_\theta(z,t,\rho_t)]\rho_t .
\]

Dynamic or unbalanced OT estimates a minimal-cost pseudo-lineage coupling between observed snapshots. SwarmLineage-OT does not redefine the OT cost as its main claim; it learns an executable control field whose finite-agent birth-death-diffusion process approximates the OT-inferred pseudo-lineage while exposing local rules and perturbation knobs.

## Testable Proposition

Assume the OT teacher couplings are consistent with a reference developmental process, hazards are bounded, interaction kernels are Lipschitz, and terminal fate sets are separated by margin \(\gamma>0\) in latent space. If snapshot, barycentric and growth losses converge to zero over observed times, then the simulated terminal fate composition converges to the teacher terminal composition. More concretely, the Wasserstein/Fisher-Rao discrepancy between simulated and teacher terminal fate mass is bounded by a constant times

\[
L_{\mathrm{snapshot}}^{1/2}
+ L_{\mathrm{barycentric}}^{1/2}
+ L_{\mathrm{growth}}^{1/2}
+ \epsilon_{\mathrm{teacher}}
+ O(N^{-1/2}),
\]

where \(\epsilon_{\mathrm{teacher}}\) is pseudo-lineage error and \(N^{-1/2}\) is finite-agent sampling error. The proof sketch couples simulated particles to teacher barycentric descendants, applies stability of continuity equations with bounded source terms, and projects terminal mass onto separated fate basins.

