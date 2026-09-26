# HubbardSuperconductorOrderParameter — manual

Source: DiagHam wiki page `HubbardSuperconductorOrderParameter`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [HubbardSuperconductorOrderParameter](../HubbardModels/HubbardSuperconductorOrderParameter.md)
HubbardSuperconductorOrderParameter computes the order parameter matrix $\hat{P}_{N+2}c^\dagger_{i,\sigma} c^\dagger_{j,\sigma'} \hat{P}_{N} c_{j,\sigma'} c_{i,\sigma}$ where  $\hat{P}_{N+2}=\frac{1}{n_{N+2}}\sum_{\alpha=1}^{n_{N+2}}|\Psi_{N+2,\alpha}\rangle\langle\Psi_{N+2,\alpha}|$ is the projector on the $n_{N+2}$ dimensional subspace (called the left subspace) of states with $N+2$ particles. Similarly,  $\hat{P}_{N}=\frac{1}{n_{N}}\sum_{\beta=1}^{n_{N}}|\Psi_{N,\beta}\rangle\langle\Psi_{N,\beta}|$ is the projector on the $n_{N}$ dimensional subspace (called the right subspace) of states with $N$ particles.

If the left and right subspaces are each made of a single state, one can just use

*\$PATHTODIAGHAM/build/FTI/src/Programs/HubbardModels/HubbardSuperconductorOrderParameter --left-state fermions_kitaev_heisenberg_gutzwiller_stripe_n_6_x_6_t_0_tK_1_j1_0_j2_0.666667.0.vec --right-state fermions_kitaev_heisenberg_gutzwiller_stripe_n_4_x_6_t_0_tK_1_j1_0_j2_0.666667.0.vec*

The code creates a text file fermions_kitaev_heisenberg_gutzwiller_stripe_n_6_x_6_t_1_tK_1_j1_0.49_j2_0.51.0.orderparam.dat whose name is deduced from the left state file name. It contains the order parameter for any possible indices i.e.

```text
# <Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R> with sigma,sigma' = 0 (down) or 1 (up)
# i j sigma sigma' |<Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R>|^2
0 1 0 0  3.4319578038925e-35
0 1 0 1 0.015539340372598
0 1 1 0 0.015539340372598
0 1 1 1 5.0494284504268e-33
...
```

Site labels depend on the geometry (for example, see here *(wiki page HubbardKitaevHeisenbergModel#Stripe configuration)* for a the stripe honeycomb model).

If the subspaces are spanned by more than one state, the --degenerate-rightstates and/or --degenerate-leftstates have to be used. In the following example

*\$PATHTODIAGHAM/build/FTI/src/Programs/HubbardModels/HubbardSuperconductorOrderParameter --left-state fermions_kitaev_heisenberg_gutzwiller_stripe_n_6_x_6_t_1_tK_0_j1_0.666667_j2_0.0.vec --degenerate-rightstates basis_n_4_x_6_t_1_tK_0_j1_0.666667_j2_0.dat*

the left subspace is made of a single state while the right subspace basis_n_4_x_6_t_1_tK_0_j1_0.666667_j2_0.dat is an ASCII file that lists the two states than span the right subspace

```text
fermions_kitaev_heisenberg_gutzwiller_stripe_n_4_x_6_t_1_tK_0_j1_0.666667_j2_0.0.vec
fermions_kitaev_heisenberg_gutzwiller_stripe_n_4_x_6_t_1_tK_0_j1_0.666667_j2_0.1.vec
```

The resulting file is 

```text
# <Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R> with sigma,sigma' = 0 (down) or 1 (up)
# i j sigma sigma' |<Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R>|^2
0 1 0 0  2.6222915549751e-35
0 1 0 1 0.016085277676507
0 1 1 0 0.016085277676507
0 1 1 1 9.6782473831841e-33
```

Beyond the fourth column are the order parameter matrix eigenvalues. The number of eigenvalues is equal to th smallest subspace dimension.
