# Package index

## Model Fitting

Fit DBN models to temporal network data

- [`dbn()`](https://netify-dev.github.io/dbn/reference/dbn.md) : Dynamic
  Bilinear Network Analysis
- [`dbn_static()`](https://netify-dev.github.io/dbn/reference/dbn_static.md)
  : Static DBN MCMC
- [`dbn_dynamic()`](https://netify-dev.github.io/dbn/reference/dbn_dynamic.md)
  : Dynamic DBN MCMC
- [`dbn_lowrank()`](https://netify-dev.github.io/dbn/reference/dbn_lowrank.md)
  : DBN Low-rank Model
- [`dbn_lowrank_accurate()`](https://netify-dev.github.io/dbn/reference/dbn_lowrank_accurate.md)
  : DBN Low-rank Model (Accurate)
- [`dbn_hmm()`](https://netify-dev.github.io/dbn/reference/dbn_hmm.md) :
  DBN HMM Model
- [`dbn-options`](https://netify-dev.github.io/dbn/reference/dbn-options.md)
  : DBN Package Options
- [`dbn-methods`](https://netify-dev.github.io/dbn/reference/dbn-methods.md)
  : S3 Methods for DBN Objects

## Simulation

Generate synthetic network data

- [`simulate`](https://netify-dev.github.io/dbn/reference/simulate.md) :
  Simulate Data from DBN Models
- [`simulate_test_data()`](https://netify-dev.github.io/dbn/reference/simulate_test_data.md)
  : Simulate Simple Test Data
- [`simulate_static_dbn()`](https://netify-dev.github.io/dbn/reference/simulate_static_dbn.md)
  : Simulate from Static DBN Model
- [`simulate_dynamic_dbn()`](https://netify-dev.github.io/dbn/reference/simulate_dynamic_dbn.md)
  : Simulate from Dynamic DBN Model
- [`simulate_piecewise_dbn()`](https://netify-dev.github.io/dbn/reference/simulate_piecewise_dbn.md)
  : Simulate Data from Piecewise DBN
- [`simulate_lowrank_dbn()`](https://netify-dev.github.io/dbn/reference/simulate_lowrank_dbn.md)
  : Simulate from Low-Rank DBN Model
- [`simulate_hmm_dbn()`](https://netify-dev.github.io/dbn/reference/simulate_hmm_dbn.md)
  : Simulate from HMM DBN Model

## Posterior Inference

Posterior prediction, draws, and forecasting

- [`posterior_predict_dbn()`](https://netify-dev.github.io/dbn/reference/posterior_predict_dbn.md)
  : Generate posterior predictive samples
- [`derive_draws()`](https://netify-dev.github.io/dbn/reference/derive_draws.md)
  : Derive new quantities from posterior draws
- [`regime_probs()`](https://netify-dev.github.io/dbn/reference/regime_probs.md)
  : Extract regime probabilities for HMM models

## Impulse Response Analysis

Compute and analyze network shocks

- [`compute_irf()`](https://netify-dev.github.io/dbn/reference/compute_irf.md)
  : Compute Network-Level Impulse Response Functions
- [`impulse_response_const()`](https://netify-dev.github.io/dbn/reference/impulse_response_const.md)
  : Compute impulse response for constant A,B matrices
- [`impulse_response_dynamic()`](https://netify-dev.github.io/dbn/reference/impulse_response_dynamic.md)
  : Compute impulse response for time-varying A,B matrices
- [`build_shock()`](https://netify-dev.github.io/dbn/reference/build_shock.md)
  : Build Shock Matrix for IRF Analysis
- [`debug_irf()`](https://netify-dev.github.io/dbn/reference/debug_irf.md)
  : Debug IRF

## Network Analysis

Summarize and analyze network structure

- [`network_summary()`](https://netify-dev.github.io/dbn/reference/network_summary.md)
  : Network-level posterior summary
- [`edge_prob()`](https://netify-dev.github.io/dbn/reference/edge_prob.md)
  : Posterior edge probability
- [`theta_summary()`](https://netify-dev.github.io/dbn/reference/theta_summary.md)
  : Summarize Theta over posterior draws
- [`theta_credible()`](https://netify-dev.github.io/dbn/reference/theta_credible.md)
  : Posterior credible intervals for Theta
- [`param_summary()`](https://netify-dev.github.io/dbn/reference/param_summary.md)
  : Summarize scalar parameters
- [`latent_summary()`](https://netify-dev.github.io/dbn/reference/latent_summary.md)
  : Summarize latent means (M arrays)
- [`tidy_dbn_lowrank()`](https://netify-dev.github.io/dbn/reference/tidy_dbn_lowrank.md)
  : Tidy Extractor for Low-Rank Factor Paths
- [`theta_slice()`](https://netify-dev.github.io/dbn/reference/theta_slice.md)
  : Extract Theta slices from posterior draws
- [`net_snapshot()`](https://netify-dev.github.io/dbn/reference/net_snapshot.md)
  : Network Snapshot
- [`dyad_path()`](https://netify-dev.github.io/dbn/reference/dyad_path.md)
  : Plot Dyad Trajectory
- [`role_trajectory()`](https://netify-dev.github.io/dbn/reference/role_trajectory.md)
  : Track Role Evolution
- [`get_group_influence()`](https://netify-dev.github.io/dbn/reference/get_group_influence.md)
  : Extract Group Influence Trajectories
- [`compare_group_influence()`](https://netify-dev.github.io/dbn/reference/compare_group_influence.md)
  : Compare Group Influences

## Network Statistics

Compute descriptive network statistics

- [`stat_density()`](https://netify-dev.github.io/dbn/reference/stat_density.md)
  : Network Statistic: Density
- [`stat_in_degree()`](https://netify-dev.github.io/dbn/reference/stat_in_degree.md)
  : Network Statistic: In-Degree
- [`stat_out_degree()`](https://netify-dev.github.io/dbn/reference/stat_out_degree.md)
  : Network Statistic: Out-Degree
- [`stat_reciprocity()`](https://netify-dev.github.io/dbn/reference/stat_reciprocity.md)
  : Network Statistic: Reciprocity
- [`stat_transitivity()`](https://netify-dev.github.io/dbn/reference/stat_transitivity.md)
  : Network Statistic: Transitivity

## Visualization

Plot models, diagnostics, and results

- [`plot_theta()`](https://netify-dev.github.io/dbn/reference/plot_theta.md)
  : Plot Theta Heatmap
- [`plot_trace()`](https://netify-dev.github.io/dbn/reference/plot_trace.md)
  : Plot Parameter Trace Plots
- [`plot_regime_probs()`](https://netify-dev.github.io/dbn/reference/plot_regime_probs.md)
  : Plot Regime Probabilities
- [`plot_ppc_density()`](https://netify-dev.github.io/dbn/reference/plot_ppc_density.md)
  : Posterior predictive density plot
- [`plot_ppc_ecdf()`](https://netify-dev.github.io/dbn/reference/plot_ppc_ecdf.md)
  : Posterior predictive ECDF plot

## Model Methods

S3 methods for dbn objects

- [`print(`*`<dbn_derived>`*`)`](https://netify-dev.github.io/dbn/reference/print.dbn_derived.md)
  : Print method for derived quantities
- [`print(`*`<dbn_irf>`*`)`](https://netify-dev.github.io/dbn/reference/print.dbn_irf.md)
  : Print IRF
- [`print(`*`<dbn_ppd>`*`)`](https://netify-dev.github.io/dbn/reference/print.dbn_ppd.md)
  : Print method for posterior predictive distribution
- [`summary_dbn()`](https://netify-dev.github.io/dbn/reference/summary_dbn.md)
  : Summary method for DBN objects
- [`plot(`*`<dbn_irf>`*`)`](https://netify-dev.github.io/dbn/reference/plot.dbn_irf.md)
  : Plot IRF
- [`plot_dbn()`](https://netify-dev.github.io/dbn/reference/plot_dbn.md)
  : Plot method for DBN objects
- [`plot_group_influence()`](https://netify-dev.github.io/dbn/reference/plot_group_influence.md)
  : Plot Group Influence Profile
- [`plot_ppc_density()`](https://netify-dev.github.io/dbn/reference/plot_ppc_density.md)
  : Posterior predictive density plot
- [`plot_ppc_ecdf()`](https://netify-dev.github.io/dbn/reference/plot_ppc_ecdf.md)
  : Posterior predictive ECDF plot
- [`plot_regime_probs()`](https://netify-dev.github.io/dbn/reference/plot_regime_probs.md)
  : Plot Regime Probabilities
- [`plot_theta()`](https://netify-dev.github.io/dbn/reference/plot_theta.md)
  : Plot Theta Heatmap
- [`plot_trace()`](https://netify-dev.github.io/dbn/reference/plot_trace.md)
  : Plot Parameter Trace Plots
- [`compare_dbn()`](https://netify-dev.github.io/dbn/reference/compare_dbn.md)
  : Compare Multiple DBN Models
- [`compare_blocks()`](https://netify-dev.github.io/dbn/reference/compare_blocks.md)
  : Compare Influence Matrices Across Blocks
- [`check_convergence()`](https://netify-dev.github.io/dbn/reference/check_convergence.md)
  : Check MCMC Convergence

## Utilities

Configuration and computational utilities

- [`estimate_memory()`](https://netify-dev.github.io/dbn/reference/estimate_memory.md)
  : Estimate Memory Usage for Dynamic DBN
- [`get_dbn_threads()`](https://netify-dev.github.io/dbn/reference/get_dbn_threads.md)
  : Get the current number of threads used by dbn
- [`set_dbn_threads()`](https://netify-dev.github.io/dbn/reference/set_dbn_threads.md)
  : Set the number of threads used by dbn
- [`ffbs_theta_struct_5arg_cpp()`](https://netify-dev.github.io/dbn/reference/ffbs_theta_struct_5arg_cpp.md)
  : Structured FFBS for Theta

## Data

Example datasets

- [`example_data`](https://netify-dev.github.io/dbn/reference/example_data.md)
  : Example Dynamic Bilinear Network Data
