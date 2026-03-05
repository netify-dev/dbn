# DBN Package Options

options controlling R vs C++ backend selection.

## Performance Options

- dbn.use_cpp_update_ab:

  Logical (default: TRUE). Use C++ for A/B updates in HMM.

- dbn.use_cpp_build_f:

  Logical (default: TRUE). Use C++ for design matrix F in lowrank.

- dbn.use_batch_ffbs:

  Logical (default: TRUE). Use batch FFBS updates.

- dbn.use_cpp_stability:

  Logical (default: TRUE). Use C++ for stability functions.

- dbn.use_ffbs_dlm_cpp:

  Logical (default: TRUE). Use C++ FFBS for DLMs.

- dbn.use_ffbs_cpp:

  Logical (default: TRUE). Use C++ time-varying FFBS.

- dbn.use_cpp_ranklik:

  Logical (default: TRUE). Use C++ rank likelihood sampling.

## Setting Options

    options(dbn.use_cpp_ranklik = FALSE)
