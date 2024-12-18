# README

## Overview

This repository provides the implementation of the TT-IPP and MC-IPP algorithms introduced in the paper **"Inexact Proximal Point Algorithms for Zeroth-Order Global Optimization"** ([arXiv:2412.11485](https://arxiv.org/abs/2412.11485)).

The **`Main_test`** function is the primary entry point to test and compare several optimization algorithms, including TT-IPP and MC-IPP, on various benchmark test functions. The folder `optimization_algorithms` contains the following key implementations:

- **TT-IPP**: To utilize the TT-IPP algorithm for a function in tensor train format, run:

  ```matlab
  TT_IPP_optimization_fun
  ```

- **MC-IPP**: To utilize the MC-IPP algorithm, run:

  ```matlab
  MC_IPP_optimization_fun
  ```

Definitions of the benchmark test functions are located in `choose_example.m`. The TT-IPP algorithm relies on the [TT-Toolbox](https://github.com/oseledets/TT-Toolbox/), which is included in the `Aux_functions` folder.

---

## Citation

If this code or approach contributes to your research, please cite the following:

```bibtex
@misc{zhang2024inexactproximalpointalgorithms,
      title={Inexact Proximal Point Algorithms for Zeroth-Order Global Optimization},
      author={Minxin Zhang and Fuqun Han and Yat Tin Chow and Stanley Osher and Hayden Schaeffer},
      year={2024},
      eprint={2412.11485},
      archivePrefix={arXiv},
      primaryClass={math.OC},
      url={https://arxiv.org/abs/2412.11485},
}
```

