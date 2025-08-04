# FLASH-D
In this work we implement in C++ hardware accelerator for attention inference following FlashAttention-2 and our proposed Flash-D approach, that computes the output vector using as weighted sum of its previous value and an incoming value vector. The weight is computed incrementally and its formulation helps eliminate max computation needed for ensuring numerical stability and hide the final vector division operation in a sigmoid function calculation. In order to evaluate power metrics we run inference using Microsoft's [PromptBench](https://github.com/microsoft/promptbench.git) framework. More specifically, we run pytorch models from [huggingface](https://huggingface.co/) and extracted inter layer results for different prompts to use as inputs in ```main.cc```.

Most of the floating-point functionality utilizes the [Fast-Float4HLS](https://github.com/ic-lab-duth/Fast-Float4HLS.git) library, publicly available on github.

## Repository Hierarchy

This repository is organized as follows:

```bash
.
├── src
│   ├── defines.h
│   ├── flash_atten.h
│   ├── file_io.h
│   ├── fp_arithm.h
│   ├── logging.h
│   ├── main.cc
│   ├── math_ops.h
│   └── reduction.h
│
├── utils
│   ├── gen_pwl_coeff.py
│   └── pack.py
│
├── LICENSE
├── README.md
└── setup.sh
```

* ```./src/``` This directory contains the C++ source files.
  * `flash_atten.h` file contains the implementation of FlashAttention-2 and Flash-D accelerators
* ```./utils/``` This directory contains Python utility scripts.
* ```./setup.sh``` A bash script to fetch all required dependencies.

## Pending Features

* Update design files for to their optimized version.
* Python scripts flows for running PromptBench outputs on C++ designs.
* Provide extention files for 8-bit reduced precision Fast-Float4HLS datatypes.

## Reference

TODO

## Contributors

Currently active: [Kosmas Alexandridis](https://github.com/kosmalex) and [Giorgos Dimitrakopoulos](https://github.com/gdimitrak)



## License

Flash-D is licensed with the MIT License. You are completely free to re-distribute your work derived from Flash-D
