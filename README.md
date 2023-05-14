# Structure 
This repository contains the code for the paper "On Perfect Linear Approximations and Differentials over Two-Round SPNs".
The code is structured into three directories.
[ciphers](./ciphers) contains basic implementations (sboxes and linear layers) of the examined ciphers.
[linear](./linear) contains our code for the linear and [differential](./differential) the code for the differential case.
[main.sage](main.sage) runs all our experiments.

# Dependencies
Our code is based on [sage](https://www.sagemath.org/).
We also make use of a sage [module for linear layers](./ciphers/linearlayer.py) which is not part of sage.
Its original version can be found [here](https://git.sagemath.org/sage.git/tree/src/sage/crypto/linearlayer.py?h=u/asante/linear_layer_module&id=9155aa93f62e2c1bab51a095137b62367c1e7fbb).
Further, we use [tabulate](https://pypi.org/project/tabulate/) to generate nice
tables and [tqdm](https://pypi.org/project/tqdm/) as a progress bar.
To install them, run `pip install tabulate` and `pip install tqdm` inside sage.

# Executing the Code
Run `sage main.sage` to run the experiments.
We used `SageMath version 9.6`.
Earlier versions might not work!
Please note: in some cases sage is not good with Gröbner basis computations. Hence, sometimes we generate magma files instead.
