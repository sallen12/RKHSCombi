# RKHSCombi

This repository contains R code to reproduce the results presented in the paper  

> Allen, S., Ginsbourger, D. and Ziegel, J. (2024).
> Efficient pooling of predictions via kernel embeddings.
> _Transactions on Machine Learning Research_.
> [openreview.net/forum?id=dji9MfONNP](https://openreview.net/forum?id=dji9MfONNP)

The paper concerns methods to efficiently combine probabilistic forecasts via kernel mean embeddings.

The results can be reproduced and saved by sourcing the `main.R` file. The saved results are stored in the `Results` folder. Functions to perform the RKHS-based forecast combination methods can be found in the `utility_funcs.R` file. 

We are extremely grateful to MeteoSwiss (The Swiss Federal Office of Meteorology and Climatology MeteoSwiss) for providing the wind speed data used in this study, and for allowing it to be shared publicly. The data is available in the `Data` folder. A separate `RDS` file is provided for the forecasts and observations at each lead time. Files containing the station info, such as longitude, latitude, and height, are also available. 
