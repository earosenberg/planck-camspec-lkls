This is the montepython likelihood for the CamSpec PR4_12.6cl and CamSpec_PR4_12.7cl likelihoods, described in [Rosenberg, Gratton, and Efstathiou (2022), MNRAS, 517, 4620](https://academic.oup.com/mnras/article/517/3/4620/6717656) ([arxiv:2205.10869](https://arxiv.org/abs/2205.10869)) and [Efstathiou, Rosenberg, and Poulin (2024)](https://arxiv.org/abs/2311.00524).
This is a high-ell likelihood using dust-cleaned TT, TE, and EE power spectra computed from the Planck 2020 (PR4 / NPIPE) maps.
Produced by George Efstathiou, Steven Gratton and Erik Rosenberg.
Python likelihood code by Antony Lewis based on Fortran code by Efstathiou, Gratton, and Lewis. Ported to montepython by Erik Rosenberg.

To use, copy the contents of ./data/ to [montepython_dir]/data and ./code to [montepython_dir]/montepython/likelihoods
Note that all likelihoods require Planck_CamSpec_NPIPE_TTTEEE/Planck_CamSpec_NPIPE.py
