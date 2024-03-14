# planck-camspec-lkls
CamSpec likelihoods using Planck data

This repository contains code and data for the CamSpec PR4_12.6cl and CamSpec_PR4_12.7cl likelihoods, described in [Rosenberg, Gratton, and Efstathiou (2022), MNRAS, 517, 4620](https://academic.oup.com/mnras/article/517/3/4620/6717656) ([arxiv:2205.10869](https://arxiv.org/abs/2205.10869)) and [Efstathiou, Rosenberg, and Poulin (2024)](https://arxiv.org/abs/2311.00524). These are based on and very similar to the PR3 likelihoods described in detail in [Efstathiou and Gratton (2021), OJA, 4, 8](https://astro.theoj.org/article/27518-a-detailed-description-of-the-camspec-likelihood-pipeline-and-a-reanalysis-of-the-planck-high-frequency-maps).

These are high-ell likelihood using dust-cleaned TT, TE, and EE power spectra computed from the Planck 2020 (PR4 / NPIPE) maps.
Produced by George Efstathiou, Steven Gratton and Erik Rosenberg.
Python likelihood code by Antony Lewis based on Fortran code by Efstathiou, Gratton, and Lewis. Ported to montepython by Erik Rosenberg.

Code for [Cobaya](https://github.com/CobayaSampler/cobaya) is provided with the main Cobaya distribution, and allows for automatic download of the data.

The appropriate code and data for CosmoMC are available for CamSpec PR3_12.5cl [here](https://people.ast.cam.ac.uk/~stg20/camspec/index.html), and for the newer data on request.
