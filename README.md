# Matched Field Processsing for complex Earth structure

Repository to accompany the paper entitled "Matched Field Processsing for complex Earth structure" by Schippkus & Hadziioannou, published as pre-print on EarthArXiv (---) and submitted to --- for peer review.

Here, we provide all scripts and data necessary to reproduce all of our results.

`/tutorial_notebooks` contains a Python Jupyter notebook that introduces standard Matched Field Processing

`/settings_files` contains yaml settings files used in the Matched Field Processing code developed for this paper (see [seismology-hamburg/matched_field_processing](https://github.com/seismology-hamburg/matched_field_processing)). These are the exact input files used to genereate the synthetic tests and real-data results. Note that for all figures, multiple `.yml` files are provided. These correspond to the subplots.

`/data_download` contains Python scripts to download the seismic waveforms and station metadata to reproduce Figures 5 & 6. For the Chino Hills earthquake, this includes data from 54 seismic stations of the CI network (estimated file size 40 MB). For the continous seismic data in February 2019, this includes data from 342 of the BE, BW, CH, CZ, DK, EI, FR, G, GB, GE, GR, GU, II, IM, IU, IV, LX, MN, NL, NO, OE, OX, PL, PM, RD, SX, TH, and WM networks (estimated file size 2.0 GB). For citations of these networks refer to our paper or [fdsn.org/networks](https://fdsn.org/networks).

`/figure_scripts` contains Python scripts to produce the figures in our paper from the output files of [seismology-hamburg/matched_field_processing](https://github.com/seismology-hamburg/matched_field_processing) using the settings files above. We also directly provide the output files.
