# ForTran2R_transaltion_SBILysis
<span style="color: green"> Fortran to R translation | Simulation Based Inference (SBI) with Markov Chain Monte Carlo (MCMC) - A. A. Ammerman | OG Fortran code : Clotting Simulation by Dr. Brittany Bannish of the University of Central Oklahoma</span >

---
This is a translation from Fortran to R of the original source code for Dr.B Bannish's Micro-scale Lysis Simulation (which can be found here: [Simulation Repository](https://github.com/UCO-OpResearch/lysis); likewise, see the Original Fortran Code (11/01/2026 - MM/DD/YYYY)[here](https://github.com/UCO-OpResearch/lysis/blob/main/src/fortran/micro_rates.f90). Please do mind that some of my descriptions on context and SBI may not be fully adequate at this time.


The purpose of translating this was firstly to become well acquainted with the logic behind simulation, secondly to make this simulation more accessible for those more versed in R or other similar, modern languages; and thirdly (the not-as-important reason), plainly just because R is something I can debug more naturally in the event of some error. So, I stubbornly spent an extra 6 or 7 hours learning FT code and deciphering it into R (regex being a major help here).

The next step from this working model is to successfully implement MCMC Sampling for generating a healthy sample for Simulation Based Inference.

## Status

Currently the code is in rough shape with scattered, unpolished comments. Sections may look shrunk since R has certain functions and built in features that supersede certain set-ups necessary for the original Fortran Code (e.g. cumsum, sample, runif). However, it is necessary for me to loop back around to remove other hard-coding that I may have persisted from the original.

Another good practice down the line would be for me to make a logic flow diagram to nail down the *intention* of the code. This is not yet in the works.

### Current Simulation speed

Notably, I am doing this on a laptop that is allegedly "cyber waste." Although the conversion to Linux Mint freed up computational power of the machine, it is fair to note that simulations will likely run faster on literally any device with more more ram and cores than mine (i have 16 gb and 2 cores). To test later: I'll attempt this on my much more powerful desktop and see how fast it can pass 1000 simulations.

At this point in time (one a 2012ish? laptop) it took 8 minutes for 1000 runs.

|          Computer     |           Spec Type       |   Spec Detail  |
|-----------------------|---------------------------|----------------|
|2012 Cyberwaste Laptop | Time per 1000 Simulations |  approx. 8 min |


## References

As is necessary, if you happened across this by chance, here are the references that either use the model and/or describe it:

Bannish, Brittany E., James P. Keener, and Aaron L. Fogelson. "Modelling fibrinolysis: a 3D stochastic multiscale model." Mathematical medicine and biology: a journal of the IMA 31.1 (2014): 17-44. https://doi.org/10.1093/imammb/dqs029

Papers that used this code (see the matching releases [here](https://github.com/UCO-OpResearch/lysis/tree/main#releases-and-publications)):

Risman, R. A., Paynter, B., Percoco, V., Shroff, M., Bannish, B. E., & Tutwiler, V. (2024). Internal fibrinolysis of fibrin clots is driven by pore expansion. Scientific Reports, 14(1), 2623. https://doi.org/10.1038/s41598-024-52844-4 (dataset: https://doi.org/10.5281/zenodo.8115180)

Bannish, B. E., Paynter, B., Risman, R. A., Shroff, M., & Tutwiler, V. (2024). The effect of plasmin-mediated degradation on fibrinolysis and tissue plasminogen activator diffusion. Biophysical Journal, 123(5), 610-621. https://doi.org/10.1016/j.bpj.2024.02.002

Risman, R. A., Percoco, V., Paynter, B., Bannish, B. E., & Tutwiler, V. (2025). Protofibril packing density of individual fibers alters fibrinolysis. Research and Practice in Thrombosis and Haemostasis, 9(2), 102708. https://doi.org/10.1016/j.rpth.2025.102708 (dataset: https://doi.org/10.5281/zenodo.15151618)
