# MIARMA
Sources files for the gap-filling algorithm based on ARMA models.
Reference for this algorithm: Pascual-Granado, J., Garrido, R., Suárez, J. C. A&A 2015, 575, A78

We used a forward-backward predictor based on autoregressive moving-average modelling (ARMA) in the time domain. The algorithm is particularly suitable for replacing invalid data such as those present in the light curves of CoRoT satellite due to the passing through the South Atlantic Anomaly, and eventually for the data gathered by the NASA planet hunter Kepler.

Files in folder \Plus are external tools developed by the same author for other purposes which are also necessary to run MIARMA.

## Version numbers
Each version has (since 1.5.11.0) four numbers x.x.x.x, e.g. 1.0.0.0
- The fourth number refers to just minor corrections in MIARMA.m or individual subroutines not involving changes in other ones, e.g. 1.0.0.1 
- The third number refers to corrections involving several subroutines, e.g. 1.0.1.0
- The second number refers to major changes in the algorithms, e.g. 1.1.0.0
- The first number refers to a complete reface of the program or significant changes in the interaction, e.g. 2.0.0.0

Code modifications aimed to improve the readibility of the code (e.g. comments, tabs, etc.) don't change the version number.

The numbering used is different for the different subroutines called.
