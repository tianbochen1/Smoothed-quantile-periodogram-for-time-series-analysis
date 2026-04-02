# Smoothed-quantile-periodogram-for-time-series-analysis

	
The non-smoothness and non-differentiability of the check loss function in quantile regression (QR) hinder the use of gradient-based optimization methods and consequently reduce computational efficiency. This issue also affects the quantile periodogram (QP), a spectral estimator based on trigonometric QR. In this paper, we smooth the check loss function using the generalized multiquadric (GMQ) function and accordingly propose the smoothed quantile periodogram (SQP). We investigate the asymptotic relationship between the SQP and its associated smoothed quantile spectrum (SQS), and further develop an SQP-based Fisher--$g$-type test for periodicity detection. Simulation studies show that the SQP substantially improves computational efficiency while retaining the ability of the original QP to detect hidden periodicities with only a modest loss of robustness. Moreover, the SQP is smooth across quantile levels. In practice, spectral estimation can therefore be obtained by smoothing only over frequency, without an additional smoothing step along the quantile levels.  Finally, we apply the proposed SQP to the analysis of bearing data, where it successfully detects a hidden frequency at 160 Hz and exhibits greater temporal stability than the PG.


'fn.R': code 
'IR007_0.mat': bearing data
'cwru.R': bearing data reproduction
