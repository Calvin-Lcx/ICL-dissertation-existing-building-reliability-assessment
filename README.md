# ICL-dissertation-existing-building-reliability-assessment
Imperial College MSc dissertation code: Existing Building Reliability Assessment: A Review and Case Study

The computational framework developed in this study was applied to a case study of a four-storey steel office building located in an urban environment, which has already been in service for thirty years. The aim was to evaluate the time-dependent reliability of the structural system over its 50-year design service life.

Three MATLAB codes were used in this research to perform the reliability analysis of two beam structures and one column structure. Monte Carlo simulation was adopted to estimate system reliability, with a sample size of 2×10⁷ to ensure statistical robustness. Reliability assessment was performed at 5-year intervals, enabling the tracking of structural reliability degradation throughout the service period.

Each code follows a consistent computational framework consisting of the following components: 
(i) Definition of time-dependent basic variables. 
(ii) Specification of potential failure modes.
(iii) Calculation of correlations between failure modes.
(iv) Evaluation of the system reliability index. 
(v) Sensitivity analyses to quantify the influence of corrosion rate, material degradation, and accumulated damage from load history on the overall structural reliability.
(vi) Predict the residual service life of each structural element.

The set of basic variables considered in the codes reflects both material and environmental influences as well as geometric and loading characteristics. 
These variables include: 
1.Corrosion rate; 2.Average metal–environment-specific time exponent; 3.Relative humidity; 4.Temperature; 5.Fatigue damage;
Material properties: 6.Yield strength and 7.Young’s modulus; 
Sectional geometry: 8.Width, 9.Height, 10.Flange thickness, 11.Web thickness, 12.Beam length, and 13.Beam bay width;
Applied loads: 14.Dead load, 15.Live load, and 16.Accumulated damage from load history.

