
Data description
----------------

Study population: adults with type 1 diabetes, who are in the Steno Diabetes Center system. It is recommended for patients to come for regular follow-up visits every 4 months. These data represent information available from the past year at the time that the predictions are to be made, which includes some measurements at 2 time points: pre and post time 0. At time 0, a decision was made to either start statins or not. The variable pid is a unique patient identifier.

Outcome: Cardiovascular disease within 5 years of time 0. The is coded as 0 for no, 1 for yes, in the variable cvd_5year

Predictors, all measured at or before time 0 but no more than 1 year before time 0: 

- age: Age in years at the time 0
- sex_male: Binary indicator of male sex (0 = female, 1 = male)
- diabetes_duration: Time in years since diagnosis of type 1 diabetes
- steno_prs: A continuous polygenic risk score for CVD in in type 1 diabetes.
- smoking: Binary indicator of current smoking status (0 = nonsmoker, 1 = smoker)
- motion: Binary indicator of regular physical activity (0 = sedentary, 1 = some moderate physical activity)
- HBA1C_pre_trt, HBA1C_post_trt: Glycated hemoglobin, also called hemoglobin A1C value at two time points in mmol/mol. Smaller values indicate better glycemic control. 
- urine_albumin_pre_trt, urine_albumin_post_trt: The urine albumin to creatinine ratio, in mg/g. A larger relative concentration of albumin in the urine is an indicator of poor kidney function
- LDL_pre_trt, LDL_post_trt: Low-density lipoprotein measurement in serum, in mmol/L. 
- SBP_pre_trt, SBP_post_trt: Systolic blood pressure in mmHg
- eGFR_pre_trt, eGFR_post_trt: Estimated glomerular filtration rate, in mL/min/1.73m^2. An estimate of kidney function, larger values = better function
- statin: indicator of initiating statin use at time 0 (0 = not on statins, 1 = on statins)


