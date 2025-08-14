Objectives for the practicals Day 1

1. Read about the data, read it into R, and examine the data and benchmark models
2. Make your own prediction model using the training dataset
3. Provide an estimate of some metric to quantify how well your model performs
4. Predict the probability of the outcome in the validation set


Data description
----------------

Study population: adults with type 1 diabetes, who are in the Steno Diabetes Center system. It is recommended for patients to come for regular follow-up visits every 4 months. These data represent information available from the past year at the time that the predictions are to be made, which includes some measurements at 2 time points. The variable pid is a unique patient identifier.

Outcome: Cardiovascular disease within 5 years of time 0. The is coded as 0 for no, 1 for yes, in the variable cvd_5year

Predictors, all measured at or before time 0 but no more than 1 year before time 0: 

- age: Age in years at the time 0
- sex_male: Binary indicator of male sex (0 = female, 1 = male)
- diabetes_duration: Time in years since diagnosis of type 1 diabetes
- smoking: Binary indicator of current smoking status (0 = nonsmoker, 1 = smoker)
- motion: Binary indicator of regular physical activity (0 = sedentary, 1 = some moderate physical activity)
- HBA1C_time1, HBA1C_time2: Glycated hemoglobin, also called hemoglobin A1C value at two time points in mmol/mol. Smaller values indicate better glycemic control. 
- urine_albumin_time1, urine_albumin_time2: The urine albumin to creatinine ratio, in mg/g. A larger relative concentration of albumin in the urine is an indicator of poor kidney function
- LDL_time1, LDL_time2: Low-density lipoprotein measurement in serum, in mmol/L. 
- SBP_time1, SBP_time2: Systolic blood pressure in mmHg
- eGFR_time1, eGFR_time2: Estimated glomerular filtration rate, in mL/min/1.73m^2. An estimate of kidney function, larger values = better function
- statin: indicator of statin use between times 1 and 2 (0 = not on statins, 1 = on statins)


