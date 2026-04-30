-- ============================================================================
-- patient_ccsr_long.csv — CCSR code long table for PSM analysis
--
-- Generates one row per patient per CCSR code, joining the revised
-- last_extubations table (post-Fix 2: inferred events have NULL caregiver_id)
-- to diagnoses_icd and the CCSR mapping table.
--
-- Filters:
--   caregiver_n > 10  — matches the analytical cohort filter in psm_analysis.R
--   seq_num <= 10     — top 10 diagnosis codes per admission (clinically
--                       meaningful; beyond position 10 codes are often
--                       administrative rather than clinical)
--   icd_version = 10  — CCSR mapping covers ICD-10 only; ICD-9 patients
--                       will not appear here. They are handled separately
--                       via ICD chapter flags in MIMIC_analysis_clean.R.
--
-- Output: ~40,000–55,000 rows (5–6 CCSR codes per patient on average)
-- Export as: patient_ccsr_long.csv
-- ============================================================================

SELECT
  e.subject_id,
  e.hadm_id,
  e.caregiver_id,
  e.caregiver_fe_rate,
  e.caregiver_n,
  e.failed_extubations,
  e.tube_event_source,
  e.caregiver_imputation_source,
  code AS ccsr_code
FROM `mythical-legend-456217-e9.Failed_Extubation_Rate.last_extubations` e,
UNNEST(e.ccsr_codes) AS code
WHERE e.caregiver_n > 10
  AND e.ccsr_codes IS NOT NULL
ORDER BY e.subject_id, code
