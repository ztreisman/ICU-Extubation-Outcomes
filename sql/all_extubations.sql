-- ============================================================================
-- All Extubation Events — for trajectory analysis
--
-- Unlike ICU_Last_Extubations.sql which selects one row per patient (the
-- last extubation), this query returns ALL extubation events per patient
-- with the per-event failed_extubation_flag as the outcome.
--
-- This is the correct design for trajectory analysis: does physiology
-- in the 72h before THIS extubation predict whether THIS extubation fails?
--
-- failed_extubation_flag = 1 if this specific extubation was followed by
--   reintubation within 72h, where both events are explicit procedure events.
--
-- Restricted to:
--   - Explicit events only (source = 'explicit')
--   - caregiver_n > 10 (FE rate estimated from > 10 extubations)
--   - vent_hours > 1 and < 1000
--
-- Output: one row per extubation event
-- Export as: all_extubations.csv
--
-- For trajectory analysis, join to RSBI and P/F measurements in R via
-- bigrquery using stay_id and event_time.
-- ============================================================================

WITH tube_events AS (
  SELECT *
  FROM (
    -- Explicit events only
    SELECT
      p.subject_id,
      p.hadm_id,
      p.stay_id,
      p.caregiver_id,
      p.starttime AS event_time,
      CASE
        WHEN p.itemid = 224385 THEN 'intubation'
        WHEN p.itemid = 227194 THEN 'extubation'
        WHEN p.itemid = 225468 THEN 'unplanned_extubation_patient'
        WHEN p.itemid = 225477 THEN 'unplanned_extubation_nonpatient'
      END AS tube_event,
      'explicit' AS source
    FROM `physionet-data.mimiciv_3_1_icu.procedureevents` p
    WHERE p.itemid IN (224385, 227194, 225468, 225477)
  )

),

ranked_events AS (
  SELECT
    *,
    LAG(tube_event)   OVER (PARTITION BY subject_id ORDER BY event_time) AS prev_tube_event,
    LAG(event_time)   OVER (PARTITION BY subject_id ORDER BY event_time) AS prev_event_time,
    LEAD(event_time)  OVER (PARTITION BY subject_id ORDER BY event_time) AS next_event_time,
    LEAD(tube_event)  OVER (PARTITION BY subject_id ORDER BY event_time) AS next_tube_event,
    LEAD(source)      OVER (PARTITION BY subject_id ORDER BY event_time) AS next_event_source
  FROM tube_events
),

-- Ventilation duration for each extubation event
ventilation AS (
  SELECT
    stay_id,
    starttime,
    endtime,
    ventilation_status
  FROM `mythical-legend-456217-e9.Failed_Extubation_Rate.ventilation`
  WHERE ventilation_status = 'InvasiveVent'
),

extubations AS (
  SELECT
    re.subject_id,
    re.hadm_id,
    re.stay_id,
    re.caregiver_id,
    re.event_time,
    re.source,
    -- vent_hours: from previous intubation event
    CASE
      WHEN re.prev_tube_event = 'intubation'
           AND TIMESTAMP_DIFF(re.event_time, re.prev_event_time, SECOND) > 0
        THEN TIMESTAMP_DIFF(re.event_time, re.prev_event_time, SECOND) / 3600.0
      WHEN v.starttime IS NOT NULL
        THEN TIMESTAMP_DIFF(v.endtime, v.starttime, SECOND) / 3600.0
      ELSE NULL
    END AS vent_hours,
    -- per-event failed extubation flag
    -- requires both extubation and reintubation to be explicit events
    IF(
      re.next_tube_event = 'intubation'
      AND TIMESTAMP_DIFF(re.next_event_time, re.event_time, HOUR) <= 72
      AND re.next_event_source = 'explicit',
      1, 0
    ) AS failed_extubation_flag
  FROM ranked_events re
  LEFT JOIN ventilation v
    ON re.stay_id = v.stay_id
    AND v.endtime = re.event_time
  WHERE re.tube_event = 'extubation'
    AND re.source = 'explicit'
),

-- Caregiver stats: leave-one-out FE rate
-- For each event, compute caregiver FE rate from all OTHER events
-- to avoid circularity (this event's outcome contributing to its own exposure)
caregiver_all AS (
  SELECT
    caregiver_id,
    COUNT(*)                                          AS caregiver_n_total,
    SUM(CAST(failed_extubation_flag AS FLOAT64))      AS caregiver_total_failures
  FROM extubations
  WHERE caregiver_id IS NOT NULL
    AND vent_hours > 1
    AND vent_hours < 1000
  GROUP BY caregiver_id
),

patient_info AS (
  SELECT subject_id, anchor_age, gender, dod
  FROM `physionet-data.mimiciv_3_1_hosp.patients`
),

mortality_info AS (
  SELECT subject_id, MAX(hospital_expire_flag) AS hospital_expire_flag
  FROM `physionet-data.mimiciv_3_1_hosp.admissions`
  GROUP BY subject_id
),

first_day_sofa AS (
  SELECT subject_id, MAX(sofa) AS sofa
  FROM `physionet-data.mimiciv_3_1_derived.first_day_sofa`
  GROUP BY subject_id
),

charlson AS (
  SELECT subject_id, MAX(charlson_comorbidity_index) AS charlson
  FROM `physionet-data.mimiciv_3_1_derived.charlson`
  GROUP BY subject_id
),

norepinephrine_equivalent AS (
  SELECT
    i.subject_id,
    MAX(n.norepinephrine_equivalent_dose) AS max_norepinephrine_equivalent_dose
  FROM `physionet-data.mimiciv_3_1_derived.norepinephrine_equivalent_dose` n
  JOIN extubations i ON n.stay_id = i.stay_id
  GROUP BY i.subject_id
)

SELECT
  e.subject_id,
  e.hadm_id,
  e.stay_id,
  icu.first_careunit,
  p.anchor_age,
  p.gender,
  sofa.sofa,
  c.charlson,
  COALESCE(n.max_norepinephrine_equivalent_dose, 0) AS norepinephrine,
  e.vent_hours,
  e.event_time,
  e.source                                          AS tube_event_source,
  e.caregiver_id,
  -- Leave-one-out caregiver FE rate
  CASE
    WHEN ca.caregiver_n_total > 1 THEN
      (ca.caregiver_total_failures - CAST(e.failed_extubation_flag AS FLOAT64))
      / (ca.caregiver_n_total - 1)
    ELSE NULL
  END                                               AS caregiver_fe_rate_loo,
  ca.caregiver_n_total                              AS caregiver_n,
  e.failed_extubation_flag,
  m.hospital_expire_flag,
  p.dod
FROM extubations e
LEFT JOIN patient_info p          ON e.subject_id = p.subject_id
LEFT JOIN mortality_info m        ON e.subject_id = m.subject_id
LEFT JOIN first_day_sofa sofa     ON e.subject_id = sofa.subject_id
LEFT JOIN charlson c              ON e.subject_id = c.subject_id
LEFT JOIN norepinephrine_equivalent n ON e.subject_id = n.subject_id
LEFT JOIN caregiver_all ca        ON e.caregiver_id = ca.caregiver_id
LEFT JOIN `physionet-data.mimiciv_3_1_icu.icustays` icu ON e.stay_id = icu.stay_id
WHERE e.vent_hours > 1
  AND e.vent_hours < 1000
  AND e.caregiver_id IS NOT NULL
  AND ca.caregiver_n_total > 10
ORDER BY e.subject_id, e.event_time
