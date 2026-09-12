# Optional validation inputs

Two small tables used only for the measurement-validation results are not part of the Mendeley Data deposit and are not distributed with this repository. `run_all.py` detects whether they are present and skips the two steps that need them, so everything else is reproduced without them.

| File | Contents | Reproduces |
|---|---|---|
| `sound_event_validation_clips.csv` | one row per 10 s validation clip (300): `clip_id`, `source_audio_file`, `day_night`, the classifier's primary and secondary event codes and the three annotators' primary and secondary codes (`model_primary_event_code`, `A_primary_event`, `B_primary_event`, `C_primary_event`, `model_secondary_event_code`, `A_secondary_event`, `B_secondary_event`, `C_secondary_event`); blank cells are unlabelled | Supplementary Table S2 (`33_sound_event_validation.py`, then `si_tables.R`) |
| `collocation_pairs.csv` | one row per collocated 10-min session (14): `session`, `time` (HH:MM), `reference_db` (Class 1 reference sound level meter LAeq) and `recorder_db` (field recorder LAeq), both in dB(A) | the recorder-versus-reference bias, RMSE and regression figures quoted in Methods (`32_collocation_check.R`) |

Event codes follow the twelve-category scheme described in the manuscript's Methods; `33_sound_event_validation.py` maps them to the four source groups (traffic, human, natural, mechanical).

Requests for these tables go to the corresponding author.
