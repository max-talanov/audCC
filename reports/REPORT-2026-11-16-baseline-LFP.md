# Baseline activity and field potentials of the model without stimulation, and comparison with the literature (16 November 2026)

Corticothalamic sleep model (NEURON, MPI, MareNostrum 5): 8 runs of 120 s on 5,031 cells
(scale 1.65, seeds 0-7, 100 ranks), with no auditory input. The runs are bit-identical to the
first 120 s of the sham arms of F2 (`reports/REPORT-2026-11-10.md`), so this is the undisturbed
network that every stimulation experiment is compared with. The recipe is in `MN5_NEURON.md`
("B120"). File paths are relative to the repository root; raw outputs (`res/`) are not in git.

## 1. Summary

- **Slow oscillation:** 0.70 ± 0.07 UP onsets per second, irregular (CV 0.44). This is inside the
  usual range for the slow oscillation (about 0.75 Hz in human recordings, 0.5-1.5 Hz as a band).
- **Spindles:** at 13.8 Hz they are slightly faster than the human mean (12.3 Hz) but inside the fast
  centroparietal range (13-15 Hz). They start at -47° of the slow oscillation, close to the preferred
  phase in the human data (about -20° to -32°). But there are 26 per minute against 1.6-2.3 per minute in the
  intracranial and scalp recordings, they are shorter (0.51 s against 0.70 s), and every one of them is
  locked to the SO (vector strength 0.99), whereas in human scalp EEG only about 10% of spindles
  closely follow an SO. See section 4 for how much of this is the definition.
- **Field potentials:** the cumulative cortical proxy has 32% of its 0.3-40 Hz power in the SO band
  (0.3-1.5 Hz) and 30% in delta (1.5-4 Hz); the thalamic proxy has 56% in the SO band. L4, TC and RE
  are the most SO-dominated traces (82-89%).
- **UP-state onset:** the UP onset is detected from L5 and L6, so those two layers lead by
  construction. After them, L2/3 interneurons (10 ms), RE (22 ms), L2/3 pyramidal cells (66 ms), TC
  (80 ms) and L4 (85-90 ms) follow.
- **Two things that do not look like sleep data:** L5 fires all the time in the DOWN state (L5 E 2.8
  Hz, L5 I 37 Hz) with a tonic 17.5 Hz rhythm, which Chauvette et al. (2010) found only in slices and
  "never ... in vivo in anesthetized or naturally sleeping cats"; and the interneuron rates are high (L5 I
  44 Hz mean, L4 I 15 Hz, L6 I 17 Hz).
- **Firing rates:** the only cortical reference in the eight papers is the mean rate of recorded
  human units, 3.1 ± 0.3 Hz (Andrillon et al. 2011); none gives rates per cell type, so the model's
  rates per group are not compared numerically.

## 2. Methods

Analysis: `neuron/lfp_baseline.py` (3 unit tests in `neuron/test_lfp_baseline.py`). First 2 s of
each run skipped; values are mean ± SD over the 8 runs.

- **LFP proxy** of a trace = its spikes convolved with a PSP kernel (2 ms rise, 100 ms decay;
  `ctx_analyze._synaptic_lfp`), z-scored, with the deep layers (L5, L6) and RE sign-flipped as in the
  earlier composite LFP. The thalamic traces are proxies, not measured fields, and the sign
  convention is a modelling choice.
- **Traces:** the 10 neuron groups (L2/3, L4, L5, L6 excitatory and inhibitory, TC, RE); the 4
  layers (excitatory and inhibitory spikes pooled); and the cumulative traces, the cortical
  composite (mean of the 4 layers), the thalamic composite (mean of TC and RE) and the network
  (mean of the two).
- **Spectra:** Welch, 8 s segments (0.125 Hz resolution), as fractions of the 0.3-40 Hz power. The
  100 ms kernel attenuates 12 Hz about 8-fold, so the sigma-band share of the LFP proxy is
  underestimated; the sigma peak and share are therefore also taken from the population rate in 5 ms
  bins.
- **UP onsets** from the L5 and L6 excitatory rate (`stage2_report.up_onsets`). UP window = 0-300 ms
  after the onset, DOWN window = -400 to -100 ms.
- **Spindles** = RE volley trains of at least 6 volleys with at most 150 ms between volleys
  (`volley_stats.py`). SO phase at the spindle start from the zero-phase SO band (0.3-1.5 Hz) of the
  cortical E rate, 0° = UP maximum.

## 3. Results

Figures: `figures/b120_traces.png` (12 s of all traces), `figures/b120_spectra.png`,
`figures/b120_rates.png`. Full table: `res/2026-11-16/b120_table.txt`.

![LFP proxies](../figures/b120_traces.png)

### 3.1 Slow oscillation and spindles

| Measure | Value (8 runs) |
|---|---|
| UP onsets per s | 0.70 ± 0.07 |
| CV of UP intervals | 0.44 ± 0.06 |
| Spindles per min | 25.8 ± 2.9 |
| Spindle duration (median) | 0.51 ± 0.04 s |
| Cycles per spindle | 8.0 ± 0.9 |
| Frequency inside a spindle | 13.8 ± 0.1 Hz |
| SO phase at spindle start | -47° (0° = UP maximum), vector strength 0.99, n = 406 |

### 3.2 Cumulative and per-layer field potentials

Fraction of the 0.3-40 Hz power:

| Trace | SO 0.3-1.5 | delta 1.5-4 | theta 4-8 | sigma 9-16 | beta 16-40 |
|---|---|---|---|---|---|
| **cortex (cumulative)** | 0.32 ± 0.02 | 0.30 ± 0.02 | 0.12 | 0.08 | 0.12 ± 0.02 |
| **thalamus (cumulative)** | 0.56 ± 0.03 | 0.20 ± 0.03 | 0.03 | 0.09 | 0.05 |
| **network (cumulative)** | 0.39 ± 0.03 | 0.30 ± 0.02 | 0.11 | 0.05 | 0.08 ± 0.01 |
| L2/3 | 0.51 | 0.22 | 0.10 | 0.09 | 0.04 |
| L4 | 0.86 | 0.07 | 0.01 | 0.02 | 0.02 |
| L5 | 0.54 ± 0.08 | 0.17 | 0.02 | 0.02 | 0.21 ± 0.07 |
| L6 | 0.40 | 0.24 | 0.12 | 0.11 | 0.06 |

![Spectra](../figures/b120_spectra.png)

Each layer and neuron group is in `res/2026-11-16/b120_table.txt`. In short:

- The SO band dominates TC (0.89), L4 I (0.87), L4 E (0.82) and RE (0.83), which follow the UP
  states with a burst and then go quiet.
- L5 has a large share at 16-40 Hz (0.21; L5 I 0.28), which comes from a tonic rhythm with a peak at
  17.5 Hz that runs all the time between UP states (visible in the traces).
- L2/3 and L6 have the most theta and sigma power (0.10-0.12 and 0.09-0.14).
- The thalamic spectrum has the sharpest SO band of all traces, and the cumulative cortical spectrum
  is the flattest, since layers with different timing partly cancel.
- Population-rate spectra (5 ms bins): the highest peak between 5 and 30 Hz is at about 16 Hz in L4,
  TC and RE (the spindle volleys), at 17.5 Hz in L5 (the tonic rhythm) and at 9-10 Hz in L2/3, L6
  I and most L6 E runs. The sigma band (9-16 Hz) holds 22-27% of the rate power in L2/3, L4, L6, TC
  and RE and only 5% in L5.

### 3.3 Firing rates (Hz per cell)

| Group | Mean | UP (0-300 ms) | DOWN (-400 to -100 ms) | UP/DOWN | Onset latency |
|---|---|---|---|---|---|
| L2/3 E | 0.18 ± 0.01 | 0.47 | 0.16 | 2.9 | 66 ± 44 ms |
| L2/3 I | 3.67 ± 0.21 | 8.04 | 3.28 | 2.4 | 10 ms |
| L4 E | 2.39 ± 0.29 | 8.19 | 1.49 | 5.5 | 85 ± 5 ms |
| L4 I | 15.1 ± 1.6 | 43.6 | 8.5 | 5.1 | 90 ms |
| L5 E | 3.53 ± 0.27 | 5.57 | 2.76 | 2.0 | -10 ms |
| L5 I | 44.4 ± 2.9 | 61.6 | 37.3 | 1.7 | 0 ms |
| L6 E | 1.68 ± 0.20 | 4.68 | 1.02 | 4.6 | 0 ms |
| L6 I | 16.6 ± 1.7 | 30.9 | 11.2 | 2.8 | 0 ms |
| TC | 1.71 ± 0.23 | 5.32 | 0.84 | 6.3 | 80 ms |
| RE | 6.33 ± 0.65 | 24.8 | 4.5 | 5.5 | 22 ± 4 ms |

![Rates](../figures/b120_rates.png)

The latency is the time at which the rate crosses half of the rise from the DOWN level to the peak,
in 10 ms bins. The L5 and L6 values are 0 by construction, because the UP onset is defined on them.

## 4. Comparison with the literature

The eight papers in `~/Downloads` were read for the numbers below (page-level values from the full
texts). Human numbers come from intracranial and scalp recordings, and depend on the detector, so
they are ranges, not targets.

| Reference | What it is |
|---|---|
| Andrillon et al. 2011, J Neurosci 31:17821 (`zns17821.pdf`) | intracranial EEG and single units in 13 patients; spindle density, frequency, SO coupling, unit firing |
| Gonzalez et al. 2022, J Neurosci 42:4517 (`zns4517.pdf`) | 20 patients, 357 cortical sEEG channels; spindle density, duration, frequency (Table 2) |
| Staresina et al. 2015, Nat Neurosci 18:1679 (`emss-64810.pdf`) | SO-spindle-ripple nesting in human intracranial recordings |
| Malerba et al. 2022, Sleep (`zsac132.pdf`) | scalp EEG, 34 subjects; spindle types and coupling to SOs |
| Coppieters 't Wallant et al. 2016, Neural Plast (`NP2016-6783812.pdf`) | review of spindle definition and detection |
| Chauvette et al. 2010, Cereb Cortex 20:2660 | origin of active states, cat, multisite recordings |
| Miyawaki and Diba 2016, Curr Biol 26:893 (`nihms-767564.pdf`) | hippocampal firing across NREM (rodent); defines DOWN states |
| Halonen et al. 2021, Front Behav Neurosci 15:712774 (`fnbeh-15-712774.pdf`) | SO-spindle coupling and memory; background only |

| Feature | Model | Literature (source) | Verdict |
|---|---|---|---|
| SO frequency | 0.70 UP/s | about 0.75 Hz (Staresina); slow waves 0.5-4 Hz and "about 1 Hz" in NREM (Miyawaki and Diba) | inside the range |
| Spindle frequency | 13.8 Hz (volley trains) | cortex mean 12.28 Hz NREM (Gonzalez, Table 2); fast centroparietal 13-15 Hz, slow frontal 9-12 Hz (Andrillon); AASM 11-16 Hz, most commonly 12-14 Hz (Coppieters) | in the fast range, about 1.5 Hz above the cortical mean |
| Spindle duration | 0.51 s | cortex 0.70 s NREM, 0.68-0.72 s per region (Gonzalez, Table 2); minimum 0.5 s "seemingly without any firm objective physiological" basis, some authors use 0.3 s (Coppieters) | shorter |
| Spindle density | 25.8 per min | cortex 1.64 per min NREM (frontal 2.12, parietal 2.30; Gonzalez); scalp 1.9 per min, frontal 1.8, parietal 1.3 (Andrillon) | 11-16 times too many |
| Intraspindle frequency | not measured | decreases within a spindle, about -0.8 Hz/s (Andrillon) and more in frontal than parietal spindles (Gonzalez) | not tested |
| Spindles and the SO up-state | start at -47°, R = 0.99, all spindles | spindles occur more often after the transition to the up-state; preferred SO phase of spindle power 328° (about -32°) at Cz and 160° with the up-state at 180° in the hippocampus, i.e. about -20° (Staresina); only about 10% of S2 spindles closely follow an SO (Malerba) | right phase, far too many coupled |
| Spindle activity of units | not measured | firing rates "are not consistently modulated" during spindles, although 19.5% of units are phase-locked (Andrillon) | not tested |
| Where an UP state starts | L5/L6 first (by construction), then L2/3 I, RE, L2/3 E, TC, L4 | deep cells fire first on average by a few ms (-1.1 ± 7.4 ms for deep and +0.5 ± 4.4 ms for middle cells in the sigmoid-fit method, superficial cells positive and close to zero), but activity "may start in any neuron" (Chauvette) | the same direction but a much larger spread in the data; the model result is circular (section 5) |
| Layer 5 in the DOWN state | L5 E 2.8 Hz, L5 I 37 Hz, tonic 17.5 Hz rhythm | L5 cells firing continuously in the silent state was seen only in slices; "never observed in vivo" in anesthetized or naturally sleeping cats (Chauvette) | disagrees |
| DOWN state | population never silent in L5 | defined as more than 50 ms without unit spikes (Miyawaki and Diba, citing refs 8 and 9) | the L5 population never satisfies it |
| Mean firing rate | per group 0.18-44 Hz (table 3.3) | recorded human units 3.1 ± 0.3 Hz on average; phase-locked units 5.9 ± 1.2 Hz (Andrillon); no rates per cell type in the eight papers | not comparable per group |

What the comparison says:

1. **The slow oscillation is right and the spindle phase is right.** The UP rate is in the human range,
   and the spindle start at -47° of the SO is close to the preferred phases found in humans (-32° and
   -20°).
2. **The spindle count is too high.** The density is 11-16 times the human values (1.6-2.3 per min), and
   only part of that is the definition. The model counts RE volley trains of at least 6 volleys at
   intervals up to 150 ms, which is a looser criterion than amplitude detection in a band-passed LFP; I
   have not run the same detector on the model's LFP proxy. A threshold of 6 cycles also puts a floor
   under the duration (6 cycles at 13.8 Hz is 0.43 s), so the 0.51 s is mostly the criterion.
3. **The coupling to the SO is far too deterministic.** In humans spindles are only partly coupled to
   SOs (about 10% of S2 spindles closely follow an SO in the scalp data, Malerba), and the model has
   one fixed phase for all of them (R = 0.99). This suggests that spindles in the model are triggered by
   the UP-state volley and do not appear on their own, as they do in recordings. This is my
   interpretation.
4. **Layer 5 behaves unlike cortex in vivo.** Its tonic firing between UP states is the pattern that
   Chauvette et al. found only in slices. The L5 I mean of 44 Hz is high against the human unit mean of
   3.1 Hz, but that mean mixes cell types, so I treat it as a flag, not a result.

## 5. Caveats

- The onset latency is circular for L5 and L6, because the UP onset is detected on them. The
  conclusion "deep layers first" cannot come from this analysis; an onset detected from a neutral
  signal (for example the sum of all cortical spikes, or the L2/3 E rate) is needed.
- LFP proxies, not fields: a spike-based PSP kernel with a sign convention, no geometry, no real
  extracellular potential. The 100 ms kernel hides the spindle band, so the sigma shares of the
  proxy are lower than in a real LFP.
- Runs are 120 s (about 85 UP states each) and 8 seeds; the SD between runs is shown, not a
  bootstrap interval. The latency of L2/3 E varies a lot between runs (66 ± 44 ms).
- The literature numbers were read in the full texts of the eight papers (in the Downloads folder). Malerba et al. 2022 was only partly readable (the PDF is mostly images; the abstract was
  read, the results were not). Spindle recurrence every 4 s, which I mentioned before, is not in any of
  the eight papers, so it is dropped. Rates per cell type are not in them either.
- The human values are from patients with epilepsy (intracranial) and scalp EEG; Chauvette is in cats;
  Miyawaki and Diba is in rodents. The comparison is order-of-magnitude, not a fit.

## 6. Next steps

1. Redo the onset latencies with an onset detected from a signal that does not use L5 or L6.
2. Detect spindles on the cortical LFP proxy with an amplitude criterion (sigma-band envelope), to
   compare density and duration with the human numbers on equal terms, and measure the spread of the
   spindle phase on that detector.
3. Find NREM firing rates per cell type (pyramidal, interneuron; UP and DOWN) in papers with unit
   recordings, and compare with the table in section 3.3. Look at the L5 tonic rhythm and the L5
   DOWN-state activity in the model (L5 recurrence and IB cells), since Chauvette et al. found that
   pattern only in slices.
3b. Measure spindle coupling per spindle (fraction of spindles that follow an UP within a window, as in
   the scalp analysis) to compare with the 10% value.
4. Repeat after any change in the thalamic excitability (K_Na, IC-to-TC weight), since the same
   network is used for all stimulation experiments.

## 7. Where things are

| What | Where |
|---|---|
| Analysis and tests | `neuron/lfp_baseline.py`, `neuron/test_lfp_baseline.py` |
| Recipe (B120) | `MN5_NEURON.md` |
| Table | `res/2026-11-16/b120_table.txt` (local) |
| Figures | `figures/b120_traces.png`, `figures/b120_spectra.png`, `figures/b120_rates.png` |
| Earlier reports | `reports/REPORT-2026-11-10.md`, `reports/REPORT-2026-10-23-stimulation-comparison.md` |
