# Calibration check for malware detectors

A web page that measures how well a malware detector's probabilities match observed
malware frequencies over time and applies the decision procedure of the article
"Temporal Calibration Failure and Delayed Recalibration in Malware Detection" by
Ali Al Bataineh, Munther Gdeisat, Ali Aburomman, Suleiman Almasri, Abdulsalam
Alarabeyyat and Akram Alkouz.

Open the page: https://aalbata.github.io/calibration-check/

All calculations run in the browser. The file you select is read on your computer and
is not uploaded, and the page loads no external scripts, fonts or other resources.

## Use

1. Open the page.
2. Select a CSV file, or load one of the two examples.
3. Check the role of each period, set the label delay, and run the check.

### Input file

One row per sample, with a header row. Column names are not case sensitive.

| Column | Content |
|---|---|
| `period` | Time period of the sample, for example `2024-03`. Periods are ordered alphabetically, so use a sortable format. |
| `label` | 1 for malware, 0 for benign. |
| `score` or `prob` | A log-odds (logit-scale) score, or a probability between 0 and 1. |
| `role` (optional) | `calibration`, `development` or `deployment`. |

The calibration window is the labeled data on which the probability map is fitted
before deployment, development periods follow it, and deployment periods are the later
periods whose labels have arrived. Roles can be assigned or changed on the page.

The page checks the file before it calculates.

- A row is read only when its period is not empty, its label is 0 or 1, and its score
  or probability is a written number. An empty cell is never read as zero. Rows that
  fail are left out and reported with their line in the file.
- The roles must follow the order calibration, development, deployment in time. The
  page does not run when a period is assigned to an earlier role than a period before
  it, and it names the periods concerned.
- When a development period or the first deployment period holds samples of one class
  only, the page reports insufficient data and gives no recommendation.
- A period with too few samples per bin for its calibration error to be read is marked
  as sparse and listed next to the recommendation. When no bin of a period can be used,
  its error is shown as n/a in every output, and averages over periods are shown only
  when every period entered them.
- With the map set to none, the values are used as given, and the page reports no
  refits and no labels used for them.

### Output

- The recommended action, with each question of the decision procedure of the article
  (its Fig. 6) and the quantities behind the answer.
- Per period: sample count, malware share, AUC, the debiased squared calibration error
  (DSCE) and its square root (RMSCE), and the calibration intercept and slope.
- Slot-averaged RMSCE and DSCE on the deployment periods when the map is never refitted
  or refitted every one, two or three periods, with the labels each schedule uses.
- A CSV export of the per-period values.

## Scope

The page applies the procedure with point estimates, with the first deployment period
standing for the first period whose labels have arrived. The procedure was derived from
the BODMAS and LAMDA corpora and has not been validated prospectively. The page
implements three of the six probability maps of the article and fixed refit schedules.
Isotonic regression, beta calibration, scaling-binning, drift-triggered refits and
bootstrap intervals are available in the Python package of the main repository only.

## Source and documentation

`index.html`, `calibration.js` and `examples.js` are copies of the `app/` folder of
https://github.com/aalbata/temporal-calibration-malware, the repository of the article.
That repository holds the analysis code and saved outputs, the full description of this
page (`app/README.md`), the script that writes the two examples, and the test that
compares the calculations of the page with the Python implementation. The button
"Run the self-check" on the page repeats part of that comparison in the browser.

To run the page without a network connection, download this repository, run
`python -m http.server 8000` in its folder, and open `http://localhost:8000`.

## License

The MIT License (see `LICENSE`) covers the code. The two examples contain model scores
and labels derived from BODMAS and LAMDA. They are provided to demonstrate the page,
are not covered by the MIT License, and remain subject to the terms of those datasets.

## Citation

If you use the page, please cite the article named above and the main repository.
