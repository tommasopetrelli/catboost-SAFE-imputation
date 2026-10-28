const fs = require('fs');
const { Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell,
  WidthType, ShadingType, AlignmentType, LevelFormat, BorderStyle, Footer, PageNumber } = require('docx');

const FONT = 'Calibri';
const W = 9026; // A4 text width with 1" margins (11906 - 2*1440)
const border = { style: BorderStyle.SINGLE, size: 4, color: 'B7B7B7' };
const borders = { top: border, bottom: border, left: border, right: border };

const r = (t, o = {}) => new TextRun({ text: t, font: FONT, ...o });
const p = (runs, o = {}) => new Paragraph({ children: (typeof runs === 'string' ? [r(runs)] : runs), spacing: { after: 120, line: 276 }, ...o });
const h1 = t => new Paragraph({ heading: HeadingLevel.HEADING_1, children: [r(t)], spacing: { before: 280, after: 140 } });
const h2 = t => new Paragraph({ heading: HeadingLevel.HEADING_2, children: [r(t)], spacing: { before: 220, after: 110 } });
const h3 = t => new Paragraph({ heading: HeadingLevel.HEADING_3, children: [r(t)], spacing: { before: 160, after: 80 } });
const b = (runs, lvl = 0) => new Paragraph({ numbering: { reference: 'bul', level: lvl }, children: (typeof runs === 'string' ? [r(runs)] : runs), spacing: { after: 80, line: 276 } });
const code = t => new Paragraph({ children: [new TextRun({ text: t, font: 'Consolas', size: 18 })], spacing: { after: 0 }, shading: { type: ShadingType.CLEAR, fill: 'F3F3F3' }, indent: { left: 284 } });
const eq = t => new Paragraph({ alignment: AlignmentType.CENTER, children: [r(t, { italics: true })], spacing: { before: 60, after: 120 } });

function table(widths, rows, header = true) {
  return new Table({
    width: { size: W, type: WidthType.DXA }, columnWidths: widths,
    rows: rows.map((row, i) => new TableRow({
      tableHeader: header && i === 0,
      children: row.map((cell, j) => new TableCell({
        borders, width: { size: widths[j], type: WidthType.DXA },
        margins: { top: 60, bottom: 60, left: 100, right: 100 },
        shading: header && i === 0 ? { type: ShadingType.CLEAR, fill: 'DCE6F0' } : undefined,
        children: [new Paragraph({ children: [r(String(cell), { bold: header && i === 0, size: 19 })] })]
      }))
    }))
  });
}
const gap = () => new Paragraph({ children: [], spacing: { after: 60 } });

const children = [
  new Paragraph({ children: [r('CatBoost imputation of SAFE Q8A', { size: 36, bold: true })], spacing: { after: 60 } }),
  new Paragraph({ children: [r('Changes to the pipeline relative to the MSc thesis', { size: 28 })], spacing: { after: 60 } }),
  new Paragraph({ children: [r('Reference: T. Petrelli (2026), “Measuring Financing Gaps with Machine Learning — Survey Evidence from SAFE in the European Union”, Bocconi MSc thesis. Repository: tommasopetrelli/catboost-SAFE-imputation.', { size: 19, color: '555555' })], spacing: { after: 60 } }),
  new Paragraph({ children: [r('Change log version 1 — 25 October 2026', { size: 19, color: '555555' })], border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: '2E5C8A', space: 6 } }, spacing: { after: 200 } }),

  h1('1. Purpose of this document'),
  p('This document records every modification made to the imputation code after the thesis was submitted, and explains for each one what the thesis did, why it was changed, what the code now does, and how the thesis results are affected. It is updated each time the pipeline is modified, so that the paper’s methodology section and any reply to referees can be traced back to specific code changes.'),
  p('Each change has an identifier (C1, C2, …), the script and block it touches, and the thesis section it relates to. The thesis setting can always be reproduced: where a change alters estimates, the original behaviour is kept behind a switch.'),

  h1('2. Summary of changes'),
  table([700, 2700, 1500, 2200, 1926], [
    ['ID', 'Change', 'Thesis ref.', 'Script / block', 'Effect on results'],
    ['C1', 'Class weights: optional, and corrected before evaluation and imputation draws', '§3.2.2, §4.1', 'Model I — Blocks 0, 1, 5, 6, 7, 8, 9', 'Changes imputed values and all holdout metrics'],
    ['C2', 'TVD and mean shift computed on expected shares (mean of predicted probabilities) instead of argmax labels', '§4.2.1, §4.3.1, Table 9, §4.5', 'Model I — Blocks 1, 8, 11', 'Changes Table 9 distributional rows; argmax version still reported'],
    ['C3', 'Correct grouping variables for “by sector” and “by size” diagnostics', 'Table 9', 'Model I — Block 8', 'Relabels / recomputes two rows of Table 9'],
  ]),
  gap(),
  p([r('Status: ', { bold: true }), r('code changes implemented and tested on synthetic data (formulas checked by hand). The pipeline has not yet been re-run on the SAFE microdata, so the new numbers for Tables 9–11 are pending. Model II is unchanged.')]),

  h1('3. Change C1 — Class weights'),
  h3('What the thesis did'),
  p('Both models weighted each training observation by w_k = median(n_1,…,n_K) / n_k, where n_k is the number of observations in bracket k (thesis §3.2.2, “Class imbalance and weighting”). Rare brackets — above all the top bracket (> EUR 1 million) — therefore received a much larger weight. The stated aim was to reduce over-allocation to the modal bracket.'),
  h3('Why this is a problem'),
  p('Weighting class k by w_k in a multinomial log-loss multiplies the odds of that class by w_k. The model’s predicted probabilities are therefore no longer estimates of P(Y = k | X): they are systematically shifted towards the rare brackets.'),
  p('This matters because the final imputation (Block 9) does not assign the most probable bracket; it draws a bracket at random from the predicted probability vector. With weighted probabilities, rare brackets are drawn too often. Because the top bracket is valued at EUR 1.5 million, this mechanically inflates the euro-valued financing gap. The weights help only when the bracket is assigned by argmax, which is not how the dataset is filled.'),
  p('Illustration on synthetic data with the thesis bracket frequencies (top bracket ≈ 3%): the average predicted probability of the top bracket rises from 3.0% without weights to 18.9% with the thesis weights, and returns exactly to 3.0% after the correction below.'),
  h3('What the code does now'),
  p('A switch CLASS_WEIGHT_MODE in Block 0 selects one of three settings:'),
  b([r('"none" (default): ', { bold: true }), r('no class weights. Predicted probabilities stay calibrated to the observed bracket frequencies, so the random draw reproduces the true bracket shares on average.')]),
  b([r('"corrected": ', { bold: true }), r('class weights as in the thesis during training, followed by a prior-shift correction of every probability vector before it is used for evaluation, hyper-parameter selection or the imputation draw:')]),
  eq('p̃_ik = (p_ik / w_k) / Σ_j (p_ij / w_j)'),
  b([r('"legacy": ', { bold: true }), r('thesis setting (weights, no correction). Used only to reproduce the thesis numbers.')]),
  p('The correction is applied consistently in the grid search (Block 6, so log-loss is computed on corrected probabilities), in the holdout evaluation (Block 8) and in the imputation (Block 9). New helper functions in Block 1: class_weight_map(), row_weights(), correct_probs().'),
  h3('Effect on results'),
  p('Under "none" or "corrected", the imputed brackets, the euro-space metrics (which use the probability-weighted expected midpoint) and the financing-gap aggregates all change. The expected direction is fewer top-bracket imputations and therefore lower euro gaps, especially in the fully imputed waves after 2022. Individual-level metrics computed on argmax labels (accuracy, QWK, macro F1, tail-risk rates) may get worse under "none", because without weights the argmax rarely selects small brackets; this is expected and is not a problem for the estimand (see C2).'),

  h1('4. Change C2 — TVD on expected shares'),
  h3('What the thesis did'),
  p('Distributional credibility (Total Variation Distance and mean shift, Table 9) was computed after assigning each holdout firm its most probable bracket (argmax). The thesis justified this choice in §4.2.1: random draws would add Monte Carlo noise and make the metrics depend on the seed.'),
  h3('Why this is a problem'),
  p('The argmax rule is not the rule used to build the dataset. TVD therefore evaluated a different object from the one that enters the financing-gap aggregates. Combined with C1, the two choices reinforce each other: the class weights make argmax-based TVD look good (they push enough firms into the rare brackets), while they distort the probabilities that the actual draw uses. The headline TVD of 0.0585 for Model I does not measure the distributional quality of the imputed data.'),
  h3('What the code does now'),
  p('The predicted share of bracket k in group g (overall, wave, country, sector, size) is the average predicted probability of that bracket across the firms in the group:'),
  eq('ŝ_gk = (1 / n_g) · Σ_{i∈g} p̂_ik        s_gk = (1 / n_g) · Σ_{i∈g} 1{k_i = k}        TVD_g = ½ · Σ_k | ŝ_gk − s_gk |'),
  p('ŝ_gk is exactly the share that the random draw produces on average, so this TVD evaluates the imputation actually used, while being deterministic and reproducible — the same property the thesis wanted from argmax. The mean shift is computed in the same way, using the expected bracket index Σ_k k · ŝ_gk.'),
  h3('Worked example'),
  p('Four holdout firms, three brackets (Low, Medium, High):'),
  table([1500, 1800, 1800, 1800, 2126], [
    ['Firm', 'P(Low)', 'P(Medium)', 'P(High)', 'Observed bracket'],
    ['1', '0.70', '0.20', '0.10', 'Low'],
    ['2', '0.60', '0.30', '0.10', 'Low'],
    ['3', '0.20', '0.60', '0.20', 'Medium'],
    ['4', '0.10', '0.30', '0.60', 'High'],
    ['Expected share (column mean)', '0.40', '0.35', '0.25', '—'],
    ['Observed share', '0.50', '0.25', '0.25', '—'],
  ]),
  gap(),
  p('TVD = ½ · (|0.40 − 0.50| + |0.35 − 0.25| + |0.25 − 0.25|) = 0.10. With argmax, the four firms would be assigned Low, Low, Medium, High, giving TVD = 0 — an apparently perfect score that ignores the model’s uncertainty, which the random draw does reproduce (firm 2 is imputed as Medium 30% of the time).'),
  h3('Implementation'),
  b('New function make_dist_prob() in Block 1; Block 8 computes both versions for all groups.'),
  b('Console output, the metrics table exported in Block 11 (CSV and LaTeX) and the exported .dta files report both: “expected_share” (primary) and “argmax” (thesis definition, kept for comparability). The .dta files have a tvd_method column.'),
  b('A table of observed vs expected vs argmax bracket shares (overall) is printed as a calibration check.'),
  b('Accuracy, within-1 accuracy, QWK, macro F1 and the tail-risk rates are still computed on argmax labels, because they measure classification of individual firms, not aggregates.'),
  h3('Effect on results and on model selection'),
  p('The distributional rows of Table 9 change. The model selection rule of §4.5 selects the model with the lowest TVD among those passing the guardrails. For Model II the snapped prediction is deterministic and is the value actually imputed, so its existing TVD is already coherent; the comparison between models should now use the expected-share TVD for Model I. The selection outcome must be re-checked once the pipeline is re-run.'),

  h1('5. Change C3 — Sector and size grouping variables'),
  h3('What the thesis did'),
  p('In the thesis script, the rows “Avg. TVD by sector” and “Avg. TVD by size” of Table 9 were computed by grouping on variables d1 and d2 respectively.'),
  h3('Why this is a problem'),
  p('In SAFE, D1 is firm size (number of employees) and D3 is the sector (thesis Chapter 2). The row labelled “by sector” was therefore computed by size class, and the row labelled “by size” by variable d2, which is not the size-class variable used elsewhere in the thesis (to be checked against the SAFE codebook). The heatmap in Block 12 already used the correct variables (d3_rec for sector, d1_rec for size).'),
  h3('What the code does now'),
  p('Sector is grouped by d3_rec (falling back to d3) and size by d1_rec (falling back to d1), consistently with the heatmap. The variable used is printed next to each row.'),
  h3('Effect on results'),
  p('The thesis value 0.1240 labelled “by sector” is actually a by-size TVD (argmax definition). Both rows must be recomputed.'),

  h1('6. What has not changed'),
  b('Sample construction, Q8A cleaning, predictor screening and feature set (Blocks 2–4). Open issues on these are listed in section 7.'),
  b('Holdout design: last 20% of observed waves (waves 26–28).'),
  b('Hyper-parameter grid, early stopping, final refit on all observed Q8A.'),
  b('Bracket midpoints, including EUR 1.5 million for the open top bracket.'),
  b('Model II (Catboost_Model_II.R), RF_Old_Method.R and Post_Imputation_Comparison.'),

  h1('7. Open issues identified (not yet implemented)'),
  p('Listed here so that future changes can refer back to them. Each will receive a change ID when implemented.'),
  table([700, 4200, 4126], [
    ['Ref.', 'Issue', 'Planned direction'],
    ['A', 'Questions discontinued or reworded after 2022 survive the >90% missingness filter and are coded “MISSING” for every post-2022 firm — a pattern never seen in training. The holdout (waves 26–28) cannot detect it.', 'Restrict predictors to questions asked identically in both periods; mask post-2022-unavailable questions in the holdout.'],
    ['B', 'wave_agg and d0 enter as predictors. Trees cannot split beyond wave 28, so all post-2022 waves receive the wave-28 effect (implicit, undocumented assumption).', 'Remove wave_agg; replace with observed macro drivers (CatboostX), with a parametric component for extrapolation.'],
    ['D', 'Holdout used for screening, grid search, early stopping and final reporting.', 'Separate validation and test blocks; rolling-origin time cross-validation.'],
    ['E', 'From 2024 the reference period of Q7A/Q7B changed from 6 to 3 months; the 2023 pilot split-sample experiment gives an application rate ratio of about 0.8.', 'Estimate a bridge factor from the pilots and apply it to 2024+.'],
    ['F', 'SAFE weights restore employment, not firm counts: Σ w·Gap is an employment-weighted mean, not a euro total.', 'Gross up with Eurostat SBS enterprise counts by country × size × sector.'],
    ['G', 'Uncertainty bands hold the model fixed (no parameter uncertainty); the script producing them is not in the repository.', 'Bootstrap refits and Rubin’s rules; add the script to the repository.'],
  ]),

  h1('8. Version history'),
  table([1300, 1700, 6026], [
    ['Version', 'Date', 'Content'],
    ['1', '25 Oct 2026', 'C1 class weights, C2 expected-share TVD, C3 sector/size grouping variables (Catboost_Model_I.R).'],
  ]),
];

const doc = new Document({
  styles: {
    default: { document: { run: { font: FONT, size: 21 } } },
    paragraphStyles: [
      { id: 'Heading1', name: 'Heading 1', basedOn: 'Normal', next: 'Normal', quickFormat: true, run: { size: 28, bold: true, color: '2E5C8A', font: FONT }, paragraph: { outlineLevel: 0 } },
      { id: 'Heading2', name: 'Heading 2', basedOn: 'Normal', next: 'Normal', quickFormat: true, run: { size: 24, bold: true, color: '2E5C8A', font: FONT }, paragraph: { outlineLevel: 1 } },
      { id: 'Heading3', name: 'Heading 3', basedOn: 'Normal', next: 'Normal', quickFormat: true, run: { size: 21, bold: true, color: '404040', font: FONT }, paragraph: { outlineLevel: 2 } },
    ]
  },
  numbering: { config: [{ reference: 'bul', levels: [{ level: 0, format: LevelFormat.BULLET, text: '•', alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 567, hanging: 283 } } } }] }] },
  sections: [{
    properties: { page: { margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [r('Changes vs thesis — page ', { size: 16, color: '777777' }), new TextRun({ children: [PageNumber.CURRENT], size: 16, color: '777777', font: FONT })] })] }) },
    children
  }]
});
Packer.toBuffer(doc).then(buf => fs.writeFileSync(process.argv[2], buf));
