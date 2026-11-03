# Haematopoietic stem and progenitor cell single-cell qRT-PCR data

## Source

Downloaded on 2026-11-02 from the Göttgens lab data page,
http://blood.stemcells.cam.ac.uk/single_cell_qpcr.html, "Profiles of primary bone marrow cells":

- `hspc_qpcr_raw.csv`: http://blood.stemcells.cam.ac.uk/qpcr_data/hspc_qpcr_raw.csv
- `hspc_qpcr_processed.csv`: http://blood.stemcells.cam.ac.uk/qpcr_data/hspc_qpcr_processed.csv

The two CSVs are unchanged copies. The same page also offers profiles of the 416B and HoxB8-FL cell lines,
which are not included here.

The data are from Hamey et al. 2017, "Reconstructing blood stem cell regulatory network models from single-cell
molecular profiles", PNAS 114:5822-5829, doi:10.1073/pnas.1610609114 (main text and SI in `documents/`). That
paper extends the data set of Wilson et al. 2015 (Cell Stem Cell). According to the paper:

- Bone marrow was taken from the femurs, tibiae and iliac crest of 8- to 12-week-old C57BL/6 mice. The paper
  does not say how many mice were used, or which mouse each cell came from.
- Cells were lineage-depleted and sorted by surface markers into 12 stem and progenitor populations.
- Expression of 48 genes was measured by single-cell qRT-PCR on the Fluidigm Biomark: 33 transcription factors,
  12 other genes involved in stem and progenitor cell biology, and 3 housekeeping genes.

Hamey et al.'s HSC -> LMPP branch (531 cells, 31 transcription factors, binarized, in pseudotime order) is a
subset of these cells. It is in `data/hamey_pseudotime_inference`, copied from
github.com/fionahamey/Pseudotime-network-inference, and all 531 cells appear in the processed file.

## Files

Both CSVs have bare carriage-return line endings, and their first two columns are in opposite orders.

`hspc_qpcr_raw.csv`: `Name,Cell_type,<48 genes>`, raw Ct values (cycles), from about 5 to 40, where 40 means
not detected. It has 2,971 rows: 381 control wells and 2,590 cells. Controls have `NTC` (no template), `10X` or
`20X` (10- or 20-cell controls) in their names. The cells include those later removed by quality control.

`hspc_qpcr_processed.csv`: `Cell_type,Name,<48 genes>`, for the 2,167 cells that passed quality control. Values
are delta Ct: each gene's raw Ct minus a per-cell offset, so values are on a log2 scale and lower means more
expression. The SI says the normalization is to the housekeeping genes *Ubc* and *Polr2a*. Checked on this file,
raw minus processed is exactly one constant per cell across genes, and roughly the mean Ct of those two genes.
Undetected values are replaced by one fixed number per gene, between 8.3 and 16.9. Detected values range from
-12.3 to 13.0. Cells per population:

| HSC1 | HSC2 | HSC3 | HSC4 | ST-HSC1 | ST-HSC2 | MPP | LMPP | CMP | GMP | PreMegE | MEP |
|------|------|------|------|---------|---------|-----|------|-----|-----|---------|-----|
| 198  | 166  | 197  | 198  | 233     | 199     | 188 | 178  | 147 | 185 | 154     | 124 |

The SI also says *Egfl7*, *Gfi1* and *Sfpi1* were removed from Hamey et al.'s analysis for technical reasons.
They are still in both files.

## Derived file

`continuous_expression_hspc.txt` is written by `single_cell_inference/normalize_hamey_qpcr.py` (run from the
repository root; it downloads the two CSVs if they are missing). The file is tab separated in the layout of
`binary_expression_LMPP.txt`: a header row of genes, then a cell name and one value per gene. It covers all 2,167
cells and 48 genes.

Each gene is scaled to [0, 1] over all 2,167 cells:

- A detected cell (raw Ct < 40) gets (max dCt - dCt) / (max dCt - min dCt), with the max and min taken over
  the gene's detected cells. So the most expressing cell is 1, and the least expressing detected cell is 0.
- An undetected cell gets 0.

The script also checks agreement with Hamey et al.'s binary LMPP matrix (531 cells x 31 genes):

- Thresholding at 0.5 agrees on 94.5% of entries.
- Detection alone (value > 0) agrees on 98.3%.

So the published binarization is essentially a detection call. *Runx1* (44%) and *Notch* (74%) disagree most at
the 0.5 threshold.
