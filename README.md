# PhD_work

## Project Description
CD4+ T helper (Th) cells dictate host defense responses against immunological challenges, while self-reactive Th progenitor cells that escape from negative selection are kept in check from causing autoimmune diseases. How Th cell immunity and tolerance is regulated remains one of the most fundamental questions in immunology. 

We are interested in studying how the TGF-β pathway and the PD-1 pathway cooperatively regulate peripheral Th cell tolerance. We found mice lacking TGF-β receptor 2 (TGFBR2) in mature CD4+ T cells on the background of PD-1 deficiency (Thpokcre Tgfbr2fl/fl Pdcd1-/-: 2KO) develop lethal autoimmunity at around 3-4 week of age, while single deficiency of each pathway did not trigger autoimmunity. 

Here, we used single-cell RNA- and TCR-sequencing to profile the Th cell population in the 2KO mice compared to the control groups (WT, Thpokcre Tgfbr2fl/fl, Pdcd1-/-) in lymph nodes and peripheral tissues. We also profiled the MHC-II+ antigen presenting cells in the tissue to characterize the effector APCs that interact with Th cells to contribute to the immunopathology.

## Conclusions
In the 2KO mice, lacking both TGF-β and the PD-1 inhibitions enhanced priming and reactivation of conventional CD4+ T cells in lymph nodes and peripheral tissues, respectively, with PD-1 deficiency exacerbating T cell clonal expansion and differentiation into a population of hyperactivated T helper 1 (Th1) and follicular helper T (Tfh) cells. Furthermore, the excessive Th cell response and the ensuing autoimmunity was dependent on CD40, which was highly expressed in a population of monocyte-derive macrophages. 

Based on these findings, we propose a model that TGF-β and the PD-1 serve as ‘tonic’ and ‘reactive’ regulators of peripheral Th cell tolerance, respectively, by inhibiting cellular immunity driven by CD40-expressing myeloid APCs. We also demonstrated that targeting the TGF-β pathway specifically on CD4+ T cells in combination with PD-L1 blockade significantly improve antitumor responses in a low mutational burden hepatocellular carcinoma mouse model, showing that such a Th cell-intrinsic control module can be targeted to revive protective immunity in cancer and chronic infections.

## Experiment design
### CD4 T cell scRNA scTCR
- library 1: samples are sorted CD4+ T cells from the liver draining lymph node of 3-week old mice from the following 5 genotypes, 1 mouse per genotype, hash-tagged

| Genotype                             | Abbreviation | Phenotype              |
|--------------------------------------|--------------|------------------------|
| wild-type                            | WT           | healthy                |
| Thpokcre Tgfbr2fl/fl                 | RII          | healthy                |
| Pdcd1-/-                             | PD1          | healthy                |
| Thpokcre Tgfbr2fl/fl Pdcd1-/-        | 2KO          | lethal autoimmunity    |
| Thpokcre Tgfbr2fl/fl Pdcd1-/- Cd40-/-| 3KO          | healthy (rescue)       |

- library 2: samples are sorted CD4+Cd1d- T cells from the matching liver of 3-week old mice from the following 5 genotypes, 1 mouse per genotype (except RII), hash-tagged

| Genotype                             | Abbreviation | Phenotype              |
|--------------------------------------|--------------|------------------------|
| wild-type                            | WT           | healthy                |
| Thpokcre Tgfbr2fl/fl                 | RII_1 (few cells)| healthy            |
| Thpokcre Tgfbr2fl/fl                 | RII_2        | healthy                |
| Pdcd1-/-                             | PD1          | healthy                |
| Thpokcre Tgfbr2fl/fl Pdcd1-/-        | 2KO          | lethal autoimmunity    |
| Thpokcre Tgfbr2fl/fl Pdcd1-/- Cd40-/-| 3KO          | healthy (rescue)       |

### MHCII APC scRNA
- library 3: samples are sorted MHCII+Lineage- (Tcrb-NK1.1-CD90-CD19-SiglecF-Ly6G-) APCs from the liver of 3-week old mice from the following 5 genotypes, 2 mice per genotype, hash-tagged

| Genotype                             | Abbreviation | Phenotype              |
|--------------------------------------|--------------|------------------------|
| wild-type                            | WT           | healthy                |
| Thpokcre Tgfbr2fl/fl                 | RII          | healthy                |
| Pdcd1-/-                             | PD1          | healthy                |
| Thpokcre Tgfbr2fl/fl Pdcd1-/-        | 2KO          | lethal autoimmunity    |
| Thpokcre Tgfbr2fl/fl Pdcd1-/- Cd40-/-| 3KO          | healthy (rescue)       |
