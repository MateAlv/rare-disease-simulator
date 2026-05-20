# MCRD-Backed Disease Candidates

## Local Evidence Checked

This note is based on the local GraPhens checkout at:

```text
/home/mate/FIUBA/trabajo-profesional/GenPhenAI/GraPhens
```

Files checked:

```text
validation_arena/MCRD.dataset.tsv
data/genes_to_phenotype.txt
data/phenotype.hpoa
data/simulation/output/simulated_patients_all_20251202_154807.json
```

MCRD is gene-labeled, not disease-labeled:

```text
ID    Gene    Phenotype
```

Therefore, "disease present in MCRD" means: at least one causal or strongly associated gene for the disease appears in MCRD. Disease names below are inferred from local HPO/GraPhens gene-to-disease annotations, not from MCRD itself.

Summary:

```text
MCRD cases: 229
Unique MCRD genes: 196
MCRD genes found in GraPhens genes_to_phenotype.txt: 196 / 196
```

## Strong Candidate Replacements

These candidates are present in MCRD by gene, have local disease annotations, and add useful clinical variety. The selected MVP set is listed in the next section.

| Candidate disease | Gene | Local disease IDs | MCRD cases | MCRD IDs | Why useful |
| --- | --- | --- | ---: | --- | --- |
| Neuronal ceroid lipofuscinosis 6 | CLN6 | OMIM:601780, OMIM:204300 | 3 | MCRD-052, MCRD-055, MCRD-056 | Progressive neurodegenerative lysosomal/storage phenotype; good contrast with static neurodevelopmental disorders. |
| Glycogen storage disease II / Pompe disease | GAA | OMIM:232300, ORPHA:308552 | 1 | MCRD-049 | Treatable metabolic-neuromuscular disease; useful for onset and severity modeling. |
| Hypophosphatasia | ALPL | OMIM:241500, OMIM:241510, OMIM:146300 | 2 | MCRD-019, MCRD-097 | Skeletal/metabolic disease with onset variation. |
| Noonan syndrome | PTPN11 | ORPHA:648, OMIM:163950 | 2 | MCRD-007, MCRD-148 | Developmental, craniofacial, cardiac, short stature; well documented. |
| Ataxia-telangiectasia | ATM | OMIM:208900, ORPHA:100 | 2 | MCRD-063, MCRD-145 | Neurologic, immune, cancer-risk phenotype; useful multisystem profile. |
| PTEN hamartoma tumor syndrome / Cowden spectrum | PTEN | ORPHA:201, ORPHA:109, OMIM:158350 | 2 | MCRD-062, MCRD-103 | Overgrowth, neurodevelopmental, dermatologic, tumor-risk phenotype. |
| Fibrodysplasia ossificans progressiva | ACVR1 | ORPHA:337, OMIM:135100 | 1 | MCRD-174 | Highly distinctive skeletal phenotype; good contrast case. |
| Mucopolysaccharidosis type VI | ARSB | OMIM:253200 | 1 | MCRD-188 | Lysosomal storage disease with skeletal, ocular, cardiac features. |
| Autosomal dominant optic atrophy | OPA1 | ORPHA:98673, OMIM:165500 | 1 | MCRD-131 | Sensory/ophthalmologic phenotype. |
| Gitelman syndrome | SLC12A3 | ORPHA:358, OMIM:263800 | 1 | MCRD-192 | Renal/electrolyte disorder; non-neurodevelopmental contrast. |
| Juvenile neuronal ceroid lipofuscinosis / CLN3 disease | CLN3 | OMIM:204200 | 1 | MCRD-163 | Neurodegenerative + visual phenotype; complements CLN6 or can replace it. |
| KID syndrome / GJB2-related disease | GJB2 | ORPHA:477, OMIM:148210 | 1 | MCRD-096 | Deafness/skin/ocular phenotype; possible sensory replacement, but gene-disease mapping is broader. |

## Selected MCRD-Backed MVP Set

The MVP now uses this MCRD-backed 8-disease set:

| Disease | Gene | MCRD evidence |
| --- | --- | --- |
| Niemann-Pick disease type C1 | NPC1 | MCRD-101 |
| Rett syndrome / atypical Rett syndrome | MECP2 | MCRD-026, MCRD-112, MCRD-198, MCRD-221 |
| Duchenne muscular dystrophy | DMD | MCRD-110 |
| Cystic fibrosis | CFTR | MCRD-213 |
| Neuronal ceroid lipofuscinosis 6 | CLN6 | MCRD-052, MCRD-055, MCRD-056 |
| Pompe disease / glycogen storage disease II | GAA | MCRD-049 |
| Hypophosphatasia | ALPL | MCRD-019, MCRD-097 |
| Noonan syndrome | PTPN11 | MCRD-007, MCRD-148 |

The set is intentionally balanced across neurodevelopmental, neuromuscular, respiratory/digestive, neurodegenerative, metabolic, skeletal, and syndromic presentations.

## Full MCRD Gene List

The local MCRD file contains these 196 unique genes:

```text
AAAS, ABCC9, ABCD1, ACO2, ACTA1, ACTB, ACTG1, ACTG2, ACVR1, ADCY5,
AEBP1, AFG3L2, AHDC1, ALDH18A1, ALG12, ALPL, AMPD2, ANKRD11, ANO3,
ARID1B, ARSB, ASH1L, ASXL1, ASXL3, ATM, ATP1A3, B4GALNT1, BCL11A,
BMP4, BRCA2, BRPF1, C12orf57, C2CD3, CACNA1A, CACNA1E, CACNA1F,
CASR, CBL, CFTR, CHAT, CHD2, CHD3, CHD7, CIC, CLN3, CLN6, CNOT1,
COL2A1, COL4A2, COL4A4, CREBBP, CTCF, CTLA4, CUL3, CUL4B, CYP11A1,
CYP7B1, DDX3X, DEAF1, DEGS1, DGAT1, DMD, DOCK3, DOCK8, EBF3,
EFTUD2, EHMT1, ENPP1, EP300, ERCC6L2, EXOSC3, EZH2, FA2H, FAR1,
FBXO11, FGFR1, FLG, FOXP1, GAA, GABBR2, GABRA1, GFAP, GJA1, GJB2,
GLMN, GNPTAB, GRIN1, HCFC1, HMGCS2, HNRNPU, HTRA1, IDH1, IFIH1,
JAK2, KANSL1, KAT6A, KCNA2, KCNB1, KDM6A, KIF11, KMT2A, KMT2B,
KMT2D, LGI1, LOX, LTBP3, LZTR1, MAP3K7, MECP2, MEFV, MPZ, MSTO1,
MUC1, MYH11, NAA15, NDUFS4, NEB, NEFL, NFKB1, NKX2-1, NONO, NPC1,
NR2F1, NR2F2, NSD1, OCA2, OPA1, PACS1, PANK2, PBX1, PGAP3, PIEZO1,
PIGV, PIK3R1, PKD1, PLA2G6, POGZ, POLR3A, POMGNT1, PPP3CA, PSEN1,
PTEN, PTPN11, PYCR2, RECQL4, RPS6KA3, RUNX1, RYR1, SCN2A, SERPINA1,
SETD5, SHANK3, SIN3A, SLC12A3, SLC20A2, SLC25A46, SLC26A3, SLC35C1,
SLC6A1, SLC6A8, SMAD4, SON, SOX2, SPG7, STAG1, SYNGAP1, TBR1, TBX4,
TCF4, TGFBR2, THAP1, THOC6, TIA1, TNFRSF1A, TRAPPC9, TRIO, TRIP12,
TRMU, TRPS1, TRPV4, TRRAP, TSPEAR, TTN, TUBA1A, TUBB4A, UBA1,
UBE2A, UBE3A, UMOD, USP9X, WAC, XPC, ZBTB18, ZC4H2, ZEB2, ZSWIM6
```
