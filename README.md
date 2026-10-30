# Muscle Contractility Lab Analysis: DMD (Dystrophin)

## 1. Metadata & Assignment Overview
* **Student Assignment Number:** 37
* **Gene Symbol:** `DMD` (Dystrophin)
* **NCBI Gene ID:** 1756
* **Taxonomy ID:** `txid9606` (*Homo sapiens*)
* **Assigned RefSeq Transcript Accession:** `NM_004006.3`
* **Retrieved RefSeq Protein Accession:** `NP_003997.2`
* **Fixed ClinVar Variation ID:** `11209`
* **Assigned Variant (HGVS):** `c.2791G>T` (`p.Glu931Ter`)
* **Repository Name:** `Muscle-Contractility-DMD`

---

## 2. Research Question
How do structural variations and nonsense mutations in the *DMD* gene impact the full-length dystrophin protein product, and what are the downstream molecular mechanisms connecting premature stop codons to sarcolemmal instability, intracellular calcium dysregulation, and muscular dystrophy pathology?

---

## 3. Step-by-Step Methods & R Code Execution

### Step 3.1: Project Directory & Workspace Setup
```R
# Install required packages (run once if needed)
install.packages(c("rentrez", "ggplot2"))
library(rentrez)
library(ggplot2)

# Define assignment variables
gene_symbol <- "DMD"
transcript_accession <- "NM_004006.3"
clinvar_id <- "11209"
variant_protein <- "p.Glu931Ter"

# Create directory structure
for (x in c("data", "results", "figures")) {
  dir.create(x, showWarnings = FALSE)
}
```

---

### Step 3.2: NCBI Gene Record Retrieval
```R
hits <- entrez_search(
  db = "gene",
  term = paste(gene_symbol, "[Gene Name] AND txid9606[Organism]"),
  retmax = 20
)

if (length(hits$ids) == 0) stop("Gene not found")

gene_id <- NA_character_
for (id in hits$ids) {
  rec <- entrez_summary(db = "gene", id = id)
  if (toupper(rec$name) == toupper(gene_symbol)) {
    gene_id <- id; gene_info <- rec; break
  }
}

if (is.na(gene_id)) stop("Exact human gene not found")

writeLines(c(
  paste("Gene symbol:", gene_info$name),
  paste("NCBI Gene ID:", gene_id),
  paste("Official name:", gene_info$description),
  paste("Chromosome:", gene_info$chromosome)
), "results/gene_summary.txt")
```

---

### Step 3.3: RefSeq Transcript & Protein FASTA Retrieval
```R
gb <- entrez_fetch(
  db = "nuccore", id = transcript_accession,
  rettype = "gb", retmode = "text"
)

pat <- '/protein_id="NP_[0-9]+[.][0-9]+"'
fragment <- regmatches(gb, regexpr(pat, gb))
if (length(fragment) == 0) stop("NP protein ID not found")
protein_accession <- gsub('^/protein_id="|"$', "", fragment)

protein_fasta <- entrez_fetch(
  db = "protein", id = protein_accession,
  rettype = "fasta", retmode = "text"
)

if (!startsWith(protein_fasta, ">")) stop("Invalid FASTA")

fasta_file <- file.path("data", paste0(gene_symbol, "_protein.fasta"))
writeLines(protein_fasta, fasta_file)
writeLines(paste("NP accession:", protein_accession), "results/protein_accession.txt")
```

---

### Step 3.4: Amino Acid Composition Calculation & Visualization
```R
fasta_lines <- readLines(fasta_file)
sequence <- paste(fasta_lines[!grepl("^>", fasta_lines)], collapse = "")
aa <- strsplit(toupper(sequence), "")[[1]]
protein_length <- length(aa)

standard <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
counts <- table(factor(aa, levels = standard))
results <- data.frame(Amino_Acid = standard, Count = as.integer(counts))

write.csv(results, "results/amino_acid_composition.csv", row.names = FALSE)

# Generate Figure 1
fig1 <- ggplot(results, aes(Amino_Acid, Count)) +
  geom_col(fill = "grey45") + theme_bw() +
  labs(title = paste(gene_symbol, "amino-acid composition"),
       x = "Amino acid", y = "Number of residues")

ggsave("figures/01_amino_acid_composition.png", fig1, width = 8, height = 5, dpi = 300)
```

---

### Step 3.5: Protein Variant Mapping
```R
if (nzchar(variant_protein)) {
  matched <- regmatches(variant_protein, regexpr("[0-9]+", variant_protein))
  position <- as.integer(matched)
  
  if (!is.na(position) && position >= 1 && position <= protein_length) {
    fig2 <- ggplot() +
      geom_segment(aes(x = 1, xend = protein_length, y = 0, yend = 0), linewidth = 2) +
      geom_point(aes(x = position, y = 0), size = 4, color = "red") +
      annotate("text", x = position, y = 0.25, label = variant_protein, size = 4) +
      xlim(0, protein_length) + ylim(-0.6, 0.6) +
      theme_void() +
      labs(title = paste(gene_symbol, "reference protein position"), 
           x = "Amino-acid position") +
      theme(axis.title.x = element_text())
    
    ggsave("figures/02_variant_position.png", fig2, width = 8, height = 2.5, dpi = 300)
  }
}
```

---

## 4. Results Summary

### 4.1 Quantitative Sequence Metrics
* **Official Gene Name:** Dystrophin (*DMD*)
* **NCBI Gene ID:** 1756 (Chromosome X)
* **Matched RefSeq Protein:** `NP_003997.2` (Dystrophin isoform Dp427m)
* **Total Protein Length:** 3,685 amino acids
* **Top 3 Most Abundant Amino Acids:**
  1. **Leucine (L):** 455 residues (12.35%)
  2. **Glutamic Acid (E):** 387 residues (10.50%)
  3. **Glutamine (Q):** 306 residues (8.30%)

### 4.2 Generated Visualizations

#### Figure 1: Amino Acid Composition Bar Chart
<img width="2400" height="1500" alt="12026" src="https://github.com/user-attachments/assets/9d90cb81-3732-4551-8a4c-e51075279252" />

* *Interpretation:* The distribution demonstrates a prominent overrepresentation of charged and polar residues (Leucine, Glutamic Acid, Lysine, and Glutamine), which structural studies show form the characteristic coiled-coil alpha-helical structural motifs in dystrophin's 24 spectrin-like repeat domain.

#### Figure 2: Linear Variant Position Plot
<img width="748" height="1822" alt="12007" src="https://github.com/user-attachments/assets/a570a00c-c8e9-4f18-ba09-710d7892918f" />

* *Interpretation:* Highlights the relative location of residue position 931 along the 3,685-residue polypeptide chain, showing that the `p.Glu931Ter` mutation occurs early in the central rod domain.

#### Figure 3: Annotated Protein Domain Architecture
<img width="1376" height="768" alt="12028" src="https://github.com/user-attachments/assets/b6276411-4155-4a62-8bf2-5b836e95aeb9" />

* *Interpretation:* Structural schematic displaying dystrophin's functional domains: N-terminal Actin-Binding Domain (ABD, aa 1–240), Central Rod Domain with 24 spectrin repeats (aa 253–3040), Cysteine-Rich Domain (CRD, aa 3080–3360), and C-terminal Domain (CTD, aa 3361–3685). The red marker at residue 931 illustrates that truncation prevents synthesis of all subsequent rod repeats, CRD, and CTD domains.

---

## 5. Answers to Laboratory Questions (Section V)

### Q1: What does your assigned gene encode, and how does this protein contribute to muscle function?
The *DMD* gene encodes **dystrophin** (isoform Dp427m), a massive 427 kDa cytoskeletal protein localized at the inner surface of the sarcolemma (cell membrane) in skeletal and cardiac muscle fibers. Dystrophin serves as a essential mechanical shock absorber: its N-terminus binds intracellular F-actin filaments, while its C-terminal region binds the transmembrane dystrophin-glycoprotein complex (DGC). This bridge transfers contractile force generated by myofibrils to the extracellular matrix, stabilizing the fragile sarcolemma during cycles of muscle contraction and relaxation.

### Q2: What was the exact protein accession and length retrieved with R?
* **Protein Accession:** `NP_003997.2`
* **Total Protein Length:** 3,685 amino acids.

### Q3: What did the amino-acid composition chart show? What can it NOT tell you?
* **What it showed:** The bar chart (`01_amino_acid_composition.png`) showed an abundance of Leucine (455), Glutamic Acid (387), Glutamine (306), and Lysine (295). This high concentration of polar/charged amino acids supports the formation of extended triple-helical coiled-coil motifs characteristic of spectrin-like repeats in the central rod domain.
* **What it CANNOT tell you:** Amino acid composition alone cannot reveal the three-dimensional folding pattern, quaternary interactions, domain spatial coordinates, or functional clinical severity resulting from specific point mutations.

### Q4: What is your fixed ClinVar variant, and what is its current classification and evidence strength?
* **Fixed Variant:** `NM_004006.3:c.2791G>T` (`p.Glu931Ter`), ClinVar Variation ID: `11209`.
* **Clinical Classification:** **Pathogenic**.
* **Evidence Strength:** High (3 stars / multiple submitters with assertion criteria provided and no conflicting classifications). Classified by major genetic clinical centers (e.g., GeneDx, OMIM) and supported by functional literature (PMID: 1549595, PMID: 25525159).

### Q5: What is the scientifically supported mechanism linking this variant to muscle function, and what remains uncertain?
* **Supported Mechanism:** The nucleotide substitution (`c.2791G>T`) introduces a premature stop codon early in the reading frame (`p.Glu931Ter`). This triggers **Nonsense-Mediated mRNA Decay (NMD)**, preventing full-length dystrophin protein translation. Without dystrophin, the physical link between F-actin and the sarcolemma is lost, causing sarcolemmal fragility, contraction-induced microtears, uncontrolled extracellular $\text{Ca}^{2+}$ influx, and progressive skeletal muscle fiber necrosis leading to Duchenne Muscular Dystrophy (DMD).
* **What Remains Uncertain:** It is uncertain whether a small fraction of mRNA evades NMD to produce a truncated peptide, or whether alternative internal translation initiation/exon skipping mechanisms can naturally modulate phenotypic severity in individual patients.

---

## 6. Limitations & Proposed Experimental Validation

### 6.1 Limitations of Computational Analysis
* **Compositional Limitations:** Global amino acid frequencies indicate physical composition (e.g., high alpha-helical structural potential) but cannot predict tertiary folding, domain stability, or clinical disease severity.
* **Linear Map Limitations:** The linear plot marks the sequence coordinate of a stop codon but does not verify whether a truncated protein is actually produced or degraded in vivo.

### 6.2 Proposed Validation Experiment
To verify whether `c.2791G>T` leads to complete loss of expression via NMD versus production of an unstable truncated fragment:
1. **Model System:** Differentiate patient-derived induced pluripotent stem cells (iPSCs) harboring the `c.2791G>T` mutation into functional skeletal myotubes.
2. **Quantitative RT-qPCR:** Measure *DMD* mRNA transcript abundance relative to wild-type controls to measure NMD efficiency.
3. **Western Blotting:** Probe myotube lysates with N-terminal and C-terminal specific anti-dystrophin antibodies to confirm total protein absence versus truncated peptide accumulation.
4. **Immunofluorescence Microscopy:** Co-stain myotubes for N-terminal dystrophin and $\beta$-dystroglycan to confirm whether sarcolemmal localization of the DGC complex is completely lost.

---

## 7. Database References & Citations

1. **NCBI Gene Database:** [Gene ID 1756 (DMD)](https://www.ncbi.nlm.nih.gov/gene/1756) (Accessed: Oct 30, 2026).
2. **NCBI RefSeq Nucleotide:** [NM_004006.3 Record](https://www.ncbi.nlm.nih.gov/nuccore/NM_004006.3) (Accessed: Oct 30, 2026).
3. **NCBI RefSeq Protein:** [NP_003997.2 Record](https://www.ncbi.nlm.nih.gov/protein/NP_003997.2) (Accessed: Oct 30, 2026).
4. **NCBI ClinVar Database:** [Variation ID 11209](https://www.ncbi.nlm.nih.gov/clinvar/variation/11209/) (Accessed: Oct 30, 2026).
5. **Roberts RG, et al. (1992).** Direct detection of point mutations in the dystrophin gene. *Proceedings of the National Academy of Sciences USA*, 89(6), 2331-2335. PMID: [1549595](https://pubmed.ncbi.nlm.nih.gov/1549595/).
6. **Xiong HY, et al. (2015).** The human splicing code reveals new insights into the genetic determinants of disease. *Science*, 347(6218), 1254806. PMID: [25525159](https://pubmed.ncbi.nlm.nih.gov/25525159/).
