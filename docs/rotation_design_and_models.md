---
title: "ACTIVATE Lentil-Wheat Rotation Experiment"
author: ""
date: "`r Sys.Date()`"
output: html_document
---

# ACTIVATE Lentil-Wheat Rotation Experiment: Design, Modeling Challenges, and Analytical Targets

## Objectives

`LegumeRotation` is an R package designed to support analysis of a two-year lentil-wheat rotation experiment. The package has four connected analytical goals:

1.  Characterize phenotypic performance of lentil and wheat as independent field trials.
2.  Test whether the lentil genotype grown in Year 1 influences the agronomic performance of wheat grown in the same plot in Year 2.
3.  Identify wheat genotypes that are better rotational partners after lentil, based on yield, biomass, nitrogen uptake, nitrogen partitioning, protein, and residue-related traits.
4.  Estimate whether specific observed lentil-wheat genotype combinations show additional rotational value beyond their additive lentil and wheat genotype effects.

The main rotation questions are therefore:

> Which lentil genotypes create positive or negative legacy effects for the following wheat crop?

> Which wheat genotypes are better able to capture, convert, or partition lentil legacy effects into useful agronomic and nutrient outcomes?

> Which observed lentil-wheat genotype pairs perform better or worse than expected from the average behavior of the lentil and wheat genotypes separately?

## Experimental Design

The field layout contains 76 rows by 30 columns, for 2,280 plots per trial layout, organized into 48 incomplete blocks. The experiment includes:

-   100 lentil genotypes in Year 1.
-   100 wheat genotypes in Year 2.
-   20 observations per genotype per location for the full-trial phenotypic layers.
-   A structured incomplete lentil-wheat factorial design.

### Structured Incomplete Factorial Design

The experiment uses a structured incomplete factorial design. A full 100 lentil × 100 wheat factorial would require 10,000 genotype combinations before replication and locations, which is not feasible. Instead, the design is divided into 10 wheat-partner facets.

Each facet contains 10 wheat genotypes and 10 lentil genotypes. Within a facet, all 10 × 10 = 100 lentil-wheat combinations are evaluated. Therefore, the design is globally incomplete relative to the full 100 × 100 factorial, but locally complete within each facet.

For example, the first wheat facet contains W001-W010. These wheats are paired with L001, L011, L021, L031, L041, L051, L061, L071, L081, and L091, giving 100 observed lentil-wheat combinations. The second wheat facet contains W011-W020 and is paired with L002, L012, L022, L032, L042, L052, L062, L072, L082, and L092. This logic continues until the final facet, which contains W091-W100 paired with L010, L020, L030, L040, L050, L060, L070, L080, L090, and L100.

Thus, `Facet` defines a local lentil-wheat network. Legacy effects should primarily be interpreted within ENV × Facet, because lentils in different facets were not evaluated against the same wheat genotypes.

### Facet Definition

The wheat file stores the local network in the `Facet` column. Facet labels such as `W001_W010` identify the 10-wheat partner set used for that local factorial. In the intended ACTIVATE design, each ENV × Facet contains:

``` text
10 lentil genotypes × 10 wheat genotypes = 100 observed combinations
```

### Expected Lentil and Wheat Membership by Facet

The lentil membership follows a modulo pattern:

``` text
Facet 1:  W001-W010 paired with L001, L011, L021, ..., L091
Facet 2:  W011-W020 paired with L002, L012, L022, ..., L092
Facet 3:  W021-W030 paired with L003, L013, L023, ..., L093
...
Facet 10: W091-W100 paired with L010, L020, L030, ..., L100
```

Wheat membership follows blocks of ten: W001-W010 are facet 1, W011-W020 are facet 2, and so on through W091-W100 in facet 10.

### Replication Terms: Rep_gen vs Rep_combo

The dataset contains two replicate identifiers that should not be used interchangeably:

-   `Rep_gen` is the genotype-level replicate. It is appropriate for independent lentil or wheat genotype performance models, where the same genotype is represented many times within an environment.
-   `Rep_combo` is the lentil × wheat combo replicate. It is the correct replicate design term for rotation analyses, because each observed `Combo` is replicated twice within an environment.

For the rotation predecessor and pair-compatibility models, use `Rep_combo` as the replicate design term. Do not use `Rep_gen` in these rotation models, because it describes genotype-level replication rather than the two replicated plots of the same observed lentil-wheat combination.

### Why Legacy Effects Are Interpreted Within ENV × Facet

Because each facet has its own wheat partner set, `Legacy_Value` is a within-facet deviation. It compares a lentil genotype against the other lentils tested with the same 10 wheat genotypes in the same environment. A cross-facet ranking of `Legacy_Value` is therefore a ranking of within-facet deviations, not a ranking from a fully connected 100 × 100 factorial.

## Data Layers

The package is designed around five major data layers.

### Lentil Full Trial

The full lentil phenotypic trial includes traits such as:

``` text
DTE, DTF, VegP, DTM, RepP, lodging, YLD, PRT, DS
```

This layer is used to estimate lentil genotype performance, spatial trends, trait distributions, BLUPs or adjusted means, and trait correlations.

### Lentil Subsample Layer

The subsample layer includes more detailed biomass and nutrient traits, such as:

``` text
biomass, straw, seed weight, NIR, LECO C/N/P, NHI, C:N ratios
```

This layer has a lower-density sampling structure ([2 reps × genotype + checks] × site), with treatment observations plus repeated check observations nested in blocks. It is analytically useful, but its design differs from the full phenotypic layer, so checks and block structure need careful handling.

### Lentil Microbiome Layer

The microbiome layer comes from the Year 1 lentil plots and includes plot-level 16S-derived predictors such as:

``` text
Observed ASV richness, Chao1, ACE, Shannon, Simpson, Pielou evenness,
Aitchison PCoA axes, microbiome legacy score, Rhizobiaceae abundance,
Bradyrhizobiaceae abundance, legume-symbiont abundance
```

For the legacy-driver analysis, these plot-level microbiome traits are aggregated to lentil genotype × environment means, then the environment labels are shifted from 2024 to the matching 2025 wheat environment. This lets the microbiome layer be treated as another Year 1 predictor of Year 2 wheat legacy.

### Wheat Full Trial

The wheat layer includes Year 2 traits such as:

``` text
HD, HT, MAT, LD, YLD, Y.ADJ, TWT, KWT, PRO
```

This layer contains the current wheat genotype, the previous lentil genotype, the observed lentil-wheat combination, the wheat facet, spatial coordinates, block information, and agronomic response traits.

### Wheat Subsample Layer

The wheat subsample layer includes detailed biomass and nutrient traits measured on Year 2 wheat plants or grain subsamples, such as:

``` text
grain biomass, straw biomass, total aboveground biomass,
seed N %, straw N %, seed C %, straw C %,
grain N content, straw N content, total aboveground N,
grain C content, straw C content, total aboveground C,
N harvest index, harvest index, C:N ratios, protein yield
```

This layer is used to estimate wheat rotational partner value and G×G rotational value for mechanistic rotation traits. It allows the workflow to move beyond whether wheat yield increased after lentil and ask whether wheat genotypes differ in their ability to capture nitrogen, partition nitrogen into grain, maintain protein, or return useful biomass and nitrogen in straw.

## Main Modeling Challenges

### 1. Spatial Heterogeneity

The field has row and column structure. Yield and related traits can vary across the field because of soil gradients, moisture, fertility, management effects, or other spatial patterns. If spatial trends are ignored, genotype effects can be confounded with field position.

The package uses SpATS spatial models where possible:

``` r
trait ~ fixed design terms + genotype effect + PSANOVA(Row, Col)
```

The spatial surface helps separate genotype signal from field heterogeneity.

### Spatial Model Scope

Spatial correction differs across the current model layers. The predecessor and wheat rotational partner models can use SpATS `PSANOVA(Row, Col)` surfaces. The pair-compatibility and G×G rotational models use the lme4 implementation with polynomial row and column terms plus their interaction, because those models require the observed `Combo` random effect and optional design random effects such as `Block`.

This difference is practical rather than biological. PSANOVA is more flexible, but it can be underpowered when `fit_scope = "env_facet"` because each local 10 × 10 network may contain only about 200 plots. Users should verify that each ENV × Facet subset spans enough rows and columns for a spatial surface to be meaningful, inspect SpATS convergence and spatial effective dimensions using `fit$eff.dim`, and compare against simpler row/column corrections. If the within-facet spatial surface is poorly supported, consider using linear row and column terms within the facet model, or using `fit_scope = "env_global"` to estimate the spatial correction from the full environment before centering corrected values by facet.

### 2. Structured Incomplete Lentil-Wheat Pairing

Not every lentil is paired with every wheat across the full experiment. This creates a globally incomplete genotype network. However, the design is locally complete within each ENV × Facet: the 10 lentils in that facet are each evaluated against the same 10 wheat genotypes.

This is why wheat genotype and facet structure must be included in the analysis. Within a facet, wheat genotype correction compares lentils against a common local wheat-partner set. Across facets, raw or corrected means should not be interpreted as if every lentil had been tested against all 100 wheat genotypes.

### 3. Facet-Specific Baselines

Because the pairing network is structured by `Facet`, comparing every lentil against the full environment mean can be misleading. A lentil in the W001-W010 facet should be compared against the expected performance of that ENV × Facet network, not against the mean of all 100 wheat genotypes.

The current preferred predecessor-effect output therefore reports:

$$
\mathrm{Legacy\_Value}
= \mathrm{Corrected\_Mean}
- \operatorname{mean}_{\mathrm{ENV}\times\mathrm{Facet}}\left(\mathrm{Corrected\_Mean}\right)
$$

It also keeps a global environment deviation for comparison:

$$
\mathrm{Legacy\_Value\_Global}
= \mathrm{Corrected\_Mean}
- \operatorname{mean}_{\mathrm{ENV}}\left(\mathrm{Corrected\_Mean}\right)
$$

The facet-aware value is the preferred interpretation for the structured incomplete design.

### 4. Pair Compatibility Is Not the Same as Lentil Legacy

Average predecessor effect and pair compatibility answer different biological questions.

Average predecessor effect asks:

``` text
Which lentil genotypes tend to leave better wheat-growing conditions?
```

Pair compatibility asks:

``` text
Which observed lentil-wheat combinations perform better than expected from their separate lentil and wheat effects?
```

Because the design is globally incomplete, pair compatibility should be interpreted within the observed local 10 × 10 network. It should not be treated as evidence for unobserved lentil-wheat combinations outside that ENV × Facet.

## Design Validation

Before fitting predecessor-effect models, pair-compatibility models, plotting functions, or GWAS-ready legacy exports, the observed data should be checked against the intended structured incomplete factorial design.

### Check Local 10 × 10 Completeness

Use:

``` r
design_summary <- check_rotation_design(
  data = wheat_treatment,
  env_col = "ENV",
  facet_col = "Facet",
  lentil_col = "Lentil",
  wheat_col = "Wheat",
  combo_col = "Combo"
)
```

For a clean treatment design, each ENV × Facet should usually report:

``` text
N_Previous_Genotypes = 10
N_Current_Genotypes  = 10
N_Combos             = 100
Expected_10x10       = TRUE
```

The helper also reports `Expected_Combos`, `Missing_Combos`, and `Is_Complete_Local_Factorial`. These diagnostics guard against accidental missing combinations, extra check plots, or incorrectly assigned facets.

### Check Lentil-to-Facet Assignment

Use:

``` r
lentil_assignment <- check_lentil_facet_assignment(
  data = wheat_treatment,
  lentil_col = "Lentil",
  facet_col = "Facet"
)
```

The expected lentil pattern is:

``` text
Facet 1:  L001, L011, L021, ..., L091
Facet 2:  L002, L012, L022, ..., L092
...
Facet 10: L010, L020, L030, ..., L100
```

The `Facet_Assignment_OK` column should be `TRUE` for all treatment genotypes.

### Check Wheat-to-Facet Assignment

Use:

``` r
wheat_assignment <- check_wheat_facet_assignment(
  data = wheat_treatment,
  wheat_col = "Wheat",
  facet_col = "Facet"
)
```

The expected wheat pattern is:

``` text
Facet 1:  W001-W010
Facet 2:  W011-W020
...
Facet 10: W091-W100
```

The `Facet_Assignment_OK` column should be `TRUE` for all treatment genotypes.

### Global Design Audit

Use:

``` r
design_audit <- audit_rotation_design(wheat_treatment)
design_audit$design_ok
design_audit$design_summary
```

The model functions now attach design diagnostics when possible. Inspect `result$design_audit` and `result$diagnostics` before ranking lentils or interpreting pair effects.

## Independent Trait Models

The generic trait model is implemented in `model_traits()`.

### SpATS Trait Model

For lentil or wheat phenotyping, the spatial model is conceptually:

``` r
trait ~ design terms + genotype random effect + PSANOVA(Row, Col)
```

Mathematically, within each environment:

$$
Y_i = \mu_e + d_{e,r(i)} + G_{e,j(i)} + f_e(R_i, C_i) + \varepsilon_i
$$

where:

-   `y_i` is the observed value of the trait for plot `i`.
-   `mu_e` is the environment-specific intercept or trial mean.
-   `d_{e,r(i)}` is the fixed design effect for the replicate, block, or other design term assigned to plot `i`.
-   `g_{e,j(i)}` is the genotype effect for genotype `j` in environment `e`; in the SpATS and lme4 implementations this is estimated as a random effect/BLUP.
-   `f_e(row_i, col_i)` is the smooth two-dimensional spatial surface fitted from plot row and column, implemented with `PSANOVA(Row, Col)` in SpATS.
-   `epsilon_i` is the residual plot-level error.

The genotype prediction reported in the BLUP table is:

$$
\mathrm{Predicted}_{e,j} = \hat{\mu}_e + \hat{G}_{e,j}
$$

after accounting for the design and spatial terms in the fitted model.

In package terms:

``` r
model_traits(
  data = lentil_treatment,
  method = "SpATS",
  trait_cols = c("YLD"),
  gen_col = "Lentil",
  env_col = "ENV",
  rep_col = "Rep_gen",
  spatial_cols = c("Row", "Col")
)
```

For wheat:

``` r
model_traits(
  data = wheat_treatment,
  method = "SpATS",
  trait_cols = c("Y.ADJ"),
  gen_col = "Wheat",
  env_col = "ENV",
  rep_col = "Rep_gen",
  spatial_cols = c("Row", "Col")
)
```

This independent wheat genotype model uses `Rep_gen` because the target is wheat genotype performance and each wheat genotype is represented across about 20 plots within an environment. Rotation models switch to `Rep_combo` because their experimental unit is the observed lentil-wheat combination, replicated twice.

### Insight

This model gives:

-   Genotype BLUPs or adjusted predictions.
-   Generalized heritability estimates from SpATS.
-   Spatial trend surfaces.
-   Trait-by-environment summaries.

These outputs answer:

``` text
Which genotypes perform well for a trait after correcting for field spatial trends?
```

## Predecessor Effect Model

The preferred function is:

``` r
model_predecessor_effect()
```

### Primary ENV × Facet Model

The primary recommended model is fit within each ENV × Facet network by default:

``` r
model_predecessor_effect(..., fit_scope = "env_facet")
```

Conceptually, the model is:

``` r
wheat_trait ~ previous_lentil + current_wheat + design terms + spatial trend
```

Depending on the implementation, current wheat genotype may be fitted as a fixed effect or random effect. The important point is that wheat genotype is included so that lentil predecessor effects are not confounded with the wheat genotypes assigned to the same facet.

In the SpATS implementation:

``` r
Y.ADJ ~ Lentil + Rep_combo + PSANOVA(Row, Col) + random(Wheat)
```

Here `Rep_combo` is the relevant replicate term because each lentil-wheat `Combo` is replicated twice within an environment. `Rep_gen` is not used in this model; it belongs to independent genotype performance analyses where the same genotype has many replicate plots.

Mathematically, within each ENV × Facet network:

$$
Y_i = \mu_{e,b} + L_{e,b,l(i)} + d_{e,b,q(i)} + W_{e,b,w(i)} + f_{e,b}(R_i, C_i) + \varepsilon_i
$$

where:

-   `y_i` is the observed Year 2 wheat trait value for plot `i`.
-   `mu_{e,b}` is the ENV × Facet mean.
-   `L_{e,b,l(i)}` is the effect of the previous lentil genotype `l` grown in that plot in Year 1. This is the predecessor effect of interest.
-   `d_{e,b,q(i)}` is the fixed combo-replicate design effect, usually `Rep_combo`, for plot `i`.
-   `W_{e,b,w(i)}` is the correction for the current wheat genotype `w` in Year 2.
-   `f_{e,b}(row_i, col_i)` is the spatial field trend for the ENV × Facet network.
-   `epsilon_i` is residual plot-level error.

The model-corrected predecessor mean for lentil genotype `l` is:

$$
\mathrm{Corrected\_Mean}_{e,b,l} = \hat{\mu}_{e,b} + \hat{L}_{e,b,l}
$$

For a lentil genotype assigned to wheat facet `b`, the facet baseline is:

$$
\mathrm{Baseline\_Mean}_{e,b}
=
\frac{
\sum_{l \in (e,b)} \omega_{e,b,l}\mathrm{Corrected\_Mean}_{e,b,l}
}{
\sum_{l \in (e,b)} \omega_{e,b,l}
}
$$

where `\omega_{e,b,l} = 1 / SE_{e,b,l}^{2}` when valid SEs are available. If SEs are unavailable, the implementation falls back to the unweighted mean.

The reported legacy quantities are:

$$
\mathrm{Legacy\_Value}_{e,b,l}
= \mathrm{Corrected\_Mean}_{e,b,l} - \mathrm{Baseline\_Mean}_{e,b}
$$

$$
\mathrm{Legacy\_Value\_Global}_{e,b,l}
= \mathrm{Corrected\_Mean}_{e,b,l}
-
\frac{
\sum_{b'}\sum_{l' \in (e,b')} \omega_{e,b',l'}\mathrm{Corrected\_Mean}_{e,b',l'}
}{
\sum_{b'}\sum_{l' \in (e,b')} \omega_{e,b',l'}
}
$$

$$
\mathrm{Total\_Correction}_{e,b,l}
= \mathrm{Corrected\_Mean}_{e,b,l} - \mathrm{Raw\_Mean}_{e,b,l}
$$

where:

-   `b(l)` is the wheat-partner facet associated with lentil genotype `l`.
-   `Raw_Mean_{e,b,l}` is the uncorrected mean wheat performance observed after lentil genotype `l` within facet `b`.
-   `Total_Correction_{e,b,l}` measures how much the model shifted the raw lentil-associated wheat mean after accounting for wheat partners, design terms, and spatial position.

`Legacy_Value_Global` is a secondary cross-facet ranking obtained by averaging facet-stratified corrected means across the environment. This comparison is only approximately valid because facets are not genetically connected by a complete 100 × 100 design; they are linked mainly through the shared environment and spatial correction. The ENV × Facet `Legacy_Value` remains the preferred interpretation.

`Network_Correction` is retained as a backward-compatible column name in the package output, but `Total_Correction` is the clearer interpretation. This value reflects all model adjustments collectively, including wheat partner correction, spatial correction, and design terms, not only the genotype network component.

Wheat genotype is included as a correction for the incomplete pairing network. If a lentil was paired with low-yielding wheat genotypes, the wheat term helps correct the lentil estimate upward. If it was paired with high-yielding wheat genotypes, the model corrects downward.

After the model estimates corrected lentil means, the preferred legacy value is calculated relative to the facet baseline:

$$
\mathrm{Legacy\_Value}
= \mathrm{Corrected\_Mean}
- \operatorname{mean}_{\mathrm{ENV}\times\mathrm{Facet}}\left(\mathrm{Corrected\_Mean}\right)
$$

The raw-vs-corrected output also reports:

$$
\mathrm{Total\_Correction}
= \mathrm{Corrected\_Mean} - \mathrm{Raw\_Mean}
$$

### Example

``` r
avg_predecessor_yadj <- model_predecessor_effect(
  data = wheat_treatment,
  trait = "Y.ADJ",
  env_col = "ENV",
  prev_gen_col = "Lentil",
  curr_gen_col = "Wheat",
  spatial_cols = c("Row", "Col"),
  baseline_col = "Facet",
  fixed_effect_cols = c("Rep_combo"),
  type_col = "Type",
  include_checks = FALSE,
  method = "SpATS"
)
```

### Main Outputs

The output table includes:

-   `Previous_Genotype`: lentil genotype.
-   `Corrected_Mean`: model-corrected wheat performance following that lentil.
-   `Raw_Mean`: uncorrected mean of wheat performance following that lentil.
-   `Baseline_Group`: wheat facet used as the comparison baseline.
-   `Baseline_Mean`: corrected mean within the ENV × Facet baseline.
-   `Legacy_Value`: corrected deviation from the ENV × Facet baseline.
-   `Legacy_Value_Global`: corrected deviation from the full ENV mean.
-   `Total_Correction`: difference between corrected and raw means after all model adjustments.
-   `Network_Correction`: backward-compatible alias for `Total_Correction`.

### Insight

This formula answers:

``` text
Which lentil genotypes are associated with higher or lower following-wheat performance relative to the other lentils tested within the same ENV × Facet network?
```

The main plot shows `Legacy_Value`, centered on the ENV × Facet mean.

The correction plot shows:

-   Hollow point: raw mean.
-   Solid point: corrected mean.
-   Dashed vertical line: facet baseline.

This helps diagnose whether the model is correcting a lentil upward or downward because of its wheat partners and spatial position.

### GWAS-ready legacy values

For lentil GWAS, the useful output is one `Legacy_Value` per lentil genotype, environment, and wheat trait. `model_predecessor_effect()` returns these values in:

``` r
avg_predecessor_yadj$legacy_values
```

The key column is:

``` text
Legacy_Value
```

This is the model-corrected wheat response after a lentil genotype, expressed as a deviation from the relevant ENV × Facet baseline.

Mathematically:

$$
\mathrm{Legacy\_Value}_{e,l,t}
= \mathrm{Corrected\_Mean}_{e,l,t} - \mathrm{Baseline\_Mean}_{e,b(l),t}
$$

where:

-   `e` is the wheat environment.
-   `l` is the previous lentil genotype.
-   `t` is the wheat trait being used as the legacy target, such as `Y.ADJ` or `PRO`.
-   `b(l)` is the wheat-partner facet associated with lentil genotype `l`.
-   `Legacy_Value_{e,l,t}` is the facet-corrected wheat legacy value estimated from the predecessor model.

This value can be interpreted as:

``` text
How much better or worse did wheat perform after this lentil genotype than expected for its wheat-partner facet?
```

Positive values indicate lentil genotypes associated with improved following-wheat performance. These values can be used as Year 2 legacy phenotypes for association with lentil genomic markers. The biological target is a lentil genomic signature associated with improving the next crop, rather than improving the lentil crop itself.

The ranked legacy-value plot is available as:

``` r
avg_predecessor_yadj$ranked_plot
```

This plot puts all lentil genotypes into one ranked list per environment, while coloring bars by `Facet` to preserve the structured incomplete design context.

The SpATS output includes prediction standard errors for the adjusted lentil predecessor means. These are kept in the output table as `SE`, but they are not drawn by default in the ranked plots because `Legacy_Value` is a deviation from a facet baseline and the SE bars can visually overwhelm the ranking.

Boxplots are better reserved for plot-level raw or residualized observations. The `legacy_values` table has one adjusted legacy value per lentil genotype × environment × trait, so a boxplot would imply a distribution that is not present in that output.

For GWAS, `Legacy_Value` is an estimated phenotype and its uncertainty should be propagated when the GWAS software supports it. Use inverse-variance weights such as `1 / SE^2`, or the available reliability column when appropriate, in weighted regression frameworks such as GEMMA weighted linear models or heteroscedastic Bayesian models. If the GWAS workflow cannot use weights, a conservative alternative is to subset to higher-confidence phenotypes, for example `Confidence_Class == "high"`, before association testing.

### Secondary Global Sensitivity Model

A global environment-level sensitivity model can be requested with:

``` r
model_predecessor_effect(..., fit_scope = "env_global")
```

This fits across the full environment and then centers corrected lentil means by facet. It is useful as a sensitivity analysis, but the primary interpretation remains the ENV × Facet model because each facet is the locally complete 10 × 10 network.

### Uncertainty and Reliability

The output keeps `SE` and `Corrected_Mean_SE` for the adjusted lentil predecessor mean before facet-centering. `Legacy_Value_SE` is retained as `NA` unless a valid centered-contrast SE is available. The output also includes `Confidence_Class`, which summarizes SE into `high`, `moderate`, `low`, or `unknown` categories for diagnostics.

## Pair Compatibility Model

The preferred function is:

``` r
model_pair_compatibility()
```

This function asks whether specific observed lentil-wheat pairs have extra compatibility beyond additive genotype effects.

The conceptual formula is:

``` r
wheat_trait ~ Lentil + Wheat + Rep_combo + spatial terms + random(Combo) + random(Block)
```

The `Combo` random effect is the pair-specific compatibility value.

Mathematically, within each ENV × Facet network:

$$
Y_i = \mu_{e,b} + L_{e,b,l(i)} + W_{e,b,w(i)} + d_{e,b,q(i)}
      + f_{e,b}(R_i, C_i) + C_{e,b,c(i)} + B_{e,b,k(i)} + \varepsilon_i
$$

where:

-   `y_i` is the observed Year 2 wheat trait value for plot `i`.
-   `mu_{e,b}` is the mean for environment `e` and facet `b`.
-   `L_{e,b,l(i)}` is the additive effect of the previous lentil genotype `l`.
-   `W_{e,b,w(i)}` is the additive effect of the current wheat genotype `w`.
-   `d_{e,b,q(i)}` is the fixed combo-replicate effect, usually `Rep_combo`.
-   `f_{e,b}(row_i, col_i)` is the spatial correction used by the lme4 implementation, represented by polynomial row and column terms plus their interaction.
-   `C_{e,b,c(i)}` is the random effect for the observed lentil-wheat pair `c`, stored in `Combo`.
-   `B_{e,b,k(i)}` is the optional random block or design effect for block `k`.
-   `epsilon_i` is residual plot-level error.

The additive expectation for pair `c = (l,w)` is:

$$
\mathrm{Expected\_Additive\_Mean}_{e,b,c}
= \frac{1}{n_c}\sum_{i:c(i)=c}
\left(
\hat{\mu}_{e,b} + \hat{L}_{e,b,l(i)} + \hat{W}_{e,b,w(i)}
+ \hat{d}_{e,b,q(i)} + \hat{f}_{e,b}(R_i, C_i)
\right)
$$

over the plots where that pair was observed. The corrected pair mean is:

$$
\mathrm{Corrected\_Pair\_Mean}_{e,b,c}
= \frac{1}{n_c}\sum_{i:c(i)=c}
\left(
\mathrm{Expected\_Additive}_i + \hat{C}_{e,b,c(i)} + \hat{B}_{e,b,k(i)}
\right)
$$

and the compatibility value is:

$$
\mathrm{Compatibility\_Value}_{e,b,c} = \hat{C}_{e,b,c}
$$

$$
\mathrm{Compatibility\_Pct}_{e,b,c}
= 100 \times
\frac{\mathrm{Compatibility\_Value}_{e,b,c}}
     {\mathrm{Expected\_Additive\_Mean}_{e,b,c}}
$$

Positive `Compatibility_Value` means the observed pair performed better than expected from the additive lentil and wheat effects in that ENV × Facet network.

### Local 10 × 10 Factorial Interpretation

The current implementation fits this model within each ENV × Facet network. This is important because the full 100 × 100 factorial is not observed. However, within each facet, the design is locally complete: 10 lentil genotypes are evaluated against 10 wheat genotypes, giving 100 observed combinations. Therefore, `Compatibility_Value` should be interpreted as a pair-specific deviation from additive lentil and wheat effects within that local 10 × 10 network.

### Compatibility_Value vs Compatibility_Pct

`Compatibility_Value` is the primary output. `Compatibility_Pct` is optional and can be unstable when `Expected_Additive_Mean` is near zero or when the response is centered, ordinal, residualized, or can be negative. Use:

``` r
model_pair_compatibility(
  ...,
  compute_compatibility_pct = TRUE,
  pct_min_denominator = 1e-6
)
```

If the denominator is too small, `Compatibility_Pct` is set to `NA` and the function warns that `Compatibility_Value` should be used as the primary interpretation.

`GxG_Rotational_Pct` is computed through the same code path as `Compatibility_Pct` in `model_pair_compatibility()`. The `compute_pct` and `pct_min_denominator` arguments in `model_gxg_rotational_value()` provide the same guard against near-zero additive expectations. Use `GxG_Rotational_Value` as the primary quantity when the additive expectation is small, negative, centered, ordinal, or otherwise difficult to interpret as a denominator.

### Corrected Pair Means

The pair output distinguishes:

-   `Conditional_Predicted_Mean`: fitted mean including observed random design effects such as block.
-   `Marginal_Corrected_Pair_Mean`: fitted additive expectation plus the pair-specific `Combo` BLUP, excluding nuisance random block effects.
-   `Compatibility_Value`: the pair-specific `Combo` BLUP and preferred ranking variable.

`Corrected_Pair_Mean` is retained as an alias for `Marginal_Corrected_Pair_Mean` for backward compatibility.

### Example

``` r
pair_compatibility_yadj <- model_pair_compatibility(
  data = wheat_treatment,
  trait = "Y.ADJ",
  env_col = "ENV",
  prev_gen_col = "Lentil",
  curr_gen_col = "Wheat",
  combo_col = "Combo",
  spatial_cols = c("Row", "Col"),
  baseline_col = "Facet",
  fixed_effect_cols = c("Rep_combo"),
  random_effect_cols = c("Block"),
  type_col = "Type",
  include_checks = FALSE
)
```

### Main Outputs

The output table includes:

-   `Combo`: observed lentil-wheat pair.
-   `Previous_Genotype`: lentil genotype.
-   `Current_Genotype`: wheat genotype.
-   `Baseline_Group`: wheat facet.
-   `Raw_Mean`: raw mean for the pair.
-   `Expected_Additive_Mean`: predicted value from lentil, wheat, and spatial effects without the pair effect.
-   `Conditional_Predicted_Mean`: fitted mean including observed design/block random effects.
-   `Marginal_Corrected_Pair_Mean`: corrected pair mean excluding nuisance random block effects.
-   `Corrected_Pair_Mean`: backward-compatible alias for `Marginal_Corrected_Pair_Mean`.
-   `Compatibility_Value`: pair-specific deviation from the additive expectation.
-   `Compatibility_Pct`: optional compatibility as a percent of the additive expectation.

### Insight

This formula answers:

``` text
Which observed lentil-wheat pairs are unusually good or bad after accounting for the average lentil effect, average wheat effect, and spatial field position?
```

Interpretation should stay inside the observed network. A positive pair effect means the observed pair performed better than expected for that specific lentil and wheat combination. It does not directly predict unobserved pairings.

The default pair-compatibility plot is now a heatmap rather than a long bar plot. Each panel is one observed ENV × Facet network, the x-axis is the current wheat genotype, the y-axis is the previous lentil genotype, and the tile color is `Compatibility_Value`.

``` r
pair_compatibility_pro$heatmap
pair_compatibility_pro$ranked_plot
pair_compatibility_pro$diagnostics
```

The heatmap avoids saturating the y-axis with long lentil-wheat combination names. The ranked plot is still available for the strongest positive and negative pair effects, but it is intended as a summary view, not as the main map of the design.

If a panel shows all pair effects as zero, that does not mean the underlying observations are missing. It means the mixed model estimated the `Combo` variance component as zero for that ENV × Facet. In that situation, the fitted model says that the additive lentil effect, additive wheat effect, spatial terms, and design random effects explain the available signal, and there is no detectable pair-specific deviation left for that facet. The pair BLUPs are therefore all shrunk to zero.

By default, pair-compatibility plots now drop ENV × Facet groups with no detectable pair-specific signal. These groups remain in `pair_compatibility_pro$diagnostics`, but they are not shown in `heatmap` or `ranked_plot`. To include them for auditing:

``` r
plot_pair_compatibility_heatmap(pair_compatibility_pro, drop_no_signal = FALSE)
plot_pair_compatibility_ranked(pair_compatibility_pro, drop_no_signal = FALSE)
```

### Diagnostics and Singular Fits

The diagnostics table should be checked before interpreting pair effects. Useful columns include:

-   `N_Plots`
-   `N_Previous_Genotypes`
-   `N_Current_Genotypes`
-   `N_Combos`
-   `Expected_Combos`
-   `Missing_Combos`
-   `Local_Factorial_Complete`
-   `Expected_10x10`
-   `Min_Reps_Per_Combo`
-   `N_Single_Rep_Combos`
-   `Status`
-   `Model_Singular`
-   `Combo_Variance`
-   `Residual_Variance`

For this design, a reliable ENV × Facet group should usually have 10 lentils, 10 wheats, 100 observed combos, and at least two plots per combo. If local completeness is broken, the diagnostic `Status` flags `incomplete_local_factorial` or `low_replication`. If a panel appears empty in an old bar plot, that is usually a visualization artifact from crossing all environment and facet labels. The heatmap uses only observed ENV × Facet panels and should not generate fake tiles for combinations outside the design.

## Formula Layer 4: Wheat Rotational Partner Value

The wheat subsample layer allows a third biological question to be addressed:

> Which wheat genotypes are better rotational partners after lentil?

This question differs from the lentil legacy question. The lentil legacy model asks which lentil genotypes leave better conditions for the following wheat crop. The wheat rotational partner model asks which wheat genotypes are better able to capture, convert, or partition those legacy effects into useful agronomic and nutrient outcomes.

A good wheat rotational partner is not necessarily the highest yielding wheat genotype overall. It may be a genotype that performs well after lentil by producing high grain yield, capturing more nitrogen, partitioning nitrogen efficiently into grain, maintaining protein, or returning useful straw biomass and straw nitrogen to the system.

Relevant wheat subsample traits include grain biomass, straw biomass, total aboveground biomass, seed N, straw N, total aboveground N, grain N content, straw N content, N harvest index, seed C:N ratio, straw C:N ratio, and protein yield.

Within each ENV × Facet network, the conceptual model is:

``` r
wheat_subsample_trait ~ Wheat + Lentil + design_terms + spatial_terms
```

The corrected wheat mean is calculated for each wheat genotype after accounting for the previous lentil genotype, design terms, and spatial field structure.

The wheat rotational value is calculated relative to the ENV × Facet wheat baseline:

$$
\mathrm{Wheat\_Rotational\_Value}_{e,b,w,t}
=
\mathrm{Corrected\_Wheat\_Mean}_{e,b,w,t}
-
\frac{
\sum_{w' \in (e,b)} \omega_{e,b,w',t}\mathrm{Corrected\_Wheat\_Mean}_{e,b,w',t}
}{
\sum_{w' \in (e,b)} \omega_{e,b,w',t}
}
$$

where `e` is the environment, `b` is the facet, `w` is the wheat genotype, `t` is the wheat subsample trait, and `\omega = 1 / SE^2` when valid SEs are available. If SEs are unavailable, the implementation falls back to the unweighted mean.

The secondary global wheat rotational value is:

$$
\mathrm{Wheat\_Rotational\_Value\_Global}_{e,b,w,t}
=
\mathrm{Corrected\_Wheat\_Mean}_{e,b,w,t}
-
\frac{
\sum_{b'}\sum_{w' \in (e,b')} \omega_{e,b',w',t}\mathrm{Corrected\_Wheat\_Mean}_{e,b',w',t}
}{
\sum_{b'}\sum_{w' \in (e,b')} \omega_{e,b',w',t}
}
$$

As with `Legacy_Value_Global`, this is a secondary cross-facet ranking. It should be interpreted cautiously because wheat genotypes are not evaluated in a complete 100 × 100 lentil × wheat design.

A positive value means that the wheat genotype performed better than the average wheat genotype in the same ENV × Facet network for that trait, after correcting for the lentil predecessor and field structure.

This output answers:

``` text
Which wheat genotypes are better following crops after lentil within their local 10 × 10 rotation network?
```

### Example

``` r
wheat_subsample_derived <- derive_wheat_rotation_traits(
  wheat_subsample,
  grain_biomass_col = "Grain_Biomass_g",
  straw_biomass_col = "Straw_Biomass_g",
  seed_n_pct_col = "Seed_N_pct",
  straw_n_pct_col = "Straw_N_pct",
  seed_c_pct_col = "Seed_C_pct",
  straw_c_pct_col = "Straw_C_pct"
)

wheat_rot_n <- model_wheat_rotational_value(
  data = wheat_subsample_derived,
  trait = "Total_Aboveground_N_g",
  env_col = "ENV",
  facet_col = "Facet",
  prev_gen_col = "Lentil",
  curr_gen_col = "Wheat",
  spatial_cols = c("Row", "Col"),
  fixed_effect_cols = c("Rep_combo"),
  type_col = "Type",
  include_checks = FALSE
)

wheat_rot_n$rotational_values
wheat_rot_n$ranked_plot
wheat_rot_n$raw_vs_corrected_plot
```

The output table includes:

-   `Current_Genotype`: wheat genotype.
-   `Corrected_Wheat_Mean`: model-corrected wheat performance after accounting for previous lentil genotype, design terms, and spatial structure.
-   `Raw_Mean`: uncorrected wheat mean.
-   `Baseline_Mean`: corrected wheat mean within ENV × Facet.
-   `Wheat_Rotational_Value`: corrected deviation from the ENV × Facet wheat baseline.
-   `Wheat_Rotational_Value_Global`: corrected deviation from the ENV-wide wheat baseline.
-   `Total_Correction`: difference between corrected and raw wheat means after all model adjustments.
-   `Network_Correction`: backward-compatible alias for `Total_Correction`.

## Formula Layer 5: G×G Rotational Value

The wheat subsample layer can also be used to estimate a rotation-specific genotype-by-genotype interaction value.

This asks:

> Which observed lentil-wheat genotype combinations perform better than expected for rotation-relevant traits such as yield, grain N, total N uptake, straw N, biomass, N harvest index, or residue quality?

Within each ENV × Facet network, the conceptual model is:

``` r
wheat_subsample_trait ~ Lentil + Wheat + Rep_combo + spatial_terms + random(Combo) + random(Block)
```

The additive expectation for each pair is based on the average lentil legacy effect and the average wheat rotational partner effect. The pair-specific random effect for `Combo` represents the G×G rotational value:

$$
\mathrm{GxG\_Rotational\_Value}_{e,b,l,w,t}
=
\hat{C}_{e,b,l,w,t}
$$

where `e` is the environment, `b` is the facet, `l` is the previous lentil genotype, `w` is the current wheat genotype, and `t` is the rotation-relevant wheat trait.

A positive G×G rotational value means that the observed lentil-wheat pair performed better than expected from the additive lentil legacy effect and the additive wheat rotational partner effect within the same local 10 × 10 network.

This output answers:

``` text
Which specific lentil-wheat combinations are unusually good or bad for rotation-relevant outcomes?
```

The G×G rotational value should only be interpreted within observed ENV × Facet networks. It should not be extrapolated to unobserved lentil-wheat combinations.

### Example

``` r
gxg_rot_n <- model_gxg_rotational_value(
  data = wheat_subsample_derived,
  trait = "Total_Aboveground_N_g",
  env_col = "ENV",
  facet_col = "Facet",
  prev_gen_col = "Lentil",
  curr_gen_col = "Wheat",
  combo_col = "Combo",
  spatial_cols = c("Row", "Col"),
  fixed_effect_cols = c("Rep_combo"),
  random_effect_cols = c("Block"),
  type_col = "Type",
  include_checks = FALSE
)

gxg_rot_n$gxg_values
gxg_rot_n$heatmap
gxg_rot_n$ranked_plot
```

The output table includes:

-   `Combo`: observed lentil-wheat pair.
-   `Previous_Genotype`: previous lentil genotype.
-   `Current_Genotype`: current wheat genotype.
-   `Trait`: rotation-relevant wheat trait.
-   `Raw_Mean`: raw mean for the observed pair.
-   `Expected_Additive_Mean`: expected value from additive lentil and wheat effects.
-   `Marginal_Corrected_Pair_Mean`: corrected pair mean excluding nuisance random block effects.
-   `GxG_Rotational_Value`: pair-specific deviation from additive lentil and wheat effects.
-   `GxG_Rotational_Pct`: optional percentage deviation from the additive expectation.
-   `N_Plots`: number of plots contributing to the pair estimate.
-   `Status`: model or diagnostic status for the ENV × Facet group.

`GxG_Rotational_Pct` uses the same denominator guard as `Compatibility_Pct`. If `Expected_Additive_Mean` is near zero, the percentage value is set to `NA`; the signed `GxG_Rotational_Value` should remain the primary interpretation.

For rotation-specific subsample traits, `GxG_Rotational_Value` is the clearer synonym for the existing `Compatibility_Value`. Both refer to the pair-specific deviation from additive lentil and wheat effects within an observed local ENV × Facet network.

### Rotation Indices

Wheat rotational partner values can be combined across traits using a transparent user-defined index:

``` r
wheat_rotation_index <- build_wheat_rotation_index(
  rotational_values = dplyr::bind_rows(
    wheat_rot_yield$rotational_values,
    wheat_rot_total_n$rotational_values,
    wheat_rot_grain_n$rotational_values,
    wheat_rot_nhi$rotational_values
  ),
  trait_weights = c(
    Y.ADJ = 1,
    Total_Aboveground_N_g = 1,
    Grain_N_Content_g = 1,
    N_Harvest_Index = 1
  )
)
```

This index is a decision-support score, not an intrinsic biological property. Different breeding goals may require different trait weights.

`build_wheat_rotation_index()` and `build_gxg_rotation_index()` standardize values internally within ENV × Facet × Trait by default through their `standardize_within` argument, so high-variance traits do not dominate the score purely because of measurement scale. Trait direction is not automatically aligned. Users should set positive or negative trait weights so that positive index values consistently represent desirable rotational performance for the specific breeding or agronomic goal.

## Correlation and Driver Analysis

Once lentil BLUPs and wheat legacy values are estimated, the package can combine them into a long table and use:

``` r
plot_legacy_correlations()
```

This creates a correlation heatmap across traits.

The correlation input can now be built from modular pieces:

``` r
lentil_full_trial_predictors <- prepare_lentil_blup_predictors(
  blups = lentil_lme4_models$blups,
  trait_prefix = "FullTrial_"
)

lentil_subsample_predictors <- prepare_lentil_subsample_predictors(
  data = lentil_subsample,
  trait_cols = lentil_subsample_traits,
  trait_prefix = "Subsample_"
)

lentil_microbiome_predictors <- prepare_lentil_microbiome_predictors(
  data = lentil_microbiome,
  trait_cols = lentil_microbiome_traits,
  trait_prefix = "Microbiome_"
)

wheat_legacy_targets <- prepare_wheat_legacy_targets(
  YADJ = avg_predecessor_yadj,
  PRO = avg_predecessor_pro
)

legacy_correlation_input <- build_legacy_correlation_input(
  predictors = list(
    lentil_full_trial_predictors,
    lentil_subsample_predictors,
    lentil_microbiome_predictors
  ),
  targets = wheat_legacy_targets
)
```

This lets the same predictor table be reused for different wheat legacy targets:

``` r
plot_legacy_correlations(
  data = legacy_correlation_input,
  target_trait = "Wheat_Legacy_PRO"
)
```

The same modular correlation framework can be extended to wheat rotational partner values and G×G rotational values. For example, wheat phenotypic traits, wheat subsample traits, or wheat genomic predictors can be correlated with `Wheat_Rotational_Value` to ask which wheat characteristics are associated with better performance after lentil. Similarly, pair-level predictors can be explored against `GxG_Rotational_Value`, although these analyses should remain within observed ENV × Facet networks.

Mathematically, for each environment, lentil predictor trait, and wheat legacy target:

$$
X_{e,l,p} = \text{lentil predictor value for genotype } l
$$

$$
Z_{e,l,t} = \text{wheat legacy target for genotype } l
$$

The environment-specific correlation is:

$$
r_{e,p,t}
=
\frac{
\sum_l \left(X_{e,l,p} - \bar{X}_{e,p}\right)
       \left(Z_{e,l,t} - \bar{Z}_{e,t}\right)
}{
\sqrt{
\sum_l \left(X_{e,l,p} - \bar{X}_{e,p}\right)^2
\sum_l \left(Z_{e,l,t} - \bar{Z}_{e,t}\right)^2
}
}
$$

The reported driver summaries use:

$$
R^2_{e,p,t} = r_{e,p,t}^{2}
$$

where:

-   `X_{e,l,p}` can come from lentil full-trial BLUPs, lentil subsample summaries, or lentil microbiome summaries.
-   `Z_{e,l,t}` is usually `Legacy_Value` from `model_predecessor_effect()`, renamed as a wheat legacy target such as `Wheat_Legacy_YADJ`.
-   `l` indexes lentil genotypes with both predictor and target values.
-   `r_{e,p,t}` is the Pearson or Spearman correlation between a Year 1 lentil predictor and a Year 2 wheat legacy target.
-   `R2_{e,p,t}` is the proportion of target variation explained by that single predictor in the environment-specific exploratory correlation.
-   p-values and FDR values summarize statistical evidence for the correlation, but they do not establish causality.

`plot_legacy_driver_correlations()` applies Benjamini-Hochberg FDR correction within each `Environment × Target_Trait` group. This controls the expected false discovery rate across the predictor traits tested for a given target in a given environment, but it does not control a single global FDR across every trait, target, and environment in the full analysis.

### Exploratory Nature

These correlations are exploratory because both lentil predictors and wheat legacy values may be estimated quantities. Correlation results should be interpreted as hypothesis-generating associations, not causal evidence. When available, SE or reliability values should be used to assess whether strong correlations are driven by uncertain legacy estimates. P-values and FDR values summarize statistical evidence for correlation, but they do not establish causality.

### Minimum Sample Size

When many subsample traits are included, the focused driver plot is easier to read than the full heatmap. The default minimum number of complete genotype pairs is now 30:

``` r
legacy_driver_pro <- plot_legacy_driver_correlations(
  data = legacy_correlation_input,
  target_trait = "Wheat_Legacy_PRO",
  top_n = 20,
  min_pairs = 30
)

legacy_driver_pro$correlations
legacy_driver_pro$plot
```

The driver table reports the environment-specific correlation, R2, p-value, Benjamini-Hochberg FDR, and number of genotype pairs used for each lentil predictor trait.

### Weighted/Spearman Options

Use `method = "spearman"` for rank-based exploratory correlations when linearity or outliers are a concern:

``` r
plot_legacy_driver_correlations(
  data = legacy_correlation_input,
  target_trait = "Wheat_Legacy_PRO",
  method = "spearman",
  min_pairs = 30
)
```

Weighted correlations are optional because predictor and legacy SEs are not always available on the same scale. `prepare_wheat_legacy_targets()` carries `SE`, `Reliability`, and a derived `Weight` column when available. To weight correlations by the wheat legacy target uncertainty, use:

``` r
plot_legacy_driver_correlations(
  data = legacy_correlation_input,
  target_trait = "Wheat_Legacy_PRO",
  method = "pearson",
  weight_col = "Weight",
  min_pairs = 30
)
```

Weighted correlations report `P_Value = NA` because the weighted correlation is treated as an exploratory effect-size screen rather than a standard unweighted correlation test.

### Microbiome Predictor Cautions

Microbiome predictors should be interpreted carefully because amplicon-derived relative abundances are compositional. Taxon-level predictors such as Rhizobiaceae or Bradyrhizobiaceae should preferably be analyzed with centred log-ratio (CLR) transformed abundances rather than raw relative abundances in regression or correlation contexts. Isometric log-ratio (ILR) balances are preferable for formal inference when a defensible balance structure is available. Aitchison PCoA axes are already based on CLR/Aitchison geometry and do not require an additional taxon-level log-ratio transformation, but their orientation should still be checked before interpreting the sign of correlations with wheat legacy values.

The microbiome aggregation step should report the number of microbiome samples per lentil × environment and should be checked for imbalance across facets. Raw relative abundance should not be interpreted as absolute abundance.

### Insight

This analysis asks:

``` text
Which Year 1 lentil traits are associated with positive or negative Year 2 wheat legacy effects?
```

Examples of possible biological interpretations include:

-   Lentil yield or biomass could correlate with wheat legacy if residue quantity matters.
-   Lentil maturity or phenology could correlate with wheat legacy if timing affects soil water or residue decomposition.
-   Lentil protein or nutrient traits could correlate with wheat legacy if nitrogen or nutrient cycling contributes to the effect.
-   Microbiome diversity, ordination axes, or symbiont-associated abundance could correlate with wheat legacy if microbial community structure contributes to residue decomposition, nutrient cycling, or soil biological carryover.

These correlations are exploratory. They can identify hypotheses, but they do not prove mechanism without follow-up modeling or experimental validation.

## Recommended Interpretation Workflow

1.  Inspect trait distributions and field structure.
2.  Run `audit_rotation_design()` and confirm local 10 × 10 completeness before modeling.
3.  Fit independent lentil and wheat trait models.
4.  Examine spatial trends to confirm field correction is needed.
5.  Estimate average lentil predecessor effects with `model_predecessor_effect(fit_scope = "env_facet")`.
6.  Use `Legacy_Value` for the ENV × Facet interpretation.
7.  Use `Legacy_Value_Global` only as a secondary sensitivity comparison.
8.  Use the raw-vs-corrected plot to show how much the model corrected for wheat partners and spatial structure.
9.  Estimate pair compatibility with `model_pair_compatibility()`.
10. Interpret pair compatibility only within observed ENV × Facet 10 × 10 networks.
11. Derive wheat subsample rotation traits with `derive_wheat_rotation_traits()`.
12. Estimate wheat rotational partner value with `model_wheat_rotational_value()`.
13. Estimate rotation-specific pair effects with `model_gxg_rotational_value()`.
14. Use `build_wheat_rotation_index()` or `build_gxg_rotation_index()` only as transparent, user-weighted decision-support summaries.
15. Correlate lentil BLUPs with wheat legacy values to generate biological hypotheses.

## Practical Notes

-   Check plots should usually be excluded from rotation-effect models unless the question specifically concerns checks.
-   `Facet` should be used as the local comparison baseline for predecessor legacy values.
-   `Facet` should not be naively added as a fixed effect if it is redundant with genotype structure.
-   Pair compatibility results are conditional on the observed local 10 × 10 design.
-   `Wheat_Rotational_Value` is the within-facet ability of a wheat genotype to perform after lentil, correcting for the previous lentil genotype.
-   `GxG_Rotational_Value` is the pair-specific deviation from additive lentil and wheat effects within an observed local 10 × 10 network.
-   Singular fits can occur in pair models when variance components are small or a facet has too little information. This can happen even inside local 10 × 10 networks and should be interpreted carefully.
-   The testing script `test_rotation_formulas_and_plots.R` is the best place to validate formulas and plots before promoting analyses into package vignettes or manuscript-ready scripts.
