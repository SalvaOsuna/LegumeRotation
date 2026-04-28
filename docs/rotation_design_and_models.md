# Legume Rotation Experiment: Design, Modeling Challenges, and Analytical Targets

## Purpose of the Package

`LegumeRotation` supports analysis of a two-year lentil-wheat rotation experiment. The package has two connected goals:

1.  Characterize phenotypic performance of lentil and wheat as independent field trials.
2.  Test whether the lentil genotype grown in a plot in Year 1 influences the agronomic performance of the wheat genotype grown in the same plot in Year 2.

The main rotation question is a predecessor or legacy-effect question:

> Does lentil genotype X leave a measurable condition in the plot that benefits or penalizes wheat grown the following season?

The package also supports a more specific compatibility question:

> Are some observed lentil-wheat genotype pairs better or worse than expected from the average behavior of the lentil and wheat genotypes separately?

These two questions are related, but they are not identical. The first asks about the average legacy of a lentil genotype. The second asks whether a specific lentil-wheat pairing has additional value beyond additive genotype effects.

## Experimental Design

The field layout contains 76 rows by 30 columns, for 2,280 plots per trial layout, organized into incomplete blocks. The experiment includes:

-   100 lentil genotypes in Year 1.
-   100 wheat genotypes in Year 2.
-   Approximately 20 observations per genotype per location for the full-trial phenotypic layers.
-   A sparse lentil-wheat pairing network rather than a complete factorial design.

A complete lentil by wheat factorial would require:

``` text
100 lentil genotypes x 100 wheat genotypes = 10,000 genotype combinations
```

before replication and locations. That is not feasible in the available field space. Instead, the design uses a sparse or incomplete factorial structure. Each lentil genotype is paired with a subset of wheat genotypes, often around 10 wheat partners. For example:

``` text
L001 paired with W001-W010
L002 paired with W011-W020
...
L100 paired with W091-W100
```

The wheat file stores this grouping in the `Facet` column. A facet is a wheat-partner network or wheat block of genotype partners. This matters because lentils in different facets are not necessarily connected through the same wheat genotype set.

## Data Layers

The package is designed around three major data layers.

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

This layer has a lower-density sampling structure, with treatment observations plus repeated check observations nested in blocks. It is analytically useful, but its design differs from the full phenotypic layer, so checks and block structure need careful handling.

### Wheat Full Trial

The wheat layer includes Year 2 traits such as:

``` text
HD, HT, MAT, LD, YLD, Y.ADJ, TWT, KWT, PRO
```

This layer contains the current wheat genotype, the previous lentil genotype, the observed lentil-wheat combination, the wheat facet, spatial coordinates, block information, and agronomic response traits.

## Main Modeling Challenges

### 1. Spatial Heterogeneity

The field has row and column structure. Yield and related traits can vary across the field because of soil gradients, moisture, fertility, management effects, or other spatial patterns. If spatial trends are ignored, genotype effects can be confounded with field position.

The package uses SpATS spatial models where possible:

``` r
trait ~ fixed design terms + genotype effect + PSANOVA(Row, Col)
```

The spatial surface helps separate genotype signal from field heterogeneity.

### 2. Incomplete Lentil-Wheat Pairing

Not every lentil is paired with every wheat. This creates an incomplete genotype network. A raw mean for a lentil can be biased if that lentil happened to be paired with unusually high-performing or low-performing wheat genotypes.

For example, if L001 is only paired with W001-W010, and W001-W010 are generally low-yielding wheat genotypes, L001's raw Year 2 wheat performance could look poor even if its true legacy effect is positive.

This is why wheat genotype must be included in the model.

### 3. Facet-Specific Baselines

Because the pairing network is structured by `Facet`, comparing every lentil against the full environment mean can be misleading. A lentil in the W001-W010 facet should be compared against the expected performance of that facet, not necessarily against the mean of all 100 wheat genotypes.

The current preferred predecessor-effect output therefore reports:

``` text
Legacy_Value = Corrected_Mean - mean(Corrected_Mean within ENV x Facet)
```

It also keeps a global environment deviation for comparison:

``` text
Legacy_Value_Global = Corrected_Mean - mean(Corrected_Mean within ENV)
```

The facet-aware value is the preferred interpretation for the sparse pairing design.

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

Because the design is sparse, pair compatibility should be interpreted within the observed network. It should not be treated as evidence for unobserved lentil-wheat combinations.

## Formula Layer 1: Independent Trait Models

The generic trait model is implemented in `model_traits()`.

### SpATS Trait Model

For lentil or wheat phenotyping, the spatial model is conceptually:

``` r
trait ~ design terms + genotype random effect + PSANOVA(Row, Col)
```

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
  rep_col = "Rep_combo",
  spatial_cols = c("Row", "Col")
)
```

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

## Formula Layer 2: Average Lentil Predecessor Effect

The preferred function is:

``` r
model_predecessor_effect()
```

Conceptually, the model is:

``` r
wheat_trait ~ previous_lentil + design terms + wheat genotype random effect + spatial trend
```

In the SpATS implementation:

``` r
Y.ADJ ~ Lentil + Rep_combo + PSANOVA(Row, Col) + random(Wheat)
```

The wheat genotype is modeled as a random effect. This is the core correction for the incomplete pairing network. If a lentil was paired with low-yielding wheat genotypes, the wheat random effects help correct the lentil estimate upward. If it was paired with high-yielding wheat genotypes, the model corrects downward.

After the model estimates corrected lentil means, the preferred legacy value is calculated relative to the facet baseline:

``` text
Legacy_Value = Corrected_Mean - mean(Corrected_Mean within ENV x Facet)
```

The raw-vs-corrected output also reports:

``` text
Network_Correction = Corrected_Mean - Raw_Mean
```

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
-   `Baseline_Mean`: corrected mean within the ENV x Facet baseline.
-   `Legacy_Value`: corrected deviation from the ENV x Facet baseline.
-   `Legacy_Value_Global`: corrected deviation from the full ENV mean.
-   `Network_Correction`: difference between corrected and raw means.

### Insight

This formula answers:

``` text
Which lentil genotypes leave better or worse wheat-growing conditions within their wheat-facet network?
```

The main plot shows `Legacy_Value`, centered on the ENV x Facet mean.

The correction plot shows:

-   Hollow point: raw mean.
-   Solid point: corrected mean.
-   Dashed vertical line: facet baseline.

This helps diagnose whether the model is correcting a lentil upward or downward because of its wheat partners and spatial position.

### GWAS-ready predecessor phenotypes

For lentil GWAS, the useful output is one value per lentil genotype, environment, and trait. `model_predecessor_effect()` now returns this table directly:

``` r
avg_predecessor_yadj$gwas_phenotypes
```

The key column is:

``` text
Predecessor_Phenotype
```

By default, this is the same biological quantity as `Legacy_Value`: the model-corrected wheat response after a lentil genotype, expressed as a deviation from the relevant ENV x Facet baseline.

This phenotype can be interpreted as:

``` text
How much better or worse did wheat perform after this lentil genotype than expected for its wheat-partner facet?
```

Positive values indicate lentil genotypes associated with improved following-wheat performance. These values can be used as Year 2 legacy phenotypes for association with lentil genomic markers. The biological target is a lentil genomic signature associated with improving the next crop, rather than improving the lentil crop itself.

The ranked GWAS-style plot is available as:

``` r
avg_predecessor_yadj$ranked_plot
plot_predecessor_gwas_phenotypes(avg_predecessor_yadj)
```

This plot puts all lentil genotypes into one ranked list per environment, while coloring bars by `Facet` to preserve the sparse-design context.

## Formula Layer 3: Lentil-Wheat Pair Compatibility

The preferred function is:

``` r
model_pair_compatibility()
```

This function asks whether specific observed lentil-wheat pairs have extra compatibility beyond additive genotype effects.

The conceptual formula is:

``` r
wheat_trait ~ Lentil + Wheat + spatial terms + random(Combo) + random(Block)
```

The `Combo` random effect is the pair-specific compatibility value.

The current implementation fits this model within each ENV x Facet network. That is important because the full 100 x 100 factorial is not observed and the facets are not fully connected through all genotype combinations.

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
-   `Corrected_Pair_Mean`: predicted value including the pair effect.
-   `Compatibility_Value`: pair-specific deviation from the additive expectation.
-   `Compatibility_Pct`: compatibility as a percent of the additive expectation.

### Insight

This formula answers:

``` text
Which observed lentil-wheat pairs are unusually good or bad after accounting for the average lentil effect, average wheat effect, and spatial field position?
```

Interpretation should stay inside the observed network. A positive pair effect means the observed pair performed better than expected for that specific lentil and wheat combination. It does not directly predict unobserved pairings.

## Formula Layer 4: Older Legacy Wrappers

The older functions:

``` r
model_legacy_spats()
model_legacy_spats2()
```

now call the newer predecessor-effect engine. They are retained for backward compatibility.

`model_legacy_spats()` returns the facet-aware predecessor plot.

`model_legacy_spats2()` returns the raw-vs-corrected correction plot as its main plot.

These wrappers are useful if older draft scripts already call them, but new analyses should prefer:

``` r
model_predecessor_effect()
```

## Formula Layer 5: Facet-Corrected Legacy Alternative

The function:

``` r
model_legacy_facet()
```

implements an older alternative approach:

1.  Estimate wheat genotype baseline performance.
2.  Estimate lentil-associated observed wheat performance.
3.  Compare lentil performance against the expected mean of its specific wheat partners.

Conceptually:

``` text
Legacy_Value = Observed_Lentil_Associated_Wheat_Performance - Facet_Mean
```

### Insight

This is a useful sensitivity analysis because it makes the facet correction explicit. However, the preferred model is currently `model_predecessor_effect()` because it directly models previous lentil genotype, current wheat genotype, and spatial trend in one model.

## Formula Layer 6: Correlating Lentil Traits with Wheat Legacy

Once lentil BLUPs and wheat legacy values are estimated, the package can combine them into a long table and use:

``` r
plot_legacy_correlations()
```

This creates a correlation heatmap across traits.

### Insight

This analysis asks:

``` text
Which Year 1 lentil traits are associated with positive or negative Year 2 wheat legacy effects?
```

Examples of possible biological interpretations include:

-   Lentil yield or biomass could correlate with wheat legacy if residue quantity matters.
-   Lentil maturity or phenology could correlate with wheat legacy if timing affects soil water or residue decomposition.
-   Lentil protein or nutrient traits could correlate with wheat legacy if nitrogen or nutrient cycling contributes to the effect.

These correlations are exploratory. They can identify hypotheses, but they do not prove mechanism without follow-up modeling or experimental validation.

## Recommended Interpretation Workflow

1.  Inspect trait distributions and field structure.
2.  Fit independent lentil and wheat trait models.
3.  Examine spatial trends to confirm field correction is needed.
4.  Estimate average lentil predecessor effects with `model_predecessor_effect()`.
5.  Use `Legacy_Value` for the facet-aware interpretation.
6.  Use `Legacy_Value_Global` only as a secondary comparison.
7.  Use the raw-vs-corrected plot to show how much the model corrected for wheat partners and spatial structure.
8.  Estimate pair compatibility with `model_pair_compatibility()`.
9.  Interpret pair compatibility only within observed ENV x Facet networks.
10. Correlate lentil BLUPs with wheat legacy values to generate biological hypotheses.

## Practical Notes

-   Check plots should usually be excluded from rotation-effect models unless the question specifically concerns checks.
-   `Facet` should be used as the comparison baseline for predecessor legacy values.
-   `Facet` should not be naively added as a fixed effect if it is redundant with genotype structure.
-   Sparse pair compatibility results are conditional on the observed design.
-   Singular fits can occur in pair models when variance components are small or a facet has too little information. This is expected in sparse networks and should be interpreted carefully.
-   The testing script `test_rotation_formulas_and_plots.R` is the best place to validate formulas and plots before promoting analyses into package vignettes or manuscript-ready scripts.
