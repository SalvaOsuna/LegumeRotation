---
---
---

# ACTIVATE Lentil-Wheat Rotation Experiment: Design, Modeling Challenges, and Analytical Targets

## Objectives

`LegumeRotation` is an R package designed to supports analysis of a two-year lentil-wheat rotation experiment. The package has two connected goals:

1.  Characterize phenotypic performance of lentil and wheat as independent field trials.
2.  Test whether the lentil genotype grown in a plot in Year 1 influences the agronomic performance of the wheat genotype grown in the same plot in Year 2.

The main rotation question is a predecessor or legacy-effect question:

> Does lentil genotype X leave a measurable condition in the plot that benefits or penalizes wheat grown the following season?

The package also supports a more specific compatibility question:

> Are some observed lentil-wheat genotype pairs better or worse than expected from the average behavior of the lentil and wheat genotypes separately?

These two questions are related, but they are not identical. The first asks about the average legacy of a lentil genotype. The second asks whether a specific lentil-wheat pairing has additional value beyond additive genotype effects.

## Experimental Design

The field layout contains 76 rows by 30 columns, for 2,280 plots per trial layout, organized into 48 incomplete blocks. The experiment includes:

-   100 lentil genotypes in Year 1.
-   100 wheat genotypes in Year 2.
-   20 observations per genotype per location for the full-trial phenotypic layers.
-   A sparse lentil-wheat pairing network rather than a complete factorial design.

A complete lentil by wheat factorial would require:

``` text
100 lentil genotypes x 100 wheat genotypes = 10,000 genotype combinations
```

before replication and locations. That is not feasible in the available field space. Instead, the design uses a sparse or incomplete factorial structure. Each lentil genotype is paired with a subset of 10 wheat partners, following this sctructure:

``` text
L001 paired with W001-W010
L002 paired with W011-W020
...
L100 paired with W091-W100
```

The wheat file stores this grouping in the `Facet` column. A facet is a wheat-partner network or wheat block of genotype partners. This matters because lentils in different facets are not necessarily connected through the same wheat genotype set.

## Data Layers

The package is designed around four major data layers.

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

This layer has a lower-density sampling structure ([2 reps x genotype + checks] x site), with treatment observations plus repeated check observations nested in blocks. It is analytically useful, but its design differs from the full phenotypic layer, so checks and block structure need careful handling.

### Lentil Microbiome Layer

The microbiome layer comes from the Year 1 lentil plots and includes plot-level 16S-derived predictors such as:

``` text
Observed ASV richness, Chao1, ACE, Shannon, Simpson, Pielou evenness,
Aitchison PCoA axes, microbiome legacy score, Rhizobiaceae abundance,
Bradyrhizobiaceae abundance, legume-symbiont abundance
```

For the legacy-driver analysis, these plot-level microbiome traits are aggregated to lentil genotype x environment means, then the environment labels are shifted from 2024 to the matching 2025 wheat environment. This lets the microbiome layer be treated as another Year 1 predictor of Year 2 wheat legacy.

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

Mathematically, within each environment:

$$
Y_i = \mu_e + L_{e,l(i)} + d_{e,r(i)} + W_{e,w(i)} + f_e(R_i, C_i) + \varepsilon_i
$$

where:

-   `y_i` is the observed Year 2 wheat trait value for plot `i`.
-   `mu_e` is the environment-specific mean.
-   `L_{e,l(i)}` is the effect of the previous lentil genotype `l` grown in that plot in Year 1. This is the predecessor effect of interest.
-   `d_{e,r(i)}` is the fixed design effect, such as `Rep_combo`, for plot `i`.
-   `W_{e,w(i)}` is the random effect of the current wheat genotype `w` in Year 2.
-   `f_e(row_i, col_i)` is the spatial field trend for the wheat trial.
-   `epsilon_i` is residual plot-level error.

The model-corrected predecessor mean for lentil genotype `l` is:

$$
\mathrm{Corrected\_Mean}_{e,l} = \hat{\mu}_e + \hat{L}_{e,l}
$$

For a lentil genotype assigned to wheat facet `b`, the facet baseline is:

$$
\mathrm{Baseline\_Mean}_{e,b}
= \frac{1}{n_{e,b}}\sum_{l \in (e,b)} \mathrm{Corrected\_Mean}_{e,l}
$$

The reported legacy quantities are:

$$
\mathrm{Legacy\_Value}_{e,l}
= \mathrm{Corrected\_Mean}_{e,l} - \mathrm{Baseline\_Mean}_{e,b(l)}
$$

$$
\mathrm{Legacy\_Value\_Global}_{e,l}
= \mathrm{Corrected\_Mean}_{e,l}
- \frac{1}{n_e}\sum_{l \in e} \mathrm{Corrected\_Mean}_{e,l}
$$

$$
\mathrm{Network\_Correction}_{e,l}
= \mathrm{Corrected\_Mean}_{e,l} - \mathrm{Raw\_Mean}_{e,l}
$$

where:

-   `b(l)` is the wheat-partner facet associated with lentil genotype `l`.
-   `Raw_Mean_{e,l}` is the uncorrected mean wheat performance observed after lentil genotype `l`.
-   `Network_Correction_{e,l}` measures how much the model shifted the raw lentil-associated wheat mean after accounting for wheat partners, design terms, and spatial position.

The wheat genotype is modeled as a random effect. This is the core correction for the incomplete pairing network. If a lentil was paired with low-yielding wheat genotypes, the wheat random effects help correct the lentil estimate upward. If it was paired with high-yielding wheat genotypes, the model corrects downward.

After the model estimates corrected lentil means, the preferred legacy value is calculated relative to the facet baseline:

$$
\mathrm{Legacy\_Value}
= \mathrm{Corrected\_Mean}
- \operatorname{mean}_{\mathrm{ENV}\times\mathrm{Facet}}\left(\mathrm{Corrected\_Mean}\right)
$$

The raw-vs-corrected output also reports:

$$
\mathrm{Network\_Correction}
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

### GWAS-ready legacy values

For lentil GWAS, the useful output is one `Legacy_Value` per lentil genotype, environment, and wheat trait. `model_predecessor_effect()` returns these values in:

``` r
avg_predecessor_yadj$legacy_values
```

The key column is:

``` text
Legacy_Value
```

This is the model-corrected wheat response after a lentil genotype, expressed as a deviation from the relevant ENV x Facet baseline.

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

This plot puts all lentil genotypes into one ranked list per environment, while coloring bars by `Facet` to preserve the sparse-design context.

The SpATS output includes prediction standard errors for the adjusted lentil predecessor means. These are kept in the output table as `SE`, but they are not drawn by default in the ranked plots because `Legacy_Value` is a deviation from a facet baseline and the SE bars can visually overwhelm the ranking.

Boxplots are better reserved for plot-level raw or residualized observations. The `legacy_values` table has one adjusted legacy value per lentil genotype x environment x trait, so a boxplot would imply a distribution that is not present in that output.

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

Mathematically, within each ENV x Facet network:

$$
Y_i = \mu_{e,b} + L_{e,b,l(i)} + W_{e,b,w(i)} + d_{e,b,r(i)}
      + f_{e,b}(R_i, C_i) + C_{e,b,c(i)} + B_{e,b,k(i)} + \varepsilon_i
$$

where:

-   `y_i` is the observed Year 2 wheat trait value for plot `i`.
-   `mu_{e,b}` is the mean for environment `e` and facet `b`.
-   `L_{e,b,l(i)}` is the additive effect of the previous lentil genotype `l`.
-   `W_{e,b,w(i)}` is the additive effect of the current wheat genotype `w`.
-   `d_{e,b,r(i)}` represents any optional fixed design terms.
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
+ \hat{d}_{e,b,r(i)} + \hat{f}_{e,b}(R_i, C_i)
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

Positive `Compatibility_Value` means the observed pair performed better than expected from the additive lentil and wheat effects in that ENV x Facet network.

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

The default pair-compatibility plot is now a heatmap rather than a long bar plot. Each panel is one observed ENV x Facet network, the x-axis is the current wheat genotype, the y-axis is the previous lentil genotype, and the tile color is `Compatibility_Value`.

``` r
pair_compatibility_pro$heatmap
pair_compatibility_pro$ranked_plot
pair_compatibility_pro$diagnostics
```

The heatmap avoids saturating the y-axis with long lentil-wheat combination names. The ranked plot is still available for the strongest positive and negative pair effects, but it is intended as a summary view, not as the main map of the design.

If a panel shows all pair effects as zero, that does not mean the underlying observations are missing. It means the mixed model estimated the `Combo` variance component as zero for that ENV x Facet. In that situation, the fitted model says that the additive lentil effect, additive wheat effect, spatial terms, and design random effects explain the available signal, and there is no detectable pair-specific deviation left for that facet. The pair BLUPs are therefore all shrunk to zero.

By default, pair-compatibility plots now drop ENV x Facet groups with no detectable pair-specific signal. These groups remain in `pair_compatibility_pro$diagnostics`, but they are not shown in `heatmap` or `ranked_plot`. To include them for auditing:

``` r
plot_pair_compatibility_heatmap(pair_compatibility_pro, drop_no_signal = FALSE)
plot_pair_compatibility_ranked(pair_compatibility_pro, drop_no_signal = FALSE)
```

The diagnostics table should be checked before interpreting pair effects. Useful columns include:

-   `N_Plots`
-   `N_Previous_Genotypes`
-   `N_Current_Genotypes`
-   `N_Combos`
-   `Min_Reps_Per_Combo`
-   `N_Single_Rep_Combos`
-   `Status`
-   `Model_Singular`
-   `Combo_Variance`
-   `Residual_Variance`

For this design, a reliable ENV x Facet group should usually have 10 lentils, 10 wheats, 100 observed combos, and at least two plots per combo. If a panel appears empty in an old bar plot, that is usually a visualization artifact from crossing all environment and facet labels. The new heatmap uses only observed ENV x Facet panels.

## Formula Layer 4: Correlating Lentil Traits with Wheat Legacy

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
-   `r_{e,p,t}` is the Pearson correlation between a Year 1 lentil predictor and a Year 2 wheat legacy target.
-   `R2_{e,p,t}` is the proportion of target variation explained by that single predictor in the environment-specific exploratory correlation.
-   p-values and FDR values summarize statistical evidence for the correlation, but they do not establish causality.

When many subsample traits are included, the focused driver plot is easier to read than the full heatmap:

``` r
legacy_driver_pro <- plot_legacy_driver_correlations(
  data = legacy_correlation_input,
  target_trait = "Wheat_Legacy_PRO",
  top_n = 20,
  min_pairs = 10
)

legacy_driver_pro$correlations
legacy_driver_pro$plot
```

The driver table reports the environment-specific correlation, R2, p-value, FDR, and number of genotype pairs used for each lentil predictor trait.

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
