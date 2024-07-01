---
title: 'CRAN Task View: Power analysis and sample size calculation'
output: 
  html_document: 
    highlight: null
    self_contained: no
    theme: null
---

<!--
<name>PowerAnalysis</name>  
<topic>Power analysis and sample size calculation</topic>  
<maintainer email="rwebster@cheo.on.ca">Richard Webster</maintainer>  
<maintainer email="ethancbrown@gmail.com">Ethan C. Brown</maintainer>  
<version>2024-07-01</version>  

<info>
-->

Power calculations and sample size estimation are an important step in
study design. It is costly to collect data, and either too little or too
much data collection may result in either unreliable findings or precise
estimates of the true effect at a higher than necessary cost. By making
educated assumptions at the study design phase, statistical power and
related design considerations can:

-   Inform study feasibility

-   Determine the point of diminishing returns for sampling, and

-   Allow decision makers to better allocate their limited resources.

There are many statistical tests and study designs, so subsequently
there is a vast array of statistical power calculations. Power
calculations may be an analytic formula approach, a simulation
approach, or a combination of the two. While the analytic formula
approach has a theoretical foundation, often for complex models there
are no analytical solutions available.

We start with several general-purpose simulation packages that are useful for 
power analysis. Although simulation-based power analysis can be performed in 
base R, several packages provide frameworks that make it easier to manage the 
full workflow of implementing an arbitrary data-generating mechanism, varying 
parameters, and gathering/displaying results.

The remainder of the Task View is devoted to power/study design packages that 
are specific to a specific type of model or analysis technique. These are
organized by overall category. We conclude with several useful links for 
performing analysis in R.

Please let us know if we have omitted something of importance,
for instance a new package that should be mentioned here, by submitting an issue or pull request on our [Github page](https://github.com/statisfactions/ctv-power).

#### General simulation study packages useful for power analysis

General-purpose simulation study packages facilitate understanding of study design 
and can accomodated arbitrary data-generating mechanisms and model-fitting
techniques. Although many packages can be useful at various stages of a 
simulation study, we only include packages here that are specifically for managing
the workflow of a simulation study. Packages that use simulation methods to
evaluate a specific type of model are integrated in to the **<u>Specific models/methods/tests</u>** section below.

-   [DeclareDesign](https://CRAN.R-Project.org/package=DeclareDesign) provides a general framework for specifying research designs, simulated based on specified designs, and assessing various properties of the designs, including but not limited to power.

-   [simpr](https://github.com/statisfactions/simpr) is a GitHub R package that has a framework to data simulation—based on tidyverse / broom syntax—and can perform power calculations

-   [MonteCarlo](https://CRAN.R-Project.org/package=MonteCarlo) provides functions for running
custom data-generating functions on a parameter grid and summarizing the results in a LaTeX 
table.

-   [simstudy](https://CRAN.R-Project.org/package=simstudy) provides a syntax for
defining datasets and referencing parameter values when generated within a for() loop, 
and includes support for 15 different types of distributions that can be used in studies, with the assumption that the model and the population of the investigated study should match.


#### Specific models/methods/tests (ordered alphabetically)

This list of packages for power calculations / sample size estimation is
ordered alphabetically by statistical test.

-   **Adaptive Study Design**

    -   [spass](https://CRAN.R-Project.org/package=spass), [esDesign](https://CRAN.R-Project.org/package=esDesign) 

-   **Algorithm Optimization**

    -   [CAISEr](https://CRAN.R-Project.org/package=CAISEr) determines the sample size required for comparing the
        performance of a set of *k* algorithms on a given problem
        instance when power inaccuracy is predefined

-   **ANCOVA**

    -   [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl) computes power for one or two factor ANCOVA with a
        single covariate

-   **ANOVA:**

    -   General options include [pwr](https://CRAN.R-Project.org/package=pwr) (one-way ANOVA), [pwr2](https://CRAN.R-Project.org/package=pwr2) &
        [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl) (up to two-way ANOVA) and [WebPower](https://CRAN.R-Project.org/package=WebPower) (up to
        two-way and repeated measures). The [easypower](https://CRAN.R-Project.org/package=easypower) package, based
        on the [pwr](https://CRAN.R-Project.org/package=pwr) package, simplifies the user input for factorial
        ANOVA.

    -   [BUCSS](https://CRAN.R-Project.org/package=BUCSS) allows any number of factors, using uncertainty and
        publication bias correction.

    -   [powerbydesign](https://CRAN.R-Project.org/package=powerbydesign) provides functions for bootstrapping the power
        of ANOVA designs based on estimated means and standard
        deviations of the conditions

    -   [powerAnalysis](https://CRAN.R-Project.org/package=powerAnalysis) only should be used with balanced one-way
        analysis of variance tests.

    -   [powerMediation](https://CRAN.R-Project.org/package=powerMediation) performs power calculation for interaction
        effect in 2x2 two-way ANOVA.

    -   [Superpower](https://CRAN.R-Project.org/package=Superpower) uses simulations and analytic power solutions for
        ANOVA designs of up to three factors to calculate power and
        average observed effect sizes.

-   **Bayesian:**

    -   [bayescount](https://CRAN.R-Project.org/package=bayescount) provides analysis and power calculations for
        count data.

    -   [BayesESS](https://CRAN.R-Project.org/package=BayesESS) determines effective sample size of a parametric
        prior distribution in Bayesian conjugate models (beta-binomial,
        gamma-exponential, gamma-Poisson, dirichlet-multinomial,
        normal-normal, inverse chi-squared-normal,
        inverse-gamma-normal), Bayesian linear and logistic regression
        models, Bayesian continual reassessment method (CRM), and
        Bayesian time to event model.

    -   [BayesianPower](https://CRAN.R-Project.org/package=BayesianPower) determine the required sample size for the
        evaluation of inequality constrained hypotheses by means of a
        Bayes factor

    -   [SampleSizeMeans](https://CRAN.R-Project.org/package=SampleSizeMeans) calculates sample size requirements for estimating means, and [SampleSizeProportions](https://CRAN.R-Project.org/package=SampleSizeProportions) for estimating proportions, using
        three different Bayesian criteria in the context of designing an
        experiment. Functions for calculation of required sample
        sizes for the Average Length Criterion, the Average Coverage
        Criterion and the Worst Outcome Criterion. Functions for both the fully Bayesian
        and the mixed Bayesian/likelihood approaches are provided.

-   **Beta Distribution:**

    -   [PASSED](https://CRAN.R-Project.org/package=PASSED) and [BetaPASS](https://CRAN.R-Project.org/package=BetaPASS) helps find the power for a given
        sample size or the sample size given power, when the
        outcome/response/dependent variable is a proportion or
        constrained between 0-1, with a beta distribution. **BetaPASS**
        has functions to plot power curves demonstrated in its vignette.

-   **Bioequivalence Study**

    -   [Power2Stage](https://CRAN.R-Project.org/package=Power2Stage) contains functions to obtain the operational
        characteristics of bioequivalence studies in Two-Stage Designs
        (TSD) via simulations.

    -   [PowerTOST](https://CRAN.R-Project.org/package=PowerTOST) calculates power and sample size for various study
        designs used in bioequivalence studies.

-   **Case-Control study:**

    -   [CoRpower](https://CRAN.R-Project.org/package=CoRpower) calculates power for assessment of intermediate biomarker
        responses as correlates of risk in the active treatment group in
        clinical efficacy trials. The methods differ from past
        approaches by accounting for the level of clinical treatment
        efficacy overall and in biomarker response subgroups, which
        enables the correlates of risk results to be interpreted in
        terms of potential correlates of efficacy/protection.

    -   [epiR](https://CRAN.R-Project.org/package=epiR) calculates sample size, power, and detectable odds ratios.

    -   [samplesizelogisticcasecontrol](https://CRAN.R-Project.org/package=samplesizelogisticcasecontrol) determines sample size for
        case-control studies to be analyzed using logistic regression.

-   **Chi-squared test:**

    -   Options include [powerAnalysis](https://CRAN.R-Project.org/package=powerAnalysis), [pwr](https://CRAN.R-Project.org/package=pwr), [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl), [bimetallic](https://CRAN.R-Project.org/package=bimetallic), [ssd](https://CRAN.R-Project.org/package=ssd).

-   **Cochran-Mantel-Haenszel Test:**

    -   [samplesizeCMH](https://CRAN.R-Project.org/package=samplesizeCMH) calculates the power and sample size for
        Cochran-Mantel-Haenszel tests, with several helper functions for
        working with probability, odds, relative risk, and odds ratio
        values.

-   **Competing Risks Analysis:**

    -   [powerCompRisk](https://CRAN.R-Project.org/package=powerCompRisk) is power analysis tool for jointly testing the
        cause-1 cause-specific hazard and the any-cause hazard with
        competing risks data.

-   **Complex Surveys**

    -   [samplesize4surveys](https://CRAN.R-Project.org/package=samplesize4surveys) computes the required sample size for
        estimation of totals, means, and proportions under complex
        sampling designs.

-   **Cost-Effectiveness Studies**

    -   [EBASS](https://CRAN.R-Project.org/package=EBASS) calculates the estimated sample size based on the
        Expected Value of Perfect Information (EVPI).

-   **Cross-sectional study**

    -   [epiR](https://CRAN.R-Project.org/package=epiR) provides functionality for sample size, power, or detectable prevalence.

-   **Diagnostic Test**

    -   [MKpower](https://CRAN.R-Project.org/package=MKpower) computes sample size, power, delta, or significance level of
        a diagnostic test for an expected sensitivity or specificity.

-   **Dirichlet-Multinomial distribution:**

    -   [HMP](https://CRAN.R-Project.org/package=HMP) uses the Dirichlet-Multinomial distribution to provide
        several functions for formal hypothesis testing, power and
        sample size calculations.

-   **Factorial Design**

    -   [BDEsize](https://CRAN.R-Project.org/package=BDEsize) calculates the sample size required to detect a
        certain standardized effect size, under a significance level
        (two-level fractional factorial, randomized complete block
        design, full factorial design, and split-plot design). This
        package also provides three graphs; detectable standardized
        effect size vs power, sample size vs detectable standardized
        effect size, and sample size vs power, which show the mutual
        relationship between the sample size, power and the detectable
        standardized effect size.

    -   [H2x2Factorial](https://CRAN.R-Project.org/package=H2x2Factorial) estimates the required number of clusters or
        the achieved power level under different types of hypothesis
        tests in a hierarchical 2x2 factorial trial with unequal cluster
        sizes and a continuous outcome.

-   **Fisher’s Test**

    -   Power calculations for differences between binomial proportions
        can be achieved with [Exact](https://CRAN.R-Project.org/package=Exact) for unconditional exact tests
        with 2x2 contingency tables. [MIDN](https://CRAN.R-Project.org/package=MIDN) computes the exact sample
        sizes required based on the Boschloo’s technique and
        Fisher-Boschloo test mid-N estimates.

    -   [ssanv](https://CRAN.R-Project.org/package=ssanv) provides a calibrated power calculation by leveraging
        the uncertainty in either nonadherence or parameter estimation.

-   **Gamma Distribution**

    -   [PASSED](https://CRAN.R-Project.org/package=PASSED) determines the power for a test of two sample means
        with Gamma distributions, or determine parameters to obtain a
        target power.

-   -   

-   **Genetic Association Studies:**

    -   [genpwr](https://CRAN.R-Project.org/package=genpwr) is for genetic association studies, considering the
        impact of misspecification of the genetic model.

    -   [pbatR](https://CRAN.R-Project.org/package=pbatR) provides pedigree/family-based genetic association tests
        analysis and power.

    -   [powerGWASinteraction](https://CRAN.R-Project.org/package=powerGWASinteraction) provides analytical power calculations for GxE
        and GxG interactions for case-control studies of candidate genes
        and genome-wide association studies (GWAS). This includes power
        calculation for four two-step screening and testing procedures.
        It can also calculate power for GxE and GxG without any
        screening.

    -   [survSNP](https://CRAN.R-Project.org/package=survSNP) conducts asymptotic and empirical power and sample
        size calculations for Single-Nucleotide Polymorphism (SNP)
        association studies with right censored time to event outcomes.

-   **Goodness of Fit Test**

    -   [pwr](https://CRAN.R-Project.org/package=pwr) supports calculation of power for goodness of fit tests.

    -   [PoweR](https://CRAN.R-Project.org/package=PoweR) provides functions for the computation of power and level
        tables for hypothesis tests, in Latex format, functions to build
        explanatory graphs for studying power of test statistics.

-   **Group Sequential Design**

    -   [PwrGSD](https://CRAN.R-Project.org/package=PwrGSD) provides tools for the evaluation of interim analysis plans
        for sequentially monitored trials on a survival endpoint;
        deriving power of a sequential design at a specified
        alternative.

-   **Hierarchical Data**

    -   [HierO](https://CRAN.R-Project.org/package=HierO) calculates statistical power for given type I error
        (alpha), effect size (Delta) and non-centrality parameter
        (ncpar) of a non-central chi-square distribution.

-   **High Dimensional Classification Study**

    -   [HDDesign](https://CRAN.R-Project.org/package=HDDesign) determines the sample size requirement to achieve
        the target probability of correct classification (PCC) for
        studies employing high-dimensional features.

-   **Human Microbiome Experiment**

    -   [HMP](https://CRAN.R-Project.org/package=HMP) uses the Dirichlet-Multinomial distribution to provide
        several functions for formal hypothesis testing, power and
        sample size calculations for human microbiome experiments.

-   **Intraclass Correlation**

    -   [ICC.Sample.Size](https://CRAN.R-Project.org/package=ICC.Sample.Size) calculates power for given value of p, the
        null hypothesis p0, number of raters (k), number of
        comparisons (N) and alpha. Can also be used to calculate the
        effect of increasing N at given intervals to a maximum N, or to
        calculate the increase in sample size to obtain increasing power
        with a given maximum N.

-   **Interobserver Agreement Studies**

    -   [kappaSize](https://CRAN.R-Project.org/package=kappaSize) providesbasic tools for sample size estimation in studies
        of interob-server/interrater agreement (reliability). Includes
        functions for both the power-based and confidence
        interval-based methods, with binary or multinomial outcomes and
        two through six raters.

-   **Likelihood Ratio test**

    -   [asypow](https://CRAN.R-Project.org/package=asypow) calculates power utilizing asymptotic likelihood
        ratio methods

-   **Linear Regression**

    -   [pwr](https://CRAN.R-Project.org/package=pwr), [powerMediation](https://CRAN.R-Project.org/package=powerMediation) and [WebPower](https://CRAN.R-Project.org/package=WebPower) provide analytical
        power calculations.

    -   [simr](https://CRAN.R-Project.org/package=simr) and [simpr](https://github.com/statisfactions/simpr) provided simulated power calculations.

    -   [BayesESS](https://CRAN.R-Project.org/package=BayesESS) provided simulated power calculations for Bayesian
        models

-   **Local Average Treatment Effect (LATE)**

    -   [powerLATE](https://CRAN.R-Project.org/package=powerLATE) is an implementation of the generalized power
        analysis for the local average treatment effect (LATE). The
        method uses standardized effect sizes to place a conservative
        bound on the power under minimal assumptions. Package allows
        users to recover power, sample size requirements, or minimum
        detectable effect sizes. Package also allows users to work with
        absolute effects rather than effect sizes, to specify an
        additional assumption to narrow the bounds, and to incorporate
        covariate adjustment.

-   **Logistic Regression**

    -   [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl), [powerMediation](https://CRAN.R-Project.org/package=powerMediation) and [WebPower](https://CRAN.R-Project.org/package=WebPower) calculates power
        for simple logistic regression with binary/continuous predictors

    -   [simr](https://CRAN.R-Project.org/package=simr) provides simulated power calculations.

    -   [BayesESS](https://CRAN.R-Project.org/package=BayesESS) provided simulated power calculations for Bayesian
        models

-   **Longitudinal Data Analysis**

    -   [LPower](https://CRAN.R-Project.org/package=LPower) computes power, or sample size or the detectable
        difference for a repeated measures model with attrition. It
        requires the variance covariance matrix of the observations but
        can compute this matrix for several common random effects
        models.

    -   [longpower](https://CRAN.R-Project.org/package=longpower) computes power and sample size for linear models of
        longitudinal data. Supported models include mixed-effects models
        and models fit by generalized least squares and generalized
        estimating equations.

    -   [powerlmm](https://github.com/rpsychologist/powerlmm) calculates power for the 'time x treatment' effect
        in two- and three-level multilevel longitudinal studies with
        missing data. Studies with partially nested designs, unequal
        cluster sizes, unequal allocation to treatment arms, and
        different dropout patterns per treatment are supported. For all
        designs power can be calculated both analytically and via
        simulations.

    -   [powerMediation](https://CRAN.R-Project.org/package=powerMediation) calculates power for testing if mean
        changes for 2 groups are the same or not for longitudinal study
        with 2 or more time points.

    -   [ssrm.logmer](https://CRAN.R-Project.org/package=ssrm.logmer) estimates required sample size for a longitudinal study with
        binary outcome in order to attain a pre-specified power while
        strictly maintaining the Type I error rate.

    -   [WebPower](https://CRAN.R-Project.org/package=WebPower) also supports estimation of power for longitudinal models.

-   **MANOVA**

    -   [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl) supports power for one factor MANOVA with up to 2 levels and 4
        measures.
		
-   **Mcnemar Test**

    -   [MIDN](https://CRAN.R-Project.org/package=MIDN) computes the exact sample sizes required in the
        randomized UMPU test and its conservative nonrandomized
        counterpart for attaining prespecified power. However, in
        contrast to the parallel group setting, the midpoint of the
        interval between these two numbers shall now be used as a nearly
        exact value of the number of pairs to be observed in the
        asymptotic test based on the score-statistic corrected for
        possible exceedances of the nominal significance level.

-   **Mediation Analysis (see also Structural Equation Modeling)**

    -   [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl) and [WebPower](https://CRAN.R-Project.org/package=WebPower) both support power analysis of mediation.  WebPower additionally supports statistical power analysis for structural equation
        modeling / mediation based on Monte Carlo simulation.
		
-   **Meta-Analysis**

    -   [metapower](https://CRAN.R-Project.org/package=metapower) calculates power for meta-analysis, including power analysis
        of main effects.

    -   [metameta](https://github.com/dsquintana/metameta) is a Github package for re-analyzing published
        meta-analysis, including calculating power for each study in a
        meta-analysis.

-   **Microarray Study**

    -   [FDRsampsize](https://CRAN.R-Project.org/package=FDRsampsize) computes average power and sample size for
        microarray studies that use the false discovery rate as the
        final measure of statistical significance.

    -   [ssize.fdr](https://CRAN.R-Project.org/package=ssize.fdr) contains a set of functions that
        calculates appropriate sample sizes for one-sample t-tests,
        two-sample t-tests, and F-tests for microarray experiments based
        on desired power while controlling for false discovery rates.

-   **Micro-randomized trials**

    -   [MRTSampleSize](https://CRAN.R-Project.org/package=MRTSampleSize).

    -   [MRTSampleSizeBinary](https://CRAN.R-Project.org/package=MRTSampleSizeBinary) is a sample size calculator for
        micro-randomized trials (MRTs) with binary outcomes.

-   **Mixed Models**

    -   [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl) and [pass.lme](https://CRAN.R-Project.org/package=pass.lme) provide power and sample size
        calculation for the fixed effect term

    -   [pamm](https://CRAN.R-Project.org/package=pamm) assesses the power of a dataset to estimate
        significant random effects (intercept or slope) in a mixed
        model. The functions are based on the `lme4` and
        `lmerTest` packages.

    -   Simulations of lmms and glmms can be performed with [simr](https://CRAN.R-Project.org/package=simr) and
        [simglm](https://CRAN.R-Project.org/package=simglm).

-   **Multiple Comparisons**

    -   [rPowerSampleSize](https://CRAN.R-Project.org/package=rPowerSampleSize) computes sample size for single-step (Bonferroni)
        and step-wise procedures (Holm and Hochberg) to control the Type-II
        Generalized Family-Wise Error Rate.

    -   [FDRsampsize](https://CRAN.R-Project.org/package=FDRsampsize) and [pwrFDR](https://CRAN.R-Project.org/package=pwrFDR) provide functions to compute
        power and sample sizes for studies using false discovery rate as
        the final measure of statistical significance.

-   **Multivariable Prediction Model**

    -   [pmsampsize](https://CRAN.R-Project.org/package=pmsampsize) computes the minimum sample size required for the
        development of a new multivariable prediction model. This can be used
        to calculate the minimum sample size for the development of
        models with continuous, binary or survival (time-to-event)
        outcomes.

-   **Pearson Correlation**

    -   [WebPower](https://CRAN.R-Project.org/package=WebPower), [pwr](https://CRAN.R-Project.org/package=pwr), and [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl) all contain methods for power
        analysis of correlations.

-   **Pharmacokinetic Study Design:**

    -   [PharmPow](https://CRAN.R-Project.org/package=PharmPow) analyzes mixed (sparse/dense sampled) pharmacokinetic study
        designs.

-   **Poisson Distributions and Poisson Regression**

    -   [PASSED](https://CRAN.R-Project.org/package=PASSED) computes the power for a test of two sample means
        with Poisson distributions, or determine parameters to obtain a
        target power.

    -   [powerMediation](https://CRAN.R-Project.org/package=powerMediation) calculates power for simple Poisson
        regression & mediation effect in poisson regression based on
        Vittinghoff, Sen and McCulloch’s (2009) method.

    -   [RSPS](https://CRAN.R-Project.org/package=RSPS) provides functions for estimating power or sample size
        for RNA-Seq studies. Empirical approach is used and the data is
        assumed to be count in nature.

    -   [WebPower](https://CRAN.R-Project.org/package=WebPower) contains methods for estimating power for Poisson
        regression.

-   **Portfolio Efficiency**

    -   [GRS.test](https://CRAN.R-Project.org/package=GRS.test) implements power analysis for the GRS test for
        portfolio efficiency.

-   **Proportion Test**

    -   [BAEssd](https://CRAN.R-Project.org/package=BAEssd) employs a Bayesian average error based approach to
        sample size determination. Several functions are included for
        sample size calculation for common designs in clinical trials
        including one- and two-sample binary and normal responses.

    -   [MKpower](https://CRAN.R-Project.org/package=MKpower) contains methods for one-sample proportion tests
        while [PASSED](https://CRAN.R-Project.org/package=PASSED) contains methods for two, while
        [powerAnalysis](https://CRAN.R-Project.org/package=powerAnalysis), [pwr](https://CRAN.R-Project.org/package=pwr), [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl), and [WebPower](https://CRAN.R-Project.org/package=WebPower) support both.

    -   [binomSamSize](https://CRAN.R-Project.org/package=binomSamSize) computes confidence intervals and necessary
        sample sizes for the parameter *p* of the Bernoulli B(p)
        distribution under simple random sampling or under pooled
        sampling. Such computations are e.g. of interest when
        investigating the incidence or prevalence in populations.

-   **Quantitative Comparative Analysis**

    -   [qcapower](https://CRAN.R-Project.org/package=qcapower) estimates power of a sufficient term using
        permutation tests. A term can be anything: A condition,
        conjunction or disjunction of any combination of these.

-   **Randomized Controlled Trials**

	-   [CP](https://CRAN.R-Project.org/package=CP) is for calculating the condition power for different models in survival time
    analysis within RCTs with two different treatments to be compared
    and survival as an endpoint.

	-   [odr](https://CRAN.R-Project.org/package=odr) performs power analyses with and without accommodating cost
    structures of sampling for experimental studies under a budget
    constraint.

	-   [powerCompRisk](https://CRAN.R-Project.org/package=powerCompRisk) provides 2-group comparisons cause-specific hazard and the
	all-cause hazard sample size calculations using an asymptotic
    chi-square joint test, accounting for censoring (e.g., lost to
    follow-up, staggered entry and administrative censoring).

	-   [RCT](https://CRAN.R-Project.org/package=RCT) computes the minimum population needed to detect difference
    between control group and each treatment, given a target minimum
    detectable effect.

	-   [Sample.Size](https://CRAN.R-Project.org/package=Sample.Size) computes the required sample size using optimal
    designs with multiple constraints. This optimal method is designed
    for two-arm, randomized phase II clinical trials, and the required
    sample size can be optimized either using fixed or flexible
    randomization allocation ratios.

	-   [SampleSize4ClinicalTrials](https://CRAN.R-Project.org/package=SampleSize4ClinicalTrials) calculates sample size when comparing
    means or proportions in Phase III clinical trials with different
    research goals.

	-   [smartsizer](https://CRAN.R-Project.org/package=smartsizer) is a set of tools for determining the necessary sample
    size in order to identify the optimal dynamic treatment regime in a
    sequential, multiple assignment, randomized trial (SMART).

	-   [ssev](https://CRAN.R-Project.org/package=ssev) computes the optimal sample size for various 2-group
    designs (e.g., when comparing the means of two groups assuming equal
    variances, unequal variances, or comparing proportions) when the aim
    is to maximize the rewards over the full decision procedure of a)
    running a trial (with the computed sample size), and b) subsequently
    administering the winning treatment to the remaining N-n units in
    the population.

	-   [odr](https://CRAN.R-Project.org/package=odr) performs power analyses with and without accommodating cost
    structures of sampling for experimental studies under a budget
    constraint for multisite RCTs. [WebPower](https://CRAN.R-Project.org/package=WebPower) also supports multisite
    RCTs.

-   **Cluster Randomized Trials:**

    -   [clusterPower](https://CRAN.R-Project.org/package=clusterPower) provides calculations for wide array of cluster design
        trials (1) Simple two-arm comparison designs, (2) Difference in
        difference designs, (3) Individually randomized group treatment
        trials, (4) Stepped wedge designs, (5) Multiarm trial designs.

    -   [CRTSize](https://CRAN.R-Project.org/package=CRTSize) supports sample size estimation in cluster (group)
        randomized trials.

    -   [odr](https://CRAN.R-Project.org/package=odr) performs power analyses with and without accommodating
        cost structures of sampling for experimental studies under a
        budget constraint.

    -   [PowerUpR](https://CRAN.R-Project.org/package=PowerUpR) calculates statistical power, minimum detectable
        effect size (MDES), MDES dif-ference (MDESD), and minimum
        required sample size for various multilevel randomized
        experi-ments with continuous outcomes. Some of the functions can
        assist with planning two- and three-level cluster-randomized
        trials (CRTs) sensitive to multilevel moderation and mediation
        (2-1-1, 2-2-1, and 3-2-1).

    -   [SteppedPower](https://CRAN.R-Project.org/package=SteppedPower) provides tools for power and sample size calculation as
        well as design diagnostics for longitudinal mixed model
        settings, with a focus on stepped wedge designs.

    -   [swdpwr](https://CRAN.R-Project.org/package=swdpwr) proivdes statistical power calculation for stepped wedge
        cluster randomized trials.

    -   [WebPower](https://CRAN.R-Project.org/package=WebPower) also provides support for cluster randomized designs.

	-   **Adaptive Study Design**

		-   [spass](https://CRAN.R-Project.org/package=spass) and [esDesign](https://CRAN.R-Project.org/package=esDesign) provide sample size calculation where eligibility criteria of the trial is adaptively updated.

		-   [rpact](https://CRAN.R-Project.org/package=rpact) can be used for adaptive clinical trial design as well as many general functions for use outside an RCT context.
		
	-   **Equivalence trial**

		-   [epiR](https://CRAN.R-Project.org/package=epiR) calculates sample size for parallel equivalence trial with binary
        or continuous outcomes.

	-   **Non-Inferiority Trial**

		-   [epiR](https://CRAN.R-Project.org/package=epiR) supports non-inferiority trials with both binary or continuous outcomes.

		-   [blindrecalc](https://CRAN.R-Project.org/package=blindrecalc) computes key characteristics and plots
        for blinded sample size recalculation, including power. Continuous and binary endpoints are supported.

	-   **Sequential Multiple Assignment Randomized Trial (SMART)**

	    -   [SMARTbayesR](https://CRAN.R-Project.org/package=SMARTbayesR) supports optimal dynamic treatment regimes and sample
        size for a SMART design in the Bayesian setting.

		-   [SMARTp](https://CRAN.R-Project.org/package=SMARTp) provides sample size calculation to detect dynamic treatment
        regime (DTR) effects based on change in clinical attachment
        level (CAL) outcomes from a non-surgical chronic periodontitis
        treatments study. The clustered tooth (sub-unit) level CAL
        outcomes are skewed, spatially-referenced, and non-randomly
        missing.
	
	-   **Superiority Trial**

	    -   [blindrecalc](https://CRAN.R-Project.org/package=blindrecalc) computes key characteristics and plots
        for blinded sample size recalculation, including power.
        Continuous and binary endpoints are supported.

		-   [epiR](https://CRAN.R-Project.org/package=epiR) calculates sample size for parallel superiority trial with binary
			or continuous outcome)

		-   [MIDN](https://CRAN.R-Project.org/package=MIDN) provides nearly exact sample size calculation for exact
        powerful nonrandomized tests for differences between binomial
        proportions.
		
-   **Rasch Model**

    -   [pwrRasch](https://CRAN.R-Project.org/package=pwrRasch) provides statistical power simulation for testing the Rasch
        Model based on a three-way analysis of variance design with
        mixed classification.

-   **Regression (also see Mixed Models and Poisson Distribution)**

    -   [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl), [pwr](https://CRAN.R-Project.org/package=pwr), and [WebPower](https://CRAN.R-Project.org/package=WebPower) all have general methods for
        various forms of regression. `pwr2ppl::pwr.f2.test` provides power for
        comparing dependent coefficients in multiple regression with two
        or three predictors.

    -   [simr](https://CRAN.R-Project.org/package=simr) and [simglm](https://CRAN.R-Project.org/package=simglm) both support simulation of regression
        models with various structures.

    -   [Hmisc](https://CRAN.R-Project.org/package=Hmisc) simulates power of a Bayesian odds ordinal logistic
        model.

    -   [PASSED](https://CRAN.R-Project.org/package=PASSED) and [MKpower](https://CRAN.R-Project.org/package=MKpower) for comparing two negative binomial
        rates

    -   [RSPS](https://CRAN.R-Project.org/package=RSPS) provides an empirical approach for RNA-seq data

-   **Regression Discontinuity (RD) Design**

    -   Both [rdpower](https://CRAN.R-Project.org/package=rdpower) and [rddapp](https://github.com/felixthoemmes/rddapp) provide tools to perform power,
        sample size and MDE calculations in RD designs.

-   **Repertory Grid Studies**

    -   [gridsampler](https://CRAN.R-Project.org/package=gridsampler) is a simulation tool to facilitate determination
        of required sample size to achieve category saturation for
        studies using multiple repertory grids in conjunction with
        content analysis.

-   **Sib Pair Design**

	-   [powerpkg](https://CRAN.R-Project.org/package=powerpkg) computes power of an affected sib pair design.
	
-   **Sign Test**

    -   [FDRsampsize](https://CRAN.R-Project.org/package=FDRsampsize) provides methods for microarray studies that use
    the false discovery rate as the final measure of statistical
    significance.

-   **Sobel Test**

    -   [powerMediation](https://CRAN.R-Project.org/package=powerMediation) calculates power for testing a mediation
        effect based on Sobel’s test.

-   **Spatial Clustering**

    -   [sparrpowR](https://CRAN.R-Project.org/package=sparrpowR) calculates the statistical power to detect
        clusters using kernel-based spatial relative risk functions.

-   **Structural Equation Modelling (SEM)**

    -   [semPower](https://CRAN.R-Project.org/package=semPower) provides a-priori, post-hoc, and compromise
        power-analyses for structural equation models (SEM))

    -   [simsem](https://CRAN.R-Project.org/package=simsem) provides an easy framework for Monte Carlo simulation
        in structural equation modeling, which can be used for various
        purposes, such as such as model fit evaluation, power analysis,
        or missing data handling and planning.

-   **Survival analysis**

    -   [powerMediation](https://CRAN.R-Project.org/package=powerMediation) calculates power for testing a mediation
        effect in Cox regression based on Vittinghoff, Sen &
        Mc-Culloch’s (2009) method.

    -   [DelayedEffect.Design](https://CRAN.R-Project.org/package=DelayedEffect.Design) uses a piecewise weighted log-rank test
        to properly and efficiently incorporate the delayed effect into
        the study design.

    -   [NPHMC](https://CRAN.R-Project.org/package=NPHMC) calculates sample size of a survival trial with or
        without cure fractions

    -   [npsurvSS](https://CRAN.R-Project.org/package=npsurvSS) provides sample size and power calculation for
        comparing two survival curves, including the difference in (or
        ratio of) t-year survival, difference in (or ratio of) p-th
        percentile survival, difference in (or ratio of) restricted mean
        survival time, and the weighted log-rank test with allowance for
        flexible accrual, censoring, and survival, and supports
        simulation methods as well.

    -   [powerSurvEpi](https://CRAN.R-Project.org/package=powerSurvEpi) provides power and sample size calculation for
        survival analysis of epidemiological studies, taking into
        account the correlation between the covariate of interest and
        other covariates

    -   [SSRMST](https://CRAN.R-Project.org/package=SSRMST) calculates the study sample size and power in
        designing clinical trials using the difference in RMST. Two
        types of one-sided tests, non-inferiority and superiority tests,
        are prepared.

-   **T-Test**

    -   Many packages can calculate power for basic one- and two-sample
        t-tests: [FDRsampsize](https://CRAN.R-Project.org/package=FDRsampsize), [WebPower](https://CRAN.R-Project.org/package=WebPower), [pwr](https://CRAN.R-Project.org/package=pwr), [PASSED](https://CRAN.R-Project.org/package=PASSED), [powerAnalysis](https://CRAN.R-Project.org/package=powerAnalysis), [samplesize](https://CRAN.R-Project.org/package=samplesize), [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl), [pwrAB](https://CRAN.R-Project.org/package=pwrAB).

    -   [PowerTOST](https://CRAN.R-Project.org/package=PowerTOST) calculates power for the two-one-sided t-test

    -   Simulate power through [MKpower](https://CRAN.R-Project.org/package=MKpower) for classical (equal
        variances), Welch and Hsu T-tests.

    -   [ssizeRNA](https://CRAN.R-Project.org/package=ssizeRNA) computes appropriate sample sizes for two-sample t-test for
        RNA-seq experiments with fixed or varied set of parameters.

    -   [MKpower](https://CRAN.R-Project.org/package=MKpower) computes power for Welch and Hsu *t-*tests for
        balanced designs; [MKmisc](https://CRAN.R-Project.org/package=MKpower) can also be used for Welch *t*-tests as well as comparing negative binomial rates.

-   **TDT Study (transmission disequilibrium test)**

    -   [powerpkg](https://CRAN.R-Project.org/package=powerpkg) calculates of how many TDT trios need to be
        studied in order to have an adequately powered TDT study.

-   **Trend-in-Trend Design**

    -   [TrendInTrend](https://CRAN.R-Project.org/package=TrendInTrend) provides estimation of causal odds ratio and
        power calculation given trends in exposure prevalence and
        outcome frequencies of stratified data.

-   **Wald Test**

    -   [SteppedPower](https://CRAN.R-Project.org/package=SteppedPower) computes power of a Wald Test given a standard
        error, an effect size, the degrees of freedom of the
        t-distribution and a significance level

-   **Wilcoxon-Mann-Whitney Test (Wilcoxon Rank Sum Test)**

    -   [FDRsampsize](https://CRAN.R-Project.org/package=FDRsampsize).

    -   [MKpower](https://CRAN.R-Project.org/package=MKpower).

    -   [WMWssp](https://CRAN.R-Project.org/package=WMWssp) calculates the minimal sample size for the
        Wilcoxon-Mann-Whitney test that is needed for a given power and
        two sided type I error rate. The method works for metric data
        with and without ties, count data, ordered categorical data, and
        even dichotomous data.

-   **Wilcoxon Sign Rank Test**

    -   [BetaPASS](https://CRAN.R-Project.org/package=BetaPASS) will help find the power for a given sample size or
        the sample size given power, when testing the null hypothesis
        that the means for the control and treatment groups are equal
        against a two-sided alternative.

    -   [MKpower](https://CRAN.R-Project.org/package=MKpower).

    -   [samplesize](https://CRAN.R-Project.org/package=samplesize) works for categorical data and allows for ties.

<!-- Links Moved to Separate document, ctv-power-links.md -->
