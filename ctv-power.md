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
<maintainer email="ecbrown@umn.edu">Ethan C. Brown</maintainer>  
<version>2021-10-20</version>  

<info><p>
Functional data analysis (FDA) deals with data that <a href="https://en.wikipedia.org/wiki/Functional_data_analysis">"provides information about curves, surfaces or anything else varying over a continuum."</a> This task view catalogues available packages in this rapidly developing field.</p>
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
calculations may be I) an analytic formula approach, II) a simulation
approach, or III) a combination of the two. While the analytic formula
approach has a theoretical foundation, often for complex models there
are no analytical solutions available. This page covers both analytical
and simulation-based methods for power analysis, as well as several more
general simulation packages that are useful for power analysis.

This Task View page hopes to make it easier to find sample size
functionality in R, and has the following sections:

**<u>Specific models/methods/tests (ordered alphabetically)</u>**

**<u>General simulation packages useful for power analysis</u>** 

**<u>Reporting power analysis</u>**

**<u>Discipline specific solutions</u>**

**<u>CRAN packages:</u>**

**<u>Related links</u>**

*<u>Please let us know</u>* if we have omitted something of importance,
for instance a new package that should be mentioned here (rwebster AT
cheo DOT on DOT ca).

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

    -   [SampleSizeMeans](https://CRAN.R-Project.org/package=SampleSizeMeans) (calculating sample size requirements using
        three different Bayesian criteria in the context of designing an
        experiment to estimate a normal mean orthe difference between
        two normal means. Functions for calculation of required sample
        sizes for the Average Length Criterion, the Average Coverage
        Criterion and the Worst Outcome Criterion in the context of
        normal means are provided. Functions for both the fully Bayesian
        and the mixed Bayesian/likelihood approaches are provided.)

    -   [SampleSizeProportions](https://CRAN.R-Project.org/package=SampleSizeProportions) (calculating sample size requirements
        using three different Bayesian criteria in the context of
        designing an experiment to estimate the difference between two
        binomial proportions. Functions for calculation of required
        sample sizes for the Average Length Criterion, the Average
        Coverage Criterion and the Worst Outcome Criterion in the
        context of binomial observations are provided. In all cases,
        estimation of the difference between two binomial proportions is
        considered. Functions for both the fully Bayesian and the mixed
        Bayesian/likelihood approaches are provided.)

-   **Beta Distribution:**

    -   [PASSED](https://CRAN.R-Project.org/package=PASSED) and [BetaPASS](https://CRAN.R-Project.org/package=BetaPASS) helps find the power for a given
        sample size or the sample size given power, when the
        outcome/response/dependent variable is a proportion or
        constrained between 0-1, with a beta distribution. **BetaPASS**
        has functions to plot power curves demonstrated in its vignette.

-   **Bioequivalence Study**

    -   [Power2Stage](https://CRAN.R-Project.org/package=Power2Stage) (Contains functions to obtain the operational
        characteristics of bioequivalence studies in Two-Stage Designs
        (TSD) via simulations.)

    -   [PowerTOST](https://CRAN.R-Project.org/package=PowerTOST) (Calculate power and sample size for various study
        designs used in bioequivalence studies)

-   **Case-Control study:**

    -   [CoRpower](https://CRAN.R-Project.org/package=CoRpower) (power for assessment of intermediate biomarker
        responses as correlates of risk in the active treatment group in
        clinical efficacy trials. The methods differ from past
        approaches by accounting for the level of clinical treatment
        efficacy overall and in biomarker response subgroups, which
        enables the correlates of risk results to be interpreted in
        terms of potential correlates of efficacy/protection)

    -   [epiR](https://CRAN.R-Project.org/package=epiR) (sample size, power, and detectable odds ratio)

    -   [samplesizelogisticcasecontrol](https://CRAN.R-Project.org/package=samplesizelogisticcasecontrol) (determine sample size for
        case-control studies to be analyzed using logistic regression)

-   **Chi-squared test:**

    -   Options include **powerAnalysis, pwr, pwr2ppl, bimetallic, ssd**

-   **Cochran-Mantel-Haenszel Test:**

    -   [sampsizeCMH](https://CRAN.R-Project.org/package=sampsizeCMH) calculates the power and sample size for
        Cochran-Mantel-Haenszel tests, with several helper functions for
        working with probability, odds, relative risk, and odds ratio
        values.

-   **Competing Risks Analysis:**

    -   [powerCompRisk](https://CRAN.R-Project.org/package=powerCompRisk) (A power analysis tool for jointly testing the
        cause-1 cause-specific hazard and the any-cause hazard with
        competing risks data.)

-   **Complex Surveys**

    -   [samplesize4surveys](https://CRAN.R-Project.org/package=samplesize4surveys) (computes the required sample size for
        estimation of totals, means, and proportions under complex
        sampling designs)

-   **Cost-Effectiveness Studies**

    -   [EBASS](https://CRAN.R-Project.org/package=EBASS) (calculate the estimated sample size based on the
        Expected Value of Perfect Information (EVPI))

-   **Cross-sectional study**

    -   [epiR](https://CRAN.R-Project.org/package=epiR) (sample size, power, or detectable prevalence)

-   **Diagnostic Test**

    -   [MKpower](https://CRAN.R-Project.org/package=MKpower) sample size, power, delta, or significance level of
        a diagnostic test for an expected sensitivity or specificity.

-   **Dirichlet-Multinomial distribution:**

    -   [HMP](https://CRAN.R-Project.org/package=HMP) using Dirichlet-Multinomial distribution to provide
        several functions for formal hypothesis testing, power and
        sample size calculations.

-   **Factorial Design**

    -   [BDEsize](https://CRAN.R-Project.org/package=BDEsize) (calculates the sample size required to detect a
        certain standardized effect size, under a significance level
        (two-level fractional factorial, randomized complete block
        design, full factorial design, and split-plot design). This
        package also provides three graphs; detectable standardized
        effect size vs power, sample size vs detectable standardized
        effect size, and sample size vs power, which show the mutual
        relationship between the sample size, power and the detectable
        standardized effect size)

    -   [H2x2Factorial](https://CRAN.R-Project.org/package=H2x2Factorial) (estimates the required number of clusters or
        the achieved power level under different types of hypothesis
        tests in a hierarchical 2x2 factorial trial with unequal cluster
        sizes and a continuous outcome)

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

    -   [genpwr](https://CRAN.R-Project.org/package=genpwr) (for genetic association studies, considering the
        impact of misspecification of the genetic model.)

    -   [pbatR](https://CRAN.R-Project.org/package=pbatR) (Pedigree/Family-Based Genetic Association Tests
        Analysis and Power)

    -   [powerGWASinteraction](https://CRAN.R-Project.org/package=powerGWASinteraction) (Analytical power calculations for GxE
        and GxG interactions for case-control studies of candidate genes
        and genome-wide association studies (GWAS). This includes power
        calculation for four two-step screening and testing procedures.
        It can also calculate power for GxE and GxG without any
        screening.)

    -   [survSNP](https://CRAN.R-Project.org/package=survSNP) (Conduct asymptotic and empirical power and sample
        size calculations for Single-Nucleotide Polymorphism (SNP)
        association studies with right censored time to event outcomes.)

-   **Goodness of Fit Test**

    -   [pwr](https://CRAN.R-Project.org/package=pwr) 

    -   [PoweR](https://CRAN.R-Project.org/package=PoweR) functions for the computation of power and level
        tables for hypothesis tests, in Latex format, functions to build
        explanatory graphs for studying power of test statistics.

-   **Group Sequential Design**

    -   [PwrGSD](https://CRAN.R-Project.org/package=PwrGSD) (Tools for the evaluation of interim analysis plans
        for sequentially monitored trials on a survival endpoint;
        deriving power of a sequential design at a specified
        alternative)

-   **Hierarchical Data**

    -   [HierO](https://CRAN.R-Project.org/package=HierO) (Calculates statistical power for given type I error
        (alpha), effect size (Delta) and non-centrality parameter
        (ncpar) of a non-central chi-square distribution)

-   **High Dimensional Classification Study**

    -   [HDDesign](https://CRAN.R-Project.org/package=HDDesign) (Determine the sample size requirement to achieve
        the target probability of correct classification (PCC) for
        studies employing high-dimensional features.)

-   **Human Microbiome Experiment**

    -   [HMP](https://CRAN.R-Project.org/package=HMP) (Using Dirichlet-Multinomial distribution to provide
        several functions for formal hypothesis testing, power and
        sample size calculations for human microbiome experiments)

-   **Interclass Correlation**

    -   [ICC.Sample.Size](https://CRAN.R-Project.org/package=ICC.Sample.Size) calculates power for given value of p, the
        null hypothesis p0, number of raters (k), number of
        comparisons (N) and alpha. Can also be used to calculate the
        effect of increasing N at given intervals to a maximum N, or to
        calculate the increase in sample size to obtain increasing power
        with a given maximum N.

-   **Interobserver Agreement Studies**

    -   [kappaSize](https://CRAN.R-Project.org/package=kappaSize) (basic tools for sample size estimation in studies
        of interob-server/interrater agreement (reliability). Includes
        functions for both the power-based and confi-dence
        interval-based methods, with binary or multinomial outcomes and
        two through six raters.)

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

    -   [powerLATE](https://CRAN.R-Project.org/package=powerLATE) (An implementation of the generalized power
        analysis for the local average treatment effect (LATE). The
        method uses standardized effect sizes to place a conservative
        bound on the power under minimal assumptions. Package allows
        users to recover power, sample size requirements, or minimum
        detectable effect sizes. Package also allows users to work with
        absolute effects rather than effect sizes, to specify an
        additional assumption to narrow the bounds, and to incorporate
        covariate adjustment.)

-   **Logistic Regression**

    -   [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl), [powerMediation](https://CRAN.R-Project.org/package=powerMediation) and [WebPower](https://CRAN.R-Project.org/package=WebPower) calculates power
        for simple logistic regression with binary/continuous predictors

    -   [simr](https://CRAN.R-Project.org/package=simr) provides simulated power calculations.

    -   [BayesESS](https://CRAN.R-Project.org/package=BayesESS) provided simulated power calculations for Bayesian
        models

-   **Longitudinal Data Analysis**

    -   [LPower](https://CRAN.R-Project.org/package=LPower) (Computes power, or sample size or the detectable
        difference for a repeated measures model with attrition. It
        requires the variance covariance matrix of the observations but
        can compute this matrix for several common random effects
        models)

    -   [longpower](https://CRAN.R-Project.org/package=longpower) (power and sample size for linear models of
        longitudinal data. Supported models include mixed-effects models
        and models fit by generalized least squares and generalized
        estimating equations)

    -   [powerlmm](https://CRAN.R-Project.org/package=powerlmm) (Calculate power for the 'time x treatment' effect
        in two- and three-level multilevel longitudinal studies with
        missing data. Studies with partially nested designs, unequal
        cluster sizes, unequal allocation to treatment arms, and
        different dropout patterns per treatment are supported. For all
        designs power can be calculated both analytically and via
        simulations)

    -   [powerMediation](https://CRAN.R-Project.org/package=powerMediation) (Power calculation for testing if mean
        changes for 2 groups are the same or not for longitudinal study
        with 2 or more time point.)

    -   [ssrm.logmer](https://CRAN.R-Project.org/package=ssrm.logmer) (sample size for a longitudinal study with
        binary outcome in order to attain a pre-specified power while
        strictly maintaining the Type I error rate.)

    -   [WebPower](https://CRAN.R-Project.org/package=WebPower) 

-   **MANOVA**

    -   [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl) one factor MANOVA with up to 2 levels and 4
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

    -   [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl) 

    -   [WebPower](https://CRAN.R-Project.org/package=WebPower) (Statistical Power Analysis for Structural Equation
        Modeling / Mediation based on Monte Carlo Simulation; and, for
        Simple Mediation)

-   **Meta-Analysis**

    -   [metapower](https://CRAN.R-Project.org/package=metapower) (Power for meta-analysis, including power analysis
        of main effects)

    -   [metameta](https://CRAN.R-Project.org/package=metameta) (Github package for re-analyzing published
        meta-analysis, including calculating power for each study in a
        meta-analysis.)

-   **Microarray Study**

    -   [FDRsampsize](https://CRAN.R-Project.org/package=FDRsampsize) (compute average power and sample size for
        microarray studies that use the false discovery rate as the
        final measure of statistical significance)

    -   [ssize.fdr](https://CRAN.R-Project.org/package=ssize.fdr) (This package contains a set of functions that
        calculates appropriate sample sizes for one-sample t-tests,
        two-sample t-tests, and F-tests for microarray experiments based
        on desired power while controlling for false discovery rates.)

-   **Micro-randomized trials**

    -   [MRTSampleSize](https://CRAN.R-Project.org/package=MRTSampleSize) 

    -   [MRTSampleSizeBinary](https://CRAN.R-Project.org/package=MRTSampleSizeBinary) (sample size calculator for
        micro-randomized trials (MRTs) with binary outcomes)

-   **Mixed Models**

    -   [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl) & [pass.lme](https://CRAN.R-Project.org/package=pass.lme) provide power and sample size
        calculation for the fixed effect term

    -   [pamm](https://CRAN.R-Project.org/package=pamm) (assess or explore the power of a dataset to estimates
        significant random effects (intercept or slope) in a mixed
        model. The functions are based on the \`\`lme4'' and
        \`\`lmerTest'' packages)

    -   Simulations of lmms and glmms can be performed with [simr](https://CRAN.R-Project.org/package=simr) and
        [simglm](https://CRAN.R-Project.org/package=simglm).

-   **Multiple Comparisons**

    -   [rPowerSampleSize](https://CRAN.R-Project.org/package=rPowerSampleSize) sample size for single-step (Bonferroni)
        and step-wise procedures (Holm and Hochberg) to control Type-II
        Generalized Family-Wise Error Rate

    -   [FDRsampsize](https://CRAN.R-Project.org/package=FDRsampsize) and [pwrFDR](https://CRAN.R-Project.org/package=pwrFDR) provides functions to compute
        power and sample sizes for studies using false discovery rate as
        the final measure of statistical significance.

-   **Multisite Randomized Trial**

-   **Multivariable Prediction Model**

    -   [pmsampsize](https://CRAN.R-Project.org/package=pmsampsize) computes the minimum sample size required for the
        development of a new multivariable prediction model. Can be used
        to calculate the minimum sample size for the development of
        models with continuous, binary or survival (time-to-event)
        outcomes

-   **Pearson Correlation**

    -   [WebPower](https://CRAN.R-Project.org/package=WebPower), [pwr](https://CRAN.R-Project.org/package=pwr) and [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl) all contain methods for power
        analysis of correlations.

-   **Pharmacokinetic Study Design:**

    -   [PharmPow](https://CRAN.R-Project.org/package=PharmPow): mixed (sparse/dense sampled) pharmacokinetic study
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
        regression

-   **Portfolio Efficiency**

    -   [GRS.test](https://CRAN.R-Project.org/package=GRS.test) implements power analysis for the GRS test for
        portfolio efficiency

-   **Proportion Test**

    -   [BAEssd](https://CRAN.R-Project.org/package=BAEssd) (Employes a Bayesian average error based approach to
        sample size determination. Several functions are included for
        sample size calculation for common designs in clinical trials
        including one- and two-sample binary and normal responses)

    -   [MKpower](https://CRAN.R-Project.org/package=MKpower) contains methods for one-sample proportion tests
        while [PASSED](https://CRAN.R-Project.org/package=PASSED) contains methods for two, while
        [PowerAnalysis](https://CRAN.R-Project.org/package=PowerAnalysis), [pwr](https://CRAN.R-Project.org/package=pwr), [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl), and [WebPower](https://CRAN.R-Project.org/package=WebPower) support both.

    -   [binomSamSize](https://CRAN.R-Project.org/package=binomSamSize) (compute confidence intervals and necessary
        sample sizes for the parameter *p* of the Bernoulli B(p)
        distribution under simple random sampling or under pooled
        sampling. Such computations are e.g. of interest when
        investigating the incidence or prevalence in populations.)

-   **Quantitative Comparative Analysis**

    -   [qcapower](https://CRAN.R-Project.org/package=qcapower) (Estimate power of a sufficient term using
        permutation tests. A term can be anything: A condition,
        conjunction or disjunction of any combination of these. The
        package further allows users to plot the estimation results and
        to estimate the number of cases required to achieve a certain
        level of power, given a prespecified null and alternative
        hypothesis.)

**Randomized Controlled Trial**

-   [CP](https://CRAN.R-Project.org/package=CP) (for condition power for different models in survival time
    analysis within RCT's with two different treatments to be compared
    and survival as an endpoint)

-   [odr](https://CRAN.R-Project.org/package=odr) (perform power analyses with and without accommodating cost
    structures of sampling for experimental studies under a budget
    constraint)

-   [powerCompRisk](https://CRAN.R-Project.org/package=powerCompRisk) (2-group comparisons cause-specific hazard and the
    all-cause hazard sample size calculations using an asymptotic
    chi-square joint test, accounting for censoring (e.g., lost to
    follow-up, staggered entry and administrative censoring))

-   [RCT](https://CRAN.R-Project.org/package=RCT) (Computes the minimum population needed to detect difference
    between control group and each treatment, given a target minimum
    detectable effect)

-   [Sample.Size](https://CRAN.R-Project.org/package=Sample.Size) (Computes the required sample size using the optimal
    designs with multiple constraints. This optimal method is designed
    for two-arm, randomized phase II clinical trials, and the required
    sample size can be optimized either using fixed or flexible
    randomization allocation ratios.)

-   [SampleSize4ClinicalTrials](https://CRAN.R-Project.org/package=SampleSize4ClinicalTrials) (calculate sample size when comparing
    means or proportions in Phase III clinical trials with different
    research goals)

-   [smartsizer](https://CRAN.R-Project.org/package=smartsizer) (Set of tools for determining the necessary sample
    size in order to identify the optimal dynamic treatment regime in a
    sequential, multiple assignment, randomized trial (SMART))

-   [ssev](https://CRAN.R-Project.org/package=ssev) (computes the optimal sample size for various 2-group
    designs (e.g., when comparing the means of two groups assuming equal
    variances, unequal variances, or comparing proportions) when the aim
    is to maximize the rewards over the full decision procedure of a)
    running a trial (with the computed sample size), and b) subsequently
    administering the winning treatment to the remaining N-n units in
    the population.)

-   [odr](https://CRAN.R-Project.org/package=odr) performs power analyses with and without accommodating cost
    structures of sampling for experimental studies under a budget
    constraint for multisite RCTs. [WebPower](https://CRAN.R-Project.org/package=WebPower) also supports multisite
    RCTs.

-   **Cluster Randomized Trials:**

    -   [clusterPower](https://CRAN.R-Project.org/package=clusterPower) (calculations for wide array of cluster design
        trials (1) Simple two-arm comparison designs, (2) Difference in
        difference designs, (3) Individually randomized group treatment
        trials, (4) Stepped wedge designs, (5) Multiarm trial designs)

    -   [CRTSize](https://CRAN.R-Project.org/package=CRTSize) (Sample size estimation in cluster (group)
        randomized trials)

    -   [odr](https://CRAN.R-Project.org/package=odr) (perform power analyses with and without accommodating
        cost structures of sampling for experimental studies under a
        budget constraint)

    -   [PowerUpR](https://CRAN.R-Project.org/package=PowerUpR) (calculate statistical power, minimum detectable
        effect size (MDES), MDES dif-ference (MDESD), and minimum
        required sample size for various multilevel randomized
        experi-ments with continuous outcomes. Some of the functions can
        assist with planning two- and three-level cluster-randomized
        trials (CRTs) sensitive to multilevel moderation and mediation
        (2-1-1, 2-2-1, and 3-2-1)

    -   [SteppedPower](https://CRAN.R-Project.org/package=SteppedPower) (Tools for power and sample size calculation as
        well as design diagnostics for longitudinal mixed model
        settings, with a focus on stepped wedge designs.)

    -   [swdpwr](https://CRAN.R-Project.org/package=swdpwr) (statistical power calculation for stepped wedge
        cluster randomized trials,)

    -   [WebPower](https://CRAN.R-Project.org/package=WebPower) 

-   **Adaptive Study Design**

    -   [Spass](https://CRAN.R-Project.org/package=Spass), [esDesign](https://CRAN.R-Project.org/package=esDesign) (sample size calculation where eligibility
        criteria of the trial is adaptively updated)

-   **Equivalence trial**

    -   [epiR](https://CRAN.R-Project.org/package=epiR) (sample size for parallel equivalence trial with binary
        or continuous outcome)

-   **Non-Inferiority Trial**

    -   [epiR](https://CRAN.R-Project.org/package=epiR) (binary or continuous outcome)

    -   [blindrecalc](https://CRAN.R-Project.org/package=blindrecalc) (computation of key characteristics and plots
        for blinded sample size recalculation, including powero.
        Continuous and binary endpoints are supported)

-   **Sequential Multiple Assignment Randomized Trial (SMART)**

    -   [SMARTbayesR](https://CRAN.R-Project.org/package=SMARTbayesR) (optimal dynamic treatment regimes and sample
        size for a SMART design in the Bayesian setting)

    -   [SMARTp](https://CRAN.R-Project.org/package=SMARTp) (Sample size calculation to detect dynamic treatment
        regime (DTR) effects based on change in clinical attachment
        level (CAL) outcomes from a non-surgical chronic periodontitis
        treatments study. The clustered tooth (sub-unit) level CAL
        outcomes are skewed, spatially-referenced, and non-randomly
        missing.)

-   **Superiority Trial**

    -   [blindrecalc](https://CRAN.R-Project.org/package=blindrecalc) (computation of key characteristics and plots
        for blinded sample size recalculation, including power.
        Continuous and binary endpoints are supported)

    -   [epiR](https://CRAN.R-Project.org/package=epiR) (sample size for parallel superiority trial with binary
        or continuous outcome)

    -   [MIDN](https://CRAN.R-Project.org/package=MIDN) (Nearly exact sample size calculation for exact
        powerful nonrandomized tests for differences between binomial
        proportions)

-   **Rasch Model Test**

    -   [pwrRasch](https://CRAN.R-Project.org/package=pwrRasch) (Statistical power simulation for testing the Rasch
        Model based on a three-way analysis of variance design with
        mixed classification)

-   **Regression (also see Mixed Models and Poisson Distribution)**

    -   [pwr2ppl](https://CRAN.R-Project.org/package=pwr2ppl), [pwr](https://CRAN.R-Project.org/package=pwr), and [WebPower](https://CRAN.R-Project.org/package=WebPower) all have general methods for
        various forms of regression. `pwr2ppl::pwr.f2.test` provides power for
        comparing dependent coefficients in multiple regression with two
        or three predictors.

    -   [simr](https://CRAN.R-Project.org/package=simr) and [simglm](https://CRAN.R-Project.org/package=simglm) both support simulation of regression
        models with various structures.

    -   [Hmisc](https://CRAN.R-Project.org/package=Hmisc) simulates power of a Bayesian odds ordinal logistic
        model.

    -   [PASSED](https://CRAN.R-Project.org/package=PASSED) & [MKpower](https://CRAN.R-Project.org/package=MKpower) for comparing two negative binomial
        rates

    -   [RSPS](https://CRAN.R-Project.org/package=RSPS) provides an empirical approach for RNA-seq data

-   **Regression Discontinuity (RD) Design**

    -   Both [rdpower](https://CRAN.R-Project.org/package=rdpower) and [rddapp](https://CRAN.R-Project.org/package=rddapp) provide tools to perform power,
        sample size and MDE calculations in RD designs.

-   **Repertory Grid Studies**

    -   [gridsampler](https://CRAN.R-Project.org/package=gridsampler) is a simulation tool to facilitate determination
        of required sample size to achieve category saturation for
        studies using multiple repertory grids in conjunction with
        content analysis.

-   **Sib Pair Design**

    -   

-   [powerpkg](https://CRAN.R-Project.org/package=powerpkg) computes power of an affected sib pair design. **Sign
    Test**

    -   [FDRsampsize](https://CRAN.R-Project.org/package=FDRsampsize) provides methods for microarray studies that use
        the false discovery rate as the final measure of statistical
        significance.

-   **Single Arm Trial that uses previous trials as historical control**

    -   [HCT](https://CRAN.R-Project.org/package=HCT) (Given a database of previous treatment/placebo
        estimates, their standard errors and sample sizes, the program
        calculates a significance criteria and power estimate that takes
        into account the among trial variation.)

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

    -   [Coprimary](https://CRAN.R-Project.org/package=Coprimary) provides sample size calculation for two primary
        time-to-event endpoints in clinical trials

    -   [powerMediation](https://CRAN.R-Project.org/package=powerMediation) calculates power for testing a mediation
        effect in cox regression based on Vittinghoff, Sen &
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

    -   [ssizeRNA](https://CRAN.R-Project.org/package=ssizeRNA) (appropriate sample sizes for two-sample t-test for
        RNA-seq experiments with fixed or varied set of parameters)

    -   [MKpower](https://CRAN.R-Project.org/package=MKpower) computes power for Welch and Hsu *t-*tests for
        balanced designs.

-   **TDT Study (transmission disequilibrium test)**

    -   [powerpkg](https://CRAN.R-Project.org/package=powerpkg) (Calculation of how many TDT trios need to be
        studied in order to have an adequately powered TDT study.)

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

    -   [samplesize](https://CRAN.R-Project.org/package=samplesize) (For categorical data. Allows for ties)

#### General simulation packages useful for power analysis

In some cases, there are either no analytical solutions (a knowledge
gap) or no programmed function of an analytical solution (an
implementation gap) in R. Empirical power (i.e., a simulation study of
the power) can be used in either of these cases, or when the assumptions
of the analytical solution are not met. General packages helpful for
simulation are shown below; method-specific packages are integrated into
the **<u>Specific models/methods/tests</u>** section above.

-   [DeclareDesign](https://CRAN.R-Project.org/package=DeclareDesign) provides a general framework for specifying research designs, simulated based on specified designs, and assessing various properties of the designs, including but not limited to power.

-   [simpr](https://github.com/simpr) is a GitHub R package that has a framework to data simulation—based on the tidyverse / broom syntax—and can perform power calculations

-   Various general-purpose MonteCarlo simulation packages exists, including the [simulator](https://CRAN.R-Project.org/package=simulator), [MonteCarlo](https://CRAN.R-Project.org/package=MonteCarlo), and others.

#### Related links

· Tutorial: [WebPower](https://webpower.psychstat.org/wiki/) is a forum
with a Stack Overflow like feature for asking questions about power
calculations in R. WebPower also has [extensive educational documents
for many tests](https://webpower.psychstat.org/wiki/models/index). Click
on a test name and then select ‘Parameter:help’.

· Tutorial: [useR! 2020: simpr: concise and readable simulations for the
tidyverse (E. Brown)](https://www.youtube.com/watch?v=MkfMSe9re2U)

· Tutorial: Precision calculation
<https://jonlefcheck.net/2015/03/31/how-much-is-enough-a-new-technique-for-quantifying-precision-of-community-surveys/>

· Tutorial: simulating power for a Longitudinal Markov Ordinal Outcome
Trials <https://hbiostat.org/R/Hmisc/markov/sim.pdf> and
<https://hbiostat.org/R/Hmisc/markov/funs.r>

· Tutorial: Power calculation for statistical interactions
<https://statmodeling.stat.columbia.edu/2018/03/15/need-16-times-sample-size-estimate-interaction-estimate-main-effect/>

· Tutorial: [Cross
Validation](https://stats.stackexchange.com/questions/tagged/statistical-power)
tags for statistical power

· Journal article: Lenth, R. V. (2001), “Some Practical Guidelines for
Effective Sample Size Determination,” *The American Statistician*, 55,
187-193.

· Journal article: Hoenig, John M. and Heisey, Dennis M. (2001), “The
Abuse of Power: The Pervasive Fallacy of Power Calculations for Data
Analysis,” *The American Statistician*, 55, 19-24.

· Book: Zhang, Z., & Yuan, K.-H. (2018). [Practical Statistical Power
Analysis Using Webpower and R](https://amzn.to/2V42nFp) (Eds). Granger,
IN: ISDSA Press.

· Book: Liu, X. S. (2013). [Statistical Power Analysis for the Social
and Behavioral Sciences: Basic and Advanced
Techniques](http://www.amazon.com/gp/product/1848729812/ref=as_li_qf_sp_asin_tl?ie=UTF8&camp=1789&creative=9325&creativeASIN=1848729812&linkCode=as2&tag=pyschstatorg-20).
Routledge.

Online Book: Blair, G., Coppock, A., & Humphreys, M. (in press).
Research Design: Declare, Diagnose, Redesign. Draft manuscript under
advance contract, Princeton University Press. Available at
<https://book.declaredesign.org/>

**RCT with multiple co-primary endpoints** (based on z-test):
mpe::power.known.var

<!-- <https://cran.r-project.org/web/packages/mpe/> -->

**TOSTER:** Power analysis for TOST (Two one-sided tests) for difference
between two proportions using Z-test (pooled)

-   Two one-sided tests procedure to test equivalence for t-tests,
    correlations, differences between proportions, and meta-analyses,
    including power analysis for t-tests and correlations. Allows you to
    specify equivalence bounds in raw scale units or in terms of effect
    sizes. See: Lakens (2017) .

**SteppedPower**: Power analysis for two sample z-test
