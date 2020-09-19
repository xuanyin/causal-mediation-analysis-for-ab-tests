# Causal Mediation Analysis for A/B Tests
## KDD 2019: Accepted and 20-Min Oral Presentation (6.43% selection rate, 45 out of 700)

## Citation
Xuan Yin and Liangjie Hong. 2019. [The Identification and Estimation of Direct and Indirect Effects in A/B Tests through Causal Mediation Analysis](https://www.kdd.org/kdd2019/accepted-papers/view/the-identification-and-estimation-of-direct-and-indirect-effects-in-online-). In *The 25th ACM SIGKDD Conference on Knowledge Discovery and Data Mining (KDD '19), August 4-8, 2019, Anchorage, AK, USA.* ACM, New York, NY, 11 pages. [https://doi.org/10.1145/3292500.3330769](https://doi.org/10.1145/3292500.3330769)

## <ins>[The Published Version of the Paper](https://www.kdd.org/kdd2019/accepted-papers/view/the-identification-and-estimation-of-direct-and-indirect-effects-in-online-)</ins>
## <ins>[Code-as-Craft Post: The Causal Analysis of Cannibalization in Online Products](https://codeascraft.com/2020/02/24/the-causal-analysis-of-cannibalization-in-online-products/)
## <ins>[LinkedIn Post: The Causal Analysis of Cannibalization in Online Products](https://www.linkedin.com/pulse/quantifying-business-impacts-cannibalization-ml-xuan-yin-ph-d-/)
## <ins>[Presentation Slides](KDD2019_short_version_Slides_paper_ads1688o_Causal_Mediation_Analysis.pdf)</ins>
## <ins>[Promotional Video on Youtube](https://youtu.be/coEpqU9HWWM)</ins>
## <ins>[Poster](KDD2019_Poster_paper_ads1688o_Causal_Mediation_Analysis.pdf)</ins>

## Implementation
The R function [cma.R](https://github.com/xuanyin/causal-mediation-analysis-for-ab-tests/blob/master/cma.R) implements causal mediation analysis via generalized method of moments.

## Abstract
E-commerce companies have a number of online products, such as organic search, sponsored search, and recommendation modules, to fulfill customer needs. Although each of these products provides a unique opportunity for users to interact with a portion of the overall inventory, they are all similar channels for users and compete for limited time and monetary budgets of users. To optimize users' overall experiences on an E-commerce platform, instead of understanding and improving different products separately, it is important to gain insights into the evidence that a change in one product would induce users to change their behaviors in others, which may be due to the fact that these products are functionally similar. In this paper, we introduce causal mediation analysis as a formal statistical tool to reveal the underlying causal mechanisms. Existing literature provides little guidance on cases where multiple unmeasured causally-dependent mediators exist, which are common in A/B tests.  We seek a novel approach to identify in those scenarios direct and indirect effects of the treatment. In the end, we demonstrate the effectiveness of the proposed method in data from Etsy's real A/B tests and shed lights on complex relationships between different products.

## Keywords
A/B test; causal inference; mediation analysis; potential outcome; experiment; online product
