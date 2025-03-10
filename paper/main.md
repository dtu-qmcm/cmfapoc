---
title: Compositional Metabolic Flux Analysis
bibliography: refs.bib
author:
  - name: Teddy Groves
    affil: [1]
    email: "tedgro@dtu.dk"
    orcid: asdfasdef
  - name: Te Chen
    affil: [1]
    orcid: asdfasdef
  - name: Sergi Muyo Abad
    affil: [1]
    orcid: asdfasdef
  - name: Nicholas Luke Cowie
    affil: [1]
    orcid: asdfasdef
  - name: Daria Volkova
    affil: [1]
    orcid: asdfasdef
  - name: Christian Brinch
    affil: [2]
    orcid: asdfasdef
  - name: Lars Keld Nielsen
    affil: [1, 3]
    orcid: asdfasdef
affiliation:
  - num: 1
    org: The Novo Nordisk Center for Biosustainability, DTU, Kongens Lyngby, Denmark
  - num: 2
    org: National Food Institute, DTU, Kongens Lyngby, Denmark
  - num: 3
    org: Australian Institute for Bioengineering and Nanotechnology (AIBN), The University of Queensland, St Lucia 4067, Australia
abstract: |
  Metabolic Flux Analysis aims to infer the values of metabolic fluxes from measurements of isotope labelling distributions. Since these distributions are positive, sum-constrained and relatively low-dimensional, we argue that they should be analysed using specialised methods that target compositional data. We illustrate our argument using a simple pedagogical example, then show how compositional analysis leads to improved results on a typical dataset.
---

<!-- email: "techen@dtu.dk" -->
<!-- email: "lakeni@dtu.dk" -->
<!-- email: "nicow@dtu.dk" -->
<!-- email: "dariav@dtu.dk" -->
<!-- email: "cbri@dtu.dk" -->

# Introduction

Metabolic flux analysis is the study of

## Previous work

This section briefly reviews previous work in labelling-based metabolic flux analysis. For more detailed review papers see XXXXX

### Experimental methods
### The forward problem
- Cumomers
- EMU

### Software
- OpenFlux
- Freeflux
- INCA
- 13CFlux
- ...

### Bayesian 13C MFA

## Problem statement

The topic of how to statistically model isotope labelling pattern measurements has received relatively little attention in the development of metabolic flux analysis. Most presentations and software applications quantify the discrepancy between a species's measured and predicted isotope labelling distribution using a Euclidean distance, and advocate choosing a flux configuration whose labelling pattern minimises this distance, possibly with per-species and/or per-isotope-equivalence-class weights. This is equivalent to using maximum likelihood estimation, where the likelihood is given by an indepedent normal distribution centered on the predictedc labelling distribution, with error standard deviations determined by the weights, i.e., for each species s,

<!-- $$ -->
<!-- y_s \sim N(\hat{y_s}, \sigma_s) -->
<!-- $$ -->

<!-- where $y_s$ is the observed labelling distribution for species $s$, $\hat{y_s}$ is the predicted labelling distribution and $sigma_s$ is a vector of standard deviations. -->

<!-- There are two key reasons why this approach is flawed in the case where $y_s$ and $\hat{y_s}$ are compositions. First, the Euclidean distance is inappropriate for measuring discrepancies between compositions. Second, the use of an independent error model neglects the fact that composition components are intrinsically correlated. This issue is especially pronounced in the case where there are relatively few composition components. -->


