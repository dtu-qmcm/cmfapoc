---
title: Compositional Metabolic Flux Analysis
bibliography: refs.bib
graphics: yes
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

# Note to Christian
Metabolism refers to all chemical reactions that occur within living cells, allowing organisms to grow, reproduce, and respond to their environment. These reactions form interconnected networks called metabolic pathways, where each step converts one molecule to another through enzymatic reactions. Just as understanding traffic flow in a city requires knowing not just the road layout but also how many vehicles travel each route, understanding cellular function requires knowing both the structure of metabolic pathways and the rates at which molecules flow through them. These rates are called "metabolic fluxes." While measuring the concentrations of molecules (metabolites) in a cell is relatively straightforward with modern techniques, determining how quickly these molecules are being produced, converted, or consumed inside the cell remains challenging. Scientists have developed approaches using isotopes—variants of elements with different atomic weights—as tracers to track the flow of atoms through metabolic pathways. These isotopes act like molecular GPS trackers, allowing researchers to follow the journey of specific atoms through the complex network of cellular metabolism. The following section explains how these isotope-based approaches work and introduces a new perspective on analyzing the resulting data.

# Introduction

The rates, or "fluxes", of metabolic reactions are key to understanding how cells behave. The rates of internal metabolic reactions, which create and destroy chemicals inside of cells, are particularly interesting, but difficult to measure directly. Isotope labelling experiments make it possible to infer these internal fluxes indirectly, thanks to the fact that metabolic reactions split and join molecules in predictable ways.

Specifically, given a known "atom mapping" describing how each reaction in a network moves potentially-labelled atoms, a specification of fluxes, and some known isotope distributions (typically the feed), one can calculate expected isotope distributions for internal metabolites. These internal isotope distributions can also be measured, for example using liquid chromatography and mass spectrometry. Comparing measured and expected isotope distributions makes it possible to say infer flux specifications.

This technique has existed in its modern form for a while, leading to impressive results, but we think that a crucial step is missing. Measured and expected isotope distributions are compositions, but they have not previously been analysed using compositional methods. 

## Isotopes, isotopomers and isotopologues

Isotopes are atoms whose nuclei have the same number of protons but different numbers of neutrons. Isotopes of the same element have similar chemical properties, but have different atomic masses and physical properties. For example, carbon has three naturally occurring isotopes: 12C, 13C and 14C, with respective atomic masses 12, 13 and 14. 14C occurs in negligible quantities, and the natural ratio of 12C to 13C is known, making carbon suitable for isotope labelling experiments where 12C is artificially replaced by 13C.

Isotopomers are forms of a compound that differ only by substitution of isotopes. For example, [1-13C] glucose and [2-13C] glucose are isotopomers of glucose that differ only in the isotopes of the carbon atoms in positions 1 and 2 (see [figure below](#glucose)). In general, for every $a$ occurrences of an atom with $I$ isotopes in a compound, there are $I^a$ corresponding isotopomers. For example, glucose has six carbon atoms: assuming only 12C and 13C isotopes are present, there are $2^6=64$ carbon isotopomers.

\begin{figure}[!h]
\centering
\includegraphics[width=0.4\textwidth, height=!]{img/1-13C-glucose.png}
\includegraphics[width=0.4\textwidth, height=!]{img/2-13C-glucose.png}
\caption{Left: [1-13C] glucose Right: [2-13C] glucose}
\label{glucose}
\end{figure}

Isotopologues are forms of a compound that differ in atomic mass. For example, the isotopomers [1-13C] glucose and [2-13C] glucose each have five 12C atoms and one 13C atom and therefore belong to the same isotopologue $M_1$ with atomic mass 181.15 g/mol. Isotopologues are important because measurements can often distinguish between isotopologues, but not between isotopomers with the same atomic mass.

## Isotope labelling based Metabolic Flux Analysis

In more detail, MFA considers a known metabolic network consisting of $M$ compounds and $N$ reactions with stoichiometric coefficients $S\in\mathbb{R}^{M\times N}$ representing the amount of each compound consumed and produced by each reaction, plus an atom transition map for each reaction. The atom transition map for a reaction specifies in what order the potentially-labelled atoms occur in each of the reaction's substrates and products.

In this study we focus on the most straightforward case, where it is safe to assume that the network is in isotopic and metabolic steady state. This means that the concentrations of internal metabolites and isotopomers are constant.

The remaining input for MFA is as follows:

- Known isotopomer proportions for some compounds, typically external ones like the feed.

- Measured fluxes for some reactions.

- Measured isotopologue proportions for some compounds, possibly with known
measurement error.

The "forward problem" of inferring the label pattern given a flux specification can be solved by writing down a balance equation describing the rate of change of each isotopomer in the network. These equations can be found by combining the atom map and stoichiometric matrix to produce a matrix $S_i\in\mathbb{R}^{I\times N}$ of stoichiometric coefficients for isotopomers. At isotopic and metabolic steady state we have $Sv=0$, giving $N$ balance equations. Given some known isotopomer proportions, other isotopomer proportions can be calculated by solving these equations.

While conceptually simple, this approach to solving the forward problem is typically unfeasible due to the prohibitively large number of isotopomers. Practical applications therefore require a method to solve the forward problem more simply: see [@daiUnderstandingMetabolismFlux2017, §3.1] for an overview of approaches to this problem. The most widely used of these is the "elementary metabolite unit" representation introduced in [@antoniewiczElementaryMetaboliteUnits2007].

The inverse problem of inferring steady state fluxes from measured isotopomer or isotopologue distributions can be solved using a statistical model that links these measurements with latent parameters representing flux configurations. In general, such a model specifies the probability $p(r_{obs}\mid r(v))$ of observing labelling pattern $r_{obs}$ given a true flux assignment $v$ and true labelling pattern $r(v)$. For example, assuming an independent linear model with measurement standard deviation vector $\sigma$ yields the following relationship:

$$
p(r_{obs}\mid r(v)) = pdf_{normal}(r_{obs} | r(v), \sigma)
\label{linear}
$$

Sometimes a statistical model is not mentioned explicitly, but $v$ is optimised by minimising a scaled Euclidean distance from $r(v)$ to $r_{obs}$, which produces the same result as optimising by maximising the likelihood [$p(r_{obs}\mid r(v))$](#linear).

## Previous work

This section briefly reviews previous work in labelling-based metabolic flux analysis. For more detailed review papers see [@daiUnderstandingMetabolismFlux2017], [@falcoMetabolicFluxAnalysis2022].

### Experimental methods

See the following references:

- [@longHighresolution13CMetabolic2019] 
- [@falcoMetabolicFluxAnalysis2022]

### Bayesian methods

See the following references:

- [@theorellBeCertainUncertainty2017]
- [@theorellReversibleJumpMCMC2020]

### Software

Existing implementations of MFA include:

- INCA
- 13CFLUX2
- Metran
- OpenFlux(2)
- FluxPyt
- mfapy
- Sysmetab
- iso2flux
- Flux-P
- WUFlux
- OpenMebius
- influx\_s

See [REFERENCE] for reviews of available software implementing MFA.

There is no previous implementation of compositional regression analysis in the context of MFA; all previous implementations apply a linear model either explicitly as in [@theorellBeCertainUncertainty2017, Eq. 3] or more commonly implicitly through the use of Euclidean distance based optimisation.

# Problem statement

Discrepancies between observed and expected isotopologue distributions have previously been modelled using linear regression models, despite the fact that both are compositions. It is well known that, in general, applying non-compositional data analysis methods to compositional data is dangerous because these methods can easily misinterpret constraint-induced correlations [@aitchisonjStatisticalAnalysisCompositional, Ch. 3]. This issue is especially pronounced in the case where there are relatively few composition components.

We therefore propose modelling augmenting existing labelling-based MFA methods with specialised compositional data analysis. In this section we draw on previous literature where compositional data analysis to discuss common arguments against the use of specialised compositional methods, explaining why we think they don't apply in the case of labelling based MFA.

## Arguments against using compositional methods

# Simple pedagogical example

# Realistic example
