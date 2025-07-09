---
title: "Options Pricing in Discrete Lévy Models"
subtitle: "Research Implementation and Analysis"
authors:
  - Joseph Loss
  - Yuchen Duan  
  - Daniel Liberman
date: 2025-01-18
abstract: |
  This research explores the implementation of two algorithms for options pricing in Lévy models, 
  specifically comparing the Normal Inverse Gaussian (NIG) process with the Inverse Transform Method. 
  Our implementation builds upon the foundational work by Feng et al. on simulating Lévy processes 
  from their characteristic functions.
---

# Options Pricing in Discrete Lévy Models

## Research Overview

This research proposal outlines our comprehensive study of options pricing methodologies within Lévy models, specifically focusing on the implementation and comparison of two distinct algorithmic approaches.

```{figure} assets/levy-process-diagram.png
:width: 600px
:align: center
:alt: Lévy Process Diagram

Illustration of Lévy process paths and their application in options pricing
```

## Key Research Components

### Normal Inverse Gaussian (NIG) Process
Implementation of advanced stochastic processes for more accurate market modeling. The NIG process captures important stylized facts of financial markets including:
- Heavy tails in return distributions
- Volatility clustering
- Asymmetric returns

### Inverse Transform Method
Comparative analysis of computational efficiency and pricing accuracy through:
- Direct simulation from characteristic functions
- Optimized numerical integration techniques
- Performance benchmarking against traditional methods

### Theoretical Foundation
Building upon the seminal work by Feng et al. on simulating Lévy processes from characteristic functions, we extend their methodology to:
- Discrete-time option pricing models
- Efficient calibration procedures
- Risk-neutral measure transformations

### Practical Applications
Real-world options pricing scenarios and performance benchmarking including:
- European option pricing
- American option approximations
- Exotic derivatives valuation
- Greeks calculation and hedging strategies

## Research Objectives

1. **Develop robust implementations** of both pricing algorithms in R
2. **Conduct comprehensive performance analysis** comparing computational efficiency
3. **Validate pricing accuracy** against market data and theoretical benchmarks
4. **Provide practical guidance** for practitioners in quantitative finance

## Methodology

Our research methodology combines theoretical analysis with empirical validation:

```{mermaid}
graph TD
    A[Theoretical Framework] --> B[Algorithm Implementation]
    B --> C[Performance Testing]
    C --> D[Market Validation]
    D --> E[Results Analysis]
    E --> F[Practical Guidelines]
```

## Expected Contributions

This research aims to contribute to the quantitative finance literature by:

- Providing open-source R implementations of advanced Lévy process pricing models
- Establishing performance benchmarks for different algorithmic approaches
- Offering practical insights for model selection in real-world applications
- Creating educational resources for students and practitioners

## Research Timeline

| Phase | Duration | Deliverables |
|-------|----------|--------------|
| Literature Review | 2 months | Comprehensive survey of Lévy models |
| Implementation | 3 months | R package with core algorithms |
| Testing & Validation | 2 months | Performance benchmarks and accuracy metrics |
| Documentation | 1 month | Research paper and user guides |

## Download Research Proposal

Access the complete research proposal with detailed methodology, mathematical formulations, and expected outcomes:

```{button-link} documentation/Project Research Proposal.pdf
:color: primary
:align: center
📄 Download Full Research Proposal PDF
```

## Collaborators

This research is conducted at the Illinois Institute of Technology by:

- **Joseph Loss** - Lead Researcher, Algorithm Implementation
- **Yuchen Duan** - Mathematical Modeling, Theoretical Analysis
- **Daniel Liberman** - Empirical Validation, Market Data Analysis

## References

1. Feng, L., Linetsky, V., & Morales, J. L. (2016). *Options Pricing in Lévy Models*. Academic Paper.
2. Barndorff-Nielsen, O. E. (1997). *Normal inverse Gaussian distributions and stochastic volatility modelling*. Scandinavian Journal of Statistics, 24(1), 1-13.
3. Cont, R., & Tankov, P. (2004). *Financial modelling with jump processes*. Chapman and Hall/CRC.
