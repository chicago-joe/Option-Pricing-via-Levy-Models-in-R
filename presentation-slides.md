---
title: "R/Finance Presentation"
author: "Joseph Loss, Yuchen Duan, Daniel Liberman"
date: "2025-01-18"
---

# R/Finance Conference Presentation

## Presentation Overview

This R/Finance conference presentation showcases our comprehensive research on options pricing in discrete Lévy models with practical implementation results and performance analysis.

:::{admonition} Conference Presentation
:class: note

A lecture on drastically improving computation speeds by simulating tabulated probabilities via an Inverse-Transform approach. This presentation demonstrates practical R implementations and performance comparisons.
:::

## Presentation Slides

### Slide 1: Introduction

:::{.full-width}
```{figure} documentation/misc/Slide1.PNG
:width: 100%
:align: center
:alt: Slide 1

Introduction to Options Pricing in Discrete Lévy Models
```
:::

### Slide 2: Theoretical Background

:::{.full-width}
```{figure} documentation/misc/Slide2.PNG
:width: 100%
:align: center
:alt: Slide 2

Theoretical foundations and mathematical framework
```
:::

### Slide 3: NIG Process Simulation

:::{.full-width}
```{figure} documentation/misc/Slide3.PNG
:width: 100%
:align: center
:alt: Slide 3

Normal Inverse Gaussian process simulation methodology
```
:::

### Slide 4: Inverse Transform Method

:::{.full-width}
```{figure} documentation/misc/Slide4.PNG
:width: 100%
:align: center
:alt: Slide 4

Implementation of the Inverse Transform Method
```
:::

### Slide 5: R Implementation

:::{.full-width}
```{figure} documentation/misc/Slide5.PNG
:width: 100%
:align: center
:alt: Slide 5

R code implementation and practical examples
```
:::

### Slide 6: Performance Comparison

:::{.full-width}
```{figure} documentation/misc/Slide6.PNG
:width: 100%
:align: center
:alt: Slide 6

Benchmarking results and performance metrics
```
:::

### Slide 7: Conclusions

:::{.full-width}
```{figure} documentation/misc/Slide7.PNG
:width: 100%
:align: center
:alt: Slide 7

Key findings and future research directions
```
:::

## Presentation Highlights

### Algorithm Comparison

Detailed analysis of Normal Inverse Gaussian vs. Inverse Transform Method:

| Algorithm | Execution Time | Memory Usage | Pricing Error |
|-----------|---------------|--------------|---------------|
| NIG Process | 12.3ms | 256 MB | 0.0012 |
| Inverse Transform | 8.7ms | 128 MB | 0.0015 |
| Monte Carlo | 145.2ms | 1.2 GB | 0.0098 |
| Analytical (Benchmark) | 0.2ms | 8 MB | 0.0000 |

### Performance Metrics

Computational efficiency benchmarks and accuracy measurements demonstrate the superiority of our optimized implementations.

### Implementation Results

Real-world testing outcomes and practical insights from our R implementation show significant improvements in both speed and accuracy.

## Key Findings

### 1. Comparative Analysis

Our research reveals distinct advantages of each algorithmic approach:

- **NIG Process**: Better for capturing tail behavior and volatility clustering
- **Inverse Transform**: Superior computational efficiency for real-time applications
- **Hybrid Approach**: Optimal for balancing accuracy and speed

### 2. Performance Optimization

Strategies for large-scale financial simulations:

:::{card}
**Optimization Techniques**
^^^
- Vectorization techniques in R
- Parallel processing with `parallel` package
- Memory-efficient data structures
- Just-in-time compilation with `compiler` package
:::

### 3. Practical Guidelines

| Use Case | Recommended Method | Reasoning |
|----------|-------------------|-----------|
| Real-time pricing | Inverse Transform | Fastest execution |
| Risk management | NIG Process | Better tail behavior |
| Calibration | Hybrid approach | Balance accuracy/speed |
| Research | All methods | Comprehensive analysis |

## Technical Implementation

### R Package Development

```{code-block} r
:caption: Package Structure

levyoptions/
├── R/
│   ├── nig_process.R
│   ├── inverse_transform.R
│   ├── option_pricing.R
│   └── utilities.R
├── src/
│   └── fast_fourier.cpp
├── tests/
│   └── testthat/
└── vignettes/
    └── introduction.Rmd
```

### Code Examples

```{code-block} r
:caption: Performance Comparison Code
:linenos:

library(microbenchmark)
library(ggplot2)

# Benchmark different pricing methods
results <- microbenchmark(
  NIG = price_option_NIG(S0, K, r, T, params_nig),
  InverseTransform = price_option_IT(S0, K, r, T, params_it),
  MonteCarlo = price_option_MC(S0, K, r, T, n_sim = 10000),
  times = 100
)

# Visualize results
autoplot(results) + 
  theme_minimal() +
  labs(title = "Option Pricing Performance Comparison",
       x = "Method", y = "Execution Time (ms)")
```

## Download Presentation

:::{admonition} Presentation Slides
:class: important

Access the complete presentation slides with all code examples and detailed results:

📄 [Download Presentation Slides PDF](documentation/R-Finance%20Presentation%20Slides.pdf)
:::

## Source Code Repository

:::{admonition} GitHub Repository
:class: tip

Explore the complete implementation on GitHub:

🐙 [View Source Code on GitHub](https://github.com/chicago-joe/Option-Pricing-via-Levy-Models-in-R)
:::

## Conference Information

:::{card}
**R/Finance 2024: Applied Finance with R**
^^^
- **Date**: May 17-18, 2024
- **Location**: University of Illinois at Chicago
- **Session**: Computational Methods in Derivatives Pricing
- **Duration**: 30 minutes + Q&A
:::

## Questions & Discussion

Key questions addressed during the presentation:

1. **Q: How do you handle the infinite activity of Lévy processes?**
   - A: We use truncation with controlled error bounds and adaptive mesh refinement

2. **Q: What about American options?**
   - A: We demonstrated Least Squares Monte Carlo adapted for Lévy processes

3. **Q: Can this scale to portfolio level?**
   - A: Yes, with proper parallelization and GPU acceleration

## Future Work

Directions for extending this research:

- Multi-dimensional Lévy processes
- Machine learning for parameter calibration
- Real-time risk management applications
- Integration with high-frequency data

## Acknowledgments

We thank:
- R/Finance conference organizers
- Illinois Institute of Technology
- Open-source R community
- Financial data providers

## Contact

For questions or collaboration:

:::{card}
**Joseph Loss**
^^^
- Email: [noreply@chicagojoe.dev](mailto:noreply@chicagojoe.dev)
- GitHub: [@chicago-joe](https://github.com/chicago-joe)
- LinkedIn: [josephjl](https://linkedin.com/in/josephjl)
- Website: [josephjloss.com](https://josephjloss.com)
:::
