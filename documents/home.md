### Welcome to Randk\uft...

**Randkluft** is a fast and intuitive tool for identifying **rare, true-positive cells** with high marker intensity in multiplexed imaging data.

In many markers, true positive signals are obscured by a large population of true negative cells with low intensity. This overlap makes it difficult to define a reliable, objective gating threshold using conventional approaches.

Randkluft addresses this challenge by analyzing the **density profile of each marker**. Using a stochastic search strategy inspired by the **Robbins–Monro algorithm**, it identifies a discriminative threshold that is often hidden within the **right-hand shoulder of the distribution**, where rare positive cells reside.

By combining automated estimation with visual inspection and manual refinement, Randkluft enables robust, reproducible gating of marker expression across samples.