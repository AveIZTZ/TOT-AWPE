# TOT-AWPE: A third-order tensor decomposition based algorithm for speech dereverberation

## 📖 Paper

> **A third-order tensor decomposition based algorithm for speech dereverberation**  
> *Signal Processing*, 239 (2026): 110210  

🔗 https://doi.org/10.1016/j.sigpro.2025.110210

---

## 🔍 Overview

Speech dereverberation is a fundamental technique for improving speech quality and intelligibility in reverberant acoustic environments.  
Among existing approaches, **adaptive weighted prediction error (AWPE)** methods based on multichannel linear prediction have shown strong performance, but they typically suffer from **high computational complexity** due to the long prediction filters involved.

This work proposes **TOT-AWPE**, a **third-order tensor decomposition based adaptive WPE algorithm**, which:

- Decomposes the global multichannel prediction filter using a **third-order Kronecker product structure**
- Replaces the estimation of one long filter with **multiple short, low-rank filters**
- Achieves **significant computational complexity reduction**
- Maintains **comparable or superior dereverberation performance** compared to AWPE and KAWPE

---

## ✨ Key Contributions

- **Third-order tensor decomposition for WPE**  
  Extends Kronecker-product-based linear prediction from second-order to third-order tensor decomposition.

- **Low-rank adaptive filter estimation**  
  Introduces a structured low-rank formulation that enables efficient recursive least-squares (RLS) adaptation.

- **Reduced computational complexity**  
  Lowers the complexity from  
  \(\mathcal{O}(L^2)\) (AWPE)  
  to  
  \(\mathcal{O}(Q^2 L_2^2 (L_{11}^2 + L_{12}^2))\) (TOT-AWPE).

- **Robustness to array geometry**  
  Demonstrates effectiveness for both uniform linear arrays (ULA) and uniform circular arrays (UCA).

---

## 🧩 Method: TOT-AWPE Framework

The key idea of **TOT-AWPE** is to represent the long multichannel linear prediction filter using a **third-order Kronecker product decomposition**:

\[
\mathbf{g}(t) =
\sum_{l=1}^{L_2} \sum_{q=1}^{Q}
\mathbf{g}_{2,l}(t)
\otimes
\mathbf{g}_{12,lq}(t)
\otimes
\mathbf{g}_{11,lq}(t)
\]

This formulation allows the original high-dimensional filter to be estimated through **three groups of much shorter filters**, which are updated iteratively using RLS-type recursions.

---

## ⚙️ Algorithm

The complete algorithm alternates between estimating:

- \(\mathbf{g}_2(t)\)
- \(\mathbf{g}_{12}(t)\)
- \(\mathbf{g}_{11}(t)\)

and updating the dereverberated speech estimate accordingly.

The full derivation and the detailed procedure are provided in **Algorithm 1** of the paper.

---

## 📊 Computational Complexity

| Method  | Computational Complexity |
|--------|--------------------------|
| AWPE   | \(\mathcal{O}(L^2)\) |
| KAWPE | Reduced compared to AWPE |
| **TOT-AWPE** | \(\mathcal{O}(Q^2 L_2^2 (L_{11}^2 + L_{12}^2))\) |

When \(Q \ll L\) and appropriate decomposition parameters are chosen, **TOT-AWPE achieves substantial computational savings**.

---

## 🧪 Experimental Results

Experiments were conducted using:

- 8-channel microphone arrays
- Reverberation times \(T_{60} \in \{400, 500, 600, 700\}\) ms
- Objective metrics:
  - Cepstral Distance (CD)
  - PESQ
  - Frequency-weighted segmental SNR (FWSNR)
  - Log-likelihood ratio (LLR)

### Key Findings

- TOT-AWPE consistently improves all objective metrics
- Performance improves as the tensor rank parameter \(Q\) increases
- With moderate \(Q\), TOT-AWPE outperforms AWPE and KAWPE
- Strong robustness across different array geometries (ULA and UCA)

---

## 📁 Code Availability

🚧 **Code will be released upon request or in future updates.**  
This repository currently serves as a **research reference implementation** accompanying the paper.

---

## 📌 Citation

If you find this work useful, please cite:

```bibtex
@article{zhu2026totawpe,
  title={A third-order tensor decomposition based algorithm for speech dereverberation},
  author={Zhu, Yujie and Huang, Gongping and Jin, Jilu and Chen, Jingdong and Benesty, Jacob},
  journal={Signal Processing},
  volume={239},
  pages={110210},
  year={2026},
  publisher={Elsevier}
}
