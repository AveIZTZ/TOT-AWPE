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
