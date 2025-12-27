# Chou–Fasman Protein Secondary Structure Predictor

This repository contains a Python implementation of the **Chou–Fasman method** for predicting
protein secondary structure elements —
**α-helices (H)** and **β-strands (S)** — from a raw amino acid sequence.

The implementation follows the classical Chou–Fasman workflow, including:

* **Seed (nucleation) identification**
* **Bidirectional region growth**
* **Resolution of overlapping helix and strand assignments**

The program reports both **predicted structural regions** and a **final residue-wise secondary structure annotation**.

---

## Key Features

* Predicts the following secondary structure types:
  * **Helix (H)**
  * **Strand (S)**
  * **Unstructured / Coil (-)**
* Fully implements Chou–Fasman decision rules:
  * Fixed-size nucleation windows for helices and strands
  * Sliding-window extension in both N- and C-terminal directions
  * Overlap handling based on average propensity comparison
* Produces detailed output including:
  * Start and end indices of predicted regions
  * Length and amino acid content of each region
  * Final structure string aligned with the input sequence  
    (e.g., `HHHHSSS---HHHSS...`)

---

## Method Overview

| Phase                     | Description                                                         |
| ------------------------- | ------------------------------------------------------------------- |
| **Helix Nucleation**         | Examine **6-residue windows**; valid if **≥4 residues have Pa > 1.0** |
| **Helix Extension**          | Extend while **ΣPa over 4 residues ≥ 4.0**                           |
| **Strand Nucleation**        | Examine **5-residue windows**; valid if **≥3 residues have Pb > 1.0** |
| **Strand Extension**         | Extend while **ΣPb over 4 residues > 4.0**                           |
| **Conflict Resolution**      | Resolve conflicts by comparing **average Pa vs Pb**                 |

---


## Running the Python Program

Run normally:
```bash
python3 chou_fasman_predictor.py
```

---

## Author

**Mritunjay Kumar**
B.Tech - CSB – IIIT Delhi

---
