# Research and Industry Applications of Molecular Simulations

This lecture examines four recent examples where molecular simulation has moved from academic methodology to engineering outcome. Together they illustrate a common thread: computation is no longer just a tool for understanding materials after the fact — it is becoming a tool for *designing* them before synthesis begins.

---

## 1. The Shift from Trial-and-Error to Computation-First

Thomas Edison tested more than 2,000 materials before finding a viable filament for the lightbulb. For most of the twentieth century, materials science worked the same way: propose a composition, synthesize it, test it, and repeat. This cycle is expensive and slow — the development of commercial lithium-ion batteries, for example, spanned nearly two decades from laboratory discovery to consumer product.

The emergence of high-performance computing and atomistic simulation changes this equation. Instead of asking "what happens when I make this material?", we can now ask "which of these 100,000 candidate materials is worth making?" Computation cannot replace synthesis, but it can dramatically narrow the search space before a single gram of material is ever processed.

The four case studies in this lecture represent this shift across different application domains: energy storage, polymer design, polymer recycling, and universal model development.

---

## 2. Case Study 1: How Computation Found a Better Battery Cathode

### The Problem

A major US battery manufacturer approached researchers at MIT with a challenge: find a better cathode material for alkaline (non-rechargeable) batteries. They provided $1 million, one year, and access to their supercomputer. The existing cathode — manganese dioxide — had been in use for decades.

The traditional approach would have been to synthesize and test candidate materials in the laboratory. With the number of possible inorganic compounds running into the hundreds of thousands, this was not feasible.

### The Computational Approach

Prof. Kristin Persson and her collaborators ran high-throughput density functional theory (DFT) calculations across the entire known space of inorganic compounds. Their initial screening pool contained approximately 130,000 potential compositions. By applying successive filters — voltage, energy capacity, energy density, and finally thermodynamic stability — they reduced this to around 200 compounds that met all the key performance criteria.

Those 200 computationally identified candidates were handed to the manufacturer. The company selected one, validated it experimentally, and in 2019 launched a new battery product with approximately twice the energy of the previous design.

### The Broader Impact

Persson's group and other collaborators have been developing high-throughput computational databases for a long time. In 2010, she proposed the Materials Genome Browser — a tool for scientists to search for materials by target properties — through a competition at DOE's Lawrence Berkeley National Laboratory. This proposal, combined with a follow-on Laboratory Directed Research and Development project focused on data mining for materials design, became the foundation for the Materials Project (materialsproject.org). The Materials Project officially launched in 2011 as one of the five inaugural software centers under DOE's Materials Genome Initiative. It is now one of the largest open-access databases of computed materials properties in the world, allowing researchers to filter across hundreds of thousands of known and predicted inorganic compounds using just a few lines of Python code.

The Materials Project uses **AiiDA** and similar workflow managers to run thousands of DFT calculations automatically, store every result with full provenance, and expose the data through a Python API (`pymatgen`). This is precisely the kind of high-throughput workflow discussed in Lecture 12.

---

## 3. Case Study 2: High-Throughput MD for Polymer Electrolyte Discovery

### The Problem

Solid polymer electrolytes could replace liquid electrolytes in lithium-ion batteries, eliminating the risk of flammable electrolyte leakage and enabling higher energy densities. However, lithium-ion conductivity in solid polymers at room temperature is too low for practical use. Finding a polymer with both high conductivity and good mechanical properties has resisted decades of experimental effort.

The challenge is the size of the design space. Polymer chemistry is combinatorially vast — small changes to monomer chemistry, chain architecture, or functional group placement can dramatically affect ion transport. No experiment can efficiently explore this space.

### The TRI HTP-MD Platform

Researchers at the Toyota Research Institute (TRI) and MIT built a closed-loop polymer discovery platform called **HTP-MD** (High-Throughput Molecular Dynamics). The workflow has three coupled components:

**Generative AI:** Models trained on known polymer structures — using transformer architectures similar to those behind large language models — propose new polymer candidates by learning the "syntax" of molecular structure. Candidates are checked for chemical validity, synthetic accessibility, and novelty before being passed forward.

**Automated MD simulation:** Each candidate polymer is automatically run through a molecular dynamics pipeline that computes ion transport properties — diffusion coefficients and conductivity — at multiple temperatures. The pipeline handles system building, force field assignment, equilibration, production, and analysis without human intervention.

**Experimental validation:** The most promising simulation-identified candidates are synthesized at MIT using a high-throughput experimental characterization platform that measures transport properties across temperatures and electrolyte concentrations.

The three components operate in a closed loop: experimental results feed back into the generative model, improving its proposals over successive iterations.

### Why MD Is Central

The key insight is that MD simulations can compute the same property (ion conductivity) that experiment would measure, but for compounds that have never been synthesized. Even if individual simulation predictions are imperfect, the *relative ranking* of candidates is reliable enough to guide synthesis toward the most promising materials. MD acts as a filter that costs compute time rather than lab time.

The platform has been published as open-source code (`htp_md` on GitHub) alongside an open database of computed polymer properties, enabling other research groups to build on the same infrastructure.

---

## 4. Case Study 3: Simulating Spatiotemporal Heating for Plastic Recycling

### The Problem

The global plastics waste crisis is partly a chemistry problem. Most commodity plastics are thermosets or semicrystalline polymers that cannot be efficiently recycled by conventional thermal processing — slow, uniform heating causes random chain scission and carbonization before useful monomers can be recovered.

Researchers including Prof. Lele and collaborators at the University of Maryland, Purdue University, and Princeton University developed an approach based on **electrified spatiotemporal heating**: passing electrical current through a thin conductive layer in contact with the polymer to generate extremely rapid, localized temperature spikes (Joule heating). The key discovery was that this spatiotemporal heating profile — heating rates of ~10⁵ K/s, temperatures exceeding 1000 K for milliseconds — could selectively break the polymer backbone into monomer units with high yield, before slower side reactions degraded the products.

### The Role of Simulation

Molecular simulation was essential for understanding the mechanism. Classical MD and reactive MD (ReaxFF) simulations were used to:

- Characterize the spatial temperature gradients at the polymer–conductor interface during Joule heating
- Track the kinetics of bond scission events as a function of the local temperature history
- Distinguish between depolymerization pathways that yield useful monomers versus those that produce char or small-molecule byproducts
- Understand how heating rate and duration — parameters that can be controlled electrically — influence product selectivity

The simulations provided atomic-level insight that could not be obtained from bulk temperature measurements. Specifically, they revealed that the extreme heating rate kinetically suppresses the unwanted side reactions that dominate under slow conventional heating — the system moves through the temperature regime where depolymerization dominates before it has time to undergo carbonization.

The results were published in *Nature* in 2023 and demonstrated selective conversion of polyethylene, polypropylene, and mixed plastic waste streams into recoverable hydrocarbon products.

### Broader Significance

This example illustrates a use case for reactive potentials (ReaxFF) that was discussed in Lecture 9: situations where bond breaking and forming are central to the phenomenon of interest and cannot be captured by classical force fields. It also shows that simulation adds value not only in the materials discovery phase but also in the mechanistic interpretation of novel processes — understanding *why* something works, not just *that* it works.

---

## 5. Case Study 4: Universal ML Potentials at Scale

The three case studies above each required a specific simulation setup. The Duracell work used DFT. The TRI polymer work used classical force fields parameterized for specific polymer chemistries. The depolymerization work used ReaxFF. Each had a setup cost — force field selection, parameterization, validation — that limited how broadly the approach could be deployed.

The most significant recent development in the field is the emergence of **universal ML potentials**: single models trained on sufficiently diverse data that they can be applied to a wide range of chemistries without system-specific parameterization.

### MatterSim (Microsoft, 2024)

MatterSim is a deep learning model trained on millions of DFT calculations covering a broad chemical space. It uses graph neural network architectures (M3GNet and Graphormer) trained with active learning to achieve near-DFT accuracy for energies, forces, and stresses across elements, temperatures (0–5000 K), and pressures (0–1000 GPa).

Key capabilities reported:

- Prediction of temperature- and pressure-dependent free energies across a wide range of inorganic solids
- Construction of phase diagrams at computational cost orders of magnitude lower than DFT
- Ten-fold improvement in prediction accuracy compared to previous state-of-the-art models on off-equilibrium structures at finite temperature and pressure
- Zero-shot generalization — reasonable predictions for systems not represented in the training data

The phase diagram capability is particularly significant. Computing a full phase diagram with DFT requires thousands of calculations and weeks of compute time. MatterSim reduces this to hours or minutes, making it practical to screen many candidate compositions before committing to expensive DFT calculations.

### UMA (Meta FAIR, 2025)

The Universal Models for Atoms (UMA) family from Meta represents a further scaling of this approach. UMA models are trained on approximately 500 million unique 3D atomic structures — the largest training run reported to date — drawn from molecular, materials, and catalysis datasets simultaneously.

The key architectural innovation is **mixture of linear experts**: despite having 1.4 billion total parameters, only about 50 million parameters are active for any given atomic structure. This allows the model to scale capacity without proportionally increasing inference cost.

The reported result is that a single UMA model, without any system-specific fine-tuning, performs comparably to or better than specialized models on diverse tasks across multiple chemical domains. The model weights, code, and associated training data have been released publicly.

### What These Models Enable

Universal ML potentials shift the cost structure of atomistic simulation:

| Before | After universal MLPs |
|:-------|:---------------------|
| Each new material needs force field parameterization (weeks–months) | Apply existing model immediately |
| Accuracy limited by empirical force field functional form | Near-DFT accuracy for supported chemistries |
| Materials discovery requires full DFT workflow | Rapid screening with MLP, DFT only for top candidates |
| Specialist knowledge needed to set up each simulation | Accessible to broader community |

The direction of travel is clear: within the next few years, universal ML potentials will likely become the first tool researchers reach for when starting any atomistic simulation project — exactly as BLAST became the first tool for sequence alignment in biology after the genomics revolution.

---

## 6. A Unifying View: The Materials Discovery Pipeline

Looking across the four case studies, a common pipeline emerges:

```
IDEA / HYPOTHESIS
        │
        ▼
COMPUTATIONAL SCREENING (DFT or ML potential)
  → Formation energies, stability, transport properties
  → High-throughput: 10³–10⁶ candidates
        │
        ▼
CANDIDATE SELECTION (top 10–200 materials)
        │
        ▼
DETAILED SIMULATION (MD, free energy, mechanism)
  → Property prediction at operating conditions
  → Mechanistic understanding
        │
        ▼
EXPERIMENTAL VALIDATION (synthesis and testing)
        │
        ▼
PRODUCT / PUBLICATION
```

In the Duracell case, the pipeline ran from 130,000 candidates to 1 commercial product. In the TRI polymer case, the pipeline runs continuously in a closed loop, with each experimental result improving the generative model. In the depolymerization case, simulation provided mechanistic understanding that guided the experimental design. In all cases, computation reduced the experimental burden by orders of magnitude.

---

## 7. Limitations and Honest Caveats

These success stories are real, but they represent the high end of what is achievable. Several important caveats apply:

**Simulation accuracy is bounded by the underlying theory.** ML potentials trained on DFT-PBE data will reproduce DFT-PBE results — including DFT-PBE errors. For properties where DFT-PBE is known to be inaccurate (band gaps, correlated electron systems, van der Waals interactions), the ML potential inherits those errors.

**Training data coverage matters.** Universal ML potentials work well for chemistries and conditions similar to their training data. For novel chemistries, extreme conditions, or rare events far from equilibrium, extrapolation is unreliable. The GNoME crystal discovery database, for example, has faced scrutiny over whether a fraction of reported novel structures were near-duplicates of known materials — a reminder that scale does not substitute for careful validation.

**Synthesis and stability under real conditions are not guaranteed.** A compound that is computationally stable (negative formation energy) may be kinetically inaccessible, require extreme synthesis conditions, or degrade rapidly in operation. Computation predicts thermodynamic stability; it does not predict synthesizability.

**The loop must close experimentally.** In each success story, computation identified candidates but experiment delivered the outcome. The Duracell project took 15 years from the first calculations to a product on store shelves. The TRI polymer platform has not yet resulted in a commercial electrolyte. Computation accelerates the process; it does not replace it.

---

## 8. Summary

| Case Study | Simulation Role | Key Method | Outcome |
|:-----------|:----------------|:-----------|:--------|
| Duracell / Materials Project | Screen 130,000 compounds | High-throughput DFT | Commercial battery (2019) |
| TRI Polymer Electrolytes | Evaluate AI-generated candidates | Automated MD pipeline | Open platform, ongoing |
| Plastic Depolymerization | Mechanistic understanding | ReaxFF, classical MD | *Nature* publication (2023); new recycling approach |
| UMA / MatterSim | Universal property prediction | ML potentials at scale | Accelerates all of the above |

Molecular simulation has crossed a threshold. It is no longer a niche interpretive tool — it is an active participant in the materials design loop. The combination of high-throughput workflows (Lecture 12), ML potentials (Lecture 10), and automated experimental platforms is producing a genuinely new paradigm for how materials are discovered and developed. For engineers who understand these tools, the opportunities are substantial.

---

## References and Further Reading

- Persson, K. A. et al., Materials Project: [materialsproject.org](https://materialsproject.org)
- Xie, T. et al., *APL Machine Learning* 1, 046108 (2023) — TRI HTP-MD platform
- Yang, Z. et al., arXiv:2312.06470 (2023) — GPT-based polymer generation at TRI
- Dong, Q., Lele, A. D. et al., *Nature* 616, 488–494 (2023) — electrified spatiotemporal depolymerization of plastics
- Wood, B. M. et al., arXiv (2025) — UMA: Universal Models for Atoms (Meta FAIR)
- Yang, H. et al., arXiv:2405.04967 (2024) — MatterSim (Microsoft)
- Merchant, A. et al., *Nature* 624, 80–85 (2023) — GNoME: scaling deep learning for materials discovery
