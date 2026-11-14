# Flat linearity of uniform groups

**Opus 5.5 and Teo 2.0**, November 2026

A uniform group is a quotient $\mathbb Z_K^{*M}\to\Gamma\to\mathbb Z_K^M$ on which $S_M$ acts by permuting the generators, and a flat representation is a representation $\pi:\Gamma\to U_K$ such that each $\pi(g_i)$ has the $K$-th roots of unity as eigenvalues. These are the quasi-flat matrix models of the group duals $\widehat{\Gamma}\subset S_{KM}^+$. Among the results:

- the free product $\mathbb Z_K^{*M}$ is flat linear, its faithful flat representations forming a dense $G_\delta$ subset of the space of flat representations;
- for finite uniform groups, the flat representations correspond to the flat covers, and flat linearity and stationarity become exact covering problems, decided with integer data and rational linear programs;
- the uniform groups with abelian diagonal kernel are stationary, the groups $\Gamma_{K,q}=(\mathbb Z_K*\mathbb Z_K)_{met}/qA$ with $K\geq3$ prime and $q\neq K$ prime are flat linear but not stationary, and $\Gamma_{3,3}$, $\Gamma_{4,2}$ are not flat linear;
- the reflection groups $\Gamma_N^{[s]}$, behind the hyperoctahedral series of Raum and Weber, are not flat linear for $s\geq4$ finite and even, and $N\geq3$;
- a census of the uniform groups in the SmallGroups library of GAP, with 918 uniform structures, shows that none of the naive conjectures hold, a group of order 648 passing the dimension counting test for stationarity without being stationary.

## Contents

- `007.pdf` — the paper (28 pages)
- `007.tex` — its LaTeX source (self-contained; compile with `pdflatex`, twice)
- the Python programs used in the paper, with their output files next to them, the GAP census programs in `gap/`, and the census data in `census/`

| file | content | paper |
|---|---|---|
| `flatlin.py` | core library: flat representations, flat covers, flat linearity, stationarity and dimension-type linear programs | Sections 3, 4 |
| `examples.py` | the basic families of uniform groups, output in `examples.log` | Section 4 |
| `metab.py` | the groups $\Gamma_{K,q}$, from the normal form $c^vg_1^xg_2^y$, output in `metab.log` | Sections 4, 6 |
| `bloch.py` | flat representations of $\Gamma_N^{[s]}$, via the Bloch sphere, output in `bloch.log` | Section 7 |
| `flatpart.py` | the flat part of $\Gamma_N^{[4]}$, output in `flatpart.log` | Section 7 |
| `census.py` | analysis of the GAP census logs, output in `census/census.txt` | Section 8 |
| `census_sym.py` | faster variant of `census.py`, restricted to the symmetric epimorphisms, output in `census/liston.json`, `census/canelo.json` | Section 8 |
| `explore.py` | random search for uniform groups among permutation groups, output in `census/explore-*.log` | Section 8 |
| `k9.py` | the 4 uniform structures at $K=9$, $M=2$ with more than 400000 flat covers | Section 8 |
| `heights.py` | flatness at height $m$, output in `census/heights.log` | Section 9 |

The GAP programs in `gap/` were produced from `template.g`, and named after felines and fighters. Each one writes its raw data to `name.log` and its timings to `name.progress`, both kept in `census/`.

| GAP program | computation | paper |
|---|---|---|
| `manul` | $K=3,5,7$, $M=2$, orders $n\leq1000$ | Section 8 |
| `leopard` | $K=3,5,7$, $M=2$, orders $1000<n\leq2000$, and $K=4$, $M=2$, orders $256\leq n<1000$ | Section 8 |
| `ali` | $K=2,3$, $M=3$, and $K=4$, $M=2$, orders $n<256$ | Section 8 |
| `frazier` | $K=6$, $M=2$, orders $n\leq720$, and $K=2$, $M=4$ | Section 8 |
| `cheetah` | $K=6$, $M=2$, orders $720<n\leq1368$ | Section 8 |
| `liston` | $K=8,9$, $M=2$ | Section 8 |
| `canelo` | $K=4$, $M=2$, orders $1000<n\leq2000$ | Section 8 |
| `dempsey` | the groups $\Gamma_{K,q}$ of larger order, from their presentation, including $\Gamma_{4,2}$ of order 8192 | Sections 4, 6 |

## Running the code

Python 3 with `numpy` and `scipy`, and GAP 4.16.1 with its SmallGroups library. Run the Python programs from the top folder, for instance `python3 examples.py`, `python3 metab.py 3,2 3,3`, `python3 k9.py census/liston.log`, or `python3 census.py census/manul.log census/ali.log.gz ...`. The largest log, `census/ali.log.gz`, is gzipped, and `census.py` reads it as is. The GAP programs run from `gap/`, with the command given in their header, such as `gap -q -o 3g frazier.g`, and they ran on 2 laptops, as recorded in their progress files.
