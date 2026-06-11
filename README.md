# SJTU-IEEE-Second-Semester-Freshman-Course

This repository contains all the notes, assignments, and materials from the second semester of my freshman year in the IEEE program at Shanghai Jiao Tong University (SJTU).

**Author:** Li Jinshuo (李谨硕)

**License:** [MIT](LICENSE)

---

## Table of Contents

- [Courses](#courses)
  - [CS1604 — Program Design](#cs1604--program-design)
  - [EST2501 — Digital Electronics](#est2501--digital-electronics)
  - [MATH1204 — Mathematical Analysis](#math1204--mathematical-analysis)
  - [MATH1409 — Linear Algebra for Artificial Intelligence](#math1409--linear-algebra-for-artificial-intelligence)
- [Repository Structure](#repository-structure)
- [Usage](#usage)

---

## Courses

### CS1604 — Program Design

An introductory C++ programming course covering fundamental programming concepts, control flow, data structures, and algorithms.

| Directory | Description |
|---|---|
| `note/` | Lecture notes (lec1–lec7) in Markdown, covering C++ syntax, types, operators, control flow, functions, arrays, pointers, classes, and more |
| `assignment0/` | Warm-up: Hello World and A+B problem, with a batch judger |
| `assignment1/` | Four problems: string length, binary conversion, prime numbers, and encryption |
| `assignment2/` | Three problems: string comparison (comp), pattern matching (match), and a Tic-Tac-Toe game |
| `assignment3/` | Three problems: Wordle game, fraction arithmetic (a+b), and regex matching |
| `assignment4/` | Three problems: browser history, maze pathfinding, and an assembly interpreter |
| `assignment5/` | Two problems: matrix operations (custom `Matrix` class) and a versioned array (`gitArray`) |
| `midterm/` | Midterm practice: 10 LeetCode-style C++ solutions and two midterm test files |

Each assignment includes C++ source files, test data (`.in`/`.out`), and Python-based auto-judger scripts.

### EST2501 — Digital Electronics

A course on digital circuits and logic design fundamentals.

| Directory / File | Description |
|---|---|
| `final.md` | Comprehensive final exam review notes covering digital concepts, logic gates, Boolean algebra, combinational/sequential circuits, flip-flops, counters, and 555 timers |
| `final.pdf` | Compiled PDF version of the review notes |
| `final_tex/` | LaTeX source for a polished digital circuits reference document, compiled with XeLaTeX. Includes circuit diagrams (SR latches, JK flip-flops, counters, 555 timers, etc.) sourced from *Digital Fundamentals* by Thomas Floyd |

### MATH1204 — Mathematical Analysis

A continuation of calculus/analysis, focusing on multivariable integration.

| Directory / File | Description |
|---|---|
| `final.md` | Final exam review: multiple integrals (double/triple integrals), area/volume computation, change of variables, line integrals, and surface integrals |
| `midterm.md` | Midterm review notes |
| `mai-ii.pdf` / `maii.md` | Supplementary reference materials |

### MATH1409 — Linear Algebra for Artificial Intelligence

An advanced linear algebra course tailored for AI applications, covering linear spaces, matrix decompositions, and computational methods.

| Directory / File | Description |
|---|---|
| `la_final.md` | Final exam review: linear spaces, linear transformations, Jordan canonical form, matrix decompositions (full-rank, QR, LU, spectral, SVD), vector norms, and gradient matrices |
| `la_update/` | LaTeX source and compiled PDF for updated linear algebra notes |
| `report/` | QR decomposition report (LaTeX + PDF) |
| `program_solution/` | Programming assignments and projects |

**Programming Solutions** — A custom pure-Python matrix library (`MatLibrary`) plus C++ implementations:

| Path | Description |
|---|---|
| `main.py` / `matlibrary.py` | Unified `Matrix` class with dual arithmetic modes (exact `fraction` and numerical `float`), supporting Gaussian elimination, determinants, inverses, LU/Cholesky/QR decompositions, eigenvalue computation (Jacobi + QR iteration), SVD, and Jordan normal form |
| `cpp/` | C++ matrix library implementation (`matlibrary.cpp`) with a driver program |
| `python_solution/` | Individual homework solutions: Cholesky decomposition, incomplete LU, SVD, and a full matrix class project |
| `more/fastmatrix_linux/` | A C++ accelerated matrix backend (`libfastmatrix.so`) with Python bindings and benchmark scripts |
| `README.md` / `README_en.md` / `README_cn.md` | MatLibrary documentation in English and Chinese |

---

## Repository Structure

```
.
├── LICENSE
├── README.md
├── CS1604_Program_Design/
│   ├── note/                  # Lecture notes (Markdown)
│   ├── assignment0/           # Warm-up assignments
│   ├── assignment1/           # Assignment 1
│   ├── assignment2/           # Assignment 2
│   ├── assignment3/           # Assignment 3
│   ├── assignment4/           # Assignment 4
│   ├── assignment5/           # Assignment 5
│   └── midterm/               # Midterm practice problems
├── EST2501_Digital_Electronics/
│   ├── final.md / final.pdf   # Final exam review
│   └── final_tex/             # LaTeX source with diagrams
├── MATH1204_Mathematical_Analysis/
│   ├── final.md               # Final exam review
│   ├── midterm.md             # Midterm review
│   └── mai-ii.pdf / maii.md   # Supplementary materials
└── MATH1409_Linear_Algebra_for_Artificial_Intelligence/
    ├── la_final.md            # Final exam review
    ├── la_update/             # LaTeX notes
    ├── report/                # QR decomposition report
    └── program_solution/      # MatLibrary and homework solutions
```

---

## Usage

- **Lecture notes** are standalone Markdown files — open with any Markdown viewer or editor.
- **C++ assignments** can be compiled with any C++ compiler (e.g., `g++`). Use the provided `judger.py` / `judger_batch.py` scripts to auto-test solutions against the supplied test data.
- **LaTeX documents** (Digital Circuits, Linear Algebra notes, QR report) should be compiled with XeLaTeX or pdfLaTeX as appropriate.
- **MatLibrary** requires Python 3.10+. See `MATH1409_Linear_Algebra_for_Artificial_Intelligence/program_solution/README.md` for full documentation.
