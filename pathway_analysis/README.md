# Pathway Analysis Shiny App  
**Repository:** `pathway_analysis` (part of `Rshiny-apps` by Alex Galaras)  
**Version:** 1.0.0  
**Date:** 2025-11-06  
**Author:** Alexandros Galaras  

---

## Table of Contents  
1. [Overview](#overview)  
2. [Features](#features)  
3. [Getting Started](#getting-started)  
   - [Prerequisites](#prerequisites)  
   - [Installation](#installation)  
   - [Running the App](#running-the-app)
4. [Getting Started](#getting-started) 

## Overview  
This Shiny application enables pathway analysis of gene expression datasets, using ClusterProfiler, ultimately providing an interactive interface for uploading data, and pathway visualization for both upregulated & downregulated genes.

---

## Features  
- Upload your own dataset (expression matrix)  
- Choose from built-in pathway databases (BP, MF, CC)  
- Visualise pathway enrichment results (bar plots, dot plots) for both upregulated and downregulated genes 
- Interactive tables for result sorting, filtering and exploration  
- Downloadable results (tables, figures) for downstream use  
- Customisation of analysis parameters: p-value cutoff, minimum gene count, etc.

---

## Getting Started

### Prerequisites  
- R version = 4.5  
- Install packages located in install.R

### Installation  
1. Clone the repository:  
   ```bash
   git clone https://github.com/alex-galaras/Rshiny-apps.git
   cd Rshiny-apps/pathway_analysis

## Data Input

Upload a gene expression file (txt or csv format) containing the headers 'gene	log2FC	Stat'. The Stat column can contain either pvalue or FDR based on your preference.
