# Statistical Characterization of Medical Images of Bone

Research conducted at the **Feil Family Brain & Mind Research Institute, Weill Cornell Medicine** (Jun 2022 - May 2023), under the supervision of [**Dr. Jonathan Victor**](https://github.com/jvlab) ([Lab Page](http://www-users.med.cornell.edu/~jdvicto/labshots.html)).

## Overview

Visual perception is evolved to process natural images, but medical images are generated through different physical processes, which may result in different statistical characteristics. This project focuses on two types of medical images: bone radiographs and scintigrams.

## Attribution

Most of the analysis code in this repository was written by **Dr. Jonathan Victor** ([github.com/jvlab](https://github.com/jvlab)). My contributions included running and tweaking the analysis scripts, processing and curating the medical images, writing portions of the VSS 2023 poster, and presenting the work at LANS 2022.

## Key Analyses

- **Spatial Frequency Content**: Comparing the spatial power spectra of medical and natural images
- **Local Image Features**: Examining correlations among image patches to assess informativeness

## Dataset

- **Medical Images**: 51 radiographs and 20 scintigrams, obtained from the public MedPix database. Images were curated to remove any labeling or artifacts.
- **Patch Extraction**: Images were divided into regions of interest (ROIs) and further subdivided into 64 x 64-pixel patches, resulting in 1918 patches from radiographs and 1405 patches from scintigrams.

## Methodology

**Spatial Power Spectrum Analysis:**
- The spectral slope of the power spectra was calculated for both radiographs and scintigrams
- Medical images had a spectral slope of -3.4, contrasting with the typical slope of -2 in natural images
- Scintigram spectra showed a decrease in slope above 0.1 cycles per pixel

**Local Feature Analysis:**
- Images were whitened based on their power spectra
- Binarization at the median was performed
- Pairwise, triplet, and quadruplet correlations were computed at several scales
- The analysis was compared to a parallel study on natural images (Hermundstad et al., eLife 2014)

## Key Findings

- **Spatial Frequency Differences**: Medical images (both radiographs and scintigrams) have distinct spatial frequency characteristics when compared to natural images
- **Pairwise correlations** were less informative in medical images than in natural images
- **Triplet correlations** were especially uninformative in bone radiographs, which may be due to a relative lack of T-junctions in these images

## Publications & Presentations

- **VSS 2023** — Poster presented at the Vision Sciences Society annual meeting ([Abstract in Journal of Vision](https://jov.arvojournals.org/article.aspx?articleid=2791798))
- **LANS 2022** — Oral presentation at the Laboratory for Auditory Neuroscience Symposium

## Repository Structure

```
├── src/             Latest MATLAB analysis scripts
├── lib/             Utility functions (ffdm_package)
├── data/            Excel databases
├── results/         Output files (.mat, .txt, .fig, .emf)
├── images/          Medical images (synpic dataset)
├── poster/          VSS 2023 poster materials
├── presentations/   LANS 2022 presentation, summary slides
├── archive/         Date-versioned script history
```

## Conclusion

The statistical properties of medical images differ significantly from those of natural images in both spatial frequency and local image statistics. This study sheds light on how medical images are processed and highlights the need for tailored approaches when applying vision models to medical imaging.
