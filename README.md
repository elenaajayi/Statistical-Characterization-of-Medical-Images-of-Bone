# Statistical-Characterization-of-Bone

**Statistical Comparison of Medical and Natural Images**
This repository contains the code and data for a study that compares the statistical characteristics of medical images and natural images, focusing on the analysis of spatial frequency content and local features. The investigation explores how the formation of medical images differs from natural images and how these differences impact their statistical properties, particularly for radiologists in the diagnostic process.


**Overview**
Visual perception is evolved to process natural images, but medical images are generated through different physical processes, which may result in different statistical characteristics. This project focuses on two types of medical images: bone radiographs and scintigrams.

**The key analyses performed include:**
Spatial Frequency Content: Comparing the spatial power spectra of medical and natural images.
Local Image Features: Examining correlations among image patches to assess informativeness.

**Dataset**
Medical Images: **51 radiographs and 20 scintigrams**, obtained from the public **MedPix** database. Images were curated to remove any labeling or artifacts.
Patch Extraction: Images were divided into regions of interest (ROIs) and further subdivided into 64 x 64-pixel patches, resulting in:
**1918 patches from radiographs
1405 patches from scintigrams**

**Methodology**
Spatial Power Spectrum Analysis:
The spectral slope of the power spectra was calculated for both radiographs and scintigrams.
Medical images had a **spectral slope of -3.4,** contrasting with the t**ypical slope of -2** in natural images.
Scintigram spectra showed a decrease in slope above 0.1 cycles per pixel.

**Local Feature Analysis:**
A custom pipeline was used to extract local image statistics:
Images were whitened based on their power spectra.
Binarization at the median was performed.
Pairwise, triplet, and quadruplet correlations were computed at several scales.
The analysis was compared to a parallel study on natural images **(Hermundstad et al., eLife 2014**).

**Key Findings**
Spatial Frequency Differences: Medical images (both radiographs and scintigrams) have distinct spatial frequency characteristics when compared to natural images.
Local Correlations:
Pairwise correlations were less informative in medical images than in natural images.
Triplet correlations were especially uninformative in bone radiographs, which may be due to a relative lack of T-junctions in these images.

**Conclusion**
The statistical properties of medical images differ significantly from those of natural images in both spatial frequency and local image statistics. This study sheds light on how medical images are processed and highlights the need for tailored approaches when applying vision models to medical imaging.


