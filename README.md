# <h1 align = "center"> Interferometry Analysis Software (IAS) (v.3.0.0)</h1>
<p align="justify">
Interferometry Analysis Software (IAS) is a dedicated tool for studying gas, vapor, and plasma jets. IAS uses a new algorithm and GUI to analyze interferograms. It can retrieve the accumulated phase shift and apply Abel inversion to estimate the density profile of targets (gas, vapor, or plasma).
</p>

<p align="center">
  <img src = '/Images/Figure1.png' width="80%" align="center">
</p>

![License](https://img.shields.io/badge/license-MIT-green)
![version](https://img.shields.io/badge/version-v.3.0.0-green)
![status](https://img.shields.io/badge/status-under%20development-yellow)

## Summary
* [Introduction](#introduction)
* [Installation](#installation)
* [How to use it](#how-to-use-it)
  * [Main Screen](#main-screen)
  * [Interferograms](#interferograms)
  * [Options](#options)
  * [Gas-Jet Profile](#gas-jet-profile)
* [How it works](#how-it-works)
* [Example](#example)
* [Authors](#authors)
* [Acknowledgment](#acknowledgment)
* [License](#license)
* [Citation](#citation)
* [Reference](#reference)

## Introduction
  The development of diagnostic tools is significant for a better understanding of laser-plasma interactions [[1]](#reference). An accurate diagnostic is crucial, as instabilities in both target and laser pulses can result in low reproducibility of processes and impair the quality of the intended interaction [[2]](#reference). Among the various non-perturbing optical methods that can be used to diagnose a gaseous target [[3-6]](#reference), interferometry is a very accurate technique capable of quantifying tiny optical path differences and, therefore, suitable for measuring density variations of gases [[7, 8]](#reference) and laser-induced plasmas [[1]](#reference). The main drawback of the technique is that returns the integrated phase along the light path, requiring deconvolution methods for retrieving the target density profile. The software IAS was developed due to the need for a new diagnostic tool to aid in the characterization of supersonic gas jets, vapor, and plasmas, quickly and reliably. It is a compilation of two other software developed by our research group: *Interferometry Analysis - Gas-Jet* and *Interferometry Analysis - LIP*. IAS as a part of two complementary works: 
  - Studies to implement a laser-plasma accelerator infrastructure laser isotopic separation at the Nuclear and Energy Research Institute (IPEN), in Brazil.
  - Studies to implement a laser isotopic separation process for nuclear medicine application in the project of Brazilian Multipurpose Reactor (RMB).    

## Installation
The *Interferometry Analysis Software (IAS)* was developed in Python 3.11. The use of this algorithm requires the installation of the following packages: [NumPy](https://numpy.org/) [[9]](#reference), [Scipy](https://scipy.org/) [[10]](#reference) and [PyAbel](https://pyabel.readthedocs.io/en/latest/index.html) [[11]](#reference) for data processing, [Pillow](https://pypi.org/project/Pillow/) [[12]](#reference) and Scikit-image [[13]](#reference) for the processing of interferogram images, [Matplotlib](https://matplotlib.org/stable/index.html) [[14]](#reference) to plot results, and [PySimpleGui](https://www.pysimplegui.org/en/latest/) to create the user's template.

Users also can create a single .exe file using the [pyinstaller](https://pyinstaller.org/en/stable/) package through the following terminal command:

<code>   pyinstaller --onefile -w IAS_V3.0.0.py                </code>

Users who do not use Python IDEs can utilize the software through the executable file available for download [here](https://drive.google.com/file/d/1K_-3wm8TOzxROLyAc22AkRXlRvD3GjjI/view?usp=sharing)

## How to use it
The “Interferometry Analysis – Gas-Jet” has a graphical interface to facilitate its use, and this section provides a simple review of the software's functions and how to employ them.

### Main Screen
The Software Main Screen (*Fig. 1*) can be divided into 3 main parts: Interferograms, Options, and LIP Profile. Each of these parts will be detailed below.

|<img src = '/Images/Figure1.PNG'> |
|:--:| 
| *Fig.1. Software Main Screen* |

### Interferograms
- ***1. [Interferogram (Gas-Jet)]*** interferogram frame.

  - ***[Open File(s)]*** Open interferogram(s) file(s) with the presence of a gaseous target. Image file extensions should preferably be *.png* or *.snp.* (Newport proprietary format) for Newport CCD. However, all image extensions (*.gif*, *.jpg*, *.bmp*, etc) can be used. The path to the opened file is shown in the text box immediately above. If more than one file has been opened, each one is analyzed individually, and the average of all results is presented to the user.
    > **Warning**   
    >  Interferometry Analysis - Gas-Jet software only works with grayscale image files. 
  
  - ***[Rotate]*** The image rotation in degrees. Positive degrees promote counterclockwise rotation.  

  - ***[Original Size]*** Original dimensions of the image file (width, height). 
    > **Note** The interferogram shown is scaled to screen size (428,342) for users' viewing only. However, all processes to determine the plasma density profile are done with the original dimensions of the image file.

- ***2. [Interferogram (Ref.)]*** Scaled reference interferogram.

  - ***[Open Ref.]*** Open an undisturbed interferogram file. Image file extensions should preferably be .png or .snp. However, all image extensions (*.gif*, *.jpg*, *.bmp*, etc) can be used. The path to open the file is shown in the textbox above. Unlike interferogram gas jet files, the algorithm allows the insertion of only one reference file.
    > **Warning**   
    >  Interferometry Analysis - Gas-Jet software only works with grayscale image files. 

- ***3. [Analyse Data]*** From this command button, the software will apply data processing to generate the accumulated phase-shift map, the radial phase-shift map, and the map of the molecular density distribution of the gaseous target.

- ***4. [Clear]*** Button to clear input and output data.


### Options
- ***5. [Select Area]*** Parameters frame for users select the interferogram area to apply the algorithm. The selected area is defined by a rectangle with edges defined by X and Y coordinates. The user can select an area using the mouse click over the image or the combo box ***[Y Coord]*** and ***[X Coord]***.
  > **Note:** The first click of the mouse defines de first value of the X and Y triangle coordinates, and the second click defines the end coordinates of the triangle. Case, the initial X (or Y) is bigger than the final X (or Y), these values will be exchanged. 
  - ***[BG Phase Fit]*** This parameter defines the values used to construct the background of the accumulated phase &Delta;&phi;. Values are set based on a percentage of smaller phase values in the selected area, and this background is defined by a 2D plane (*default is 1%*).
  - 
- ***6. [Input Parameteres]*** Frame to set the experimental parameters used to obtain the interferogram. These parameters are:
  - ***[Scaling Factor]*** and ***[Uncertainty Scaling Factor]*** Interferogram scale in micrometers/pixel (*default is 1.000&plusmn;0.001 &nu;m/pixel*).
  - ***[Fringes Orientation]*** Definition of the interferogram fringes orientation (*default is horizontal*).
  
  - ***[Laser Wavelength]*** (&lambda;) and ***[Laser FHWM]*** (&Delta;&lambda;) in nm (*default is 395&plusmn;0 nm, respectively*).
  - ***[Gas type]*** list box of some types of gases: *H<sub>2</sub>*, *N<sub>2</sub>*, *He* and *Ar* (*default is N<sub>2</sub>*).
  - ***[Polarizability]*** (&alpha;) in angstrom³. This parameter usually refers to the tendency of matter to acquire an electric dipole moment when subjected to an electric field (*default is 1.710 A³*).
    > **Note**
    > The polarizability value is automatically filled in after selecting the gas type. If the user wants to use gases not yet listed, the gas polarizability value can be entered manually.

- ***7. [Analysis Parameters]*** Parameters frame to analyze the interferogram.
  - ***[Filter Frequency]*** ((&nu;<sub>0</sub>) This parameter is set automatically by the algorithm and this position defines which frequency will be used to apply the Inverse Fourier Transform and build the phase map of the plasma. This parameter is given in pixels.
  - ***[Filter Range]*** (&Delta;&nu;) frequency spread of the Gaussian frequency filter in pixel. The initial &Delta;&nu; depends on the image dimension but can changed by the user.
    > **Note:** The algorithm sets the frequency that generates a negative phase map. Because the refractive index of the plasma is less than 1. This is an intrinsic characteristic of plasmas, and it is considered in the calculations of the algorithm.
  - ***[Gaussian Blur]*** (&sigma;<sub>blur</sub>) Spread of the bi-dimensional Gaussian image filter. The standard deviation of the Gaussian filter ($\sigma$) defined by the user is equal for all axes. The ***[Gaussian Blur]*** is used to improve the target symmetry (*default is 5 pixels*).  
  
  - ***[Axisymmetric Orientation]*** Definition of the axis of symmetry (or axisymmetric) to apply the Inverse Abel Transform. The axisymmetric can be horizontal or vertical (*default is vertical*).
  - ***[Axisymmetric Position]*** Axisymmetric position on the accumulated phase map to apply the Abel inversion.

### Gas-Jet Profile
- ***8. [Stages]:*** Stages frame allows the visualization of each result of the algorithm.
  - ***[Fourier Transform]*** This frequency map (Fig. 3.a) is built from the Fourier Transform of the plasma interferogram. Among the frequency positions highlighted in the 1D frequency graph (Fig. 3.b), the software automatically selects the frequency that generates a negative phase shift map. The selected frequency is marked with a red line over a pixel line (or column) identifying the ***[Filter Frequency]*** (&nu;<sub>0</sub>). If the ***[Filter Frequency]***  is equal to zero, the software will set the new valor value automatically.   
    > **Note:** The user can change this ***[Filter Frequency]*** manually.
  
  - ***[Gaussian Filter]*** This image is the Gaussian filter map (Fig. 3.c) applied to generate the phase map. This filter is built from the selected frequency &nu;<sub>0</sub> and the range frequency &Delta;&nu;(Fig. 3.d).  

|<img src = '/Images/Figure3_A.PNG' width="80%">|
|:--:| 
| *Fig. 3. (a) 2D and (b) 1D frequency domain obtained by the interferogram Fourier Transform with the selected frequency to be filtered; (c) Gaussian filter to be applied on the selected frequency; (d) 1D frequency domain with range frequency to apply the Gaussian filter.* |

For the next three steps, users have the option of viewing the 2D maps or 1D curves with standard deviation using the ***[Standard Deviation]*** checkbox.
 
  - ***[Acc. Phase-shift]*** Accumulated phase-shift (&Delta;&phi;) of the plasma (in rad) retrieved from the interferograms.

|<img src = '/Images/Figure4.PNG' width="80%">|
|:--:| 
|*Fig. 4. (a) 2D accumulated phase-shift map and (b) 2D standard deviation map; (c) 1D accumulated phase curves and (d) standard deviation of one curve. All phase values are given in rad.*|   
    
  - ***[Radial Phase-shift]*** Radial phase-shift (&Delta;&phi;<sup>r</sup>) map in $rad/\mu m$ obtained after applying an Inverse Abel Transform from Accumulated Phase-shift map (&Delta;&phi;).

|<img src='/Images/Figure5.PNG' width="80%">|
|:--:| 
|*Fig. 5. (a) 2D radial phase-shift map and (b) 2D standard deviation map; (c) and (d) accuracy between 1D radial phase-shift and normalized phase-shift curves. All radial phase values are given in rad/&mu;m.*|  

  - ***[Density Profile]*** molecular density distribution ($N$) of the Gas-Jet in $cm^{−3}$ built from the radial phase-shift ($\Delta\phi_r$) and ***[Laser Wavelength]*** ($\lambda$).
    
|<img src='/Images/Figure6.PNG' width="80%">|
|:--:| 
|*Fig. 5. (a) 2D plasma density map and (b) 2D standard deviation map; (c) 1D plasma density curves and (d) standard deviation of one density curve. All density values are given in cm&oline;³ .*|

- ***9. [1D Profile]*** This button enables a 1D frame (*Item 15 in Fig. 1*) with options for the user to visualize the curves of each selected stage for different positions on the chosen symmetry axis.
- ***10. [2D Profile]*** This button enables the visualization of each ***[Stage]*** in 2D images.
- ***11. [Uncertainty of Measurement]*** This checkbox enables the visualization of the uncertainty of the accumulated phase, radial phase and density for 1D and 2D profiles.
- ***12. [Result Frame]*** In this frame user can verify the results of each stage of the data processing. The results can be seen in 1D or 2D. Results in 2D have options for ***[Colormap distribution]*** With this list box the user can choose between three colormaps distributions: linear (*Fig 7.a*), quadratic (*Fig 7.a*), or cubic(*Fig 7.c*).
- ***13. [Save Plot]*** This button allows the user to save the visualized plot as an image file (*.png*, *.jpg*, *.bmp*, etc).
- ***14. [Save Data]*** This button allows the user to save the 2D array that generated the visualized plot as a *.dat* or *.txt* file.

|<img src = '/Images/Figure7.PNG' width="80%">|
|:--:| 
|*Fig. 7. Examples with the colormaps distributions: (a) Linear distribution, (b) Quadratic distribution, (c) Cubic Distribution.* |

- ***15. [1D Profile]*** This Frame allows the visualization of 1D profiles (accumulated phase, radial phase and density) for different positions over the Axisymmetry.

## How it works
A detailed description of the algorithm will be presented in a future article. However, the summarized data processing by the software algorithm is described by the flowchart shown in *Fig. 8*:

|<img src = '/Images/Figure8.PNG' width="80%">|
|:--:| 
| *Fig. 8. Scheme of the algorithm data processing.* |

In the scheme of the algorithm data processing (*Fig. 8*): *I<sub>Target</sub>* and *I<sub>Reference</sub>* are the intensity functions of the bi-dimensional fringes fields obtained from the gas and the reference interferograms, respectively. The hats denote the Fourier transform of the intensities, *N<sup>gas</sup>* is the gas jet density calculated from its refractive index, *n*, and polarizability, &alpha;, using the Lorentz-Lorenz relation [[14,15]](#reference).

## Example
In the Example folder of this repository, the user will find the two interferogram images shown in *Fig. 9*. These images were obtained using a Mach-Zehnder-like interferometer, as discussed in [[16]](#reference).

|<img src = '/Example/Interferogram_reference.png' width='40%'> <img src = '/Example/Interferogram_gas.png' width='40%'> |
|:--:| 
| *Fig. 9. Examples of Interferogram images: Reference image (on the left), and gas-jet image (on the right).* |

The gas jet in this example was generated by a homemade supersonic nozzle [[17]](#reference) fixed inside a vacuum chamber. During data collection, the nozzle backpressure was maintained at $~50$ $bar$ with nitrogen gas. The interferograms were generated by the second harmonic (395 nm) of ultrashort pulses from a Ti:sapphire multipass CPA system [[1]](#reference).
All input parameters used for gas-jet characterization are shown in *Fig. 10*.

|<img src = '/Example/MainScreen_Example.png'> |
|:--:| 
| *Fig. 10. Software main screen: input parameters and result of the example.* |

## Authors
Interferometry Analysis - Gas-Jet software was developed by researchers of the High-Power Ultrashort Pulse Lasers Group from the Center for Lasers and
Applications (CLA) from the Instituto de Pesquisas Energéticas e Nucleares ([IPEN](https://www.ipen.br/portal_por/portal/default.php)), and of the Instituto Tecnológico de Aeronáutica ([ITA](http://www.ita.br/)).

* Jhonatha Ricardo dos Santos [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0001-7877-0580)
* Armando Valter Felicio Zuffi [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0001-5705-1499)
* Nilson Dias Vieira Junior [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0003-0092-9357)
* Edison Puig Maldonado [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0002-9462-8151)
* Ricardo Elgul Samad [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0001-7762-8961)

## Acknowledgment
Interferogram Analysis – Gas-Jet was developed to help with the analysis of the gas-jets generated by micro nozzles. This software is part of the work performed by the High-Power Ultrashort Pulse Lasers Group of the CLA/IPEN in the [Extreme Light Laboratory](https://www.unl.edu/diocles/home) (ELL)  at the University of Nebraska - Lincoln (UNL). This partnership was able due to funding provided by the São Paulo Research Foundation (FAPESP) and by the US Department of Energy, through the [LaserNetUS](https://lasernetus.org/) program.

The authors thank all members of the ELL team, especially Prof. Dr. Donald Umstadter director of ELL.

The author Jhonatha Ricardo dos Santos also acknowledges the FAPESP for doctoral fellowship 2021/13737-8. 

## License
Interferogram Analysis - Gas-Jet Profile is licensed under the [MIT license](/LICENSE).

Copyright (c) 2023 Jhonatha Ricardo dos Santos

## Citation
You can find the DOI for the latest version at [Zenodo](https://zenodo.org/badge/latestdoi/601841375).

## Reference
- [1] A. V. F. Zuffi, J. R. dos Santos, E. P. Maldonado, N D. Vieira, and R. E. Samad, "Femtosecond laser-plasma dynamics study by a time-resolved Mach–Zehnder-like interferometer," Appl. Opt. 62, C128-C134 (2023) [DOI: 10.1364/AO.477395](https://doi.org/10.1364/AO.477395).
- [2]  P. Sprangle, B. Hafizi, and J. R. Peñano, “Laser pulse modulation instabilities in plasma channels,” Phys. Rev. E 61, 4381–4393 (2000).[DOI: 10.1103/PhysRevE.61.4381](https://doi.org/10.1103/PhysRevE.61.4381).
- [3] G. Costa, M. P. Anania, F. Bisesto, E. Chiadroni, A. Cianchi, A. Curcio,M. Ferrario, F. Filippi, A. Marocchino, F. Mira, R. Pompili, and A. Zigler,“Characterization of self-injected electron beams from LWFA experiments at SPARC_LAB,” Nucl. Instrum. Methods A 909, 118–122 (2018).[DOI 10.1016/j.nima.2018.02.008](https://doi.org/10.1016/j.nima.2018.02.008).
- [4] G. S. Settles, Schlieren and shadowgraph techniques: visualizing phenomena in transparent media, in Experimental Fluid Mechanics (Springer, 2001), pp. xviii.
- [5] S. Shiraishi, C. Benedetti, A. J. Gonsalves, K. Nakamura, B. H. Shaw, T. Sokollik, J. van Tilborg, C. G. R. Geddes, C. B. Schroeder, C. Toth, E. Esarey, and W. P. Leemans, “Laser red shifting based characterization of wakefield excitation in a laser-plasma accelerator,” Phys. Plasmas 20, 063103 (2013).[DOI 10.1063/1.4810802](https://doi.org/10.1063/1.4810802).
- [6] A. J. Goers, G. A. Hine, L. Feder, B. Miao, F. Salehi, J. K. Wahlstrand, and H. M. Milchberg, “Multi-MeV electron acceleration by Subterawatt laser pulses,” Phys. Rev. Lett. 115, 194802 (2015).[DOI 10.1103/PhysRevLett.115.194802](https://doi.org/10.1103/PhysRevLett.115.194802).
- [7] F. Brandi and L. A. Gizzi, “Optical diagnostics for density measurement in high-quality laser-plasma electron accelerators,” High Power Laser Sci. Eng. 7, e26 (2019).[DOI 10.1017/hpl.2019.11](https://doi.org/10.1017/hpl.2019.11).
- [8] A. K. Arunachalam, “Investigation of laser-plasma interactions at near-critical densities,” Dissertation (University of Jena, 2017)..
- [9] Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). [DOI: 10.1038/s41586-020-2649-2](https://www.nature.com/articles/s41586-020-2649-2). 
- [10] Pauli Virtanen, et. al. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272. [DOI: 10.1038/s41592-019-0686-2](https://www.nature.com/articles/s41592-019-0686-2).
- [11] Gibson, Stephen; Hickstein, Daniel D.; Yurchak, Roman; Ryazanov, Mikhail; Das, Dhrubajyoti; Shih, Gilbert.(2022) PyAbel, PyAbel: v0.9.0, Zenodo,  [DOI: 10.5281/zenodo.7438595](https://doi.org/10.5281/zenodo.7438595).
- [12] Clark, A. (2015). Pillow (PIL Fork) Documentation. readthedocs. Retrieved from [https://buildmedia.readthedocs.org/media/pdf/pillow/latest/pillow.pdf](https://buildmedia.readthedocs.org/media/pdf/pillow/latest/pillow.pdf).
- [13] Stéfan van der Walt, Johannes L. Schönberger, Juan Nunez-Iglesias, François Boulogne, Joshua D. Warner, Neil Yager, Emmanuelle Gouillart,
Tony Yu and the scikit-image contributors. scikit-image: Image processing in Python. PeerJ 2:e453 (2014). [DOI: 10.7717/peerj.453](https://doi.org/10.7717/peerj.453)
- [13] J. D. Hunter, Matplotlib: A 2D Graphics Environment. Computing in Science & Engineering, 9 (3), 90-95 (2007). [
DOI: 10.1109/MCSE.2007.55] (https://ieeexplore.ieee.org/document/4160265)
- [14] H. A. Lorentz, "Über die Beziehungzwischen der Fortpflanzungsgeschwindigkeit des Lichtes derKörperdichte", Ann. Phys. 9, 41-665, (1880). [DOI: 10.1002/andp.18802450406]( https://doi.org/10.1002/andp.18802450406)
- [15] L. Lorenz, "Über die Refractionsconstante", Ann. Phys. 11, 70-103  (1880). [DOI: 10.1002/andp.18802470905](https://doi.org/10.1002/andp.18802470905)
- [16] A. V. F. Zuffi, E. P. Maldonado, N. D. Vieira, and R. E. Samad, “Development of a modified Mach-Zehnder interferometer for time and space density measurements for laser wakefield acceleration,” in SBFoton International Optics and Photonics Conference (IEEE, 2021).
- [17] B. B. Chiomento, A. V. F. Zuffi, N. D. V. Junior, F. B. D. Tabacow, E. P. Maldonado, and R. E. Samad, “Development of dielectric de Laval
nozzles for laser electron acceleration by ultrashort pulses micromachining,” in SBFoton International Optics and Photonics Conference (IEEE, 2021).



