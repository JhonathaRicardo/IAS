# Example 1 - Analysis of Laser-induced Plasma

To demonstrate the software we will display one of the data obtained in the experiments carried out in the High Power Laser Laboratory at IPEN. The interferograms were recorded using a home-made time-resolved Mach-Zehnder-like interferometer (MZI), a detailed description of which can be found in [[16]](#Reference). The laser-induced plasma was produced in the atmosphere by focusing pump pulses with energy of 200 𝜇𝐽, duration of 25 𝑓𝑠 (FWHM), and 𝑀2 ≈ 1.2, to beam waist of ≈ 4 𝜇𝑚, reaching intensities above 10<sup>16</sup> 𝑊/𝑐𝑚2[[1]](#Reference). The CCD used in the setup saved images with dimensions of 1280 × 1024. The interference pattern fringes were adjusted to be perpendicular to the laser propagation direction with a high spatial frequency to make the plasma-shifted fringes visible.

## How to use
Detailed usage of the software will be described below:

### 1. Open Files
The user must open a Target and Reference interferogram files, using [***Open File***] and [***Open Ref.***] buttons, respectively.
In this directory (*Interferogram Data*), the users can find two interferogram files to download. A Reference *Interferogram(reference).png* (in *Fig. 1. a*) and a Target *Interferogram(plasma).png* in (*Fig. 1. b*). These image files are used in this example. 
|<img src = '/Examples/Example 2/Interferogram Data/Interferogram(reference).png' width='40%'> <img src = '/Examples/Example 1/Interferogram Data/Interferogram(plasma).png' width='40%'> |
|:--:| 
| *Fig. 1. Examples of Interferogram images: Reference image (on the left), and Plasma image (on the right).* |
  >Note: The users can use *n* files as Target interferogram and *n* Reference interferogram. Or only one Reference file. 

### 2. Options - Target Type and Select Area
The user must select *Plasma* as ***[Target Type]***.
  >Note: *Gas/Vapor* type is the default.

The User must select the rectangle area of ​​the image for IAS analysis, using the [***mouse click***] or [***X/Y Coord.***]. The coordinates values depend on the monitor screen size and do not correspond with the original image size.
  > Note: The first mouse click sets the first value of coordinates [***X/Y Coord.***] and the second click sets the second [***X/Y Coord.***].

The values of the coordinates are [***X Coord.***] = (150, 280) and [***Y Coord.***] = (160,340). We suggest a BG Phase fit of 10% (*Fig. 2*). 

> **Suggestion**: We suggested to users that these ***[Input Parameters]*** can be defined after all analysis. A ***[Scaling Factor]*** equal to 1.00 $\mu m/ pixel$ facilitates the prior analysis of users.

|<img src = '/Examples/Example 1/Images/Figure2.png'> |
|:--:| 
| *Fig. 2. IAS Main Screen for example 1 - Prior result* |

### 3. Analyse Data - Prior Analysis
Use the ***[Analyse Data]*** button to apply a retrieve phase, radial phase, and density of the target. The ***[Filter Freq.]*** parameters are defined automatically by IAS code. After the first Accumulated Phase-shift $\phi$ result (*Fig. 2*), the user will be able to refine the results by modifying the ***[Options Parameters]***. 
>Note: After each modification, the users need to press again ***[Analyse Data]*** button.  

### 4. Refining the Acc. Phase-shift ($\phi$)
To refine the Acc. Phase-shift ($\phi$) the user must verify the ***[Frequency Domain]*** 2D and 1D. In this example, the $\nu_y$ position defined by code is $46$ and the $\nu_x$ is null because the fringes orientation of interferograms is horizontal. The Fig. 3.a and 3.b show that $\nu_y = 46$ is the right position to extract the phase-shift and $\Delta\nu = 5 pixel$ is a good range to apply a Gaussian filter. The Fig. 3.c is a phase-shift map retrieved from these parameters.
|<img src = '/Examples/Example 1/Images/Fig3.png'> |
|:--:| 
| *Fig. 3. a) 2D Frequency Domain with* $\nu_x = 0$ *and* $\nu_y = 46$*; b) Identification of frequency position $\nu_x and filter range $\Delta\nu$*|

From these parameters, IAS code builds the acc. phase-shift $\phi$ and uncertainty of phase-shift $\sigma_{\phi}$ map.
|<img src = '/Examples/Example 1/Images/Fig4.png'> |
|:--:| 
| *Fig. 4. a) Phase-shift and b) uncertainty phase-shift map of the gas target using 5 pixels as Gaussian blur filter.*|

### 5. Improving The Target Symmetry to Apply Abel-inversion
The accuracy of applying the Abel-inversion is associated with the symmetry of the accumulated phase map. For some cases like this example, the symmetry failure of targets can promote great error in retrieving of density profile. To improve the cylindrical symmetry of these targets, the IAS code allows to users apply a Gaussian blur filter, which able for an improvement of symmetry. So, improve an Abel-inversion result. Fig. 5 shows two examples of Gaussian blur application.
  > Note: The Gaussian blur application can mask results. however, it is an important tool for constructing symmetrical images.
|<img src = '/Examples/Example 1/Images/Fig5.png'> |
|:--:| 
| *Fig. 5. Phase-shift map of the gas target using a Gaussian blur of a) 10 pixels and b) 20 pixels.* |

### 6. Radial Phase-shift Analysis
Applying a Gaussian blur filter of 20 pixels and setting the axis of symmetry in 105th-pixel position, the user retrieves a radial phase-shift shown in Fig. 7.


### 7. Retrieve a Density Distribution of Gaseous Target
During data collection, we used a $500\mu m$ diameter dental probe as a calibration parameter. For this example, the user can set ***[Scaling Factor]*** as 1.81 $\mu m/ pixel$ and N<sub>2</sub> as Gas type. The wavelength was indicated before and is equal to ($395 \pm 5 nm$). Applying these experimental parameters the users obtain the wished density distribution of the target (Fig. 8).
The profiles of this density map can be visualized using a ***[1D Plot]*** and the slide control in ***[1D Frame]***. 
