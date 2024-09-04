# Example 2 - Analysis of Laser-induced Plasma

To demonstrate the software we will display one of the data obtained in the experiments carried out in the High Power Laser Laboratory at IPEN. The interferograms were recorded using a home-made time-resolved Mach-Zehnder-like interferometer (MZI), a detailed description of which can be found in [[16]](#Reference). The laser-induced plasma was produced in the atmosphere by focusing pump pulses with energy of 200 𝜇𝐽, duration of 25 𝑓𝑠 (FWHM), and 𝑀2 ≈ 1.2, to beam waist of ≈ 4 𝜇𝑚, reaching intensities above 10<sup>16</sup> 𝑊/𝑐𝑚2[[1]](#Reference). The CCD used in the setup saved images with dimensions of 1280 × 1024. The interference pattern fringes were adjusted to be perpendicular to the laser propagation direction with a high spatial frequency to make the plasma-shifted fringes visible.

## How to use
Detailed usage of the software will be described below:

### 1. Open Files
The user must open a Target and Reference interferogram files, using [***Open File***] and [***Open Ref.***] buttons, respectively.
In this directory (*Interferogram Data*), the users can find two interferogram files to download. A Reference *Interferogram(reference).png* (in *Fig. 1. a*) and a Target *Interferogram(plasma).png* in (*Fig. 1. b*). These image files are used in this example. 
|<img src = '/Examples/Example 2/Interferogram Data/interferogram (reference).png' width='40%'> <img src = '/Examples/Example 2/Interferogram Data/interferogram (plasma).png' width='40%'> |
|:--:| 
| *Fig. 1. Examples of Interferogram images: Reference image (on the left), and Plasma image (on the right).* |
  >Note: The users can use *n* files as Target interferogram and *n* Reference interferogram. Or only one Reference file. 

### 2. Options - Target Type and Select Area
The user must select *Plasma* as ***[Target Type]***.
  >Note: *Gas/Vapor* type is the default.

The User must select the rectangle area of ​​the image for IAS analysis, using the [***mouse click***] or [***X/Y Coord.***]. The coordinates values depend on the monitor screen size and do not correspond with the original image size.
  > Note: The first mouse click sets the first value of coordinates [***X/Y Coord.***] and the second click sets the second [***X/Y Coord.***].

The values of the coordinates are [***X Coord.***] = (115, 235) and [***Y Coord.***] = (55,95). We suggest a BG Phase fit of 10% (*Fig. 2*). 

> **Suggestion**: We suggested to users that these ***[Input Parameters]*** can be defined after all analysis. A ***[Scaling Factor]*** equal to 1.00 $\mu m/ pixel$ facilitates the prior analysis of users.

|<img src = '/Examples/Example 2/Images/Figure2Ex2.png'> |
|:--:| 
| *Fig. 2. IAS Main Screen for example 2 - Prior result* |

### 3. Analyse Data - Prior Analysis
Use the ***[Analyse Data]*** button to apply a retrieve phase, radial phase, and density of the target. The ***[Filter Freq.]*** parameters are defined automatically by IAS code. After the first Accumulated Phase-shift $\phi$ result (*Fig. 2*), the user will be able to refine the results by modifying the ***[Options Parameters]***. 
>Note: After each modification, the users need to press again ***[Analyse Data]*** button.  

### 4. Refining the Acc. Phase-shift ($\phi$)
To refine the Acc. Phase-shift ($\phi$) the user must verify the ***[Frequency Domain]*** 2D and 1D. In this example, the $\nu_x$ position defined by code is $267$ and the $\nu_y$ is null because the fringes orientation of interferograms is vertical. The Fig. 3.a and 3.b show that $\nu_x = 267$ is the right position to extract the phase-shift and $\Delta\nu = 7 pixel$ is a good range to apply a Gaussian filter.
|<img src = '/Examples/Example 2/Images/Figure3Ex2.png'> |
|:--:| 
| *Fig. 3. a) 2D Frequency Domain with* $\nu_x = 267$ *and* $\nu_y = 0$*; b) Identification of frequency position $\nu_x and filter range $\Delta\nu$*|
>Note: For plasmas analysis the IAS code defines the frequency $-\nu_x$ and $-\nu_y$ as default. But, these parameters can depend on the interferogram data. 

From these parameters, IAS code builds the acc. phase-shift $\phi$ and uncertainty of phase-shift $\sigma_{\phi}$ map using a Gaussian blur of 5 pixels (*Default value*).
|<img src = '/Examples/Example 2/Images/Figure4Ex2.png'> |
|:--:| 
| *Fig. 4. a) Phase-shift and b) uncertainty phase-shift map of the gas target using 5 pixels as Gaussian blur filter.*|

### 5. Improving The Target Symmetry to Apply Abel-inversion
The accuracy of applying the Abel-inversion is associated with the symmetry of the accumulated phase map. In general, the phase-shift of laser-induced plasmas presents a good cylindrical symmetry. This way, the default value of the Gaussian blur (5 pixels) can be reduced. Fig. 5.b shows a radial phase-shift with failure in symmetry. It was retrieved from Abel-inversion in the phase-shift map without Gaussian blur (Fig. 5.a). To improve the cylindrical symmetry of this target, we apply a Gaussian blur of 2 pixels in phase-shift (Fig. 5.c) to obtain a better result of radial phase-shift (Fig. 5.d).
  > Note: The Gaussian blur application can mask results. however, it is an important tool for constructing symmetrical images.
|<img src = '/Examples/Example 2/Images/Figure5Ex2.png'> |
|:--:| 
| *Fig. 5. a) Phase-shift, and b) radial phase-shift map of plasma without Gaussian blur. c) Phase-shift, and d) radial phase-shift map of plasma using a Gaussian blur of 2 pixels.* |

### 6. Radial Phase-shift Analysis
Applying a Gaussian blur filter of 2 pixels and setting the axis of symmetry in 46th-pixel position, the user retrieves a radial phase-shift and uncertainty map, shown in Fig. 7.a and 7.b, respectively.
|<img src = '/Examples/Example 2/Images/Figure6Ex2.png'> |
|:--:| 
| *Fig. 6. a) Radial phase-shift, and b) uncertainty radial phase-shift maps of the plasma using a Gaussian blur of 20 pixels.* |

### 7. Retrieve a Density Distribution of Gaseous Target
During data collection, we used a $500\mu m$ diameter dental probe as a calibration parameter. For this example, the user can set ***[Scaling Factor]*** as 1.81 $\mu m/ pixel$ and N<sub>2</sub> as Gas type. The wavelength was indicated before and is equal to ($395 \pm 5 nm$). Applying these experimental parameters the users obtain the wished density distribution of the target (Fig. 8).
The profiles of this density map can be visualized using a ***[1D Plot]*** and the slide control in ***[1D Frame]***. 
|<img src = '/Examples/Example 2/Images/Figure7Ex2.png'> |
|:--:| 
| *Fig. 7. a) Electron density, and b) uncertainty density distribution of plasma* |
