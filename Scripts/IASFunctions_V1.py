#################################################################################
# FUNCTIONS OF INTERFEROGRAM ANALYSIS V.3.X.X
################################################################################
# Authors: Jhonatha Ricardo dos Santos, Armando Zuffi, Ricardo Edgul Samad, Nilson Dias Vieira Junior
# Python 3.11
# Last update: 2024_09_06

import os
import io

import numpy as np
import shutil

from io import BytesIO
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.ndimage import rotate
from scipy.signal import peak_widths, find_peaks
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.optimize import curve_fit
from PIL import Image, ImageDraw


# GET BINARY DATA
def getBinaryData(filename):
    '''
    :path file name to binary value
    :param filename: path of file
    :return: binary value
    '''
    binary_values = []
    with open(filename, 'rb') as f:
        data = f.read(1)
        while data != b'':
            binary_values.append(ord(data))
            data = f.read(1)
        return binary_values

# DRAW FIGURE FROM FILES
def draw_figure(canvas, figure):
    '''
    Drawing rectangle figure on canvas
    :param canvas: image canvas
    :param figure: original interferogram
    :return: rectangle drawn on figure
    '''
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg

# GET VALUES FROM INPUTBOX
def get_value(key, values):
    '''
    convert string labels to float values.
    :param key: labels
    :param values: labels value
    :return: float of label values
    '''
    value = values[key]
    return float(value)

# GET IMAGE FILE FROM PATH FILE
def image_to_data(im):
    '''
    convert image to data image
    :param im: image
    :return: data of image
    '''
    with BytesIO() as output:
        im.save(output, format="PNG")
        data = output.getvalue()
    return data

# DRAW RECTANGLE ON INTERFEROMETER FIGURE
def apply_drawing(values, window, tmp_file, size):
    '''
    :param values: x and y labels of rectangle
    :param window: main window
    :return: rectangle drown on temp image file
    '''

    image_file = values["file1"]
    begin_x = get_value("-BEGIN_X-", values)
    begin_y = get_value("-BEGIN_Y-", values)
    end_x = get_value("-END_X-", values)
    end_y = get_value("-END_Y-", values)
    rotate_degree = get_value("-DEGREE-", values)

    if begin_x > end_x:
        begin_x = get_value("-END_X-", values)
        end_x = get_value("-BEGIN_X-", values)
    if begin_y > end_y:
        begin_y = get_value("-END_Y-", values)
        end_y = get_value("-BEGIN_Y-", values)

    if os.path.exists(image_file):
        shutil.copy(image_file, tmp_file)
        imagetmp = Image.open(tmp_file)
        imagetmp = imagetmp.resize(size)
        imagetmp = imagetmp.rotate(rotate_degree, resample=Image.Resampling.BICUBIC)
        draw = ImageDraw.Draw(imagetmp)
        draw.rectangle((begin_x, begin_y, end_x, end_y), width=2, outline='white')  ##DCDCDC
        imagetmp.save(tmp_file)
        bio = io.BytesIO()
        imagetmp.save(bio, format='PNG')
        window["image1"].update(data=bio.getvalue(), size=size)

def rotate_img(values, window):
    '''
    :param values: x and y labels of rectangle
    :param window: main window
    :return: rectangle drown on temp image file
    '''
    image_file = values["file1"]
    rotate_degree = get_value("-DEGREE-", values)

    if os.path.exists(image_file):
        shutil.copy(image_file, tmp_file)
        imagetmp = Image.open(tmp_file)
        imagetmp = imagetmp.resize(size)
        imagetmp = imagetmp.rotate(rotate_degree, resample=Image.Resampling.BICUBIC)
        imagetmp.save(tmp_file)
        bio = io.BytesIO()
        imagetmp.save(bio, format='PNG')
        window["image1"].update(data=bio.getvalue(), size=size)

# CREATING MEAN MAPS/ARRAY AND STD ARRAY
def mean_maps(data):
    '''
    2D Array mean
    :param n: group of 2D arrays
    :return: 2D array
    '''
    mean_data = data[0] / len(data)
    for i in range(1, len(data)):
        mean_data = mean_data + data[i] / len(data)
    return mean_data

def std_maps(data, mean_data):
    '''
    standard deviation of 2D Array maps
    :param n: group of 2D arrays and mean 2D array
    :return: std of 2D array
    '''
    desv = (data[0] - mean_data) * (data[0] - mean_data)
    for i in range(1, len(data)):
        desv = desv + (data[i] - mean_data) * (data[i] - mean_data)

    return np.sqrt(desv / len(data))

# CREATING FRINGES WIDTHS
def fringes_width(data, fang_deg):
    '''
    Calculate 2D array shifts and widths fringes distribution
    :param n: 2D array, 2D array of ref. image.
    :return: mean fringe width
    '''
    data1 = rotate(data, 90 - fang_deg, reshape = False)
    nl, nr = np.shape(data1)
    f_width = np.zeros(np.shape(data1))
    data2 = np.zeros(np.shape(data1))
    for r in range(0, nr):
        col = data1[:, r]
        data2[:, r] = np.where(col == 0, np.max(col), col)

    for l in range(0, nl):
        li = (data2[l, :])
        try:
            ypeaks1, _ = find_peaks(li, height=0.5 * np.max(li), width=1)

            x = np.arange(0,nr,1)
            y = np.interp(x, ypeaks1[0:-1], np.diff(ypeaks1))

            f_width[l, :] = y
        except:
            f_width[l, :] = np.zeros(np.shape(li))

    result = np.asfarray(rotate(f_width, -(90 - fang_deg), reshape = False))
    result = np.where(result == 0, np.mean(f_width), result)

    return result, np.std(f_width)

def baseline2D_gas(data, base_ref):
    '''
    path file name to binary value
    :param 2D array and base ref (%)
    :return: 2D array BG map and 2D array std of BG map
    '''
    data_sort = np.sort(data.flatten())
    bg_map = np.mean(data_sort[0:int(len(data_sort)* base_ref)]) * np.ones(np.shape(data))
    bg_std = np.std(data_sort[0:int(len(data_sort)* base_ref)]) * np.ones(np.shape(data))

    return bg_map, bg_std

def baseline2D_plasma(data, base_ref):
    nlines, nrows = np.shape(data)
    X = np.arange(0, nrows, 1)
    Y = np.arange(0, nlines, 1)
    XY = np.meshgrid(X, Y)

    # creating vector from data
    dataxyz = []
    for i in range(0, nrows):
        for j in range(0, nlines):
            if (i <= int(nrows * base_ref) or i >= int(nrows - nrows * base_ref)) or (
                j <= int(nlines * base_ref) or j >= int(nlines - nlines * base_ref)):
                dataxyz.append([i, j, data[j, i]])

    # extracting x, y and z info
    xy, z = [], []
    for i in range(0, np.shape(dataxyz)[0]):
        xy.append(dataxyz[i][0:2])
        z.append(dataxyz[i][2])

    z = np.transpose(z)
    xy = np.transpose(xy)

    popt, pcov = curve_fit(func_baseline2D, xy, z)

    base_map = func_baseline2D(XY, *popt)
    std_bmap = func_baseline2D(XY, *(np.diag(pcov)))

    return base_map, np.abs(std_bmap)

def func_baseline2D(xy, A, B, C, E, F, G, H, I, J, K, L):
    y = np.asfarray(xy[1])
    x = np.asfarray(xy[0])
    p_func = A + B * x * y + C * x + E * x ** 2 + F * y ** 2 + G * x ** 3 + H * y ** 3 + I * x ** 4 + J * y ** 4 + K * x ** 5 + L * y ** 5
    return p_func

def func_gfilter(data, centerfh, centerfv, f_range, sigma_gfilter):
    '''
    2D gaussian filter
    :param 2D array, filter position, filter range and sigma of Gfilter
    :return: 2D arrayfilter and sigma of Gfilter
    '''
    gfilterv = np.zeros(np.shape(data))
    gfilterh = np.zeros(np.shape(data))
    nl, nr = np.shape(data)
    X = np.arange(0, nr, 1)
    Y = np.arange(0, nl, 1)
    # Creating Filter for Horizontal/vertical fringes orientation
    if sigma_gfilter == 0:
        # sigma filter is a func of image dimensions and f_rqnge
        sigma_gfilter = (2 * f_range)
    for i in X:
        gfilterh[:, i] = np.exp(-np.square(Y - centerfh) / (2 * np.square(sigma_gfilter)))

    for i in Y:
        gfilterv[i, :] = np.exp(-np.square(X - centerfv) / (2 * np.square(sigma_gfilter)))

    if centerfh == 0:
        gfilter = gfilterv

    elif centerfv == 0:
        gfilter = gfilterh

    else:
        gfilter = gfilterv * gfilterh

    return gfilter, int(sigma_gfilter)

#Function to detect fringes orientation
def func_cfilter(data, centerfh, centerfv, oppfx, oppfy):

    if oppfx == False: fx = 0
    else: fx = -1
    if oppfy == False: fy = 0
    else: fy = -1

    nl, nr = np.shape(data)

    summapv = np.sum(data, axis=0) - als(np.sum(data, axis=0))
    summaph = np.sum(data, axis=1) - als(np.sum(data, axis=1))


    try:
        fpv = find_peaks(summapv, height=0.5 * np.max(summapv), width=1)[0]
        # fwhm of peaks for sigma of gaussian filter
        pwv = (peak_widths(summapv, fpv, rel_height=0.5)[0])
    except:
        fpv = np.array([0,0])
        pwv = np.array([0,0])
    if len(fpv) == 0:
        fpv = np.array([0, 0])
        pwv = np.array([0, 0])

    try:
        fph = find_peaks(summaph, height=0.5 * np.max(summaph), width=1)[0]
        # Range of gaussian filter
        pwh = (peak_widths(summaph, fph, rel_height=0.5)[0])

    except:
        fph = np.array([0,0])
        pwh = np.array([0,0])
    if len(fph) == 0:
        fph = np.array([0, 0])
        pwh = np.array([0, 0])

    if centerfh == 0:
        centerfh = fph[fx]
    if centerfv == 0:
        centerfv = fpv[fy]

    if centerfh == 0: fang_rad = np.pi/2
    else: fang_rad = np.tan(centerfv / centerfh)

    fang_deg = (fang_rad) * 180 / np.pi


    fh_range = pwh[fx]
    fv_range = pwv[fy]
    f_range = np.max([fh_range, fv_range,2])

    return summaph,fph, summapv, fpv, centerfh, centerfv, f_range, fang_deg

def als(data, lam=5e2, p=0.1, itermax=10):
    r"""
    Implements an Asymmetric Least Squares Smoothing
    baseline correction algorithm (P. Eilers, H. Boelens 2005)
    Inputs:
        y:
            input data (i.e. chromatogram of spectrum)
        lam:
            parameter that can be adjusted by user. The larger lambda is,
            the smoother the resulting background, z
        p:
            wheighting deviations. 0.5 = symmetric, <0.5: negative
            deviations are stronger suppressed
        itermax:
            number of iterations to perform
    Output:
        the fitted background vector

    """
    L = len(data)
#  D = sparse.csc_matrix(np.diff(np.eye(L), 2))
    D = sparse.eye(L, format='csc')
    D = D[1:] - D[:-1]  # numpy.diff( ,2) does not work with sparse matrix. This is a workaround.
    D = D[1:] - D[:-1]
    D = D.T
    w = np.ones(L)
    for i in range(itermax):
        W = sparse.diags(w, 0, shape=(L, L))
        Z = W + lam * D.dot(D.T)
        z = spsolve(Z, w * data)
        w = p * (data > z) + (1 - p) * (data < z)
    return z
