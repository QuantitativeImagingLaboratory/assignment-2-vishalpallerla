# For this part of the assignment, You can use inbuilt functions to compute the fourier transform
# You are welcome to use fft that are available in numpy and opencv

import numpy as np
import math
import cv2


class Filtering:
    image = None
    filter = None
    cutoff = None
    order = None

    def __init__(self, image, filter_name, cutoff, order=0):
        """initializes the variables frequency filtering on an input image
        takes as input:
        image: the input image
        filter_name: the name of the mask to use
        cutoff: the cutoff frequency of the filter
        order: the order of the filter (only for butterworth
        returns"""
        self.image = image
        if filter_name == 'ideal_l':
            self.filter = self.get_ideal_low_pass_filter
        elif filter_name == 'ideal_h':
            self.filter = self.get_ideal_high_pass_filter
        elif filter_name == 'butterworth_l':
            self.filter = self.get_butterworth_low_pass_filter
        elif filter_name == 'butterworth_h':
            self.filter = self.get_butterworth_high_pass_filter
        elif filter_name == 'gaussian_l':
            self.filter = self.get_gaussian_low_pass_filter
        elif filter_name == 'gaussian_h':
            self.filter = self.get_gaussian_high_pass_filter

        self.cutoff = cutoff
        self.order = order

    def get_ideal_low_pass_filter(self, shape, cutoff):
        """Computes a Ideal low pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the ideal filter
        returns a ideal low pass mask"""

        [h, w] = shape

        mask_image = np.zeros((h, w))
        for i in range(h):
            for j in range(w):
                distance = math.sqrt((i - (h / 2)) * (i - (h / 2)) + (j - (w / 2)) * (j - (w / 2)))
                if distance <= cutoff:
                    mask_image[i][j] = 1
                else:
                    mask_image[i][j] = 0
        # cv2.imshow("IdealLowpass", mask_image)
        # cv2.waitKey(0)
        return mask_image

    def get_ideal_high_pass_filter(self, shape, cutoff):
        """Computes a Ideal high pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the ideal filter
        returns a ideal high pass mask"""

        # Hint: May be one can use the low pass filter function to get a high pass mask
        mask_image = 1 - self.get_ideal_low_pass_filter(shape, cutoff)
        return mask_image

    def get_butterworth_low_pass_filter(self, shape, cutoff):
        """Computes a butterworth low pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the butterworth filter
        order: the order of the butterworth filter
        returns a butterworth low pass mask"""

        [h, w] = shape
        mask_image = np.zeros((h, w))
        for i in range(h):
            for j in range(w):
                distance = math.sqrt((i - (h / 2)) ** 2 + ((j - (w / 2)) ** 2))
                mask_image[i][j] = 1 / (1 + ((distance / cutoff) ** (2 * self.order)))

        # cv2.imshow("ButterLow", mask_image)
        # cv2.waitKey(0)
        return mask_image

    def get_butterworth_high_pass_filter(self, shape, cutoff):
        """Computes a butterworth high pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the butterworth filter
        order: the order of the butterworth filter
        returns a butterworth high pass mask"""

        # Hint: May be one can use the low pass filter function to get a high pass mask
        mask_image = 1 - self.get_butterworth_low_pass_filter(shape, cutoff)
        return mask_image

    def get_gaussian_low_pass_filter(self, shape, cutoff):
        """Computes a gaussian low pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the gaussian filter (sigma)
        returns a gaussian low pass mask"""

        [h, w] = shape

        mask_image = np.zeros((h, w))
        for i in range(h):
            for j in range(w):
                distance = math.sqrt((i - (h/ 2)) ** 2 + ((j - (w / 2)) ** 2))
                mask_image[i][j] = 1 / math.exp((distance * distance) / (2 * cutoff * cutoff))

        # cv2.imshow("Gauss Low", mask_image)
        # cv2.waitKey(0)
        return mask_image

    def get_gaussian_high_pass_filter(self, shape, cutoff):
        """Computes a gaussian high pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the gaussian filter (sigma)
        returns a gaussian high pass mask"""

        # Hint: May be one can use the low pass filter function to get a high pass mask
        mask_image = 1 - self.get_gaussian_low_pass_filter(shape, cutoff)
        return mask_image

    def post_process_image(self, image):
        """Post process the image to create a full contrast stretch of the image
        takes as input:
        image: the image obtained from the inverse fourier transform
        return an image with full contrast stretch
        -----------------------------------------------------
        1. Full contrast stretch (fsimage)
        2. take negative (255 - fsimage)
        """
        c_min = np.min(image)
        c_max = np.max(image)
        new_min = 0
        new_max = 255
        stretch_image = np.zeros((np.shape(image)), dtype=np.uint8)
        for i in range(0, image.shape[0]):
            for j in range(0, image.shape[1]):
                stretch_image[i][j] = (image[i][j] - c_min) * ((new_max - new_min) / (c_max - c_min)) + new_min

        return stretch_image

    def filtering(self):
        """Performs frequency filtering on an input image
        returns a filtered image, magnitude of DFT, magnitude of filtered DFT        
        ----------------------------------------------------------
        You are allowed to used inbuilt functions to compute fft
        There are packages available in numpy as well as in opencv
        Steps:
        1. Compute the fft of the image
        2. shift the fft to center the low frequencies
        3. get the mask (write your code in functions provided above) the functions can be called by self.filter(shape, cutoff, order)
        4. filter the image frequency based on the mask (Convolution theorem)
        5. compute the inverse shift
        6. compute the inverse fourier transform
        7. compute the magnitude
        8. You will need to do a full contrast stretch on the magnitude and depending on the algorithm you may also need to
        take negative of the image to be able to view it (use post_process_image to write this code)
        Note: You do not have to do zero padding as discussed in class, the inbuilt functions takes care of that
        filtered image, magnitude of DFT, magnitude of filtered DFT: Make sure all images being returned have grey scale full contrast stretch and dtype=uint8 
        """




        given_image = self.image
        s = given_image.shape

        fft_image = np.fft.fft2(self.image)
        shift_image = np.fft.fftshift(fft_image)
        dft_image = np.uint8(np.log(np.absolute(shift_image))*10)

        mask = self.filter(s, self.cutoff)
        filter_image = shift_image * mask
        filter_finalimg = np.log(np.absolute(filter_image))*10

        ishift_image = np.fft.ifftshift(filter_image)
        ifft_image = np.fft.ifft2(ishift_image)
        mag_image = np.absolute(ifft_image)
        f = self.post_process_image(mag_image)

        return [f, dft_image, filter_finalimg]





