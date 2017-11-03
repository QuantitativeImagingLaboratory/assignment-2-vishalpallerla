# Report

1. DFT: 
- Forward Fourier Transform: 
- Height and width of the image are taken using the shape property.
- A result matrix with complex data type is generated using all zeroes with new height and weight so that it we can use it to store the calculated forward fourier transform values.
- Running through each pixel in the result matrix from left to right and top to bottom calculate the FFT using the input matrix values according to the below mentioned formula and assign them resultant matrix.
- verified the resultant matrix with the inbuilt fft2 function of numpy and it matched.

- Inverse Fourier Transform:
- Height and width of the image are taken using the shape property.
- A result matrix with complex data type is generated using all zeroes with new height and weight so that it we can use it to store the calculated inverse fourier transform values.
- Running through each pixel in the result matrix from left to right and top to bottom calculate the IFT using the input matrix values according to the below mentioned formula and assign them to the resultant matrix.
- Checked the resultant matrix with inbuilt ifft2 function of numpy and it did not match with the result (searched about the reason and got to know that inbuilt functions normalize the values where as I did not and hence the difference)

- Discrete Cosine Transform:
- Height and width of the image are taken using the shape property.
- A result matrix is generated using all zeroes with new height and weight so that it we can use it to store the calculated discrete cosine transform values.
- Running through each pixel in the result matrix from left to right and top to bottom calculate the DCT using the input matrix values according to the below mentioned formula and assign them to the resultant matrix.
- No complex part, so used normal matrix

- Magnitude:
- Height and width of the image are taken using the shape property.
- A result matrix is generated using all zeroes with new height and weight so that it we can use it to store the calculated magnitude values.
- Running through each pixel in the result matrix from left to right and top to bottom calculate the magnitude using the input matrix values according to the formula and assign them to the resultant matrix.

2. Filtering:
- Ideal Low Pass Filter:
- Allows low frequencies than cut off and discards high frequencies than cut off.
- Height and width of the image are taken using the shape property.
- A result matrix is generated using all zeroes with new height and weight so that it we can use it to store the calculated values.
- Running through each pixel in the result matrix from left to right and top to bottom calculate the low pass using the input matrix values according to the formula and assign them to the resultant matrix.

- Observations:
- Output of Low pass filter was observed to be smoothened(blurred)
- As we increase the mask size, image gets more smoothened.
- Hence it is called blurring filter.
- Ringing effect is observed

- Ideal High Pass Filter:
- Opposite to Ideal Low pass filter and hence 1 � Ideal Low Pass Filter output.
- Observations: 
- High pass filter output is observed to have edge content visible.
- As the size of the mask grows, more edge content is increased
- Hence, can be used for edge detection.

- Butterworth Low Pass Filter:
- Height and width of the image are taken using the shape property.
- A result matrix is generated using all zeroes with new height and weight so that it we can use it to store the calculated values.
- Running through each pixel in the result matrix from left to right and top to bottom calculate the butterworth  low pass using the input matrix values according to the formula and assign them to the resultant matrix.
- mask_image[i][j] = 1 / (1 + ((distance / cutoff) ** (2 * self.order)))
- Observations:
- If we keep the order constant and increase the cut off, the output gets clearer.
- Ringing effect is avoided.

- Butterworth High Pass Filter:
- Opposite to Butterworth LowPass filter

- Gaussian Low Pass Filter:
- Transition is more smooth compared to Ideal Low Pass filter
- A low pass Gaussian filter is used to connect broken text

- Gaussian High Pass Filter:
- Transition is more smooth compared to Ideal High Pass filter
- Edges and fine detail in images are associated with high frequency components 

- Performing Frequency Filtering:

- Frequency filtering is done according to the given steps as follows:

-  1. Compute the fft of the image
2. shift the fft to center the low frequencies
3. get the mask by self.filter(shape, cutoff, order)
4. filter the image frequency based on the mask
5. compute the inverse shift
6. compute the inverse fourier transform
7. compute the magnitude
8. full contrast stretch on the magnitude is done using  post_process_image function
- 

