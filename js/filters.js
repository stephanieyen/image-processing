"use strict";

const Filters = {};

////////////////////////////////////////////////////////////////////////////////
// General utility functions
////////////////////////////////////////////////////////////////////////////////

// Hardcoded Pi value
// const pi = 3.14159265359;
const pi = Math.PI;

// Constrain val to the range [min, max]
function clamp(val, min, max) {
    /* Shorthand for:
    * if (val < min) {
    *   return min;
    * } else if (val > max) {
    *   return max;
    * } else {
    *   return val;
    * }
    */
    return val < min ? min : val > max ? max : val;
}

// Extract vertex coordinates from a URL string
function stringToCoords(vertsString) {
    const centers = [];
    const coordStrings = vertsString.split("x");
    for (let i = 0; i < coordStrings.length; i++) {
        const coords = coordStrings[i].split("y");
        const x = parseInt(coords[0]);
        const y = parseInt(coords[1]);
        if (!isNaN(x) && !isNaN(y)) {
            centers.push({ x: x, y: y });
        }
    }

    return centers;
}

// Blend scalar start with scalar end. Note that for image blending,
// end would be the upper layer, and start would be the background
function blend(start, end, alpha) {
    return start * (1 - alpha) + end * alpha;
}

// ----------- STUDENT CODE BEGIN ------------
// ----------- Our reference solution uses 72 lines of code.
// ----------- STUDENT CODE END ------------

////////////////////////////////////////////////////////////////////////////////
// Filters
////////////////////////////////////////////////////////////////////////////////

// You've already implemented this in A0! Feel free to copy your code into here
Filters.fillFilter = function(image, color) {
    image.fill(color);

    return image;
};

// You've already implemented this in A0! Feel free to copy your code into here
Filters.brushFilter = function(image, radius, color, vertsString) {
    // centers is an array of (x, y) coordinates that each defines a circle center
    const centers = stringToCoords(vertsString);

    // draw a filled circle centered at every location in centers[].
    // radius and color are specified in function arguments.
    for (let i = 0; i < centers.length; i++) {
        let centerX = centers[i].x;
        let centerY = centers[i].y;
        for (let x = centerX - radius; x <= centerX + radius; x++) {
          // Use circle equation: x^2 + y^2 = r^2
          let yDiff = Math.sqrt(Math.abs( Math.pow(radius, 2) - Math.pow(Math.abs(centerX - x), 2) ));
          for (let y = centerY - yDiff; y <= centerY + yDiff; y++) {
            // Change color of pixel at (x, y)
            image.setPixel(x, y, color);
          }
        }
    }

    return image;
};

// You've already implemented this in A0! Feel free to copy your code into here
Filters.softBrushFilter = function(image, radius, color, alpha_at_center, vertsString) {
    // centers is an array of (x, y) coordinates that each defines a circle center
    const centers = stringToCoords(vertsString);

    // draw a filled circle with opacity equals to alpha_at_center at the center of each circle
    // the opacity decreases linearly along the radius and becomes zero at the edge of the circle
    // Save RGB values to overlay
    let overlayR = color.data[0];
    let overlayG = color.data[1];
    let overlayB = color.data[2];

    for (let i = 0; i < centers.length; i++) {
        let centerX = centers[i].x;
        let centerY = centers[i].y;
        for (let x = centerX - radius; x <= centerX + radius; x++) {
        // Use circle equation: x^2 + y^2 = r^2
        let yDiff = Math.sqrt(Math.abs( Math.pow(radius, 2) - Math.pow(Math.abs(centerX - x), 2) ));
        for (let y = centerY - yDiff; y <= centerY + yDiff; y++) {        
            // Save current pixel at (x, y)
            let currPixel = image.getPixel(x, y);

            // At pixel (x, y), calculate the alpha value to overlay based on its distance from center of circle
            let distance = Math.sqrt((x-centerX)*(x-centerX) + (y-centerY)*(y-centerY));
            let overlayAlpha = -(alpha_at_center/radius)*distance + alpha_at_center;

            // Use the calculated alpha to calculate the pixel's new RGB values
            color.data[0] = overlayAlpha*overlayR + currPixel.data[0]*(1 - overlayAlpha);       
            color.data[1] = overlayAlpha*overlayG + currPixel.data[1]*(1 - overlayAlpha);    
            color.data[2] = overlayAlpha*overlayB + currPixel.data[2]*(1 - overlayAlpha);   

            image.setPixel(x, y, color);
        }
        }
    }

    return image;
};

// Ratio is a value in the domain [-1, 1]. When ratio is < 0, linearly blend the image
// with black. When ratio is > 0, linearly blend the image with white. At the extremes
// of -1 and 1, the image should be completely black and completely white, respectively.
Filters.brightnessFilter = function(image, ratio) {
    let alpha, dirLuminance;
    if (ratio < 0.0) {
        alpha = 1 + ratio;
        dirLuminance = 0; // blend with black (decreasing brightness)
    } else {
        alpha = 1 - ratio;
        dirLuminance = 1; // blend with white (increasing brightness)
    }

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            pixel.data[0] = alpha * pixel.data[0] + (1 - alpha) * dirLuminance;
            pixel.data[1] = alpha * pixel.data[1] + (1 - alpha) * dirLuminance;
            pixel.data[2] = alpha * pixel.data[2] + (1 - alpha) * dirLuminance;

            image.setPixel(x, y, pixel);
        }
    }

    return image;
};

// Reference at this:
//      https://en.wikipedia.org/wiki/Image_editing#Contrast_change_and_brightening
// value = (value - 0.5) * (tan ((contrast + 1) * PI/4) ) + 0.5;
// Note that ratio is in the domain [-1, 1]
Filters.contrastFilter = function(image, ratio) {
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            // For each color channel, apply formula and clamp pixel
            pixel.data[0] = (pixel.data[0] - 0.5) * (Math.tan ((ratio + 1) * pi/4) ) + 0.5;
            pixel.data[0] = clamp(pixel.data[0], 0, 1);
            pixel.data[1] = (pixel.data[1] - 0.5) * (Math.tan ((ratio + 1) * pi/4) ) + 0.5;
            pixel.data[1] = clamp(pixel.data[1], 0, 1);
            pixel.data[2] = (pixel.data[2] - 0.5) * (Math.tan ((ratio + 1) * pi/4) ) + 0.5;
            pixel.data[2] = clamp(pixel.data[2], 0, 1);

            image.setPixel(x, y, pixel);
        }
    }

    return image;
};

// Note that the argument here is log(gamma)
Filters.gammaFilter = function(image, logOfGamma) {
    const gamma = Math.exp(logOfGamma);

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            pixel.data[0] = Math.pow(pixel.data[0], gamma);
            pixel.data[1] = Math.pow(pixel.data[1], gamma);
            pixel.data[2] = Math.pow(pixel.data[2], gamma);

            image.setPixel(x, y, pixel);
        }
    }

    return image;
};

/*
* The image should be perfectly clear up to innerRadius, perfectly dark
* (black) at outerRadius and beyond, and smoothly increase darkness in the
* circular ring in between. Both are specified as multiples of half the length
* of the image diagonal (so 1.0 is the distance from the image center to the
* corner).
*
* Note that the vignette should still form a perfect circle!
*/
Filters.vignetteFilter = function(image, innerR, outerR) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 17 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('vignetteFilter is not implemented yet');
};

/*
* You will want to build a normalized CDF of the L channel in the image.
*/
Filters.histogramEqualizationFilter = function(image) {
    const L = 256; // number of gray levels ("bins")
    const N = image.width * image.height; // number of pixels

    // Initialize arrays, indices are each i-th gray level 
    let pdf = []; // probability density function
    let cdf = []; // cumulative distribution function
    let numPixels = []; // number of pixels in each grey level
    for (let i = 0; i < L; i++) {
        pdf[i] = 0;
        cdf[i] = 0;
        numPixels[i] = 0;
    }

    // Build PDF
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            // Get the pixel's luminance value by converting into HSL space
            const pixelHsl = pixel.rgbToHsl();
            const luminance = pixelHsl.data[2];
            const i = Math.round(luminance*255); // the i-th grey level that corresponds to the luminance value
            numPixels[i]++; // increment number of pixels in this grey level

            // Make PDF 
            pdf[i] = numPixels[i];
        }
    }

    // Build CDF (normalized)
    cdf[0] = pdf[0];
    for (let j = 1; j <= L-1; j++) {
        cdf[j] = pdf[j] + cdf[j-1];
    }

    for (let i = 0; i < cdf.length; i++) {
        cdf[i] /= N;
    }

    // Iterate through image and adjust luminance value
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            // Find new luminance value for pixel by mapping original luminance to corresponding value in CDF
            const pixelHsl = pixel.rgbToHsl();
            const luminance = pixelHsl.data[2];
            const i = Math.round(luminance*255); // the i-th grey level that corresponds to the luminance value

            pixelHsl.data[2] = cdf[i]; // new luminance value

            const newPixel = pixelHsl.hslToRgb();
            image.setPixel(x, y, newPixel);
        }
    }

    return image;
};

///////////////////////////////////////////
// Color Operations
///////////////////////////////////////////

// Set each pixel in the image to its luminance
Filters.grayscaleFilter = function(image) {
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);
            const luminance = 0.2126 * pixel.data[0] + 0.7152 * pixel.data[1] + 0.0722 * pixel.data[2];
            pixel.data[0] = luminance;
            pixel.data[1] = luminance;
            pixel.data[2] = luminance;

            image.setPixel(x, y, pixel);
        }
    }

    return image;
};

// Adjust each channel in each pixel by a fraction of its distance from the average
// value of the pixel (luminance).
// See: http://www.graficaobscura.com/interp/index.html
Filters.saturationFilter = function(image, ratio) {
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);
            const luminance = 0.2126 * pixel.data[0] + 0.7152 * pixel.data[1] + 0.0722 * pixel.data[2];

            // For each color channel, apply formula and clamp pixel
            pixel.data[0] = pixel.data[0] + (pixel.data[0] - luminance) * ratio;
            pixel.data[0] = clamp(pixel.data[0], 0, 1);
            pixel.data[1] = pixel.data[1] + (pixel.data[1] - luminance) * ratio;
            pixel.data[1] = clamp(pixel.data[1], 0, 1);
            pixel.data[2] = pixel.data[2] + (pixel.data[2] - luminance) * ratio;
            pixel.data[2] = clamp(pixel.data[2], 0, 1);

            image.setPixel(x, y, pixel);
        }
    }

    return image;
};

// Apply the Von Kries method: convert the image from RGB to LMS, divide by
// the LMS coordinates of the white point color, and convert back to RGB.
Filters.whiteBalanceFilter = function(image, white) {
    // Save the color of "white" in LMS
    let whiteXyz = white.rgbToXyz();
    let whiteLms = whiteXyz.xyzToLms();

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            // Save each pixel in the LMS color space
            const pixel = image.getPixel(x, y);
            let pixelXyz = pixel.rgbToXyz();
            let pixelLms = pixelXyz.xyzToLms();

            // For each color channel, divide by the color of "white" in LMS and clamp pixel
            pixelLms.data[0] = pixelLms.data[0] / whiteLms.data[0];
            pixelLms.data[0] = clamp(pixelLms.data[0], 0, 1);
            pixelLms.data[1] = pixelLms.data[1] / whiteLms.data[1];
            pixelLms.data[1] = clamp(pixelLms.data[1], 0, 1);
            pixelLms.data[2] = pixelLms.data[2] / whiteLms.data[2];
            pixelLms.data[2] = clamp(pixelLms.data[2], 0, 1);

            // Convert back to RGB color space
            pixelXyz = pixelLms.lmsToXyz();
            image.setPixel(x, y, pixelXyz.xyzToRgb());
        }
    }

    return image;
};

// This is similar to the histogram filter, except here you should take the
// the CDF of the L channel in one image and map it to another
Filters.histogramMatchFilter = function(image, refImg) {
    const L = 256; // number of gray levels ("bins")

    // --- Compute lightness histograms of the two images ---

    // Initialize arrays, indices are each i-th gray level 
    // Probability density functions
    let refImgPdf = []; 
    let imagePdf = [];
    // Cumulative distribution functions
    let refImgCdf = []; 
    let imageCdf = [];
    // Number of pixels in each grey level
    let refImgNumPixels = [];
    let imageNumPixels = [];
    // Histogram mapping function 
    
    for (let i = 0; i < L; i++) {
        refImgPdf[i] = 0;
        refImgCdf[i] = 0;
        refImgNumPixels[i] = 0;
        imagePdf[i] = 0;
        imageCdf[i] = 0;
        imageNumPixels[i] = 0;
    }

    // --- Compute F1, lightness histogram of refImg --- 
    const refImgN = refImg.width * refImg.height; // number of pixels

    // Build PDF
    for (let x = 0; x < refImg.width; x++) {
        for (let y = 0; y < refImg.height; y++) {
            const pixel = refImg.getPixel(x, y);

            // Get the pixel's luminance value by converting into HSL space
            const pixelHsl = pixel.rgbToHsl();
            const luminance = pixelHsl.data[2];
            const i = Math.round(luminance*255); // the i-th grey level that corresponds to the luminance value
            refImgNumPixels[i]++; // increment number of pixels in this grey level

            // Make PDF 
            refImgPdf[i] = refImgNumPixels[i];
        }
    }

    // Build CDF (normalized)
    refImgCdf[0] = refImgPdf[0];
    for (let j = 1; j <= L-1; j++) {
        refImgCdf[j] = refImgPdf[j] + refImgCdf[j-1];
    }

    for (let i = 0; i < refImgCdf.length; i++) {
        refImgCdf[i] /= refImgN;
    }

    // --- Compute F2, lightness histogram of image --- 
    const imageN = image.width * image.height; // number of pixels

    // Build PDF
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            // Get the pixel's luminance value by converting into HSL space
            const pixelHsl = pixel.rgbToHsl();
            const luminance = pixelHsl.data[2];
            const i = Math.round(luminance*255); // the i-th grey level that corresponds to the luminance value
            imageNumPixels[i]++; // increment number of pixels in this grey level

            // Make PDF 
            imagePdf[i] = imageNumPixels[i];
        }
    }

    // Build CDF (normalized)
    imageCdf[0] = imagePdf[0];
    for (let j = 1; j <= L-1; j++) {
        imageCdf[j] = imagePdf[j] + imageCdf[j-1];
    }

    for (let i = 0; i < imageCdf.length; i++) {
        imageCdf[i] /= imageN;
    }

    // --- For each i-th gray level G1 (in ref image), find G2 for which F1(G1) = F2(G2) ---

    let minDiff = L;
    for (let x = 0; x < L; x++) { // all gray levels of original image
        for (let i = 0; i < L; i++) { // all gray levels of reference image
            const diff = Math.abs(imageCdf[x] - refImgCdf[i]);

            // Set each pixel's gray level from x to their corresponding i (for which CDF values are minimized)
            if (diff < minDiff) {
                minDiff = Math.min(minDiff, diff);
                imageCdf[x] = refImgCdf[i];
            }
        }
    }
    
    // --- Apply the mapping of G2 --> G1 to match image's histogram to refImg's histogram ---

    // Iterate through source image and adjust luminance value
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            // Find new luminance value for pixel by mapping original luminance to corresponding value in CDF
            const pixelHsl = pixel.rgbToHsl();
            const luminance = pixelHsl.data[2];
            const i = Math.round(luminance*255); // the i-th grey level that corresponds to the luminance value

            pixelHsl.data[2] = imageCdf[i]; // new luminance value

            const newPixel = pixelHsl.hslToRgb();
            image.setPixel(x, y, newPixel);
        }
    }

    return image;
};

///////////////////////////////////////////
// Filter Operations
///////////////////////////////////////////

// Convolve the image with a gaussian filter.
// NB: Implement this as a seperable gaussian filter
Filters.gaussianFilter = function(image, sigma) {
    // note: this function needs to work in a new copy of the image
    //       to avoid overwriting original pixel values needed later
    // create a new image with the same size as the input image
    let newImg = image.createImg(image.width, image.height);
    // the filter window will be [-winR, winR] for a total diameter of roughly Math.round(3*sigma)*2+1;
    const winR = Math.round(sigma * 3);

    // ----------- STUDENT CODE BEGIN ------------

    // Pre-compute 1D Gaussian kernel 
    let weights = []; // indices reflect distance from center of kernel
    let weightsSum = 0;

    // Calculate the weight at each location in the kernel based on its distance from the center
    for (let i = -winR; i <= winR; i++) {
        weights[i] = Math.exp((-1*i*i)/(2*sigma*sigma));
        weightsSum += weights[i];
    }

    // Normalize kernel values (so they sum to 1) by dividing each value by the sum
    for (let i = -winR; i <= winR; i++) {
        weights[i] /= weightsSum;
    }

    // Linear separation optimization
    // 1. Convolve original image with 1D vertical filter
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);
            pixel.data[0] = 0;
            pixel.data[1] = 0;
            pixel.data[2] = 0;

            // Iterate through the neighboring vertical pixels
            for (let i = -winR; i <= winR; i++) {
                const neighborPixel = image.getPixel(x, clamp(y+i, 0, image.height-1));

                // For each color channel, sum the product of the kernel weights and corresponding neighboring pixel values
                pixel.data[0] += weights[i] * neighborPixel.data[0];
                pixel.data[1] += weights[i] * neighborPixel.data[1];
                pixel.data[2] += weights[i] * neighborPixel.data[2];
            }

            newImg.setPixel(x, y, pixel);
        }
    }

    image = newImg;

    // 2. Convolve resulting image with the 1D horizontal filter to get final image
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);
            pixel.data[0] = 0;
            pixel.data[1] = 0;
            pixel.data[2] = 0;

            // Iterate through the neighboring vertical pixels
            for (let i = -winR; i <= winR; i++) {
                const neighborPixel = image.getPixel(clamp(x+i, 0, image.width-1), y);

                // For each color channel, sum the product of the kernel weights and corresponding neighboring pixel values
                pixel.data[0] += weights[i] * neighborPixel.data[0];
                pixel.data[1] += weights[i] * neighborPixel.data[1];
                pixel.data[2] += weights[i] * neighborPixel.data[2];
            }

            newImg.setPixel(x, y, pixel);
        }
    }

    return newImg;
};

/*
* First the image with the edge kernel and then add the result back onto the
* original image.
*/
Filters.sharpenFilter = function(image) {
    // Make a copy of the image to avoid overwriting original pixel values needed later
    let newImg = image.createImg(image.width, image.height);

    // Pre-compute sharpen kernel
    let weights = new Array(3);
    for (let i = 0; i < weights.length; i++) {
        weights[i] = new Array(3).fill(-1);
        if (i == 1) 
            weights[i][1] = 9;
    }

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            pixel.data[0] = 0;
            pixel.data[1] = 0;
            pixel.data[2] = 0;

            // Each pixel becomes the sum of the product of the kernel weights and corresponding neighboring pixel values
            for (let i = -1; i <= 1; i++) {
                for (let j = -1; j <= 1; j++) {
                    const neighborPixel = image.getPixel(clamp(x+i, 0, image.width-1), clamp(y+j, 0, image.height-1));

                    // For each color channel, sum the product of the kernel weights and corresponding neighboring pixel values
                    pixel.data[0] += weights[i+1][j+1] * neighborPixel.data[0];
                    pixel.data[1] += weights[i+1][j+1] * neighborPixel.data[1];
                    pixel.data[2] += weights[i+1][j+1] * neighborPixel.data[2];
                }
            }

            newImg.setPixel(x, y, pixel);
        }
    }

    return newImg;
};

/*
* Convolve the image with the edge kernel from class. You might want to define
* a convolution utility that convolves an image with some arbitrary input kernel
*
* For this filter, we recommend inverting pixel values to enhance edge visualization
*/
Filters.edgeFilter = function(image) {
    // Make a copy of the image to avoid overwriting original pixel values needed later
    let newImg = image.createImg(image.width, image.height);

    // Pre-compute edge kernel
    let weights = new Array(3);
    for (let i = 0; i < weights.length; i++) {
        weights[i] = new Array(3).fill(-1);
        if (i == 1) 
            weights[i][1] = 8;
    }

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            pixel.data[0] = 0;
            pixel.data[1] = 0;
            pixel.data[2] = 0;

            // Each pixel becomes the sum of the product of the kernel weights and corresponding neighboring pixel values
            for (let i = -1; i <= 1; i++) {
                for (let j = -1; j <= 1; j++) {
                    const neighborPixel = image.getPixel(clamp(x+i, 0, image.width-1), clamp(y+j, 0, image.height-1));

                    // For each color channel, sum the product of the kernel weights and corresponding neighboring pixel values
                    pixel.data[0] += weights[i+1][j+1] * neighborPixel.data[0];
                    pixel.data[1] += weights[i+1][j+1] * neighborPixel.data[1];
                    pixel.data[2] += weights[i+1][j+1] * neighborPixel.data[2];
                }
            }

            // Invert pixel
            pixel.data[0] = 1 - pixel.data[0];
            pixel.data[1] = 1 - pixel.data[1];
            pixel.data[2] = 1 - pixel.data[2];

            newImg.setPixel(x, y, pixel);
        }
    }

    return newImg;
};

// Set a pixel to the median value in its local neighbor hood. You might want to
// apply this seperately to each channel.
Filters.medianFilter = function(image, winR) {
    // winR: the window will be  [-winR, winR];
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 36 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('medianFilter is not implemented yet');
    return image;
};

// Apply a bilateral filter to the image. You will likely want to reference
// precept slides, lecture slides, and the assignments/examples page for help.
Filters.bilateralFilter = function(image, sigmaR, sigmaS) {
    // reference: https://en.wikipedia.org/wiki/Bilateral_filter
    // we first compute window size and preprocess sigmaR
    const winR = Math.round((sigmaR + sigmaS) * 1.5);
    sigmaR = sigmaR * (Math.sqrt(2) * winR);

    // ----------- STUDENT CODE BEGIN ------------

    // Make a copy of the image to avoid overwriting original pixel values needed later
    let newImg = image.createImg(image.width, image.height);

    // Initialize kernel
    let weights = new Array(winR*2+1);
    for (let i = 0; i < weights.length; i++) {
        weights[i] = new Array(winR*2+1);
    }

    // Iterate through all pixels
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            let weightsSum = 0;

            // Compute bilateral kernel at this pixel based on both spatial distance and color distance
            for (let i = -winR; i <= winR; i++) {
                for (let j = -winR; j <= winR; j++) {
                    const pixelToCompare = image.getPixel(x+i, y+j);

                    const spatialDist = (-1 * (Math.pow(i, 2) + Math.pow(j, 2))) / (2*sigmaS*sigmaS);
                    const colorDist = (-1 * (Math.pow(pixelToCompare.data[0]*255 - pixel.data[0]*255, 2) 
                                    + Math.pow(pixelToCompare.data[1]*255 - pixel.data[1]*255, 2) 
                                    + Math.pow(pixelToCompare.data[2]*255 - pixel.data[2]*255, 2))) / (2*sigmaR*sigmaR);

                    weights[i+winR][j+winR] = Math.exp(spatialDist + colorDist);
                    weightsSum += weights[i+winR][j+winR];
                }
            }

            pixel.data[0] = 0;
            pixel.data[1] = 0;
            pixel.data[2] = 0;

            // Each pixel becomes the sum of the product of the kernel weights and corresponding neighboring pixel values
            for (let i = -winR; i <= winR; i++) {
                for (let j = -winR; j <= winR; j++) {
                    // Normalize pixel weights
                    weights[i+winR][j+winR] /= weightsSum;

                    const neighborPixel = image.getPixel(clamp(x+i, 0, image.width-1), clamp(y+j, 0, image.height-1));

                    // For each color channel, sum the product of the kernel weights and corresponding neighboring pixel values
                    pixel.data[0] += weights[i+winR][j+winR] * neighborPixel.data[0];
                    pixel.data[1] += weights[i+winR][j+winR] * neighborPixel.data[1];
                    pixel.data[2] += weights[i+winR][j+winR] * neighborPixel.data[2];
                }
            }

            newImg.setPixel(x, y, pixel);
        }
    }

    return newImg;
};

///////////////////////////////////////////
// Dithering Operations
///////////////////////////////////////////

// Conver the image to binary
Filters.quantizeFilter = function(image) {
    // convert to grayscale
    image = Filters.grayscaleFilter(image);

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            for (let c = 0; c < 3; c++) {
                pixel.data[c] = Math.round(pixel.data[c]);
            }

            pixel.clamp();
            image.setPixel(x, y, pixel);
        }
    }

    return image;
};

// To apply random dithering, first convert the image to grayscale, then apply
// random noise, and finally quantize
Filters.randomFilter = function(image) {
    // convert to grayscale
    image = Filters.grayscaleFilter(image);

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            // Apply random noise before quantization
            const val = pixel.data[0] + Math.random() - 0.5;

            // Quantize pixel
            for (let c = 0; c < 3; c++) {
                pixel.data[c] = clamp(val, 0, 1);
            }
            
            image.setPixel(x, y, pixel); 
        }
    }

    return image;
};

// Apply the Floyd-Steinberg dither with error diffusion
Filters.floydFilter = function(image) {
    // convert to grayscale
    image = Filters.grayscaleFilter(image);

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            // Save neighboring pixels
            const pixelRight = image.getPixel(x+1, y);
            const pixelBottomLeft = image.getPixel(x-1, y+1);
            const pixelBottom = image.getPixel(x, y+1);
            const pixelBottomRight = image.getPixel(x+1, y+1);

            // For each color channel, quantize and perform error diffusion
            for (let c = 0; c < 3; c++) {
                // Compute quantization error (diff b/w original pixel and quantized pixel)
                const err = pixel.data[c] - Math.round(pixel.data[c]);

                // Quantize pixel
                pixel.data[c] = clamp(Math.round(pixel.data[c], 0, 1));

                // Spread quantization error over four neighboring pixels 
                pixelRight.data[c] += (7/16)*err;
                pixelBottomLeft.data[c] += (3/16)*err;
                pixelBottom.data[c] += (5/16)*err;
                pixelBottomRight.data[c] += (1/16)*err;
            }

            image.setPixel(x, y, pixel);
            image.setPixel(x+1, y, pixelRight); 
            image.setPixel(x-1, y+1, pixelBottomLeft); 
            image.setPixel(x, y+1, pixelBottom); 
            image.setPixel(x+1, y+1, pixelBottomRight); 
        }
    }
    
    return image;
};

// Apply ordered dithering to the image. We recommend using the pattern from the
// examples page and precept slides.
Filters.orderedFilter = function(image) {
    // convert to gray scale
    image = Filters.grayscaleFilter(image);

    // Pre-compute "Bayer 4" filter
    const M = 4;
    const weights = [
        [ 15, 7, 13, 5 ],
        [ 3, 11, 1, 9 ],
        [ 12, 4, 14, 6 ],
        [ 0, 8, 2, 10 ]
    ];

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            const i = x % M;
            const j = y % M;

            // Compute quantization error (using any color channel) and threshold values
            const err = pixel.data[0] - Math.floor(pixel.data[0]);
            const threshold = (weights[i][j] + 1) / (M*M + 1);

            // For each color channel, quantize and perform error diffusion
            for (let c = 0; c < 3; c++) {
                if (err > threshold) {
                    pixel.data[c] = Math.ceil(pixel.data[c]);
                } else {
                    pixel.data[c] = Math.floor(pixel.data[c]);
                }
            }

            image.setPixel(x, y, pixel);
        }
    }

    return image;
};

///////////////////////////////////////////
// Resampling Operations
///////////////////////////////////////////

// Implement bilinear and Gaussian sampling (in addition to the basic point sampling).
// This operation doesn't appear on GUI and should be used as a utility function.
// Call this function from filters that require sampling (e.g. scale, rotate)
Filters.samplePixel = function(image, x, y, mode) {
    if (mode === "bilinear") {
        
        x = clamp(x, 0.5, image.width - 1.5);
        y = clamp(y, 0.5, image.height - 1.5);

        const pixel = image.getPixel(x, y);

        // Get closest pixel coordinates in input image 
        const x1 = Math.floor(x);
        const x2 = Math.floor(x) + 1;
        const y1 = Math.floor(y) + 1;
        const y2 = Math.floor(y);
        const q11 = image.getPixel(x1, y1);
        const q12 = image.getPixel(x1, y2);
        const q21 = image.getPixel(x2, y1);
        const q22 = image.getPixel(x2, y2);

        // For each color channel of the current pixel, get the weighted average of the input neighboring pixels' values
        for (let c = 0; c < 3; c++) {
            // pixel.data[c] = clamp(Math.round(pixel.data[c], 0, 1));
            pixel.data[c] = (q11.data[c] * (x2 - x) * (y2 - y)
                            + q21.data[c] * (x - x1) * (y2 - y)
                            + q12.data[c] * (x2 - x) * (y - y1)
                            + q22.data[c] * (x - x1) * (y - y1)) /
                            ((x2 - x1) * (y2 - y1));
        }

        pixel.clamp();
        return pixel;

    } else if (mode === "gaussian") {

        x = clamp(x, 0.5, image.width - 1.5);
        y = clamp(y, 0.5, image.height - 1.5);

        const winR = 3; // equal to 3*sigma, assume sigma = 1
        
        // Pre-compute 1D Gaussian kernel
        let weights = []; // indices reflect distance from center of kernel
        let weightsSum = 0;
        
        // Calculate the weight at each location in the kernel based on its distance from the center
        for (let i = -winR; i <= winR; i++) {
            weights[i] = Math.exp((-1*i*i)/2);
            weightsSum += weights[i];
        }
        
        // Normalize kernel values (so they sum to 1) by dividing each value by the sum
        for (let i = -winR; i <= winR; i++) {
            weights[i] /= weightsSum;
        }
        
        // Linear separation optimization
        // 1. Convolve original image with 1D vertical filter
        const pixel = image.getPixel(x, y);
        pixel.data[0] = 0;
        pixel.data[1] = 0;
        pixel.data[2] = 0;

        // Iterate through the neighboring vertical pixels
        for (let i = -winR; i <= winR; i++) {
            const neighborPixel = image.getPixel(x, clamp(y+i, 0, image.height-1));

            // For each color channel, sum the product of the kernel weights and corresponding neighboring pixel values
            pixel.data[0] += weights[i] * neighborPixel.data[0];
            pixel.data[1] += weights[i] * neighborPixel.data[1];
            pixel.data[2] += weights[i] * neighborPixel.data[2];
        }
        
        // 2. Convolve resulting image with the 1D horizontal filter to get final image
        pixel.data[0] = 0;
        pixel.data[1] = 0;
        pixel.data[2] = 0;

        // Iterate through the neighboring vertical pixels
        for (let i = -winR; i <= winR; i++) {
            const neighborPixel = image.getPixel(clamp(x+i, 0, image.width-1), y);

            // For each color channel, sum the product of the kernel weights and corresponding neighboring pixel values
            pixel.data[0] += weights[i] * neighborPixel.data[0];
            pixel.data[1] += weights[i] * neighborPixel.data[1];
            pixel.data[2] += weights[i] * neighborPixel.data[2];
        }
        
        return pixel;

    } else {
        // point sampling
        y = Math.max(0, Math.min(Math.round(y), image.height - 1));
        x = Math.max(0, Math.min(Math.round(x), image.width - 1));
        return image.getPixel(x, y);
    }
};

// Translate the image by some x, y and using a requested method of sampling/resampling
Filters.translateFilter = function(image, x, y, sampleMode) {
    // Note: set pixels outside the image to RGBA(0,0,0,0)

    // Define start and end values for the final image
    const iStart = x;
    const iEnd = image.width + x;
    const jStart = y;
    const jEnd = image.height + y;

    // Make a copy of the image to avoid overwriting original pixel values needed later
    let newImg = image.createImg(image.width, image.height);

    for (let i = iStart; i < iEnd; i++) {
        for (let j = jStart; j < jEnd; j++) {
            // Look up the pixel you want in the input image (with resampling technique) and get its value
            const pixel = this.samplePixel(image, i - x, j - y, sampleMode);

            newImg.setPixel(i, j, pixel);
        }
    }

    return newImg;
};

// Scale the image by some ratio and using a requested method of sampling/resampling
Filters.scaleFilter = function(image, ratio, sampleMode) {
    // Save center of input image
    const xCenter = Math.round(image.width / 2);
    const yCenter = Math.round(image.height / 2);

    // Make a copy of the image to avoid overwriting original pixel values needed later
    const newWidth = Math.ceil(image.width * ratio);
    const newHeight = Math.ceil(image.height * ratio);
    let newImg = image.createImg(newWidth, newHeight);
    const newXCenter = Math.round(newImg.width / 2);
    const newYCenter = Math.round(newImg.height / 2);

    for (let i = 0; i < newImg.width; i++) {
        for (let j = 0; j < newImg.height; j++) {
            // Find destination pixel in source image
            // Find distance from new center, scale by ratio, add original center back
            const u = (i - newXCenter) / ratio + xCenter;
            const v = (j - newYCenter) / ratio + yCenter;

            // Look up the pixel you want in the input image (with resampling technique) and get its value
            const pixel = this.samplePixel(image, u, v, sampleMode);

            newImg.setPixel(i, j, pixel);
        }
    }

    return newImg;
};

// Rotate the image by some angle and using a requested method of sampling/resampling
Filters.rotateFilter = function(image, radians, sampleMode) {
    // Note: set pixels outside the image to RGBA(0,0,0,0)

    // Save center of input image
    const xCenter = Math.round(image.width / 2);
    const yCenter = Math.round(image.height / 2);

    // Make a copy of the image to avoid overwriting original pixel values needed later
    // Size the canvas with the diagonal, the greatest possible dimension needed to fit the rotated image
    const diagonal = Math.ceil( Math.sqrt(image.width*image.width + image.height*image.height) ); 
    let newImg = image.createImg(diagonal, diagonal);
    const newXCenter = Math.round(newImg.width / 2);
    const newYCenter = Math.round(newImg.height / 2);

    for (let i = 0; i < newImg.width; i++) {
        for (let j = 0; j < newImg.height; j++) {
            // Find destination pixel in source image
            // Find distance from new center, rotate by radians, add original center back
            const u = ( (i - newXCenter)*Math.cos(-1*radians) - (j - newYCenter)*Math.sin(-1*radians) ) + xCenter;
            const v = ( (i - newXCenter)*Math.sin(-1*radians) + (j - newYCenter)*Math.cos(-1*radians) ) + yCenter;

            let pixel = new Pixel(0,0,0,0,"rgb");
            if (u < 0 || u > image.width || v < 0 || v > image.height) {
                newImg.setPixel(i, j, pixel);
            } else {
                // Look up the pixel you want in the input image (with resampling technique) and get its value
                pixel = this.samplePixel(image, u, v, sampleMode);
                newImg.setPixel(i, j, pixel);
            }
        }
    }

    return newImg;
};

// Swirl the filter about its center. The rotation of the swirl should be in linear increase
// along the radial axis up to radians
Filters.swirlFilter = function(image, radians, sampleMode) {
    const SCALE_FACTOR = 125;

    // Save center of input image
    const xCenter = Math.round(image.width / 2);
    const yCenter = Math.round(image.height / 2);

    // Make a copy of the image to avoid overwriting original pixel values needed later
    let newImg = image.createImg(image.width, image.height);

    for (let i = 0; i < newImg.width; i++) {
        for (let j = 0; j < newImg.height; j++) {
            // Find this pixel's distance from center
            // The farther it is, the more theta increases
            const dist = Math.sqrt((i-xCenter)*(i-xCenter) + (j-yCenter)*(j-yCenter));
            const theta = -1 * (dist / SCALE_FACTOR) * radians;

            // Find destination pixel in source image
            // Find distance from new center, rotate by radians, add original center back
            const u = ( (i - xCenter)*Math.cos(theta) - (j - yCenter)*Math.sin(theta) ) + xCenter;
            const v = ( (i - xCenter)*Math.sin(theta) + (j - yCenter)*Math.cos(theta) ) + yCenter;

            const pixel = this.samplePixel(image, u, v, sampleMode);
            newImg.setPixel(i, j, pixel);
        }
    }

    return newImg;
};

///////////////////////////////////////////
// Composite Operations
///////////////////////////////////////////

// Set alpha from luminance (sets the alpha value of the image from the luminance of the opacity map)
Filters.getAlphaFilter = function(backgroundImg, foregroundImg) {
    for (let i = 0; i < backgroundImg.height; i++) {
        for (let j = 0; j < backgroundImg.width; j++) {
            const pixelBg = backgroundImg.getPixel(j, i);
            const pixelFg = foregroundImg.getPixel(j, i);
            const luminance =
            0.2126 * pixelFg.data[0] + 0.7152 * pixelFg.data[1] + 0.0722 * pixelFg.data[2];
            pixelBg.a = luminance;
            backgroundImg.setPixel(j, i, pixelBg);
        }
    }

    return backgroundImg;
};

// Composites the foreground image over the background image, using the alpha
// channel of the foreground image to blend two images.
Filters.compositeFilter = function(backgroundImg, foregroundImg) {
    // Assume the input images are of the same sizes.

    for (let x = 0; x < backgroundImg.width; x++) {
        for (let y = 0; y < backgroundImg.height; y++) {
            const pixelBg = backgroundImg.getPixel(x, y);
            const pixelFg = foregroundImg.getPixel(x, y);

            // Use alpha channel of foreground image to blend images together
            pixelBg.data[0] = pixelFg.a * pixelFg.data[0] + (1 - pixelFg.a) * pixelBg.data[0];
            pixelBg.data[1] = pixelFg.a * pixelFg.data[1] + (1 - pixelFg.a) * pixelBg.data[1];
            pixelBg.data[2] = pixelFg.a * pixelFg.data[2] + (1 - pixelFg.a) * pixelBg.data[2];

            backgroundImg.setPixel(x, y, pixelBg);
        }
    }

    return backgroundImg;
};

// Helper function to warp a single line in an image by mapping X (in destination image) to X' (in source image)
// i, j are the coordinates of X: the pixel that we want to map
Filters.warpLine = function(i, j, lineSrc, lineDst) {
    // Save P, Q as the endpoints of lineDst
    const pI = lineDst.x0;
    const pJ = lineDst.y0;
    const qI = lineDst.x1;
    const qJ = lineDst.y1;

    // Save P', Q' as the endpoints of lineSrc
    const pPrimeI = lineSrc.x0;
    const pPrimeJ = lineSrc.y0;
    const qPrimeI = lineSrc.x1;
    const qPrimeJ = lineSrc.y1;

    // --- 1. Calculate u, v based on lineDst endpoints ---
 
    // Vector PX (P to pixel)
    const pToPixelXDir = i - pI;
    const pToPixelYDir = j - pJ;
    // Vector PQ
    const pToQXDir = qI - pI;
    const pToQYDir = qJ - pJ;
    const pToQMag = Math.sqrt(pToQXDir*pToQXDir + pToQYDir*pToQYDir);
    
    // u = scalar projection of PX onto PQ
    const u = ( (pToPixelXDir * pToQXDir) + (pToPixelYDir * pToQYDir) ) / ( pToQXDir*pToQXDir + pToQYDir*pToQYDir );
    
    // Perpendicular vector of PQ
    const pToQPerpXDir = pToQYDir;
    const pToQPerpYDir = -1 * pToQXDir;
    
    // v = pixel's signed distance to PQ
    const v = ( (pToPixelXDir * pToQPerpXDir) + (pToPixelYDir * pToQPerpYDir) ) / pToQMag;

    // --- 2. Calculate coordinates of pixel X' based on u, v and lineSrc endpoints ---

    // Vector P'Q'
    const pPrimeToQPrimeXDir = qPrimeI - pPrimeI;
    const pPrimeToQPrimeYDir = qPrimeJ - pPrimeJ;

    // Perpendicular vector of P'Q'
    const pPrimeToQPrimePerpXDir = pPrimeToQPrimeYDir;
    const pPrimeToQPrimePerpYDir = -1 * pPrimeToQPrimeXDir;

    // Unit vector of P'Q'
    const pPrimeToQPrimeMag = Math.sqrt(pPrimeToQPrimeXDir*pPrimeToQPrimeXDir + pPrimeToQPrimeYDir*pPrimeToQPrimeYDir);
    const pPrimeToQPrimeUnitXDir = pPrimeToQPrimePerpXDir / pPrimeToQPrimeMag;
    const pPrimeToQPrimeUnitYDir = pPrimeToQPrimePerpYDir / pPrimeToQPrimeMag;

    // Get coordinates of X'
    const iPrime = pPrimeI + (u * pPrimeToQPrimeXDir) + (v * pPrimeToQPrimeUnitXDir);
    const jPrime = pPrimeJ + (u * pPrimeToQPrimeYDir) + (v * pPrimeToQPrimeUnitYDir);

    // 3. Calculate displacement between X and X' for this line
    const displacementI = iPrime - i;
    const displacementJ = jPrime - j;

    // 4. Calculate shortest distance from X to lineDst
    const pToPixelMag = Math.sqrt(pToPixelXDir*pToPixelXDir + pToPixelYDir*pToPixelYDir);
    const qToPixelMag = Math.sqrt((i-qI)*(i-qI) + (j-qJ)*(j-qJ));
    let dist = 0;
    if (u < 0) 
        dist = pToPixelMag;
    else if (u > 1)
        dist = qToPixelMag;
    else 
        dist = Math.abs(v);

    // 5. Calculate the weight that line segment PQ contributes to the warping of X's location
    const p = 0.5;
    const a = 0.01;
    const b = 2;
    const weight = Math.pow( (Math.pow(pPrimeToQPrimeMag, p) / (a + dist)), b); 

    return { displacementI, displacementJ, weight };
}

// Helper function to warp many line pairs in an image 
Filters.warpLines = function(initialImg, linesSrc, linesDst, sampleMode) {
    // Create a new image with the same size as the input image
    let newImg = initialImg.createImg(initialImg.width, initialImg.height);

    // Iterate through all pixels in the destination image
    for (let i = 0; i < newImg.width; i++) {
        for (let j = 0; j < newImg.height; j++) {
            let dSumI = 0;
            let dSumJ = 0;

            let weightSum = 0; 

            // Iterate through all lines
            for (let l = 0; l < linesDst.length; l++) {
                // Algorithm to warp a single line
                let { displacementI, displacementJ, weight } = this.warpLine(i, j, linesSrc[l], linesDst[l]);

                // Save dSum and weightSum values
                dSumI += displacementI * weight;
                dSumJ += displacementJ * weight;

                weightSum += weight; 
            }

            // Finalize coordinates for X' by averaging based on line weights
            const iPrime = i + dSumI / weightSum;
            const jPrime = j + dSumJ / weightSum;

            // Look up the pixel you want in the source image (with resampling technique) and get its value
            const pixel = this.samplePixel(initialImg, iPrime, jPrime, sampleMode);
            newImg.setPixel(i, j, pixel);
        }
    }

    return newImg;
}

// Morph two images according to a set of correspondence lines
Filters.morphFilter = function(initialImg, finalImg, alpha, sampleMode, linesFile) {
    // Generate corresponding line segments to be aligned between the before & after images
    const lines = Parser.parseJson("images/" + linesFile);

    // The provided linesFile represents lines in a flipped x, y coordinate system
    //  (i.e. x for vertical direction, y for horizontal direction).
    // Therefore we first fix the flipped x, y coordinates here.
    for (let i = 0; i < lines.initial.length; i++) {
        [lines.initial[i].x0, lines.initial[i].y0] = [lines.initial[i].y0, lines.initial[i].x0]; // foreground lines
        [lines.initial[i].x1, lines.initial[i].y1] = [lines.initial[i].y1, lines.initial[i].x1];
        [lines.final[i].x0, lines.final[i].y0] = [lines.final[i].y0, lines.final[i].x0]; // background lines
        [lines.final[i].x1, lines.final[i].y1] = [lines.final[i].y1, lines.final[i].x1];
    }

    // ----------- STUDENT CODE BEGIN ------------
    let image = finalImg.createImg(finalImg.width, finalImg.height);

    // Make a deep copy of lines.initial to store the values for the intermediate line coordinates
    let linesInt = JSON.parse(JSON.stringify(lines.initial));

    // Iterate through t (time frames) to get a morphing animation
    // for (let t = 0; t < alpha; t++) { 

    const t = alpha; // get one snapshot 
    
    // Interpolate all line correspondences to create an intermediate image
    for (let i = 0; i < lines.initial.length; i++) {
        // Current line coordinates = linear combination of initial lines and final lines
        // L[i] = (1-t) * L0[i] + t * L1[i];
        linesInt[i].x0 = (1-t) * lines.initial[i].x0 + t * lines.final[i].x0;
        linesInt[i].x1 = (1-t) * lines.initial[i].x1 + t * lines.final[i].x1;

        linesInt[i].y0 = (1-t) * lines.initial[i].y0 + t * lines.final[i].y0;
        linesInt[i].y1 = (1-t) * lines.initial[i].y1 + t * lines.final[i].y1;
    }

    // Warp the initial image to the intermediate correspondence
    const fgToInt = this.warpLines(initialImg, lines.initial, linesInt, sampleMode);

    // Warp the final image to the intermediate correspondence
    const bgToInt = this.warpLines(finalImg, lines.final, linesInt, sampleMode);

    // Alpha blend fgToInt and bgToInt
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            const fgMorphPixel = fgToInt.getPixel(x, y);
            const bgMorphPixel = bgToInt.getPixel(x, y);

            pixel.data[0] = (1-t) * fgMorphPixel.data[0] + t * bgMorphPixel.data[0];
            pixel.data[1] = (1-t) * fgMorphPixel.data[1] + t * bgMorphPixel.data[1];
            pixel.data[2] = (1-t) * fgMorphPixel.data[2] + t * bgMorphPixel.data[2];

            image.setPixel(x, y, pixel);
        }
    }

    // }

    return image;
};

///////////////////////////////////////////
// Miscellaneous Operations
///////////////////////////////////////////

// Use k-means to extract a pallete from an image
Filters.paletteFilter = function(image, colorNum) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 89 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('paletteFilter is not implemented yet');
    return image;
};

// Read the following paper and implement your own "painter":
//      http://mrl.nyu.edu/publications/painterly98/hertzmann-siggraph98.pdf
Filters.paintFilter = function(image, value) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 59 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('paintFilter is not implemented yet');
    return image;
};

/*
* Read this paper for background on eXtended Difference-of-Gaussians:
*      http://www.cs.princeton.edu/courses/archive/spring19/cos426/papers/Winnemoeller12.pdf
* Read this paper for an approach that develops a flow field based on a bilateral filter
*      http://www.cs.princeton.edu/courses/archive/spring19/cos426/papers/Kang09.pdf
*/
Filters.xDoGFilter = function(image, value) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 70 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('xDoGFilter is not implemented yet');
    return image;
};

// You can use this filter to do whatever you want, for example
// trying out some new idea or implementing something for the
// art contest.
// Currently the 'value' argument will be 1 or whatever else you set
// it to in the URL. You could use this value to switch between
// a bunch of different versions of your code if you want to
// code up a bunch of different things for the art contest.
Filters.customFilter = function(image, value) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 0 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('customFilter is not implemented yet');
    return image;
};