function patch = gaborPatch(patchSize, gratingFrequency, gratingRotation, ...
    gratingType, filterSigma, filterAspect, filterRotation, foreColour, ...
    backColour, contrast, style, phase)
%gaborPatch     Generate gabor patch
%   patch = gaborPatch(patchSize, gratingFrequency, gratingRotation, gratingType,
%                      filterSigma, filterAspect, filterRotation, foreColour,
%                      backColour, contrast, style)
%
%   Draws a grating filtered by a 2D Gaussian function, both centred on the image.
%   The spatial frequency, type, and orientation of the grating are adjustable.
%   The standard deviation, aspect ratio, and orientation of the filter are also
%   adjustable. So too are the foreground and background colour and contrast.
%
%   Only the first argument is required. Any number of the others may be provided
%   and individual arguments may be left unspecified by replacing with []. Any
%   unspecified arguments will be replaced by a default value.
%
%   patchSize is the width and height of the image in pixels. This parameter must
%       be provided. The resultant image will always have an odd number of rows
%       and columns with the grating and filter centred on the centre pixel. When
%       an even size is specified, one row and one column will be added.
%   gratingFrequency is the spatial frequency of the grating in pixels per cycle. 
%       default is patchSize / 10.
%   gratingRotation is the orientation of the grating clockwise from vertical in
%       degrees. 
%       default is 0.
%   gratingType is the of grating and can be cos (cosine wave), sin (sine wave), 
%       square (square wave), sawtooth (sawtooth wave). square and sawtooth
%       require the signal processing toolbox.
%       default is 'cos' (patch centred on a peak in grating).
%   filterSigma is the standard deviation of the Gaussian filter in pixels. 
%       default is patchSize / 10.
%   filterAspect dictates the shape of the Gaussian filter. when filterAspect 
%       = 1, the filter is circular. otherwise, the filter is elliptical. 
%       default is 1 (circular).
%   filterRotation adjusts the rotation of the filter clockwise in degrees when
%       the filter is elliptical. 
%       default is 0.
%   foreColour is an rgb triplet for the colour of the grating. 
%       default is white [1 1 1].
%   backColour is an rgb triplet for the background colour. 
%       default is grey [.5 .5 .5].
%   contrast is the contrast of the grating against the background. 
%       default is 1 (maximum).
%   style dictates whether the grating is simply added to the background ('uni').
%       or both adds and substracts from the background ('bi'). default is 'uni'.
%   phase adjusts the phase of the grating which might be useful if you want to 
%       produce animated gratings. The values should be in the range 0 to 1 and
%       is the proportion of a whole cycle by which the grating is shifted. 
%       Obviously, a cos grating at phase 0 is equivalent to a sin grating at 
%       phase 0.25 so there is a certain amount of redundancy here.
%       defalut is 0.

%   D.George 2016
%check that neither too few or too many arguments have been given
minargin = 1; maxargin = 12;
narginchk(minargin, maxargin);
%deal with any unspecified arguments and make sure that some values are valid
if nargin < 12 || isempty(phase) || phase < 0, phase = 0; end
if nargin < 11 || isempty(style) || ~any(ismember({'uni', 'bi'}, style)), style = 'uni'; end
if nargin < 10 || isempty(contrast) || ~(contrast >= 0 && contrast <= 1), contrast = 1; end
if nargin < 9 || isempty(backColour) || size(backColour, 2) < 3, backColour = [.5 .5 .5]; end
if nargin < 8 || isempty(foreColour) || size(foreColour, 2) < 3, foreColour = [1 1 1]; end
if nargin < 7 || isempty(filterRotation), filterRotation = 0; end
if nargin < 6 || isempty(filterAspect) || filterAspect <= 0, filterAspect = 1; end
if nargin < 5 || isempty(filterSigma) || filterSigma <= 0, filterSigma = patchSize / 10; end
if nargin < 4 || isempty(gratingType) || ~((exist(gratingType, 'file') || exist(gratingType, 'builtin')) && ...
        any(ismember({'cos', 'sin', 'square', 'sawtooth'}, gratingType))), gratingType = 'cos'; end
if nargin < 3 || isempty(gratingRotation), gratingRotation = 0; end
if nargin < 2 || isempty(gratingFrequency) || gratingFrequency <= 0, gratingFrequency = patchSize / 10; end

%generate grating
xArray = -floor(patchSize / 2):floor(patchSize / 2); %distance from midpoint along x-axis
yArray = -floor(patchSize / 2):floor(patchSize / 2); %distance from midpoint along y-axis
[x, y] = meshgrid(xArray, yArray); %create meshgrids
rotationInRadian = gratingRotation * pi / 180; %convert orientation of grating into radians so that we can use square and sawtooth functions for which there is no equivalent to sind or cosd
spatialFrequency = 1 / gratingFrequency; %convert grating frequency from pixels per cycle to cycles per pixel
radsPerPixel = spatialFrequency * (2 * pi); %and then convert that into radians per pixel
a = cos(rotationInRadian) * radsPerPixel; %work out how much horizontal and vertical is needed to create desired orientation
b = sin(rotationInRadian) * radsPerPixel;
switch gratingType
    case 'sin'
        grating = sin((a * x) + (b * y) + (phase * 2 * pi));
    case 'cos'
        grating = cos((a * x) + (b * y) + (phase * 2 * pi));
    case 'square' %requires signal processing toolbox - if not installed, will revert to default
        grating = square((a * x) + (b * y) + (phase * 2 * pi));
    case 'sawtooth' %requires signal processing toolbox - if not installed, will revert to default
        grating = sawtooth((a * x) + (b * y) + (phase * 2 * pi));
    otherwise
        return
end

%generate Gaussian filter
if filterAspect == 1 %circular filter, so just quickly multiply a couple of univariate normal distributions
    gaussFilter = exp(-(yArray.^2) / (2 * filterSigma^2))' * exp(-(xArray.^2) / (2 * filterSigma^2));
else %otherwise we must draw a 2D elliptical distribution which is a little more involved
    sigmaX = filterSigma; %stated sigma applies to the x-axis
    sigmaY = filterSigma * filterAspect; %aspect ratio affects just the y-axis sigma
    theta = -filterRotation; %makes the grating and filter rotate in the same direction
    A = 1; %A, a, b, c, are parameters in the 2D Gaussian function
    a = (cosd(theta)^2 / (2 * sigmaX^2)) + (sind(theta)^2 / (2 * sigmaY^2));
    b = -(sind(2 * theta) / (4 * sigmaX^2)) + (sind(2 * theta) / (4 * sigmaY^2));
    c = (sind(theta)^2 / (2 * sigmaX^2)) + (cosd(theta)^2 / (2 * sigmaY^2));
    gaussFilter = A * exp(-((a * x.^2) + (2 * b * x.*y) + (c * y.^2)));
end

%apply filter to grating for each colour channel as specified by foreColour and backColour
patch = zeros(size(grating, 1), size(grating, 2), 3); %empty matrix containing three colour channels
for c = 1:1:3
	if strcmp(style, 'bi') %grating sits on top of background
        patch(:, :, c) = backColour(c) + ((foreColour(c) - backColour(c)) * contrast * (grating .* gaussFilter));
    else %grating is centred on background and extends in both colour directions - clipping may occur when backColour is not [.5 .5 .5]
        patch(:, :, c) = backColour(c) + ((foreColour(c) - backColour(c)) * contrast * ((0.5 + (grating * 0.5)) .* gaussFilter));  
	end
end
end