%create object of transducer class requires a defined kgrid object
%and structure field structure field consists of two tables what can and
%cannot be changed postition is set and cannot be changed but steering
%angle and focus distance can be changed

%define an input signal, any 1d array 
%we can alter input signal as needed 

%will run simulation in single precision and hopefully speed up process 
%DATA_CAST = 'gpuArray-single';% need parallel computing toolbox
DATA_CAST = 'single';
%define Kgrid 
%kgrid.t_array
%specific parameters can be defined by markers

%desired grid size based on phantom limb size in m 
x_size = 0.1;
y_size = 0.1;
z_size = 0.05;
%min speed of sound varies in tissue, lowest is 343 in air 
sound_min = 343;
%max FUS frq is 20 kH, so max is set to 5 higher 
fmax = 100000;
%points per wavelength is generally 2-3 but for het medium is 4-5
ppw = 4;
%grid parameter calculation
dx = sound_min/(ppw*fmax);
Nx = round(x_size/dx);

dy = sound_min/(ppw*fmax);
Ny = round(y_size/dy);

dz = sound_min/(ppw*fmax);
Nz = round(z_size/dz);

%create computational grid 
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

%define acoustic medium 
%five material properties: sound speed, mass density, non linearity, alpha
%coefficient, alpha exponent (must be scalar) 


%CT data must be in houndsfield units

%load image (will need to change to the dicom data we will use)
trialImage = dicomread('case1_008'); %load in DICOM file

info = dicominfo('case1_008');%extract DiCOM information

rescale = info.RescaleSlope; %get slope 
ct_data = zeros(size(trialImage)); %create blank mat to hold hf data

%convert dicom data to houndsfield units
for x = 1 : size(trialImage,1)
    for y = 1: size(trialImage,2)
        ct_data(x,y) = trialImage(x,y)*rescale;
    end 
end

[density, sound_speed] = hounsfield2density(ct_data); %specific point by point density 

medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
%medium.density = density;
%medium.sound_speed = sound_speed;
medium.density = 1500;
medium.sound_speed = 1000;

%create t array using the medium
t_end = 40e-6;                  % [s]
kgrid.makeTime(medium.sound_speed, [], t_end);

%Setting both source and sensor to the transducer 

%define input signal 
source_strength = 1e6;         % [MPa] set to generic for now 
burst_freq = 2e6;        % [Hz] high frequency ultrasound frequency
burst_cycles = 5;         %set high for demonstation purposes
input_signal = toneBurst(1/kgrid.dt, burst_freq, burst_cycles); %creates an input signal
input_signal = (source_strength ./ (medium.sound_speed * medium.density)) .* input_signal;

%create transducer object 

tr.number_elements = 72;    % total number of transducer elements
tr.element_width = 1;       % width of each element [grid points]
tr.element_length = 12;     % length of each element [grid points]
tr.element_spacing = 0;     % spacing (kerf width) between the elements [grid points]
tr.radius = inf;        % radius of curvature of the transducer in m (only inf is supported)
tr.sound_speed = 1500;      % sound speed 
%}
%these inputs will need to be taken from the real time transducer 
tr.focus_distance = 0.006529;           % focus distance [m] from specs
tr.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
tr.steering_angle = 0;                  % steering angle [degrees]

%transducer width 
t_width = tr.number_elements*tr.element_width+(tr.number_elements - 1) * tr.element_spacing;

%transducer position 
tr.position = round([1, Ny/2 - t_width/2, Nz/2 - tr.element_length/2]);%generically sets position to 1/2 

% apodization
tr.transmit_apodization = 'Rectangular';%transducer is source
tr.receive_apodization = 'Rectangular';%transducer is sensor 

% define the transducer elements that are currently active
tr.active_elements = zeros(tr.number_elements, 1);
tr.active_elements(21:52) = 1;

% set input signal as created signal
tr.input_signal = input_signal;

% create the transducer
transducer = kWaveTransducer(kgrid, tr);

%for now, a sensor mask 
% create a binary sensor mask with four detection positions
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask([round(Nx/4), round(Nx/2), round(3*Nx/4)], round(Ny/2), round(Nz/2)) = 1;

% run the simulation
input_args = {'DisplayMask', transducer.active_elements_mask,'DataCast', DATA_CAST,};

%with sensor
%[sensor_data] = kspaceFirstOrder3D(kgrid, medium, transducer, sensor, input_args{:});

%tr as sensor
[sensor_data] = kspaceFirstOrder3D(kgrid, medium, transducer, transducer, input_args{:});