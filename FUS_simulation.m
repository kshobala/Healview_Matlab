%create object of transducer class requires a defined kgrid object
%and structure field structure field consists of two tables what can and
%cannot be changed postition is set and cannot be changed but steering
%angle and focus distance can be changed

%define an input signal, any 1d array 
%we can alter input signal as needed 

%will run simulation in single precision and hopefully speed up process 
%DATA_CAST = 'gpuArray-single'; % need parallel computing toolbox
DATA_CAST = 'single';

%%%%%%%%%%%%%%%%%%%%%%%%
%define Kgrid 
%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%
%define acoustic medium 
%%%%%%%%%%%%%%%%%%%%%%%%%%

%five material properties: sound speed, mass density, non linearity, alpha
%coefficient, alpha exponent (must be scalar) 


%CT data must be in houndsfield units

%load image (will need to change to the dicom data we will use)
trialImage = dicomread('case1_008'); %load in DICOM file

info = dicominfo('case1_008');%extract DiCOM information

rescale = info.RescaleSlope; %get slope 
ct_data = zeros(Nx,Ny,Nz); %create blank mat to hold hf data

%convert dicom data to houndsfield units
for x = 1 : size(ct_data,1)
    for y = 1: size(ct_data,2)
        ct_data(x,y) = trialImage(x,y)*rescale;
    end 
end

[density, sound_speed] = hounsfield2density(ct_data); %specific point by point density

medium.alpha_coeff = 0.75;% [dB/(MHz^y cm)]
medium.alpha_power = 1.5;

%
%using DICOM
%medium.density = density;
%medium.sound_speed = sound_speed;
%}

%{
%Scalar
medium.density = 1500;
medium.sound_speed = 1000;
%}

%
%calculated
medium.sound_speed = 1500 * ones(Nx, Ny,Nz);   % [m/s]
medium.sound_speed(1:Nx/2, :) = 1800;       % [m/s]
medium.density = 1000 * ones(Nx, Ny,Nz);       % [kg/m^3]
medium.density(:, Ny/4:Ny) = 1200;          % [kg/m^3]
%}
%create t array using the medium
t_end = 40e-6;                  % [s]
kgrid.makeTime(medium.sound_speed, [], t_end);

%Setting both source and sensor to the transducer 

%define input signal 
source_strength = 1e6;         % [MPa] set to generic for now 
burst_freq = 0.5e6;        % [Hz] high frequency ultrasound frequency
burst_cycles = 5;         %set high for demonstation purposes
input_signal = toneBurst(1/kgrid.dt, burst_freq, burst_cycles); %creates an input signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create transducer object 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
tr.number_elements = 52;    % total number of transducer elements
tr.element_width = 1;       % width of each element [grid points]
tr.element_length = 1;     % length of each element [grid points]
tr.element_spacing = 0;     % spacing (kerf width) between the elements [grid points]
tr.radius = inf;        % radius of curvature of the transducer in m (only inf is supported)
tr.sound_speed = 1500;      % sound speed 



%these inputs will need to be taken from the real time transducer 
tr.focus_distance = 0.006529;           % focus distance [m] from specs
tr.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
tr.steering_angle = 0;                  % steering angle [degrees]

%transducer width 
t_width = tr.number_elements*tr.element_width+(tr.number_elements - 1) * tr.element_spacing;

%transducer position 
tr.position = round([1, 1, 1]);%generically sets position to 1/2 

%apodization
tr.transmit_apodization = 'Rectangular';%transducer is source
tr.receive_apodization = 'Rectangular';%transducer is sensor 

%define the transducer elements that are currently active
tr.active_elements = zeros(tr.number_elements, 1);
tr.active_elements(21:52) = 1;

%set input signal as created signal
tr.input_signal = input_signal;

%create the transducer
transducer = kWaveTransducer(kgrid, tr);

%sensor mask 
%sensor.mask = zeros(Nx, Ny, Nz);
%sensor.mask([round(Nx/4), round(Nx/2), round(3*Nx/4)], round(Ny/2), round(Nz/2)) = 1;

% run the simulation
input_args = {'DisplayMask', transducer.active_elements_mask,'DataCast', DATA_CAST,'RecordMovie', true, 'MovieName', 'tr_source_het_med', ...
        'MovieProfile', 'MPEG-4', 'MovieArgs', {'FrameRate', 10}};

%with sensor
%[sensor_data] = kspaceFirstOrder3D(kgrid, medium, transducer, sensor, input_args{:});

%tr as sensor
[sensor_data] = kspaceFirstOrder3D(kgrid, medium, transducer, transducer, input_args{:});
transducer.plot;
%}
%end transducer sim 
%
%trial with a curved transducer 
% define a curved transducer element
bowl_pos  = [round(Nx/2),round(Ny/2),Nz];
radius    = 30;
diameter  = 11;
focus_pos = [round(Nx/2),round(Ny/2), 32];

% create bowl
source.p_mask = makeBowl([Nx,Ny,Nz],bowl_pos, radius, diameter, focus_pos, 'Plot', true);

% set input sigal 
source.p = input_signal;

% define a time varying sinusoidal source
%source_freq = 20e6;  % [Hz]
%source_mag = 1;     % [Pa]
%source.p = source_mag * sin(2 * pi * source_freq * kgrid.t_array);

% create a display mask to display the transducer
display_mask = source.p_mask;

sensor.mask = [1, 1,1, Nx, Ny,Nz].';

% assign the input options
%input_args = {'DisplayMask', display_mask,'PMLInside', false, 'PlotPML', false,'RecordMovie', true, 'MovieName', 'circ_source_het_med', ...'MovieProfile', 'MPEG-4', 'MovieArgs', {'FrameRate', 10}};
input_args = {'DisplayMask', display_mask,'PMLInside', false, 'PlotPML',false}; % trial without movie


[sensor_data] = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
%}
