# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# The LED has a fake white light obtained by a blue monochromatic source and a phosphorescent tail

# <codecell>

from scipy.stats import norm, maxwell
import random as rnd
import matplotlib.cm as cm

# Here we define the simulation parameters
n_leds = 3 # number of different light sources (different distributions)
n_samples = 30 # number of differently colored objects
object_size = 2.0 # object size (adimensional)
object_h_mean = 0.5
object_h_spread = 0.01 # width of the normal curve describing the height distribution of the means for various machines
object_h_subspread = 0.01 # width of the normal curve describing the height distribution around each machine mean value
black_camera_absorbance = 5 # absorbance level of the black camera

# define wavelength domain from 0 to 1000nm with 1000 bins
x=linspace(0,1000,1000)

##############################
#### LED
led_samples = []
for _ in range(n_leds):
    position = 1.*(random.random()-0.5)
    width = (random.random()-0.5)*2.
    # define the emission spectrum of the white LED
    led_blue = norm.pdf(x, loc=460+position, scale=10.+width) # blue emission
    led_blue /= led_blue.max()
    led_red = norm.pdf(x, loc=560+position, scale=50.+width) # phosphorescence
    led_red /= led_red.max()
    # the amplitude due to phosphorescence is ~40% of the main blue source
    led_total = led_blue + 0.4*led_red
    led_samples.append(led_total)

subplot(411)
grid(True)
axis([350,750,0,1.5])
xlabel('Wavelength $\lambda$ (nm)')
ylabel('White LED Relative Intensity')
title(r'White LED Spectrum')
for led in led_samples:
    plot(led)
###############################
    
#### Object reflectivity (set of differently colored objects)
obj_samples = []
for _ in range(n_samples):
    position = 500+200*(random.random()-0.5)
    width = 90*random.random()
    obj = norm.pdf(x, loc=position, scale=width)
    obj /= obj.max()
    obj_samples.append(obj)

subplot(412)
grid(True)
axis([350,750,0,1.5])
xlabel('Wavelength $\lambda$ (nm)')
ylabel('Reflectance coefficient')
title(r'Sample reflectivity behaviour')
for obj in obj_samples:
    plot(obj)
###############################
    
#### Black camera reflectivity

black_camera = maxwell.pdf(x, loc=350, scale=200,size=1000)
black_camera /= black_camera.max() * black_camera_absorbance

subplot(413)
grid(True)
axis([350,750,0,1.5])
xlabel('Wavelength $\lambda$ (nm)')
ylabel('Reflectance coefficient')
title(r'Black camera reflectivity behaviour')
plot(black_camera)
###############################

#### Reflected light from object and camera

subplot(414)
grid(True)
axis([350,750,0,1.5])
xlabel('Wavelength $\lambda$ (nm)')
ylabel('Intensity')
title(r'Reflected intensity arriving at the sensor')

for led in led_samples:
    for obj in obj_samples:
        plot(led * obj)
    plot(led * black_camera)

subplots_adjust(hspace=0.3)
###############################

#### Sensor characteristics (TAOS 3210)

figure(figsize=(8, 16))
subplot(211)
grid(True)
axis([350,1000,0,1.5])
xlabel('Wavelength $\lambda$ (nm)')
ylabel('Sensitivity')
title(r'Sensor sensitivity for RGB')

sens_red = norm.pdf(x, loc=730, scale=100) #red
sens_red_uv = norm.pdf(x, loc=400, scale=10) #uv
sens_red_uv /= sens_red_uv.max() / 0.2
sens_red /= sens_red.max() / 1.0
sens_red += sens_red_uv

sens_green = norm.pdf(x, loc=524, scale=50) #green
sens_green_ir = norm.pdf(x, loc=850, scale=10) #infrared
sens_green_ir /= sens_green_ir.max() / 0.75
sens_green /= sens_green.max() / 0.55
sens_green += sens_green_ir

sens_blue = norm.pdf(x, loc=470, scale=75) #blue
sens_blue_ir = norm.pdf(x, loc=840, scale=13) #ir
sens_blue_ir /= sens_blue_ir.max() / 0.85
sens_blue /= sens_blue.max() / 0.47
sens_blue += sens_blue_ir

plot(sens_blue)
plot(sens_green)
plot(sens_red)

subplot(212)
grid(True)
xlabel('Incoming angle')
ylabel('Intensity')
title(r'Angular response of the sensor')
angles=linspace(-3.14/2.0,3.14/2.0,1000)
axis([-3.14/2,3.14/2,0,1.1])
plot(angles,cos(angles)) # plot the lambertian cosine law for the sensor

subplots_adjust(hspace=0.15)
###############################

#### Computation of the total flux for every channel due to object and camera

grid(True)
xlabel('Photodiode Channel')
ylabel('Output')
title(r'Numerical output of the sensor')

fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

colors = cm.rainbow(linspace(0, 1, n_samples))
symbols = ['o','+','x','o','o','o']
for machine in range(3): # loop over some number of distinct machines
    loc_machine = random.normal(loc=object_h_mean, scale=object_h_spread, size=None)
    for led in led_samples: # loop over various illuminations
        for _ in range(15): # loop randomizes for object height within the chamber
            object_h = random.normal(loc=loc_machine, scale=object_h_subspread, size=None)
            min_ang_obj = -math.atan(object_size/(2.*object_h))
            max_ang_obj = math.atan(object_size/(2.*object_h))
            # define the integration angles for object and camera (the camera is modelled as an infinite sphere)
            angles_obj = [val if val>min_ang_obj and val<max_ang_obj else pi/2. for val in angles]
            angles_camera = [val if val>max_ang_obj or val<min_ang_obj else pi/2. for val in angles]
            
            for i, light in enumerate(obj_samples*led):
                obj_red = (sens_red*light).sum()
                obj_red *= cos(angles_obj).sum()
                camera_red = (sens_red*black_camera*led).sum()
                camera_red *= cos(angles_camera).sum()
                tot_red = obj_red+camera_red
                
                obj_green = (sens_green*light).sum()
                obj_green *= cos(angles_obj).sum()
                camera_green = (sens_green*black_camera*led).sum()
                camera_green *= cos(angles_camera).sum()
                tot_green = obj_green+camera_green
                
                obj_blue = (sens_blue*light).sum()
                obj_blue *= cos(angles_obj).sum()
                camera_blue = (sens_blue*black_camera*led).sum()
                camera_blue *= cos(angles_camera).sum()
                tot_blue = obj_blue+camera_blue
            
                scatter(tot_red/(tot_red+tot_green+tot_blue), tot_green/(tot_red+tot_green+tot_blue),marker=symbols[machine], s=25, c=colors[i])
                #ax.scatter (tot_red, tot_green, tot_blue, marker=symbols[machine], s=5, c=i)

# <codecell>


