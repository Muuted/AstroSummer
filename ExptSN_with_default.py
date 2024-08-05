import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from astropy import units as u
import pytz

def calculate_airmass(observatory_coords, star_coords):
    location = EarthLocation(lat=observatory_coords[0], lon=observatory_coords[1], height=observatory_coords[2])
    timezone = pytz.timezone(observatory_coords[3])
    time = Time.now().to_datetime(timezone)  # Current time in local time zone
    time = Time(time)
    altaz_frame = AltAz(obstime=time, location=location)

    star_altaz = star_coords.transform_to(altaz_frame)
    alt = star_altaz.alt.deg
    
    # Compute the airmass scaled to be 1 at zenith
    if alt > 90:
        alt = 90  # Limit to 90 degrees as zenith

    airmass = 1 / np.cos(np.radians(90 - alt))  # Airmass calculation, secant of the zenith angle
    return airmass

def signal_to_noise_ratio(flux, exposure_time, aperture_diameter, efficiency, airmass):
    extinction_coefficient = 0.1  # Magnitudes per airmass
    
    # Telescope area
    area = np.pi * (aperture_diameter / 2) ** 2
    
    # Total photons collected (corrected for efficiency and atmospheric extinction)
    #total_photons = flux * exposure_time * area * efficiency * 10**(-0.4 * extinction_coefficient * (airmass - 1))
    total_photons = flux * exposure_time* area

    # Assuming shot noise is dominant, SNR = sqrt(signal)
    snr = np.sqrt(total_photons)
    return snr

# Default values
default_observatory_name = "Mt. Kent"
default_obs_lat = -27.0
default_obs_lon = 151
default_obs_alt = 682
default_obs_timezone = 'Australia/Brisbane'
default_aperture_diameter = 0.6

default_star_name = "WASP-74"
default_star_ra = "20:18:09.32"
default_star_dec = "-01:04:33.61"
default_star_magnitude = 9.2

# Input prompts with defaults
observatory_name = input(f"Enter the observatory name (default: {default_observatory_name}): ") or default_observatory_name
obs_lat = float(input(f"Enter the observatory latitude (in degrees, default: {default_obs_lat}): ") or default_obs_lat)
obs_lon = float(input(f"Enter the observatory longitude (in degrees, default: {default_obs_lon}): ") or default_obs_lon)
obs_alt = float(input(f"Enter the observatory altitude (in meters, default: {default_obs_alt}): ") or default_obs_alt)
obs_timezone = input(f"Enter the observatory timezone (e.g., 'Europe/Copenhagen', default: {default_obs_timezone}): ") or default_obs_timezone
aperture_diameter = float(input(f"Enter the telescope's aperture diameter in meters (default: {default_aperture_diameter}): ") or default_aperture_diameter)

star_name = input(f"Enter the star name (default: {default_star_name}): ") or default_star_name
star_ra = input(f"Enter the star's right ascension (RA) in hh:mm:ss format (default: {default_star_ra}): ") or default_star_ra
star_dec = input(f"Enter the star's declination (Dec) in dd:mm:ss format (default: {default_star_dec}): ") or default_star_dec
star_magnitude = float(input(f"Enter the star's magnitude (default: {default_star_magnitude}): ") or default_star_magnitude)

# Observing conditions
flux = 1e2  # photons/s/m^2
efficiency = 0.8  # rough estimate

# Calculate airmass
observatory_coords = (obs_lat, obs_lon, obs_alt, obs_timezone)
star_coords = SkyCoord(f"{star_ra} {star_dec}", unit=(u.hourangle, u.deg))
airmass = calculate_airmass(observatory_coords, star_coords)
print(f'Airmass: {airmass:.2f}')

# Exposure times range (in seconds)
exposure_times = np.linspace(1, 3600, 100)  # From 1 second to 1 hour

# Calculate SNR for each exposure time
snr_values = [signal_to_noise_ratio(flux, t, aperture_diameter, efficiency, airmass) for t in exposure_times]

# Find exposure times corresponding to target SNR values
target_snrs = [4, 5, 6, 10]
exposure_times_for_snrs = []

for target_snr in target_snrs:
    # Find the exposure time where SNR is closest to the target SNR
    idx = np.argmin(np.abs(np.array(snr_values) - target_snr))
    exposure_times_for_snrs.append(exposure_times[idx])

# Print exposure times for each target SNR
for snr, exp_time in zip(target_snrs, exposure_times_for_snrs):
    print(f'Exposure time for SNR {snr}: {exp_time:.2f} seconds')

# Plot the results
fig, ax1 = plt.subplots(figsize=(10, 6))

color = 'tab:blue'
ax1.set_xlabel('Exposure Time (s)')
ax1.set_ylabel('Signal-to-Noise Ratio (SNR)', color=color)
ax1.plot(exposure_times, snr_values, color=color, label='SNR')
ax1.tick_params(axis='y', labelcolor=color)

# Add text to the plot for observatory and star information
plt.text(0.05, 0.95, f'Observatory: {observatory_name} {obs_lat};{obs_lon}@{obs_alt}m', transform=plt.gcf().transFigure)
plt.text(0.05, 0.90, f'Star: {star_name}: {star_ra}:{star_dec} Mag: {star_magnitude}', transform=plt.gcf().transFigure)

plt.title('SNR vs Exposure Time')
fig.tight_layout()  
plt.grid(True)
plt.show()
