import os
import pickle


from tudatpy import constants
from tudatpy.interface import spice
from tudatpy.astro.time_conversion import DateTime
from tudatpy.numerical_simulation import environment_setup
from tudatpy.trajectory_design import transfer_trajectory, shape_based_thrust
from tudatpy.trajectory_design.porkchop import porkchop, plot_porkchop

# Load SPICE kernels
spice.load_standard_kernels()

# Set frame and bodies
global_frame_orientation = 'ECLIPJ2000'
global_frame_origin = 'Sun'
bodies_to_create = ['Sun', 'Earth', 'Mars']
body_settings = environment_setup.get_default_body_settings(
    bodies_to_create, global_frame_origin, global_frame_orientation)

bodies = environment_setup.create_system_of_bodies(body_settings)

departure_body = 'Earth'
target_body = 'Mars'

earliest_departure_time = DateTime(2041,  1,  1)
latest_departure_time   = DateTime(2041, 12,  31)

earliest_arrival_time   = DateTime(2041, 5,  1)
latest_arrival_time     = DateTime(2042, 4,  31)

# Set time resolution IN DAYS as 0.5% of the smallest window (be it departure, or arrival)
# This will result in fairly good time resolution, at a runtime of approximately 10 seconds
# Tune the time resolution to obtain results to your liking!
time_window_percentage = 0.5
time_resolution = time_resolution = min(
        latest_departure_time.epoch() - earliest_departure_time.epoch(),
        latest_arrival_time.epoch()   - earliest_arrival_time.epoch()
) / constants.JULIAN_DAY * time_window_percentage / 100

# File
data_file = 'porkchop.pkl'

# Whether to recalculate the porkchop plot or use saved data
RECALCULATE_delta_v = input(
    '\n    Recalculate ΔV for porkchop plot? [y/N] '
).strip().lower() == 'y'
print()

if not os.path.isfile(data_file) or RECALCULATE_delta_v:
    # Regenerate plot
    [departure_epochs, arrival_epochs, ΔV] = porkchop(
        bodies,
        departure_body,
        target_body,
        earliest_departure_time,
        latest_departure_time,
        earliest_arrival_time,
        latest_arrival_time,
        time_resolution
    )
    # Save data
    pickle.dump(
        [departure_epochs, arrival_epochs, ΔV],
        open(data_file, 'wb')
    )
else:
    # Read saved data
    [departure_epochs, arrival_epochs, ΔV] = pickle.load(
        open(data_file, 'rb')
    )
    # Plot saved data
    plot_porkchop(
        departure_body   = departure_body,
        target_body      = target_body,
        departure_epochs = departure_epochs,
        arrival_epochs   = arrival_epochs,
        delta_v          = ΔV,
        threshold        = 15
    )
