import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
from astropy import units as u
from astropy.time import Time, TimeDelta
from poliastro.bodies import Mars
from poliastro.twobody import Orbit
from poliastro.util import time_range
from astropy.coordinates import get_body_barycentric

# def plot_shadow_cones(ax, sun_direction, mars_radius, R_sun, sun_distance, shadow_length=50000):
#     alpha_umbra = np.arctan((mars_radius - R_sun) / sun_distance)
#     alpha_penumbra = np.arctan((mars_radius + R_sun) / sun_distance)
#     L = shadow_length

#     def make_cone(angle, color, label):
#         h = np.linspace(0, L, 30)
#         theta = np.linspace(0, 2 * np.pi, 30)
#         H, Theta = np.meshgrid(h, theta)
#         R_cone = H * np.tan(angle) * 100
#         X = R_cone * np.cos(Theta)
#         Y = R_cone * np.sin(Theta)
#         Z = H

#         cone_pts = np.stack([X, Y, Z], axis=-1).reshape(-1, 3)
#         z_axis = np.array([0, 0, 1])
#         rotvec = np.cross(z_axis, sun_direction)
#         if np.linalg.norm(rotvec) != 0:
#             angle_rot = np.arccos(np.dot(z_axis, sun_direction))
#             rotation = R.from_rotvec(rotvec / np.linalg.norm(rotvec) * angle_rot)
#             cone_pts = rotation.apply(cone_pts)

#         cone_pts = cone_pts.reshape(X.shape + (3,))
#         Xr, Yr, Zr = cone_pts[..., 0], cone_pts[..., 1], cone_pts[..., 2]
#         ax.plot_surface(Xr, Yr, Zr, alpha=0.15, color=color, label=label)

#     make_cone(alpha_umbra, "gray", "Umbra")
#     make_cone(alpha_penumbra, "orange", "Penumbra")
#     ax.scatter(0, 0, 0, color='brown', s=100, label='Mars')
#     ax.set_xlabel("X (km)")
#     ax.set_ylabel("Y (km)")
#     ax.set_zlabel("Z (km)")
#     ax.legend()
#     ax.set_box_aspect([1, 1, 1])


def compute_eclipse_duration():
    start_time = Time("2025-01-01 00:00:00", scale="utc")
    mars_radius = Mars.R.to(u.km).value
    R_sun = 695700  # km

    alt = 200 * u.km
    a = mars_radius * u.km + alt
    orb = Orbit.circular(Mars, alt, inc=92 * u.deg, epoch=start_time)

    period = orb.period
    times = time_range(start_time, spacing=TimeDelta((period / 500).to(u.s)), periods=500)

    sun_positions = np.array([
        (get_body_barycentric("sun", t) - get_body_barycentric("mars", t)).xyz.to(u.km).value
        for t in times
    ])
    sc_positions = np.array([
        orb.propagate(t - orb.epoch).rv()[0].to(u.km).value
        for t in times
    ])

    penumbra_flags = []
    for sc_pos, sun_pos in zip(sc_positions, sun_positions):
        v_sun = sun_pos
        v_sat = sc_pos
        angle = np.arccos(np.dot(v_sun, v_sat) / (np.linalg.norm(v_sun) * np.linalg.norm(v_sat)))
        r_sc = np.linalg.norm(v_sat)
        r_sun = np.linalg.norm(v_sun)
        alpha = np.arcsin(mars_radius / r_sc)
        beta = np.arcsin(R_sun / r_sun)
        penumbra = (angle > alpha - beta) and (angle < alpha + beta)
        penumbra_flags.append(penumbra)

    penumbra_flags = np.array(penumbra_flags)
    dt = (period / len(times)).to(u.s)
    penumbra_duration = np.sum(penumbra_flags) * dt
    print(f"Penumbra eclipse duration: {penumbra_duration.to(u.min):.2f}")

    # 2D Time Plot
    plt.figure()
    plt.plot(times.datetime, penumbra_flags.astype(int), drawstyle="steps-post")
    plt.ylabel("In Penumbra")
    plt.xlabel("Time (UTC)")
    plt.title("Penumbra Eclipse Over One Orbit")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # 3D Visualization
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # sun_direction = -sun_positions[0] / np.linalg.norm(sun_positions[0])
    # sun_distance = np.linalg.norm(sun_positions[0])
    # shadow_direction = -sun_direction 
    # plot_shadow_cones(ax, shadow_direction, mars_radius, R_sun, sun_distance)
    # ax.plot(sc_positions[:, 0], sc_positions[:, 1], sc_positions[:, 2], label='Orbit', color='blue')
    # ax.scatter(sc_positions[penumbra_flags, 0], sc_positions[penumbra_flags, 1], sc_positions[penumbra_flags, 2], color='red', s=10, label='Penumbra')
    # ax.legend()
    # plt.tight_layout()
    # plt.show()

if __name__ == "__main__":
    compute_eclipse_duration()
