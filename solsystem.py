import csv, os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button
from astropy.time import Time
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.image as mpimg

#Initialize the 3D plot
fig = plt.figure(figsize=(100,100))
fig.subplots_adjust(top=1.5, bottom=-.5)
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=25, azim=25)
anim_running=True

def reading_data():
    filepath = r"C:\Users\Regis\Desktop\Python\soldata.csv"
    object_data = []
    with open(filepath, mode='r', encoding='utf-8-sig') as file:
        reader = csv.DictReader(file, skipinitialspace=True, delimiter=";")
        for row in reader:
            row["semi_major_axis"] = safe_float_conversion(row["semi_major_axis"])
            row["semi_major_axis_dot"] = safe_float_conversion(row["semi_major_axis_dot"])
            row["eccentricity"] = (safe_float_conversion(row["eccentricity"]))
            row["inclination"] = np.radians(safe_float_conversion(row["inclination"]))
            row["p"] = np.radians((safe_float_conversion(row["longitude_perihelion"])))
            row["p_delta"] = np.radians((safe_float_conversion(row["longitude_perihelion_change"])))
            row["W"] = np.radians((safe_float_conversion(row["ascending_node"])))
            row["W_delta"] = np.radians((safe_float_conversion(row["node_change"])))
            row["L0"] = ((safe_float_conversion(row["mean_longitude"])))
            row["L_change"] = ((safe_float_conversion(row["mean_longitude_change"])))
            row["M0"] = np.radians((safe_float_conversion(row["M0"])))
            row["b"] = ((safe_float_conversion(row["b"])))
            row["c"] = ((safe_float_conversion(row["c"])))
            row["s"] = ((safe_float_conversion(row["s"])))
            row["f"] = ((safe_float_conversion(row["f"])))
            object_data.append(row)
    return {object['name']: CelestialObjects(object['name'], object) for object in object_data}

def safe_float_conversion(value, default=0.0):
    try:
        return float(value)
    except (ValueError, TypeError):
        return default

class interface:
    t_delta_text = None
    date_text = None
    
    def create_date_text(fig):
        interface.date_text = fig.text(0.015, 0.95, "", fontsize=12, color='white')
    def create_t_delta_text(fig):
        interface.t_delta_text = fig.text(0.015, 0.93, f'Earth days per second: {delta_increment}', fontsize=12, color='white')
    def max_zoom():
        try:
            mng = plt.get_current_fig_manager()
            mng.window.state('zoomed')  # Works on some Windows systems
        except AttributeError:
            # Fallback if the above does not work
            try:
                mng.window.showMaximized()  # An alternative method
            except AttributeError:
                # If maximizing isn't supported, we just pass
                pass
        
        zoom_factor = 0.45  # 200% means reducing the limits by 
        
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        zlim = ax.get_zlim()

        # Calculate the center of the plot
        x_center = (xlim[1] + xlim[0]) / 2
        y_center = (ylim[1] + ylim[0]) / 2
        z_center = (zlim[1] + zlim[0]) / 2

        # Define the zoom factor (2x zoom)
        zoom_factor = 0.45  # 200% means reducing the limits by half

        # Set new limits based on zoom factor
        ax.set_xlim([x_center - (xlim[1] - x_center) * zoom_factor,
                    x_center + (xlim[1] - x_center) * zoom_factor])
        ax.set_ylim([y_center - (ylim[1] - y_center) * zoom_factor,
                    y_center + (ylim[1] - y_center) * zoom_factor])
        ax.set_zlim([z_center - (zlim[1] - z_center) * zoom_factor,
                    z_center + (zlim[1] - z_center) * zoom_factor])         
    def hide_axes():
        ax.set_axis_off()
    def background():
        num_stars = 2000
        x_stars = np.random.uniform(-70, 70, num_stars)
        y_stars = np.random.uniform(-70, 70, num_stars)
        z_stars = np.random.uniform(-70, 70, num_stars)
        
        ax.set_facecolor('black')
        fig.patch.set_facecolor('black')
        
        ax.scatter(x_stars, y_stars, z_stars, color='white', s=0.7)
    def start_animation(event):
        global anim_running
        if not anim_running:
            ani.event_source.start()
            anim_running = True       
    def stop_animation(event):
        global anim_running
        if anim_running:
            ani.event_source.stop()
            anim_running = False
    def increase_speed(event):
        global delta_increment
        delta_increment += 1
        interface.update_t_delta_text()
    def decrease_speed(event):
        global delta_increment
        delta_increment -= 1
        interface.update_t_delta_text()
    def update_t_delta_text():
        interface.t_delta_text.set_text(f'Earth days per second: {delta_increment}')
    def update_date_annotation(t_julian, t_delta):
        t_julian_current = t_julian + t_delta
        time_obj = Time(t_julian_current, format='jd')
        current_date = time_obj.to_value('iso', 'date')
        interface.date_text.set_text(f"Date: {current_date}")

class belts:
    @staticmethod
    def gen_kuiper_belt():
        KUIPER_BELT_INNER = 30  # Inner edge of the Kuiper Belt in AU
        KUIPER_BELT_OUTER = 50  # Outer edge of the Kuiper Belt in AU
        KUIPER_BELT_CENTER = (KUIPER_BELT_INNER + KUIPER_BELT_OUTER) / 2  # Center of the belt
        KUIPER_BELT_STD_DEV = 3  # Standard deviation for the belt distribution
        KUIPER_BELT_POINTS = 10000  # Number of points to represent the Kuiper Belt
        KUIPER_BELT_HEIGHT = 3  # Height of the Kuiper Belt in AU

        # Generate random points within the Kuiper Belt with normal distribution
        kuiper_belt_r = np.random.normal(KUIPER_BELT_CENTER, KUIPER_BELT_STD_DEV, KUIPER_BELT_POINTS)

        # Ensure the values stay within the defined inner and outer edges
        kuiper_belt_r = np.clip(kuiper_belt_r, KUIPER_BELT_INNER, KUIPER_BELT_OUTER)

        kuiper_belt_theta = np.random.uniform(0, 2 * np.pi, KUIPER_BELT_POINTS)
        kuiper_belt_x = kuiper_belt_r * np.cos(kuiper_belt_theta)
        kuiper_belt_y = kuiper_belt_r * np.sin(kuiper_belt_theta)

        # Generate random z-coordinates to give the belt some height
        kuiper_belt_z = np.random.uniform(-KUIPER_BELT_HEIGHT / 2, KUIPER_BELT_HEIGHT / 2, KUIPER_BELT_POINTS)

        # Set the inclination in degrees
        inclination_degrees = 1.86
        inclination_radians = np.radians(inclination_degrees)

        # Rotate the points around the x-axis by the inclination angle
        kuiper_belt_y_rot = kuiper_belt_y * np.cos(inclination_radians) - kuiper_belt_z * np.sin(inclination_radians)
        kuiper_belt_z_rot = kuiper_belt_y * np.sin(inclination_radians) + kuiper_belt_z * np.cos(inclination_radians)

        # Define the color of the Kuiper Belt
        belt_color = (0.5, 0.5, 0.5)

        # Plot the Kuiper Belt with the applied inclination and height
        ax.scatter(kuiper_belt_x, kuiper_belt_y_rot, kuiper_belt_z_rot, color=belt_color, s=1)
        
    def gen_main_belt():
        MAIN_BELT_INNER = 1.70  # Inner edge of the Main Belt in AU
        MAIN_BELT_OUTER = 4.00  # Outer edge of the Main Belt in AU
        MAIN_BELT_CENTER = (MAIN_BELT_INNER + MAIN_BELT_OUTER) / 2  # Center of the belt
        MAIN_BELT_STD_DEV = 0.35  # Standard deviation for the belt distribution
        MAIN_BELT_POINTS = 3000  # Number of points to represent the Main Belt
        MAIN_BELT_HEIGHT = 0.8  # Height of the Main Belt in AU

        # Generate random points within the Main Belt with normal distribution
        main_belt_r = np.random.normal(MAIN_BELT_CENTER, MAIN_BELT_STD_DEV, MAIN_BELT_POINTS)
        main_belt_r = np.clip(main_belt_r, MAIN_BELT_INNER, MAIN_BELT_OUTER)
        
        main_belt_theta = np.random.uniform(0, 2 * np.pi, MAIN_BELT_POINTS)
        main_belt_x = main_belt_r * np.cos(main_belt_theta)
        
        # Generate random z-coordinates with a normal distribution centered at 0
        main_belt_z = np.random.normal(0, MAIN_BELT_HEIGHT / 4, MAIN_BELT_POINTS)  # Adjust standard deviation as needed
        
        # Create a "diamond-shaped" distribution in the y-axis based on z
        main_belt_y = main_belt_r * np.sin(main_belt_theta) * (1 - np.abs(main_belt_z) / MAIN_BELT_HEIGHT)

        # Set the inclination in degrees
        inclination_degrees = 2
        inclination_radians = np.radians(inclination_degrees)

        # Rotate the points around the x-axis by the inclination angle
        main_belt_y_rot = main_belt_y * np.cos(inclination_radians) - main_belt_z * np.sin(inclination_radians)
        main_belt_z_rot = main_belt_y * np.sin(inclination_radians) + main_belt_z * np.cos(inclination_radians)

        # Define the color of the Main Belt
        belt_color = (0.4, 0.4, 0.4)

        # Plot the Main Belt with the applied inclination and height
        ax.scatter(main_belt_x, main_belt_y_rot, main_belt_z_rot, color=belt_color, s=0.5)

class CelestialObjects:
    def __init__ (self, name, data):
        self.name = name
        self.ecc = data["eccentricity"]
        self.semi_major_axis = data["semi_major_axis"]
        self.semi_major_axis_dot = data["semi_major_axis_dot"]
        self.inclination = data["inclination"]
        self.p = data["p"]
        self.p_delta = data["p_delta"]
        self.W = data["W"]
        self.W_delta = data["W_delta"]
        self.L0 = data["L0"]
        self.Ldot = data["L_change"]
        self.M = data["M0"]
        self.b = data["b"]
        self.c = data["c"]
        self.s = data["s"]
        self.f = data["f"]
        self.cen_semi_major_axis = self.semi_major_axis
        self.cen_p = self.p
        self.cen_W = self.W
        self.cen_L0 = self.L0
        self.cen_b = 0
        self.cen_c = 0
        self.cen_s = 0
        
    def time(self, t_delta):
        sim_start_time = Time.now().jd
        t_delta_start_now = sim_start_time - 2451545.000000000 + t_delta
        t_delta_cen = t_delta_start_now / 36525.0
        return sim_start_time, t_delta_cen
        
    def update_stats_delta(self, t_delta):
        _, t_delta_cen = CelestialObjects.time(self, t_delta)
        self.cen_semi_major_axis = self.semi_major_axis + self.semi_major_axis_dot * t_delta_cen
        self.cen_p = self.p + self.p_delta * t_delta_cen
        self.cen_W = self.W + self.W_delta * t_delta_cen
        self.cen_L0 = self.L0 + self.Ldot * t_delta_cen
        self.cen_b = self.b * t_delta_cen**2
        self.cen_c = self.c * np.cos(self.f * t_delta_cen)
        self.cen_s = self.s * np.sin(self.f * t_delta_cen)
    
    def __str__ (self):
        attributes = ', '.join(f"{key}={value}" for key, value in vars(self).items())
        return f"CelestialObject({attributes})"
    
    def solve_kepler_equation (self):
        L = (self.cen_L0
            + self.cen_b
            + self.cen_c
            + self.cen_s)
        
        L = np.mod(L, 360)  # Normalize L to within 0-360 degrees
        
        # Calculate Mean Anomaly (M) and Argument of Perihelion (w)
        M = np.mod(L - np.degrees(self.cen_p), 360)  # Mean anomaly, normalized to 0-360
        w = np.degrees(self.cen_p) - np.degrees(self.cen_W)  # Argument of periapsis
        w = np.mod(w, 360)
        
        # Modulate Mean Anomaly M to fall within [-180°, +180°]
        M = np.mod(M + 180, 360) - 180
        
        # Convert M and w to radians
        self.M_radians = np.radians(M)
        self.w_radians = np.radians(w)
        
        # Solve Kepler's equation to find Eccentric Anomaly (E)
        E = self.M_radians  # Initial guess for E is M
        tol = 1e-6
        while True:
            dE = (E - self.ecc * np.sin(E) - self.M_radians) / (1 - self.ecc * np.cos(E))
            E -= dE
            if abs(dE) < tol:
                break
        self.E = E
    
    def draw_orbit(self, theta):
            x_orb = []
            y_orb = []
            z_orb = []

            for t in theta:
                # Calculate Eccentric Anomaly for each angle t (true anomaly)
                CelestialObjects.solve_kepler_equation(self)
                E = 2 * np.arctan(np.sqrt((1 - self.ecc) / (1 + self.ecc)) * np.tan(t / 2))
                
                # Calculate x, y, z coordinates
                x_orbit, y_orbit, z_orbit = self.xyz_calculation(E)

                # Append the calculated positions to the lists
                x_orb.append(x_orbit)
                y_orb.append(y_orbit)
                z_orb.append(z_orbit)

            # Plot the orbit for this planet
            ax.plot(x_orb, y_orb, z_orb, label=f'{self} Orbit')
            ax.set_aspect('equal')
            
    def xyz_calculation(self, E):
        # Compute the heliocentric coordinates in the orbital plane
        x_prime = self.cen_semi_major_axis * (np.cos(E) - self.ecc)
        y_prime = self.cen_semi_major_axis* np.sqrt(1 - self.ecc**2) * np.sin(E)
        z_prime = 0  # By definition in the orbital plane
        
        # Compute the coordinates in the ecliptic plane
        cos_w = np.cos(self.w_radians)
        sin_w = np.sin(self.w_radians)
        cos_W = np.cos(self.cen_W)
        sin_W = np.sin(self.cen_W)
        cos_i = np.cos(self.inclination)
        sin_i = np.sin(self.inclination)
        
        self.x_orbit = (cos_w * cos_W - sin_w * sin_W * cos_i) * x_prime + (-sin_w * cos_W - cos_w * sin_W * cos_i) * y_prime
        self.y_orbit = (cos_w * sin_W + sin_w * cos_W * cos_i) * x_prime + (-sin_w * sin_W + cos_w * cos_W * cos_i) * y_prime
        self.z_orbit = (sin_w * sin_i) * x_prime + (cos_w * sin_i) * y_prime
        
        self.d = np.sqrt(self.x_orbit**2 + self.y_orbit**2 + self.z_orbit**2)
        
        return self.x_orbit, self.y_orbit, self.z_orbit
    
    def current_position(self):
        self.solve_kepler_equation()
        x_orbit, y_orbit, z_orbit = self.xyz_calculation(self.E)
        #print (f"Current position, {self.name}: X: {x_orbit}   Y: {y_orbit}   Z: {z_orbit} Distance: {self.d} Time: {t_delta}")
        return x_orbit, y_orbit, z_orbit
    
delta_increment=1
t_delta = 0
object_data = reading_data()
num_points = 500  # Number of points to simulate
theta = np.linspace(0, 2 * np.pi, num_points)

def animate(i):
    global t_delta
    if t_delta==0:
        for planet_name, planet in object_data.items():
            planet.update_stats_delta(0)
            
    for planet_name, planet in object_data.items():
        x_next, y_next, z_next = planet.current_position()
        planet_objects[planet_name].set_data(x_next, y_next)
        planet_objects[planet_name].set_3d_properties(z_next)
            
    t_delta += (delta_increment*0.5)
    for planet_name, planet in object_data.items():    
        planet.update_stats_delta(t_delta)
        
    interface.update_date_annotation(Time.now().jd, t_delta)

planet_objects = {}
for planet_name, planet in object_data.items():
    x_initial, y_initial, z_initial = 0,0,0
    planet_objects[planet_name], = ax.plot([x_initial], [y_initial], [z_initial], 'o') # Initialize but don't show label

for planet_name, planet in object_data.items():
    planet.draw_orbit(theta)
    
# Create buttons
ax_start = plt.axes([0.015, 0.85, 0.05, 0.035])
ax_stop = plt.axes([0.065, 0.85, 0.05, 0.035])
ax_increase = plt.axes([0.015, 0.80, 0.05, 0.035])
ax_decrease = plt.axes([0.065, 0.80, 0.05, 0.035])
btn_start = Button(ax_start, 'Start')
btn_stop = Button(ax_stop, 'Stop')
btn_increase = Button(ax_increase, '+ Speed')
btn_decrease = Button(ax_decrease, '- Speed')

# Connect the buttons to their respective functions
btn_start.on_clicked(interface.start_animation)
btn_stop.on_clicked(interface.stop_animation)
btn_increase.on_clicked(interface.increase_speed)
btn_decrease.on_clicked(interface.decrease_speed)

sun, = ax.plot([0], [0], [0], 'o', color='yellow', label='Sun', markersize=13)
ani = FuncAnimation(fig, animate, frames=num_points, interval=500, blit=False)


interface.hide_axes()
interface.background()
interface.max_zoom()
interface.create_t_delta_text(fig)
interface.create_date_text(fig)
belts.gen_kuiper_belt()
belts.gen_main_belt()

plt.show()