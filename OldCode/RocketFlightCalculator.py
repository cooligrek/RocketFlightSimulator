# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
# Press Ctrl+F8 to toggle the breakpoint.

import matplotlib.pyplot as plt

def flight_calc():
    print('Calculating...\n')

    timeArray: list = []
    accArray: list = []
    velArray: list = []
    heightArray: list = []
    fuelMassArray: list = []
    oxidizerMassArray: list = []
    totalMassArray: list = []
    desnityArray: list = []
    dragArray: list = []

    maxAcc: float = 0
    maxVel: float = 0
    maxHeight: float = 0

    timeAtMaxAcc: float = 0
    timeAtMaxVel: float = 0
    timeAtMaxHeight: float = 0

    dry_mass: float = 200  # lbs
    fuel_mass: float = 16.39  # lbs
    oxidizer_mass: float = 29.48  # lbs

    fuel_mass_flowrate: float = 1.49  # lbs/s
    oxidizer_mass_flowrate: float = 2.68  # lbs/s
    burn_time: float = 11  # s

    gravity: float = 32.17  # ft/s^2
    atmospheric_density: float = 0.0765  # lbs/ft^3

    frontal_area_Rocket: float = 0.3490  # ft^2
    frontal_area_Parachute_Pilot: float = 78.53  # ft^2
    frontal_area_Parachute_Droug1: float = 78.53  # ft^2
    frontal_area_Parachute_Droug2: float = 78.53  # ft^2
    frontal_area_Parachute_Main: float = 314.16  # ft^2

    coeff_drag_Rocket: float = 0.4
    coeff_drag_Parachute_Pilot: float = 0.97
    coeff_drag_Parachute_Droug1: float = 2.2
    coeff_drag_Parachute_Droug2: float = 2.2
    coeff_drag_Parachute_Main: float = 2.2

    drag_rocket: float = 0  # N
    drag_pilot: float = 0  # N
    drag_drog1: float = 0  # N
    drag_drog2: float = 0  # N
    drag_main: float = 0  # N

    force_thrust: float = 999  # lbs
    force_gravity: float = -(dry_mass + fuel_mass + oxidizer_mass)  # lbs
    force_drag: float = 0  # lbs

    acceleration: float = 0  # ft/s^2
    velocity: float = 0  # ft/s
    pastHeight: float = 0  # ft
    height: float = 0  # ft
    time: float = 0  # ms

    while(height >= 0):

        burn_time = loss_calculator(burn_time, 1/1000)
        #print(burn_time)

        fuel_mass = loss_calculator(fuel_mass, fuel_mass_flowrate/1000)
        oxidizer_mass = loss_calculator(oxidizer_mass, oxidizer_mass_flowrate/1000)
        masses = [dry_mass, fuel_mass, oxidizer_mass]
        #print(masses)

        # rough approximation equation of atmospheric density good up to 50k feet
        atmospheric_density = -0.05440907 + ((0.07646961 - -0.05440907) / (1 + (height / 50621.47)**1.055246))
        #print(atmospheric_density)

        if (maxAcc < acceleration):
            maxAcc = acceleration
            timeAtMaxAcc = time / 1000

        if (maxVel < velocity):
            maxVel = velocity
            timeAtMaxVel = time / 1000

        if (maxHeight < height):
            maxHeight = height
            timeAtMaxHeight = time / 1000

        drag_rocket = drag_calculator(velocity, atmospheric_density, gravity, coeff_drag_Rocket, frontal_area_Rocket)

        if(velocity < 0):
            drag_pilot = drag_calculator(velocity, atmospheric_density, gravity, coeff_drag_Parachute_Pilot, frontal_area_Parachute_Pilot)

            if(time / 1000 > timeAtMaxHeight + 30):
                drag_drog1 = drag_calculator(velocity, atmospheric_density, gravity, coeff_drag_Parachute_Droug1, frontal_area_Parachute_Droug1)
                drag_drog2 = drag_calculator(velocity, atmospheric_density, gravity, coeff_drag_Parachute_Droug2, frontal_area_Parachute_Droug2)

            if (height < 1000):
                drag_main = drag_calculator(velocity, atmospheric_density, gravity, coeff_drag_Parachute_Main, frontal_area_Parachute_Main)

        drags = [drag_rocket, drag_pilot, drag_drog1, drag_drog2, drag_main]
        print(str(drags) + "\n")

        #print(coeff_dragArray)
        #print(frontal_areaArray)

        if (burn_time <= 0):
            force_thrust = 0
        force_gravity = -1 * sum(masses)
        force_drag = sum(drags)
        forces = [force_gravity, force_drag, force_thrust]
        #print(forces)


        acceleration = sum(forces) / (sum(masses)/gravity)
        velocity = velocity + acceleration/1000
        pastHeight = height
        height = height + velocity/1000
        time = time + 1

        print('Time (s): ' + str(time/1000) +
              ' | Acceleration (ft/s^2): ' + str(acceleration) +
              ' | Velocity (ft/s): ' + str(velocity) +
              ' | Hight (ft): ' + str(height))

        timeArray.append(time/1000)
        accArray.append(acceleration)
        velArray.append(velocity)
        heightArray.append(height)
        fuelMassArray.append(fuel_mass)
        oxidizerMassArray.append(oxidizer_mass)
        totalMassArray.append(force_gravity)
        desnityArray.append(atmospheric_density)
        dragArray.append(force_drag)
    
    print('Max Acceleration (ft/s^2): ' + str(maxAcc) + ' | Time at max Acceleration (s): ' + str(timeAtMaxAcc))
    print('Max Velocity (ft/s^2): ' + str(maxVel) + ' | Time at max Velocity (s): ' + str(timeAtMaxVel))
    print('Max Height (ft/s^2): ' + str(maxHeight) + ' | Time at max Height (s): ' + str(timeAtMaxHeight))
    
    figure, grid = plt.subplots(2,5)

    grid[0, 0].plot(timeArray, accArray)
    grid[0, 0].set_title('Acceleration vs Time')

    grid[0, 1].plot(timeArray, velArray)
    grid[0, 1].set_title('Velocity vs Time')

    grid[0, 2].plot(timeArray, heightArray)
    grid[0, 2].set_title('Height vs Time')

    grid[1, 1].plot(timeArray, fuelMassArray)
    grid[1, 1].set_title('Fuel Mass vs Time')

    grid[1, 2].plot(timeArray, oxidizerMassArray)
    grid[1, 2].set_title('Oxidizer Mass vs Time')

    grid[1, 0].plot(timeArray, totalMassArray)
    grid[1, 0].set_title('Total Mass vs Time')

    grid[0, 3].plot(heightArray, desnityArray)
    grid[0, 3].set_title('Atmospheric Density vs Height')

    grid[1, 3].plot(timeArray, desnityArray)
    grid[1, 3].set_title('Atmospheric Density vs Time')

    grid[0, 4].plot(timeArray, dragArray)
    grid[0, 4].set_title('Drag vs Time')

    plt.show()

def loss_calculator(input, delta):

    if(input > 0):
        output = input - delta
    else:
        output = 0

    return output

def drag_calculator(velocity, density, gravity, area, coef):
    return -0.5 * density * velocity * abs(velocity) * area * coef / gravity


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    flight_calc()
