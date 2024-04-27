#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

// Constants
const double g = 9.81; // Acceleration due to gravity (m/s^2)

// Function to calculate the forces acting on the sphere and net force
void calculateForces(double mass, double volume, double velocity, double dragCoefficient,
                     double fluidDensity, double airDensity, double crossSectionalArea,
                     double& airDragForce, double& viscousDragForce, double& buoyantForce, double& netForce) {
    // Gravitational force
    double gravitationalForce = mass * g;

    // Buoyant force
    buoyantForce = fluidDensity * volume * g;

    // Viscous drag force (proportional to velocity)
    viscousDragForce = dragCoefficient * velocity;

    // Air drag force (proportional to velocity squared)
    airDragForce = 0.5 * airDensity * crossSectionalArea * dragCoefficient * velocity * velocity;

    // Net force acting on the sphere
    netForce = gravitationalForce - buoyantForce - viscousDragForce - airDragForce;
}

// Function to simulate the motion of a falling spherical body using Euler's method
void simulateFallingSphere(double mass, double volume, double dragCoefficient, double fluidDensity,
                           double airDensity, double crossSectionalArea, double timeStep, double duration) {
    // Initial conditions
    double y = 0.0; // Initial position (m)
    double v = 0.0; // Initial velocity (m/s)
    double a = 0.0; // Initial acceleration (m/s^2)

    // Open files to save the data for each type of force
    ofstream airDragFile("air_drag_data.txt");
    ofstream viscousDragFile("viscous_drag_data.txt");
    ofstream buoyantForceFile("buoyant_force_data.txt");

    // Write headers for each file
    airDragFile << "Time Position Velocity Acceleration\n";
    viscousDragFile << "Time Position Velocity Acceleration\n";
    buoyantForceFile << "Time Position Velocity Acceleration\n";

    // Simulate the motion over the specified duration
    for (double t = 0.0; t <= duration; t += timeStep) {
        // Variables to hold the forces
        double airDragForce, viscousDragForce, buoyantForce, netForce;

        // Calculate the forces and net force acting on the sphere
        calculateForces(mass, volume, v, dragCoefficient, fluidDensity, airDensity, crossSectionalArea,
                        airDragForce, viscousDragForce, buoyantForce, netForce);

        // Calculate acceleration using Newton's second law
        a = netForce / mass;

        // Update position and velocity using Euler's method
        y += v * timeStep;
        v += a * timeStep;

        // Save the current time, position, velocity, and acceleration to each file
        airDragFile << t << " " << y << " " << v << " " << a << "\n";
        viscousDragFile << t << " " << y << " " << v << " " << a << "\n";
        buoyantForceFile << t << " " << y << " " << v << " " << a << "\n";

        // Optional: Output the current time, position, velocity, and acceleration
        cout << "Time: " << t << " s, Position: " << y << " m, Velocity: " << v << " m/s, Acceleration: " << a << " m/s^2" << endl;
    }

    // Close the files
    airDragFile.close();
    viscousDragFile.close();
    buoyantForceFile.close();
}

// Main function
int main() {
    // Variables to hold user inputs
    double mass, volume, dragCoefficient, fluidDensity, airDensity, crossSectionalArea;
    double timeStep, duration;

    // Prompt user for various parameters
  cout << "Enter mass of the sphere (kg): ";
  cin >> mass;

    cout << "Enter volume of the sphere (m^3): ";
    cin >> volume;

   cout << "Enter drag coefficient of the sphere: ";
   cin >> dragCoefficient;

    cout << "Enter fluid density (kg/m^3): ";
    cin >> fluidDensity;

    cout << "Enter air density (kg/m^3): ";
    cin >> airDensity;

    cout << "Enter cross-sectional area of the sphere (m^2): ";
    cin >> crossSectionalArea;

    cout << "Enter time step for the simulation (s): ";
    cin >> timeStep;

    cout << "Enter duration of the simulation (s): ";
    cin >> duration;

    // Run the simulation
    simulateFallingSphere(mass, volume, dragCoefficient, fluidDensity, airDensity, crossSectionalArea, timeStep, duration);

    return 0;
}

