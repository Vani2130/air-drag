#include <iostream>
#include <cmath>
#include <fstream>

// Constants
const double g = 9.81; // Acceleration due to gravity (m/s^2)

// Function to calculate the net force acting on the sphere
double calculateNetForce(double mass, double volume, double velocity, double dragCoefficient,
                         double fluidDensity, double airDensity, double crossSectionalArea) {
    // Gravitational force
    double gravitationalForce = mass * g;

    // Buoyant force
    double buoyantForce = fluidDensity * volume * g;

    // Viscous drag force (proportional to velocity)
    double viscousDragForce = dragCoefficient * velocity;

    // Air drag force (proportional to velocity squared)
    double airDragForce = 0.5 * airDensity * crossSectionalArea * dragCoefficient * velocity * velocity;

    // Net force acting on the sphere
    double netForce = gravitationalForce - buoyantForce - viscousDragForce - airDragForce;
    return netForce;
}

// Function to simulate the motion of a falling spherical body using Euler's method
void simulateFallingSphere(double mass, double volume, double dragCoefficient, double fluidDensity,
                           double airDensity, double crossSectionalArea, double timeStep, double duration) {
    // Initial conditions
    double y = 0.0; // Initial position (m)
    double v = 0.0; // Initial velocity (m/s)

    // Open a file to save the data (space-separated)
    std::ofstream dataFile("motion_data.txt");
    dataFile << "Time Position Velocity\n"; // Write headers for better readability

    // Simulate the motion over the specified duration
    for (double t = 0.0; t <= duration; t += timeStep) {
        // Calculate net force acting on the sphere
        double netForce = calculateNetForce(mass, volume, v, dragCoefficient, fluidDensity, airDensity, crossSectionalArea);

        // Calculate acceleration using Newton's second law
        double acceleration = netForce / mass;

        // Update position and velocity using Euler's method
        y += v * timeStep;
        v += acceleration * timeStep;

        // Save the current time, position, and velocity to the file (space-separated)
        dataFile << t << " " << y << " " << v << "\n";

        // Output the current time, position, and velocity
        std::cout << "Time: " << t << " s, Position: " << y << " m, Velocity: " << v << " m/s" << std::endl;
    }

    // Close the file
    dataFile.close();
}

// Main function
int main() {
    // Variables to hold user inputs
    double mass, volume, dragCoefficient, fluidDensity, airDensity, crossSectionalArea;
    double timeStep, duration;

    // Prompt user for various parameters
    std::cout << "Enter mass of the sphere (kg): ";
    std::cin >> mass;

    std::cout << "Enter volume of the sphere (m^3): ";
    std::cin >> volume;

    std::cout << "Enter drag coefficient of the sphere: ";
    std::cin >> dragCoefficient;

    std::cout << "Enter fluid density (kg/m^3): ";
    std::cin >> fluidDensity;

    std::cout << "Enter air density (kg/m^3): ";
    std::cin >> airDensity;

    std::cout << "Enter cross-sectional area of the sphere (m^2): ";
    std::cin >> crossSectionalArea;

    std::cout << "Enter time step for the simulation (s): ";
    std::cin >> timeStep;

    std::cout << "Enter duration of the simulation (s): ";
    std::cin >> duration;

    // Run the simulation
    simulateFallingSphere(mass, volume, dragCoefficient, fluidDensity, airDensity, crossSectionalArea, timeStep, duration);

    return 0;
}
