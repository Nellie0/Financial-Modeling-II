This is a Monte Carlo simulator for European options that accepts using Antithetic Variance and/or Van der Corput sequences.
The user will be asked for the number of steps and simulations they want for the Monte Carlo simulator.
They will then be asked if they wish to use Antithetic Variance Reduction.
And finally, they will be asked if they wish to use Van der Corput sequences for quasi random numbers. If yes, then they will be asked for the base and index values.

Start the simulation by typing "dotnet run".
The simulation will return the option price, greeks, and standard error (the standard error for Van der Corput will be zero).
