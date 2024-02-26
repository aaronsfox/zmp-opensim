TODO: add any details

This code attempts to implement the same ZMP process outlined in the Python implementation, but using C++. It leverages similar ideas from the OpenSim code (e.g. InverseDynamics tool etc.).

**Notes:**

- **<u>Issues with resolving for udot in the static manner vs. using accelerations generated from a Moco solution (i.e. different udot resulting in different ZMP calculations).</u>**
  - This issue gets fixed when using the data from the Moco trajectory vs. realizing to dynamics and getting udot
- Still an issue with the COP allocation to right or left body
  - Using a simpler approach of looking at which body is closer to the ground (y-axis lower) seems to work better
  - Seems like there is then an issue with the distance to CoP calculations
- There are some random 1's added to certain columns in ZMP forces for some reason?