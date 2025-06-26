TODO: add any details

This code attempts to implement the same ZMP process outlined in the Python implementation, but using C++. It leverages similar ideas from the OpenSim code (e.g. InverseDynamics tool etc.).

**Notes:**

- **<u>Issues with resolving for udot in the static manner vs. using accelerations generated from a Moco solution (i.e. different udot resulting in different ZMP calculations).</u>**
  - This issue gets fixed when using the data from the Moco trajectory vs. realizing to accelerations and getting udot
- Distance calculation fixed
- Random 1's in outputs fixed



**TODO:**

- Add moments in outputs?