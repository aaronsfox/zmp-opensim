TODO: add any details

The likely easiest way to implement the ZMP calculatgions on an existing motion would be through the analyze tool, so this code attempts to leverage similar ideas to other pilot tests around using the Inverse Dynamics solver on an existing motion. There are some easier steps here given that the entire motion is prescribed, so it should be easier.

- The `test.cpp` file is an initial test of the code to be implemented in an analysis plugin
  - Note that it currently doesn't work in allocating forces appropriately to the feet, as the code to position the model at the appropriate time-stamp isn't written - so each time it checks for force allocation the model remains in the default pose. The theory is that this will work better when used in the analyze tool. The forces look like they are coming out appropriately though.