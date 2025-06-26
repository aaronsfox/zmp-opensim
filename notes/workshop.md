**Notes:**

- How to structure model outputs, components etc.?
  - A component is probably better
  - <u>Component Structure</u>
    - A list of contact bodies that are reviewed for ground contact
      - Body that is assessed for contact with the ground (socket, body) - name property instead of a socket?
      - A series of points on the body that are checked against a distance threshold to determine if they are in contact with the ground (list of stations) - NOPE, list of Vec3's associated with body?
        - Contact points may also be useful in identifying a polygon to constraint ZMP COP within
    - The distance to check against each of the points for relevant contact to the ground plane (integer, in metres)
    - The free body joint that connects the body to the ground (socket, joint)
    - Threshold for vertical force to be considered true contact (optional, integer)
- Determine how to manage contact identification?
  - Points relative to the ground (i.e. distance)
  - Do these need to be contact geometries or can they simply be points with a radius?
- Determine cause of difference with getUdot() vs. Moco trajectory accelerations
  - Probably due to understanding the forces being applied at that point of the state?

**To Do:**

- Create model components, outputs etc.
- Test out initial components in state by state estimates - possible to do?
- Create example scripts estimating ground reactions from kinematics
- Create force connected component and test creating external loads from these
- Create Moco goal (and constraint?)
- Test Moco goal