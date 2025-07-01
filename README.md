# SLPConstructor
A package for Constructing Subsystem Lifted Product Codes from a Base Matrix

In our GNARSIL paper,  https://ieeexplore.ieee.org/document/10821364, we introduced the concept of Subsystem Lifted Product Codes or SLP codes. A More in-depth view SLP codes can be found here:

https://wp.optics.arizona.edu/alumni/wp-content/uploads/sites/113/2024/12/onovakthesis.pdf

These are exciting codes due to their QLDPC nature*. They can be seen as natural extensions of Subsystem Hypergraph Product (SHP) codes.

Our package is very simple to use:

  1. Build a Base Matrix using np arrays and Galois Polynomial Objects. The Field should always be set to GF2 inside of each object instance.

  2. Call SLP with the base matrix and Lift size L as shown. This will output the entire construction, including  Logical Generator Matrices (not in standard form), Gauge Operators, and Stabilizers.
     
  4. !!!You will need to use the logOps script provided (or another similar script) to get logical generators that form symplectic pairs from the given generator matrices LX, LZ!!!
     
  6. !!!The stabilizers will have a few linearly dependent rows necessary to ensure full rank in binary!!!
  7. !!!These past two conditions are an artifact of the method used for computing the Base Generator Matrix over the polynomial ring, which may not always be simplified to standard form over the ring!!!

  8. Have Fun!!

!!!Please Cite the GNARSIL paper if you use this package for publication!!!

* The gauges are QLDPC, Stabilizers are not.
