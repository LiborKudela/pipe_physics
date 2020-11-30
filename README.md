# Efficient integration of machine learning into district heating predictive models
Suplementary materials

Usage:
```console
# generate mesh etc.
python3 GenerateTopology.py

# generate training data
mpirun -n 40 python3 SimulationTraining.py

# generate test data
mpirun -n 40 python3 SimulationTest.py

# run optimization using lofi
mpirun -n 40 python3 optimize.py
```
