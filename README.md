Version 1.0
Created 11.02.2022
by Jari Vepsäläinen

# High-speed delivery robot powertrain design using synthetic driving cycles
Electric vehicle powertrain dimensioning approach for a high-speed delivery robot.
The powertrain components of interest are the motor and the battery.
Synthetic driving cycles are generated based on a [Wolt dataset](https://github.com/woltapp/data-science-summer-intern-2021) of deliveries .
The Wolt dataset contains start and end coordinates, and the shortest pedestrian routes are acquired with Open Steertmap tool [OSMnx](https://github.com/gboeing/osmnx).
Elevation data is also added, which is provided by the [Open Topo Data API](https://www.opentopodata.org/).

## Required packages
```
pip install -r requirements.txt
```

## Please reference if you use it in your work:

Vepsäläinen, Jari. 2022. "[Energy Demand Analysis and Powertrain Design of a High-Speed Delivery Robot Using Synthetic Driving Cycles](https://www.mdpi.com/1996-1073/15/6/2198)". Energies 15, no. 6: 2198. https://doi.org/10.3390/en15062198

### BibTex

```
@article{vepsalainen2022,
  title={Energy Demand Analysis and Powertrain Design of a High-Speed Delivery Robot Using Synthetic Driving Cycles},
  author={Veps{\"a}l{\"a}inen, Jari},
  journal={Energies},
  volume={15},
  number={6},
  pages={21},
  year={2022},
  publisher={Multidisciplinary Digital Publishing Institute}
}
```



