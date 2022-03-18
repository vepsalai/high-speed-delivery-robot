Version 1.0
11.02.2022
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
[Paper submitted, not published]



