# Unsynchronized-Pseudolite-PL-Positioning-Dataset-Kalman-Filter-Scripts
This repository provides public-release datasets and reference EKF scripts for reproducing the vehicle-borne and UAV-borne positioning results in our manuscript. The released data are organized so that reviewers can run the scripts end-to-end and verify the reported positioning performance.
## 1) Folder Structure
```text
root/
├─ observations/
│  ├─ car_observations_data.txt     # vehicle: time + code pseudorange + pseudorange-rate
│  ├─ uav_observations.txt          # UAV: time + code pseudorange
│  └─ uav_base_stations.txt         # UAV: base station positions in LLH (lat, lon, h)
├─ trajectory/
│  ├─ CAR.txt                       # vehicle trajectory (SBG/RTK log)
│  └─ UAV.txt                       # UAV trajectory (SBG/RTK log)
└─ Code/
   ├─ CarEKF.m                      # public EKF script for vehicle dataset (PR + PRR + height)
   ├─ UavEKF.m                      # public EKF script for UAV dataset (PR + height)
   └─ common/                       # shared MATLAB helpers (e.g., geo2utm_batch.m)
```
## 2) Data Files
### `observations/car_observations_data.txt`

Vehicle observations (tab-separated):

- Column 1: `timeSec[s]`
- Columns 2–6: `codePseudorange1..5 [m]`
- Columns 7–11: `pseudorangeRate1..5 [m/s]`

### `observations/uav_observations.txt`

UAV observations (tab-separated):

- Column 1: `timeSec[s]`
- Columns 2–6: `codePseudorange1..5 [m]`

### `observations/uav_base_stations.txt`

UAV base stations (tab-separated):

- Column 1: `ID`
- Column 2: `lat [deg]`
- Column 3: `lon [deg]`
- Column 4: `h [m]`

### `trajectory/CAR.txt` and `trajectory/UAV.txt`

Vehicle/UAV trajectories (SBG/RTK log files used by the scripts).
