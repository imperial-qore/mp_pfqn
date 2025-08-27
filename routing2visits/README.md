# Routing2Visits Utility

A standalone utility that converts routing probability matrices and service times into QN format for queueing network analysis.

## Features

- Exact rational arithmetic using GMP library
- Calculates mean number of visits from routing matrices
- Outputs in standard QN format compatible with queueing network tools
- Supports multiple classes and stations

## Building

```bash
make all
```

This creates:
- `routing2visits` - standalone utility
- `librouting2visits.a` - static library for integration

## Usage

```bash
./routing2visits <num_stations> <num_classes> <routing_prefix> <service_times_file> [population1] [population2] ...
```

### Arguments

- `num_stations` - Total number of stations (including delay station)
- `num_classes` - Number of customer classes  
- `routing_prefix` - Prefix for routing matrix files (e.g., 'P' for P1.txt, P2.txt, ...)
- `service_times_file` - File containing service times matrix
- `population1, ...` - Population for each class (optional, defaults to 0)

### Input File Formats

#### Routing Matrix Files (P1.txt, P2.txt, ...)

Each file contains the routing probability matrix for one class as GMP rationals:

```
0/1 1/3 2/3
1/2 0/1 1/2  
1/2 1/2 0/1
```

Matrix dimensions: (num_stations × num_stations)
- Row/column 0: delay station
- Rows/columns 1 to M: queue stations

#### Service Times File

Contains service times for each station and class as GMP rationals:

```
10/1 20/1
5/1  8/1
3/1  12/1  
```

Matrix dimensions: (num_stations × num_classes)
- Row 0: think times at delay station
- Rows 1 to M: service times at queue stations

### Output Format

The utility outputs in QN format:

```
R                    # Number of classes
N1 N2 ... NR        # Population for each class
Z1 Z2 ... ZR        # Think times (service_demand × visits)
M                    # Number of queues (excluding delay)
mi1 L11 L12 ... L1R  # Server count and service demands for queue 1
mi2 L21 L22 ... L2R  # Server count and service demands for queue 2
...
miM LM1 LM2 ... LMR  # Server count and service demands for queue M
```

Where service demands L[m][r] = service_time[m][r] × visits[r][m]

## Example

```bash
# Example: 3 stations (1 delay + 2 queues), 2 classes
./routing2visits 3 2 P service_times.txt 0 0
```

This reads:
- P1.txt, P2.txt (routing matrices for classes 1,2)
- service_times.txt (service times)
- Sets population to 0 for both classes

## Mathematical Background

The utility solves the visit equation system:
- v[r][0] = 1 (delay station reference)  
- v[r][m] = Σ(k=0 to M) v[r][k] × P[r][k][m] for m = 1 to M

Where:
- v[r][m] = mean visits of class r to station m
- P[r][k][m] = routing probability from station k to m for class r

## Dependencies

- GMP library for exact rational arithmetic
- Standard C library

## Integration

Include `routing2visits.h` and link with `librouting2visits.a` and `-lgmp`.