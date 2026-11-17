# Crazyflie Autonomous Navigation

A simulation-to-hardware project exploring path planning and autonomous navigation, beginning with a MATLAB grid-based pathfinding simulation and intended to progress to a Bitcraze Crazyflie nano-quadcopter.

> **Project status:** MATLAB pathfinding simulation implemented. Crazyflie integration is planned and has not yet been implemented.

## Current implementation: MATLAB pathfinding simulation

The current MATLAB implementation is contained in one file:

- `MATLAB_Simulation/Assignment3_Pathfinding.m`

It uses **Breadth-First Search (BFS)** on a 2D occupancy grid.

### Environment and movement model

- The map is a 2D grid: `0` represents a free cell and `1` represents a blocked cell.
- The agent moves in four directions: up, down, left, and right.
- Each move has a cost of 1.
- Diagonal movement is not allowed.
- With equal step costs, BFS finds a minimum-move route when a route exists.

### Implemented features

- Original 16 × 16 corridor map.
- Three prescribed obstacle scenarios.
- Dynamic obstacle insertion and route validation.
- Replanning only when the remaining route becomes invalid.
- Retaining the existing route when a newly added obstacle does not affect it.
- Handling of an unreachable goal.
- Animated agent movement.
- Interactive MATLAB GUI with random map generation, start/goal selection, obstacle editing, and simulation controls.
- Automated checks covering map construction, obstacle updates, route validity, replanning, route retention, and unreachable-goal handling.

## Requirements

- MATLAB.
- The source header describes the implementation as a single-file simulation with no additional toolboxes.

## Run the simulation

Open MATLAB, set the current folder to `MATLAB_Simulation`, and run one of the following commands:

| Command | Behaviour |
|---|---|
| `Assignment3_Pathfinding` | Run the original map and the three obstacle cases with animation |
| `Assignment3_Pathfinding('fast')` | Run the same scenarios without animation pauses |
| `Assignment3_Pathfinding('retain')` | Demonstrate retaining a still-valid route |
| `Assignment3_Pathfinding('unreachable')` | Run the synthetic unreachable-goal test |
| `Assignment3_Pathfinding('test')` | Run the automated verification checks |
| `Assignment3_Pathfinding('gui')` | Open the interactive GUI |

The optional agent-step argument can be supplied for the three scripted obstacle cases:

```matlab
Assignment3_Pathfinding('run', [4 4 4])
```

The values specify how many moves the agent makes before the obstacle is introduced in Cases 1–3.

## Repository structure

```text
crazyflie-autonomous-navigation/
├── README.md
└── MATLAB_Simulation/
    └── Assignment3_Pathfinding.m
```

As the project develops, additional folders can be added for Crazyflie path planning, localization, obstacle detection, flight control, hardware notes, and real-world test results.

## Planned work

The following are future work, not part of the current implementation:

- Select and configure the Crazyflie platform and required hardware.
- Choose a localization and positioning approach.
- Adapt path planning from grid coordinates to physical waypoints.
- Integrate obstacle sensing and navigation.
- Validate the system in simulation and controlled real-world tests.

## Notes

The current MATLAB code implements BFS. It does not currently implement A*, Dijkstra, 3D planning, or Crazyflie flight control.

The MATLAB simulation is based on an academic assignment. Before publishing, check that the course brief, supplied map, and any third-party materials are permitted to be shared publicly. This repository currently includes the implementation file and this project description only.
