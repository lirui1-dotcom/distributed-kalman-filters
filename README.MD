## High-pass / Low-pass Distributed Kalman Filter

Implementation of distributed Kalman Filter algorithms for multi-agent systems and sensor networks using MATLAB/Simulink.

![MATLAB](https://img.shields.io/badge/MATLAB-R2022b%2B-blue?logo=mathworks) ![Simulink](https://img.shields.io/badge/Simulink-v2022b%2B-orange?logo=simulink)

### Overview


The Distributed Kalman Filter (DKF) builds on the information form of the Kalman filter (IKF), incorporating two **Dynamic Average Consensus (DAC)** filters to share information across sensors nodes in a fully distributed way. There are two main DAC approaches:
- **High-pass DAC**: More robust and theoretically grounded
- **Low-pass DAC**: Ideal for noise reduction but network-density dependent

### How It Works

**Dynamic Average Consensus (DAC)** is an extension of the Static Average Consensus (SAC). While SAC drives states to the average of initial values, DAC extends this to track the average of N time-varying signals. By premultiplying the information form Kalman filter update equations by 1/N, DAC can be seamlessly integrated into IKF for distributed state estimation.


### References

**Static Average Consensus:**
- Consensus and Cooperation in Networked Multi-Agent Systems
- Consensus Protocols for Networks of Dynamic Agents

**Low-pass Dynamic Average Consensus:**
- Consensus Filters for Sensor Networks and Distributed Sensor Fusion

**High-pass Dynamic Average Consensus:**
- DYNAMIC CONSENSUS FOR MOBILE NETWORKS
- Tutorial on Dynamic Average Consensus

### Testing Network Density Effects for LP-DAC

To observe how network density affects performance, open `LowPass.slx` and modify the number of nodes in:  
`Modelling → Model Settings → Model properties → Callbacks → InitFcn`