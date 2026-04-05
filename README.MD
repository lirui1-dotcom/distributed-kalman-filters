## High-pass / Low-pass Distributed Kalman Filter

Implementation of distributed Kalman Filter algorithms for multi-agent systems and sensor networks using MATLAB/Simulink.

![MATLAB](https://img.shields.io/badge/MATLAB-R2022b%2B-blue?logo=mathworks) ![Simulink](https://img.shields.io/badge/Simulink-v2022b%2B-orange?logo=simulink)

### Overview

The Distributed Kalman Filter extends the information form of the Kalman filter by incorporating a **Dynamic Average Consensus (DAC)** algorithm to share information across distributed agents. There are two main variants:
- **High-pass DAC**: More robust and theoretically grounded
- **Low-pass DAC**: Ideal for noise reduction but network-density dependent

### How It Works

**Dynamic Average Consensus** is an extension of Static Average Consensus (SAC). While SAC drives states to the average of initial values, DAC extends this to time-varying signals. By reformulating the information form Kalman filter equations, DAC can be integrated to enable distributed estimation across multiple agents.

**High-pass DAC** is more practical and theoretically robust, with a solid foundation independent of network density ([Tutorial on Dynamic Average Consensus](#references)).
 
**Low-pass DAC** is simpler to implement but relies on network density (number of nodes), which can limit its applicability in sparse networks.

### References

**Static Average Consensus:**
- Consensus and Cooperation in Networked Multi-Agent Systems
- Consensus Protocols for Networks of Dynamic Agents

**Low-pass Dynamic Average Consensus:**
- Consensus Filters for Sensor Networks and Distributed Sensor Fusion

**High-pass Dynamic Average Consensus:**
- DYNAMIC CONSENSUS FOR MOBILE NETWORKS
- Tutorial on Dynamic Average Consensus

### Testing Network Density Effects

To observe how network density affects performance, open `LowPass.slx` and modify the number of nodes in:  
`Modelling → Model Settings → Model properties → Callbacks → InitFcn`