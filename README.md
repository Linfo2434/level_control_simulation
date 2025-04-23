# Level Control Simulation

Provides a modular Python simulation for control systems. Integrates nonlinear dynamics with a 4th‐order Runge–Kutta solver under random inlet disturbances, plotting level, rate, control action, and flow. Ideal for education and research.

This project simulates a tank level control system using five different strategies:
- **P (Proportional)**
- **PI (Proportional–Integral)**
- **Super-Twisting Sliding Mode Control**
- **Predefined-Time Control**
- **Sliding Mode Control (SMC)**

Each controller is tested under random inlet‐flow disturbances, integrated with a 4th-order Runge–Kutta solver, and visualized via plots of:
1. Tank level \(h(t)\)  
2. Level rate \(\dot h(t)\)  
3. Control signal \(u(t)\)  
4. Inlet/outlet flows \(F_{\text{in}}(t)\) & \(F_{\text{out}}(t)\)
