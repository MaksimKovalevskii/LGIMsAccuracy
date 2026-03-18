# LGIMsAccuracy
This study conducts a systematic comparison of six methods: conventional second-order differential equations and Newton-Euler angular velocity-based equations, each implemented with both Euler parameters and Cartesian rotation vectors, alongside corresponding Lie group integration (LGIM) variants using RKMK integration.

Three problems are being analyzed (each one is in separate folder):
1. **Pendulum** - constrained pendulum exhibiting smooth, large rotations
<br> Read *READMEPendulum.txt* to make comparison on Pendulum problem

2. **Tennis racket phenomena** also known as the intermediate axis theorem or Dzhanibekov effect,
describes the unstable rotational behavior of a rigid body about its intermediate principal axis.
<br> Read *READMETennis.txt* to make comparison on Tennis racket phenomena

3. **Heavy Top** - classical benchmark problem in multibody system dynamics that combines constrained motion with fast rotations and rapid rotational changes.
<br> Read *READMEHeavyTop.txt* to make comparison on Heavy Top problem

System used for development and testing:
- Operating system: Windows 10 (64‑bit)
- Machine type: Work laptop, multi‑core CPU
- MATLAB: R2024b (64‑bit), standard desktop installation
- Disk: At least 25 GB free space recommended for full batch
