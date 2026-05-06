# Presentation Script — RIS-Assisted CubeSat-to-Ground Communication

**Total estimated time:** 12–14 minutes (excluding Q&A).

Below is a slide-by-slide script with what to say and what to show. Each slide has a target time and key visual.

---

## Slide 1 — Title (30 seconds)

**Visual:** Project title, your name, institution, advisor (if any), date.

**What to say:**

> Good morning everyone. My project is titled "RIS-Assisted CubeSat-to-Ground Communication: Algorithmic Optimization for Low-Elevation Satellite Links." I'm APache from MCTE. Today I'll walk you through why low-elevation satellite communication is a hard problem, how a Reconfigurable Intelligent Surface — basically a programmable smart mirror for radio waves — can help solve it, and what my MATLAB simulations show about the performance gain we can achieve through smart algorithms alone, with no hardware changes required.

---

## Slide 2 — The CubeSat Communication Challenge (1 minute)

**Visual:** Diagram of a CubeSat orbiting Earth, with a ground station antenna. Show the satellite's arc across the sky, marking zenith and horizon.

**What to say:**

> CubeSats are small satellites, roughly the size of a shoebox, used for Earth observation, IoT networks, and scientific missions. They orbit in Low Earth Orbit, meaning they zip across the sky in just a few minutes per pass — typically 5 to 12 minutes of visibility from any one ground station.
>
> The problem: link quality is not constant during the pass. When the satellite is directly overhead — at zenith — the link is excellent. But as the satellite approaches the horizon, communication quality collapses. We can lose 10 to 20 dB of signal strength at low elevation angles compared to zenith. Effectively, we waste several minutes of every pass.

---

## Slide 3 — Why Low Elevation Is Hard (1 minute)

**Visual:** Annotated diagram showing four problems at low elevation: longer atmospheric path, terrain obstruction, multipath reflections, higher Doppler.

**What to say:**

> Four physical reasons make low elevation difficult:
>
> First, the radio wave travels through much more atmosphere when the satellite is near the horizon, so atmospheric absorption is significantly worse.
>
> Second, the line-of-sight is often blocked by buildings, trees, or terrain.
>
> Third, the signal bounces off ground clutter, arriving from multiple directions with different delays — this is multipath fading, and it causes destructive interference at the receiver.
>
> Fourth, the Doppler shift is largest near the horizon, complicating receiver synchronization.
>
> The combined effect is a high bit error rate that wastes the contact window. Our project addresses this purely through algorithms — without launching new satellites or building bigger antennas.

---

## Slide 4 — What Is an RIS? (1.5 minutes)

**Visual:** Image of an RIS panel — a flat surface with a grid of small reflecting elements. Side-by-side: a passive metal mirror vs an RIS with electronically tuned elements.

**What to say:**

> The Reconfigurable Intelligent Surface — RIS for short — is the central idea. Picture a flat panel covered with hundreds of small elements. Each element can be electronically tuned to apply a specific phase shift to incoming radio waves. By carefully choosing the phase on every element, the panel acts as a beam-steerable mirror. It can take an incoming radio wave from one direction and reflect it precisely toward any chosen direction, with all reflected components arriving in phase at the target.
>
> Three things make the RIS powerful:
>
> One, it's passive. It doesn't generate signals or amplify them, so it consumes very little power.
>
> Two, it's reconfigurable in microseconds. As the satellite moves across the sky, the phase pattern can be updated in real time.
>
> Three, performance scales with the number of elements. Adding more elements gives you more beamforming gain.
>
> In our setup, the RIS is mounted near the ground station — say, on a nearby building. The CubeSat's signal reaches the receiver via two paths now: a direct path, which is weak at low elevation, and a reflected path through the RIS. When the RIS phases are correctly tuned, the reflected path adds coherently to the direct path, dramatically improving received signal strength.

---

## Slide 5 — System Model (1 minute)

**Visual:** Block diagram. CubeSat (transmitter) → two arrows: one direct to GS, one via RIS panel with N elements → receiver. Label all channels: h_d, h_sr, h_rg.

**What to say:**

> Here's the system model. The CubeSat transmits a symbol x. It reaches the ground station through three channels: a direct channel h_d; a satellite-to-RIS channel h_sr, which is a vector of N coefficients, one per element; and an RIS-to-ground channel h_rg, also N coefficients.
>
> The RIS applies a diagonal matrix of phase shifts — one phase per element. The received signal is the direct path plus the reflected path summed across all N elements, multiplied by the transmitted symbol, plus thermal noise.
>
> Our optimization variable is the vector of phase angles. Our objective is to maximize the magnitude of the effective channel — because a larger effective channel gain means a stronger received signal, which means a lower bit error rate.

---

## Slide 6 — Channel Model: Rician Fading (1 minute)

**Visual:** Two panels showing Rician distribution. Left: high K-factor — narrow distribution around the LOS peak. Right: low K-factor — wide distribution.

**What to say:**

> We use the Rician fading model, which captures channels that have both a strong line-of-sight component and a weaker scattered component. The Rician K-factor is the ratio of LOS power to scattered power. High K means deterministic line-of-sight dominates; low K means the channel is heavily affected by multipath.
>
> Crucially, we make the K-factor depend on elevation angle: at zenith, K is high — clean line-of-sight; at the horizon, K is low — multipath dominates. This single modeling choice is what mechanically captures why low-elevation links are harder, and it's what lets us quantify the RIS benefit specifically in that regime.

---

## Slide 7 — Optimization Problem (1 minute)

**Visual:** Equation: maximize |h_d + sum(h_sr_n · h_rg_n · exp(j·θ_n))|² over θ ∈ [0, 2π)^N. Then show the closed-form solution: θ_n* = angle(h_d) − angle(h_sr_n · h_rg_n).

**What to say:**

> The optimization problem is to find the N phase angles that maximize the effective channel magnitude.
>
> For this idealized single-input single-output case, there's actually a closed-form solution: set each element's phase so that its reflected contribution is aligned with the direct path. We use this as a theoretical upper bound — the absolute best any algorithm can do under perfect channel knowledge.
>
> But in practice, you don't always have perfect knowledge, you might have hardware that only supports a handful of discrete phase levels, and you may have multiple users to serve. So we need an iterative algorithm that generalizes. That's where Particle Swarm Optimization comes in.

---

## Slide 8 — Particle Swarm Optimization (1.5 minutes)

**Visual:** Animation or schematic of a swarm of particles converging in a 2D landscape toward a peak.

**What to say:**

> PSO is a metaheuristic inspired by bird flocking. We initialize a swarm of, say, 20 particles, each one a random candidate phase vector. Each iteration, every particle remembers its own best-found position and the swarm's globally best position. Each particle updates its velocity based on inertia, attraction toward its personal best, and attraction toward the global best. After enough iterations, the swarm converges on a high-quality solution.
>
> Why PSO instead of the closed-form solution? Because PSO generalizes naturally to discrete phase shifts, imperfect channel state information, and multi-user constraints — all the things real hardware actually deals with. The closed-form is just our reference ceiling.

---

## Slide 9 — Methodology (1 minute)

**Visual:** Flowchart: parameters → channel generation → four phase configurations → BPSK transmission → demodulation → error counting → BER calculation.

**What to say:**

> The simulation methodology is Monte Carlo. For each SNR value, we run hundreds of independent trials. In each trial, we generate fresh random channel realizations from the Rician model. We compute four phase configurations: no RIS (baseline); random RIS, which is our lower bound; PSO-optimized RIS; and the closed-form optimal RIS, our upper bound.
>
> For each configuration, we transmit a thousand BPSK symbols, add Gaussian noise scaled to the target SNR, demodulate coherently, and count errors. Averaging across all trials gives us the bit error rate for that configuration at that SNR. Repeat across the SNR range, and we get our headline plot.

---

## Slide 10 — Results: BER vs SNR (1.5 minutes)

**Visual:** The BER vs SNR plot from `main.m`. Four curves clearly separated. Annotate the SNR gap at BER = 1e-3.

**What to say:**

> Here are the results at 10 degrees elevation with 64 RIS elements. The red curve is no-RIS — it improves slowly with SNR because the direct channel is weak.
>
> The black dashed curve is random RIS — barely better than no RIS. Random reflections don't add coherently.
>
> The blue curve is PSO-optimized RIS. It tracks the green curve very closely. The green curve is the closed-form upper bound.
>
> At a target BER of 10⁻³, the PSO-RIS scheme requires roughly 12 to 15 dB less SNR than no-RIS to achieve the same error rate. That's a massive practical improvement — equivalent to having a 30-times-stronger transmitter, but achieved purely through smart software control of a passive surface.

---

## Slide 11 — PSO Convergence (1 minute)

**Visual:** PSO convergence plot from `pso_convergence.m`. Multiple trials converging toward the dashed optimum line.

**What to say:**

> One question we need to answer: how fast does PSO actually converge, and how close does it get to the theoretical optimum?
>
> This figure shows five independent PSO runs on the same channel realization. Within roughly 30 iterations, the swarm reaches over 95 percent of the theoretical optimum fitness. By 50 iterations, the gap is negligible. So PSO is computationally affordable: a 30-iteration run with a 20-particle swarm is feasible to execute within a single channel coherence interval, even for fast-moving LEO satellites.

---

## Slide 12 — Performance vs Elevation (1 minute)

**Visual:** BER vs elevation plot from `elevation_sweep.m`. The gap between no-RIS and PSO-RIS curves is largest at low elevation and shrinks at high elevation.

**What to say:**

> This is arguably the most important result. We sweep elevation from 5 to 85 degrees at fixed SNR. Without RIS, BER is high at low elevation and improves as the satellite rises. With PSO-RIS, BER stays low across the entire elevation range — including the low-elevation regime where the no-RIS link is essentially unusable.
>
> The headline takeaway: RIS doesn't just improve average performance, it specifically rescues the low-elevation portion of every satellite pass. That extends the usable contact window — which in CubeSat operations directly translates to more downloaded data per orbit.

---

## Slide 13 — Limitations and Future Work (45 seconds)

**Visual:** Bullet list.

**What to say:**

> A few honest limitations. We assume continuous phase shifts; real hardware uses 1-bit or 2-bit quantization. We assume perfect channel state information, which is optimistic. And we treat the channel as quasi-static within each Monte Carlo run, ignoring Doppler dynamics.
>
> Future work: I want to add quantized-phase analysis to quantify the hardware-realistic gap, replace PSO with Deep Reinforcement Learning using DDPG or SAC for online phase control, and incorporate the Lutz Markov state model for shadowing transitions.

---

## Slide 14 — Conclusions (45 seconds)

**Visual:** Three-bullet summary.

**What to say:**

> Three takeaways.
>
> First, low-elevation CubeSat links are a real bottleneck and an algorithmic problem worth solving.
>
> Second, an RIS controlled by Particle Swarm Optimization can deliver 10 to 15 dB of SNR gain at low elevation, recovering most of the contact window that's currently wasted.
>
> Third, this is achievable with no satellite-side hardware modification and a single passive panel near the ground station — a deployable, scalable approach to a hard problem in non-terrestrial networks.
>
> Thank you. I'm happy to take questions.

---

## Anticipated Q&A

**Q: Why not just use a bigger antenna at the ground station?**
> Bigger antennas help on the direct path, but they don't address multipath or terrain blockage. RIS adds a controllable second path and works alongside any existing antenna.

**Q: Does the RIS need to know the channel?**
> Yes. The phase optimization assumes channel state information. In practice, this is obtained through pilot signals between the satellite and ground station with the RIS in known reference configurations. Imperfect CSI is one of my future-work items.

**Q: How does this compare to active relays?**
> An active relay amplifies signals and consumes power proportional to its gain. RIS is passive, much lower power, simpler hardware, and doesn't add receiver-side noise. The tradeoff is that RIS only redirects existing energy; it doesn't add power. For our scenario, the tradeoff favors RIS.

**Q: What about the Doppler shift?**
> Within a single coherence interval — typically a few milliseconds — the channel is approximately static and PSO can converge. The phase configuration is then updated for each new interval. Doppler-aware optimization is a future extension.

**Q: Why PSO and not gradient descent?**
> The optimization landscape is non-convex with many local optima, especially with quantized phases. PSO handles this robustly without requiring gradient information. Gradient methods do work but need careful initialization.

**Q: Have you compared with discrete-phase RIS?**
> Not yet — that's the immediate next experiment. Literature suggests roughly 1 to 2 dB loss for 2-bit quantization compared to continuous, which is acceptable.
