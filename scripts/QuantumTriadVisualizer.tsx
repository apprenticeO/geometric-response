import React, { useState, useEffect, useRef } from 'react';
import { Play, Pause, RotateCcw, Zap, Droplet, Link, Maximize2, Minimize2 } from 'lucide-react';

const QuantumTriadVisualizer = () => {
  const canvasRef = useRef(null);
  const [isRunning, setIsRunning] = useState(true);
  // Internal time (ref to avoid re-render/reset jitter)
  const timeRef = useRef(0);
  const [params, setParams] = useState({
    N: 8,
    coupling: 1.0,
    speed: 1.0
  });
  // Gravity/response parameters
  const [resp, setResp] = useState({
    tauG: 2.0,       // s
    G0: 1.0,         // dimensionless DC scale (visual proxy)
    wxx: 1.0,
    wyy: 1.0,
    wxy: 1.0,
  });
  const [metrics, setMetrics] = useState({
    fluctuation: 0,
    entropy: 0,
    correlation: 0,
    triad: 0
  });
  // Gravity diagnostics
  const [grav, setGrav] = useState({
    omegaG: 0,
    lG: 0,
    mG: 0,
    residual: 0,
    residualEMA: 0,
    slowBandOK: true,
  });
  const [fullscreen, setFullscreen] = useState(true);
  const [panelOpen, setPanelOpen] = useState(true);
  // Visualization scale factor (multiplies nominal size)
  const [vizScale, setVizScale] = useState(1.2);
  const fitToScreen = () => {
    // Choose a comfortable scale that keeps margins and room for the panel
    setPanelOpen(true);
    setVizScale(1.4);
  };
  const setNumElements = (n: number) => {
    const nInt = Math.max(3, Math.min(24, Math.round(n)));
    setParams({ ...params, N: nInt });
    // force re-init of particle layout next frame
    particles.current = [];
  };

  // Particle system for quantum state
  const particles = useRef([]);
  const connections = useRef([]);
  // Debye filter state for tensor components
  const debyeState = useRef<{ y_xx: number; y_yy: number; y_xy: number }>({ y_xx: 0, y_yy: 0, y_xy: 0 });
  // Export buffer (ring) for tensor time series
  const exportBuf = useRef<
    Array<{ t: number; Rinfo_xx: number; Rinfo_yy: number; Rinfo_xy: number; R_xx: number; R_yy: number; R_xy: number; triad: number }>
  >([]);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    
    const ctx = canvas.getContext('2d');
    const width = canvas.width = canvas.offsetWidth * window.devicePixelRatio;
    const height = canvas.height = canvas.offsetHeight * window.devicePixelRatio;
    ctx.scale(window.devicePixelRatio, window.devicePixelRatio);
    
    const displayWidth = canvas.offsetWidth;
    const displayHeight = canvas.offsetHeight;
    // Reserve space on the right when the panel is open so the viz doesn't hide behind it
    const rightPanelReserved = panelOpen ? 340 : 16; // px

    // Initialize particles (subsystems)
      if (particles.current.length === 0) {
      for (let i = 0; i < params.N; i++) {
        const angle = (i / params.N) * Math.PI * 2;
          const base = Math.min(displayWidth - rightPanelReserved, displayHeight);
          const radius = Math.min(base * 0.48 * vizScale, base * 0.9);
        particles.current.push({
          x: displayWidth / 2 + Math.cos(angle) * radius,
          y: displayHeight / 2 + Math.sin(angle) * radius,
          vx: 0,
          vy: 0,
          phase: Math.random() * Math.PI * 2,
          energy: 0,
          entangled: false
        });
      }
    }

    const animate = () => {
      if (!isRunning) return;

      ctx.fillStyle = 'rgba(15, 23, 42, 0.2)';
      ctx.fillRect(0, 0, displayWidth, displayHeight);

      const t = timeRef.current;
      const omega = 2 * Math.PI * 0.3;
      const decay = Math.exp(-0.1 * t);

      // Calculate quantum metrics
      const fluctuation = 2.0 * (1 + 0.5 * Math.sin(omega * t)) * decay + 0.3 * Math.random();
      const entropy = Math.min(Math.log(params.N), 
                               1.5 * (1 - Math.exp(-0.5 * t)) * Math.log(params.N) * 
                               (1 + 0.2 * Math.sin(1.5 * omega * t)));
      const correlation = Math.min(2 * Math.log(params.N), 
                                    2 * entropy * (0.8 + 0.2 * Math.cos(omega * t * 0.7)));
      
      const triad = fluctuation * entropy * correlation;

      setMetrics({
        fluctuation: fluctuation,
        entropy: entropy,
        correlation: correlation,
        triad: triad
      });

      // Update particle behavior based on metrics
      particles.current.forEach((p, i) => {
        // Fluctuation creates jitter
        const jitterStrength = fluctuation * 0.5;
        p.vx += (Math.random() - 0.5) * jitterStrength;
        p.vy += (Math.random() - 0.5) * jitterStrength;

        // Damping
        p.vx *= 0.95;
        p.vy *= 0.95;

        p.x += p.vx;
        p.y += p.vy;

        // Keep in bounds with soft constraint
        const centerX = displayWidth / 2;
        const centerY = displayHeight / 2;
        const dx = p.x - centerX;
        const dy = p.y - centerY;
        const dist = Math.sqrt(dx * dx + dy * dy);
        const maxDist = Math.min(displayWidth, displayHeight) * 0.35;
        
        if (dist > maxDist) {
          p.x = centerX + (dx / dist) * maxDist;
          p.y = centerY + (dy / dist) * maxDist;
          p.vx *= -0.5;
          p.vy *= -0.5;
        }

        // Update phase based on entropy
        p.phase += entropy * 0.1;
        p.energy = fluctuation;
        p.entangled = correlation > 1.0;
      });

      // === Informational curvature tensor (2D) ===
      const cx = displayWidth / 2;
      const cy = displayHeight / 2;
      let Ixx = 0, Iyy = 0, Ixy = 0, Ew = 0;
      let vxMean = 0, vyMean = 0;
      particles.current.forEach((p) => {
        vxMean += p.vx;
        vyMean += p.vy;
      });
      vxMean /= Math.max(1, particles.current.length);
      vyMean /= Math.max(1, particles.current.length);
      particles.current.forEach((p) => {
        const dx = p.x - cx;
        const dy = p.y - cy;
        const e = Math.max(0, p.energy);
        Ixx += e * dx * dx;
        Iyy += e * dy * dy;
        Ixy += e * dx * dy;
        Ew += e;
      });
      if (Ew > 0) {
        Ixx /= Ew; Iyy /= Ew; Ixy /= Ew;
      }
      const trI = Math.max(1e-12, Ixx + Iyy);
      Ixx /= trI; Iyy /= trI; Ixy /= trI; // unit-trace normalization
      // Projector weights and scalar triad amplitude
      const Rinfo_xx = triad * (resp.wxx * Ixx);
      const Rinfo_yy = triad * (resp.wyy * Iyy);
      const Rinfo_xy = triad * (resp.wxy * Ixy);

      // === Debye filter (single-pole) per component ===
      const dt = 0.05 * params.speed; // simulation timestep (same as timeStep increment)
      const alpha = Math.exp(-Math.max(0, dt) / Math.max(resp.tauG, 1e-6));
      debyeState.current.y_xx = alpha * debyeState.current.y_xx + (1 - alpha) * Rinfo_xx;
      debyeState.current.y_yy = alpha * debyeState.current.y_yy + (1 - alpha) * Rinfo_yy;
      debyeState.current.y_xy = alpha * debyeState.current.y_xy + (1 - alpha) * Rinfo_xy;
      const R_xx = debyeState.current.y_xx;
      const R_yy = debyeState.current.y_yy;
      const R_xy = debyeState.current.y_xy;

      // === Einstein proxy in 2D and residual ===
      const trR = R_xx + R_yy;
      const Gxx = R_xx - 0.5 * trR;
      const Gyy = R_yy - 0.5 * trR;
      const Gxy = R_xy; // off-diagonal unchanged
      // Matter proxy: rho ~ triad; p ~ velocity variance
      let varV = 0;
      particles.current.forEach((p) => {
        const dvx = p.vx - vxMean;
        const dvy = p.vy - vyMean;
        varV += 0.5 * (dvx * dvx + dvy * dvy);
      });
      varV /= Math.max(1, particles.current.length);
      const rho = triad;
      const pEff = varV;
      const eightPi = 8.0 * Math.PI;
      const Txx = rho;
      const Tyy = pEff;
      const Txy = 0.0;
      const Exx = Gxx - eightPi * resp.G0 * Txx;
      const Eyy = Gyy - eightPi * resp.G0 * Tyy;
      const Exy = Gxy - eightPi * resp.G0 * Txy;
      const normE = Math.sqrt(Exx * Exx + Eyy * Eyy + 2 * Exy * Exy);
      const normG = Math.max(1e-12, Math.sqrt(Gxx * Gxx + Gyy * Gyy + 2 * Gxy * Gxy));
      const nrmse = normE / normG;
      // EMA of residual
      const emaAlpha = 0.02;
      const newEMA = (1 - emaAlpha) * grav.residualEMA + emaAlpha * nrmse;

      // === Update gravity scales from tauG ===
      const C_LIGHT = 299792458.0; // m/s
      const HBAR = 1.054571817e-34; // J*s
      const omegaG = 1.0 / Math.max(resp.tauG, 1e-12);
      const lG = C_LIGHT * resp.tauG;
      const mG = HBAR / (C_LIGHT * C_LIGHT * Math.max(resp.tauG, 1e-12));
      const slowBandOK = (dt / Math.max(resp.tauG, 1e-6)) <= 0.05;
      setGrav({
        omegaG,
        lG,
        mG,
        residual: nrmse,
        residualEMA: newEMA,
        slowBandOK,
      });

      // === Draw tensor glyph for R ===
      // Eigen decomposition for 2x2 symmetric [R_xx R_xy; R_xy R_yy]
      const a = R_xx, b = R_xy, c = R_yy;
      const tr = a + c;
      const det = a * c - b * b;
      const disc = Math.max(0, tr * tr / 4 - det);
      const lam1 = tr / 2 + Math.sqrt(disc);
      const lam2 = tr / 2 - Math.sqrt(disc);
      // Principal axis angle
      const angle = 0.5 * Math.atan2(2 * b, a - c);
      const base = Math.min(displayWidth - rightPanelReserved, displayHeight);
      const scale = Math.min(base * 0.48 * vizScale, base * 0.95);
      const r1 = Math.max(0, lam1);
      const r2 = Math.max(0, lam2);
      ctx.save();
      ctx.translate(cx, cy);
      ctx.rotate(angle);
      ctx.strokeStyle = 'rgba(255, 255, 255, 0.8)';
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.ellipse(0, 0, Math.sqrt(1 + r1) * (scale * 0.9), Math.sqrt(1 + r2) * (scale * 0.9), 0, 0, Math.PI * 2);
      ctx.stroke();
      ctx.restore();

      // Residual heat ring (smaller baseline radius, thinner ring, amplified alpha)
      const residualGain = 4.0; // visual amplification of residual
      const ringAlpha = Math.min(0.9, residualGain * nrmse);
      const ringRadius = scale * (0.65 + 0.10 * Math.tanh(6 * nrmse)); // smaller base + pulsation
      const gradRing = ctx.createRadialGradient(cx, cy, ringRadius * 0.95, cx, cy, ringRadius); // thin ring
      gradRing.addColorStop(0, `rgba(244, 63, 94, ${0.0})`);
      gradRing.addColorStop(1, `rgba(244, 63, 94, ${ringAlpha})`);
      ctx.fillStyle = gradRing;
      ctx.beginPath();
      ctx.arc(cx, cy, ringRadius, 0, Math.PI * 2);
      ctx.fill();

      // Push to export buffer (cap length)
      exportBuf.current.push({
        t,
        Rinfo_xx,
        Rinfo_yy,
        Rinfo_xy,
        R_xx,
        R_yy,
        R_xy,
        triad,
      });
      if (exportBuf.current.length > 5000) {
        exportBuf.current.shift();
      }

      // Draw entanglement connections (correlation)
      connections.current = [];
      if (correlation > 0.5) {
        for (let i = 0; i < particles.current.length; i++) {
          for (let j = i + 1; j < particles.current.length; j++) {
            const p1 = particles.current[i];
            const p2 = particles.current[j];
            const dx = p2.x - p1.x;
            const dy = p2.y - p1.y;
            const dist = Math.sqrt(dx * dx + dy * dy);
            
            // Probability of entanglement decreases with distance
            const entangleProb = correlation / (2 * Math.log(params.N)) * Math.exp(-dist / 100);
            
            if (Math.random() < entangleProb * 0.3) {
              connections.current.push({ i, j, strength: entangleProb });
              
              // Draw connection
              const alpha = entangleProb * 0.6;
              const gradient = ctx.createLinearGradient(p1.x, p1.y, p2.x, p2.y);
              gradient.addColorStop(0, `rgba(139, 92, 246, ${alpha})`);
              gradient.addColorStop(0.5, `rgba(59, 130, 246, ${alpha})`);
              gradient.addColorStop(1, `rgba(139, 92, 246, ${alpha})`);
              
              ctx.strokeStyle = gradient;
              ctx.lineWidth = 1 + entangleProb * 2;
              ctx.beginPath();
              ctx.moveTo(p1.x, p1.y);
              
              // Wavy connection based on phase
              const midX = (p1.x + p2.x) / 2;
              const midY = (p1.y + p2.y) / 2;
              const wave = Math.sin(t * 2 + i + j) * 10 * correlation;
              const perpX = -(p2.y - p1.y) / dist;
              const perpY = (p2.x - p1.x) / dist;
              
              ctx.quadraticCurveTo(
                midX + perpX * wave,
                midY + perpY * wave,
                p2.x, p2.y
              );
              ctx.stroke();
            }
          }
        }
      }

      // Draw particles (subsystems)
      particles.current.forEach((p, i) => {
        const baseRadius = 8;
        const fluctuationRadius = baseRadius + fluctuation * 2;
        
        // Outer glow (fluctuation)
        const glowRadius = fluctuationRadius * 2;
        const gradient = ctx.createRadialGradient(p.x, p.y, 0, p.x, p.y, glowRadius);
        gradient.addColorStop(0, `rgba(59, 130, 246, ${fluctuation * 0.3})`);
        gradient.addColorStop(1, 'rgba(59, 130, 246, 0)');
        
        ctx.fillStyle = gradient;
        ctx.beginPath();
        ctx.arc(p.x, p.y, glowRadius, 0, Math.PI * 2);
        ctx.fill();

        // Middle ring (entropy - mixedness)
        const entropyAlpha = entropy / Math.log(params.N);
        ctx.strokeStyle = `rgba(16, 185, 129, ${entropyAlpha})`;
        ctx.lineWidth = 2 + entropy * 0.5;
        ctx.beginPath();
        ctx.arc(p.x, p.y, fluctuationRadius * 1.5, 0, Math.PI * 2);
        ctx.stroke();

        // Core particle
        const coreGradient = ctx.createRadialGradient(p.x, p.y, 0, p.x, p.y, fluctuationRadius);
        coreGradient.addColorStop(0, p.entangled ? '#a855f7' : '#3b82f6');
        coreGradient.addColorStop(1, p.entangled ? '#6366f1' : '#1e40af');
        
        ctx.fillStyle = coreGradient;
        ctx.beginPath();
        ctx.arc(p.x, p.y, fluctuationRadius, 0, Math.PI * 2);
        ctx.fill();

        // Phase indicator (rotating dot)
        const phaseX = p.x + Math.cos(p.phase) * fluctuationRadius * 0.6;
        const phaseY = p.y + Math.sin(p.phase) * fluctuationRadius * 0.6;
        ctx.fillStyle = '#fbbf24';
        ctx.beginPath();
        ctx.arc(phaseX, phaseY, 2, 0, Math.PI * 2);
        ctx.fill();
      });

      // Draw central triad indicator
      const centerX = displayWidth / 2;
      const centerY = displayHeight / 2;
      const triadSize = Math.min(50, triad * 5);
      
      const triadGradient = ctx.createRadialGradient(centerX, centerY, 0, centerX, centerY, triadSize);
      triadGradient.addColorStop(0, `rgba(251, 191, 36, ${Math.min(1, triad * 0.1)})`);
      triadGradient.addColorStop(1, 'rgba(251, 191, 36, 0)');
      
      ctx.fillStyle = triadGradient;
      ctx.beginPath();
      ctx.arc(centerX, centerY, triadSize, 0, Math.PI * 2);
      ctx.fill();

      timeRef.current = t + 0.05 * params.speed;
    };

    let animationId;
    const loop = () => {
      animate();
      animationId = requestAnimationFrame(loop);
    };

    if (isRunning) {
      loop();
    }

    return () => {
      if (animationId) cancelAnimationFrame(animationId);
    };
  }, [isRunning, params, panelOpen, vizScale]);

  const reset = () => {
    setIsRunning(false);
    timeRef.current = 0;
    particles.current = [];
    connections.current = [];
    setMetrics({ fluctuation: 0, entropy: 0, correlation: 0, triad: 0 });
    debyeState.current = { y_xx: 0, y_yy: 0, y_xy: 0 };
    exportBuf.current = [];
    setGrav({ omegaG: 0, lG: 0, mG: 0, residual: 0, residualEMA: 0, slowBandOK: true });
  };

  const getIntensityColor = (value, max) => {
    const intensity = Math.min(1, value / max);
    if (intensity < 0.3) return 'bg-red-500';
    if (intensity < 0.7) return 'bg-yellow-500';
    return 'bg-green-500';
  };

  const exportCSV = () => {
    const rows = exportBuf.current;
    if (!rows.length) return;
    const header = ['t', 'R_info_xx', 'R_info_yy', 'R_info_xy', 'R_xx', 'R_yy', 'R_xy', 'triad'].join(',');
    const body = rows.map(r =>
      [r.t, r.Rinfo_xx, r.Rinfo_yy, r.Rinfo_xy, r.R_xx, r.R_yy, r.R_xy, r.triad].join(',')
    ).join('\n');
    const csv = header + '\n' + body;
    const blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'tensor_transfer_visualizer.csv';
    a.click();
    URL.revokeObjectURL(url);
  };

  return (
    <div className="w-full h-screen bg-gradient-to-br from-slate-950 via-slate-900 to-slate-950 text-white flex flex-col">
      {/* Header hidden in fullscreen */}
      {!fullscreen && (
        <div className="p-4 border-b border-slate-800 bg-slate-900/50 backdrop-blur">
          <h1 className="text-2xl font-bold bg-gradient-to-r from-blue-400 via-purple-400 to-yellow-400 bg-clip-text text-transparent">
            Quantum Structural Triad - Visual Representation
          </h1>
          <p className="text-sm text-slate-400 mt-1">
            Observing the interdependence of fluctuation, entropy, and correlation
          </p>
        </div>
      )}

      <div className={fullscreen ? "fixed inset-0" : "flex-1 flex"}>
        {/* Main Canvas */}
        <div className={fullscreen ? "absolute inset-0" : "flex-1 relative"} style={fullscreen ? { position: 'fixed', inset: 0 } : undefined}>
          <canvas 
            ref={canvasRef}
            className="w-full h-full"
            style={{ width: '100vw', height: '100vh', display: 'block' }}
          />
          
          {/* Compact Controls Overlay */}
          <div className="absolute bottom-4 left-1/2 -translate-x-1/2 bg-slate-900/90 backdrop-blur rounded-lg px-3 py-2 border border-slate-700 flex items-center gap-2 shadow-xl" style={{ zIndex: 100 }}>
            <button
              onClick={() => setIsRunning(!isRunning)}
              className="bg-blue-600 hover:bg-blue-700 px-3 py-2 rounded-md text-sm font-semibold"
            >
              <div className="flex items-center gap-2">{isRunning ? <Pause size={16}/> : <Play size={16}/>} {isRunning ? 'Pause' : 'Start'}</div>
            </button>
            <button
              onClick={fitToScreen}
              className="bg-emerald-600 hover:bg-emerald-700 px-3 py-2 rounded-md text-sm font-semibold"
              title="Fit viz to screen"
            >
              Fit
            </button>
            <button
              onClick={reset}
              className="bg-slate-700 hover:bg-slate-600 px-3 py-2 rounded-md text-sm"
            >
              <div className="flex items-center gap-2"><RotateCcw size={16}/> Reset</div>
            </button>
            <button
              onClick={() => setFullscreen(!fullscreen)}
              className="bg-amber-600 hover:bg-amber-700 px-3 py-2 rounded-md text-sm font-semibold"
              title="Toggle fullscreen layout"
            >
              <div className="flex items-center gap-2">{fullscreen ? <Minimize2 size={16}/> : <Maximize2 size={16}/>} {fullscreen ? 'Window' : 'Full'}</div>
            </button>
            <button
              onClick={() => setPanelOpen(!panelOpen)}
              className="bg-slate-600 hover:bg-slate-500 px-3 py-2 rounded-md text-sm font-semibold"
              title="Toggle info panel"
            >
              {panelOpen ? 'Hide info' : 'Show info'}
            </button>
          </div>
          
          {/* Legend Overlay (always visible) */}
          <div className="absolute top-4 left-4 bg-slate-900/95 backdrop-blur rounded-lg p-3 border border-slate-700 shadow-xl" style={{ zIndex: 100 }}>
            <div className="text-xs font-semibold text-slate-300 mb-2 border-b border-slate-600 pb-1">Visual Guide</div>
            <div className="text-xs space-y-1.5">
              <div className="flex items-center gap-2">
                <div className="w-3 h-3 rounded-full bg-blue-500"></div>
                <span>Subsystems (particles)</span>
              </div>
              <div className="flex items-center gap-2">
                <div className="w-3 h-0.5 bg-gradient-to-r from-purple-500 to-blue-500"></div>
                <span>Entanglement links</span>
              </div>
              <div className="flex items-center gap-2">
                <div className="w-3 h-3 rounded-full border-2 border-green-500"></div>
                <span>Entropy ring</span>
              </div>
              <div className="flex items-center gap-2">
                <div className="w-4 h-2 border border-white/80 rounded-sm"></div>
                <span>Response tensor</span>
              </div>
              <div className="flex items-center gap-2">
                <div className="w-3 h-3 rounded-full border-2 border-red-400/60"></div>
                <span>Einstein residual</span>
              </div>
              <div className="flex items-center gap-2">
                <div className="w-2 h-2 rounded-full bg-yellow-400"></div>
                <span>Phase (rotating dot)</span>
              </div>
            </div>
          </div>

          {/* Triad Value Display - only when panel is closed */}
          {!panelOpen && (
            <div className="absolute top-4 right-4 bg-slate-900/90 backdrop-blur rounded-lg p-4 border border-yellow-600/50 shadow-xl" style={{ zIndex: 100 }}>
              <div className="text-center">
                <div className="text-xs text-slate-400 mb-1">Triad Product Π</div>
                <div className="text-3xl font-bold text-yellow-400 font-mono">
                  {metrics.triad.toFixed(2)}
                </div>
                <div className="text-xs text-slate-400 mt-1">
                  {metrics.triad > 1.0 ? '✓ Viable' : '⋯ Evolving'}
                </div>
              </div>
            </div>
          )}

          {/* Info/Controls Panel (overlay, pinned inside the viewport) */}
          {panelOpen && (
            <div className="absolute top-4 right-4 w-80 bg-slate-900/95 border border-slate-700 rounded-lg p-4 space-y-3 overflow-y-auto max-h-[88vh] shadow-xl" style={{ zIndex: 100 }}>
              
              {/* Basic Controls */}
              <div className="border-b border-slate-700 pb-3">
                <h3 className="text-sm font-semibold text-slate-300 mb-2">Display Settings</h3>
                <div className="space-y-2">
                  <div>
                    <label className="block text-xs text-slate-400 mb-1">Number of subsystems N: <span className="font-mono">{params.N}</span></label>
                    <input
                      type="range"
                      min="3"
                      max="24"
                      step="1"
                      value={params.N}
                      onChange={(e) => setNumElements(parseInt(e.target.value))}
                      className="w-full"
                    />
                  </div>
                  <div>
                    <label className="block text-xs text-slate-400 mb-1">Visualization size: {vizScale.toFixed(1)}×</label>
                    <input type="range" min="0.8" max="2.5" step="0.1" value={vizScale} onChange={(e) => setVizScale(parseFloat(e.target.value))} className="w-full" />
                  </div>
                </div>
              </div>

              {/* Response Parameters */}
              <div className="border-b border-slate-700 pb-3">
                <h3 className="text-sm font-semibold text-slate-300 mb-2">Response Parameters</h3>
                <div className="space-y-2">
                  <div>
                    <label className="block text-xs text-slate-400 mb-1">Response time τ_G: <span className="font-mono">{resp.tauG.toFixed(2)} s</span></label>
                    <input type="range" min="0.1" max="10" step="0.1" value={resp.tauG} onChange={(e) => setResp({...resp, tauG: parseFloat(e.target.value)})} className="w-full" />
                  </div>
                  <div>
                    <label className="block text-xs text-slate-400 mb-1">Coupling strength G₀: <span className="font-mono">{resp.G0.toFixed(2)}</span></label>
                    <input type="range" min="0.1" max="5" step="0.1" value={resp.G0} onChange={(e) => setResp({...resp, G0: parseFloat(e.target.value)})} className="w-full" />
                  </div>
                </div>
              </div>

              {/* Tensor Weights */}
              <div className="border-b border-slate-700 pb-3">
                <h3 className="text-sm font-semibold text-slate-300 mb-2">Tensor Projector Weights</h3>
                <div className="grid grid-cols-3 gap-2">
                  <div>
                    <label className="block text-xs text-slate-400 mb-1">w_xx</label>
                    <input type="range" min="0" max="2" step="0.1" value={resp.wxx} onChange={(e) => setResp({...resp, wxx: parseFloat(e.target.value)})} className="w-full" />
                    <div className="text-xs text-center font-mono text-slate-300">{resp.wxx.toFixed(1)}</div>
                  </div>
                  <div>
                    <label className="block text-xs text-slate-400 mb-1">w_yy</label>
                    <input type="range" min="0" max="2" step="0.1" value={resp.wyy} onChange={(e) => setResp({...resp, wyy: parseFloat(e.target.value)})} className="w-full" />
                    <div className="text-xs text-center font-mono text-slate-300">{resp.wyy.toFixed(1)}</div>
                  </div>
                  <div>
                    <label className="block text-xs text-slate-400 mb-1">w_xy</label>
                    <input type="range" min="0" max="2" step="0.1" value={resp.wxy} onChange={(e) => setResp({...resp, wxy: parseFloat(e.target.value)})} className="w-full" />
                    <div className="text-xs text-center font-mono text-slate-300">{resp.wxy.toFixed(1)}</div>
                  </div>
                </div>
              </div>

              {/* Live Metrics */}
              <div className="border-b border-slate-700 pb-3">
                <h3 className="text-sm font-semibold text-slate-300 mb-2">Live Triad Metrics</h3>
                <div className="text-xs space-y-1 bg-slate-800/50 rounded p-2">
                  <div className="flex justify-between">
                    <span className="text-slate-400">Fluctuation √F_Q:</span>
                    <span className="font-mono text-blue-400">{metrics.fluctuation.toFixed(2)}</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-slate-400">Entropy S_A:</span>
                    <span className="font-mono text-green-400">{metrics.entropy.toFixed(2)}</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-slate-400">Correlation I(A:Ā):</span>
                    <span className="font-mono text-purple-400">{metrics.correlation.toFixed(2)}</span>
                  </div>
                  <div className="flex justify-between border-t border-slate-700 pt-1 mt-1">
                    <span className="text-slate-300 font-semibold">Product Π:</span>
                    <span className="font-mono text-yellow-400 font-semibold">{metrics.triad.toFixed(2)}</span>
                  </div>
                </div>
              </div>

              {/* Gravity-Scale Metrics */}
              <div className="pb-2">
                <h3 className="text-sm font-semibold text-slate-300 mb-2">Derived Scales</h3>
                <div className="text-xs space-y-1 text-slate-400">
                  <div className="flex justify-between">
                    <span>Frequency ω_G:</span>
                    <span className="font-mono">{grav.omegaG.toFixed(3)} s⁻¹</span>
                  </div>
                  <div className="flex justify-between">
                    <span>Length ℓ_G:</span>
                    <span className="font-mono">{grav.lG.toExponential(2)} m</span>
                  </div>
                  <div className="flex justify-between">
                    <span>Mass m_G:</span>
                    <span className="font-mono">{grav.mG.toExponential(2)} kg</span>
                  </div>
                  <div className="flex justify-between">
                    <span>Einstein residual:</span>
                    <span className="font-mono">{grav.residual.toFixed(3)}</span>
                  </div>
                  <div className="flex justify-between items-center mt-2 pt-2 border-t border-slate-700">
                    <span>Slow-band check:</span>
                    <span className={grav.slowBandOK ? 'text-green-400 font-semibold' : 'text-red-400 font-semibold'}>
                      {grav.slowBandOK ? '✓ OK' : '⚠ Fast'}
                    </span>
                  </div>
                </div>
              </div>

              <button onClick={exportCSV} className="w-full bg-yellow-600 hover:bg-yellow-700 px-3 py-2 rounded-md text-sm font-semibold transition-colors">
                Export Tensor Time Series (CSV)
              </button>
            </div>
          )}
        </div>

        {/* Right control panel removed in fullscreen mode */}
      </div>
    </div>
  );
};

export default QuantumTriadVisualizer;