import { useEffect, useRef } from "react";

export function BackgroundCanvas() {
  const canvasRef = useRef<HTMLCanvasElement | null>(null);

  useEffect(() => {
    const canvasEl = canvasRef.current;
    if (!canvasEl) return;
    const ctx = canvasEl.getContext("2d");
    if (!ctx) return;
    const canvas: HTMLCanvasElement = canvasEl;
    const context: CanvasRenderingContext2D = ctx;

    const reduced = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
    let width = 0;
    let height = 0;
    let frame = 0;
    let animation = 0;
    let particles: Array<{ x: number; y: number; r: number; speed: number; alpha: number }> = [];

    function resize() {
      const ratio = Math.min(window.devicePixelRatio || 1, 2);
      width = window.innerWidth;
      height = window.innerHeight;
      canvas.width = Math.floor(width * ratio);
      canvas.height = Math.floor(height * ratio);
      canvas.style.width = `${width}px`;
      canvas.style.height = `${height}px`;
      context.setTransform(ratio, 0, 0, ratio, 0, 0);
      particles = Array.from({ length: width < 760 ? 54 : 120 }, () => ({
        x: Math.random() * width,
        y: Math.random() * height,
        r: Math.random() * 1.5 + 0.3,
        speed: Math.random() * 0.12 + 0.03,
        alpha: Math.random() * 0.4 + 0.18,
      }));
    }

    function draw() {
      frame += 1;
      context.clearRect(0, 0, width, height);

      const glow = context.createRadialGradient(width * 0.68, height * 0.16, 0, width * 0.68, height * 0.16, width * 0.66);
      glow.addColorStop(0, "rgba(110, 231, 255, 0.22)");
      glow.addColorStop(0.38, "rgba(93, 255, 189, 0.09)");
      glow.addColorStop(1, "rgba(1, 4, 12, 0)");
      context.fillStyle = glow;
      context.fillRect(0, 0, width, height);

      for (const particle of particles) {
        if (!reduced) {
          particle.x += particle.speed;
          particle.y += Math.sin(frame * 0.006 + particle.x * 0.01) * 0.035;
          if (particle.x > width + 8) particle.x = -8;
        }
        context.globalAlpha = particle.alpha;
        context.fillStyle = "#f7fbff";
        context.beginPath();
        context.arc(particle.x, particle.y, particle.r, 0, Math.PI * 2);
        context.fill();
      }

      context.globalAlpha = 0.35;
      context.lineWidth = 1;
      for (let band = 0; band < 6; band += 1) {
        const offset = band * 72 + Math.sin(frame * 0.01 + band) * 18;
        const gradient = context.createLinearGradient(0, offset, width, offset + 160);
        gradient.addColorStop(0, "rgba(59, 130, 246, 0)");
        gradient.addColorStop(0.5, band % 2 ? "rgba(94, 234, 212, 0.22)" : "rgba(168, 85, 247, 0.18)");
        gradient.addColorStop(1, "rgba(59, 130, 246, 0)");
        context.strokeStyle = gradient;
        context.beginPath();
        for (let x = -80; x <= width + 80; x += 40) {
          const y = height * 0.2 + offset + Math.sin(x * 0.008 + frame * 0.008 + band) * 34;
          if (x === -80) context.moveTo(x, y);
          else context.lineTo(x, y);
        }
        context.stroke();
      }

      context.globalAlpha = 1;
      if (!reduced) animation = requestAnimationFrame(draw);
    }

    resize();
    draw();
    window.addEventListener("resize", resize);
    return () => {
      window.removeEventListener("resize", resize);
      cancelAnimationFrame(animation);
    };
  }, []);

  return <canvas className="background-canvas" ref={canvasRef} aria-hidden="true" />;
}
