import { useEffect, useRef } from "react";
import type { Observatory } from "./types";

interface EarthGlobeProps {
  observatory?: Observatory;
}

export function EarthGlobe({ observatory }: EarthGlobeProps) {
  const canvasRef = useRef<HTMLCanvasElement | null>(null);

  useEffect(() => {
    const canvasEl = canvasRef.current;
    const activeObservatory = observatory;
    if (!canvasEl || !activeObservatory) return;
    const ctx = canvasEl.getContext("2d");
    if (!ctx) return;
    const canvas: HTMLCanvasElement = canvasEl;
    const context: CanvasRenderingContext2D = ctx;
    const siteObservatory: Observatory = activeObservatory;

    const reduced = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
    const lights = Array.from({ length: 420 }, (_, index) => ({
      lon: ((index * 137.508) % 360) - 180,
      lat: Math.sin(index * 0.77) * 62,
      size: Math.random() * 1.4 + 0.35,
      pulse: Math.random() * Math.PI * 2,
    }));
    let animation = 0;

    function draw(time = 0) {
      const size = canvas.clientWidth || 320;
      const ratio = Math.min(window.devicePixelRatio || 1, 2);
      canvas.width = Math.floor(size * ratio);
      canvas.height = Math.floor(size * ratio);
      context.setTransform(ratio, 0, 0, ratio, 0, 0);
      context.clearRect(0, 0, size, size);

      const center = size / 2;
      const radius = size * 0.38;
      const rotation = -1.6 + (reduced ? 0 : time * 0.000045);
      const halo = context.createRadialGradient(center, center, radius * 0.7, center, center, radius * 1.8);
      halo.addColorStop(0, "rgba(125, 244, 255, 0.2)");
      halo.addColorStop(0.6, "rgba(85, 255, 187, 0.08)");
      halo.addColorStop(1, "rgba(0, 0, 0, 0)");
      context.fillStyle = halo;
      context.beginPath();
      context.arc(center, center, radius * 1.8, 0, Math.PI * 2);
      context.fill();

      const globe = context.createRadialGradient(center - radius * 0.42, center - radius * 0.38, 0, center, center, radius);
      globe.addColorStop(0, "#364b5f");
      globe.addColorStop(0.45, "#102035");
      globe.addColorStop(1, "#020611");
      context.fillStyle = globe;
      context.beginPath();
      context.arc(center, center, radius, 0, Math.PI * 2);
      context.fill();

      context.save();
      context.beginPath();
      context.arc(center, center, radius, 0, Math.PI * 2);
      context.clip();

      for (const light of lights) {
        const point = project(light.lon, light.lat, rotation, center, radius);
        if (!point.visible) continue;
        const alpha = 0.18 + point.depth * 0.5 + Math.sin(time * 0.003 + light.pulse) * 0.05;
        context.globalAlpha = Math.max(0.08, alpha);
        context.fillStyle = "#ffd98a";
        context.beginPath();
        context.arc(point.x, point.y, light.size * (0.6 + point.depth), 0, Math.PI * 2);
        context.fill();
      }

      context.globalAlpha = 0.18;
      context.strokeStyle = "#a9f7ff";
      context.lineWidth = 1;
      for (let lat = -60; lat <= 60; lat += 30) {
        context.beginPath();
        for (let lon = -180; lon <= 180; lon += 5) {
          const point = project(lon, lat, rotation, center, radius);
          if (!point.visible) continue;
          if (lon === -180) context.moveTo(point.x, point.y);
          else context.lineTo(point.x, point.y);
        }
        context.stroke();
      }
      for (let lon = -150; lon <= 180; lon += 30) {
        context.beginPath();
        let started = false;
        for (let lat = -80; lat <= 80; lat += 4) {
          const point = project(lon, lat, rotation, center, radius);
          if (!point.visible) continue;
          if (!started) {
            context.moveTo(point.x, point.y);
            started = true;
          } else {
            context.lineTo(point.x, point.y);
          }
        }
        context.stroke();
      }

      const site = project(siteObservatory.longitude_deg, siteObservatory.latitude_deg, rotation, center, radius);
      if (site.visible) {
        const pulse = 16 + Math.sin(time * 0.004) * 5;
        context.globalAlpha = 0.8;
        context.strokeStyle = "#78ffc6";
        context.lineWidth = 1.6;
        context.beginPath();
        context.arc(site.x, site.y, pulse, 0, Math.PI * 2);
        context.stroke();
        context.fillStyle = "#edfff7";
        context.beginPath();
        context.arc(site.x, site.y, 4.2, 0, Math.PI * 2);
        context.fill();
      }

      context.restore();
      context.globalAlpha = 1;
      context.strokeStyle = "rgba(157, 242, 255, 0.42)";
      context.lineWidth = 1.3;
      context.beginPath();
      context.arc(center, center, radius, 0, Math.PI * 2);
      context.stroke();

      if (!reduced) animation = requestAnimationFrame(draw);
    }

    const handleResize = () => draw(0);
    draw();
    window.addEventListener("resize", handleResize);
    return () => {
      window.removeEventListener("resize", handleResize);
      cancelAnimationFrame(animation);
    };
  }, [observatory]);

  return <canvas className="earth-globe" ref={canvasRef} aria-label="Illuminated Earth with selected observatory marker" />;
}

function project(lon: number, lat: number, rotation: number, center: number, radius: number) {
  const phi = (lat * Math.PI) / 180;
  const lambda = (lon * Math.PI) / 180 + rotation;
  const x3 = Math.cos(phi) * Math.sin(lambda);
  const y3 = Math.sin(phi);
  const z3 = Math.cos(phi) * Math.cos(lambda);
  return {
    x: center + x3 * radius,
    y: center - y3 * radius,
    visible: z3 > -0.18,
    depth: Math.max(0, z3),
  };
}
