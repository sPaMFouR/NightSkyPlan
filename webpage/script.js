const canvas = document.querySelector("#starfield");
const context = canvas.getContext("2d");
const prefersReducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

let width = 0;
let height = 0;
let stars = [];
let pointer = { x: 0, y: 0 };

function resizeCanvas() {
  const ratio = Math.min(window.devicePixelRatio || 1, 2);
  width = window.innerWidth;
  height = window.innerHeight;
  canvas.width = Math.floor(width * ratio);
  canvas.height = Math.floor(height * ratio);
  canvas.style.width = `${width}px`;
  canvas.style.height = `${height}px`;
  context.setTransform(ratio, 0, 0, ratio, 0, 0);

  const starCount = width < 720 ? 70 : 140;
  stars = Array.from({ length: starCount }, () => ({
    x: Math.random() * width,
    y: Math.random() * height,
    radius: Math.random() * 1.8 + 0.25,
    drift: Math.random() * 0.18 + 0.04,
    alpha: Math.random() * 0.45 + 0.32,
  }));
}

function drawStars() {
  if (prefersReducedMotion) return;

  context.clearRect(0, 0, width, height);
  context.fillStyle = "rgba(255, 255, 255, 0.88)";

  for (const star of stars) {
    star.x += star.drift;
    star.y += pointer.y * 0.00016;
    if (star.x > width + 4) star.x = -4;
    if (star.y > height + 4) star.y = -4;
    if (star.y < -4) star.y = height + 4;

    context.globalAlpha = star.alpha;
    context.beginPath();
    context.arc(star.x, star.y, star.radius, 0, Math.PI * 2);
    context.fill();
  }

  context.globalAlpha = 0.16;
  context.strokeStyle = "#69e1ff";
  context.lineWidth = 1;
  for (let index = 0; index < stars.length - 1; index += 7) {
    const a = stars[index];
    const b = stars[index + 1];
    if (!b) continue;
    const distance = Math.hypot(a.x - b.x, a.y - b.y);
    if (distance < 180) {
      context.beginPath();
      context.moveTo(a.x, a.y);
      context.lineTo(b.x, b.y);
      context.stroke();
    }
  }

  context.globalAlpha = 1;
  requestAnimationFrame(drawStars);
}

function setupReveals() {
  const revealItems = document.querySelectorAll(".reveal");
  const observer = new IntersectionObserver(
    entries => {
      for (const entry of entries) {
        if (entry.isIntersecting) {
          entry.target.classList.add("is-visible");
          observer.unobserve(entry.target);
        }
      }
    },
    { threshold: 0.22 }
  );

  revealItems.forEach(item => observer.observe(item));
}

function setupPlannerControls() {
  const airmass = document.querySelector("#airmass");
  const moon = document.querySelector("#moon");
  const ha = document.querySelector("#ha");
  const airmassValue = document.querySelector("#airmass-value");
  const moonValue = document.querySelector("#moon-value");
  const haValue = document.querySelector("#ha-value");
  const score = document.querySelector("#score");
  const scheduled = document.querySelector("#scheduled");
  const risk = document.querySelector("#risk");
  const quality = document.querySelector("#quality-pill");
  const bars = document.querySelectorAll(".chart__bar");

  function update() {
    const airmassNum = Number(airmass.value);
    const moonNum = Number(moon.value);
    const haNum = Number(ha.value);

    airmassValue.value = airmassNum.toFixed(1);
    moonValue.value = `${moonNum} deg`;
    haValue.value = `${haNum.toFixed(1)} h`;

    const constraintPressure =
      (3.8 - airmassNum) * 8 + moonNum * 0.18 + Math.max(0, 5 - haNum) * 7;
    const computedScore = Math.max(52, Math.min(98, Math.round(99 - constraintPressure)));
    const scheduledCount = Math.max(6, Math.min(14, Math.round(computedScore / 8.1)));

    score.textContent = computedScore;
    scheduled.textContent = scheduledCount;
    risk.textContent = moonNum > 74 ? "Low" : moonNum > 34 ? "Med" : "High";
    quality.textContent = computedScore > 82 ? "Stable" : computedScore > 68 ? "Tight" : "Review";

    bars.forEach((bar, index) => {
      const widthValue = Math.max(23, Math.min(82, computedScore - index * 11 + haNum * 2));
      const startValue = Math.max(4, Math.min(36, 39 - widthValue / 3 + index * 8));
      bar.style.setProperty("--bar-width", `${widthValue}%`);
      bar.style.setProperty("--bar-start", `${startValue}%`);
    });
  }

  [airmass, moon, ha].forEach(control => control.addEventListener("input", update));
  update();
}

function setupHeroParallax() {
  const heroImage = document.querySelector(".hero__image");
  const panel = document.querySelector(".mission-panel");
  window.addEventListener("pointermove", event => {
    pointer = {
      x: event.clientX - width / 2,
      y: event.clientY - height / 2,
    };

    if (prefersReducedMotion || width < 900) return;
    const moveX = pointer.x * -0.006;
    const moveY = pointer.y * -0.006;
    heroImage.style.transform = `scale(1.04) translate(${moveX}px, ${moveY}px)`;
    panel.style.transform = `translate(${pointer.x * 0.005}px, ${pointer.y * 0.005}px)`;
  });
}

function setupEarthGlobe() {
  const globe = document.querySelector("#earth-globe");
  if (!globe) return;

  const globeContext = globe.getContext("2d");
  let currentSize = 0;
  let currentRatio = 0;
  const lights = Array.from({ length: 120 }, (_, index) => ({
    lon: (index * 137.508) % 360 - 180,
    lat: Math.sin(index * 1.91) * 58,
    pulse: Math.random() * Math.PI * 2,
    size: Math.random() * 1.4 + 0.5,
  }));

  function configureCanvas() {
    const parentWidth = globe.clientWidth || 260;
    const ratio = Math.min(window.devicePixelRatio || 1, 2);
    if (parentWidth !== currentSize || ratio !== currentRatio) {
      currentSize = parentWidth;
      currentRatio = ratio;
      globe.width = Math.floor(parentWidth * ratio);
      globe.height = Math.floor(parentWidth * ratio);
      globeContext.setTransform(ratio, 0, 0, ratio, 0, 0);
    }
    return parentWidth;
  }

  function draw(time = 0) {
    const size = configureCanvas();
    const center = size / 2;
    const radius = size * 0.42;
    const rotation = prefersReducedMotion ? 0 : time * 0.000035;

    globeContext.clearRect(0, 0, size, size);

    const halo = globeContext.createRadialGradient(center, center, radius * 0.55, center, center, radius * 1.35);
    halo.addColorStop(0, "rgba(105, 225, 255, 0.18)");
    halo.addColorStop(0.58, "rgba(105, 225, 255, 0.08)");
    halo.addColorStop(1, "rgba(105, 225, 255, 0)");
    globeContext.fillStyle = halo;
    globeContext.beginPath();
    globeContext.arc(center, center, radius * 1.36, 0, Math.PI * 2);
    globeContext.fill();

    const body = globeContext.createRadialGradient(center - radius * 0.35, center - radius * 0.35, radius * 0.1, center, center, radius);
    body.addColorStop(0, "#29313a");
    body.addColorStop(0.42, "#121a24");
    body.addColorStop(1, "#04070d");
    globeContext.fillStyle = body;
    globeContext.beginPath();
    globeContext.arc(center, center, radius, 0, Math.PI * 2);
    globeContext.fill();

    globeContext.save();
    globeContext.beginPath();
    globeContext.arc(center, center, radius, 0, Math.PI * 2);
    globeContext.clip();

    globeContext.strokeStyle = "rgba(245, 242, 234, 0.11)";
    globeContext.lineWidth = 1;
    for (let i = -3; i <= 3; i += 1) {
      globeContext.beginPath();
      globeContext.ellipse(center, center, radius * (0.2 + Math.abs(i) * 0.13), radius, rotation + i * 0.04, 0, Math.PI * 2);
      globeContext.stroke();
    }
    for (let i = -2; i <= 2; i += 1) {
      globeContext.beginPath();
      globeContext.ellipse(center, center + i * radius * 0.22, radius, radius * 0.22, 0, 0, Math.PI * 2);
      globeContext.stroke();
    }

    lights.forEach(light => {
      const lon = ((light.lon * Math.PI) / 180) + rotation;
      const lat = (light.lat * Math.PI) / 180;
      const depth = Math.cos(lon);
      if (depth < -0.08) return;
      const x = center + Math.sin(lon) * Math.cos(lat) * radius;
      const y = center - Math.sin(lat) * radius;
      const alpha = Math.max(0, depth) * (0.42 + Math.sin(time * 0.002 + light.pulse) * 0.18);
      globeContext.fillStyle = `rgba(242, 177, 95, ${alpha})`;
      globeContext.beginPath();
      globeContext.arc(x, y, light.size, 0, Math.PI * 2);
      globeContext.fill();
    });

    const hanle = projectGlobe(78.96, 32.78, rotation, center, radius);
    if (hanle.visible) {
      globeContext.strokeStyle = "rgba(153, 244, 199, 0.55)";
      globeContext.lineWidth = 1.5;
      globeContext.beginPath();
      globeContext.arc(hanle.x, hanle.y, 20 + Math.sin(time * 0.004) * 4, 0, Math.PI * 2);
      globeContext.stroke();
      globeContext.fillStyle = "#99f4c7";
      globeContext.beginPath();
      globeContext.arc(hanle.x, hanle.y, 4.8, 0, Math.PI * 2);
      globeContext.fill();
    }

    globeContext.restore();

    if (!prefersReducedMotion) {
      requestAnimationFrame(draw);
    }
  }

  draw();
  if (prefersReducedMotion) {
    window.addEventListener("resize", draw);
  }
}

function projectGlobe(lonDeg, latDeg, rotation, center, radius) {
  const lon = (lonDeg * Math.PI) / 180 + rotation;
  const lat = (latDeg * Math.PI) / 180;
  const depth = Math.cos(lon);
  return {
    visible: depth > -0.08,
    x: center + Math.sin(lon) * Math.cos(lat) * radius,
    y: center - Math.sin(lat) * radius,
  };
}

function setupCopyCommand() {
  const button = document.querySelector("#copy-command");
  const status = document.querySelector("#copy-status");
  const command = "streamlit run application/app.py";
  if (!button || !status) return;

  button.addEventListener("click", async () => {
    try {
      await navigator.clipboard.writeText(command);
      status.textContent = "Command copied.";
    } catch {
      status.textContent = command;
    }
  });
}

resizeCanvas();
setupReveals();
setupPlannerControls();
setupHeroParallax();
setupEarthGlobe();
setupCopyCommand();

if (!prefersReducedMotion) {
  drawStars();
}

window.addEventListener("resize", resizeCanvas);
