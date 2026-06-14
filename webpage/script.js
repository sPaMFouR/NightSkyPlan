const canvas = document.querySelector("#starfield");
const context = canvas.getContext("2d");
const prefersReducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

const observatories = {
  HCT: { name: "Hanle, India", lat: 32.78, lon: 78.96, horizon: 25, zenith: 85 },
  DOT: { name: "Devasthal, India", lat: 29.36, lon: 79.68, horizon: 15, zenith: 87.5 },
  KT: { name: "Higashi-Hiroshima, Japan", lat: 34.38, lon: 132.78, horizon: 10, zenith: 90 },
};

const targets = [
  { name: "SN2022jli", ra: 8.688, dec: -8.39, exposure: 45, priority: 3 },
  { name: "SN2020tlf", ra: 220.042, dec: 42.778, exposure: 30, priority: 2 },
  { name: "SN2018zd", ra: 94.513, dec: 78.367, exposure: 20, priority: 1 },
];

let width = 0;
let height = 0;
let stars = [];
let pointer = { x: 0, y: 0 };
let currentPlan = null;

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
    radius: Math.random() * 1.6 + 0.2,
    drift: Math.random() * 0.14 + 0.03,
    alpha: Math.random() * 0.4 + 0.28,
  }));
}

function drawStars() {
  if (prefersReducedMotion) return;

  context.clearRect(0, 0, width, height);
  context.fillStyle = "rgba(255, 255, 255, 0.88)";

  for (const star of stars) {
    star.x += star.drift;
    star.y += pointer.y * 0.00012;
    if (star.x > width + 4) star.x = -4;
    if (star.y > height + 4) star.y = -4;
    if (star.y < -4) star.y = height + 4;

    context.globalAlpha = star.alpha;
    context.beginPath();
    context.arc(star.x, star.y, star.radius, 0, Math.PI * 2);
    context.fill();
  }

  context.globalAlpha = 0.12;
  context.strokeStyle = "#69e1ff";
  context.lineWidth = 1;
  for (let index = 0; index < stars.length - 1; index += 8) {
    const a = stars[index];
    const b = stars[index + 1];
    if (!b) continue;
    const distance = Math.hypot(a.x - b.x, a.y - b.y);
    if (distance < 170) {
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
    { threshold: 0.18 }
  );

  revealItems.forEach(item => observer.observe(item));
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
    heroImage.style.transform = `scale(1.035) translate(${pointer.x * -0.005}px, ${pointer.y * -0.005}px)`;
    panel.style.transform = `translate(${pointer.x * 0.004}px, ${pointer.y * 0.004}px)`;
  });
}

function setupPlanner() {
  renderTargetList();
  updateControlLabels();
  runPlan();

  document.querySelector("#planner-form").addEventListener("submit", event => {
    event.preventDefault();
    runPlan();
  });

  ["observatory-select", "plan-date", "airmass", "moon", "ha"].forEach(id => {
    document.querySelector(`#${id}`).addEventListener("input", () => {
      updateControlLabels();
      runPlan();
    });
  });
}

function updateControlLabels() {
  document.querySelector("#airmass-value").value = Number(document.querySelector("#airmass").value).toFixed(1);
  document.querySelector("#moon-value").value = `${document.querySelector("#moon").value} deg`;
  document.querySelector("#ha-value").value = `${Number(document.querySelector("#ha").value).toFixed(1)} h`;
}

function renderTargetList() {
  const targetList = document.querySelector("#target-list");
  targetList.innerHTML = targets
    .map(
      target => `
        <article class="target-item">
          <div>
            <strong>${target.name}</strong>
            <span>RA ${target.ra.toFixed(2)} deg / Dec ${target.dec.toFixed(2)} deg</span>
          </div>
          <dl>
            <div><dt>Exp</dt><dd>${target.exposure}m</dd></div>
            <div><dt>P</dt><dd>${target.priority}</dd></div>
          </dl>
        </article>
      `
    )
    .join("");
  document.querySelector("#target-count").textContent = `${targets.length} targets`;
}

function runPlan() {
  const obs = observatories[document.querySelector("#observatory-select").value];
  const dateValue = document.querySelector("#plan-date").value;
  const maxAirmass = Number(document.querySelector("#airmass").value);
  const minMoonSep = Number(document.querySelector("#moon").value);
  const haLimit = Number(document.querySelector("#ha").value);
  const moon = approximateMoon(dateValue);
  const timeSlots = Array.from({ length: 25 }, (_, index) => 18 + index * 0.5);

  const tracks = targets.map(target => {
    const samples = timeSlots.map(hour => sampleTarget(target, obs, dateValue, hour, moon));
    const valid = samples.map(
      sample =>
        sample.alt >= obs.horizon &&
        sample.alt <= obs.zenith &&
        sample.airmass <= maxAirmass &&
        Math.abs(sample.ha) <= haLimit &&
        sample.moonSep >= minMoonSep
    );
    const windows = contiguousWindows(samples, valid);
    return { target, samples, valid, windows };
  });

  const schedule = buildSchedule(tracks);
  const meanAirmass =
    schedule.length > 0
      ? schedule.reduce((sum, row) => sum + row.airmass, 0) / schedule.length
      : Math.min(maxAirmass, 2.8);
  const pressure = (3.8 - maxAirmass) * 7 + minMoonSep * 0.16 + Math.max(0, 5 - haLimit) * 7;
  const score = Math.max(48, Math.min(98, Math.round(100 - pressure + schedule.length * 2)));
  const risk = minMoonSep > 74 ? "Low" : minMoonSep > 34 ? "Med" : "High";

  currentPlan = { obs, dateValue, tracks, schedule, score, risk, meanAirmass };
  renderPlan(currentPlan);
  updateEarthSites(obs);
}

function sampleTarget(target, obs, dateValue, hourLocal, moon) {
  const utcHour = (hourLocal - timezoneOffsetHours(obs.lon) + 24) % 24;
  const lst = localSiderealDegrees(dateValue, utcHour, obs.lon);
  let haDeg = normalizeDegrees(lst - target.ra);
  if (haDeg > 180) haDeg -= 360;
  const haRad = degToRad(haDeg);
  const latRad = degToRad(obs.lat);
  const decRad = degToRad(target.dec);
  const sinAlt = Math.sin(decRad) * Math.sin(latRad) + Math.cos(decRad) * Math.cos(latRad) * Math.cos(haRad);
  const alt = radToDeg(Math.asin(sinAlt));
  const airmass = alt > 0 ? 1 / Math.max(0.18, Math.cos(degToRad(90 - alt))) : 9.9;
  const moonSep = angularSeparation(target.ra, target.dec, moon.ra, moon.dec);
  return { hour: hourLocal, alt, airmass, ha: haDeg / 15, moonSep };
}

function buildSchedule(tracks) {
  return tracks
    .map(track => {
      const best = track.samples.reduce(
        (winner, sample, index) => {
          if (!track.valid[index]) return winner;
          const score =
            track.target.priority * 120 +
            sample.alt -
            sample.airmass * 18 +
            sample.moonSep * 0.2 -
            Math.abs(sample.ha) * 4;
          if (!winner || score > winner.score) return { ...sample, score };
          return winner;
        },
        null
      );
      if (!best) return null;
      return {
        target: track.target.name,
        start: formatHour(best.hour),
        end: formatHour(best.hour + track.target.exposure / 60),
        airmass: best.airmass,
        moonSep: best.moonSep,
        priority: track.target.priority,
      };
    })
    .filter(Boolean)
    .sort((a, b) => b.priority - a.priority || Number(a.start.slice(0, 2)) - Number(b.start.slice(0, 2)));
}

function renderPlan(plan) {
  document.querySelector("#quality-pill").textContent = plan.score > 84 ? "Stable" : plan.score > 68 ? "Tight" : "Review";
  document.querySelector("#score").textContent = plan.score;
  document.querySelector("#scheduled").textContent = plan.schedule.length;
  document.querySelector("#risk").textContent = plan.risk;
  document.querySelector("#schedule-summary").textContent =
    plan.schedule.length > 0 ? `${plan.schedule.length} scheduled / ${targets.length} targets` : "No valid windows";
  document.querySelector("#hero-observatory").textContent = plan.obs.name;
  document.querySelector("#hero-date").textContent = plan.dateValue;
  document.querySelector("#hero-coordinates").textContent = formatCoordinates(plan.obs);
  document.querySelector("#hero-scheduled").textContent = plan.schedule.length;
  document.querySelector("#hero-airmass").textContent = plan.meanAirmass.toFixed(2);
  document.querySelector("#hero-dark-window").textContent = darkWindowLabel(plan.obs.lat);
  document.querySelector("#site-name").textContent = plan.obs.name;
  document.querySelector("#site-coordinates").textContent = formatCoordinates(plan.obs);
  document.querySelector("#site-horizon").textContent = `${plan.obs.horizon} deg`;

  document.querySelector("#visibility-rows").innerHTML = plan.tracks
    .map(track => {
      const bestWindow = track.windows[0];
      const widthPercent = bestWindow ? Math.max(9, ((bestWindow.end - bestWindow.start) / 12) * 100) : 6;
      const startPercent = bestWindow ? ((bestWindow.start - 18) / 12) * 100 : 2;
      const state = bestWindow ? "" : " chart__bar--blocked";
      return `
        <div class="chart__row">
          <span>${track.target.name}</span>
          <i class="chart__bar${state}" style="--bar-start:${startPercent}%; --bar-width:${widthPercent}%"></i>
        </div>
      `;
    })
    .join("");

  document.querySelector("#schedule-list").innerHTML =
    plan.schedule.length > 0
      ? plan.schedule
          .map(
            row => `
              <article class="schedule-item">
                <strong>${row.target}</strong>
                <span>${row.start} - ${row.end}</span>
                <small>X ${row.airmass.toFixed(2)} / Moon ${Math.round(row.moonSep)} deg</small>
              </article>
            `
          )
          .join("")
      : `<article class="schedule-item schedule-item--empty">No target clears the current constraints.</article>`;
}

function contiguousWindows(samples, valid) {
  const windows = [];
  let start = null;
  valid.forEach((isValid, index) => {
    if (isValid && start === null) start = samples[index].hour;
    if (start !== null && (!isValid || index === valid.length - 1)) {
      const endIndex = isValid && index === valid.length - 1 ? index : index - 1;
      windows.push({ start, end: samples[endIndex].hour + 0.5 });
      start = null;
    }
  });
  return windows.sort((a, b) => b.end - b.start - (a.end - a.start));
}

function setupEarthGlobes() {
  const canvases = [document.querySelector("#earth-globe"), document.querySelector("#earth-globe-large")].filter(Boolean);
  const globeStates = canvases.map(canvasEl => createGlobeState(canvasEl));

  function draw(time = 0) {
    const obs = currentPlan?.obs || observatories.HCT;
    globeStates.forEach(state => drawGlobe(state, obs, time));
    if (!prefersReducedMotion) requestAnimationFrame(draw);
  }

  draw();
  if (prefersReducedMotion) window.addEventListener("resize", () => draw(0));
}

function createGlobeState(canvasEl) {
  const lights = Array.from({ length: 360 }, (_, index) => ({
    lon: (index * 137.508) % 360 - 180,
    lat: Math.sin(index * 0.83) * 64,
    pulse: Math.random() * Math.PI * 2,
    size: Math.random() * 1.8 + 0.35,
  }));
  return { canvas: canvasEl, context: canvasEl.getContext("2d"), lights, size: 0, ratio: 0 };
}

function drawGlobe(state, obs, time) {
  const size = state.canvas.clientWidth || 260;
  const ratio = Math.min(window.devicePixelRatio || 1, 2);
  if (size !== state.size || ratio !== state.ratio) {
    state.size = size;
    state.ratio = ratio;
    state.canvas.width = Math.floor(size * ratio);
    state.canvas.height = Math.floor(size * ratio);
    state.context.setTransform(ratio, 0, 0, ratio, 0, 0);
  }

  const ctx = state.context;
  const center = size / 2;
  const radius = size * 0.39;
  const rotation = prefersReducedMotion ? -1.55 : -1.55 + time * 0.00005;

  ctx.clearRect(0, 0, size, size);

  const outerHalo = ctx.createRadialGradient(center, center, radius * 0.7, center, center, radius * 1.75);
  outerHalo.addColorStop(0, "rgba(105, 225, 255, 0.2)");
  outerHalo.addColorStop(0.62, "rgba(105, 225, 255, 0.08)");
  outerHalo.addColorStop(1, "rgba(105, 225, 255, 0)");
  ctx.fillStyle = outerHalo;
  ctx.beginPath();
  ctx.arc(center, center, radius * 1.75, 0, Math.PI * 2);
  ctx.fill();

  const planet = ctx.createRadialGradient(center - radius * 0.42, center - radius * 0.36, radius * 0.08, center, center, radius);
  planet.addColorStop(0, "#334453");
  planet.addColorStop(0.28, "#162332");
  planet.addColorStop(0.72, "#08111d");
  planet.addColorStop(1, "#01040a");
  ctx.fillStyle = planet;
  ctx.beginPath();
  ctx.arc(center, center, radius, 0, Math.PI * 2);
  ctx.fill();

  ctx.save();
  ctx.beginPath();
  ctx.arc(center, center, radius, 0, Math.PI * 2);
  ctx.clip();

  drawNightLights(ctx, state.lights, center, radius, rotation, time);
  drawContinentHints(ctx, center, radius, rotation);
  drawGrid(ctx, center, radius, rotation);

  const site = projectGlobe(obs.lon, obs.lat, rotation, center, radius);
  if (site.visible) {
    const pulse = 18 + Math.sin(time * 0.004) * 5;
    ctx.strokeStyle = "rgba(153, 244, 199, 0.65)";
    ctx.lineWidth = 1.4;
    ctx.beginPath();
    ctx.arc(site.x, site.y, pulse, 0, Math.PI * 2);
    ctx.stroke();
    ctx.fillStyle = "#99f4c7";
    ctx.shadowColor = "#99f4c7";
    ctx.shadowBlur = 18;
    ctx.beginPath();
    ctx.arc(site.x, site.y, 4.2, 0, Math.PI * 2);
    ctx.fill();
    ctx.shadowBlur = 0;
  }
  ctx.restore();

  ctx.strokeStyle = "rgba(105, 225, 255, 0.42)";
  ctx.lineWidth = 1.2;
  ctx.beginPath();
  ctx.arc(center, center, radius + 1, 0, Math.PI * 2);
  ctx.stroke();
}

function drawNightLights(ctx, lights, center, radius, rotation, time) {
  lights.forEach(light => {
    const projected = projectGlobe(light.lon, light.lat, rotation, center, radius);
    if (!projected.visible) return;
    const alpha = projected.depth * (0.38 + Math.sin(time * 0.002 + light.pulse) * 0.16);
    ctx.fillStyle = `rgba(255, 198, 103, ${alpha})`;
    ctx.shadowColor = "rgba(255, 198, 103, 0.65)";
    ctx.shadowBlur = 5;
    ctx.beginPath();
    ctx.arc(projected.x, projected.y, light.size, 0, Math.PI * 2);
    ctx.fill();
    ctx.shadowBlur = 0;
  });
}

function drawContinentHints(ctx, center, radius, rotation) {
  ctx.strokeStyle = "rgba(105, 225, 255, 0.13)";
  ctx.lineWidth = 1.1;
  const bands = [
    [-70, 20, 80, -10, 130, 22],
    [-20, 45, 32, 18, 90, 8],
    [70, 28, 108, -6, 146, 10],
    [-130, -18, -82, -36, -48, -12],
  ];
  bands.forEach(points => {
    ctx.beginPath();
    for (let i = 0; i < points.length; i += 2) {
      const p = projectGlobe(points[i], points[i + 1], rotation, center, radius);
      if (!p.visible) continue;
      if (i === 0) ctx.moveTo(p.x, p.y);
      else ctx.lineTo(p.x, p.y);
    }
    ctx.stroke();
  });
}

function drawGrid(ctx, center, radius, rotation) {
  ctx.strokeStyle = "rgba(245, 242, 234, 0.08)";
  ctx.lineWidth = 1;
  for (let lat = -45; lat <= 45; lat += 30) {
    ctx.beginPath();
    for (let lon = -180; lon <= 180; lon += 8) {
      const p = projectGlobe(lon, lat, rotation, center, radius);
      if (!p.visible) continue;
      if (lon === -180) ctx.moveTo(p.x, p.y);
      else ctx.lineTo(p.x, p.y);
    }
    ctx.stroke();
  }
  for (let lon = -120; lon <= 120; lon += 60) {
    ctx.beginPath();
    for (let lat = -80; lat <= 80; lat += 5) {
      const p = projectGlobe(lon, lat, rotation, center, radius);
      if (!p.visible) continue;
      if (lat === -80) ctx.moveTo(p.x, p.y);
      else ctx.lineTo(p.x, p.y);
    }
    ctx.stroke();
  }
}

function updateEarthSites(obs) {
  document.querySelector("#hero-coordinates").textContent = formatCoordinates(obs);
}

function projectGlobe(lonDeg, latDeg, rotation, center, radius) {
  const lon = degToRad(lonDeg) + rotation;
  const lat = degToRad(latDeg);
  const depth = Math.cos(lon) * Math.cos(lat);
  return {
    visible: depth > -0.12,
    depth: Math.max(0, depth),
    x: center + Math.sin(lon) * Math.cos(lat) * radius,
    y: center - Math.sin(lat) * radius,
  };
}

function approximateMoon(dateValue) {
  const days = julianDays(dateValue);
  return {
    ra: normalizeDegrees(218.32 + 13.176396 * days),
    dec: 5.1 * Math.sin(degToRad(134.9 + 13.064993 * days)),
  };
}

function localSiderealDegrees(dateValue, utcHour, lon) {
  const days = julianDays(dateValue);
  return normalizeDegrees(100.46 + 0.985647 * days + lon + 15 * utcHour);
}

function julianDays(dateValue) {
  const date = new Date(`${dateValue}T00:00:00Z`);
  return date.getTime() / 86400000 - Date.UTC(2000, 0, 1, 12) / 86400000;
}

function angularSeparation(ra1, dec1, ra2, dec2) {
  const r1 = degToRad(ra1);
  const d1 = degToRad(dec1);
  const r2 = degToRad(ra2);
  const d2 = degToRad(dec2);
  const cosSep = Math.sin(d1) * Math.sin(d2) + Math.cos(d1) * Math.cos(d2) * Math.cos(r1 - r2);
  return radToDeg(Math.acos(Math.max(-1, Math.min(1, cosSep))));
}

function timezoneOffsetHours(lon) {
  return Math.round(lon / 15);
}

function normalizeDegrees(value) {
  return ((value % 360) + 360) % 360;
}

function degToRad(value) {
  return (value * Math.PI) / 180;
}

function radToDeg(value) {
  return (value * 180) / Math.PI;
}

function formatHour(hour) {
  const wrapped = ((hour % 24) + 24) % 24;
  const h = Math.floor(wrapped);
  const m = Math.round((wrapped - h) * 60);
  return `${String(h).padStart(2, "0")}:${String(m).padStart(2, "0")}`;
}

function formatCoordinates(obs) {
  return `${Math.abs(obs.lat).toFixed(2)} ${obs.lat >= 0 ? "N" : "S"}, ${Math.abs(obs.lon).toFixed(2)} ${
    obs.lon >= 0 ? "E" : "W"
  }`;
}

function darkWindowLabel(lat) {
  const hours = Math.max(6.4, Math.min(10.2, 8.5 + Math.abs(lat - 30) * 0.03));
  const h = Math.floor(hours);
  const m = Math.round((hours - h) * 60);
  return `${h}h ${m}m`;
}

resizeCanvas();
setupReveals();
setupHeroParallax();
setupPlanner();
setupEarthGlobes();

if (!prefersReducedMotion) {
  drawStars();
}

window.addEventListener("resize", resizeCanvas);
