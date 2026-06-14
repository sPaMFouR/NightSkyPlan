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

function setupPrototypeControls() {
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

resizeCanvas();
setupReveals();
setupPrototypeControls();
setupHeroParallax();

if (!prefersReducedMotion) {
  drawStars();
}

window.addEventListener("resize", resizeCanvas);
