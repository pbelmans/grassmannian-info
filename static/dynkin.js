/**
 * Client-side Dynkin diagram renderer.
 *
 * Replicates the SVG layout from the Jinja2 `dynkin()` macro in macros.html.
 * Supports interactive node toggling for the "build your own G/P" tool.
 *
 * Usage:
 *   renderDynkin(container, "A", 5, [2, 4], { interactive: true, onToggle: fn })
 */

const SVG_NS = "http://www.w3.org/2000/svg";

function svgEl(tag, attrs = {}) {
  const el = document.createElementNS(SVG_NS, tag);
  for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v);
  return el;
}

/**
 * Render a Dynkin diagram into the given container element.
 *
 * @param {HTMLElement} container - DOM element to render into (will be cleared)
 * @param {string} type - Dynkin type letter: A, B, C, D, E, F, G
 * @param {number} rank - Rank of the Dynkin diagram
 * @param {number[]} marked - 1-based indices of marked (active) nodes
 * @param {object} [options]
 * @param {boolean} [options.interactive=false] - Whether nodes are clickable
 * @param {function} [options.onToggle] - Callback(markedArray) when a node is toggled
 * @param {number} [options.height=40] - SVG height
 * @param {boolean} [options.showLabels=false] - Show node number labels
 */
export function renderDynkin(container, type, rank, marked = [], options = {}) {
  const { interactive = false, onToggle = null, height = 40, showLabels = false } = options;
  const markedSet = new Set(marked);

  // Normalise B1/C1 → A1
  if ((type === "B" || type === "C") && rank === 1) {
    type = "A";
    rank = 1;
  }

  container.innerHTML = "";

  // Compute viewBox width
  let vbWidth;
  if (type === "D") vbWidth = rank * 20 - 20;
  else if (type === "E") vbWidth = rank * 20 - 10;
  else vbWidth = rank * 20;

  const vbHeight = 40;
  const svg = svgEl("svg", {
    class: "dynkin",
    height: String(height),
    viewBox: `-10 -10 ${vbWidth} ${vbHeight}`,
  });

  const g = svgEl("g", {
    style: "stroke-width:1; stroke-linejoin:bevel; stroke-linecap:butt; stroke:black; fill:white",
  });
  svg.appendChild(g);

  const circles = []; // { element, node (1-based) }

  function addPolyline(points) {
    g.appendChild(svgEl("polyline", { points }));
  }

  function addCircle(cx, cy, node) {
    const circle = svgEl("circle", {
      cx: String(cx),
      cy: String(cy),
      r: "5",
      "data-label": `${type}${rank}-${node}`,
    });
    if (markedSet.has(node)) circle.classList.add("active");
    g.appendChild(circle);
    circles.push({ element: circle, node });
    return circle;
  }

  function addLabel(cx, cy, node) {
    if (!showLabels) return;
    const text = svgEl("text", {
      x: String(cx),
      y: String(cy + 15),
      "text-anchor": "middle",
      "font-size": "8",
      fill: "black",
      stroke: "none",
    });
    text.textContent = String(node);
    g.appendChild(text);
  }

  // --- Type-specific rendering ---

  if (type === "A") {
    if (rank >= 2) addPolyline(`0,10 ${(rank - 1) * 20},10`);
    for (let i = 0; i < rank; i++) {
      addCircle(i * 20, 10, i + 1);
      addLabel(i * 20, 10, i + 1);
    }
  } else if (type === "B" && rank >= 2) {
    // Arrow pointing right: >
    const ax = (rank - 2) * 20;
    addPolyline(`${ax + 7},17 ${ax + 13},10 ${ax + 7},3`);
    // Single bond for first n-2 nodes
    if (rank >= 3) addPolyline(`0,10 ${(rank - 2) * 20},10`);
    // Double bond for last two
    addPolyline(`${(rank - 2) * 20},8 ${(rank - 1) * 20},8`);
    addPolyline(`${(rank - 2) * 20},12 ${(rank - 1) * 20},12`);
    for (let i = 0; i < rank; i++) {
      addCircle(i * 20, 10, i + 1);
      addLabel(i * 20, 10, i + 1);
    }
  } else if (type === "C") {
    // Arrow pointing left: <
    const ax = (rank - 2) * 20;
    addPolyline(`${ax + 13},17 ${ax + 7},10 ${ax + 13},3`);
    // Single bond for first n-2 nodes
    if (rank >= 3) addPolyline(`0,10 ${(rank - 2) * 20},10`);
    // Double bond for last two
    addPolyline(`${(rank - 2) * 20},8 ${(rank - 1) * 20},8`);
    addPolyline(`${(rank - 2) * 20},12 ${(rank - 1) * 20},12`);
    for (let i = 0; i < rank; i++) {
      addCircle(i * 20, 10, i + 1);
      addLabel(i * 20, 10, i + 1);
    }
  } else if (type === "D") {
    // Linear part: nodes 1..(n-2) on the line y=10
    if (rank >= 4) addPolyline(`0,10 ${(rank - 3) * 20},10`);
    // Fork from node (n-2) to (n-1) up and n down
    const forkX = (rank - 3) * 20;
    const tipX = (rank - 2) * 20;
    addPolyline(`${forkX},10 ${tipX},0`);
    addPolyline(`${forkX},10 ${tipX},20`);
    // Nodes on the line
    for (let i = 0; i < rank - 2; i++) {
      addCircle(i * 20, 10, i + 1);
      addLabel(i * 20, 10, i + 1);
    }
    // Forked nodes
    addCircle(tipX, 0, rank - 1);
    addLabel(tipX, -5, rank - 1);
    addCircle(tipX, 20, rank);
    addLabel(tipX, 20, rank);
  } else if (type === "E") {
    // Horizontal chain: nodes 1, 3, 4, 5, ..., n at y=20
    // Node 2 branches down from position x=40 at y=0
    addPolyline(`0,20 ${(rank - 2) * 20},20`);
    addPolyline(`40,20 40,0`);
    // Node 1 at x=0
    addCircle(0, 20, 1);
    addLabel(0, 20, 1);
    // Nodes 3..n along the line
    for (let i = 1; i < rank - 1; i++) {
      addCircle(i * 20, 20, i + 2);
      addLabel(i * 20, 20, i + 2);
    }
    // Node 2 branches up
    addCircle(40, 0, 2);
    addLabel(40, -5, 2);
  } else if (type === "F") {
    // Arrow pointing right between nodes 2 and 3
    addPolyline("27,17 33,10 27,3");
    addPolyline("0,10 20,10");
    addPolyline("20,8 40,8");
    addPolyline("20,12 40,12");
    addPolyline("40,10 60,10");
    for (let i = 0; i < 4; i++) {
      addCircle(i * 20, 10, i + 1);
      addLabel(i * 20, 10, i + 1);
    }
  } else if (type === "G") {
    // Arrow pointing left, triple bond
    addPolyline("13,17 7,10 13,3");
    addPolyline("0,10 20,10");
    addPolyline("0,7 20,7");
    addPolyline("0,13 20,13");
    for (let i = 0; i < 2; i++) {
      addCircle(i * 20, 10, i + 1);
      addLabel(i * 20, 10, i + 1);
    }
  }

  // Interactive click handling
  if (interactive) {
    svg.style.cursor = "pointer";
    for (const { element, node } of circles) {
      element.style.cursor = "pointer";
      element.addEventListener("click", () => {
        if (markedSet.has(node)) {
          markedSet.delete(node);
          element.classList.remove("active");
        } else {
          markedSet.add(node);
          element.classList.add("active");
        }
        if (onToggle) {
          onToggle(Array.from(markedSet).sort((a, b) => a - b));
        }
      });
    }
  }

  container.appendChild(svg);
  return svg;
}
