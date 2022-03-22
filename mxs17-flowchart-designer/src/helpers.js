import { LINE_COLOR, STORAGE_KEY } from "./constants";

export function drawBezierLink(ctx, x0, y0, x1, y1, fromOutput = true) {
  const cp1 = fromOutput
    ? { x: (x1 - x0) / 2 + x0, y: y0 }
    : { x: x0, y: y1 - (y1 - y0) / 5 };
  const cp2 = fromOutput
    ? { x: x1, y: (y1 - y0) / 5 + y0 }
    : { x: (x1 - x0) / 2 + x0, y: y1 };

  ctx.strokeStyle = LINE_COLOR;
  ctx.beginPath();
  ctx.moveTo(x0, y0);
  ctx.bezierCurveTo(cp1.x, cp1.y, cp2.x, cp2.y, x1, y1);
  ctx.stroke();
  drawArrow(ctx, fromOutput ? x1 : x0, fromOutput ? y1 : y0);
}

function drawArrow(ctx, x, y) {
  ctx.beginPath();
  ctx.moveTo(x, y);
  ctx.lineTo(x - 3, y - 8);
  ctx.lineTo(x + 3, y - 8);
  ctx.lineTo(x, y);
  ctx.stroke();
}

export function serializeBlocks(blocks) {
  const fields = [
    "id",
    "blockType",
    "opType",
    "symbol",
    "value",
    "shape",
    "x",
    "y",
    "w",
    "h",
    "inputs",
    "outputs",
  ];
  return blocks?.map((b) =>
    fields.reduce((acc, f) => {
      acc[f] = b[f];
      return acc;
    }, {})
  );
}

export function serializeLinks(links) {
  return JSON.parse(JSON.stringify(links));
}

export function persistChartData({ blocks, links }) {
  localStorage.setItem(
    STORAGE_KEY,
    JSON.stringify({
      blocks: serializeBlocks(blocks),
      links: serializeLinks(links),
    })
  );
}

export function restoreChartData() {
  const stored = localStorage.getItem(STORAGE_KEY);
  return stored ? JSON.parse(stored) : null;
}

export function resetPersistedData() {
  localStorage.removeItem(STORAGE_KEY);
}
