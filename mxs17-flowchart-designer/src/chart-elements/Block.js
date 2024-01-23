import {
  BLOCK_CONN_R,
  BLOCK_H,
  BLOCK_TYPE,
  BLOCK_W,
  LINE_COLOR,
  SHAPE_DIAM,
  SHAPE_RECT,
  SHAPE_ROUND,
} from "../constants";

export default class Block {
  #x;
  #y;
  #w;
  #h;
  #value;
  #isIcon = false;
  #inputs;
  #outputs;

  constructor({
    id,
    blockType,
    opType,
    symbol,
    value,
    x,
    y,
    w = BLOCK_W,
    h = BLOCK_H,
    shape = SHAPE_RECT,
    inputs = [],
    outputs = [],
  }) {
    this.id = id;
    this.blockType = blockType;
    this.opType = opType;
    this.symbol = symbol;
    this.#value = value;
    this.shape = shape;
    this.#x = x;
    this.#y = y;
    this.#w = w;
    this.#h = h;
    this.#inputs = inputs;
    this.#outputs = outputs;
  }

  get x() {
    return this.#x;
  }

  set x(value) {
    this.#x = value;
  }

  get y() {
    return this.#y;
  }

  set y(value) {
    this.#y = value;
  }

  get isIcon() {
    return this.#isIcon;
  }

  set isIcon(value) {
    this.#isIcon = value;
  }

  get inputs() {
    return this.#inputs;
  }

  get outputs() {
    return this.#outputs;
  }

  set value(v) {
    if (!v) {
      this.#value = null;
    } else if (
      this.blockType === BLOCK_TYPE.NUM &&
      !isNaN(parseInt(v, 10)) &&
      parseInt(v, 10) < 100
    ) {
      this.#value = parseInt(v, 10);
    } else if (
      this.blockType === BLOCK_TYPE.BOOL &&
      ["T", "F"].includes(v.toUpperCase())
    ) {
      this.#value = v.toUpperCase();
    }
  }

  get value() {
    return this.#value;
  }

  mouseIsOver(x, y) {
    return (
      !this.getConnectionUnderMouse(x, y) &&
      x >= this.#x &&
      x <= this.#x + this.#w &&
      y >= this.#y &&
      y <= this.#y + this.#h
    );
  }

  getConnectionCenter(index) {
    if (!this.#inputs.includes(index) && !this.#outputs.includes(index)) {
      return null;
    }
    const x0 = [0.5 * this.#w, this.#w, 0.5 * this.#w, 0][index] + this.#x;
    const y0 = [0, 0.5 * this.#h, this.#h, 0.5 * this.#h][index] + this.#y;
    return { x0, y0 };
  }

  getConnectionUnderMouse(x, y) {
    const mouseInsideCircle = (idx) => {
      const { x0, y0 } = this.getConnectionCenter(idx);
      return Math.sqrt((x - x0) ** 2 + (y - y0) ** 2) <= BLOCK_CONN_R * 2;
    };

    const inputIdx = this.#inputs.find((idx) => mouseInsideCircle(idx));
    if (inputIdx >= 0)
      return {
        type: "input",
        index: inputIdx,
        ...this.getConnectionCenter(inputIdx),
      };

    const outputIdx = this.#outputs.find((idx) => mouseInsideCircle(idx));
    return outputIdx >= 0
      ? {
          type: "output",
          index: outputIdx,
          ...this.getConnectionCenter(outputIdx),
        }
      : null;
  }

  #drawConnection(ctx, index, isInput) {
    const { x0, y0 } = this.getConnectionCenter(index);
    ctx.beginPath();
    ctx.arc(x0, y0, BLOCK_CONN_R, 0, 2 * Math.PI);
    if (isInput) {
      ctx.fillStyle = LINE_COLOR;
      ctx.fill();
    } else {
      ctx.strokeStyle = LINE_COLOR;
      ctx.stroke();
    }
  }

  draw(ctx, dragX, dragY) {
    const x = dragX || this.#x;
    const y = dragY || this.#y;
    const w = this.#w;
    const h = this.#h;

    ctx.beginPath();
    ctx.strokeStyle = LINE_COLOR;
    if (this.shape === SHAPE_DIAM) {
      ctx.moveTo(x + w / 2, y);
      ctx.lineTo(x + w, y + h / 2);
      ctx.lineTo(x + w / 2, y + h);
      ctx.lineTo(x, y + h / 2);
      ctx.lineTo(x + w / 2, y);
      ctx.stroke();
    } else if (this.shape === SHAPE_ROUND) {
      ctx.arc(x + w / 2, y + h / 2, Math.min(w, h) / 2, 0, 2 * Math.PI);
      ctx.stroke();
    } else {
      ctx.strokeRect(x, y, w, h);
    }
    ctx.font = "16px sans-serif";
    ctx.textAlign = "center";
    ctx.fillStyle = "black";
    ctx.fillText(this.#value || this.symbol, x + w / 2, y + h / 2 + 5);

    this.#inputs.forEach((idx) => this.#drawConnection(ctx, idx, true));
    this.#outputs.forEach((idx) => this.#drawConnection(ctx, idx, false));
  }

  clone({ id }) {
    return new Block({
      id: id || this.id,
      blockType: this.blockType,
      opType: this.opType,
      symbol: this.symbol,
      value: this.value,
      x: this.#x,
      y: this.#y,
      w: this.#w,
      h: this.#h,
      inputs: this.#inputs,
      outputs: this.#outputs,
      shape: this.shape,
    });
  }
}
