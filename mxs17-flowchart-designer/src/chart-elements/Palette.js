import { PALETTE_W, BLOCK_W, BLOCK_H, PALETTE_GAP } from "../constants";

export default class Palette {
  #icons;

  constructor(icons) {
    this.#icons = icons.map((icon, i) => {
      icon.x = (PALETTE_W - BLOCK_W) / 2;
      icon.y = PALETTE_GAP * 1.5 + i * (BLOCK_H + PALETTE_GAP);
      icon.isIcon = true;
      return icon;
    });
  }

  get icons() {
    return this.#icons;
  }

  mouseIsOver(x) {
    return x <= PALETTE_W;
  }

  draw({ ctx, cnv }) {
    ctx.fillStyle = "#f2f2f2";
    ctx.fillRect(0, 0, PALETTE_W, cnv.clientHeight);
    ctx.globalAlpha = 0.6;
    this.#icons.forEach((icon) => icon.draw(ctx));
    ctx.globalAlpha = 1;
  }
}
