import { PALETTE_W } from "../constants";
import { drawBezierLink } from "../helpers";

export default class Chart {
  #idCounter;
  #blocks;
  #links;

  constructor({ blocks = [], links = [], lastId = -1 } = {}) {
    this.#blocks = blocks;
    this.#links = links;
    this.#idCounter = lastId;
  }

  get blocks() {
    return this.#blocks;
  }

  set blocks(bs) {
    this.#blocks = bs;
  }

  get links() {
    return this.#links;
  }

  set links(lks) {
    this.#links = lks;
  }

  addBlock(block) {
    this.#idCounter += 1;
    block.id = `C${this.#idCounter}`;
    this.#blocks.push(block);
  }

  addLink(link) {
    this.#links.push(link);
  }

  mouseIsOver(x) {
    return x > PALETTE_W;
  }

  #drawLink(ctx, link) {
    const getPoint = ({ blockId, connIdx }) => {
      const block = this.#blocks.find((b) => b.id === blockId);
      return block.getConnectionCenter(connIdx);
    };

    const { x0, y0 } = getPoint(link.start);
    const { x0: x1, y0: y1 } = getPoint(link.end);

    drawBezierLink(ctx, x0, y0, x1, y1, link.start.connType === "output");
  }

  draw({ cnv, ctx }) {
    const x0 = PALETTE_W;
    const y0 = 0;
    const h = cnv.clientHeight;
    const w = cnv.clientWidth - x0;
    ctx.fillStyle = "#ffffe0";
    ctx.fillRect(x0, y0, w, h);

    this.#blocks.forEach((b) => b.draw(ctx));
    this.#links.forEach((lnk) => this.#drawLink(ctx, lnk));
  }
}
