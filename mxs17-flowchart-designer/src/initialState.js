import { BLOCK_H, BLOCK_TYPE, BLOCK_W, OP, SHAPE_ROUND } from "./constants";
import Block from "./chart-elements/Block";
import Chart from "./chart-elements/Chart";
import Palette from "./chart-elements/Palette";
import { restoreChartData } from "./helpers";
import Link from "./chart-elements/Link";

const STORED = restoreChartData();

const INITIAL_STATE = {
  palette: new Palette(
    [
      {
        blockType: BLOCK_TYPE.NUM,
        opType: null,
        symbol: "0-9",
        isOperand: true,
      },
      {
        blockType: BLOCK_TYPE.BOOL,
        opType: null,
        symbol: "t / f",
        isOperand: true,
      },
      { blockType: BLOCK_TYPE.OP, opType: OP.PLUS, symbol: OP.PLUS },
      { blockType: BLOCK_TYPE.OP, opType: OP.MINUS, symbol: OP.MINUS },
      { blockType: BLOCK_TYPE.OP, opType: OP.MULT, symbol: OP.MULT },
      { blockType: BLOCK_TYPE.OP, opType: OP.DIV, symbol: OP.DIV },
      { blockType: BLOCK_TYPE.OP, opType: OP.LESS, symbol: OP.LESS },
      { blockType: BLOCK_TYPE.OP, opType: OP.MORE, symbol: OP.MORE },
      { blockType: BLOCK_TYPE.OP, opType: OP.EQUALS, symbol: OP.EQUALS },
      { blockType: BLOCK_TYPE.OP, opType: OP.AND, symbol: OP.AND },
      { blockType: BLOCK_TYPE.OP, opType: OP.OR, symbol: OP.OR },
    ].map(
      ({ blockType, opType, symbol, isOperand }, i) =>
        new Block({
          id: `I${i}`,
          blockType,
          opType,
          symbol,
          inputs: [0],
          outputs: isOperand ? [] : [1, 3],
          w: BLOCK_W,
          h: BLOCK_H,
          shape: SHAPE_ROUND,
        })
    )
  ),
  chart: new Chart(
    STORED
      ? {
          blocks: STORED.blocks.map((b) => new Block({ ...b })),
          links: STORED.links.map((l) => new Link({ ...l })),
          lastId: Math.max(
            ...STORED.blocks.map(({ id }) => parseInt(id.slice(1), 10))
          ),
        }
      : { blocks: [], links: [] }
  ),
};

export default INITIAL_STATE;
