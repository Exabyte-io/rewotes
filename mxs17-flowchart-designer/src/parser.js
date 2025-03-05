import { BLOCK_TYPE, OP } from "./constants";

export function parseChart({ blocks, links }) {
  if (blocks.length < 3 && links.length < 2) return [null, {}];

  const outgoingLinksMap = links.reduce((acc, lnk) => {
    const { start, end } = lnk;
    if (!acc[start.blockId]) {
      acc[start.blockId] = {};
    }
    if (start.connIdx === 3) {
      acc[start.blockId].left = end.blockId;
    } else if (start.connIdx === 1) {
      acc[start.blockId].right = end.blockId;
    }
    return acc;
  }, {});

  const linkedBlocks = blocks.map(({ id, blockType, opType, value }) => ({
    id,
    blockType,
    opType,
    value,
  }));

  linkedBlocks.forEach((b, _, arr) => {
    const bLinks = outgoingLinksMap[b.id];
    if (!bLinks) return b;

    const leftChild = bLinks.left
      ? arr.find(({ id }) => id === bLinks.left)
      : null;
    if (leftChild) {
      b.left = leftChild;
      leftChild.parent = b;
    }

    const rightChild = bLinks.right
      ? arr.find(({ id }) => id === bLinks.right)
      : null;
    if (rightChild) {
      b.right = rightChild;
      rightChild.parent = b;
    }
  });

  const orphans = findOrphans(linkedBlocks);
  if (orphans.length) {
    return [new Error("Invalid tree: orphaned blocks were found")];
  }

  const [rootErr, root] = findRoot(linkedBlocks);
  if (rootErr) {
    return [rootErr];
  }

  const ast = traverse(root, linkedBlocks);
  return [null, ast];
}

function traverse(startNode, linkedNodesList) {
  const { blockType: type, opType: operation, value, left, right } = startNode;

  return {
    type,
    ...(operation ? { operation } : {}),
    ...(value ? { value } : {}),
    ...(left ? { left: traverse(left, linkedNodesList) } : {}),
    ...(right ? { right: traverse(right, linkedNodesList) } : {}),
  };
}

function findOrphans(linkedBlocks) {
  return linkedBlocks.filter(
    (b) =>
      (b.blockType === BLOCK_TYPE.OP && !(b.left || b.right)) ||
      (b.blockType !== BLOCK_TYPE.OP && !b.parent)
  );
}

function findRoot(linkedBlocks) {
  const roots = linkedBlocks.filter(
    (b) => b.blockType === BLOCK_TYPE.OP && !b.parent && b.left && b.right
  );

  if (!roots.length) {
    return [new Error("Invalid tree: no eligible root elements")];
  }
  if (roots.length > 1) {
    return [new Error("Invalid tree: more than one eligible root elements")];
  }
  return [null, roots[0]];
}

export function evalAst(ast) {
  const { type, operation: op, left, right, value } = ast;

  if (type === BLOCK_TYPE.NUM) {
    return value;
  } else if (type === BLOCK_TYPE.BOOL) {
    return value === "T";
  } else {
    return {
      [OP.PLUS]: evalAst(left) + evalAst(right),
      [OP.MINUS]: evalAst(left) - evalAst(right),
      [OP.MULT]: evalAst(left) * evalAst(right),
      [OP.DIV]: evalAst(left) / evalAst(right),
      [OP.LESS]: evalAst(left) < evalAst(right),
      [OP.MORE]: evalAst(left) > evalAst(right),
      [OP.EQUALS]: evalAst(left) === evalAst(right),
      [OP.AND]: evalAst(left) && evalAst(right),
      [OP.OR]: evalAst(left) || evalAst(right),
    }[op];
  }
}
