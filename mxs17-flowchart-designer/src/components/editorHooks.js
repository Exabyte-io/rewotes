import { useCallback, useEffect, useRef, useState } from "react";

import Link from "../chart-elements/Link";
import { BLOCK_TYPE } from "../constants";
import { drawBezierLink } from "../helpers";

export default function useEditorHooks(panels, updatePanels) {
  const [editedBlock, setEditedBlock] = useState(null);

  const canvasRef = useRef();
  const contextRef = useRef();
  const created = useRef();
  const selected = useRef();
  const dragStart = useRef();
  const activeLine = useRef();

  const draw = useCallback(() => {
    const cnv = canvasRef.current;
    const ctx = contextRef.current;
    ctx.clearRect(0, 0, cnv.clientWidth, cnv.clientHeight);

    panels.palette.draw({ cnv, ctx });
    panels.chart.draw({ cnv, ctx });

    if (created.current) {
      created.current.draw(ctx);
    }

    ctx.strokeStyle = "gray";
    ctx.strokeRect(0, 0, cnv.clientWidth, cnv.clientHeight);

    if (activeLine.current) {
      const { x0, y0, x1, y1, type } = activeLine.current.connection;
      drawBezierLink(ctx, x0, y0, x1, y1, type === "output");
    }
  }, [panels]);

  const getSelected = useCallback(
    (x, y) => {
      return [...panels.palette.icons, ...panels.chart.blocks].find((b) =>
        b.mouseIsOver(x, y)
      );
    },
    [panels]
  );

  const getSelectedConnection = useCallback(
    (x, y) => {
      for (let b of panels.chart.blocks) {
        const conn = b.getConnectionUnderMouse(x, y);
        if (conn) {
          return { block: b, connection: conn };
        }
      }

      return null;
    },
    [panels]
  );

  const resetDragging = useCallback(() => {
    selected.current = null;
    created.current = null;
    dragStart.current = null;
    activeLine.current = null;
  }, []);

  const handleBlockValueInput = useCallback(
    (v) => {
      if (!editedBlock) return;
      editedBlock.value = v;
      setEditedBlock(null);
      updatePanels({ ...panels });
    },
    [editedBlock, setEditedBlock, panels, updatePanels]
  );

  const handleDoubleClick = useCallback(
    (e) => {
      const cnv = canvasRef.current;
      const mouseX = e.nativeEvent.offsetX - cnv.clientLeft;
      const mouseY = e.nativeEvent.offsetY - cnv.clientTop;
      const sb = getSelected(mouseX, mouseY);
      if (sb.blockType === BLOCK_TYPE.NUM || sb.blockType === BLOCK_TYPE.BOOL) {
        setEditedBlock(sb);
      }
    },
    [getSelected]
  );

  const handleMouseDown = useCallback(
    (e) => {
      const cnv = canvasRef.current;
      dragStart.current = {};
      dragStart.current.x = e.nativeEvent.offsetX - cnv.clientLeft;
      dragStart.current.y = e.nativeEvent.offsetY - cnv.clientTop;

      const sb = getSelected(dragStart.current.x, dragStart.current.y);

      if (sb?.isIcon) {
        created.current = sb.clone({ id: "new" });
        selected.current = created.current;
        return;
      }
      if (sb) {
        selected.current = sb;
        return;
      }

      const conn = getSelectedConnection(
        dragStart.current.x,
        dragStart.current.y
      );
      if (conn) {
        activeLine.current = conn;
        return;
      }

      resetDragging();
    },
    [resetDragging, getSelected, getSelectedConnection]
  );

  const handleMouseMove = useCallback(
    (e) => {
      if (!(selected.current || dragStart.current || activeLine.current))
        return;

      const cnv = canvasRef.current;
      const mouseX = e.nativeEvent.offsetX - cnv.clientLeft;
      const mouseY = e.nativeEvent.offsetY - cnv.clientTop;

      if (activeLine.current) {
        activeLine.current.connection.x1 = mouseX;
        activeLine.current.connection.y1 = mouseY;
        return draw();
      }

      const dx = mouseX - dragStart.current.x;
      const dy = mouseY - dragStart.current.y;
      dragStart.current.x = mouseX;
      dragStart.current.y = mouseY;
      selected.current.x += dx;
      selected.current.y += dy;
      draw();
    },
    [draw]
  );

  const handleMouseUp = useCallback(
    (e) => {
      const cnv = canvasRef.current;
      const mouseX = e.nativeEvent.offsetX - cnv.clientLeft;
      const mouseY = e.nativeEvent.offsetY - cnv.clientTop;

      if (created.current && panels.chart.mouseIsOver(mouseX)) {
        created.current.x = selected.current.x;
        created.current.y = selected.current.y;
        panels.chart.addBlock(created.current);
        updatePanels({ ...panels });
      } else if (
        !created.current &&
        selected.current &&
        panels.palette.mouseIsOver(mouseX)
      ) {
        const delId = selected.current.id;
        panels.chart.blocks = panels.chart.blocks.filter(
          ({ id }) => id !== delId
        );
        panels.chart.links = panels.chart.links.filter(
          ({ start, end }) => start.blockId !== delId && end.blockId !== delId
        );
        updatePanels({ ...panels });
      } else if (!created.current && selected.current) {
        const { id, x, y } = selected.current;
        const selBlock = panels.chart.blocks.find((b) => b.id === id);
        selBlock.x = x;
        selBlock.y = y;
        updatePanels({ ...panels });
      } else if (activeLine.current) {
        const { block: startBlock, connection: startConn } = activeLine.current;
        const endConnData = getSelectedConnection(mouseX, mouseY);

        if (!endConnData) {
          resetDragging();
          return draw();
        }

        const { block: endBlock, connection: endConn } = endConnData;
        if (startBlock.id === endBlock.id) return;

        panels.chart.addLink(
          new Link({
            start: {
              blockId: startBlock.id,
              connIdx: startConn.index,
              connType: startConn.type,
            },
            end: {
              blockId: endBlock.id,
              connIdx: endConn.index,
              connType: endConn.type,
            },
          })
        );
        updatePanels({ ...panels });
      }

      resetDragging();
      return draw();
    },
    [draw, panels, updatePanels, getSelectedConnection, resetDragging]
  );

  useEffect(() => {
    const canvasElement = canvasRef.current;
    contextRef.current = canvasElement.getContext("2d");
    draw();
  }, [draw, panels]);

  return {
    canvasRef,
    handleMouseDown,
    handleMouseUp,
    handleMouseMove,
    handleDoubleClick,
    editedBlock,
    handleBlockValueInput,
  };
}
