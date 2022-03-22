import BlockInput from "./BlockInput";

import {
  BLOCK_H,
  BLOCK_W,
  EDITOR_H,
  EDITOR_W,
  INPUT_H,
  INPUT_W,
} from "../constants";

import useEditorHooks from "./editorHooks";

export default function Editor({ panels, updatePanels }) {
  const {
    canvasRef,
    handleMouseDown,
    handleMouseUp,
    handleMouseMove,
    handleDoubleClick,
    editedBlock,
    handleBlockValueInput,
  } = useEditorHooks(panels, updatePanels);

  return (
    <div className="editor-container">
      <canvas
        ref={canvasRef}
        width={EDITOR_W}
        height={EDITOR_H}
        onMouseDown={handleMouseDown}
        onMouseUp={handleMouseUp}
        onMouseMove={handleMouseMove}
        onDoubleClick={handleDoubleClick}
      ></canvas>

      {editedBlock && (
        <BlockInput
          top={editedBlock.y + (BLOCK_H - INPUT_H) / 2 - 2}
          left={editedBlock.x + (BLOCK_W - INPUT_W) / 2 - 2}
          width={`${INPUT_W}px`}
          height={`${INPUT_H}px`}
          submitValue={handleBlockValueInput}
        />
      )}
    </div>
  );
}
