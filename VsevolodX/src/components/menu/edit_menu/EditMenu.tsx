import React from 'react';
import { MenuItem } from '@blueprintjs/core';
import Menu from '../Menu';
interface EditMenuProps {
  handleUndo: () => void;
  handleRedo: () => void;
}

export const EditMenu: React.FC<EditMenuProps> = ({ handleUndo, handleRedo }) => {
  return (
    <Menu>
      <MenuItem icon="undo" text="Undo" onClick={handleUndo} />
      <MenuItem icon="redo" text="Redo" onClick={handleRedo} />
      {/* Add other Edit menu items here */}
    </Menu>
  );
};
