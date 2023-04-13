import React from 'react';
import { MenuItem } from '@blueprintjs/core';
import Menu from '../menu/Menu';

interface FileMenuProps {
  openDialog: () => void;
}

export const FileMenu: React.FC<FileMenuProps> = ({ openDialog }) => {
  return (
    <Menu>
      <MenuItem icon="insert" text="Import file" onClick={openDialog} />
      <MenuItem icon="import" text="Export file" />
    </Menu>
  );
};
