import React from 'react';
import { MenuItem } from '@blueprintjs/core';
import Menu from '../menu/Menu';
export const ViewMenu: React.FC = () => {
  return (
    <Menu>
          {[...Array(5)].map((_, i) => (
          <>
          <MenuItem className='skeleton' icon='disable' text={`Option ${i}`} disabled={true}/>
          </>
        ))}
    </Menu>
  );
};
