import React, { ReactNode } from 'react';
import { IconName } from '@blueprintjs/icons';
import { Button, Menu as BPMenu, MenuItem } from '@blueprintjs/core';
import styles from './Menu.module.scss';

interface MenuItemProps {
  icon?: IconName;
  text: string;
  onClick?: () => void;
}

const MenuItemComponent = ({ icon, text, onClick }: MenuItemProps) => {
  return (
    <MenuItem icon={icon} text={text} onClick={onClick} />
  );
};

interface MenuProps {
  children: ReactNode;
}

const Menu = ({ children }: MenuProps) => {
  return (
    <BPMenu className={styles.Menu}>
      {React.Children.map(children, (child, index) => (
        <div key={index}>
          {child}
        </div>
      ))}
    </BPMenu>
  );
};

export default Menu;
export { MenuItemComponent };
