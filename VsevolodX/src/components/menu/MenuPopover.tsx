import React, { ReactNode } from 'react';
import { Popover2 } from '@blueprintjs/popover2';
import styles from './MenuPopover.module.scss';

interface MenuPopoverProps {
  isOpen: boolean;
  onInteraction: (willOpen: boolean) => void;
  position?: 'auto' | 'top' | 'bottom-left' | 'left' | 'right';
  content: string | JSX.Element | undefined;
  children: ReactNode;
}

export const MenuPopover: React.FC<MenuPopoverProps> = ({
  isOpen,
  onInteraction,
  position = 'bottom-left',
  content,
  children,
}) => {
  return (
    <Popover2 className={styles.MenuPopover}
      isOpen={isOpen}
      onInteraction={onInteraction}
      position={position}
      content={content}
      transitionDuration={100}
    >
      {children}
    </Popover2>
  );
};

export default MenuPopover;
