import React, { ReactNode } from 'react';
import { Popover2 } from '@blueprintjs/popover2';

interface MyPopoverProps {
  isOpen: boolean;
  onInteraction: (willOpen: boolean) => void;
  position?: 'auto' | 'top' | 'bottom' | 'left' | 'right';
  content: string | JSX.Element | undefined;
  children: ReactNode;
}

export const MyPopover: React.FC<MyPopoverProps> = ({
  isOpen,
  onInteraction,
  position = 'bottom',
  content,
  children,
}) => {
  return (
    <Popover2
      isOpen={isOpen}
      onInteraction={onInteraction}
      position={position}
      hoverCloseDelay={300}
      content={content}
      transitionDuration={100}
    >
      {children}
    </Popover2>
  );
};

export default MyPopover;
