import { ReactNode } from 'react';
import './Utils.scss';
interface StackProps {
  children: ReactNode;
}

const HStack = ({ children }: StackProps) => {
  return (
    <div className='HStack'>
          {children}
    </div>
  );
};

export default HStack;
