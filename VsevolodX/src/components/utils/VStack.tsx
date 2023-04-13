import { ReactNode } from 'react';
import './Stack.scss';
interface StackProps {
  children: ReactNode;
}

const VStack = ({ children }: StackProps) => {
  return (
    <div className='VStack'>
          {children}
    </div>
  );
};

export default VStack;
