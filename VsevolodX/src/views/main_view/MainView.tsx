import { ReactNode } from 'react';
import './MainView.scss';

interface MainViewProps {
  children: ReactNode;
}

const MainView = ({ children }: MainViewProps) => {
  return (
    <div className='MainView'>
          {children}
    </div>
  );
};

export default MainView;
