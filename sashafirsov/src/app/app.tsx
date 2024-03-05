import styles from './app.module.scss';
import { ResizeObserver } from '@juggle/resize-observer'
import useMeasure from 'react-use-measure';

import { Editor } from '../xyz/Editor';
import View3d from '../view3d/View3d';

export function App() {
    const [ref, bounds] = useMeasure({ polyfill: ResizeObserver })
  return (
    <div className={styles['views-container']} ref={ref} >
        <Editor />
        <button data-testid='resize-slider'>⬌</button>
        <View3d />
    </div>
  );
}

export default App;
