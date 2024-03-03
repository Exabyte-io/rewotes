import styles from './app.module.scss';

import { Editor } from '../xyz/Editor';
import View3d from '../view3d/View3d';

export function App() {
  return (
    <div className={styles['views-container']}>
        <Editor />
        <button></button>
        <View3d />
    </div>
  );
}

export default App;
