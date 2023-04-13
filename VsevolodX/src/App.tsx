import './App.scss';
import Toolbar from './components/toolbar/Toolbar';
import SourceEditor from './views/sourceEditor/SourceEditor'
import StructureViewer from './views/structureViewer/StructureViewer';
import SettingsEditor from './views/settingsEditor/SettingsEditor';
import { useSettings } from './context/SettingsContext';
import { SourceProvider } from './context/SourceContext';
import { useState } from 'react';
import HStack from './components/utils/HStack';
import MainView from './views/main_view/MainView';

function App() {
  const settings = useSettings();
  const theme = settings.settings.theme;
  const [settingsShown, setSettingsShown] = useState(false);
  const toggleSettings = () => setSettingsShown(!settingsShown);
  
  return (
    <>
      <SourceProvider>
        <div className={`App bp4-${theme}`}>
          <header className="App-header">
            <h2>For Mat3ra. Materials Designer PoC</h2>
          </header>
          <main className="App-main">
            <Toolbar toggleSettings={toggleSettings}/>
            <MainView>
              <SourceEditor />
              <StructureViewer />
              {settingsShown && <SettingsEditor />}
            </MainView>
        </main>
          <footer className='App-footer'>
            <p>Version: 0.1.1</p>
          </footer>
        </div>
      </SourceProvider>
    </>
  );
}

export default App;
