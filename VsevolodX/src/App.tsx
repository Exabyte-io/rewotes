import React, { useEffect, useState } from 'react';

import './App.scss';
import Toolbar from './components/Toolbar';
import SourceEditor from './components/SourceEditor';
import StructureViewer from './components/StructureViewer';
import SettingsEditor from './components/SettingsEditor';
import { useSettings } from './components/SettingsContext';
import SourceProvider from './components/SourceProvider';


function App() {
  const settings = useSettings();
  const theme = settings.settings.theme;

  const [version, setVersion] = useState(null);

  useEffect(() => {
    fetch('../package.json')
      .then(response => response.json())
      .then(data => {
        setVersion(data.version);
      });
  }, []);

  return (
    <SourceProvider>
      <div className={`App bp4-${theme}`}>
        <header className="App-header">
      </header>

      <main className="App-main">
        <Toolbar />
        <div className="HStack">
          <SourceEditor />
          <StructureViewer />
          <SettingsEditor />
        </div>

      </main>
      <footer className='App-footer'>
          {version ? (
          <p>Version: {version}</p>
        ) : (
          <p>Loading version...</p>
        )}
      </footer>
        <SourceEditor />
      </div>
    </SourceProvider>
  );
}

export default App;
