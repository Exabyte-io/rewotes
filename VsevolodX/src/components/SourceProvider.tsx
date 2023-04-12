import React, { useState, useEffect, ReactNode } from 'react';
import SourceContext from './SourceContext';

interface SourceProviderProps {
  children: ReactNode;
}

const SourceProvider: React.FC<SourceProviderProps> = ({ children }) => {
  const [source, setSource] = useState<string>('');
  const [sourceName, setSourceName] = useState<string>('');
  const [isValidXYZFormat, setIsValidXYZFormat] = useState(true);
  
  const loadSourceFromLocalStorage = () => {
    const storedSource = localStorage.getItem('sourceText');
    if (storedSource) {
      setSource(storedSource);
    }
  };

  const saveSourceToLocalStorage = (newSource: string) => {
    localStorage.setItem('sourceText', newSource);
  };

  const loadSourceNameFromLocalStorage = () => {
    const storedSourceName = localStorage.getItem('sourceName');
    if (storedSourceName) {
      setSourceName(storedSourceName);
    }
  };

  const saveSourceNameToLocalStorage = (newSourceName: string) => {
    localStorage.setItem('sourceName', newSourceName);
  };

  useEffect(() => {
    loadSourceFromLocalStorage();
    loadSourceNameFromLocalStorage();
  }, []);

  useEffect(() => {
    saveSourceToLocalStorage(source);
    saveSourceNameToLocalStorage(sourceName)
  }, [source, sourceName]);

  const importSource = async (file: File) => {
    const fileReader = new FileReader();
    fileReader.onload = (e) => {
      if (typeof e.target?.result === 'string') {
        setSource(e.target.result);
      }
    };
    fileReader.readAsText(file);
    setSourceName(file.name);
  };

  return (
    <SourceContext.Provider value={{ source, setSource, importSource, sourceName, setSourceName, isValidXYZFormat, setIsValidXYZFormat  }}>
      {children}
    </SourceContext.Provider>
  );
};

export default SourceProvider;
