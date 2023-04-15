import { createContext, useContext, useState, useEffect, ReactNode } from 'react';

type SourceContextType = {
  source: string;
  setSource: (source: string) => void;
  importSource: (file: File) => Promise<void>;
  sourceName: string;
  setSourceName: (sourceName: string) => void;
  saveSourceToLocalStorage: (source: string) => void;
  isValidXYZFormat: boolean;
  setIsValidXYZFormat: (isValidXYZFormat: boolean) => void;
};

const SourceContext = createContext<SourceContextType>({
  source: '',
  setSource: () => {},
  importSource: async () => {},
  sourceName: '',
  setSourceName: () => {},
  saveSourceToLocalStorage: () => {},
  isValidXYZFormat: false,
  setIsValidXYZFormat: () => {}
});

export const useSourceContext = () => useContext(SourceContext);

export default SourceContext;


interface SourceProviderProps {
  children: ReactNode;
  initialSource?: string;
}

export const SourceProvider: React.FC<SourceProviderProps> = ({ children, initialSource = ''}) => {
  const [source, setSource] = useState<string>(initialSource);
  const [sourceName, setSourceName] = useState<string>('');
  const [isValidXYZFormat, setIsValidXYZFormat] = useState(false);

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
    <SourceContext.Provider value={{ source, setSource, importSource, sourceName, setSourceName, saveSourceToLocalStorage, isValidXYZFormat, setIsValidXYZFormat  }}>
      {children}
    </SourceContext.Provider>
  );
};
