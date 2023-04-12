import { createContext, useContext } from 'react';

type SourceContextType = {
  source: string;
  setSource: (source: string) => void;
  importSource: (file: File) => Promise<void>;
  sourceName: string;
  setSourceName: (sourceName: string) => void;
  isValidXYZFormat: boolean;
  setIsValidXYZFormat: (isValidXYZFormat: boolean) => void;
};

const SourceContext = createContext<SourceContextType>({
  source: '',
  setSource: () => {},
  importSource: async () => {},
  sourceName: '',
  setSourceName: () => {},
  isValidXYZFormat: false,
  setIsValidXYZFormat: () => {}
});

export const useSourceContext = () => useContext(SourceContext);

export default SourceContext;
